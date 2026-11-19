"""ANA-068: exact graph diagnostics of frozen GB1/TrpB model fits."""
from pathlib import Path
import ctypes, hashlib, json, os, subprocess, sys, tempfile, time, functools, shutil, platform
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import rankdata
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import yaml
ANA=Path(__file__).resolve().parents[1];ROOT=ANA.parents[1]
RUN=ROOT/'data/processed/production_landscapes/RUN-172'
DATA=ROOT/'data/processed/production_landscapes/ANA-068'
TABLES=ANA/'tables';FIGURES=ANA/'figures'
for p in [DATA,DATA/'inputs',DATA/'landscapes',DATA/'walks',TABLES,FIGURES]:p.mkdir(parents=True,exist_ok=True)
SOURCE=RUN/'provenance/source';sys.dont_write_bytecode=True
sys.path.insert(0,str(SOURCE/'external/Adapt-Env'));sys.path.insert(0,str(SOURCE/'scripts'))
from adaptenv import ThermodynamicConfig,ThermodynamicLandscape
from replicate_observation_exp065 import calibrate_noise,observe
IDS={'GB1':'gb1_wu2016_four_site_graphfla','TrpB':'trpb_johnston2024_four_site_graphfla'}
RES=[0.,.001,.01];NW=2048;STEPS=128
COLS=['degree','beneficial','neutral','uniform_steps','greedy_steps','uniform_terminal_rank','greedy_terminal_rank','uniform_hit_top','greedy_hit_top','uniform_hit_best','greedy_hit_best','accessible_top','accessible_best','nondecreasing_top','nondecreasing_best','rank','uniform_terminal_fitness','greedy_terminal_fitness']
CASES=['empirical','r4_none','r0_none','r4_independent_intrinsic','r4_independent_observed','r4_shared_intrinsic','r4_shared_observed']
LABELS={'empirical':'Empirical','r4_none':'r=4, no added noise','r0_none':'r=0, no added noise','r4_independent_intrinsic':'Independent: intrinsic','r4_independent_observed':'Independent: observed','r4_shared_intrinsic':'Shared-stop: intrinsic','r4_shared_observed':'Shared-stop: observed'}
COLORS={'empirical':'#252525','r4_none':'#0072B2','r0_none':'#D55E00','r4_independent_intrinsic':'#56B4E9','r4_independent_observed':'#009E73','r4_shared_intrinsic':'#CC79A7','r4_shared_observed':'#7B3294'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,o):Path(p).write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')

def kernel():
    source=ANA/'notebooks/graph_kernel.cpp';key=sha(source)[:16]
    lib=Path(tempfile.gettempdir())/('ana068-'+key+('.dylib' if sys.platform=='darwin' else '.so'))
    if not lib.exists():
        cmd=['clang++' if sys.platform=='darwin' else 'g++','-O3','-std=c++17','-shared','-fPIC',str(source),'-o',str(lib)]
        if sys.platform=='darwin':
            sdk=subprocess.check_output(['xcrun','--show-sdk-path'],text=True).strip()
            cmd+=['-isysroot',sdk,'-isystem',str(Path(sdk)/'usr/include/c++/v1')]
        result=subprocess.run(cmd,capture_output=True,text=True)
        if result.returncode:raise RuntimeError('Graph-kernel compilation failed: '+result.stderr)
    dll=ctypes.CDLL(str(lib));D=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');B=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
    dll.analyse.argtypes=[ctypes.c_int]*3+[D,B]+[ctypes.c_double]*3+[D,I,I]
    dll.walks.argtypes=[ctypes.c_int]*3+[D,B,I]+[ctypes.c_int]*2+[ctypes.c_double,ctypes.c_int,ctypes.c_uint64,I]
    return dll
KERNEL=None

def graph(f,valid,eps,top,floor,q=20,dim=4):
    global KERNEL
    if KERNEL is None:KERNEL=kernel()
    f=np.ascontiguousarray(f,dtype=np.float64);valid=np.ascontiguousarray(valid,dtype=np.uint8)
    out=np.zeros((len(f),18));components=np.full(len(f),-1,dtype=np.int32);active=components.copy()
    assert KERNEL.analyse(len(f),q,dim,f,valid,float(eps),float(top),float(floor),out,components,active)==0
    return out,components,active

def sample_walks(f,valid,eps,starts,greedy,seed,q=20,dim=4,steps=STEPS):
    starts=np.ascontiguousarray(starts,dtype=np.int32);paths=np.zeros((len(starts),steps+1),dtype=np.int32)
    assert KERNEL.walks(len(f),q,dim,np.ascontiguousarray(f),np.ascontiguousarray(valid,dtype=np.uint8),starts,len(starts),steps,float(eps),int(greedy),int(seed),paths)==0
    return paths

def validate_kernel():
    """Independent recursive enumeration and sparse components on small graphs."""
    checks=0
    for q,dim in [(2,3),(3,2)]:
        n=q**dim;codes=np.array(np.unravel_index(np.arange(n),(q,)*dim)).T
        neighbours=[np.flatnonzero((codes!=codes[u]).sum(axis=1)==1) for u in range(n)]
        for trial in range(5):
            rng=np.random.default_rng(6800+trial);f=np.round(rng.normal(size=n),1);valid=np.ones(n,dtype=bool)
            if trial%2:valid[1]=False
            for eps in [0.,.15]:
                target=float(np.quantile(f[valid],.8));floor=0.;out,comp,ac=graph(f,valid,eps,target,floor,q,dim)
                ranks=np.zeros(n);ranks[valid]=(rankdata(f[valid],method='average')-.5)/valid.sum();best=max(f[valid])
                @functools.lru_cache(None)
                def solve(u,greedy):
                    vv=[v for v in neighbours[u] if valid[v] and f[v]-f[u]>eps]
                    if greedy and vv:vv=[v for v in vv if f[v]==max(f[w] for w in vv)]
                    if not vv:return np.array([0,ranks[u],float(f[u]>=target),float(f[u]==best),f[u]])
                    z=np.mean([solve(v,greedy) for v in vv],axis=0);z[0]+=1
                    if f[u]>=target:z[2]=1
                    if f[u]==best:z[3]=1
                    return z
                for u in np.flatnonzero(valid):
                    for greedy,cols in [(False,[3,5,7,9,16]),(True,[4,6,8,10,17])]:np.testing.assert_allclose(out[u,cols],solve(int(u),greedy),atol=1e-12)
                    np.testing.assert_equal(out[u,0],sum(valid[neighbours[u]]))
                    for tol,col in [(eps,11),(eps,12),(None,13),(None,14)]:
                        targets=set(np.flatnonzero(valid & (f>=target if col in [11,13] else f==best)));seen={int(u)};queue=[int(u)]
                        for x in queue:
                            for v in neighbours[x]:
                                if valid[v] and v not in seen and ((f[v]>=f[x]) if tol is None else (f[v]-f[x]>tol)):seen.add(int(v));queue.append(int(v))
                        assert out[u,col]==bool(seen&targets)
                for selection,cc in [(valid,comp),(valid&(f>floor),ac)]:
                    ix=np.flatnonzero(selection)
                    adj=np.array([[int(v in neighbours[u] and abs(f[v]-f[u])<=eps) for v in ix] for u in ix])
                    if len(ix):
                        _,labels=connected_components(csr_matrix(adj),directed=False)
                        np.testing.assert_array_equal(labels[:,None]==labels[None,:],cc[ix,None]==cc[ix][None,:])
                checks+=1
    result={'small_graph_cases':checks,'independent_checks':['recursive walk moments','exhaustive reachability','SciPy components','degrees'],'kernel_sha256':sha(ANA/'notebooks/graph_kernel.cpp')}
    dump(TABLES/'kernel_validation.json',result);return result

def prepare():
    manifest=json.loads((RUN/'manifest.json').read_text());inventory=pd.read_csv(RUN/'dataset_inventory.csv').set_index('dataset_id')
    fits=json.loads((RUN/'selected_parameters.json').read_text());audit={'inputs':{},'source':{},'frozen_fits_sha256':sha(RUN/'selected_parameters.json')}
    sm=json.loads((SOURCE/'experiments/2026-10-08_EXP-066_production-landscape-aggregate-realism/source_manifest.json').read_text())
    for relative,expected in sm['files'].items():assert sha(SOURCE/relative)==expected;audit['source'][relative]=expected
    arrays={};allcodes=np.array(np.unravel_index(np.arange(20**4),(20,)*4)).T.astype(np.uint8)
    for protein,dataset in IDS.items():
        original=ROOT/inventory.loc[dataset,'parquet_path'];p=DATA/'inputs'/f'{protein}.parquet'
        expected=inventory.loc[dataset,'sha256_parquet']
        if not p.exists():
            assert sha(original)==expected;shutil.copyfile(original,p)
        assert sha(p)==expected
        table=pd.read_parquet(p);geo=json.loads((RUN/dataset/'geometry.json').read_text())
        codes=np.column_stack([pd.Categorical(table[f'pos{j+1}'].astype(str),categories=geo['alphabets'][j]).codes for j in range(4)])
        assert (codes>=0).all();indices=np.ravel_multi_index(codes.T,(20,)*4);assert len(np.unique(indices))==len(indices)
        f=np.full(20**4,np.nan);f[indices]=table.fitness.to_numpy(float);valid=np.isfinite(f)
        assert valid.sum()==geo['unique_genotypes'];arrays[protein]=(f,valid,allcodes)
        np.savez_compressed(DATA/'inputs'/f'{protein}.npz',fitness=f,valid=valid,genotypes=allcodes)
        audit['inputs'][protein]={'parquet_path':str(original.relative_to(ROOT)),'bundled_input':str(p.relative_to(DATA)),'sha256':expected,'observed_nodes':int(valid.sum()),'possible_nodes':len(f),'input_sha_in_run':geo['input_sha256']}
    donor_path=RUN/IDS['TrpB']/'donor_assignments.csv.gz';assert sha(donor_path)==manifest['files'][str(donor_path.relative_to(RUN))]
    cfg=yaml.safe_load((RUN/'config.yaml').read_text());donors=pd.read_csv(donor_path);excluded=set(json.loads((donor_path.parent/'excluded_noise_donors.json').read_text()))
    with tempfile.TemporaryDirectory(prefix='ana068-noise-') as td:noise,_=calibrate_noise(donors,excluded,cfg,Path(td))
    audit['noise_donors_sha256']=sha(donor_path);audit['analysis_kernel_sha256']=sha(ANA/'notebooks/graph_kernel.cpp')
    audit['calculation_source_sha256']=sha(__file__);audit['protocol_sha256']=sha(ANA/'protocol.md')
    audit['python_version']=platform.python_version();audit['numpy_version']=np.__version__;audit['pandas_version']=pd.__version__
    dump(DATA/'input_audit.json',audit);dump(DATA/'frozen_settings.json',[s for s in fits if s['dataset_id'] in IDS.values()])
    return arrays,fits,noise,audit

def component_summary(f,components,valid,floor):
    selected=(components>=0);ids,counts=np.unique(components[selected],return_counts=True)
    if not len(ids):return {'n':0,'largest_fraction':None,'largest_span':None,'largest_below_floor_fraction':None,'largest_size':0,'mean_size_seen_by_node':None}
    biggest=ids[np.argmax(counts)];values=f[components==biggest]
    return {'n':len(counts),'largest_fraction':float(counts.max()/len(f[valid])),'largest_span':float(np.ptp(values)),
      'largest_below_floor_fraction':float((values<=floor).mean()),'largest_size':int(counts.max()),'mean_size_seen_by_node':float((counts.astype(float)**2).sum()/counts.sum())}

def calculate():
    start=time.perf_counter();validation=validate_kernel();arrays,fits,noise,audit=prepare()
    records=[];curves=[];walk_records=[];coverage=[];mask_rows=[];neutral_sizes=[];replay=0
    for protein,(emp,valid,codes) in arrays.items():
        dataset=IDS[protein];settings=[s for s in fits if s['dataset_id']==dataset and s['arm'] in ['r4_sigmoid','r0_sigmoid']]
        base=next(s for s in settings if s['arm']=='r4_sigmoid' and s['observation']=='none')
        scale=base['calibration_scale'];origin=base['calibration_origin'];floor=base['near_floor_raw']
        ix=pd.read_csv(RUN/dataset/'latent_confirmation_index.csv');seeds=sorted(ix.seed.unique().tolist())
        def consume(case,seed,f):
            nonlocal records,curves,walk_records,neutral_sizes
            f=np.ascontiguousarray(f);top=float(np.quantile(f[valid],.99));best=float(f[valid].max())
            landscape=DATA/'landscapes'/f'{protein}_{case}_{seed}.npz';np.savez_compressed(landscape,fitness=f)
            for rel in RES:
                out,comp,active=graph(f,valid,rel*scale,top,floor)
                vals=out[valid];r=dict(protein=protein,case=case,seed=int(seed),relative_resolution=rel,epsilon=rel*scale,
                  nodes=int(valid.sum()),edges=int(vals[:,0].sum()/2),top_threshold=top,top_fraction=float((f[valid]>=top).mean()),best_fraction=float((f[valid]==best).mean()),above_floor_fraction=float((f[valid]>floor).mean()),no_beneficial_fraction=float((vals[:,1]==0).mean()),
                  neutral_edge_fraction=float(vals[:,2].sum()/vals[:,0].sum()),mean_beneficial_fraction=float(np.divide(vals[:,1],vals[:,0],out=np.zeros(len(vals)),where=vals[:,0]>0).mean()))
                r.update({COLS[j]:float(vals[:,j].mean()) for j in [3,4,5,6,7,8,9,10,11,12,13,14,16,17]})
                for name,cc,vv in [('all',comp,valid),('above_floor',active,valid&(f>floor))]:
                    st=component_summary(f,cc,vv,floor);r.update({name+'_'+k:v for k,v in st.items()})
                    _,sizes=np.unique(cc[cc>=0],return_counts=True)
                    for bound in [1,2,5,10,30,100,300,1000,3000,10000,30000,100000]:
                        neutral_sizes.append(dict(protein=protein,case=case,seed=int(seed),relative_resolution=rel,subset=name,size_at_least=bound,node_fraction=float(sizes[sizes>=bound].sum()/max(1,sizes.sum()))))
                records.append(r)
                if rel==.001:
                    start_ids=np.random.default_rng(681001+int(seed)).choice(np.flatnonzero(valid),NW,replace=True).astype(np.int32)
                    for greedy in [False,True]:
                        policy='greedy' if greedy else 'uniform';paths=sample_walks(f,valid,rel*scale,start_ids,greedy,682001+int(seed))
                        assert (paths[:,-1]==paths[:,-2]).all(),'walk truncation: increase STEPS'
                        changes=paths[:,1:]!=paths[:,:-1];lengths=changes.sum(axis=1);ends=paths[:,-1]
                        # Monte Carlo walk means must agree with exact graph expectations,
                        # allowing sampling variation in both starts and walk choices.
                        col=4 if greedy else 3
                        assert abs(lengths.mean()-out[valid,col].mean()) < 6*lengths.std()/np.sqrt(NW)+.05
                        rank_sample=out[ends,15];rank_col=6 if greedy else 5
                        assert abs(rank_sample.mean()-out[valid,rank_col].mean()) < 6*rank_sample.std()/np.sqrt(NW)+.003
                        np.savez_compressed(DATA/'walks'/f'{protein}_{case}_{seed}_{policy}.npz',starts=start_ids,length=lengths,ends=ends,terminal_rank=rank_sample)
                        for t in range(31):
                            y=(f[paths[:,t]]-origin)/scale;rank=out[paths[:,t],15]
                            curves.append(dict(protein=protein,case=case,seed=int(seed),policy=policy,step=t,mean_fitness=float(y.mean()),q10=float(np.quantile(y,.1)),median=float(np.median(y)),q90=float(np.quantile(y,.9)),mean_rank=float(rank.mean()),continuing_fraction=float((lengths>t).mean())))
                        for rankcut in np.linspace(0,1,101):walk_records.append(dict(protein=protein,case=case,seed=int(seed),policy=policy,rank=rankcut,cdf=float((out[ends,15]<=rankcut).mean())))
                    if case=='r4_none':
                        full=np.ones(len(f),dtype=bool);o,_,_=graph(f,full,rel*scale,top,floor)
                        # Same target threshold and measured start nodes isolates missing routes/neighbours.
                        for metric in ['accessible_top','uniform_hit_top','uniform_steps','no_beneficial_fraction']:
                            if metric=='no_beneficial_fraction':masked=float((out[valid,1]==0).mean());fullval=float((o[valid,1]==0).mean())
                            else:col=COLS.index(metric);masked=float(out[valid,col].mean());fullval=float(o[valid,col].mean())
                            mask_rows.append(dict(protein=protein,seed=int(seed),metric=metric,masked=masked,full=fullval,difference=masked-fullval))
            print(protein,case,seed,'finished',flush=True)
        consume('empirical',-1,np.nan_to_num(emp,nan=0.))
        empout,_,_=graph(np.nan_to_num(emp,nan=0.),valid,0.,np.quantile(emp[valid],.99),floor)
        coverage.append(dict(protein=protein,nodes=int(valid.sum()),possible_nodes=len(emp),edges=int(empout[valid,0].sum()/2),complete_neighbourhood_fraction=float((empout[valid,0]==76).mean()),mean_degree=float(empout[valid,0].mean()),calibration_origin=origin,calibration_scale=scale,near_floor_cutoff=floor))
        for seed in seeds:
            cfg=dict(base['thermodynamic_config']);cfg['seed']=int(seed);model=ThermodynamicLandscape(ThermodynamicConfig(**cfg))
            a,g=model.latent_components(codes)
            # Replay exact stored confirmation cubes for this seed.
            latent=np.load(RUN/dataset/'latent_confirmation.npz');saved_index=np.flatnonzero(ix.seed.eq(seed))
            subset=latent['genotypes'][saved_index];idx=np.ravel_multi_index(subset.reshape(-1,4).T,(20,)*4)
            np.testing.assert_allclose(a[idx],latent['a'][saved_index].ravel(),rtol=1e-10,atol=1e-10)
            np.testing.assert_allclose(g[idx],latent['g'][saved_index].ravel(),rtol=1e-10,atol=1e-10);replay+=len(idx)
            uniforms=np.random.default_rng(683001+int(seed)).random((len(codes),6))
            for s in settings:
                if s['arm']=='r0_sigmoid' and s['observation']!='none':continue
                raw=s['offset']+s['scale']*expit(s['position']+s['steepness']*(a+s['interaction_ratio']*g))
                if s['observation']=='none':
                    case='r4_none' if s['arm']=='r4_sigmoid' else 'r0_none'
                    if case=='r4_none':np.testing.assert_allclose(model.evaluate_batch(codes[::211]),raw[::211],rtol=1e-10,atol=1e-10)
                    consume(case,seed,raw)
                else:
                    prefix='r4_independent' if s['observation']=='independent_difference' else 'r4_shared'
                    consume(prefix+'_intrinsic',seed,raw)
                    r1,r2,_,_=observe(raw,noise,uniforms,shared=s['observation']=='shared_stop_sensitivity')
                    consume(prefix+'_observed',seed,(r1+r2)/2)
            pd.DataFrame(records).to_csv(DATA/'graph_summary.csv',index=False)
    for name,rows in [('graph_summary',records),('trajectories',curves),('endpoint_cdf',walk_records),('coverage',coverage),('mask_sensitivity',mask_rows),('neutral_component_sizes',neutral_sizes)]:pd.DataFrame(rows).to_csv(DATA/(name+'.csv'),index=False)
    audit.update(elapsed_seconds=time.perf_counter()-start,landscapes=len(records)//len(RES),graph_diagnostics=len(records),replayed_latent_nodes=replay,kernel_validation=validation,walks_per_policy_landscape=NW,walk_step_cap=STEPS,all_walks_terminated=True)
    audit['manifest']={str(p.relative_to(DATA)):sha(p) for p in sorted(DATA.rglob('*')) if p.is_file() and p.name!='audit.json'}
    dump(DATA/'audit.json',audit);dump(TABLES/'audit.json',audit)
    return load()

def load():
    names=['graph_summary','trajectories','endpoint_cdf','coverage','mask_sensitivity','neutral_component_sizes']
    frames={name:pd.read_csv(DATA/(name+'.csv')) for name in names}
    for name,frame in frames.items():frame.to_csv(TABLES/(name+'.csv'),index=False)
    return frames
