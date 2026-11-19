#!/usr/bin/env python3
"""Paired EXP-065 interaction-strength sweep; global calibration, unfitted dynamics."""
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ThreadPoolExecutor
import argparse,copy,hashlib,importlib.metadata,json,platform,shutil,sys,tempfile,time
import numpy as np
import pandas as pd
import yaml
from scipy.optimize import differential_evolution
from scipy.special import expit
ROOT=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'external/Adapt-Env'))
from adaptenv import ThermodynamicConfig,ThermodynamicLandscape,PUBLIC_MODEL_SYSTEM
from production_landscapes_exp066 import cube_statistics,seed_for
from replicate_observation_exp065 import calibrate_noise,observe
from observation_calibration_exp065 import observation_bank,observed
import landscape_graph as graph
IDS={'GB1':'gb1_wu2016_four_site_graphfla','TrpB':'trpb_johnston2024_four_site_graphfla'}
COLS=['degree','beneficial','neutral','uniform_steps','greedy_steps','uniform_terminal_rank','greedy_terminal_rank','uniform_hit_top','greedy_hit_top','uniform_hit_best','greedy_hit_best','accessible_top','accessible_best','nondecreasing_top','nondecreasing_best','rank','uniform_terminal_fitness','greedy_terminal_fitness']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def predict(p,a,g,ratio):return p[0]+np.exp(p[1])*expit(p[2]+np.exp(p[3])*(a+ratio*g))

def fit_readout(ratio,mode,ca,cg,target,bank,origin,scale,floor,fit,seed,anchors,quick=False):
    """Original calibration objective; topology and epistasis never enter fitting."""
    probs=np.linspace(0,1,fit['quantiles']);q=np.quantile(target,probs)
    variance=float(target.var());mean_sd=float(target.std(axis=-1).mean())
    floor_affine=(floor-origin)/scale;mass=float((target<=floor_affine).mean())
    bounds=[fit['offset_bounds'],np.log(fit['scale_bounds']).tolist(),fit['position_bounds'],np.log(fit['steepness_bounds']).tolist()]
    trace=[]
    def parts(p):
        raw=origin+scale*predict(p,ca,cg,ratio)
        if bank is not None:r1,r2,_,_=observed(raw,bank);y=((r1+r2)/2-origin)/scale
        else:y=(raw-origin)/scale
        quantile=float(np.mean((np.quantile(y,probs)-q)**2)/variance)
        sd=float((y.std(axis=-1).mean()/mean_sd-1)**2)
        floor_loss=float(((y<=floor_affine).mean()-mass)**2)
        return dict(quantile=quantile,mean_sd=sd,floor_cdf=floor_loss,total=quantile+sd+floor_loss)
    def evaluate(p,candidate):
        losses=parts(p);trace.append(dict(evaluation=len(trace),candidate=candidate,**{f'p{i}':float(v) for i,v in enumerate(p)},**losses));return losses['total']
    if ratio in [0,4]:
        anchor=anchors[int(ratio)];p=np.array(anchor['parameters']);evaluate(p,'frozen_RUN172_anchor')
        np.testing.assert_allclose(parts(p)['total'],anchor['objective']['total'],rtol=1e-9,atol=1e-10)
        converged=anchor['optimizer_converged'];message=anchor['optimizer_message'];selected_source='frozen_RUN172_anchor'
    else:
        result=differential_evolution(lambda p:evaluate(p,'differential_evolution'),bounds,seed=seed,maxiter=2 if quick else fit['maxiter'],popsize=4 if quick else fit['popsize'],tol=1e-5,polish=False,workers=1)
        for anchor_ratio in [0,4]:evaluate(anchors[anchor_ratio]['parameters'],f'RUN172_r{anchor_ratio}_candidate')
        best=min(trace,key=lambda x:x['total']);p=np.array([best[f'p{i}'] for i in range(4)])
        selected_source=best['candidate'];converged=bool(result.success);message=str(result.message)
    setting=dict(ratio=float(ratio),observation=mode,parameters=p.tolist(),offset=origin+scale*float(p[0]),scale=scale*float(np.exp(p[1])),position=float(p[2]),steepness=float(np.exp(p[3])),objective=parts(p),evaluations=len(trace),optimizer_converged=converged,optimizer_message=message,selected_source=selected_source,bounds=bounds,relative_boundary_distance=[float(min(v-lo,hi-v)/(hi-lo)) for v,(lo,hi) in zip(p,bounds)],calibration_origin=origin,calibration_scale=scale,near_floor_raw=floor)
    return setting,pd.DataFrame(trace)

def graph_record(f,valid,eps,top,floor,base):
    out,comp,active=graph.graph(f,valid,eps,top,floor);v=out[valid]
    record=dict(**base,epsilon=float(eps),nodes=int(valid.sum()),edges=int(v[:,0].sum()/2),top_fraction=float((f[valid]>=top).mean()),best_fraction=float((f[valid]==f[valid].max()).mean()),above_floor_fraction=float((f[valid]>floor).mean()),no_beneficial_fraction=float((v[:,1]==0).mean()),neutral_edge_fraction=float(v[:,2].sum()/v[:,0].sum()))
    record.update({COLS[j]:float(v[:,j].mean()) for j in [3,4,5,6,7,8,9,10,11,12,13,14,16,17]})
    for name,cc,vv in [('all',comp,valid),('above_floor',active,valid&(f>floor))]:record.update({name+'_'+k:x for k,x in graph.component_summary(f,cc,vv,floor).items()})
    return record

def cube_records(values,metadata,relative,scale,floor,base):
    stats,mut,_=cube_statistics(values,metadata,np.array(relative)*scale)
    stats['near_floor_fraction']=np.tile((values<=floor).mean(axis=1),len(relative))
    rows=[]
    for rel in relative:
        s=stats[stats.resolution.eq(rel*scale)];m=mut[mut.resolution.eq(rel*scale)]
        d=dict(**base,relative_resolution=rel,n_cubes=len(s))
        for key in ['epistatic_fraction','ho_fraction','slope','beneficial_slope','sign_epistasis_fraction','reciprocal_sign_epistasis_fraction','sd_fitness','near_floor_fraction','power_available']:
            d[key]=float(s[key].mean()) if s[key].notna().any() else np.nan
        d['slope_median']=float(m.slope.median()) if m.slope.notna().any() else np.nan
        d['slope_availability']=float(m.eligible.mean());d['power_availability']=float(s.power_available.mean())
        den=s.total_power.sum();minimum=len(s)*max((rel*scale)**2,1e-24)
        for k in range(1,5):d[f'pooled_power_{k}']=float(s[f'power_{k}'].sum()/den) if den>minimum else np.nan
        d['pooled_epistatic_fraction']=sum(d[f'pooled_power_{k}'] for k in [2,3,4]);d['pooled_ho_fraction']=sum(d[f'pooled_power_{k}'] for k in [3,4])
        rows.append(d)
    return rows

def run(config,output,data,workers=4,quick=False):
    start=time.perf_counter();output.mkdir(parents=True,exist_ok=True)
    cfg=yaml.safe_load(config.read_text());shutil.copyfile(config,output/'config.yaml')
    inputs=json.loads((data/'manifest.json').read_text())
    for rel,digest in inputs['files'].items():assert sha(data/rel)==digest,rel
    shutil.copyfile(data/'manifest.json',output/'input_manifest.json')
    source_manifest=json.loads((config.parent/'interaction_strength_source.json').read_text())
    for rel,digest in source_manifest['files'].items():assert sha(ROOT/rel)==digest,rel
    shutil.copyfile(config.parent/'interaction_strength_source.json',output/'source_manifest.json')
    environment=dict(python=platform.python_version(),packages={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','pyarrow','PyYAML']})
    dump(output/'worker_environment.json',environment)
    graph.validate_kernel(output);graph.KERNEL=graph.kernel()
    oldcfg=yaml.safe_load((data/'config.yaml').read_text());oldsettings=json.loads((data/'selected_parameters.json').read_text())
    ratios=[0.,1.,4.] if quick else cfg['interaction_ratios'];relative=cfg['measurement']['relative_resolutions']
    fields=output/'latent_fields';fields.mkdir(exist_ok=True);cubesdir=output/'cube_values';cubesdir.mkdir(exist_ok=True)
    fitsdir=output/'fits';fitsdir.mkdir(exist_ok=True)
    allfits=[];graphs=[];spectra=[];latent_graphs=[];replays=[];coverage=[]
    for protein in cfg['proteins']:
        print('START',protein,flush=True);sub=data/protein;dataset=IDS[protein]
        settings=[s for s in oldsettings if s['dataset_id']==dataset and s['arm'] in ['r0_sigmoid','r4_sigmoid']]
        anchor=next(s for s in settings if s['arm']=='r4_sigmoid' and s['observation']=='none')
        origin=anchor['calibration_origin'];scale=anchor['calibration_scale'];floor=anchor['near_floor_raw']
        geometry=json.loads((sub/'geometry.json').read_text());table=pd.read_parquet(sub/'variants.parquet')
        observed_codes=np.column_stack([pd.Categorical(table[f'pos{j+1}'].astype(str),categories=geometry['alphabets'][j]).codes for j in range(4)])
        assert (observed_codes>=0).all();observed_ix=np.ravel_multi_index(observed_codes.T,(20,)*4);assert len(np.unique(observed_ix))==len(observed_ix)
        codes=np.array(np.unravel_index(np.arange(20**4),(20,)*4)).T.astype(np.uint8)
        f=np.full(len(codes),np.nan);f[observed_ix]=table.fitness.to_numpy(float);valid=np.isfinite(f);assert valid.sum()==geometry['unique_genotypes']
        np.savez_compressed(output/f'{protein}_empirical.npz',fitness=f,valid=valid,genotypes=codes)
        panel=np.load(sub/'empirical_panel.npz');index=pd.read_csv(sub/'empirical_panel_index.csv');calmask=index.split.eq('calibration').to_numpy();confmask=~calmask
        cal=panel['fitness'][calmask];conf=panel['fitness'][confmask];conf_genotypes=panel['genotypes'][confmask]
        np.savez_compressed(output/f'{protein}_empirical_cubes.npz',fitness=conf,genotypes=conf_genotypes)
        cl=np.load(sub/'latent_calibration.npz');ca,cg=cl['a'],cl['g']
        noise=None;banks={'none':None}
        if protein=='TrpB':
            nd=output/'noise';nd.mkdir(exist_ok=True)
            noise,_=calibrate_noise(pd.read_csv(sub/'donor_assignments.csv.gz'),set(json.loads((sub/'excluded_noise_donors.json').read_text())),oldcfg,nd)
            old_noise=json.loads((sub/'noise_model.json').read_text());new_noise=json.loads((nd/'noise_model.json').read_text())
            for k in ['depth_breaks','fitness_breaks','shared_covariance_raw']:np.testing.assert_allclose(new_noise[k],old_noise[k],rtol=0,atol=1e-14)
            uniforms=np.load(sub/'fit_observation_uniforms.npz')['uniforms']
            for mode in cfg['observation_modes'][protein][1:]:banks[mode]=observation_bank(noise,uniforms,mode=='shared_stop_sensitivity')
        fitted=[]
        for mi,mode in enumerate(cfg['observation_modes'][protein]):
            anchors={int(s['interaction_ratio']):s for s in settings if s['observation']==mode}
            for ri,ratio in enumerate(ratios):
                setting,trace=fit_readout(ratio,mode,ca,cg,(cal-origin)/scale,banks[mode],origin,scale,floor,cfg['fit'],seed_for(dataset)+cfg['fit']['seed']+1000+mi*100+ri,anchors,quick)
                setting.update(protein=protein,dataset_id=dataset);fitted.append(setting);allfits.append(setting)
                trace.to_csv(fitsdir/f'{protein}_{mode}_r{ratio:g}.csv.gz',index=False)
                print('FIT',protein,mode,ratio,round(setting['objective']['total'],6),setting['selected_source'],flush=True)
        dump(output/'selected_parameters.json',allfits)
        base=dict(protein=protein,seed=-1,cohort='empirical',ratio=-1.,observation='measured',layer='empirical')
        for rel in relative:graphs.append(graph_record(np.nan_to_num(f),valid,rel*scale,float(np.quantile(f[valid],.99)),floor,dict(**base,relative_resolution=rel)))
        spectra+=cube_records(conf,pd.DataFrame({'cube':np.arange(len(conf))}),relative,scale,floor,base)
        coverage.append(dict(protein=protein,nodes=int(valid.sum()),possible_nodes=len(f),origin=origin,scale=scale,floor=floor,calibration_cubes=len(cal),confirmation_cubes=len(conf)))
        oldindex=pd.read_csv(sub/'latent_confirmation_index.csv');oldlatent=np.load(sub/'latent_confirmation.npz')
        oldseeds=sorted(oldindex.seed.unique().tolist());fresh=[seed_for(dataset)+n for n in cfg['seed_design']['fresh_seed_offsets']]
        assert not set(oldseeds)&set(fresh)
        calseeds=set(pd.read_csv(sub/'latent_calibration_index.csv').seed);assert not calseeds&set(fresh)
        if quick:oldseeds=oldseeds[:1];fresh=fresh[:1]
        tasks=[(seed,'reused',i) for i,seed in enumerate(oldseeds)]+[(seed,'fresh',i) for i,seed in enumerate(fresh)]
        def one_seed(task):
            seed,cohort,seed_index=task;model=ThermodynamicLandscape(ThermodynamicConfig(length=4,seed=int(seed),**cfg['model']))
            a,g=model.latent_components(codes)
            np.savez_compressed(fields/f'{protein}_{seed}.npz',a=a,g=g)
            if cohort=='reused':
                selected=np.flatnonzero(oldindex.seed.eq(seed));cgtypes=oldlatent['genotypes'][selected]
                ci=np.ravel_multi_index(cgtypes.reshape(-1,4).T,(20,)*4)
                np.testing.assert_allclose(a[ci],oldlatent['a'][selected].ravel(),rtol=1e-10,atol=1e-10)
                np.testing.assert_allclose(g[ci],oldlatent['g'][selected].ravel(),rtol=1e-10,atol=1e-10)
            else:
                selected=(np.arange(cfg['sampling']['synthetic_cubes_per_seed'])+seed_index*cfg['sampling']['synthetic_cubes_per_seed'])%len(conf_genotypes)
                cgtypes=conf_genotypes[selected];ci=np.ravel_multi_index(cgtypes.reshape(-1,4).T,(20,)*4)
            noise_seed=cfg['seed_design']['noise_seed_offset']+int(seed)
            uniforms=np.random.default_rng(noise_seed).random((len(codes),6))
            grow=[];srow=[];lrow=[];saved={'genotypes':cgtypes};audit=dict(protein=protein,seed=int(seed),cohort=cohort,noise_seed=noise_seed,replayed_nodes=len(ci) if cohort=='reused' else 0)
            for ratio in ratios:
                latent=a+ratio*g
                lb=dict(protein=protein,seed=int(seed),cohort=cohort,ratio=float(ratio),observation='none',layer='latent',relative_resolution=0.)
                lrow.append(graph_record(latent,valid,0.,float(np.quantile(latent[valid],.99)),float(latent[valid].min()-1),lb))
            for setting in fitted:
                ratio=setting['ratio'];mode=setting['observation'];raw=origin+scale*predict(setting['parameters'],a,g,ratio)
                public_cfg=ThermodynamicConfig(length=4,seed=int(seed),**cfg['model'],interaction_ratio=ratio,reference_occupancy=float(expit(setting['position'])),steepness=setting['steepness'],fitness_min=setting['offset'],fitness_max=setting['offset']+setting['scale'])
                np.testing.assert_allclose(ThermodynamicLandscape(public_cfg).evaluate_batch(codes[::211]),raw[::211],rtol=1e-10,atol=1e-10*scale)
                maps=[('intrinsic',raw)]
                if mode!='none':
                    r1,r2,_,_=observe(raw,noise,uniforms,shared=mode=='shared_stop_sensitivity');maps.append(('observed',(r1+r2)/2))
                for layer,y in maps:
                    b=dict(protein=protein,seed=int(seed),cohort=cohort,ratio=ratio,observation=mode,layer=layer)
                    for rel in relative:
                        rec=graph_record(y,valid,rel*scale,float(np.quantile(y[valid],.99)),floor,dict(**b,relative_resolution=rel))
                        rec['unique_fitness_fraction']=len(np.unique(y[valid]))/valid.sum();grow.append(rec)
                    cv=y[ci].reshape(cgtypes.shape[:2]);saved[f'{mode}_{layer}_r{ratio:g}']=cv
                    srow+=cube_records(cv,pd.DataFrame({'cube':np.arange(len(cv))}),relative,scale,floor,b)
            np.savez_compressed(cubesdir/f'{protein}_{seed}.npz',**saved)
            return grow,srow,lrow,audit
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for grow,srow,lrow,audit in pool.map(one_seed,tasks):
                graphs+=grow;spectra+=srow;latent_graphs+=lrow;replays.append(audit)
                pd.DataFrame(graphs).to_csv(output/'graph_summary.csv',index=False)
                print('GRAPHS',protein,audit['cohort'],audit['seed'],'complete',flush=True)
        pd.DataFrame(spectra).to_csv(output/'cube_summary.csv',index=False)
        pd.DataFrame(latent_graphs).to_csv(output/'latent_graph_summary.csv',index=False)
        pd.DataFrame(replays).to_csv(output/'seed_audit.csv',index=False)
        pd.DataFrame(coverage).to_csv(output/'coverage.csv',index=False)
    # The r=0/r=4 reused noiseless and observed graph outcomes must reproduce ANA-068;
    # this comparison is performed independently after collection, not a fitted target.
    files={str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    dump(output/'manifest.json',dict(protocol=cfg['protocol'],mode='SMOKE_ONLY' if quick else 'PAIRED_SWEEP',public_model=PUBLIC_MODEL_SYSTEM,elapsed_seconds=time.perf_counter()-start,input_manifest_sha256=sha(data/'manifest.json'),config_sha256=sha(config),source_manifest_sha256=sha(config.parent/'interaction_strength_source.json'),settings=len(allfits),graph_rows=len(graphs),latent_graph_rows=len(latent_graphs),cube_summary_rows=len(spectra),seed_audits=len(replays),files=files,interpretation=cfg['interpretation']))
    print('COMPLETE',json.dumps({'settings':len(allfits),'graph_rows':len(graphs),'seconds':time.perf_counter()-start}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--data-root',type=Path,required=True);p.add_argument('--workers',type=int,default=4);p.add_argument('--quick',action='store_true');a=p.parse_args()
    run(a.config,a.output_dir,a.data_root,a.workers,a.quick)
