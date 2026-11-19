"""ANA-065: independent numerical audit and figure-led descriptive analysis.

Assertions verify identities and provenance, never scientific success criteria.
No import or execution of the submitted runner or original external notebooks.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import math
import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ANALYSIS=Path(__file__).resolve().parents[1]
ROOT=ANALYSIS.parents[1]
BUNDLE=ROOT/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-168'
INPUT=ROOT/'data/processed/trpb-replicate-observation-exp065'
OUT=BUNDLE/'outputs'
TABLES,FIGURES=ANALYSIS/'tables',ANALYSIS/'figures'
METHODS=['average','cross','stop_adjusted']
MLABELS=['Replicate average','Cross-replicate','Stop-adjusted sensitivity']
MCOLORS=['#0072B2','#E69F00','#009E73']
MODES=['intrinsic','independent_difference','shared_stop_sensitivity']
NLABELS=['Intrinsic','Independent error','Shared-stop sensitivity']
NCOLORS=['#555555','#0072B2','#D55E00']
STRATA=['all','no_active','some_active','all_active']
SLABELS=['Full panel','No active vertices','Some active vertices','All active vertices']
RES=[1e-12,.001,.01]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})


def dump(path,obj):
    def clean(v):
        if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
        if isinstance(v,(list,tuple,np.ndarray)):return [clean(x) for x in v]
        if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
        if isinstance(v,np.integer):return int(v)
        if isinstance(v,np.bool_):return bool(v)
        return v
    Path(path).write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def mean_finite(x,axis=-1):
    ok=np.isfinite(x)
    return np.divide(np.where(ok,x,0).sum(axis),ok.sum(axis),out=np.full(ok.sum(axis).shape,np.nan),where=ok.sum(axis)>0)


def ratio(a,b,threshold):
    return np.divide(a,b,out=np.full(np.broadcast_shapes(np.shape(a),np.shape(b)),np.nan),where=np.asarray(b)>threshold)


def fwht(x):
    """Butterfly contrast transform, independent of the submitted matrix code."""
    a=np.array(x,float,copy=True)
    for stride in (1,2,4,8):
        view=a.reshape(len(a),-1,2,stride)
        lo,hi=view[:,:,0].copy(),view[:,:,1].copy()
        view[:,:,0]=lo+hi;view[:,:,1]=lo-hi
    return a/16


def components(a,b,q):
    c1,c2,cm=fwht(a),fwht(b),fwht((a+b)/2)
    order=np.array([i.bit_count() for i in range(16)])
    p=np.array([np.column_stack([(cm[:,order==k]**2).sum(1) for k in range(1,5)]),
                np.column_stack([(c1[:,order==k]*c2[:,order==k]).sum(1) for k in range(1,5)])]).transpose(1,0,2)
    corr=q.sum(1)[:,None]*np.array([math.comb(4,k) for k in range(1,5)])/256
    power=np.concatenate([p,(p[:,1]-corr)[:,None,:]],axis=1)
    num,den,corrs,effects,edges=[],[],[],[],[]
    for site in range(4):
        low=np.flatnonzero((np.arange(16)&(1<<site))==0);high=low+(1<<site)
        x1,x2=a[:,low],b[:,low];d1,d2=a[:,high]-x1,b[:,high]-x2
        xm,dm=(x1+x2)/2,(d1+d2)/2
        xc1,xc2=x1-x1.mean(1)[:,None],x2-x2.mean(1)[:,None]
        xcm=xm-xm.mean(1)[:,None]
        nc=((xc1*d2).sum(1)+(xc2*d1).sum(1))/16
        dc=(xc1*xc2).sum(1)/8
        correction=q[:,low].sum(1)*7/64
        num.append(np.column_stack([(xcm*dm).sum(1)/8,nc,nc+correction]))
        den.append(np.column_stack([(xcm*xcm).sum(1)/8,dc,dc-correction]))
        corrs.append(correction);effects.append(dm.mean(1));edges.append(dm)
    return dict(power=power,numerator=np.array(num).transpose(1,2,0),denominator=np.array(den).transpose(1,2,0),
                correction=np.array(corrs).T,effect=np.array(effects).T,edges=np.array(edges).transpose(1,0,2),
                mean=((a+b)/2).mean(1),sd=((a+b)/2).std(1))


def statistic_columns(c,res):
    threshold=max(res**2,1e-18)
    cols={'mean_fitness':c['mean'],'sd_fitness':c['sd'],'neutral_edge_fraction':(np.abs(c['edges'])<=res).mean((1,2))}
    for j,m in enumerate(METHODS):
        p=c['power'][:,j];t=p.sum(1)
        cols.update({f'{m}_power_{k+1}':p[:,k] for k in range(4)})
        cols[f'{m}_total_power']=t;cols[f'{m}_power_available']=t>threshold
        cols[f'{m}_ho_fraction']=ratio(p[:,2:].sum(1),t,threshold)
        cols[f'{m}_epistatic_fraction']=ratio(p[:,1:].sum(1),t,threshold)
        valid=c['denominator'][:,j]>threshold
        slopes=ratio(c['numerator'][:,j],c['denominator'][:,j],threshold)
        benef=c['effect']>res
        cols[f'{m}_slope']=mean_finite(slopes);cols[f'{m}_slope_n']=valid.sum(1)
        cols[f'{m}_beneficial_slope']=mean_finite(np.where(benef,slopes,np.nan))
        cols[f'{m}_beneficial_slope_n']=(valid&benef).sum(1)
    return cols


def close(a,b,rtol=2e-7,atol=2e-10):
    np.testing.assert_allclose(a,b,rtol=rtol,atol=atol,equal_nan=True)
    x,y=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float));both=np.isfinite(x)&np.isfinite(y)
    return float(np.max(np.abs(x[both]-y[both]))) if both.any() else 0.


def verify_components(z,c):
    for saved,new in [('power','power'),('slope_numerator','numerator'),('slope_denominator','denominator'),('slope_correction','correction'),('mean_effect','effect')]:
        close(z[saved],c[new],rtol=1e-10,atol=2e-12)
    for k,r in enumerate(RES):
        np.testing.assert_array_equal(z['eligible'][k],c['denominator']>max(r*r,1e-18))
        np.testing.assert_array_equal(z['beneficial'][k],c['effect']>r)


def audit_and_load():
    FIGURES.mkdir(exist_ok=True);TABLES.mkdir(exist_ok=True)
    for key,h in json.loads((BUNDLE/'remote_sha256.json').read_text()).items():assert sha(BUNDLE/key)==h,key
    manifest=json.loads((OUT/'manifest.json').read_text());assert manifest['mode']=='PRODUCTION'
    assert (BUNDLE/'logs/labproj-exit-code.txt').read_text().strip()=='0'
    for name,h in manifest['output_sha256'].items():assert sha(OUT/name)==h,name
    submitted=json.loads((BUNDLE/'submission-verification.json').read_text())
    for name,h in submitted['sha256'].items():assert sha(BUNDLE/'source'/name)==h,name
    assert sha(BUNDLE/'source/scripts/replicate_observation_exp065.py')==manifest['runner_sha256']
    assert sha(OUT/'config.yaml')==manifest['config_sha256']
    assert sha(INPUT/'input_manifest.json')==manifest['input_manifest_sha256']
    inp=json.loads((INPUT/'input_manifest.json').read_text())
    for n,h in inp['input_sha256'].items():assert sha(INPUT/n)==h,n
    for item in json.loads((INPUT/'source_manifest.json').read_text())['files']:assert sha(INPUT/item['file'])==item['sha256']
    assert sha(OUT/'frozen_parameters.json')==sha(INPUT/'prior/selected_parameters.json')
    cfg=yaml.safe_load((OUT/'config.yaml').read_text())
    offset,scale=manifest['raw_offset'],manifest['raw_scale']
    table=pd.read_csv(INPUT/'replicates.csv').sort_values('genotype').reset_index(drop=True)
    assignments=pd.read_csv(OUT/'donor_assignments.csv.gz')
    assert table.genotype.tolist()==assignments.genotype.tolist()
    for name in table.columns:
        if name in ('genotype','published_active'):assert table[name].tolist()==assignments[name].tolist()
        else:close(table[name],assignments[name],rtol=0,atol=1e-14)
    cutoffs=cfg['measurement']['activity_cutoffs']
    active=(table.replicate_1>cutoffs[0])&(table.replicate_2>cutoffs[1])
    np.testing.assert_array_equal(active,table.published_active)
    excluded=set(json.loads((INPUT/'excluded_donor_genotypes.json').read_text()))
    rng=np.random.default_rng(cfg['sampling']['donor_split_seed']);split=np.full(len(table),'excluded',dtype=object)
    for mask in (table.stop_count>0,(table.stop_count==0)&~table.genotype.isin(excluded)):
        ix=rng.permutation(table.index[mask]);n=int(.8*len(ix));split[ix[:n]]='training';split[ix[n:]]='holdout'
    np.testing.assert_array_equal(split,assignments.split)
    assert not (set(assignments.loc[(assignments.stop_count==0)&(assignments.split!='excluded'),'genotype'])&excluded)
    noise=json.loads((OUT/'noise_model.json').read_text())
    trainstop=assignments[(assignments.stop_count>0)&(assignments.split=='training')]
    traindonor=assignments[(assignments.stop_count==0)&(assignments.split=='training')]
    depth=np.quantile(trainstop.input_count,[.25,.5,.75]);close(depth,noise['depth_breaks'])
    cutoff=np.mean(cutoffs);fitness=np.array([0,cutoff,traindonor.loc[traindonor.published_mean>cutoff,'published_mean'].median()]);close(fitness,noise['fitness_breaks'])
    np.testing.assert_array_equal(assignments.depth_bin,np.searchsorted(depth,assignments.input_count,side='right'))
    np.testing.assert_array_equal(assignments.fitness_bin,np.searchsorted(fitness,assignments.published_mean,side='right'))
    close(assignments.difference,assignments.replicate_1-assignments.replicate_2)
    controls=pd.read_csv(OUT/'control_diagnostics.csv');cal=pd.read_csv(OUT/'noise_calibration.csv')
    q=[];common=[]
    for _,r in controls.iterrows():
        g=assignments[(assignments.stop_count>0)&(assignments.split==r.split)&(assignments.depth_bin==r.depth_bin)]
        a,b=g.replicate_1.to_numpy(),g.replicate_2.to_numpy();qc=((a-a.mean())*(b-b.mean())).mean()
        close([r.mean_1,r.mean_2,r.variance_1,r.variance_2,r.covariance,r.difference_variance],[a.mean(),b.mean(),a.var(),b.var(),qc,(a-b).var()])
        assert len(g)==r.n
        if r.split=='training':
            q.append(qc);v=(a+b)/2;v-=v.mean();common.append(v*np.sqrt(max(0,qc)/v.var()))
    close(q,noise['shared_covariance_raw']);q=np.array(q)
    donors,depthpools={},{}
    for f in range(4):
        fs=traindonor[traindonor.fitness_bin==f];depthpools[f]=fs.input_count.to_numpy()
        for d in range(4):
            cell=fs[fs.depth_bin==d];r=cal[(cal.fitness_bin==f)&(cal.depth_bin==d)].iloc[0]
            assert len(cell)==r.donor_n and r.fallback=='cell'  # Actual run used no fallback.
            values=cell.difference.to_numpy();donors[f,d]=values-values.mean()
            hold=assignments[(assignments.stop_count==0)&(assignments.split=='holdout')&(assignments.fitness_bin==f)&(assignments.depth_bin==d)]
            close([r.training_difference_mean,r.training_difference_variance,r.holdout_difference_mean,r.holdout_difference_variance],[values.mean(),values.var(),hold.difference.mean(),hold.difference.var(ddof=0)])
    empirical=pd.read_csv(OUT/'empirical_statistics.csv');nodes=pd.read_csv(OUT/'empirical_nodes.csv')
    lookup=table.set_index('genotype');source=lookup.loc[nodes.genotype]
    close(nodes.replicate_1,source.replicate_1,rtol=0,atol=1e-14);close(nodes.replicate_2,source.replicate_2,rtol=0,atol=1e-14)
    close(nodes.input_count,source.input_count)
    close(nodes.prior_imported_mean,(nodes.replicate_1+nodes.replicate_2)/2,rtol=0,atol=1e-7)
    prior_emp=pd.read_csv(INPUT/'prior/empirical_cubes.csv')
    expected_genotypes=[ ''.join(v) for s in prior_emp.node_keys for v in json.loads(s)]
    assert nodes.genotype.tolist()==expected_genotypes
    z=np.load(OUT/'empirical_components.npz');c=components((z['replicate_1_raw']-offset)/scale,(z['replicate_2_raw']-offset)/scale,z['shared_covariance_raw']/scale**2)
    close(z['shared_covariance_raw'],q[np.searchsorted(depth,z['input_depth'],side='right')]);verify_components(z,c)
    for res in RES:
        g=empirical[empirical.resolution==res]
        for name,val in statistic_columns(c,res).items():close(g[name],val)
    mutations=pd.read_csv(OUT/'empirical_mutations.csv')
    for res in RES:
        for j,m in enumerate(METHODS):
            g=mutations[(mutations.resolution==res)&(mutations.method==m)]
            close(g.numerator,c['numerator'][:,j].ravel());close(g.denominator,c['denominator'][:,j].ravel())
            close(g.shared_correction,c['correction'].ravel());close(g.mean_effect,c['effect'].ravel())
            close(g.slope,ratio(c['numerator'][:,j],c['denominator'][:,j],max(res*res,1e-18)).ravel())
    # Recompute all synthetic observation draws from donor pools and fixed uniforms.
    uniforms=np.random.default_rng(cfg['sampling']['observation_seed']).random((128,64,16,6))
    prior=pd.read_csv(INPUT/'prior/cube_values.csv');prior=prior[(prior.phase=='active_four_site')&(prior.split=='confirmation')]
    params=json.loads((OUT/'frozen_parameters.json').read_text())['within_ratio_pairs']
    reader=pd.read_csv(OUT/'synthetic_statistics.csv.gz',iterator=True)
    stats_parts=[];max_error=0.;n_cubes=0;n_rows=0
    for (readout,r),frame in prior.groupby(['readout','ratio'],sort=True):
        frame=frame.sort_values(['seed','cube']).reset_index(drop=True)
        par=next(x for x in params if x['readout']==readout and x['ratio']==r)
        for key in ['offset','scale','beta','p0']:
            if par[key] is not None:close(frame[key],par[key])
        raw=offset+scale*np.array([json.loads(x) for x in frame['values']])
        fbin=np.searchsorted(fitness,raw,side='right');fbin=np.broadcast_to(fbin,(128,64,16))
        dep=np.empty((128,64,16));e1=np.empty_like(dep);e2=np.empty_like(dep)
        def draw(pool,u):return pool[(u*len(pool)).astype(int)]
        for f in range(4):
            mask=fbin==f;dep[mask]=draw(depthpools[f],uniforms[...,0][mask])
        dbin=np.searchsorted(depth,dep,side='right');shared=np.empty_like(dep)
        for f in range(4):
            for d in range(4):
                mask=(fbin==f)&(dbin==d)
                e1[mask]=draw(donors[f,d],uniforms[...,1][mask])*np.where(uniforms[...,2][mask]<.5,-1,1)/np.sqrt(2)
                e2[mask]=draw(donors[f,d],uniforms[...,3][mask])*np.where(uniforms[...,4][mask]<.5,-1,1)/np.sqrt(2)
        for d in range(4):
            mask=dbin==d;shared[mask]=draw(common[d],uniforms[...,5][mask])
        for mode in MODES:
            stem=f'{readout}_r{r:g}_{mode}';z=np.load(OUT/'arrays'/f'{stem}.npz');index=pd.read_csv(OUT/'arrays'/f'{stem}_index.csv.gz')
            if mode=='intrinsic':a,b=raw,raw;expected_depth=np.full_like(raw,np.nan);expected_q=np.zeros_like(raw)
            else:
                common_error=shared if mode=='shared_stop_sensitivity' else 0
                a=(raw+e1+common_error).reshape(-1,16);b=(raw+e2+common_error).reshape(-1,16)
                expected_depth=dep.reshape(-1,16);expected_q=q[dbin].reshape(-1,16)
            close(z['replicate_1_raw'],a,rtol=1e-10,atol=2e-12);close(z['replicate_2_raw'],b,rtol=1e-10,atol=2e-12)
            close(z['input_depth'],expected_depth);close(z['shared_covariance_raw'],expected_q)
            expected_seed=np.tile(frame.seed.to_numpy(),1 if mode=='intrinsic' else 128)
            np.testing.assert_array_equal(index.seed,expected_seed)
            active_count=((a>cutoffs[0])&(b>cutoffs[1])).sum(1);np.testing.assert_array_equal(index.active_nodes,active_count)
            c=components((a-offset)/scale,(b-offset)/scale,expected_q/scale**2);verify_components(z,c)
            g=reader.get_chunk(len(index)*3)
            for res in RES:
                part=g[g.resolution==res]
                assert part.cube_id.tolist()==index.cube_id.tolist()
                for name,val in statistic_columns(c,res).items():max_error=max(max_error,close(part[name],val))
            stats_parts.append(g);n_cubes+=len(a);n_rows+=len(g)
    try:
        extra=reader.get_chunk(1);assert len(extra)==0
    except StopIteration:pass
    synthetic=pd.concat(stats_parts,ignore_index=True)
    audit=dict(remote_hashes=len(json.loads((BUNDLE/'remote_sha256.json').read_text())),output_hashes=len(manifest['output_sha256']),submitted_source_hashes=len(submitted['sha256']),input_hashes=14,empirical_measurements=len(nodes),activity_flags=len(table),excluded_genotypes=len(excluded),synthetic_pair_cubes=n_cubes,synthetic_statistic_rows=n_rows,empirical_statistic_rows=len(empirical),empirical_mutation_rows=len(mutations),replayed_observation_values=n_cubes*16*2,maximum_statistic_absolute_error=max_error,comparison_rtol=2e-7,comparison_atol=2e-10,component_rtol=1e-10,component_atol=2e-12,versions={p:importlib.metadata.version(p) for p in ['numpy','pandas','matplotlib','nbformat','nbclient','ipykernel']})
    dump(TABLES/'audit.json',audit)
    print(f'Audit complete: {n_cubes:,} replicate-pair cubes, {n_rows:,} synthetic statistic rows; all checks agree.')
    return dict(empirical=empirical,synthetic=synthetic,nodes=nodes,mutations=mutations,controls=controls,calibration=cal,assignments=assignments,audit=audit,manifest=manifest,config=cfg)


def summary_values(g,m,res):
    total=g[f'{m}_total_power'].sum();high=g[[f'{m}_power_{k}' for k in (3,4)]].to_numpy().sum()
    return dict(n_cubes=len(g),mean_total_power=total/len(g) if len(g) else np.nan,
                aggregate_ho_fraction=high/total if len(g) and total>len(g)*max(res*res,1e-18) else np.nan,
                mean_cube_ho_fraction=g[f'{m}_ho_fraction'].mean(),
                power_available_n=g[f'{m}_power_available'].sum(),mean_slope=g[f'{m}_slope'].mean(),
                slope_mutations_n=g[f'{m}_slope_n'].sum(),mean_beneficial_slope=g[f'{m}_beneficial_slope'].mean(),
                beneficial_mutations_n=g[f'{m}_beneficial_slope_n'].sum(),mean_sd=g.sd_fitness.mean(),
                mean_neutral_edge_fraction=g.neutral_edge_fraction.mean())


def summarize(c):
    empirical=c['empirical'];rows=[];bootstrap=[]
    rng=np.random.default_rng(c['config']['sampling']['bootstrap_seed'])
    old=pd.read_csv(OUT/'empirical_summary.csv');oldboot=pd.read_csv(OUT/'empirical_bootstrap.csv')
    for res in RES:
        for split in ['calibration','confirmation']:
            g=empirical[(empirical.resolution==res)&(empirical.split==split)]
            for st in STRATA:
                sub=g if st=='all' else g[g.activity_stratum==st]
                for m in METHODS:
                    vals=summary_values(sub,m,res)
                    saved=old[(old.resolution==res)&(old.split==split)&(old.stratum==st)&(old.method==m)].iloc[0]
                    for key,value in vals.items():close(saved[key],value)
                    rows.append(dict(split=split,resolution=res,stratum=st,method=m,**vals))
        g=empirical[(empirical.resolution==res)&(empirical.split=='confirmation')]
        indices=rng.integers(len(g),size=(2000,len(g)))
        for st in STRATA:
            member=np.ones(len(g),bool) if st=='all' else (g.activity_stratum==st).to_numpy()
            counts=member[indices].sum(1)
            for m in METHODS:
                total=g[f'{m}_total_power'].to_numpy();high=g[[f'{m}_power_{k}' for k in (3,4)]].sum(1).to_numpy()
                ts=np.where(member[indices],total[indices],0).sum(1);hs=np.where(member[indices],high[indices],0).sum(1)
                frac=ratio(hs,ts,np.maximum(counts,1)*max(res*res,1e-18))
                values={'aggregate_ho_fraction':frac}
                for metric,col in [('mean_cube_ho_fraction',f'{m}_ho_fraction'),('mean_slope',f'{m}_slope'),('mean_beneficial_slope',f'{m}_beneficial_slope'),('mean_total_power',f'{m}_total_power')]+[(f'mean_power_{k}',f'{m}_power_{k}') for k in range(1,5)]:
                    values[metric]=mean_finite(np.where(member[indices],g[col].to_numpy()[indices],np.nan))
                for metric,v in values.items():
                    finite=v[np.isfinite(v)]
                    lo,hi=np.quantile(finite,[.025,.975]) if len(finite) else (np.nan,np.nan)
                    rec=dict(resolution=res,stratum=st,method=m,metric=metric,n_cubes=member.sum(),valid_draws=len(finite),draws=2000,lower=lo,upper=hi)
                    bootstrap.append(rec)
                    saved=oldboot[(oldboot.resolution==res)&(oldboot.stratum==st)&(oldboot.method==m)&(oldboot.metric==metric)]
                    if len(saved):
                        close(saved[['lower','upper']].to_numpy()[0],[lo,hi]);assert saved.valid_draws.iloc[0]==len(finite)
    es=pd.DataFrame(rows);eb=pd.DataFrame(bootstrap)
    es.to_csv(TABLES/'empirical_summary.csv',index=False);eb.to_csv(TABLES/'empirical_bootstrap.csv',index=False)
    s=pd.read_csv(OUT/'synthetic_summary.csv');checks=0
    # All per-cube rows were audited above. Additionally independently aggregate
    # the first and last draw of each noisy arm and every intrinsic arm.
    for (readout,r,mode,draw,res),g in c['synthetic'][c['synthetic'].observation_draw.isin([-1,0,127])].groupby(['readout','ratio','mode','observation_draw','resolution'],sort=False):
        saved=s[(s.readout==readout)&(s.ratio==r)&(s['mode']==mode)&(s.observation_draw==draw)&(s.resolution==res)]
        for st in STRATA:
            sub=g if st=='all' else g[g.activity_stratum==st]
            for m in METHODS:
                vals=summary_values(sub,m,res);one=saved[(saved.stratum==st)&(saved.method==m)].iloc[0]
                for key,value in vals.items():close(one[key],value)
                checks+=1
    metrics=['mean_cube_ho_fraction','aggregate_ho_fraction','mean_slope','mean_beneficial_slope','mean_sd','mean_total_power','n_cubes','power_available_n','slope_mutations_n','beneficial_mutations_n','mean_neutral_edge_fraction']
    model=[]
    for keys,g in s.groupby(['readout','ratio','mode','resolution','stratum','method'],sort=False):
        base=dict(zip(['readout','ratio','mode','resolution','stratum','method'],keys))
        for metric in metrics:
            v=g[metric].to_numpy();finite=v[np.isfinite(v)];lo,hi=np.quantile(finite,[.025,.975]) if len(finite) else (np.nan,np.nan)
            model.append({**base,'metric':metric,'mean':finite.mean() if len(finite) else np.nan,'mc_lower':lo,'mc_upper':hi,'available_draws':len(finite),'draws':len(g)})
    ms=pd.DataFrame(model);ms.to_csv(TABLES/'frozen_model_summary.csv',index=False)
    contrasts=[]
    for keys,g in s.groupby(['readout','ratio','resolution','stratum','method'],sort=False):
        base=dict(zip(['readout','ratio','resolution','stratum','method'],keys));intrinsic=g[g['mode']=='intrinsic'].iloc[0]
        for mode in MODES[1:]:
            noisy=g[g['mode']==mode]
            for metric in ['mean_cube_ho_fraction','aggregate_ho_fraction','mean_slope','mean_beneficial_slope','mean_sd']:
                values=noisy[metric].to_numpy()-intrinsic[metric];finite=values[np.isfinite(values)]
                lo,hi=np.quantile(finite,[.025,.975]) if len(finite) else (np.nan,np.nan)
                contrasts.append({**base,'mode':mode,'metric':metric,'paired_change':finite.mean() if len(finite) else np.nan,'mc_lower':lo,'mc_upper':hi,'available_draws':len(finite)})
    pd.DataFrame(contrasts).to_csv(TABLES/'observation_contrasts.csv',index=False)
    pd.read_csv(OUT/'seed_summary.csv').to_csv(TABLES/'seed_summaries.csv',index=False)
    c.update(esummary=es,bootstrap=eb,model=ms,synthetic_summary=s)
    conf=empirical[(empirical.resolution==.001)&(empirical.split=='confirmation')]
    donor=c['calibration'].copy();donor['holdout_to_training_variance']=donor.holdout_difference_variance/donor.training_difference_variance
    donor.to_csv(TABLES/'donor_holdout_diagnostics.csv',index=False)
    c['controls'].to_csv(TABLES/'control_diagnostics.csv',index=False)
    counts=c['assignments'].groupby(['stop_count','split']).size()  # Retain direct total counts below.
    findings=dict(confirmation_cubes=len(conf),confirmation_active_vertices=int(conf.active_nodes.sum()),confirmation_vertices=len(conf)*16,
                  cube_strata=conf.activity_stratum.value_counts().to_dict(),top_three_cube_total_power_share=conf.average_total_power.nlargest(3).sum()/conf.average_total_power.sum(),
                  stop_adjusted_negative_total_cubes=int((conf.stop_adjusted_total_power<0).sum()),
                  stop_adjusted_outside_unit_fraction_cubes=int(((conf.stop_adjusted_ho_fraction<0)|(conf.stop_adjusted_ho_fraction>1)).sum()),
                  stop_training=int(((c['assignments'].stop_count>0)&(c['assignments'].split=='training')).sum()),
                  stop_holdout=int(((c['assignments'].stop_count>0)&(c['assignments'].split=='holdout')).sum()),
                  donor_training=int(((c['assignments'].stop_count==0)&(c['assignments'].split=='training')).sum()),
                  donor_holdout=int(((c['assignments'].stop_count==0)&(c['assignments'].split=='holdout')).sum()),
                  donor_variance_ratio_range=[donor.holdout_to_training_variance.min(),donor.holdout_to_training_variance.max()],
                  donor_fallback_cells=int((donor.fallback!='cell').sum()),
                  primary_empirical=es[(es.split=='confirmation')&(es.resolution==.001)&(es.stratum=='all')].to_dict('records'),
                  r4_model=ms[(ms.resolution==.001)&(ms.stratum=='all')&(ms.method=='average')&(ms.ratio==4)].to_dict('records'))
    dump(TABLES/'metrics.json',findings)
    c['audit']['independently_reaggregated_synthetic_summary_rows']=checks
    c['audit']['empirical_bootstrap_rows_reproduced']=len(oldboot)
    c['audit']['analysis_source_sha256']=sha(__file__)
    dump(TABLES/'audit.json',c['audit']);c['findings']=findings
    return c


def finish(fig,name):
    fig.savefig(FIGURES/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight')
    return fig


def point_interval(ax,x,mean,lower,upper,color,marker='o',label=None):
    ax.vlines(x,lower,upper,color=color,lw=1.5)
    ax.scatter(x,mean,color=color,marker=marker,s=35,label=label,zorder=3)


def empirical_point(c,method,metric,stratum='all',res=.001):
    g=c['esummary'];g=g[(g.split=='confirmation')&(g.resolution==res)&(g.stratum==stratum)&(g.method==method)]
    point=float(g.iloc[0][metric])
    b=c['bootstrap'];b=b[(b.resolution==res)&(b.stratum==stratum)&(b.method==method)&(b.metric==metric)]
    return point,*b.iloc[0][['lower','upper']].to_numpy(float)


def figure_measurement(c):
    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    nodes=c['nodes'];nodes=nodes[nodes.split=='confirmation'];active=nodes.published_active
    for sel,color,label in [(~active,'#777777','Below activity criterion'),(active,'#009E73','Active in both replicates')]:
        ax[0,0].scatter(nodes.loc[sel,'replicate_1'],nodes.loc[sel,'replicate_2'],s=15,alpha=.65,c=color,label=label)
    lim=[min(nodes.replicate_1.min(),nodes.replicate_2.min())-.02,max(nodes.replicate_1.max(),nodes.replicate_2.max())+.02]
    ax[0,0].plot(lim,lim,color='.6',lw=1);ax[0,0].axvline(c['config']['measurement']['activity_cutoffs'][0],color='.7',ls=':');ax[0,0].axhline(c['config']['measurement']['activity_cutoffs'][1],color='.7',ls=':')
    ax[0,0].set(xlabel='Replicate 1 (published fitness)',ylabel='Replicate 2 (published fitness)',title=f'A  Same 384 confirmation vertices; {int(active.sum())} active')
    ax[0,0].legend(fontsize=8)
    for split,marker in [('training','o'),('holdout','s')]:
        g=c['controls'][c['controls'].split==split]
        ax[0,1].plot(g.depth_bin+1,g.covariance,marker+'-',label=split)
        ax[1,0].plot(g.depth_bin+1,g.difference_variance/2,marker+'-',label=split)
    ax[0,1].set(xlabel='Input-depth quartile',ylabel='Stop cross covariance (raw fitness²)',title='B  Shared covariance transfers to held-out controls')
    ax[1,0].set(xlabel='Input-depth quartile',ylabel='Stop difference variance / 2 (raw fitness²)',title='C  Replicate disagreement also depends on depth')
    for a in (ax[0,1],ax[1,0]):a.set_xticks([1,2,3,4]);a.ticklabel_format(axis='y',style='sci',scilimits=(0,0));a.legend(fontsize=8)
    cal=c['calibration'];colors=['#555555','#0072B2','#E69F00','#009E73']
    for f,col in enumerate(colors):
        g=cal[cal.fitness_bin==f];ax[1,1].scatter(g.training_difference_variance,g.holdout_difference_variance,s=35,c=col,label=f'Fitness stratum {f+1}')
    lo,hi=4e-5,2e-3;ax[1,1].plot([lo,hi],[lo,hi],color='.6',lw=1)
    ax[1,1].set(xscale='log',yscale='log',xlabel='Training difference variance',ylabel='Held-out difference variance',title='D  Non-stop donor noise: 16 fixed cells')
    ax[1,1].legend(fontsize=8)
    fig.suptitle('Replicates constrain differential error; shared error remains an assumption',fontsize=14)
    return finish(fig,'01_measurement_and_controls')


def figure_wht(c):
    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    e=c['empirical'];e=e[(e.split=='confirmation')&(e.resolution==.001)]
    for j,(m,label,color) in enumerate(zip(METHODS,MLABELS,MCOLORS)):
        for k in range(1,5):
            b=c['bootstrap'];b=b[(b.resolution==.001)&(b.stratum=='all')&(b.method==m)&(b.metric==f'mean_power_{k}')].iloc[0]
            point_interval(ax[0,0],k+(j-1)*.16,e[f'{m}_power_{k}'].mean(),b.lower,b.upper,color,label=label if k==1 else None)
        for a,metric in [(ax[0,1],'mean_cube_ho_fraction'),(ax[1,0],'aggregate_ho_fraction')]:
            pt,lo,hi=empirical_point(c,m,metric);point_interval(a,j,pt,lo,hi,color)
    ax[0,0].set(xticks=[1,2,3,4],xlabel='WHT order',ylabel='Mean absolute power (RUN-166 units²)',title='A  Order-specific signed power')
    ax[0,0].ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax[0,0].legend(fontsize=8)
    for a,title in [(ax[0,1],'B  Equal weight per resolved cube'),(ax[1,0],'C  Pool powers, then form the fraction')]:
        a.set_xticks(range(3),['Average','Cross','Stop-adjusted']);a.set(ylabel='High-order fraction (orders 3–4)',title=title);a.axhline(0,color='.7',lw=.8)
    ax[0,1].set_xticks(range(3),['Average\n24/24 resolved','Cross\n24/24 resolved','Stop-adjusted\n12/24 resolved'])
    ax[1,0].set_xlabel('All 24 signed power estimates contribute to each pooled ratio')
    ax[0,1].axhline(1,color='.75',ls=':');ax[0,1].set_yscale('symlog',linthresh=.3)
    for st,color in [('no_active','#777777'),('some_active','#0072B2')]:
        g=e[e.activity_stratum==st];ax[1,1].scatter(g.average_total_power,g.stop_adjusted_total_power,c=color,s=35,label=st.replace('_',' '))
    lo,hi=e.average_total_power.min()/1.5,e.average_total_power.max()*1.5
    ax[1,1].plot([lo,hi],[lo,hi],color='.6',lw=.8);ax[1,1].axhline(0,color='black',lw=1)
    ax[1,1].set(xscale='log',yscale='symlog',xlabel='Ordinary total power',ylabel='Stop-adjusted total power',title='D  All 12 no-active cubes have negative adjusted totals')
    ax[1,1].set_yscale('symlog',linthresh=.0005);ax[1,1].legend(fontsize=8)
    fig.suptitle('WHT: estimator assumptions and weighting change the empirical target',fontsize=14)
    return finish(fig,'02_wht_weighting_and_resolution')


def figure_slopes(c):
    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for j,(m,label,color) in enumerate(zip(METHODS,MLABELS,MCOLORS)):
        for a,metric in [(ax[0,0],'mean_slope'),(ax[0,1],'mean_beneficial_slope')]:
            pt,lo,hi=empirical_point(c,m,metric);point_interval(a,j,pt,lo,hi,color)
        for i,res in enumerate(RES):
            g=c['esummary'];g=g[(g.split=='confirmation')&(g.resolution==res)&(g.stratum=='all')&(g.method==m)].iloc[0]
            ax[1,0].plot(i+(j-1)*.07,g.slope_mutations_n,'o',color=color,label=label if i==0 else None)
            ax[1,0].plot(i+(j-1)*.07,g.beneficial_mutations_n,'x',color=color)
    for a,title in [(ax[0,0],'A  All fixed 0→1 mutations'),(ax[0,1],'B  Mutations with positive mean effect')]:
        a.set_xticks(range(3),['Average','Cross','Stop-adjusted']);a.set(ylabel='Mean within-cube slope',title=title);a.axhline(0,color='.7',lw=.8)
    ax[1,0].set_xticks(range(3),['10⁻¹²','0.001','0.01']);ax[1,0].set(xlabel='Resolution in RUN-166 units',ylabel='Resolved mutation count',ylim=(-2,101),title='C  Circles: all mutations; crosses: beneficial subset');ax[1,0].legend(fontsize=8)
    g=c['mutations'];g=g[(g.split=='confirmation')&(g.resolution==.001)&(g.method=='stop_adjusted')]
    for flag,col,label in [(False,'#999999','Unresolved denominator'),(True,'#009E73','Resolved denominator')]:
        h=g[g.eligible==flag];ax[1,1].scatter(h.denominator,h.numerator,s=20,c=col,label=label)
    ax[1,1].axvline(0,color='black',lw=1);ax[1,1].axhline(0,color='.7',lw=.8)
    ax[1,1].set_xscale('symlog',linthresh=.0005);ax[1,1].set_yscale('symlog',linthresh=.0005)
    ax[1,1].set_xticks([-.001,0,.001,.01,.1],['−0.001','0','0.001','0.01','0.1'])
    ax[1,1].set(xlabel='Stop-adjusted background covariance',ylabel='Stop-adjusted effect/background covariance',title='D  Retain the numerator and denominator');ax[1,1].legend(fontsize=8)
    fig.suptitle('Diminishing returns remains sensitive to error assumptions and availability',fontsize=14)
    return finish(fig,'03_slopes_and_availability')


def figure_models(c):
    fig,axes=plt.subplots(2,3,figsize=(15,9),layout='constrained')
    ms=c['model'];metrics=['mean_cube_ho_fraction','mean_beneficial_slope','mean_sd'];labels=['Mean per-cube WHT fraction ≥3','Mean beneficial slope','Mean within-cube SD']
    for row,readout in enumerate(['linear','sigmoid']):
        for col,(metric,label) in enumerate(zip(metrics,labels)):
            ax=axes[row,col];point,lo,hi=empirical_point(c,'average',metric) if metric!='mean_sd' else (float(c['esummary'].query("split=='confirmation' and resolution==.001 and stratum=='all' and method=='average'").mean_sd.iloc[0]),np.nan,np.nan)
            ax.axhline(point,color='black',lw=1.3,label='Empirical average')
            if np.isfinite(lo):ax.axhspan(lo,hi,color='.6',alpha=.15)
            for mode,name,color in zip(MODES,NLABELS,NCOLORS):
                g=ms[(ms.readout==readout)&(ms['mode']==mode)&(ms.resolution==.001)&(ms.stratum=='all')&(ms.method=='average')&(ms.metric==metric)].sort_values('ratio')
                x=np.arange(4);ax.plot(x,g['mean'],'o-',color=color,label=name,ms=4);ax.fill_between(x,g.mc_lower,g.mc_upper,color=color,alpha=.12)
            ax.set_xticks(range(4),['0','0.25','1','4']);ax.set(xlabel='Frozen interaction/additive ratio',ylabel=label,title=f'{chr(65+row*3+col)}  {readout.capitalize()} readout')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Observation noise moves several summaries together; the generator is frozen',fontsize=14)
    return finish(fig,'04_frozen_models_with_observation')


def figure_strata(c):
    fig,axes=plt.subplots(2,3,figsize=(15,9),layout='constrained')
    ms=c['model'];metrics=['mean_cube_ho_fraction','mean_beneficial_slope']
    for col,st in enumerate(STRATA[1:]):
        for row,metric in enumerate(metrics):
            ax=axes[row,col];eg=c['esummary'];eg=eg[(eg.resolution==.001)&(eg.split=='confirmation')&(eg.stratum==st)&(eg.method=='average')].iloc[0]
            if eg.n_cubes:
                point,lo,hi=empirical_point(c,'average',metric,st);ax.axhline(point,color='black',label=f'Empirical ({int(eg.n_cubes)} cubes)');ax.axhspan(lo,hi,color='.6',alpha=.15)
            else:ax.text(.5,.95,'No empirical cubes in this stratum',ha='center',va='top',transform=ax.transAxes,fontsize=9)
            for j,(mode,name,color) in enumerate(zip(MODES,NLABELS,NCOLORS)):
                for i,readout in enumerate(['linear','sigmoid']):
                    g=ms[(ms.readout==readout)&(ms.ratio==4)&(ms['mode']==mode)&(ms.resolution==.001)&(ms.stratum==st)&(ms.method=='average')]
                    value=g[g.metric==metric].iloc[0];count=g[g.metric=='n_cubes'].iloc[0]['mean']
                    point_interval(ax,j+(i-.5)*.2,value['mean'],value.mc_lower,value.mc_upper,color,marker='o' if i==0 else 's',label=name if i==0 else None)
                    if row==0:ax.text(j+(i-.5)*.2,-.18,f'{count:.1f}',transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
            ax.set_xticks(range(3),['Intrinsic','Independent','Shared-stop']);ax.tick_params(axis='x',rotation=15)
            ax.set(ylabel='Mean WHT fraction ≥3' if row==0 else 'Mean beneficial slope',title=f'{chr(65+row*3+col)}  {SLABELS[col+1]}')
            if row==0:ax.set_xlabel('Mean synthetic cube count printed below ticks',labelpad=23,fontsize=8)
            ax.axhline(0,color='.7',lw=.7)
    handles=[Line2D([],[],marker='o',color='black',linestyle='',label='Linear'),Line2D([],[],marker='s',color='black',linestyle='',label='Sigmoid')]
    axes[1,0].legend(handles=handles,fontsize=9)
    fig.suptitle('At r=4, activity-conditioned populations differ between model and assay',fontsize=14)
    return finish(fig,'05_activity_strata')
