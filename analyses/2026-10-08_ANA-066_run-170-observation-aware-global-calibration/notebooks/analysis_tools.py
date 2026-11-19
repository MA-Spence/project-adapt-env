"""ANA-066: independent numerical audit and figure-led descriptive analysis.

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
BUNDLE=ROOT/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-170'
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
RES=[1e-12,.00066408129,.0066408129]
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


ARMS=['r4_linear','r4_sigmoid','r0_sigmoid']
ALABELS=['r=4 linear','r=4 sigmoid','r=0 sigmoid']
COLORS=['#777777','#0072B2','#D55E00']
OBS=['independent_difference','shared_stop_sensitivity']
OLABELS=['Independent error','Shared-stop sensitivity']


def raw_map(setting,a,g):
    from scipy.special import expit
    p=setting['parameters'];latent=a+setting['interaction_ratio']*g
    return p[0]+np.exp(p[1])*(latent if setting['arm'].endswith('linear') else expit(p[2]+np.exp(p[3])*latent))


def rebuild_noise(assignments,cfg):
    train=assignments[assignments.split=='training'];st=train[train.stop_count>0];dn=train[train.stop_count==0]
    depth=np.quantile(st.input_count,[.25,.5,.75]);cut=np.mean(cfg['measurement']['activity_cutoffs'])
    fitness=np.array([0,cut,dn.loc[dn.published_mean>cut,'published_mean'].median()])
    donors={};pools={};common=[];q=[];minimum=cfg['observation']['minimum_donors_per_cell']
    for f in range(4):
        fs=dn[dn.fitness_bin==f];pools[f]=(fs if len(fs)>=minimum else dn).input_count.to_numpy()
        for d in range(4):
            cell=fs[fs.depth_bin==d];selected=cell if len(cell)>=minimum else fs if len(fs)>=minimum else dn
            difference=(selected.replicate_1-selected.replicate_2).to_numpy();donors[f,d]=difference-difference.mean()
    for d in range(4):
        g=st[st.depth_bin==d];a,b=g.replicate_1.to_numpy(),g.replicate_2.to_numpy()
        qc=((a-a.mean())*(b-b.mean())).mean();q.append(qc)
        means=(a+b)/2;means-=means.mean();common.append(means*np.sqrt(max(0,qc)/means.var()))
    return depth,fitness,donors,pools,common,np.array(q)


def replay(raw,u,noise,shared):
    depth,fitness,donors,pools,common,q=noise;shape=u.shape[:-1]
    fbin=np.broadcast_to(np.searchsorted(fitness,raw,side='right'),shape)
    dval=np.empty(shape);e1=np.empty(shape);e2=np.empty(shape)
    def draw(pool,uu):return pool[(uu*len(pool)).astype(int)]
    for f in range(4):
        mask=fbin==f;dval[mask]=draw(pools[f],u[...,0][mask])
    dbin=np.searchsorted(depth,dval,side='right')
    for f in range(4):
        for d in range(4):
            mask=(fbin==f)&(dbin==d)
            e1[mask]=draw(donors[f,d],u[...,1][mask])*np.where(u[...,2][mask]<.5,-1,1)/np.sqrt(2)
            e2[mask]=draw(donors[f,d],u[...,3][mask])*np.where(u[...,4][mask]<.5,-1,1)/np.sqrt(2)
    sh=np.zeros(shape)
    if shared:
        for d in range(4):
            mask=dbin==d;sh[mask]=draw(common[d],u[...,5][mask])
    return raw+e1+sh,raw+e2+sh,dval,q[dbin]


def audit_and_load():
    FIGURES.mkdir(exist_ok=True);TABLES.mkdir(exist_ok=True)
    manifest=json.loads((OUT/'manifest.json').read_text());assert manifest['mode']=='PRODUCTION'
    assert (BUNDLE/'logs/labproj-exit-code.txt').read_text().strip()=='0'
    remote=json.loads((BUNDLE/'remote_sha256.json').read_text())
    for n,h in remote.items():assert sha(BUNDLE/n)==h,n
    for n,h in manifest['files'].items():assert sha(OUT/n)==h,n
    for n,h in manifest['inputs'].items():assert sha(INPUT/n)==h,n
    submitted=json.loads((BUNDLE/'submission-verification.json').read_text())
    for n,h in submitted['sha256'].items():assert sha(BUNDLE/'source'/n)==h,n
    cfg=yaml.safe_load((OUT/'config.yaml').read_text());settings=json.loads((OUT/'selected_parameters.json').read_text())
    nodes=pd.read_csv(OUT/'empirical_nodes.csv');original=pd.read_csv(INPUT/'replicates.csv').set_index('genotype')
    for col in ['replicate_1','replicate_2','input_count','published_mean']:close(nodes[col],original.loc[nodes.genotype,col])
    assert nodes.genotype.nunique()==2048
    old=set(json.loads((INPUT/'excluded_donor_genotypes.json').read_text()));assert not old.intersection(nodes.genotype)
    excluded=set(json.loads((OUT/'excluded_donor_genotypes.json').read_text()));assert excluded==old|set(nodes.genotype)
    assignments=pd.read_csv(OUT/'donor_assignments.csv.gz')
    assert not set(assignments.loc[(assignments.stop_count==0)&(assignments.split!='excluded'),'genotype']).intersection(excluded)
    ordered=original.reset_index().sort_values('genotype').reset_index(drop=True)
    assert ordered.genotype.tolist()==assignments.genotype.tolist()
    for col in ['replicate_1','replicate_2','input_count']:close(ordered[col],assignments[col])
    rng=np.random.default_rng(cfg['sampling']['donor_split_seed']);labels=np.full(len(ordered),'excluded',dtype=object)
    for mask in [ordered.stop_count>0,(ordered.stop_count==0)&~ordered.genotype.isin(excluded)]:
        ix=rng.permutation(ordered.index[mask]);n=int(.8*len(ix));labels[ix[:n]]='training';labels[ix[n:]]='holdout'
    np.testing.assert_array_equal(assignments.split,labels)
    noise=rebuild_noise(assignments,cfg);nm=json.loads((OUT/'noise_model.json').read_text())
    close(noise[0],nm['depth_breaks']);close(noise[1],nm['fitness_breaks']);close(noise[-1],nm['shared_covariance_raw'])
    empirical=pd.read_csv(OUT/'empirical_statistics.csv');synthetic=pd.read_csv(OUT/'synthetic_statistics.csv.gz')
    rows=0;arrays=0;replayed=0
    for indexpath in sorted(OUT.glob('*_index.csv')):
        stem=indexpath.name[:-10]
        if stem.startswith('latent_'):continue
        z=np.load(OUT/f'{stem}.npz');meta=pd.read_csv(indexpath);a,b,q=z['replicate_1'],z['replicate_2'],z['q'];c=components(a,b,q)
        for key,new in [('power','power'),('numerator','numerator'),('denominator','denominator'),('slope_correction','correction'),('mean_effect','effect')]:close(z[key],c[new],rtol=1e-10,atol=2e-12)
        if stem=='empirical':
            frame=empirical;close(a,nodes.replicate_1.to_numpy().reshape(-1,16));close(b,nodes.replicate_2.to_numpy().reshape(-1,16))
        else:
            m=meta.iloc[0];frame=synthetic[(synthetic.split==m.split)&(synthetic.observation==m.observation)&(synthetic.arm==m.arm)]
            mode=m.observation.replace('_intrinsic','');setting=next(s for s in settings if s['arm']==m.arm and s['observation']==mode)
            lat=np.load(OUT/f'latent_{m.split}.npz');raw=raw_map(setting,lat['a'],lat['g'])
            if m.observation.endswith('_intrinsic'):
                close(a,raw);close(b,raw)
            else:
                u=np.load(OUT/f'evaluation_uniforms_{m.split}.npz')['uniforms'];r1,r2,dep,qq=replay(raw,u,noise,mode==OBS[1])
                for v,vv in [(a,r1),(b,r2),(z['depth'],dep),(q,qq)]:close(v,vv.reshape(-1,16),rtol=1e-9,atol=2e-12)
                replayed+=a.size+b.size
        for res in RES:
            f=frame[np.isclose(frame.resolution,res,rtol=1e-9,atol=1e-18)]
            assert len(f)==len(a)
            for col,val in statistic_columns(c,res).items():close(f[col],val)
            rows+=len(f)
        cut=cfg['measurement']['activity_cutoffs'];ac=((a>cut[0])&(b>cut[1])).sum(1)
        np.testing.assert_array_equal(meta.active_vertices,ac)
        arrays+=len(a)
    target=json.loads((OUT/'calibration_target.json').read_text())
    cal=nodes[nodes.split=='calibration'];values=((cal.replicate_1+cal.replicate_2)/2).to_numpy().reshape(-1,16)
    close(target['quantiles'],np.quantile(values,target['probabilities']));close(target['variance'],values.var());close(target['mean_sd'],values.std(1).mean())
    lat=np.load(OUT/'latent_calibration.npz');u=np.load(OUT/'fit_uniforms.npz')['uniforms']
    for setting in settings:
        raw=raw_map(setting,lat['a'],lat['g']);a,b,_,_=replay(raw,u,noise,setting['observation']==OBS[1]);y=(a+b)/2
        terms=[np.mean((np.quantile(y,target['probabilities'])-target['quantiles'])**2)/target['variance'],(y.std(-1).mean()/target['mean_sd']-1)**2,((y<=target['cutoff']).mean()-target['floor_cdf'])**2]
        close(list(setting['objective'].values()),terms+[sum(terms)],rtol=2e-6)
        trace=pd.read_csv(OUT/f"fit_{setting['observation']}_{setting['arm']}.csv.gz")
        close(trace.total.min(),sum(terms));assert len(trace)==setting['evaluations']
    audit=dict(output_hashes=len(manifest['files']),remote_hashes=len(remote),source_hashes=len(submitted['sha256']),
        empirical_genotypes=len(nodes),excluded_genotypes=len(excluded),independent_statistic_rows=rows,
        replicate_pair_cubes=arrays,replayed_observations=replayed,fits_recomputed=len(settings),checks='Numerical/provenance identities, not scientific gates.')
    dump(TABLES/'audit.json',audit)
    pd.DataFrame([{k:v for k,v in s.items() if k not in ('objective','thermodynamic_config','bounds','parameters','relative_boundary_distance')}|s['objective'] for s in settings]).to_csv(TABLES/'selected_parameters.csv',index=False)
    dump(TABLES/'portable_parameters.json',settings)
    return dict(cfg=cfg,settings=settings,empirical=empirical,synthetic=synthetic,nodes=nodes,audit=audit)


def estimate(frame,method='average'):
    p=np.array([frame[f'{method}_power_{k}'].sum() for k in range(1,5)])
    return dict(mean_sd=frame.sd_fitness.mean(),mean_fitness=frame.mean_fitness.mean(),
        equal_ho=frame[f'{method}_ho_fraction'].mean(),pooled_ho=p[2:].sum()/p.sum(),
        equal_epistasis=frame[f'{method}_epistatic_fraction'].mean(),pooled_epistasis=p[1:].sum()/p.sum(),
        slope=frame[f'{method}_slope'].mean(),beneficial_slope=frame[f'{method}_beneficial_slope'].mean(),
        slope_availability=frame[f'{method}_slope_n'].sum()/(4*len(frame)),
        beneficial_count=frame[f'{method}_beneficial_slope_n'].sum()/len(frame),
        no_active=(frame.active_vertices==0).mean(),active_vertices=frame.active_vertices.mean()/16,
        **{f'power_{k}':p[k-1]/len(frame) for k in range(1,5)})


def bootstrap(frame,unit,method,rng,draws=1000):
    # Cluster sums/counts preserve equal-cube estimands and unavailable ratios.
    means={'mean_sd':'sd_fitness','mean_fitness':'mean_fitness','equal_ho':f'{method}_ho_fraction',
        'equal_epistasis':f'{method}_epistatic_fraction','slope':f'{method}_slope','beneficial_slope':f'{method}_beneficial_slope'}
    group=frame.groupby(unit,sort=True);n=group.size().to_numpy();ix=rng.integers(len(n),size=(draws,len(n)));total=n[ix].sum(1)
    out={}
    for metric,col in means.items():
        sums=group[col].sum().to_numpy();counts=group[col].count().to_numpy();den=counts[ix].sum(1)
        out[metric]=np.divide(sums[ix].sum(1),den,out=np.full(draws,np.nan),where=den>0)
    powers=[]
    for k in range(1,5):
        s=group[f'{method}_power_{k}'].sum().to_numpy()[ix].sum(1);powers.append(s);out[f'power_{k}']=s/total
    powers=np.array(powers);den=powers.sum(0);threshold=RES[1]**2*total
    out['pooled_ho']=ratio(powers[2:].sum(0),den,threshold);out['pooled_epistasis']=ratio(powers[1:].sum(0),den,threshold)
    for metric,vals in [('slope_availability',frame[f'{method}_slope_n']/4),('beneficial_count',frame[f'{method}_beneficial_slope_n']),('no_active',(frame.active_vertices==0).astype(float)),('active_vertices',frame.active_vertices/16)]:
        s=frame.assign(value=vals).groupby(unit).value.sum().to_numpy();out[metric]=s[ix].sum(1)/total
    return {metric:(float(np.nanquantile(v,.025)),float(np.nanquantile(v,.975)),int(np.isfinite(v).sum())) for metric,v in out.items()}


def summarize(data):
    rng=np.random.default_rng(652601);res=RES[1]
    e=data['empirical'];e=e[(e.split=='confirmation')&np.isclose(e.resolution,res,atol=1e-18)]
    s=data['synthetic'];s=s[(s.split=='confirmation')&np.isclose(s.resolution,res,atol=1e-18)]
    rows=[]
    groups=[('empirical','empirical',e,'cube')]+[(mode,arm,g,'seed') for (mode,arm),g in s.groupby(['observation','arm'],sort=False)]
    for mode,arm,g,unit in groups:
        for method in METHODS:
            est=estimate(g,method);ci=bootstrap(g,unit,method,rng)
            for metric,value in est.items():
                lo,hi,n=ci[metric];rows.append(dict(observation=mode,arm=arm,method=method,metric=metric,value=value,lower=lo,upper=hi,bootstrap_available=n,clusters=g[unit].nunique()))
    summary=pd.DataFrame(rows);summary.to_csv(TABLES/'confirmation_estimates.csv',index=False)
    data['summary']=summary;data['e']=e;data['s']=s
    return summary


def save(fig,name):
    fig.savefig(FIGURES/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(FIGURES/f'{name}.pdf',bbox_inches='tight')
    return fig


def figure_fit(data):
    fig,axes=plt.subplots(1,2,figsize=(11,4.4))
    for j,mode in enumerate(OBS):
        ax=axes[j]
        for arm,label,color in zip(ARMS,ALABELS,COLORS):
            tr=pd.read_csv(OUT/f'fit_{mode}_{arm}.csv.gz')
            ax.plot(tr.evaluation,np.minimum.accumulate(tr.total),label=label,color=color)
        ax.set(xlabel='Objective evaluations',ylabel='Best calibration loss',yscale='log',title=OLABELS[j])
    axes[0].legend(frameon=False);fig.suptitle('1 · Fit only measured marginals and average local spread',fontsize=14);fig.tight_layout()
    return save(fig,'01_calibration')


def figure_marginals(data):
    fig,axes=plt.subplots(2,2,figsize=(11,8));e=data['e'];nodes=data['nodes'];nodes=nodes[nodes.split=='confirmation']
    y=(nodes.replicate_1+nodes.replicate_2).to_numpy()/2
    def cdf(ax,v,**kw):v=np.sort(np.ravel(v));ax.plot(v,np.arange(1,len(v)+1)/len(v),**kw)
    for j,mode in enumerate(OBS):
        cdf(axes[0,j],y,color='black',lw=2,label='Empirical confirmation')
        cdf(axes[1,j],e.sd_fitness,color='black',lw=2,label='Empirical confirmation')
        for arm,label,color in zip(ARMS,ALABELS,COLORS):
            z=np.load(OUT/f'confirmation_{mode}_{arm}.npz');v=(z['replicate_1']+z['replicate_2'])/2
            cdf(axes[0,j],v,color=color,label=label)
            cdf(axes[1,j],v.std(1),color=color,label=label)
        axes[0,j].set(title=OLABELS[j],xlabel='Measured fitness',ylabel='Cumulative fraction',xlim=(-.08,.7))
        axes[1,j].set(xlabel='Within-cube SD (raw fitness)',ylabel='Cumulative fraction',xlim=(0,.3))
        axes[0,j].axvline(np.mean(data['cfg']['measurement']['activity_cutoffs']),color='grey',ls=':',lw=1)
    axes[0,0].legend(frameon=False,fontsize=9);fig.suptitle('2 · Fresh confirmation: fitness shape and local spread',fontsize=14);fig.tight_layout()
    return save(fig,'02_heldout_marginals')


def points(ax,data,metrics,mode,method='average',ylabel=None):
    table=data['summary'];sources=[('empirical','empirical','Empirical','black')]+[(mode,a,l,c) for a,l,c in zip(ARMS,ALABELS,COLORS)]
    for j,(obs,arm,label,color) in enumerate(sources):
        g=table[(table.observation==obs)&(table.arm==arm)&(table.method==method)].set_index('metric')
        for k,metric in enumerate(metrics):
            r=g.loc[metric];x=k+(j-1.5)*.16
            ax.plot([x,x],[r.lower,r.upper],color=color,lw=1.5)
            ax.scatter([x],[r.value],color=color,s=30,label=label if k==0 else None,zorder=3)
    ax.set_xticks(range(len(metrics)));ax.set_ylabel(ylabel)


def figure_wht(data):
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for j,mode in enumerate(OBS):
        points(axes[0,j],data,[f'power_{k}' for k in range(1,5)],mode,ylabel='Mean absolute WHT power')
        axes[0,j].set(title=OLABELS[j],yscale='log');axes[0,j].set_xticklabels(['Order 1','Order 2','Order 3','Order 4'])
        points(axes[1,j],data,['equal_epistasis','pooled_epistasis','equal_ho','pooled_ho'],mode,ylabel='Fraction of nonconstant power')
        axes[1,j].set_xticklabels(['Epistasis\nequal cube','Epistasis\npooled','Orders 3–4\nequal cube','Orders 3–4\npooled'],fontsize=9)
        axes[1,j].set_ylim(0,1)
    axes[0,0].legend(frameon=False,fontsize=9);fig.suptitle('3 · Withheld WHT: absolute power and two distinct weightings',fontsize=14);fig.tight_layout()
    return save(fig,'03_withheld_wht')


def figure_slopes(data):
    fig,axes=plt.subplots(3,2,figsize=(11,11))
    for j,mode in enumerate(OBS):
        points(axes[0,j],data,['slope','beneficial_slope'],mode,ylabel='Mean fixed-mutation slope')
        axes[0,j].set(title=OLABELS[j]);axes[0,j].set_xticklabels(['All mutations','Beneficial only']);axes[0,j].axhline(0,color='grey',lw=.8)
        points(axes[1,j],data,['slope_availability'],mode,ylabel='Fraction of mutations with resolved denominator')
        axes[1,j].set_xticklabels(['All-mutation availability']);axes[1,j].set_ylim(0,1.05)
        table=data['summary']
        for i,(arm,label,color) in enumerate(zip(ARMS,ALABELS,COLORS)):
            g=table[(table.observation==mode+'_intrinsic')&(table.arm==arm)&(table.method=='average')].set_index('metric')
            for k,metric in enumerate(['slope','beneficial_slope']):
                r=g.loc[metric];x=k+(i-1)*.18
                axes[2,j].plot([x,x],[r.lower,r.upper],color=color)
                axes[2,j].scatter(x,r.value,color=color,s=30)
        axes[2,j].set_xticks([0,1],['All mutations','Beneficial only'])
        axes[2,j].set_ylabel('Intrinsic map slope (before noise)')
        axes[2,j].axhline(0,color='grey',lw=.8)
    axes[0,0].legend(frameon=False,fontsize=9);fig.suptitle('4 · Withheld diminishing returns and estimator availability',fontsize=14);fig.tight_layout()
    return save(fig,'04_withheld_slopes')


def figure_activity(data):
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for j,mode in enumerate(OBS):
        points(axes[0,j],data,['no_active','active_vertices'],mode,ylabel='Fraction')
        axes[0,j].set(title=OLABELS[j],ylim=(0,1));axes[0,j].set_xticklabels(['Cubes with no active vertices','Active vertices'],fontsize=9)
        e=data['e'];axes[1,j].scatter(e.mean_fitness,e.average_ho_fraction,color='black',s=23,alpha=.65,label='Empirical')
        g=data['s'];g=g[(g.observation==mode)&(g.arm=='r4_sigmoid')]
        # Average repeated observations per original cube, retaining 16 seed clusters.
        cube=g.groupby(['seed','cube'])[['mean_fitness','average_ho_fraction']].mean()
        axes[1,j].scatter(cube.mean_fitness,cube.average_ho_fraction,color=COLORS[1],s=16,alpha=.45,label='r=4 sigmoid')
        axes[1,j].set(xlabel='Cube mean measured fitness',ylabel='Orders 3–4 fraction',ylim=(0,1))
    axes[0,0].legend(frameon=False,fontsize=9);axes[1,0].legend(frameon=False,fontsize=9)
    fig.suptitle('5 · Activity composition and where high-order fractions occur',fontsize=14);fig.tight_layout()
    return save(fig,'05_activity_and_structure')
