"""Reproducible ANA-067 calculations and figures; no fitting or outcome-based gates."""
from pathlib import Path
import hashlib, json, os, sys, tempfile, importlib.metadata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.special import expit
import yaml

ANA=Path(__file__).resolve().parents[1]
ROOT=ANA.parents[1]
DATA=ROOT/'data/processed/production_landscapes/RUN-172'
SOURCE=DATA/'provenance/source'
# Import frozen execution sources without adding bytecode to the DVC bundle.
sys.dont_write_bytecode=True
EXP=SOURCE/'experiments/2026-10-08_EXP-066_production-landscape-aggregate-realism'
TABLES=ANA/'tables';FIGURES=ANA/'figures'
for p in (TABLES,FIGURES):p.mkdir(exist_ok=True)
sys.path.insert(0,str(SOURCE/'scripts'))
sys.path.insert(0,str(SOURCE/'external/Adapt-Env'))
from adaptenv import ThermodynamicConfig,ThermodynamicLandscape
from observation_calibration_exp065 import observation_bank, observed
from replicate_observation_exp065 import calibrate_noise

COLORS={'empirical':'#242830','r4_sigmoid':'#0072B2','r0_sigmoid':'#D55E00','r4_linear':'#009E73'}
LABELS={'empirical':'Measured','r4_sigmoid':'r=4 sigmoid','r0_sigmoid':'r=0 sigmoid','r4_linear':'r=4 linear'}
FUNCTION_NAMES={'gb1_wu2016_four_site_graphfla':'GB1', 'trpb_johnston2024_four_site_graphfla':'TrpB',
 'cr9114_h1_phillips2021_graphfla':'CR9114 H1','cr9114_h3_phillips2021_graphfla':'CR9114 H3',
 'cr9114_flub_phillips2021_graphfla':'CR9114 FluB','cr6261_h1_phillips2021_graphfla':'CR6261 H1',
 'cr6261_h9_phillips2021_graphfla':'CR6261 H9','pard2_lite2020_graphfla':'ParD2','pard3_lite2020_graphfla':'ParD3',
 'proteingym_f7ybw8_mesow_aakre_2015':'Aakre ParD3','proteingym_gfp_aeqvi_sarkisyan_2016':'GFP',
 'proteingym_his7_yeast_pokusaeva_2019':'HIS7','proteingym_phot_chlre_chen_2023':'PHOT',
 'proteingym_spg1_strsg_olson_2014':'Olson GB1'}
ARMS=['r4_sigmoid','r0_sigmoid','r4_linear']


def name(dataset):
    if dataset in FUNCTION_NAMES:return FUNCTION_NAMES[dataset]
    p=dataset.split('_');return p[1].upper()+' / '+p[-1].upper()


def family(dataset):return 'Stability' if 'tsuboyama' in dataset else 'Function / binding'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, obj):
    def default(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,np.ndarray):return x.tolist()
        raise TypeError(type(x))
    Path(path).write_text(json.dumps(obj,indent=2,default=default,allow_nan=False)+'\n')


def seed(label):return int.from_bytes(hashlib.sha256(label.encode()).digest()[:4],'little')


def read_fitness(dataset,split,observation,arm):
    path=DATA/dataset
    if arm=='empirical':
        y=np.load(path/'empirical_panel.npz')['fitness'];meta=pd.read_csv(path/'empirical_panel_index.csv')
        return y[meta.split.eq(split)]
    return np.load(path/f'{split}_{observation}_{arm}.npz')['fitness']


def independent_statistics(y,res):
    """Butterfly FWHT; direct within-direction covariance and sign-face counts."""
    y=np.asarray(y);width=y.shape[-1];dim=int(np.log2(width))
    # Subtract a represented vertex before transforming to avoid cancellation
    # of a large common offset; nonconstant Walsh coefficients are invariant.
    coef=y-y[:,:1];block=1
    while block<width:
        for j in range(0,width,2*block):
            a,b=coef[:,j:j+block].copy(),coef[:,j+block:j+2*block].copy()
            coef[:,j:j+block],coef[:,j+block:j+2*block]=a+b,a-b
        block*=2
    coef/=width;order=np.array([i.bit_count() for i in range(width)])
    powers=np.array([(coef[:,order==k]**2).sum(axis=1) for k in range(1,dim+1)]).T
    total=powers.sum(axis=1);threshold=max(res**2,1e-24)
    low=np.array([[x for x in range(width) if x&(1<<j)==0] for j in range(dim)])
    high=low^(1<<np.arange(dim))[:,None];x=y[:,low];effects=y[:,high]-x
    centered=x-x.mean(axis=-1,keepdims=True);den=(centered**2).mean(axis=-1)
    num=(centered*(effects-effects.mean(axis=-1,keepdims=True))).mean(axis=-1)
    slopes=np.divide(num,den,out=np.full_like(num,np.nan),where=den>threshold)
    beneficial=effects.mean(axis=-1)>res
    def avg(a):
        finite=np.isfinite(a);return np.divide(np.where(finite,a,0).sum(axis=-1),finite.sum(axis=-1),out=np.full(len(y),np.nan),where=finite.sum(axis=-1)>0)
    epi=np.divide(powers[:,1:].sum(axis=1),total,out=np.full(len(y),np.nan),where=total>threshold)
    ho=np.divide(powers[:,2:].sum(axis=1),total,out=np.full(len(y),np.nan),where=total>threshold) if dim>2 else np.full(len(y),np.nan)
    signs=[];rsigns=[]
    for j in range(dim):
        for k in range(j+1,dim):
            base=np.array([i for i in range(width) if i&(1<<j)==0 and i&(1<<k)==0])
            ej0=y[:,base^(1<<j)]-y[:,base];ej1=y[:,base^(1<<j)^(1<<k)]-y[:,base^(1<<k)]
            ek0=y[:,base^(1<<k)]-y[:,base];ek1=y[:,base^(1<<j)^(1<<k)]-y[:,base^(1<<j)]
            jflip=((ej0>res)&(ej1 < -res))|((ej1>res)&(ej0 < -res))
            kflip=((ek0>res)&(ek1 < -res))|((ek1>res)&(ek0 < -res))
            signs.append(jflip|kflip);rsigns.append(jflip&kflip)
    answer=dict(total_power=total,mean_fitness=y.mean(axis=1),sd_fitness=y.std(axis=1),
       epistatic_fraction=epi,ho_fraction=ho,slope=avg(slopes),beneficial_slope=avg(np.where(beneficial,slopes,np.nan)),
       slope_n=(den>threshold).sum(axis=1),beneficial_slope_n=((den>threshold)&beneficial).sum(axis=1),
       neutral_edge_fraction=(abs(effects)<=res).mean(axis=(1,2)),
       sign_epistasis_fraction=np.concatenate(signs,axis=1).mean(axis=1),
       reciprocal_sign_epistasis_fraction=np.concatenate(rsigns,axis=1).mean(axis=1))
    answer.update({f'power_{k}':powers[:,k-1] if k<=dim else np.full(len(y),np.nan) for k in range(1,5)})
    return answer


def audit():
    manifest=json.loads((DATA/'manifest.json').read_text())
    assert manifest['evaluated']==63 and manifest['fitted_settings']==195 and not manifest['execution_errors']
    for relative,expected in manifest['files'].items():assert sha(DATA/relative)==expected,relative
    source_manifest=json.loads((EXP/'source_manifest.json').read_text())
    for relative,expected in source_manifest['files'].items():assert sha(SOURCE/relative)==expected,relative
    assert source_manifest['adaptenv_commit']==manifest['adaptenv_commit']
    inventory=pd.read_csv(DATA/'dataset_inventory.csv');suitability=pd.read_csv(DATA/'suitability.csv')
    missing=sorted(set(inventory.dataset_id)-set(suitability.dataset_id))
    assert missing==['flip_aav_full_data','flip_gb1_four_mutations_full_data','flip_meltome_stability_summary']
    reasons={'flip_aav_full_data':'mixed insertion/deletion library; fixed-length substitution model out of scope',
     'flip_gb1_four_mutations_full_data':'duplicate GB1 mirror; primary representation included',
     'flip_meltome_stability_summary':'cross-protein melting temperatures, not a mutational landscape'}
    for dataset in missing:
        suitability=pd.concat([suitability,pd.DataFrame([dict(dataset_id=dataset,status='out_of_scope_or_duplicate',reason=reasons[dataset],analysis_restored=True)])],ignore_index=True)
    suitability['family']=suitability.dataset_id.map(family)
    suitability.to_csv(TABLES/'suitability_complete.csv',index=False)
    fits=json.loads((DATA/'selected_parameters.json').read_text());cfg=yaml.safe_load((DATA/'config.yaml').read_text())
    frames=[];audited_rows=0;replayed_nodes=0;public_nodes=0;panel_nodes=0;loss_rows=[];roundoff=[]
    ids=sorted({s['dataset_id'] for s in fits})
    for di,dataset in enumerate(ids):
        path=DATA/dataset;geometry=json.loads((path/'geometry.json').read_text())
        expected=inventory.set_index('dataset_id').loc[dataset,'sha256_parquet'];assert geometry['input_sha256']==expected
        panel=np.load(path/'empirical_panel.npz');nodes=panel['genotypes'];values=panel['fitness'];meta=pd.read_csv(path/'empirical_panel_index.csv')
        cal=meta.split.eq('calibration').to_numpy();conf=~cal;n=cal.sum();panel_nodes+=nodes.shape[0]*nodes.shape[1]
        keys=[[v.tobytes() for v in cube] for cube in nodes];flat=[k for cube in keys for k in cube]
        assert len(set(flat))==len(flat),dataset
        dim=int(np.log2(values.shape[1]));bits=((np.arange(2**dim)[:,None]>>np.arange(dim))&1).astype(bool)
        for cube in nodes:
            sites=np.flatnonzero(np.any(cube!=cube[0],axis=0));assert len(sites)==dim
            np.testing.assert_array_equal(cube[:,sites],np.where(bits,cube[-1,sites],cube[0,sites]))
        stats=pd.read_csv(path/'statistics.csv.gz');stats['relative_resolution']=stats.resolution/next(s['calibration_scale'] for s in fits if s['dataset_id']==dataset)
        frames.append(stats)
        latents={};indices={}
        for split,mask in [('calibration',cal),('confirmation',conf)]:
            z=np.load(path/f'latent_{split}.npz');index=pd.read_csv(path/f'latent_{split}_index.csv');latents[split]=z;indices[split]=index
            np.testing.assert_array_equal(z['genotypes'],nodes[mask][index.geometry_cube.to_numpy()])
        assert set(indices['calibration'].seed).isdisjoint(indices['confirmation'].seed)
        assert indices['calibration'].seed.nunique()==8 and indices['confirmation'].seed.nunique()==16
        noise=None
        if dataset=='trpb_johnston2024_four_site_graphfla':
            donors=pd.read_csv(path/'donor_assignments.csv.gz');excluded=set(json.loads((path/'excluded_noise_donors.json').read_text()))
            assert set(donors.loc[donors.split.eq('training') & donors.stop_count.eq(0),'genotype']).isdisjoint(excluded)
            with tempfile.TemporaryDirectory(prefix='ana067-noise-') as td:noise,_=calibrate_noise(donors,excluded,cfg,Path(td))
        for setting in [s for s in fits if s['dataset_id']==dataset]:
            arm,mode=setting['arm'],setting['observation'];params=np.array(setting['parameters'])
            def raw_readout(split):
                z=latents[split];v=z['a']+setting['interaction_ratio']*z['g']
                return setting['offset']+setting['scale']*(v if arm.endswith('linear') else expit(setting['position']+setting['steepness']*v))
            for split in ['calibration','confirmation']:
                raw=raw_readout(split);label=mode+'_intrinsic' if mode!='none' else mode
                saved=read_fitness(dataset,split,label,arm)
                np.testing.assert_allclose(raw,saved,rtol=1e-10,atol=1e-10*setting['calibration_scale'])
                replayed_nodes+=raw.size
                if mode!='none':
                    reps=np.load(path/f'{split}_{mode}_{arm}_replicates.npz')
                    np.testing.assert_array_equal(((reps['replicate_1']+reps['replicate_2'])/2).reshape(-1,2**dim),read_fitness(dataset,split,mode,arm))
            if mode=='none' and arm=='r4_sigmoid':
                config=ThermodynamicConfig(**setting['thermodynamic_config']);model=ThermodynamicLandscape(config)
                # One full cube from each assay provides a source/API replay in addition to all cached-readout identities.
                test=latents['confirmation']['genotypes'][0]
                np.testing.assert_allclose(model.evaluate_batch(test),raw_readout('confirmation')[0],rtol=1e-9,atol=1e-9*setting['calibration_scale']);public_nodes+=len(test)
            target=(values[cal]-setting['calibration_origin'])/setting['calibration_scale'];raw=raw_readout('calibration')
            if mode!='none':
                uniforms=np.load(path/'fit_observation_uniforms.npz')['uniforms'];bank=observation_bank(noise,uniforms,mode=='shared_stop_sensitivity');r1,r2,_,_=observed(raw,bank);raw=(r1+r2)/2
            y=(raw-setting['calibration_origin'])/setting['calibration_scale'];probs=np.linspace(0,1,cfg['fit']['quantiles'])
            floor=(setting['near_floor_raw']-setting['calibration_origin'])/setting['calibration_scale']
            parts=dict(quantile=np.mean((np.quantile(y,probs)-np.quantile(target,probs))**2)/target.var(),mean_sd=(y.std(axis=-1).mean()/target.std(axis=-1).mean()-1)**2,floor_cdf=((y<=floor).mean()-(target<=floor).mean())**2)
            loss=sum(parts.values());np.testing.assert_allclose(loss,setting['objective']['total'],rtol=1e-8,atol=1e-9)
            trace=pd.read_csv(path/f'fit_{mode}_{arm}.csv.gz');assert len(trace)==setting['evaluations'];np.testing.assert_allclose(trace.total.min(),loss,rtol=1e-8,atol=1e-9)
            loss_rows.append(dict(dataset_id=dataset,arm=arm,observation=mode,recomputed_loss=loss,reported_loss=setting['objective']['total']))
        for (split,obs,arm),group in stats.groupby(['split','observation','arm'],sort=False):
            y=read_fitness(dataset,split,obs,arm)
            for res,rows in group.groupby('resolution',sort=False):
                recomputed=independent_statistics(y,res)
                for key,expected in recomputed.items():
                    actual=rows[key].to_numpy()
                    if key in ['epistatic_fraction','ho_fraction']:
                        np.testing.assert_array_equal(np.isnan(actual),np.isnan(expected))
                        # Relative power fractions become ill-conditioned on
                        # saturated cubes. Bound propagation of floating-point
                        # coefficient roundoff, and preserve every discrepancy.
                        bound=np.maximum(1e-9+2e-7*np.abs(expected),32*y.shape[1]*np.finfo(float).eps*np.max(abs(y),axis=1)/np.maximum(y.std(axis=1),1e-300))
                        finite=np.isfinite(expected);delta=abs(actual-expected)
                        assert np.all(delta[finite]<=bound[finite]),(dataset,split,obs,arm,key)
                        for i in np.flatnonzero(finite & (delta>1e-9+2e-7*abs(expected))):
                            roundoff.append(dict(dataset_id=dataset,split=split,observation=obs,arm=arm,metric=key,resolution=res,cube_row=int(i),absolute_difference=delta[i],roundoff_bound=bound[i],cube_sd=y[i].std()))
                    else:
                        np.testing.assert_allclose(actual,expected,rtol=2e-7,atol=1e-9,equal_nan=True,err_msg=f'{dataset}/{split}/{obs}/{arm}/{key}')
                audited_rows+=len(rows)
        if (di+1)%10==0:print('Audited',di+1,'of',len(ids),'assays',flush=True)
    stats=pd.concat(frames,ignore_index=True);stats['family']=stats.dataset_id.map(family)
    # Production summary is audited against per-assay records, not used as an independent input.
    confirmation=pd.read_csv(DATA/'confirmation_statistics.csv');assert len(confirmation)==int(stats.split.eq('confirmation').sum())
    pd.DataFrame(loss_rows).to_csv(TABLES/'fit_loss_audit.csv',index=False)
    pd.DataFrame(roundoff).to_csv(TABLES/'low_power_roundoff_audit.csv',index=False)
    audit_info=dict(output_hashes=len(manifest['files']),source_hashes=len(source_manifest['files']),restored_inventory_rows=missing,
        evaluated_assays=len(ids),fitted_settings=len(fits),statistic_rows=audited_rows,empirical_nodes=panel_nodes,
        replayed_intrinsic_nodes=replayed_nodes,public_API_replayed_nodes=public_nodes,losses_recomputed=len(loss_rows),roundoff_sensitive_fraction_rows=len(roundoff),
        inference='Numerical/provenance checks, not scientific acceptance criteria.',
        runtime_versions={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','matplotlib','nbformat','nbclient']})
    dump(TABLES/'audit.json',audit_info)
    return stats,fits,suitability,audit_info


METRICS=['mean_fitness','sd_fitness','epistatic_fraction','ho_fraction','slope','beneficial_slope',
 'neutral_edge_fraction','sign_epistasis_fraction','reciprocal_sign_epistasis_fraction','strict_peak_fraction','near_floor_fraction']


def summarize(stats,fits,draws=1500):
    """Cluster bootstrap: empirical cube versus synthetic seed; never noise draw."""
    rows=[];bootstrap={};keys=['dataset_id','split','observation','arm','relative_resolution']
    stats=stats.copy();stats['relative_resolution']=stats.relative_resolution.map(lambda x:min([1e-12,.001,.01],key=lambda r:abs(x-r)))
    for key,g in stats.groupby(keys,sort=False):
        dataset,split,obs,arm,relative=key;dim=int(g.dimension.iloc[0]);res=float(g.resolution.iloc[0])
        unit=g.cube if arm=='empirical' else g.seed
        groups=[x for _,x in g.groupby(unit,sort=True)]
        idx=np.random.default_rng(seed(str(key))).integers(len(groups),size=(draws,len(groups)))
        base=dict(zip(keys,key));base.update(family=family(dataset),dimension=dim,n_units=len(groups),n_cube_observations=len(g),resolution=res)
        for metric in METRICS:
            sums=np.array([x[metric].sum() for x in groups]);counts=np.array([x[metric].count() for x in groups])
            den=counts[idx].sum(axis=1);boot=np.divide(sums[idx].sum(axis=1),den,out=np.full(draws,np.nan),where=den>0)
            point=sums.sum()/counts.sum() if counts.sum() else np.nan
            finite=boot[np.isfinite(boot)];lo,hi=np.quantile(finite,[.025,.975]) if len(finite) else (np.nan,np.nan)
            rows.append(dict(**base,metric=metric,estimate=point,lower=lo,upper=hi,n_available=int(counts.sum())))
            bootstrap[(key,metric)]=boot
        for metric,orders in [('pooled_epistatic_fraction',range(2,dim+1)),('pooled_ho_fraction',range(3,dim+1))]:
            numerator=np.array([sum(x[f'power_{k}'].sum() for k in orders) for x in groups]);denominator=np.array([x.total_power.sum() for x in groups]);nobs=np.array([len(x) for x in groups])
            den=denominator[idx].sum(axis=1);floor=nobs[idx].sum(axis=1)*max(res**2,1e-24)
            boot=np.divide(numerator[idx].sum(axis=1),den,out=np.full(draws,np.nan),where=(den>floor)&(len(orders)>0))
            point=numerator.sum()/denominator.sum() if len(orders)>0 and denominator.sum()>len(g)*max(res**2,1e-24) else np.nan
            finite=boot[np.isfinite(boot)];lo,hi=np.quantile(finite,[.025,.975]) if len(finite) else (np.nan,np.nan)
            rows.append(dict(**base,metric=metric,estimate=point,lower=lo,upper=hi,n_available=int(g.power_available.sum())))
            bootstrap[(key,metric)]=boot
        for metric,values in [('median_cube_slope',g.slope),('power_available_fraction',g.power_available.astype(float)),('slope_available_fraction',g.slope_n/dim)]:
            rows.append(dict(**base,metric=metric,estimate=values.median() if metric.startswith('median') else values.mean(),lower=np.nan,upper=np.nan,n_available=int(values.count())))
    estimates=pd.DataFrame(rows);estimates.to_csv(TABLES/'aggregate_estimates.csv',index=False)
    primary=estimates[(estimates.split=='confirmation')&estimates.relative_resolution.eq(.001)]
    comparisons=[]
    for dataset in sorted(primary.dataset_id.unique()):
        empirical=primary[(primary.dataset_id==dataset)&primary.arm.eq('empirical')].set_index('metric')
        for setting in [s for s in fits if s['dataset_id']==dataset]:
            mode,arm=setting['observation'],setting['arm']
            labels=[mode] if mode=='none' else [mode,mode+'_intrinsic']
            for label in labels:
                syn=primary[(primary.dataset_id==dataset)&primary.arm.eq(arm)&primary.observation.eq(label)].set_index('metric')
                for metric in syn.index:
                    a,b=syn.loc[metric],empirical.loc[metric]
                    comparisons.append(dict(dataset_id=dataset,label=name(dataset),family=family(dataset),dimension=int(a.dimension),observation=label,arm=arm,metric=metric,
                        empirical=b.estimate,synthetic=a.estimate,empirical_lower=b.lower,empirical_upper=b.upper,synthetic_lower=a.lower,synthetic_upper=a.upper,
                        difference=a.estimate-b.estimate,ratio=a.estimate/b.estimate if b.estimate!=0 else np.nan))
    comparisons=pd.DataFrame(comparisons);comparisons.to_csv(TABLES/'primary_comparisons.csv',index=False)
    marginal=[]
    for setting in fits:
        dataset=setting['dataset_id'];mode=setting['observation'];arm=setting['arm']
        for split in ['calibration','confirmation']:
            emp=read_fitness(dataset,split,'measured','empirical');syn=read_fitness(dataset,split,mode,arm)
            probs=np.linspace(0,1,101);qerr=np.sqrt(np.mean((np.quantile(syn,probs)-np.quantile(emp,probs))**2)/emp.var())
            marginal.append(dict(dataset_id=dataset,family=family(dataset),arm=arm,observation=mode,split=split,quantile_nrmse=qerr,
                sd_ratio=syn.std(axis=1).mean()/emp.std(axis=1).mean(),calibration_loss=setting['objective']['total'],
                min_boundary_distance=min(setting['relative_boundary_distance']),optimizer_converged=setting['optimizer_converged']))
    marginal=pd.DataFrame(marginal);marginal.to_csv(TABLES/'marginal_recovery.csv',index=False)
    pd.DataFrame([{k:v for k,v in s.items() if not isinstance(v,(list,dict))} for s in fits]).to_csv(TABLES/'selected_parameters.csv',index=False)
    dump(TABLES/'portable_sigmoid_parameters.json',[s for s in fits if not s['arm'].endswith('linear')])
    return estimates,comparisons,marginal,bootstrap


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
        'axes.titlesize':12,'axes.labelsize':10,'figure.titlesize':15,'savefig.dpi':170,'pdf.fonttype':42})


def save(fig,stem):
    # Atomic replacement avoids partial cloud-synced exports.
    for suffix in ('.png','.pdf'):
        with tempfile.NamedTemporaryFile(suffix=suffix,delete=False) as temporary:
            path=Path(temporary.name)
        try:
            fig.savefig(path,bbox_inches='tight',facecolor='white')
            os.replace(path,FIGURES/(stem+suffix))
        finally:
            path.unlink(missing_ok=True)
    plt.show()


def caption(fig,text):fig.text(.01,-.025,text,fontsize=9,va='top',color='#555555')


def scope_calibration(suitability,marginal,fits):
    style();fig,axes=plt.subplots(1,3,figsize=(14,4.1),layout='constrained')
    counts=suitability.status.value_counts();labels=['Evaluated','Constant calibration','Insufficient cube coverage','Duplicate / out of scope'];statuses=['evaluated','constant_calibration','insufficient_complete_cube_coverage','out_of_scope_or_duplicate']
    axes[0].barh(labels,[counts.get(x,0) for x in statuses],color=['#0072B2','#999999','#c7ced7','#e1e5eb']);axes[0].invert_yaxis()
    for i,s in enumerate(statuses):axes[0].text(counts.get(s,0)+.7,i,str(counts.get(s,0)),va='center')
    axes[0].set(xlim=(0,70),xlabel='Inventory entries (86 total)',title='A  Coverage determines the estimand')
    data=marginal[(marginal.observation=='none')&(marginal.split=='confirmation')]
    for arm in ARMS:
        z=data[data.arm==arm];axes[1].scatter(z.calibration_loss,z.quantile_nrmse,s=24,alpha=.7,color=COLORS[arm],label=LABELS[arm])
    axes[1].set(xscale='log',yscale='log',xlabel='Calibration objective (lower is closer)',ylabel='Held-out quantile RMSE / empirical SD',title='B  Distributional fit is variable')
    axes[1].legend(fontsize=8)
    for i,group in enumerate(['Function / binding','Stability']):
        for j,arm in enumerate(ARMS):
            x=data[(data.family==group)&(data.arm==arm)].sd_ratio.to_numpy();pos=i*4+j
            axes[2].scatter(pos+np.linspace(-.12,.12,len(x)),x,s=20,alpha=.65,color=COLORS[arm])
            axes[2].plot([pos-.23,pos+.23],[np.median(x)]*2,color='black',lw=2)
    axes[2].axhline(1,color='#555555',ls='--');axes[2].set(yscale='log',xticks=[1,5],xticklabels=['Function / binding\n13 assays','Stability\n50 assays'],ylabel='Synthetic / empirical mean cube SD',title='C  Local spread remains underestimated')
    fig.suptitle('Production calibration gives useful regimes, with substantial residual differences')
    caption(fig,'Points are assays; all 63 evaluated assays are retained. Primary comparison: no added observation model.\n52 assays supply only two-site cubes; nine supply four-site and two supply three-site cubes. No pass/fail threshold is applied.')
    save(fig,'01_scope_and_calibration')


def global_statistics(comparisons):
    style();fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    specs=[('epistatic_fraction','Equal-cube epistatic fraction'),('pooled_epistatic_fraction','Variance-pooled epistatic fraction'),('ho_fraction','Equal-cube order ≥3 fraction'),('sign_epistasis_fraction','Sign epistasis on faces'),('neutral_edge_fraction','Neutral edge fraction'),('slope','Mean fixed-direction DR slope')]
    for ax,(metric,label) in zip(axes.ravel(),specs):
        for arm in ARMS:
            z=comparisons[(comparisons.observation=='none')&(comparisons.arm==arm)&(comparisons.metric==metric)].dropna(subset=['empirical','synthetic'])
            if arm=='r4_sigmoid':
                ax.hlines(z.synthetic,z.empirical_lower,z.empirical_upper,color=COLORS[arm],alpha=.18,lw=.6)
                ax.vlines(z.empirical,z.synthetic_lower,z.synthetic_upper,color=COLORS[arm],alpha=.18,lw=.6)
            for group,marker in [('Function / binding','o'),('Stability','^')]:
                d=z[z.family==group];ax.scatter(d.empirical,d.synthetic,s=28,marker=marker,color=COLORS[arm],alpha=.72,label=LABELS[arm] if group=='Function / binding' else None)
        if metric=='slope':
            ax.set_xscale('symlog',linthresh=1);ax.set_yscale('symlog',linthresh=1)
            z=comparisons[(comparisons.observation=='none')&(comparisons.metric==metric)];lo=min(z.empirical.min(),z.synthetic.min(),-1);hi=max(z.empirical.max(),z.synthetic.max(),1)
            ax.plot([lo,hi],[lo,hi],'--',color='#888888',lw=1)
        else:
            ax.plot([0,1],[0,1],'--',color='#888888',lw=1);ax.set(xlim=(-.025,1.025),ylim=(-.025,1.025))
        ax.set(title=label,xlabel='Empirical',ylabel='Synthetic')
    axes[0,0].legend(fontsize=8);fig.suptitle('Withheld global statistics: resemblance is assay- and descriptor-dependent')
    caption(fig,'Circles: function/binding; triangles: stability. Faint bars: 95% cube/seed bootstrap intervals for r=4 sigmoid. Diagonal denotes equality, not a gate.\nResolution = 0.001 × calibration fitness range. Order ≥3 is unavailable for two-site cubes. DR axes are symmetric-log to retain extreme ratios.')
    save(fig,'02_withheld_global_statistics')


def residual_atlas(comparisons):
    style();x=comparisons[(comparisons.observation=='none')&comparisons.arm.eq('r4_sigmoid')]
    datasets=sorted(x.dataset_id.unique(),key=lambda d:(family(d)=='Stability',name(d)))
    fig,axes=plt.subplots(1,6,figsize=(13,19),sharey=True,layout='constrained')
    specs=[('sd_fitness','log₂ SD ratio',2.,'ratio'),('epistatic_fraction','Δ epistasis',.6,'difference'),('ho_fraction','Δ order ≥3',.3,'difference'),('slope','Δ mean DR',2.,'difference'),('sign_epistasis_fraction','Δ sign epistasis',.6,'difference'),('neutral_edge_fraction','Δ neutral edges',.6,'difference')]
    for ax,(metric,label,limit,field) in zip(axes,specs):
        z=x[x.metric.eq(metric)].set_index('dataset_id').reindex(datasets);values=z[field].to_numpy()
        if field=='ratio':values=np.log2(values)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e6e6e6')
        im=ax.imshow(values[:,None],cmap=cmap,vmin=-limit,vmax=limit,aspect='auto',interpolation='nearest')
        for i,v in enumerate(values):
            if np.isfinite(v) and abs(v)>limit:ax.text(0,i,'+' if v>0 else '−',ha='center',va='center',fontsize=9,color='white',weight='bold')
        ax.set(title=label,xticks=[]);fig.colorbar(im,ax=ax,orientation='horizontal',shrink=.9,pad=.006,extend='both')
        ax.axhline(12.5,color='black',lw=.8)
    labels=[name(d)+f'  ({int(x[x.dataset_id.eq(d)].dimension.iloc[0])}D)' for d in datasets]
    axes[0].set_yticks(np.arange(len(datasets)),labels,fontsize=8)
    fig.suptitle('Every evaluated assay: r=4 sigmoid minus empirical (no added observation model)')
    caption(fig,'Blue = lower, red = higher; grey = unavailable. “+ / −” marks values outside the displayed colour range; exact values remain in primary_comparisons.csv.\nDifferent columns use different units. Rows are grouped by assay family, not ranked by goodness of fit. The 50 stability panels cannot establish higher-order recovery.')
    save(fig,'03_all_assay_residuals')


def examples(fits):
    # Five systems were chosen for the protocol smoke test before production;
    # GB1 is added transparently to show its seed/split instability, not as a best case.
    datasets=['gb1_wu2016_four_site_graphfla','trpb_johnston2024_four_site_graphfla','pard2_lite2020_graphfla','proteingym_amfr_human_tsuboyama_2023_4g3o','proteingym_phot_chlre_chen_2023','proteingym_spg1_strsg_olson_2014']
    style();fig,axes=plt.subplots(6,2,figsize=(11,17),layout='constrained')
    for row,dataset in enumerate(datasets):
        setting=next(s for s in fits if s['dataset_id']==dataset and s['observation']=='none' and s['arm']=='r4_sigmoid')
        origin,scale=setting['calibration_origin'],setting['calibration_scale']
        for arm in ['empirical']+ARMS:
            values=read_fitness(dataset,'confirmation','measured' if arm=='empirical' else 'none',arm)
            x=np.sort(((values-origin)/scale).ravel());p=np.arange(1,len(x)+1)/len(x)
            axes[row,0].plot(x,p,color=COLORS[arm],lw=1.5,label=LABELS[arm]);sd=np.sort(values.std(axis=1)/scale)
            axes[row,1].plot(sd,np.arange(1,len(sd)+1)/len(sd),color=COLORS[arm],lw=1.5)
        axes[row,0].set(xscale='symlog',title=name(dataset)+' — fitness',ylabel='Cumulative fraction',xlabel='Fitness in calibration affine units')
        axes[row,0].set_xscale('symlog',linthresh=.1)
        axes[row,1].set(title='Distribution of local spread',xlabel='Cube SD / calibration range',ylabel='Cumulative fraction');axes[row,1].set_xscale('symlog',linthresh=.01)
    axes[0,0].legend(fontsize=8);fig.suptitle('Measured and synthetic distributions in named assay contexts')
    caption(fig,'All confirmation vertices/cubes shown; no favorable seed or cube is selected. These are sampled-library distributions, not full protein-space distributions.\nThe five protocol smoke systems are shown alongside GB1, whose larger held-out spread illustrates unstable transfer from calibration.')
    save(fig,'04_distribution_examples')


def trpb_observation(estimates):
    style();fig,axes=plt.subplots(2,4,figsize=(17,8.6))
    fig.subplots_adjust(left=.055,right=.99,bottom=.13,top=.85,wspace=.29,hspace=.48)
    t=estimates[estimates.dataset_id.eq('trpb_johnston2024_four_site_graphfla')&estimates.split.eq('confirmation')&estimates.relative_resolution.eq(.001)]
    specs=[('epistatic_fraction','Equal-cube epistatic fraction'),('ho_fraction','Equal-cube order ≥3 fraction'),('pooled_ho_fraction','Variance-pooled order ≥3 fraction'),('slope','All-mutation DR slope'),('sd_fitness','Mean cube SD'),('sign_epistasis_fraction','Sign epistasis on faces'),('near_floor_fraction','Fraction below activity threshold'),('beneficial_slope','Beneficial-mutation DR slope')]
    modes=['none','independent_difference','shared_stop_sensitivity']
    for ax,(metric,label) in zip(axes.ravel(),specs):
        e=t[t.arm.eq('empirical')&t.metric.eq(metric)].iloc[0];ax.axhspan(e.lower,e.upper,color='#999999',alpha=.17);ax.axhline(e.estimate,color=COLORS['empirical'],lw=1.2)
        for mi,mode in enumerate(modes):
            for ai,arm in enumerate(ARMS):
                pos=mi*4+ai;r=t[t.observation.eq(mode)&t.arm.eq(arm)&t.metric.eq(metric)].iloc[0]
                ax.errorbar(pos,r.estimate,yerr=[[max(0,r.estimate-r.lower)],[max(0,r.upper-r.estimate)]],fmt='o',ms=5,color=COLORS[arm],capsize=2)
                if mode!='none':
                    intrinsic=t[t.observation.eq(mode+'_intrinsic')&t.arm.eq(arm)&t.metric.eq(metric)].iloc[0]
                    ax.plot([pos,pos],[intrinsic.estimate,r.estimate],color=COLORS[arm],lw=1)
                    ax.scatter(pos,intrinsic.estimate,facecolors='white',edgecolors=COLORS[arm],s=28,zorder=4)
        if metric in ['slope','beneficial_slope']:ax.set_yscale('symlog',linthresh=1)
        ax.set(title=label,xticks=[1,5,9],xticklabels=['No added\nnoise','Independent\ndifference','Shared-stop\nsensitivity'])
    from matplotlib.lines import Line2D
    handles=[Line2D([0],[0],marker='o',ls='',color=COLORS[a],label=LABELS[a]) for a in ARMS]
    fig.legend(handles=handles,loc='upper center',ncol=3,bbox_to_anchor=(.5,.955),fontsize=9)
    fig.suptitle('TrpB: observation assumptions materially change apparent global resemblance',y=.995)
    caption(fig,'Filled = observed; hollow = intrinsic map under the same fitted parameters. Black line/grey band = empirical estimate / 95% cube bootstrap.\nError bars resample 16 synthetic seeds, not repeated observation draws. Shared-stop covariance is a sensitivity assumption, not an identified biological mechanism.')
    save(fig,'05_trpb_observation_sensitivity')


def resolution_diagnostics(estimates,stats):
    style();fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    e=estimates[estimates.split.eq('confirmation')&estimates.arm.eq('empirical')]
    a=e[e.metric.eq('slope')].pivot(index='dataset_id',columns='relative_resolution',values='estimate')
    axes[0].scatter(a[1e-12],a[.001],s=28,c=[COLORS['r4_sigmoid'] if family(d)=='Stability' else '#555555' for d in a.index])
    axes[0].set_xscale('symlog',linthresh=1);axes[0].set_yscale('symlog',linthresh=1)
    lo=min(a[1e-12].min(),a[.001].min());hi=max(a[1e-12].max(),a[.001].max());axes[0].plot([lo,hi],[lo,hi],'--',color='#888888')
    axes[0].set(title='A  Mean DR can be denominator-sensitive',xlabel='Mean slope at relative resolution 10⁻¹²',ylabel='Mean slope at relative resolution 10⁻³')
    x=e[e.relative_resolution.eq(.001)];means=x[x.metric.eq('slope')].set_index('dataset_id').estimate;medians=x[x.metric.eq('median_cube_slope')].set_index('dataset_id').estimate
    axes[1].scatter(means,medians,s=28,color='#0072B2');axes[1].set_xscale('symlog',linthresh=1);axes[1].set_yscale('symlog',linthresh=1)
    lo=min(means.min(),medians.min());hi=max(means.max(),medians.max());axes[1].plot([lo,hi],[lo,hi],'--',color='#888888')
    axes[1].set(title='B  A few extreme ratios shift the mean',xlabel='Mean of cube-level slopes',ylabel='Median of cube-level slopes')
    z=x[x.metric.eq('slope_available_fraction')]
    for i,d in enumerate([2,3,4]):
        v=z[z.dimension.eq(d)].estimate.to_numpy();axes[2].scatter(np.full(len(v),i)+np.linspace(-.12,.12,len(v)),v,color='#0072B2',s=28)
    axes[2].set(xticks=[0,1,2],xticklabels=['2 sites\n2 backgrounds\nper mutation','3 sites\n4 backgrounds\nper mutation','4 sites\n8 backgrounds\nper mutation'],ylim=(-.03,1.03),ylabel='Fraction of mutation slopes available',title='C  Availability is part of the evidence')
    axes[2].tick_params(axis='x',labelsize=9)
    fig.suptitle('Diminishing returns is a distribution of uncertain ratios, not one universal number')
    caption(fig,'All empirical confirmation panels retained. Thresholds are numerical sensitivity choices, not known assay detection limits.\nTwo-site slopes have only two backgrounds per mutation; large or unstable means must not be interpreted as precise biological rates.')
    save(fig,'06_resolution_and_slope_diagnostics')
