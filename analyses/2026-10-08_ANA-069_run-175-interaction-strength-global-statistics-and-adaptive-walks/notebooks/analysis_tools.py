"""Independent checks and figure tables for RUN-175."""
from pathlib import Path
from dataclasses import asdict
import hashlib,json,sys,tempfile
import numpy as np
import pandas as pd
import yaml
from scipy.special import expit
ANA=Path(__file__).resolve().parents[1];ROOT=ANA.parents[1]
BUNDLE=ROOT/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-175'
DATA=BUNDLE/'outputs';INPUT=BUNDLE/'inputs';SOURCE=BUNDLE/'source'
TABLES=ANA/'tables';FIGURES=ANA/'figures'
sys.dont_write_bytecode=True
sys.path.insert(0,str(SOURCE/'external/Adapt-Env'));sys.path.insert(0,str(SOURCE/'scripts'))
from adaptenv import ThermodynamicConfig,ThermodynamicLandscape
from replicate_observation_exp065 import calibrate_noise
from observation_calibration_exp065 import observation_bank,observed
IDS={'GB1':'gb1_wu2016_four_site_graphfla','TrpB':'trpb_johnston2024_four_site_graphfla'}
KEYS=['protein','seed','cohort','ratio','observation','layer','relative_resolution']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def independent_wht(y):
    coeff=y.copy();n=y.shape[1]
    for stride in [1,2,4,8]:
        shaped=coeff.reshape(len(y),-1,2*stride);left=shaped[:,:,:stride].copy();right=shaped[:,:,stride:].copy()
        shaped[:,:,:stride]=left+right;shaped[:,:,stride:]=left-right
    coeff/=n
    return np.stack([(coeff[:,[i for i in range(n) if i.bit_count()==k]]**2).sum(axis=1) for k in range(1,5)],axis=1)

def verify():
    manifest=json.loads((DATA/'manifest.json').read_text());assert manifest['mode']=='PAIRED_SWEEP'
    assert manifest['settings']==24 and manifest['graph_rows']==2598 and manifest['seed_audits']==48
    for rel,digest in manifest['files'].items():assert sha(DATA/rel)==digest,rel
    im=json.loads((INPUT/'manifest.json').read_text())
    for rel,digest in im['files'].items():assert sha(INPUT/rel)==digest,rel
    sm=json.loads((DATA/'source_manifest.json').read_text())
    for rel,digest in sm['files'].items():assert sha(SOURCE/rel)==digest,rel
    cfg=yaml.safe_load((DATA/'config.yaml').read_text());oldcfg=yaml.safe_load((INPUT/'config.yaml').read_text());settings=json.loads((DATA/'selected_parameters.json').read_text())
    oldsettings=json.loads((INPUT/'selected_parameters.json').read_text())
    g=pd.read_csv(DATA/'graph_summary.csv');c=pd.read_csv(DATA/'cube_summary.csv');seed_table=pd.read_csv(DATA/'seed_audit.csv')
    assert len(g)==len(c)==2598 and not g.duplicated(KEYS).any() and not c.duplicated(KEYS).any()
    for protein in IDS:
        s=seed_table[seed_table.protein.eq(protein)];assert s.groupby('cohort').size().to_dict()=={'fresh':8,'reused':16}
    for policy in ['uniform','greedy']:
        assert (g[policy+'_steps']>=0).all();assert g[policy+'_hit_top'].between(-1e-12,1+1e-12).all()
        assert (g[policy+'_hit_top']<=g.accessible_top+1e-12).all()
        assert (g[policy+'_hit_best']<=g[policy+'_hit_top']+1e-12).all()
    # Exact reproduction of pre-existing graph anchors on the same 16 seeds.
    prior=pd.read_csv(ROOT/'analyses/2026-10-08_ANA-068_gb1-trpb-adaptive-walks-accessibility-neutral-networks/tables/graph_summary.csv')
    replay=[]
    cols=['uniform_steps','greedy_steps','uniform_hit_top','greedy_hit_top','uniform_hit_best','accessible_top','neutral_edge_fraction','all_largest_fraction']
    for row in g[g.cohort.isin(['reused','empirical'])].itertuples():
        if row.layer=='empirical':case='empirical'
        elif row.ratio in [0,4] and row.observation=='none':case=f'r{int(row.ratio)}_none'
        elif row.ratio==4 and row.observation!='none':case='r4_'+('independent' if row.observation=='independent_difference' else 'shared')+'_'+row.layer
        else:continue
        old=prior[prior.protein.eq(row.protein)&prior.case.eq(case)&prior.seed.eq(row.seed)&prior.relative_resolution.eq(row.relative_resolution)]
        assert len(old)==1
        v=np.array([getattr(row,k) for k in cols]);w=old[cols].iloc[0].to_numpy(float)
        np.testing.assert_allclose(v,w,rtol=1e-8,atol=1e-8)
        replay.append(dict(protein=row.protein,case=case,seed=row.seed,relative_resolution=row.relative_resolution,max_absolute_difference=float(abs(v-w).max())))
    pd.DataFrame(replay).to_csv(TABLES/'anchor_replay.csv',index=False)
    with tempfile.TemporaryDirectory(prefix='strength-noise-audit-') as td:
        noise,_=calibrate_noise(pd.read_csv(INPUT/'TrpB/donor_assignments.csv.gz'),set(json.loads((INPUT/'TrpB/excluded_noise_donors.json').read_text())),oldcfg,Path(td))
    fit_checks=[];field_nodes=0;cube_checks=[]
    for protein,dataset in IDS.items():
        latent=np.load(INPUT/protein/'latent_calibration.npz');ca,cg=latent['a'],latent['g']
        panel=np.load(INPUT/protein/'empirical_panel.npz');idx=pd.read_csv(INPUT/protein/'empirical_panel_index.csv');cal=panel['fitness'][idx.split.eq('calibration')]
        if protein=='TrpB':u=np.load(INPUT/protein/'fit_observation_uniforms.npz')['uniforms'];banks={m:observation_bank(noise,u,m=='shared_stop_sensitivity') for m in ['independent_difference','shared_stop_sensitivity']}
        else:banks={}
        for s in [s for s in settings if s['protein']==protein]:
            trace=pd.read_csv(DATA/'fits'/f'{protein}_{s["observation"]}_r{s["ratio"]:g}.csv.gz')
            assert np.isclose(trace.total.min(),s['objective']['total'],rtol=1e-10,atol=1e-10)
            p=np.array(s['parameters']);origin=s['calibration_origin'];scale=s['calibration_scale'];target=(cal-origin)/scale
            raw=origin+scale*(p[0]+np.exp(p[1])*expit(p[2]+np.exp(p[3])*(ca+s['ratio']*cg)))
            if s['observation']=='none':y=(raw-origin)/scale
            else:r1,r2,_,_=observed(raw,banks[s['observation']]);y=((r1+r2)/2-origin)/scale
            probs=np.linspace(0,1,cfg['fit']['quantiles']);floor=(s['near_floor_raw']-origin)/scale
            parts={'quantile':float(np.mean((np.quantile(y,probs)-np.quantile(target,probs))**2)/target.var()),'mean_sd':float((y.std(axis=-1).mean()/target.std(axis=-1).mean()-1)**2),'floor_cdf':float(((y<=floor).mean()-(target<=floor).mean())**2)};parts['total']=sum(parts.values())
            for k,v in parts.items():np.testing.assert_allclose(v,s['objective'][k],rtol=1e-9,atol=1e-10)
            if s['ratio'] in [0,4]:
                old=next(x for x in oldsettings if x['dataset_id']==dataset and x['observation']==s['observation'] and x['arm']==f'r{int(s["ratio"])}_sigmoid')
                np.testing.assert_array_equal(p,old['parameters'])
            fit_checks.append(dict(protein=protein,ratio=s['ratio'],observation=s['observation'],objective_recomputed=parts['total'],mean_sd_ratio=float(y.std(axis=-1).mean()/target.std(axis=-1).mean()),floor_mass_model=float((y<=floor).mean()),floor_mass_empirical=float((target<=floor).mean()),optimizer_converged=s['optimizer_converged'],selected_source=s['selected_source']))
        emp=np.load(DATA/f'{protein}_empirical.npz');codes=emp['genotypes'];scale=next(s['calibration_scale'] for s in settings if s['protein']==protein)
        def check_cube_values(values,base):
            power=independent_wht(values);total=power.sum(axis=1)
            for rel in [.001,.01]:
                mask=total>max((rel*scale)**2,1e-24)
                rows=c.copy()
                for key,value in base.items():rows=rows[rows[key].eq(value)]
                rows=rows[rows.relative_resolution.eq(rel)];assert len(rows)==1
                for metric,numerator in [('epistatic_fraction',power[:,1:].sum(axis=1)),('ho_fraction',power[:,2:].sum(axis=1))]:
                    v=float((numerator[mask]/total[mask]).mean()) if mask.any() else np.nan;expected=rows[metric].iloc[0]
                    np.testing.assert_allclose(v,expected,rtol=1e-7,atol=1e-8,equal_nan=True)
                for k in range(4):
                    value=power[:,k].sum()/total.sum() if total.sum()>len(values)*max((rel*scale)**2,1e-24) else np.nan
                    np.testing.assert_allclose(value,rows[f'pooled_power_{k+1}'].iloc[0],rtol=1e-7,atol=1e-8,equal_nan=True)
                cube_checks.append(dict(**base,relative_resolution=rel,cubes=len(values)))
        cp=np.load(DATA/f'{protein}_empirical_cubes.npz');check_cube_values(cp['fitness'],dict(protein=protein,seed=-1,cohort='empirical',ratio=-1,observation='measured',layer='empirical'))
        for seedrow in seed_table[seed_table.protein.eq(protein)].itertuples():
            field=np.load(DATA/'latent_fields'/f'{protein}_{seedrow.seed}.npz');model=ThermodynamicLandscape(ThermodynamicConfig(length=4,seed=int(seedrow.seed),**cfg['model']))
            a0,g0=model.latent_components(codes[::211]);np.testing.assert_allclose(a0,field['a'][::211],rtol=1e-10,atol=1e-10);np.testing.assert_allclose(g0,field['g'][::211],rtol=1e-10,atol=1e-10);field_nodes+=len(a0)
            cv=np.load(DATA/'cube_values'/f'{protein}_{seedrow.seed}.npz')
            for name in cv.files:
                if name=='genotypes':continue
                label,r=name.rsplit('_r',1);mode,layer=label.rsplit('_',1)
                check_cube_values(cv[name],dict(protein=protein,seed=int(seedrow.seed),cohort=seedrow.cohort,ratio=float(r),observation=mode,layer=layer))
    pd.DataFrame(fit_checks).to_csv(TABLES/'fit_audit.csv',index=False);pd.DataFrame(cube_checks).to_csv(TABLES/'wht_audit.csv',index=False)
    audit={'output_hashes':len(manifest['files']),'input_hashes':len(im['files']),'source_hashes':len(sm['files']),'anchor_graph_rows_replayed':len(replay),'maximum_anchor_difference':max(x['max_absolute_difference'] for x in replay),'calibration_objectives_recomputed':len(fit_checks),'field_nodes_replayed':field_nodes,'independent_wht_summary_rows':len(cube_checks),'small_graph_cases':json.loads((DATA/'kernel_validation.json').read_text())['small_graph_cases'],'settings':24,'graph_rows':2598,'latent_graph_rows':288,'model_maps':864,'elapsed_worker_seconds':manifest['elapsed_seconds'],'scientific_gates':False}
    dump(TABLES/'audit.json',audit);return audit

def load():
    g=pd.read_csv(DATA/'graph_summary.csv');c=pd.read_csv(DATA/'cube_summary.csv');latent=pd.read_csv(DATA/'latent_graph_summary.csv')
    fits=pd.json_normalize(json.loads((DATA/'selected_parameters.json').read_text()));fits.columns=[k.replace('.','_') for k in fits.columns]
    joined=g.merge(c,on=KEYS,validate='one_to_one');joined.to_csv(TABLES/'per_seed_results.csv',index=False)
    metrics=['uniform_steps','greedy_steps','uniform_hit_top','greedy_hit_top','uniform_hit_best','accessible_top','epistatic_fraction','ho_fraction','pooled_epistatic_fraction','pooled_ho_fraction','slope','slope_median','sign_epistasis_fraction','sd_fitness','power_availability','slope_availability','above_floor_largest_fraction']
    records=[]
    for key,w in joined.groupby(['protein','cohort','ratio','observation','layer','relative_resolution'],dropna=False):
        base=dict(zip(['protein','cohort','ratio','observation','layer','relative_resolution'],key))
        for metric in metrics:
            values=w[metric].dropna();records.append(dict(**base,metric=metric,n_available=len(values),median=float(values.median()) if len(values) else np.nan,q10=float(values.quantile(.1)) if len(values) else np.nan,q90=float(values.quantile(.9)) if len(values) else np.nan))
    pd.DataFrame(records).to_csv(TABLES/'aggregate_estimates.csv',index=False)
    fits.to_csv(TABLES/'selected_parameters.csv',index=False)
    return g,c,fits,latent
