#!/usr/bin/env python3
"""Matched TrpB Gamma concentration calibration; frozen Gaussian/r0 controls."""
from __future__ import annotations
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import time
import numpy as np
import pandas as pd
import yaml
import interaction_ablation_exp065 as base
from interaction_ablation_exp065 import sha, dump, MODE, EXP_REL, ROOT
from interaction_ablation_metrics_exp065 import replicate_metrics, summarize_replicates, static_graph

PROTOCOL = 'trpb_shared_stop_gamma_concentration_v1'
BASELINE_MANIFEST = 'c40c6d47a7814bd149d07735254a2ce69cd1d70c6ae33058988a0470e9d8c306'
RUNTIME_MANIFEST = 'b2f9eccfa6b19446259c9915c23a822e28a37d86b430cce336cb85f3c10b66c1'


def verify_control(root):
    report = {}
    for folder, name, pin in [('outputs','manifest.json',BASELINE_MANIFEST),
                              ('frozen_runtime','prelaunch_manifest.json',RUNTIME_MANIFEST)]:
        directory = root/folder
        if sha(directory/name) != pin:
            raise ValueError('RUN-188 manifest pin mismatch')
        files = json.loads((directory/name).read_text())['files']
        for relative, digest in files.items():
            p = directory/relative
            if Path(relative).is_absolute() or '..' in Path(relative).parts or not p.resolve().is_relative_to(directory.resolve()):
                raise ValueError('Invalid authenticated path')
            if not p.is_file() or sha(p) != digest:
                raise ValueError(f'RUN-188 checksum mismatch: {relative}')
        report[folder] = dict(files=len(files),sha256=pin)
    return report


def validate_config(config, old):
    if set(config) != set(old)|{'concentration_bounds'}:
        raise ValueError('Unexpected Gamma configuration keys')
    expected=dict(experiment_id='EXP-065',protocol=PROTOCOL,baseline_run='RUN-188',
        baseline_manifest_sha256=BASELINE_MANIFEST,motivating_analysis='ANA-073',
        arms=['r0','full','gamma'],interaction_concentration='fitted',concentration_bounds=[.05,100.])
    for k,v in expected.items():
        if config[k] != v: raise ValueError(f'Unexpected {k}')
    for k in ['fit','model','measurement','observation','ratio_bounds']:
        if config[k] != old[k]: raise ValueError(f'Changed matched setting: {k}')
    if config['optimizer'] != dict(old['optimizer'],quick_trials=24):
        raise ValueError('Changed fixed-budget optimiser')


class Inputs(base.Inputs):
    def __init__(self, baseline, config):
        self.control_root=Path(baseline).resolve()
        self.control_authentication=verify_control(self.control_root)
        old=yaml.safe_load((self.control_root/'outputs/config.yaml').read_text())
        validate_config(config,old)
        super().__init__(self.control_root/'frozen_runtime/baseline/RUN-186',old)
        self.config=config
        self.controls=json.loads((self.control_root/'outputs/selected_parameters.json').read_text())
        if [s['arm'] for s in self.controls] != ['r0','full']:
            raise ValueError('Missing Gaussian/r0 controls')
        self.matrix=np.load(self.output/'additive_matrix.npy')
        self.calibration_index=pd.read_csv(self.output/'latent_calibration_index.csv')
        with np.load(self.output/'latent_calibration.npz') as z:
            self.synthetic_codes=z['genotypes'].copy()
        # Prove the exact archived backend reproduces both additive and Gaussian fields.
        self.field_audit=[]
        np.testing.assert_allclose(self.calibration_g(None),self.cg,rtol=1e-12,atol=1e-12)
        for s in self.controls:
            loss=self.loss(s['ratio'],s['parameters'])
            for k,v in loss.items():
                np.testing.assert_allclose(v,s['objective'][k],rtol=1e-9,atol=1e-10)
        for alpha in [.05,1.,100.]:
            g=self.calibration_g(alpha)
            if not np.isfinite(g).all(): raise ValueError('Nonfinite Gamma field')
        self.field_audit=[]

    def model(self, seed, alpha):
        cfg=self.sweep.ThermodynamicConfig(length=4,seed=int(seed),
            **self.config['model'],interaction_concentration=alpha)
        model=self.sweep.ThermodynamicLandscape(cfg,additive_matrix=self.matrix)
        np.testing.assert_allclose(model.field.expected_mutation_energy()/model.sigma_a**2,1.,rtol=1e-12,atol=1e-12)
        return model

    def calibration_g(self, alpha):
        result=np.empty_like(self.cg)
        for seed,rows in self.calibration_index.groupby('seed',sort=False):
            ix=rows.index.to_numpy()
            model=self.model(seed,alpha)
            a,g=model.latent_components(self.synthetic_codes[ix].reshape(-1,4))
            np.testing.assert_allclose(a.reshape(-1,16),self.ca[ix],rtol=0,atol=1e-12)
            result[ix]=g.reshape(-1,16)
        return result

    def gamma_loss(self, ratio, parameters, alpha):
        g=self.calibration_g(alpha)
        return self.fitting.make_objective(ratio,self.ca,g,self.target,self.bank,
            self.origin,self.scale,self.floor,self.config['fit'])(parameters)

    def confirmation_g(self,seed,alpha,additive):
        a,g=self.model(seed,alpha).latent_components(self.codes)
        np.testing.assert_allclose(a,additive,rtol=0,atol=1e-12)
        return g


def fit_gamma(inputs,directory,quick=False):
    import optuna
    if importlib.metadata.version('optuna') != '4.9.0': raise ValueError('Require Optuna 4.9.0')
    cfg=inputs.config; count=cfg['optimizer']['quick_trials' if quick else 'trials']
    names=inputs.fitting.PARAMETER_NAMES; bounds=inputs.fitting.parameter_bounds(cfg['fit'])
    warm=[dict(r=s['ratio'],alpha=1.,parameters=s['parameters'],source=f"historical_r{s['ratio']:g}_alpha1") for s in inputs.historical]
    full=next(s for s in inputs.controls if s['arm']=='full')
    warm += [dict(r=full['ratio'],alpha=alpha,parameters=full['parameters'],source=f'gaussian_winner_alpha{alpha:g}') for alpha in cfg['concentration_bounds']]
    null=next(s for s in inputs.controls if s['arm']=='r0')
    warm.append(dict(r=0.,alpha=1.,parameters=null['parameters'],source='RUN188_refitted_r0'))
    directory.mkdir(parents=True,exist_ok=True);db=directory/'gamma.sqlite3'
    if db.exists(): raise ValueError('Refuse to overwrite a scientific study')
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study=optuna.create_study(study_name='TrpB_gamma_shared_stop',direction='minimize',
        storage=f'sqlite:///{db.resolve()}',load_if_exists=False,
        sampler=optuna.samplers.TPESampler(seed=cfg['optimizer']['seed'],n_startup_trials=cfg['optimizer']['n_startup_trials']))
    for k,v in dict(protocol=PROTOCOL,coordinates=['alpha','r',*names],planned_trials=count,
        noise_strength_fitted=False,warm_starts_consume_budget=True,alpha_log=True,
        gaussian_is_separate_control=True,controls_run='RUN-188').items():study.set_user_attr(k,v)
    for w in warm:study.enqueue_trial(dict(alpha=w['alpha'],r=w['r'],**dict(zip(names,w['parameters']))),user_attrs={'candidate':w['source']})
    trace=[]
    def objective(t):
        alpha=t.suggest_float('alpha',*cfg['concentration_bounds'],log=True)
        r=t.suggest_float('r',*cfg['ratio_bounds'])
        p=[t.suggest_float(n,*b) for n,b in zip(names,bounds)]
        loss=inputs.gamma_loss(r,p,alpha);t.set_user_attr('loss_components',loss)
        trace.append(dict(trial=t.number,alpha=alpha,r=r,**dict(zip(names,p)),**loss,candidate=t.user_attrs.get('candidate','optuna_tpe')))
        return loss['total']
    def progress(study,t):
        if (t.number+1)%200==0 or t.number+1==count:
            print('FIT_PROGRESS','gamma',f'{t.number+1}/{count}',f'best_loss={study.best_value:.10g}',study.best_params,flush=True)
    try:study.optimize(objective,n_trials=count,n_jobs=1,callbacks=[progress],show_progress_bar=False)
    finally:pd.DataFrame(trace).to_csv(directory/'gamma_trials.csv.gz',index=False)
    if len(study.trials)!=count or any(t.state!=optuna.trial.TrialState.COMPLETE for t in study.trials):raise ValueError('Incomplete calibration')
    b=study.best_params;p=[b[n] for n in names];loss=inputs.gamma_loss(b['r'],p,b['alpha'])
    if loss['total']>null['objective']['total']+1e-12:raise AssertionError('Lost included null candidate')
    result=dict(arm='gamma',observation=MODE,ratio=b['r'],parameters=p,interaction_concentration=b['alpha'],
        gamma_shape=b['alpha'],gamma_rate=b['alpha'],objective=loss,selected_trial=study.best_trial.number,
        selected_source=study.best_trial.user_attrs.get('candidate','optuna_tpe'),trials=count,
        calibration_origin=inputs.origin,calibration_scale=inputs.scale,near_floor_raw=inputs.floor,
        optimizer_seed=cfg['optimizer']['seed'],noise_strength_fitted=False,optimizer_completed=True,
        optimizer_converged=False,warm_starts=warm,fit_scope='Fixed budget; no global optimality claim')
    dump(directory/'gamma_selected.json',result)
    return result


def diagnostics(inputs, settings, output, quick=False):
    epsilon = np.asarray(inputs.config['measurement']['relative_resolutions'])*inputs.scale
    cube_tables,face_tables,summaries,graphs = [],[],[],[]
    standard,distributions = [],[]
    def append(a,b,q,label):
        cubes,faces = replicate_metrics(a,b,q,inputs.cutoffs,epsilon)
        summary = summarize_replicates(cubes,faces)
        values=(a+b)/2
        stats,_,_=inputs.sweep.cube_statistics(values,pd.DataFrame({'cube':np.arange(len(values))}),epsilon)
        stats['near_floor_fraction']=np.tile((values<=inputs.floor).mean(axis=1),len(epsilon))
        for key,value in label.items():stats[key]=value
        standard.append(stats)
        distributions.append(dict(**label,mean=float(values.mean()),mean_cube_sd=float(values.std(axis=1).mean()),
            near_floor_fraction=float((values<=inputs.floor).mean()),
            **{f'q{q:03d}':float(v) for q,v in zip(range(101),np.quantile(values,np.linspace(0,1,101)))}))
        for table,collection in [(cubes,cube_tables),(faces,face_tables),(summary,summaries)]:
            for key,value in label.items():table[key]=value
            collection.append(table)
    a,b,q = inputs.empirical_replicates()
    append(a,b,q,dict(arm='empirical',seed=-1,cohort='empirical',layer='observed',ratio=-1.))
    np.savez_compressed(output/'empirical_replicates.npz',replicate_1=a,replicate_2=b,
                        shared_covariance=q,genotypes=inputs.confirmation_codes)
    for eps in epsilon:
        graphs.append(dict(arm='empirical',seed=-1,cohort='empirical',layer='observed',ratio=-1.,
                           **static_graph(inputs.empirical,inputs.valid,eps)))
    seeds = inputs.audit.groupby('cohort',sort=False).head(1) if quick else inputs.audit
    saved_dir = output/'confirmation';saved_dir.mkdir()
    (output/'latent_fields').mkdir()
    comparison = []
    for row in seeds.itertuples(index=False):
        with np.load(inputs.output/'latent_fields'/f'TrpB_{row.seed}.npz') as z:
            latent_a,latent_g = z['a'].copy(),z['g'].copy()
        with np.load(inputs.output/'cube_values'/f'TrpB_{row.seed}.npz') as z:
            codes = z['genotypes'].copy()
            archived = {key:z[key].copy() for key in z.files if key != 'genotypes'}
        ix = np.ravel_multi_index(codes.reshape(-1,4).T,(20,)*4)
        uniforms = np.random.default_rng(int(row.noise_seed)).random((len(latent_a),6))
        prior = next(s for s in inputs.previous if s['observation'] == MODE)
        raw = inputs.origin+inputs.scale*inputs.sweep.predict(prior['parameters'],latent_a,latent_g,prior['ratio'])
        r1,r2,_,_ = inputs.sweep.observe(raw,inputs.noise,uniforms,shared=True)
        for layer,values in [('intrinsic',raw),('observed',(r1+r2)/2)]:
            key=f"{MODE}_{layer}_r{prior['ratio']:g}"
            diff=float(np.max(np.abs(values[ix].reshape(-1,16)-archived[key])))
            np.testing.assert_allclose(values[ix].reshape(-1,16),archived[key],rtol=1e-11,atol=1e-12)
            comparison.append(dict(seed=int(row.seed),layer=layer,maximum_difference=diff))
        saved = dict(genotypes=codes,noise_seed=np.asarray(row.noise_seed))
        for setting in settings:
            arm,r=setting['arm'],setting['ratio']
            alpha=setting.get('interaction_concentration')
            interaction=latent_g if alpha is None else inputs.confirmation_g(row.seed,alpha,latent_a)
            raw=inputs.origin+inputs.scale*inputs.sweep.predict(setting['parameters'],latent_a,interaction,r)
            first,second,depth,shared_q=inputs.sweep.observe(raw,inputs.noise,uniforms,shared=True)
            if arm in ('r0','full'):
                with np.load(inputs.control_root/'outputs/confirmation'/f'TrpB_{row.seed}.npz') as control:
                    for key,values in [('intrinsic',raw),('replicate_1',first),('replicate_2',second),('shared_covariance',shared_q),('input_depth',depth)]:
                        np.testing.assert_allclose(values[ix].reshape(-1,16),control[f'{arm}_{key}'],rtol=1e-11,atol=1e-12)
            if alpha is not None:
                np.savez_compressed(output/'latent_fields'/f'TrpB_{row.seed}.npz',a=latent_a,g=interaction)

            saved[f'{arm}_intrinsic']=raw[ix].reshape(-1,16)
            saved[f'{arm}_replicate_1']=first[ix].reshape(-1,16)
            saved[f'{arm}_replicate_2']=second[ix].reshape(-1,16)
            saved[f'{arm}_shared_covariance']=shared_q[ix].reshape(-1,16)
            saved[f'{arm}_input_depth']=depth[ix].reshape(-1,16)
            for layer,v1,v2,covariance in [('intrinsic',raw,raw,np.zeros_like(raw)),('observed',first,second,shared_q)]:
                label=dict(arm=arm,seed=int(row.seed),cohort=row.cohort,layer=layer,ratio=r)
                append(v1[ix].reshape(-1,16),v2[ix].reshape(-1,16),covariance[ix].reshape(-1,16),label)
                for eps in epsilon:
                    graphs.append(dict(**label,**static_graph((v1+v2)/2,inputs.valid,eps)))
        np.savez_compressed(saved_dir/f'TrpB_{row.seed}.npz',**saved)
        print('CONFIRMATION',row.cohort,row.seed,'complete',flush=True)
    for name,tables in [('cube_powers.csv.gz',cube_tables),('face_reproducibility.csv.gz',face_tables),
                        ('replicate_summary.csv',summaries)]:
        pd.concat(tables,ignore_index=True).to_csv(output/name,index=False)
    pd.DataFrame(graphs).to_csv(output/'static_graph_summary.csv',index=False)
    pd.concat(standard,ignore_index=True).to_csv(output/'cube_summary.csv.gz',index=False)
    pd.DataFrame(distributions).to_csv(output/'fitness_distribution.csv',index=False)
    seeds.to_csv(output/'seed_audit.csv',index=False)
    dump(output/'run186_confirmation_replay.json',comparison)
    return dict(seeds=len(seeds),cube_rows=sum(map(len,cube_tables)),face_rows=sum(map(len,face_tables)),
                summary_rows=sum(map(len,summaries)),graph_rows=len(graphs))



def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=ROOT/EXP_REL/'gamma_concentration.yaml')
    parser.add_argument('--baseline-root',type=Path,default=ROOT/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-188')
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args();start=time.perf_counter()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):raise ValueError('Output directory must be new or empty')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    config=yaml.safe_load(args.config.read_text());inputs=Inputs(args.baseline_root,config)
    inputs.empirical_replicates()
    # Timed exact backend evaluation at a declared starting value, not model selection.
    full=next(s for s in inputs.controls if s['arm']=='full')
    t=time.perf_counter();loss=inputs.gamma_loss(full['ratio'],full['parameters'],1.)
    dump(args.output_dir/'input_audit.json',dict(control_authentication=inputs.control_authentication,
        baseline_authentication=inputs.authentication,historical_objective_replays=inputs.replay,
        gaussian_calibration_fields_reproduced=True,controls_objectives_reproduced=True,
        normalized_gamma_endpoints_verified=[.05,1.,100.],
        benchmark_alpha1_loss=loss,benchmark_seconds=time.perf_counter()-t,
        calibration_cubes=len(inputs.target),synthetic_calibration_cubes=len(inputs.ca),
        confirmation_seeds=len(inputs.audit),calibration_confirmation_disjoint=True,noise_refitted=False))
    shutil.copyfile(args.config,args.output_dir/'config.yaml')
    dump(args.output_dir/'environment.json',dict(python=platform.python_version(),platform=platform.platform(),
        packages={n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','optuna','SQLAlchemy','PyYAML']}))
    if args.preflight:
        print('PREFLIGHT_PASS frozen inputs, Gaussian replay and Gamma normalization',flush=True);return
    gamma=fit_gamma(inputs,args.output_dir/'fits',args.quick)
    settings=[dict(s,reused_from='RUN-188',refitted=False) for s in inputs.controls]+[gamma]
    dump(args.output_dir/'selected_parameters.json',settings)
    np.savez_compressed(args.output_dir/'latent_calibration.npz',a=inputs.ca,g=inputs.calibration_g(gamma['interaction_concentration']),genotypes=inputs.synthetic_codes)
    counts=diagnostics(inputs,settings,args.output_dir,args.quick)
    files={str(p.relative_to(args.output_dir)):sha(p) for p in sorted(args.output_dir.rglob('*')) if p.is_file()}
    dump(args.output_dir/'manifest.json',dict(protocol=PROTOCOL,mode='SMOKE_ONLY' if args.quick else 'MATCHED_GAMMA_CONCENTRATION',
        elapsed_seconds=time.perf_counter()-start,baseline_run='RUN-188',baseline_manifest_sha256=BASELINE_MANIFEST,
        fitted_arms=1,reused_controls=['r0','full'],trials=config['optimizer']['quick_trials' if args.quick else 'trials'],
        noise_strength_fitted=False,interaction_concentration=gamma['interaction_concentration'],counts=counts,files=files,
        interpretation=config['interpretation']))
    print('COMPLETE',counts,flush=True)


if __name__=='__main__':main()
