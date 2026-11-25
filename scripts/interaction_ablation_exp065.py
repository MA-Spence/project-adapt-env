#!/usr/bin/env python3
"""Matched TrpB pairwise-versus-r0 fitting under the frozen shared-stop model.

RUN-186 fields, noise pools and code are immutable inputs. Only r and the four
readout coordinates are fitted. The r0 arm refits all four readout coordinates.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import shutil
import sys
import time

import numpy as np
import pandas as pd
import yaml

from interaction_ablation_metrics_exp065 import replicate_metrics, summarize_replicates, static_graph

ROOT = Path(__file__).resolve().parents[1]
EXP_REL = Path('experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration')
PROTOCOL = 'trpb_shared_stop_pairwise_ablation_v1'
MODE = 'shared_stop_sensitivity'
ARMS = ('r0', 'full')
BASELINE_MANIFEST = 'f02b8331489060b28c658ef3aa919a95fca734f81c7b73c5f3a420610b5d2e4c'
RUNTIME_MANIFEST = '39960cb6a6cc1ff99302ae1e6ba5ca1b7b38fe14a59b8b617fdeb72545c022d3'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def validate_config(config, baseline_config):
    allowed = {'experiment_id','protocol','baseline_run','baseline_manifest_sha256',
               'motivating_analysis','arms','observation','ratio_bounds','interaction_concentration',
               'fit','model','measurement','optimizer','interpretation'}
    if set(config) != allowed or config['baseline_run'] != 'RUN-186' or config['experiment_id'] != 'EXP-065':
        raise ValueError('Unexpected configuration keys or record identifiers')
    if config.get('protocol') != PROTOCOL or config.get('arms') != list(ARMS):
        raise ValueError('Require the declared paired r0/full protocol')
    if config.get('observation') != MODE or config.get('ratio_bounds') != [0., 8.]:
        raise ValueError('Keep shared-stop fixed and full r continuous on [0,8]')
    if config.get('optimizer') != dict(sampler='tpe', trials=11712, n_startup_trials=20, quick_trials=8, seed=58459495):
        raise ValueError('Unexpected ablation optimiser budget or seed')
    for key in ['fit', 'model', 'measurement']:
        if config[key] != baseline_config[key]:
            raise ValueError(f'RUN-186 {key} changed')
    if config.get('baseline_manifest_sha256') != BASELINE_MANIFEST or config.get('interaction_concentration', 'missing') is not None:
        raise ValueError('Baseline or Gaussian coefficient distribution changed')
    return config


def verify_baseline(baseline):
    """Authenticate all original scientific outputs and frozen source/input files."""
    baseline = Path(baseline)
    output, frozen = baseline/'outputs', baseline/'frozen_runtime'
    if sha(output/'manifest.json') != BASELINE_MANIFEST or sha(frozen/'prelaunch_manifest.json') != RUNTIME_MANIFEST:
        raise ValueError('RUN-186 manifest pin mismatch')
    report = {}
    for root, name in [(output, 'manifest.json'), (frozen, 'prelaunch_manifest.json')]:
        files = json.loads((root/name).read_text())['files']
        for relative, expected in files.items():
            p = root/relative
            if Path(relative).is_absolute() or '..' in Path(relative).parts or not p.resolve().is_relative_to(root.resolve()):
                raise ValueError('Invalid authenticated path')
            if not p.is_file() or sha(p) != expected:
                raise ValueError(f'Baseline checksum mismatch: {p}')
        report[name] = dict(files=len(files), sha256=sha(root/name))
    return report


def restore_noise(output):
    """Restore the exact archived donor pools; no covariance or strength refit."""
    info = json.loads((output/'noise/noise_model.json').read_text())
    table = pd.read_csv(output/'noise/donor_assignments.csv.gz')
    train = table[(table.stop_count == 0) & (table.split == 'training')]
    donors, depths, common = {}, {}, []
    for entry in info['cell_audit']:
        f, d = entry['fitness_bin'], entry['depth_bin']
        fs = train[train.fitness_bin == f]
        cell = fs[fs.depth_bin == d]
        selected = {'cell': cell, 'fitness': fs, 'all_training': train}[entry['fallback']]
        values = selected.difference.to_numpy()
        donors[f,d] = values-values.mean()
        depths[f] = (fs if entry['depth_fallback'] == 'fitness' else train).input_count.to_numpy()
        if len(values) != entry['donor_n'] or len(depths[f]) != entry['depth_pool_n']:
            raise ValueError('Archived donor membership changed')
    q = np.asarray(info['shared_covariance_raw'])
    for d in range(4):
        stop = table[(table.stop_count > 0) & (table.split == 'training') & (table.depth_bin == d)]
        means = (stop.replicate_1.to_numpy()+stop.replicate_2.to_numpy())/2
        means -= means.mean()
        common.append(means*math.sqrt(max(0.,q[d])/means.var()) if means.var() > 0 else means*0)
    return dict(donors=donors, depths=depths, common=common, q=q,
                depth_breaks=np.asarray(info['depth_breaks']), fitness_breaks=np.asarray(info['fitness_breaks'])), table


class Inputs:
    def __init__(self, baseline, config):
        self.baseline = Path(baseline).resolve()
        self.authentication = verify_baseline(self.baseline)
        self.output = self.baseline/'outputs'
        self.source = self.baseline/'frozen_runtime/source'
        self.config = validate_config(config, yaml.safe_load((self.output/'config.yaml').read_text()))
        # All reused scientific functions and the model come from the immutable run.
        sys.path.insert(0, str(self.source/'scripts'))
        self.fitting = importlib.import_module('optuna_readout_exp065')
        self.sweep = self.fitting.sweep
        for module in [self.fitting, self.sweep]:
            if not Path(module.__file__).resolve().is_relative_to(self.source):
                raise ValueError('Editable scientific module leaked into frozen replay')
        self.sub = self.source/'data/processed/interaction-strength-exp065-inputs/TrpB'
        original = yaml.safe_load((self.sub.parent/'config.yaml').read_text())
        self.cutoffs = np.asarray(original['measurement']['activity_cutoffs'])
        self.noise, self.donors = restore_noise(self.output)
        self.previous = json.loads((self.output/'selected_parameters.json').read_text())
        self.historical = [s for s in json.loads((self.output/'warm_start_parameters.json').read_text()) if s['observation'] == MODE]
        self.origin, self.scale, self.floor = [self.previous[0][k] for k in ['calibration_origin','calibration_scale','near_floor_raw']]
        with np.load(self.output/'latent_calibration.npz') as z:
            self.ca, self.cg = z['a'].copy(), z['g'].copy()
        split = pd.read_csv(self.sub/'empirical_panel_index.csv').split.eq('calibration').to_numpy()
        with np.load(self.sub/'empirical_panel.npz') as z:
            self.target = (z['fitness'][split]-self.origin)/self.scale
            self.calibration_codes = z['genotypes'][split].copy()
        with np.load(self.sub/'fit_observation_uniforms.npz') as z:
            self.bank = self.sweep.observation_bank(self.noise, z['uniforms'], True)
        with np.load(self.output/'TrpB_empirical.npz') as z:
            self.empirical, self.valid, self.codes = z['fitness'].copy(), z['valid'].copy(), z['genotypes'].copy()
        with np.load(self.output/'TrpB_empirical_cubes.npz') as z:
            self.confirmation_codes = z['genotypes'].copy()
            self.confirmation_fitness = z['fitness'].copy()
        self.audit = pd.read_csv(self.output/'seed_audit.csv')
        calibration_seeds = set(pd.read_csv(self.output/'latent_calibration_index.csv').seed)
        if calibration_seeds & set(self.audit.seed) or len(self.audit) != 24:
            raise ValueError('Invalid seed partitions')
        cal_ids = set(np.ravel_multi_index(self.calibration_codes.reshape(-1,4).T,(20,)*4))
        conf_ids = set(np.ravel_multi_index(self.confirmation_codes.reshape(-1,4).T,(20,)*4))
        if cal_ids & conf_ids:
            raise ValueError('Empirical calibration and confirmation overlap')
        self.alphabet = json.loads((self.sub/'geometry.json').read_text())['alphabets'][0]
        if self.alphabet != list('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Unexpected q20 genotype coding')
        bounds = self.fitting.parameter_bounds(config['fit'])
        self.replay = []
        for setting in self.historical:
            p = self.fitting._validate_warm_start(setting,MODE,bounds,self.origin,self.scale,self.floor)
            loss = self.loss(setting['ratio'],p)
            for key in loss:
                np.testing.assert_allclose(loss[key],setting['objective'][key],rtol=1e-9,atol=1e-10)
            self.replay.append(dict(r=setting['ratio'],objective=loss))
        if len(self.historical) != 6 or set(s['ratio'] for s in self.historical) != {0.,.5,1.,2.,4.,8.}:
            raise ValueError('Historical warm starts incomplete')

    def loss(self, ratio, parameters):
        return self.fitting.make_objective(ratio,self.ca,self.cg,self.target,self.bank,
            self.origin,self.scale,self.floor,self.config['fit'])(parameters)

    def empirical_replicates(self):
        names = [''.join(self.alphabet[int(i)] for i in code) for code in self.confirmation_codes.reshape(-1,4)]
        table = self.donors.set_index('genotype').loc[names]
        if table.split.ne('excluded').any():
            raise ValueError('Confirmation genotype used as noise donor')
        a,b = [table[name].to_numpy().reshape(-1,16) for name in ['replicate_1','replicate_2']]
        np.testing.assert_allclose((a+b)/2,self.confirmation_fitness,rtol=0,atol=1e-7)
        q = self.noise['q'][table.depth_bin.to_numpy()].reshape(-1,16)
        return a,b,q


def fit_arm(inputs, arm, directory, quick=False, null_setting=None):
    import optuna
    if arm not in ARMS or (arm == 'full') != (null_setting is not None):
        raise ValueError('Fit null first, then include it in the nested full model')
    if importlib.metadata.version('optuna') != '4.9.0':
        raise ValueError('Use frozen Optuna 4.9.0')
    cfg = inputs.config
    count = cfg['optimizer']['quick_trials' if quick else 'trials']
    names = inputs.fitting.PARAMETER_NAMES
    bounds = inputs.fitting.parameter_bounds(cfg['fit'])
    warm = [dict(r=0. if arm == 'r0' else s['ratio'],parameters=s['parameters'],
                 source=f"RUN180_readout_r{s['ratio']:g}") for s in inputs.historical]
    if arm == 'full':
        warm.append(dict(r=0.,parameters=null_setting['parameters'],source='paired_refitted_r0'))
    if count < len(warm):
        raise ValueError('Budget must include all warm starts')
    directory = Path(directory); directory.mkdir(parents=True,exist_ok=True)
    db = directory/f'{arm}.sqlite3'
    if db.exists():
        raise ValueError('Refuse to overwrite a scientific study')
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(study_name=f'TrpB_{arm}_shared_stop', direction='minimize',
        storage=f'sqlite:///{db.resolve()}',load_if_exists=False,
        sampler=optuna.samplers.TPESampler(seed=cfg['optimizer']['seed'],n_startup_trials=cfg['optimizer']['n_startup_trials']))
    for key,value in dict(protocol=PROTOCOL,arm=arm,observation=MODE,planned_trials=count,
        fixed_noise=True,noise_strength_fitted=False,interaction_concentration=None,
        coordinates=list(names) if arm == 'r0' else ['r',*names],
        warm_starts_consume_budget=True,matched_generator_seeds=True).items():
        study.set_user_attr(key,value)
    trace = []
    for item in warm:
        params = dict(zip(names,item['parameters']))
        if arm == 'full': params['r'] = item['r']
        study.enqueue_trial(params,user_attrs={'candidate':item['source']})
    def objective(trial):
        ratio = 0. if arm == 'r0' else trial.suggest_float('r',*cfg['ratio_bounds'])
        p = [trial.suggest_float(name,*bound) for name,bound in zip(names,bounds)]
        loss = inputs.loss(ratio,p)
        trial.set_user_attr('loss_components',loss)
        trace.append(dict(trial=trial.number,r=ratio,candidate=trial.user_attrs.get('candidate','optuna_tpe'),
                          **dict(zip(names,p)),**loss))
        return loss['total']
    def progress(study,trial):
        n = trial.number+1
        if n%200 == 0 or n == count:
            print('FIT_PROGRESS',arm,f'{n}/{count}',f'best_loss={study.best_value:.10g}',
                  f"best_r={study.best_params.get('r',0.):.10g}",flush=True)
    try:
        study.optimize(objective,n_trials=count,n_jobs=1,callbacks=[progress],show_progress_bar=False)
    finally:
        pd.DataFrame(trace).to_csv(directory/f'{arm}_trials.csv.gz',index=False)
    if len(study.trials) != count or any(t.state != optuna.trial.TrialState.COMPLETE for t in study.trials):
        raise ValueError('Incomplete or failed calibration study')
    p = [study.best_params[name] for name in names]
    r = study.best_params.get('r',0.)
    loss = inputs.loss(r,p)
    if arm == 'full' and loss['total'] > null_setting['objective']['total']+1e-14:
        raise AssertionError('Nested full fit worse than included null candidate')
    result = dict(arm=arm,observation=MODE,ratio=float(r),parameters=p,objective=loss,
        selected_trial=study.best_trial.number,selected_source=study.best_trial.user_attrs.get('candidate','optuna_tpe'),
        trials=count,optimizer_seed=cfg['optimizer']['seed'],bounds=bounds,
        calibration_origin=inputs.origin,calibration_scale=inputs.scale,near_floor_raw=inputs.floor,
        interaction_concentration=None,noise_strength_fitted=False,optimizer_completed=True,
        optimizer_converged=False,fit_scope='Fixed-budget calibration; no claim of global optimality',
        warm_starts=warm)
    dump(directory/f'{arm}_selected.json',result)
    print('FIT_COMPLETE',arm,result['ratio'],loss['total'],flush=True)
    return result


def diagnostics(inputs, settings, output, quick=False):
    epsilon = np.asarray(inputs.config['measurement']['relative_resolutions'])*inputs.scale
    cube_tables,face_tables,summaries,graphs = [],[],[],[]
    def append(a,b,q,label):
        cubes,faces = replicate_metrics(a,b,q,inputs.cutoffs,epsilon)
        summary = summarize_replicates(cubes,faces)
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
            raw=inputs.origin+inputs.scale*inputs.sweep.predict(setting['parameters'],latent_a,latent_g,r)
            first,second,depth,shared_q=inputs.sweep.observe(raw,inputs.noise,uniforms,shared=True)
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
    seeds.to_csv(output/'seed_audit.csv',index=False)
    dump(output/'run186_confirmation_replay.json',comparison)
    return dict(seeds=len(seeds),cube_rows=sum(map(len,cube_tables)),face_rows=sum(map(len,face_tables)),
                summary_rows=sum(map(len,summaries)),graph_rows=len(graphs))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=ROOT/EXP_REL/'interaction_ablation.yaml')
    parser.add_argument('--baseline-root',type=Path,default=ROOT/'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-186')
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args()
    start=time.perf_counter()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError('Output directory must be new or empty')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    config=yaml.safe_load(args.config.read_text())
    inputs=Inputs(args.baseline_root,config)
    inputs.empirical_replicates()
    dump(args.output_dir/'input_audit.json',dict(authentication=inputs.authentication,
        historical_objective_replays=inputs.replay,calibration_cubes=len(inputs.target),
        synthetic_calibration_cubes=len(inputs.ca),confirmation_cubes=len(inputs.confirmation_codes),
        confirmation_seeds=len(inputs.audit),noise_refitted=False,calibration_confirmation_disjoint=True))
    shutil.copyfile(args.config,args.output_dir/'config.yaml')
    dump(args.output_dir/'environment.json',dict(python=platform.python_version(),platform=platform.platform(),
        packages={name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','optuna','SQLAlchemy','PyYAML']}))
    if args.preflight:
        print('PREFLIGHT_PASS',len(inputs.replay),'historical objectives; frozen inputs and replicate mapping',flush=True)
        return
    null=fit_arm(inputs,'r0',args.output_dir/'fits',args.quick)
    full=fit_arm(inputs,'full',args.output_dir/'fits',args.quick,null)
    settings=[null,full]
    dump(args.output_dir/'selected_parameters.json',settings)
    counts=diagnostics(inputs,settings,args.output_dir,args.quick)
    files={str(p.relative_to(args.output_dir)):sha(p) for p in sorted(args.output_dir.rglob('*')) if p.is_file()}
    dump(args.output_dir/'manifest.json',dict(protocol=PROTOCOL,mode='SMOKE_ONLY' if args.quick else 'MATCHED_INTERACTION_ABLATION',
        elapsed_seconds=time.perf_counter()-start,baseline_run='RUN-186',baseline_manifest_sha256=BASELINE_MANIFEST,
        source_sha256={Path(__file__).name:sha(__file__),'interaction_ablation_metrics_exp065.py':sha(Path(__file__).with_name('interaction_ablation_metrics_exp065.py'))},
        fitted_arms=2,trials_per_arm=config['optimizer']['quick_trials' if args.quick else 'trials'],
        noise_strength_fitted=False,interaction_concentration=None,counts=counts,files=files,
        interpretation=config['interpretation']))
    print('COMPLETE',counts,flush=True)


if __name__=='__main__':
    main()
