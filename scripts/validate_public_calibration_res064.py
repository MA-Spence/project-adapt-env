#!/usr/bin/env python3
"""Replay frozen RUN-186/188/190 calibration through installed public APIs.

Run with AdaptEnv, NumPy, pandas and PyYAML available. No fitting, new draws or
changes to immutable RUN archives. Pass --output to save the validation report.
"""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from adaptenv import (aggregate_objective, LatentCalibrationPanel, ThermodynamicConfig,
                      observe_replicates, observation_bank)

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/processed/aggregate-thermodynamic-gnk-calibration'


def functions(path,names,namespace):
    tree=ast.parse(path.read_text())
    selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    if len(selected)!=len(names): raise ValueError('Frozen functions missing')
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),namespace)
    return namespace


def validate():
    run=BASE/'RUN-186';output=run/'outputs';source=run/'frozen_runtime/source'
    sub=source/'data/processed/interaction-strength-exp065-inputs/TrpB'
    # Use the archived restoration routine, without importing its runner.
    restoration=BASE/'RUN-190/frozen_runtime/source/scripts/interaction_ablation_exp065.py'
    if not restoration.is_file():
        restoration=BASE/'RUN-190/frozen_runtime/baseline/RUN-188/frozen_runtime/source/scripts/interaction_ablation_exp065.py'
    ns=functions(restoration,['restore_noise'],dict(np=np,pd=pd,json=json,math=math))
    noise,_=ns['restore_noise'](output)
    old=functions(source/'scripts/replicate_observation_exp065.py',['sample_array','observe'],dict(np=np))
    old=functions(source/'scripts/observation_calibration_exp065.py',['observation_bank'],old)
    uniforms=np.load(sub/'fit_observation_uniforms.npz')['uniforms']
    bank=observation_bank(noise,uniforms,True)
    frozen=old['observation_bank'](noise,uniforms,True)
    for x,y in zip(bank,frozen): np.testing.assert_array_equal(x,y)
    selected=json.loads((output/'selected_parameters.json').read_text())
    units={k:selected[0][k] for k in ('calibration_origin','calibration_scale','near_floor_raw')}
    origin,scale,floor=units.values()
    latent=np.load(output/'latent_calibration.npz')
    split=pd.read_csv(sub/'empirical_panel_index.csv').split.eq('calibration').to_numpy()
    target=(np.load(sub/'empirical_panel.npz')['fitness'][split]-origin)/scale
    cfg=yaml.safe_load((BASE/'RUN-190/outputs/config.yaml').read_text())
    panel=LatentCalibrationPanel(latent['genotypes'],pd.read_csv(output/'latent_calibration_index.csv').seed.to_numpy(),
        np.load(output/'additive_matrix.npy'),ThermodynamicConfig(length=4,**cfg['model']))
    np.testing.assert_allclose(panel.additive,latent['a'],rtol=0,atol=1e-12)
    np.testing.assert_allclose(panel(None),latent['g'],rtol=1e-12,atol=1e-12)
    raw=origin+scale*latent['a']
    for shared in (False,True):
        for x,y in zip(observe_replicates(raw,noise,uniforms,shared),old['observe'](raw,noise,uniforms,shared)):
            np.testing.assert_array_equal(x,y)
    records=[('RUN180_warm',s) for s in json.loads((output/'warm_start_parameters.json').read_text())
             if s['observation']=='shared_stop_sensitivity']
    records += [('RUN186',s) for s in selected if s['observation']=='shared_stop_sensitivity']
    records += [('RUN188',s) for s in json.loads((BASE/'RUN-188/outputs/selected_parameters.json').read_text())]
    records += [('RUN190',s) for s in json.loads((BASE/'RUN-190/outputs/selected_parameters.json').read_text()) if s['arm']=='gamma']
    comparisons=[]
    for label,s in records:
        alpha=s.get('interaction_concentration')
        g=panel(alpha)
        loss=aggregate_objective(s['ratio'],panel.additive,g,target,bank,origin,scale,floor,cfg['fit'])(s['parameters'])
        for k,v in loss.items(): np.testing.assert_allclose(v,s['objective'][k],rtol=1e-9,atol=1e-10)
        comparisons.append(dict(source=label,ratio=s['ratio'],alpha=alpha,loss=loss,
            max_absolute_difference=max(abs(loss[k]-s['objective'][k]) for k in loss)))
    import adaptenv.aggregate_calibration as public
    report=dict(passed=True,checks=['frozen independent/shared replicate draws exact',
        'frozen observation bank exact','Gaussian additive and interaction fields',
        'all four loss components for six warm starts and selected Gaussian/r0/Gamma settings'],
        comparisons=comparisons,public_module=str(Path(public.__file__).resolve()),
        public_module_sha256=hashlib.sha256(Path(public.__file__).read_bytes()).hexdigest(),
        source_archive='RUN-186; RUN-188; RUN-190',no_new_fit=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    report=validate()
    if args.output: args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
