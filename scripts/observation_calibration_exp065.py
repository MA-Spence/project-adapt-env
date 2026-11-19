#!/usr/bin/env python3
"""Global observation-aware calibration, with WHT and slopes strictly withheld."""
from __future__ import annotations
import argparse
import json
import platform
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from scipy.optimize import differential_evolution
from scipy.special import expit
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'external/Adapt-Env'))
from adaptenv import ThermodynamicConfig, ThermodynamicLandscape
from replicate_observation_exp065 import (sha, dump, calibrate_noise, observe, activity,
    pair_components, statistics, aggregate_table)
BITS = (np.arange(16)[:, None] >> np.arange(4)) & 1


def empirical_panel(table, excluded, cfg, quick=False, data=None):
    """Uniform measured anchors and alternative residues, conditional on completeness.

    Rejection only for missing/nonfinite measurements or previously used genotypes.
    No fitness/activity-based selection. Both panels precede noise estimation.
    """
    source = table[(table.stop_count == 0) & np.isfinite(table.replicate_1) &
                   np.isfinite(table.replicate_2)].sort_values('genotype').set_index('genotype')
    if source.index.has_duplicates:
        raise ValueError('Duplicate genotype records')
    s = cfg['sampling']; ncal = 2 if quick else s['calibration_cubes']
    nconf = 2 if quick else s['confirmation_cubes']
    states = source.index.tolist(); available = set(states); used = set(excluded)
    rng = np.random.default_rng(s['empirical_seed'])
    rows, keys = [], []
    old = pd.read_csv(data / 'prior/empirical_cubes.csv') if quick else None
    alphabet = [sorted({g[j] for g in states}) for j in range(4)]
    for attempt in range(100000):
        if quick:
            nodes = [''.join(n) for n in json.loads(old.iloc[len(keys)].node_keys)]
        else:
            anchor = states[int(rng.integers(len(states)))]
            alt = [rng.choice([a for a in alphabet[j] if a != anchor[j]]) for j in range(4)]
            nodes = [''.join(alt[j] if bits[j] else anchor[j] for j in range(4)) for bits in BITS]
            if any(g not in available or g in used for g in nodes):
                continue
        used.update(nodes); cube = len(keys); keys.append(nodes)
        split = 'calibration' if cube < ncal else 'confirmation'
        for node, genotype in enumerate(nodes):
            r = source.loc[genotype]
            rows.append(dict(split=split, cube=cube, node=node, genotype=genotype,
                replicate_1=r.replicate_1, replicate_2=r.replicate_2,
                input_count=r.input_count, published_mean=r.published_mean, attempt=attempt+1))
        if len(keys) == ncal+nconf:
            break
    if len(keys) != ncal+nconf:
        raise ValueError('Insufficient complete disjoint cubes')
    return pd.DataFrame(rows), used


def latent_panel(cfg, seeds, split, quick=False):
    records, av, gv = [], [], []
    for seed in seeds:
        m = cfg['model']
        model = ThermodynamicLandscape(ThermodynamicConfig(length=4, seed=seed,
            reference_rank=m['reference_rank'], representation_rank=m['representation_rank'],
            order_weights=tuple(m['order_weights'])))
        rng = np.random.default_rng(seed + 87111)
        for cube in range(2 if quick else cfg['sampling']['synthetic_cubes_per_seed']):
            anchor = rng.integers(20, size=4); alt = (anchor+rng.integers(1,20,size=4))%20
            nodes = np.where(BITS, alt, anchor)
            a, g = model.latent_components(nodes); av.append(a); gv.append(g)
            records.append(dict(split=split, seed=seed, cube=cube, anchor=json.dumps(anchor.tolist()),
                                alternatives=json.dumps(alt.tolist()), fingerprint=model.fingerprint))
    return pd.DataFrame(records), np.array(av), np.array(gv)


def readout(arm, params, a, g):
    b, c = params[0], np.exp(params[1])
    latent = a + (4*g if arm.startswith('r4') else 0)
    return b+c*(latent if arm.endswith('linear') else expit(params[2]+np.exp(params[3])*latent))


def observation_bank(noise, uniforms, shared):
    """Cache fixed error draws for each fitness bin; exact equivalent of observe."""
    breaks = noise['fitness_breaks']
    reps = [breaks[0]-1, (breaks[0]+breaks[1])/2, (breaks[1]+breaks[2])/2, breaks[2]+1]
    banks = []
    for value in reps:
        a, b, d, q = observe(np.full(uniforms.shape[1:-1],value),noise,uniforms,shared)
        banks.append(np.stack([a-value,b-value,d,q]))
    return np.stack(banks), breaks


def observed(raw, bank):
    arrays, breaks = bank
    f = np.searchsorted(breaks, raw, side='right')
    # arrays: fitness bin, quantity, draw, cube, node
    selected = np.take_along_axis(arrays, np.broadcast_to(f[None,None,None],
        (1,arrays.shape[1],arrays.shape[2],*f.shape)),axis=0)[0]
    return raw+selected[0], raw+selected[1], selected[2], selected[3]


def target_from(values, cfg):
    probabilities = np.linspace(0,1,cfg['fit']['quantiles'])
    cutoff = np.mean(cfg['measurement']['activity_cutoffs'])
    return dict(probabilities=probabilities, quantiles=np.quantile(values,probabilities),
        variance=float(values.var()), mean_sd=float(values.std(axis=-1).mean()),
        floor_cdf=float((values<=cutoff).mean()), cutoff=float(cutoff))


def loss_components(raw, bank, target):
    a,b,_,_ = observed(raw,bank); y=(a+b)/2
    quantile = float(np.mean((np.quantile(y,target['probabilities'])-target['quantiles'])**2)/target['variance'])
    spread = float((y.std(axis=-1).mean()/target['mean_sd']-1)**2)
    floor = float(((y<=target['cutoff']).mean()-target['floor_cdf'])**2)
    return dict(quantile=quantile,mean_sd=spread,floor_cdf=floor,total=quantile+spread+floor)


def fit_arm(arm, mode, a, g, bank, target, cfg, index, quick=False):
    f=cfg['fit']; bounds=[f['offset_bounds'], np.log(f['scale_bounds']).tolist()]
    if not arm.endswith('linear'):
        bounds += [f['position_bounds'], np.log(f['steepness_bounds']).tolist()]
    trace=[]
    def objective(p):
        parts=loss_components(readout(arm,p,a,g),bank,target)
        trace.append(dict(evaluation=len(trace), **{f'p{i}':float(v) for i,v in enumerate(p)}, **parts))
        return parts['total']
    opt=differential_evolution(objective,bounds,seed=f['seed']+index,
        maxiter=2 if quick else f['maxiter'],popsize=4 if quick else f['popsize'],
        polish=False,tol=1e-5,atol=0,workers=1,updating='immediate')
    p=opt.x
    result=dict(arm=arm,observation=mode,parameters=p.tolist(),offset=float(p[0]),scale=float(np.exp(p[1])),
        position=None if len(p)==2 else float(p[2]),steepness=None if len(p)==2 else float(np.exp(p[3])),
        occupancy=None if len(p)==2 else float(expit(p[2])),interaction_ratio=4 if arm.startswith('r4') else 0,
        objective=loss_components(readout(arm,p,a,g),bank,target),evaluations=opt.nfev,
        optimizer_converged=bool(opt.success),optimizer_message=str(opt.message),
        bounds=bounds,relative_boundary_distance=[float(min(v-lo,hi-v)/(hi-lo)) for v,(lo,hi) in zip(p,bounds)])
    if len(p)==4:
        result['thermodynamic_config']=dict(length=4,seed=cfg['sampling']['confirmation_seeds'][0],
            reference_rank=cfg['model']['reference_rank'],representation_rank=cfg['model']['representation_rank'],
            order_weights=cfg['model']['order_weights'],reference_occupancy=float(expit(p[2])),
            steepness=float(np.exp(p[3])),interaction_ratio=result['interaction_ratio'],
            fitness_min=float(p[0]),fitness_max=float(p[0]+np.exp(p[1])))
    return result,pd.DataFrame(trace)


def save_panel(output, name, rep1,rep2,depth,q,meta,cfg):
    a,b=rep1.reshape(-1,16),rep2.reshape(-1,16)
    d,qq=depth.reshape(-1,16),q.reshape(-1,16)
    count,stratum=activity(a,b,cfg['measurement']['activity_cutoffs'])
    meta=meta.copy();meta['active_vertices']=count;meta['activity_stratum']=stratum
    comp=pair_components(a,b,qq)
    stats,mutations=statistics(comp,meta,cfg['measurement']['sensitivity_resolutions'],mutations=True)
    np.savez_compressed(output/f'{name}.npz',replicate_1=a,replicate_2=b,depth=d,q=qq,**comp)
    meta.to_csv(output/f'{name}_index.csv',index=False)
    mutations.to_csv(output/f'{name}_mutations.csv.gz',index=False)
    return stats


def run(config,output,data,quick=False):
    start=time.perf_counter(); output.mkdir(parents=True,exist_ok=True)
    cfg=yaml.safe_load(config.read_text()); (output/'config.yaml').write_text(config.read_text())
    table=pd.read_csv(data/'replicates.csv'); excluded=set(json.loads((data/'excluded_donor_genotypes.json').read_text()))
    # Verify recorded inputs before any fitting.
    im=json.loads((data/'input_manifest.json').read_text())
    for name, expected in im['input_sha256'].items():
        if sha(data/name) != expected:
            raise ValueError(f'Input checksum mismatch: {name}')
    # Snapshot all small input manifests, and explicit hashes for the numerical inputs.
    inputs={p:sha(data/p) for p in ['replicates.csv','excluded_donor_genotypes.json','input_manifest.json','source_manifest.json']}
    dump(output/'input_hashes.json',inputs)
    nodes,used=empirical_panel(table,excluded,cfg,quick,data)
    nodes.to_csv(output/'empirical_nodes.csv',index=False);dump(output/'excluded_donor_genotypes.json',sorted(used))
    noise,assignments=calibrate_noise(table,used,cfg,output)
    emp1=nodes.replicate_1.to_numpy().reshape(-1,16);emp2=nodes.replicate_2.to_numpy().reshape(-1,16)
    emeta=nodes[nodes.node==0][['split','cube']].reset_index(drop=True)
    depth=nodes.input_count.to_numpy().reshape(-1,16)
    q=noise['q'][np.searchsorted(noise['depth_breaks'],depth,side='right')]
    empirical=save_panel(output,'empirical',emp1,emp2,depth,q,emeta,cfg)
    empirical.to_csv(output/'empirical_statistics.csv',index=False)
    aggregate_table(empirical,['split','resolution']).to_csv(output/'empirical_summary.csv',index=False)
    target=target_from(((emp1+emp2)/2)[emeta.split=='calibration'],cfg)
    dump(output/'calibration_target.json',{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in target.items()})
    s=cfg['sampling']; panels={};rng=np.random.default_rng(659201 if quick else s['observation_seed'])
    for split in ['calibration','confirmation']:
        seeds=([659301] if split=='calibration' else [659401]) if quick else s[f'{split}_seeds']
        meta,a,g=latent_panel(cfg,seeds,split,quick);panels[split]=(meta,a,g)
        np.savez_compressed(output/f'latent_{split}.npz',a=a,g=g)
        meta.to_csv(output/f'latent_{split}_index.csv',index=False)
    cmeta,ca,cg=panels['calibration']
    fitu=rng.random((2 if quick else s['fit_observation_draws'],len(cmeta),16,6))
    np.savez_compressed(output/'fit_uniforms.npz',uniforms=fitu)
    selected=[]
    # No confirmation arrays or statistics are passed to any optimizer.
    for index,(mode,arm) in enumerate((m,a) for m in cfg['observation']['modes'] for a in cfg['model']['arms']):
        bank=observation_bank(noise,fitu,mode=='shared_stop_sensitivity')
        result,trace=fit_arm(arm,mode,ca,cg,bank,target,cfg,index,quick)
        trace.to_csv(output/f'fit_{mode}_{arm}.csv.gz',index=False);selected.append(result)
        print('FIT',mode,arm,json.dumps(result['objective']),flush=True)
    dump(output/'selected_parameters.json',selected)
    frames=[]
    for split,(meta,a,g) in panels.items():
        u=rng.random((2 if quick else s['confirmation_observation_draws'],len(meta),16,6))
        np.savez_compressed(output/f'evaluation_uniforms_{split}.npz',uniforms=u)
        for setting in selected:
            mode,arm=setting['observation'],setting['arm'];raw=readout(arm,setting['parameters'],a,g)
            rep1,rep2,depth,q=observe(raw,noise,u,mode=='shared_stop_sensitivity')
            sm=pd.concat([meta.assign(draw=i) for i in range(len(u))],ignore_index=True)
            sm['arm'],sm['observation']=arm,mode
            name=f'{split}_{mode}_{arm}'
            frame=save_panel(output,name,rep1,rep2,depth,q,sm,cfg);frames.append(frame)
            # Intrinsic map is deterministic and retained once per parameter setting/panel.
            intrinsic=meta.assign(draw=-1,arm=arm,observation=mode+'_intrinsic')
            frames.append(save_panel(output,name+'_intrinsic',raw,raw,np.zeros_like(raw),np.zeros_like(raw),intrinsic,cfg))
    allstats=pd.concat(frames,ignore_index=True);allstats.to_csv(output/'synthetic_statistics.csv.gz',index=False)
    aggregate_table(allstats,['split','observation','arm','resolution']).to_csv(output/'synthetic_summary.csv',index=False)
    aggregate_table(allstats,['split','observation','arm','resolution','seed']).to_csv(output/'seed_summary.csv',index=False)
    files={p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    dump(output/'manifest.json',dict(protocol=cfg['protocol'],mode='SMOKE_ONLY' if quick else 'PRODUCTION',
        elapsed_seconds=time.perf_counter()-start,python=platform.python_version(),numpy=np.__version__,
        empirical_cubes=len(emeta),calibration_cubes=int((emeta.split=='calibration').sum()),
        confirmation_cubes=int((emeta.split=='confirmation').sum()),excluded_genotypes=len(used),
        fitted_settings=len(selected),files=files,inputs=inputs,
        interpretation='Aggregate resemblance on four-site TrpB cubes; no interaction recovery, unique thermodynamics, or all-length transfer claim.'))
    print('COMPLETE',json.dumps({k:v for k,v in json.loads((output/'manifest.json').read_text()).items() if k not in ('files','inputs')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--data-root',type=Path,required=True);p.add_argument('--quick',action='store_true');args=p.parse_args()
    run(args.config,args.output_dir,args.data_root,args.quick)
