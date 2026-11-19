#!/usr/bin/env python3
"""EXP-066: geometry-matched, aggregate-only production landscape comparison.

Fit marginal quantiles, mean local SD and near-floor mass; withhold WHT, slopes,
DFEs and topology. No empirical interaction identities enter a synthetic field.
"""
from __future__ import annotations
import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import platform
import shutil
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
from adaptenv import ThermodynamicConfig, ThermodynamicLandscape, PUBLIC_MODEL_SYSTEM
from replicate_observation_exp065 import calibrate_noise, observe, pair_components, statistics as replicate_statistics
from observation_calibration_exp065 import observation_bank, observed


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def seed_for(name):
    return int.from_bytes(hashlib.sha256(name.encode()).digest()[:4], 'little') % 1000000000


def inventory_reason(row):
    name = row.dataset_id
    if name in {'proteingym_spg1_strsg_wu_2016', 'flip_gb1_four_mutations_full_data',
                'trpb_johnston2024_full_sequence_hf'}:
        return 'duplicate assay mirror; primary compact representation included'
    if name.startswith('dhfr_papkou'):
        return 'nine-nucleotide/codon landscape; current public model uses amino-acid substitution geometry'
    if name == 'flip_aav_full_data':
        return 'mixed insertion/deletion library; outside fixed-length substitution geometry'
    if name == 'flip_meltome_stability_summary':
        return 'cross-protein melting-temperature collection, not a single mutational landscape'
    return None


def encode_table(table):
    """Compress invariant coordinates, preserving categorical substitution geometry."""
    positions = sorted((c for c in table if c.startswith('pos') and c[3:].isdigit()), key=lambda x: int(x[3:]))
    sparse = not positions
    if positions:
        symbols = table[positions].astype(str).to_numpy()
        source_positions = positions
        reference_symbols = None
    else:
        seq = table.sequence.astype(str)
        lengths = seq.str.len().unique()
        if len(lengths) != 1:
            raise ValueError('variable-length sequence collection')
        # ASCII array avoids a millions-of-Python-character intermediate.
        symbols = np.frombuffer(''.join(seq).encode('ascii'), dtype='S1').reshape(len(seq), int(lengths[0]))
        source_positions = list(range(1, int(lengths[0]) + 1))
        refs = table.target_sequence.dropna().astype(str).unique() if 'target_sequence' in table else []
        if len(refs) != 1 or len(refs[0]) != lengths[0]:
            raise ValueError('one consistent target sequence is required for sparse library sampling')
        reference_symbols = np.frombuffer(refs[0].encode('ascii'), dtype='S1')
    variable = np.flatnonzero(np.any(symbols != symbols[0], axis=0))
    codes = np.empty((len(table), len(variable)), dtype=np.uint8)
    reference, alphabets = [], []
    for j, original in enumerate(variable):
        alphabet, inverse = np.unique(symbols[:, original], return_inverse=True)
        if len(alphabet) > 20:
            raise ValueError('more than 20 categorical alleles at a site')
        codes[:, j] = inverse
        labels = [x.decode() if isinstance(x, bytes) else str(x) for x in alphabet]
        alphabets.append(labels)
        ref = reference_symbols[original] if reference_symbols is not None else alphabet[0]
        reference.append(int(np.flatnonzero(alphabet == ref)[0]) if ref in alphabet else -1)
    if any(v < 0 for v in reference):
        raise ValueError('target residue absent from measured site alphabet')
    keys = [row.tobytes() for row in codes]
    finite = np.isfinite(table.fitness.to_numpy(float))
    # If duplicate records exist, retain their mean and record the count.
    values = pd.DataFrame({'key': keys, 'fitness': table.fitness.to_numpy(float)})[finite].groupby('key', sort=True).fitness.mean()
    keys = values.index.tolist()
    unique_codes = np.array([np.frombuffer(k, dtype=np.uint8) for k in keys], dtype=np.uint8)
    info = dict(input_rows=len(table), finite_rows=int(finite.sum()), unique_genotypes=len(keys),
                duplicate_finite_rows=int(finite.sum()-len(keys)), n_variable_sites=len(variable),
                source_positions=[source_positions[i] for i in variable], alphabets=alphabets,
                sparse=sparse, reference_codes=reference,
                geometry='Only variable library sites are modeled; no whole-protein length-transfer inference.')
    return unique_codes, values.to_numpy(), np.array(reference, dtype=np.uint8), info


def sample_cubes(codes, values, reference, sparse, cfg, rng):
    """Genotype-disjoint complete cubes, sampled without consulting fitness values."""
    ncal, nconf = cfg['calibration_cubes'], cfg['confirmation_cubes']
    lookup = {r.tobytes(): i for i, r in enumerate(codes)}
    alphabets = [np.unique(codes[:, j]) for j in range(codes.shape[1])]
    audits = []
    for dim in range(min(cfg['max_dimension'], codes.shape[1]), 1, -1):
        bits = ((np.arange(2**dim)[:, None] >> np.arange(dim)) & 1).astype(bool)
        used, panels = set(), []
        rejects = Counter()
        for attempt in range(cfg['max_attempts_per_dimension']):
            anchor = codes[int(rng.integers(len(codes)))].copy()
            candidates = np.flatnonzero(anchor != reference) if sparse else np.arange(codes.shape[1])
            if len(candidates) < dim:
                rejects['too_few_mutated_sites'] += 1; continue
            sites = np.sort(rng.choice(candidates, dim, replace=False))
            alt = []
            for site in sites:
                choices = alphabets[site][alphabets[site] != anchor[site]]
                alt.append(reference[site] if sparse and rng.random() < .5 else rng.choice(choices))
            nodes = np.tile(anchor, (2**dim, 1))
            nodes[:, sites] = np.where(bits, alt, anchor[sites])
            keys = [r.tobytes() for r in nodes]
            if any(k not in lookup for k in keys):
                rejects['missing_vertex'] += 1; continue
            if any(k in used for k in keys):
                rejects['reused_vertex'] += 1; continue
            used.update(keys); panels.append([lookup[k] for k in keys])
            if len(panels) >= ncal + nconf:
                break
        audits.append(dict(dimension=dim, attempts=attempt+1, cubes=len(panels), rejected=dict(rejects)))
        if len(panels) >= 2 * cfg['minimum_cubes_per_split']:
            # Random assignment prevents rejection/occupancy order from becoming split order.
            indices = np.array(panels)[rng.permutation(len(panels))]
            n = min(ncal, nconf, len(indices)//2)
            indices = indices[:2*n]
            return codes[indices], values[indices], dim, audits
    return None, None, None, audits


def cube_statistics(values, metadata, resolutions):
    """Exact WHT, fixed-direction DR slopes and topology on complete binary cubes."""
    y = np.asarray(values, float)
    dim = int(np.log2(y.shape[1])); width = 2**dim
    if y.ndim != 2 or y.shape[1] != width or not np.isfinite(y).all():
        raise ValueError('finite complete binary cubes required')
    H = np.array([[(-1.)**((i & j).bit_count()) for j in range(width)] for i in range(width)])
    order = np.array([i.bit_count() for i in range(width)])
    coef = y @ H.T / width
    powers = np.stack([(coef[:, order == k]**2).sum(axis=1) for k in range(1, dim+1)], axis=1)
    total = powers.sum(axis=1)
    low = np.array([[i for i in range(width) if not i & (1 << j)] for j in range(dim)])
    high = low | (1 << np.arange(dim))[:, None]
    x, effects = y[:, low], y[:, high] - y[:, low]
    xc = x-x.mean(axis=-1, keepdims=True)
    denominator = (xc**2).mean(axis=-1)
    numerator = (xc*(effects-effects.mean(axis=-1, keepdims=True))).mean(axis=-1)
    neighbors = np.arange(width)[:, None] ^ (1 << np.arange(dim))
    neighbor_delta = y[:, neighbors] - y[:, :, None]
    faces = [(j, k, np.array([i for i in range(width) if not (i & (1<<j) or i & (1<<k))]))
             for j in range(dim) for k in range(j+1, dim)]
    rows, mutation_rows = [], []
    for resolution in resolutions:
        threshold = max(resolution**2, 1e-24)
        eligible = denominator > threshold
        slopes = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=eligible)
        beneficial = effects.mean(axis=-1) > resolution
        frame = metadata.reset_index(drop=True).copy()
        frame['dimension'], frame['resolution'] = dim, resolution
        frame['mean_fitness'], frame['sd_fitness'] = y.mean(axis=1), y.std(axis=1)
        frame['total_power'], frame['power_available'] = total, total > threshold
        for k in range(1, 5):
            frame[f'power_{k}'] = powers[:, k-1] if k <= dim else np.nan
        frame['epistatic_fraction'] = np.divide(powers[:, 1:].sum(axis=1), total, out=np.full(len(y), np.nan), where=total > threshold)
        frame['ho_fraction'] = (np.divide(powers[:, 2:].sum(axis=1), total, out=np.full(len(y), np.nan), where=total > threshold) if dim >= 3 else np.nan)
        for name, mask in [('slope', eligible), ('beneficial_slope', eligible & beneficial)]:
            frame[name] = np.divide(np.where(mask, slopes, 0).sum(axis=1), mask.sum(axis=1), out=np.full(len(y), np.nan), where=mask.sum(axis=1)>0)
            frame[name+'_n'] = mask.sum(axis=1)
        frame['neutral_edge_fraction'] = (abs(effects) <= resolution).mean(axis=(1,2))
        frame['beneficial_edge_fraction'] = (effects > resolution).mean(axis=(1,2))
        frame['strict_peak_fraction'] = (neighbor_delta < -resolution).all(axis=-1).mean(axis=1)
        frame['weak_peak_fraction'] = (neighbor_delta <= resolution).all(axis=-1).mean(axis=1)
        switches, reciprocal = [], []
        for j, k, base in faces:
            ej0 = y[:, base^(1<<j)]-y[:, base]; ej1 = y[:, base^(1<<j)^(1<<k)]-y[:, base^(1<<k)]
            ek0 = y[:, base^(1<<k)]-y[:, base]; ek1 = y[:, base^(1<<j)^(1<<k)]-y[:, base^(1<<j)]
            sj = ((ej0>resolution)&(ej1 < -resolution))|((ej1>resolution)&(ej0 < -resolution))
            sk = ((ek0>resolution)&(ek1 < -resolution))|((ek1>resolution)&(ek0 < -resolution))
            switches.append(sj|sk); reciprocal.append(sj&sk)
        frame['sign_epistasis_fraction'] = np.concatenate(switches,axis=1).mean(axis=1)
        frame['reciprocal_sign_epistasis_fraction'] = np.concatenate(reciprocal,axis=1).mean(axis=1)
        rows.append(frame)
        mf = metadata.loc[metadata.index.repeat(dim)].reset_index(drop=True).copy()
        mf['mutation'] = np.tile(np.arange(dim),len(y)); mf['dimension'],mf['resolution'] = dim,resolution
        mf['numerator'],mf['denominator'] = numerator.ravel(),denominator.ravel()
        mf['slope'],mf['eligible'],mf['beneficial'] = slopes.ravel(),eligible.ravel(),beneficial.ravel()
        mf['mean_effect'] = effects.mean(axis=-1).ravel(); mutation_rows.append(mf)
    return pd.concat(rows,ignore_index=True),pd.concat(mutation_rows,ignore_index=True),effects


def latent_panel(nodes, seeds, model_cfg, cubes_per_seed):
    meta, av, gv, geometries = [], [], [], []
    for i, seed in enumerate(seeds):
        selected = (np.arange(cubes_per_seed)+i*cubes_per_seed) % len(nodes)
        geometry = nodes[selected]
        model = ThermodynamicLandscape(ThermodynamicConfig(length=nodes.shape[-1], seed=seed, **model_cfg))
        a,g = model.latent_components(geometry.reshape(-1,nodes.shape[-1]))
        av.append(a.reshape(len(selected),-1));gv.append(g.reshape(len(selected),-1));geometries.append(geometry)
        meta.extend(dict(seed=seed,cube=j,geometry_cube=int(v),fingerprint=model.fingerprint) for j,v in enumerate(selected))
    return pd.DataFrame(meta),np.concatenate(av),np.concatenate(gv),np.concatenate(geometries)


def readout(arm, p, a, g):
    latent = a + (4*g if arm.startswith('r4') else 0)
    return p[0]+np.exp(p[1])*(latent if arm.endswith('linear') else expit(p[2]+np.exp(p[3])*latent))


def fit_map(arm, a, g, target, bank, origin, scale, fit, seed, quick=False, floor_raw=None):
    # Fit in calibration-only affine coordinates; observation errors remain in raw units.
    probs = np.linspace(0,1,fit['quantiles']); q = np.quantile(target,probs)
    variance = float(target.var()); mean_sd = float(target.std(axis=-1).mean())
    floor = float(np.quantile(target,.05)) if floor_raw is None else (floor_raw-origin)/scale
    floor_mass = float((target<=floor).mean())
    if variance<=1e-24 or mean_sd<=1e-12:
        return None,None
    bounds = [fit['offset_bounds'],np.log(fit['scale_bounds']).tolist()]
    if not arm.endswith('linear'):
        bounds += [fit['position_bounds'],np.log(fit['steepness_bounds']).tolist()]
    trace = []
    def parts(p):
        raw = origin+scale*readout(arm,p,a,g)
        if bank is not None:
            r1,r2,_,_ = observed(raw,bank); y=((r1+r2)/2-origin)/scale
        else: y=(raw-origin)/scale
        quantile=float(np.mean((np.quantile(y,probs)-q)**2)/variance)
        sd=float((y.std(axis=-1).mean()/mean_sd-1)**2)
        floor_loss=float(((y<=floor).mean()-floor_mass)**2)
        return dict(quantile=quantile,mean_sd=sd,floor_cdf=floor_loss,total=quantile+sd+floor_loss)
    def objective(p):
        losses=parts(p);trace.append(dict(evaluation=len(trace),**{f'p{i}':float(v) for i,v in enumerate(p)},**losses))
        return losses['total']
    result=differential_evolution(objective,bounds,seed=seed,maxiter=2 if quick else fit['maxiter'],
        popsize=4 if quick else fit['popsize'],tol=1e-5,polish=False,workers=1)
    p=result.x
    selected=dict(arm=arm,parameters=p.tolist(),offset=origin+scale*float(p[0]),
        scale=scale*float(np.exp(p[1])),position=None if len(p)==2 else float(p[2]),
        steepness=None if len(p)==2 else float(np.exp(p[3])),interaction_ratio=4 if arm.startswith('r4') else 0,
        objective=parts(p),evaluations=result.nfev,optimizer_converged=bool(result.success),
        optimizer_message=str(result.message),bounds=bounds,
        relative_boundary_distance=[float(min(v-lo,hi-v)/(hi-lo)) for v,(lo,hi) in zip(p,bounds)],
        calibration_origin=origin,calibration_scale=scale,near_floor_raw=origin+scale*floor)
    return selected,pd.DataFrame(trace)


def config_for(setting,length,seed,model_cfg):
    return ThermodynamicConfig(length=length,seed=seed,**model_cfg,
        reference_occupancy=float(expit(setting['position'])),steepness=setting['steepness'],
        interaction_ratio=setting['interaction_ratio'],fitness_min=setting['offset'],
        fitness_max=setting['offset']+setting['scale'])


def figure_table(stats):
    """One row per metric/empirical panel or synthetic seed, retaining denominators."""
    rows=[]
    keys=['dataset_id','split','source','observation','arm','dimension','resolution','seed']
    metrics=['mean_fitness','sd_fitness','epistatic_fraction','ho_fraction','slope','beneficial_slope',
             'neutral_edge_fraction','beneficial_edge_fraction','strict_peak_fraction','weak_peak_fraction',
             'sign_epistasis_fraction','reciprocal_sign_epistasis_fraction','near_floor_fraction']
    for key,frame in stats.groupby(keys,dropna=False,sort=False):
        base=dict(zip(keys,key));den=frame.total_power.sum()
        for metric in metrics:
            available=frame[metric].dropna()
            rows.append(dict(**base,metric='mean_cube_'+metric,value=available.mean(),n_available=len(available),n_cubes=len(frame)))
        for name,orders in [('pooled_epistatic_fraction',range(2,int(base['dimension'])+1)),('pooled_ho_fraction',range(3,int(base['dimension'])+1))]:
            value=sum(frame[f'power_{k}'].sum() for k in orders)/den if den>len(frame)*max(base['resolution']**2,1e-24) and len(orders)>0 else np.nan
            rows.append(dict(**base,metric=name,value=value,n_available=int(frame.power_available.sum()),n_cubes=len(frame)))
    return pd.DataFrame(rows)


def run_assay(row,data,replicate_data,cfg,output,quick=False):
    name=row.dataset_id; sub=output/name;sub.mkdir(parents=True,exist_ok=True)
    path=data/'datasets'/name/'variants.parquet'
    if sha(path)!=row.sha256_parquet: raise ValueError(f'Input checksum mismatch: {name}')
    table=pd.read_parquet(path)
    if 'fitness_is_higher_better' in table and not table.fitness_is_higher_better.astype(bool).all():
        raise ValueError('Fitness direction must be explicitly higher-is-better')
    codes,values,reference,info=encode_table(table)
    sample_cfg=cfg['sampling'].copy()
    if quick: sample_cfg.update(calibration_cubes=4,confirmation_cubes=4,minimum_cubes_per_split=2,max_attempts_per_dimension=2000)
    nodes,y,dim,audit=sample_cubes(codes,values,reference,info['sparse'],sample_cfg,np.random.default_rng(seed_for(name)+cfg['sampling']['empirical_seed']))
    dump(sub/'geometry.json',dict(**info,sampling_audit=audit,input_sha256=sha(path)))
    if nodes is None:
        return dict(dataset_id=name,status='insufficient_complete_cube_coverage',reason='No tested dimension 2–4 supplied enough disjoint complete measured cubes within the fixed sampling budget.',**{k:info[k] for k in ('input_rows','unique_genotypes','n_variable_sites')}),[],[]
    n=len(nodes)//2
    emeta=pd.DataFrame(dict(dataset_id=name,split=['calibration']*n+['confirmation']*n,source='empirical',observation='measured',arm='empirical',seed=-1,cube=np.arange(2*n)))
    np.savez_compressed(sub/'empirical_panel.npz',genotypes=nodes,fitness=y,reference=reference)
    emeta.to_csv(sub/'empirical_panel_index.csv',index=False)
    # The normalization and numerical thresholds depend only on calibration.
    origin=float(np.quantile(y[:n],.01));scale=float(np.quantile(y[:n],.99)-origin)
    if scale<=1e-12: scale=float(y[:n].std())
    if scale<=1e-12:
        return dict(dataset_id=name,status='constant_calibration',reason='No calibration variation to identify readout parameters; retained geometry and raw measurements.',dimension=dim,n_variable_sites=nodes.shape[-1]),[],[]
    resolutions=np.array(cfg['measurement']['relative_resolutions'])*scale
    stats,mutations,edges=cube_statistics(y,emeta,resolutions)
    mutations.to_csv(sub/'empirical_mutations.csv.gz',index=False);np.savez_compressed(sub/'empirical_dfe.npz',effects=edges)
    floor_raw = float(np.mean(cfg['measurement']['activity_cutoffs'])) if name=='trpb_johnston2024_four_site_graphfla' else float(np.quantile(y[:n],.05))
    stats['near_floor_fraction']=(y<=floor_raw).mean(axis=1).tolist()*len(resolutions)
    allstats=[stats];settings=[]
    model_cfg=dict(cfg['model']['config']);model_cfg['order_weights']=tuple(model_cfg['order_weights'])
    panels={}
    for split,geometry in [('calibration',nodes[:n]),('confirmation',nodes[n:])]:
        seeds=[seed_for(name)+s for s in cfg['sampling'][split+'_seeds']]
        if quick: seeds=seeds[:1]
        meta,a,g,genotypes=latent_panel(geometry,seeds,model_cfg,2 if quick else cfg['sampling']['synthetic_cubes_per_seed'])
        meta=meta.assign(dataset_id=name,split=split,source='synthetic')
        panels[split]=(meta,a,g,genotypes)
        np.savez_compressed(sub/f'latent_{split}.npz',a=a,g=g,genotypes=genotypes)
        meta.to_csv(sub/f'latent_{split}_index.csv',index=False)
    modes=['none'];noise=None
    if name=='trpb_johnston2024_four_site_graphfla':
        table_rep=pd.read_csv(replicate_data/'replicates.csv')
        manifest=json.loads((replicate_data/'input_manifest.json').read_text())
        for filename,expected in manifest['input_sha256'].items():
            if sha(replicate_data/filename)!=expected:raise ValueError(f'Replicate input checksum mismatch: {filename}')
        genotype_keys=[''.join(info['alphabets'][j][v] for j,v in enumerate(genotype)) for genotype in nodes.reshape(-1,nodes.shape[-1])]
        excluded=set(genotype_keys)|set(json.loads((replicate_data/'excluded_donor_genotypes.json').read_text()))
        noise,_=calibrate_noise(table_rep,excluded,cfg,sub)
        dump(sub/'excluded_noise_donors.json',sorted(excluded));dump(sub/'replicate_input_hashes.json',manifest['input_sha256'])
        modes+=cfg['observation']['modes']
        # Retain paired measured values and corrected diagnostics as a separate sensitivity analysis.
        replicates=table_rep.set_index('genotype').loc[genotype_keys]
        r1=replicates.replicate_1.to_numpy().reshape(y.shape);r2=replicates.replicate_2.to_numpy().reshape(y.shape)
        # The compact assay mirror is rounded at float32 precision. Check that
        # the disagreement is only representation error, including near zero.
        rounding = float(np.finfo(np.float32).eps * max(1., np.max(np.abs(y))))
        np.testing.assert_allclose(y,(r1+r2)/2,rtol=0,atol=rounding)
        dump(sub/'replicate_mirror_precision.json',dict(
            max_absolute_difference=float(np.max(abs(y-(r1+r2)/2))),
            allowed_float32_rounding=rounding,
            primary_target='Unchanged compact-mirror fitness; paired diagnostics use replicate means.'))
        depth=replicates.input_count.to_numpy().reshape(y.shape)
        q=noise['q'][np.searchsorted(noise['depth_breaks'],depth,side='right')]
        np.savez_compressed(sub/'empirical_replicates.npz',replicate_1=r1,replicate_2=r2,depth=depth,q=q)
        if dim==4:
            corrected,cm=replicate_statistics(pair_components(r1,r2,q),emeta,resolutions,mutations=True)
            corrected.to_csv(sub/'empirical_replicate_statistics.csv',index=False);cm.to_csv(sub/'empirical_replicate_mutations.csv.gz',index=False)
    cmeta,ca,cg,_=panels['calibration'];rng=np.random.default_rng(seed_for(name)+cfg['sampling']['observation_seed'])
    fitu=rng.random((2 if quick else cfg['sampling']['fit_observation_draws'],len(cmeta),2**dim,6))
    if noise is not None:np.savez_compressed(sub/'fit_observation_uniforms.npz',uniforms=fitu)
    for mi,mode in enumerate(modes):
        bank=None if mode=='none' else observation_bank(noise,fitu,mode=='shared_stop_sensitivity')
        for ai,arm in enumerate(cfg['model']['arms']):
            setting,trace=fit_map(arm,ca,cg,(y[:n]-origin)/scale,bank,origin,scale,cfg['fit'],seed_for(name)+cfg['fit']['seed']+mi*10+ai,quick,floor_raw)
            if setting is None:continue
            setting.update(dataset_id=name,observation=mode,dimension=dim,n_variable_sites=nodes.shape[-1])
            if not arm.endswith('linear'):
                setting['thermodynamic_config']=asdict(config_for(setting,nodes.shape[-1],int(panels['confirmation'][0].seed.iloc[0]),model_cfg))
            settings.append(setting);trace.to_csv(sub/f'fit_{mode}_{arm}.csv.gz',index=False)
            print('FIT',name,mode,arm,setting['objective']['total'],flush=True)
    dump(sub/'selected_parameters.json',settings)
    for split,(meta,a,g,genotypes) in panels.items():
        uniforms=rng.random((2 if quick else cfg['sampling']['confirmation_observation_draws'],len(meta),2**dim,6))
        if noise is not None:np.savez_compressed(sub/f'observation_uniforms_{split}.npz',uniforms=uniforms)
        for setting in settings:
            mode,arm=setting['observation'],setting['arm']
            raw=origin+scale*readout(arm,setting['parameters'],a,g)
            if not arm.endswith('linear'):
                for seed,indices in meta.groupby('seed').indices.items():
                    model=ThermodynamicLandscape(config_for(setting,nodes.shape[-1],int(seed),model_cfg))
                    public=model.evaluate_batch(genotypes[indices].reshape(-1,nodes.shape[-1])).reshape(len(indices),-1)
                    np.testing.assert_allclose(raw[indices],public,rtol=1e-10,atol=1e-10*scale)
            outputs=[(mode+'_intrinsic' if mode!='none' else 'none',raw,meta)]
            if mode!='none':
                r1,r2,depth,q=observe(raw,noise,uniforms,mode=='shared_stop_sensitivity')
                obs=(r1+r2)/2;ometa=pd.concat([meta.assign(draw=i) for i in range(len(uniforms))],ignore_index=True)
                outputs.append((mode,obs.reshape(-1,2**dim),ometa))
                np.savez_compressed(sub/f'{split}_{mode}_{arm}_replicates.npz',replicate_1=r1,replicate_2=r2,depth=depth,q=q)
                if dim==4:
                    cm=ometa.assign(observation=mode,arm=arm)
                    cs,mut=replicate_statistics(pair_components(r1.reshape(-1,16),r2.reshape(-1,16),q.reshape(-1,16)),cm,resolutions,mutations=True)
                    cs.to_csv(sub/f'{split}_{mode}_{arm}_replicate_statistics.csv.gz',index=False)
            for label,fitness,metadata in outputs:
                smeta=metadata.assign(observation=label,arm=arm)
                ss,mm,ee=cube_statistics(fitness,smeta,resolutions)
                ss['near_floor_fraction']=(fitness<=setting['near_floor_raw']).mean(axis=1).tolist()*len(resolutions)
                allstats.append(ss)
                stem=f'{split}_{label}_{arm}'
                np.savez_compressed(sub/f'{stem}.npz',fitness=fitness,effects=ee)
                mm.to_csv(sub/f'{stem}_mutations.csv.gz',index=False)
    combined=pd.concat(allstats,ignore_index=True);combined.to_csv(sub/'statistics.csv.gz',index=False)
    figure_table(combined).to_csv(sub/'figure_data.csv',index=False)
    return dict(dataset_id=name,status='evaluated',reason='Complete geometry-matched held-out comparison; no scientific pass/fail threshold.',dimension=dim,n_variable_sites=nodes.shape[-1],calibration_cubes=n,confirmation_cubes=n,settings=len(settings),input_rows=info['input_rows']),settings,[combined]


def run(config,output,data,replicate_data,quick=False,only=None):
    start=time.perf_counter();output.mkdir(parents=True,exist_ok=True)
    cfg=yaml.safe_load(config.read_text());shutil.copyfile(config,output/'config.yaml')
    expected_env=json.loads((config.parent/'worker_environment.json').read_text())
    actual_env=dict(python=platform.python_version(),packages={n:importlib.metadata.version(n) for n in expected_env['packages']})
    if actual_env!=expected_env:raise ValueError(f'Worker environment differs from validated version: {actual_env}')
    dump(output/'worker_environment.json',actual_env)
    source_manifest=json.loads((config.parent/'source_manifest.json').read_text())
    for relative,expected in source_manifest['files'].items():
        if sha(ROOT/relative)!=expected:raise ValueError(f'Source differs from frozen public release/protocol: {relative}')
    shutil.copyfile(config.parent/'source_manifest.json',output/'source_manifest.json')
    index=pd.read_csv(data/'dataset_index.csv')
    if sha(data/'dataset_index.csv')!=cfg['data']['index_sha256']:raise ValueError('Dataset inventory checksum mismatch')
    index.to_csv(output/'dataset_inventory.csv',index=False)
    suitability=[];selected=[];stats=[];failures=[]
    for row in index.itertuples(index=False):
        reason=inventory_reason(row)
        if reason:
            suitability.append(dict(dataset_id=row.dataset_id,status='out_of_scope_or_duplicate',reason=reason));continue
        if only and row.dataset_id not in only:continue
        if quick and row.dataset_id not in cfg['smoke_datasets']:continue
        print('START',row.dataset_id,flush=True)
        try:
            audit,settings,frames=run_assay(row,data,replicate_data,cfg,output,quick)
            suitability.append(audit);selected+=settings;stats+=frames
        except Exception as exc:
            import traceback
            traceback.print_exc();failures.append(dict(dataset_id=row.dataset_id,error=repr(exc)))
            suitability.append(dict(dataset_id=row.dataset_id,status='execution_error',reason=repr(exc)))
        pd.DataFrame(suitability).to_csv(output/'suitability.csv',index=False)
        dump(output/'selected_parameters.json',selected)
        print('DONE',row.dataset_id,suitability[-1]['status'],flush=True)
    if stats:
        allstats=pd.concat(stats,ignore_index=True)
        allstats[allstats.split=='confirmation'].to_csv(output/'confirmation_statistics.csv',index=False)
        figure_table(allstats).to_csv(output/'figure_data.csv',index=False)
    dump(output/'execution_errors.json',failures)
    files={p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    dump(output/'manifest.json',dict(protocol=cfg['protocol'],mode='SMOKE_ONLY' if quick else 'PRODUCTION',
        public_model=PUBLIC_MODEL_SYSTEM,adaptenv_commit=source_manifest['adaptenv_commit'],
        elapsed_seconds=time.perf_counter()-start,python=platform.python_version(),numpy=np.__version__,
        config_sha256=sha(config),inventory_sha256=sha(data/'dataset_index.csv'),
        evaluated=sum(r['status']=='evaluated' for r in suitability),fitted_settings=len(selected),
        execution_errors=failures,files=files,interpretation=cfg['interpretation']))
    print('COMPLETE',json.dumps(dict(elapsed_seconds=time.perf_counter()-start,evaluated=sum(r['status']=='evaluated' for r in suitability),errors=len(failures))),flush=True)
    if failures:raise RuntimeError(f'{len(failures)} assay execution errors; see execution_errors.json. Partial results retained.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--data-root',type=Path,required=True);p.add_argument('--replicate-data-root',type=Path,required=True)
    p.add_argument('--quick',action='store_true');p.add_argument('--only',nargs='+');a=p.parse_args()
    run(a.config,a.output_dir,a.data_root,a.replicate_data_root,a.quick,a.only)
