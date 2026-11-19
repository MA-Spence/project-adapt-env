#!/usr/bin/env python3
"""Frozen RUN-166 comparison with explicitly conditional observation models.

All fitness calculations use the RUN-166 empirical affine coordinates. Neither
latent parameters nor measurement variances are fitted to aggregate discrepancies.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

H = np.array([[(-1.) ** ((i & j).bit_count()) for j in range(16)] for i in range(16)])
ORDER = np.array([i.bit_count() for i in range(16)])
LOW = np.array([[i for i in range(16) if not i & (1 << j)] for j in range(4)])
HIGH = LOW | (1 << np.arange(4))[:, None]
METHODS = ('average', 'cross', 'stop_adjusted')
STRATA = ('all', 'no_active', 'some_active', 'all_active')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def cov(x, y):
    return np.mean((x - x.mean(axis=-1, keepdims=True)) *
                   (y - y.mean(axis=-1, keepdims=True)), axis=-1)


def divide(n, d, threshold):
    out = np.full(np.broadcast_shapes(np.shape(n), np.shape(d)), np.nan)
    np.divide(n, d, out=out, where=np.asarray(d) > threshold)
    return out


def finite_mean(x, axis=-1):
    ok = np.isfinite(x)
    return np.divide(np.where(ok, x, 0).sum(axis=axis), ok.sum(axis=axis),
                     out=np.full(ok.sum(axis=axis).shape, np.nan), where=ok.sum(axis=axis) > 0)


def pair_components(rep1, rep2, shared_covariance):
    """Return signed powers and fixed-direction slope components; no clipping.

    Replicates: (cubes,16). Shared covariance is per vertex in squared fitness
    units. Correction assumes errors independent between distinct vertices.
    """
    a, b = np.asarray(rep1, float), np.asarray(rep2, float)
    if a.ndim != 2 or a.shape != b.shape or a.shape[1] != 16:
        raise ValueError('Expected two equally shaped complete four-site cubes')
    q = np.broadcast_to(shared_covariance, a.shape)
    if not all(np.isfinite(v).all() for v in (a, b, q)):
        raise ValueError('Nonfinite cube or covariance')
    mean = (a + b) / 2
    c1, c2, cm = a @ H.T / 16, b @ H.T / 16, mean @ H.T / 16
    ordinary = np.stack([(cm[:, ORDER == k] ** 2).sum(axis=1) for k in range(1, 5)], axis=1)
    cross = np.stack([(c1[:, ORDER == k] * c2[:, ORDER == k]).sum(axis=1) for k in range(1, 5)], axis=1)
    correction = q.sum(axis=1)[:, None] / 256 * np.array([4, 6, 4, 1])[None, :]
    power = np.stack([ordinary, cross, cross - correction], axis=1)
    x1, x2, xm = a[:, LOW], b[:, LOW], mean[:, LOW]
    e1, e2, em = a[:, HIGH] - x1, b[:, HIGH] - x2, mean[:, HIGH] - xm
    nc = (cov(x1, e2) + cov(x2, e1)) / 2
    dc = cov(x1, x2)
    qc = (7 / 8) * q[:, LOW].mean(axis=-1)
    numerator = np.stack([cov(xm, em), nc, nc + qc], axis=1)
    denominator = np.stack([cov(xm, xm), dc, dc - qc], axis=1)
    return dict(power=power, numerator=numerator, denominator=denominator,
                slope_correction=qc, mean_effect=em.mean(axis=-1),
                mean_fitness=mean.mean(axis=1), sd_fitness=mean.std(axis=1),
                edges=em)


def statistics(comp, metadata, resolutions, mutations=False):
    frames, mutation_frames = [], []
    n = len(metadata)
    for res in resolutions:
        threshold = max(res * res, 1e-18)
        d = metadata.copy().reset_index(drop=True)
        d['resolution'] = res
        d['mean_fitness'] = comp['mean_fitness']
        d['sd_fitness'] = comp['sd_fitness']
        d['neutral_edge_fraction'] = (np.abs(comp['edges']) <= res).mean(axis=(1, 2))
        for j, method in enumerate(METHODS):
            p = comp['power'][:, j]
            total = p.sum(axis=1)
            d[f'{method}_total_power'] = total
            d[f'{method}_power_available'] = total > threshold
            for k in range(4):
                d[f'{method}_power_{k+1}'] = p[:, k]
            d[f'{method}_epistatic_fraction'] = divide(p[:, 1:].sum(axis=1), total, threshold)
            d[f'{method}_ho_fraction'] = divide(p[:, 2:].sum(axis=1), total, threshold)
            den, num = comp['denominator'][:, j], comp['numerator'][:, j]
            slopes = divide(num, den, threshold)
            eligible = den > threshold
            beneficial = comp['mean_effect'] > res
            d[f'{method}_slope'] = finite_mean(slopes)
            d[f'{method}_slope_n'] = eligible.sum(axis=1)
            d[f'{method}_beneficial_slope'] = finite_mean(np.where(beneficial, slopes, np.nan))
            d[f'{method}_beneficial_slope_n'] = (eligible & beneficial).sum(axis=1)
            if mutations:
                m = metadata.loc[metadata.index.repeat(4)].reset_index(drop=True)
                m['mutation'] = np.tile(np.arange(4), n)
                m['resolution'], m['method'] = res, method
                m['numerator'], m['denominator'] = num.ravel(), den.ravel()
                m['shared_correction'] = comp['slope_correction'].ravel()
                m['mean_effect'] = comp['mean_effect'].ravel()
                m['eligible'], m['beneficial'] = eligible.ravel(), beneficial.ravel()
                m['slope'] = slopes.ravel()
                mutation_frames.append(m)
        frames.append(d)
    return pd.concat(frames, ignore_index=True), (pd.concat(mutation_frames, ignore_index=True) if mutations else None)


def activity(a, b, cutoffs):
    count = ((a > cutoffs[0]) & (b > cutoffs[1])).sum(axis=-1)
    return count, np.where(count == 0, 'no_active', np.where(count == 16, 'all_active', 'some_active'))


def sample_array(values, uniforms):
    return np.asarray(values)[np.minimum((uniforms * len(values)).astype(int), len(values) - 1)]


def calibrate_noise(table, excluded, config, output):
    """Split before estimating any noise parameter, retaining all assignments."""
    rng = np.random.default_rng(config['sampling']['donor_split_seed'])
    table = table.sort_values('genotype').reset_index(drop=True).copy()
    table['split'] = 'excluded'
    for mask in (table.stop_count > 0,
                 (table.stop_count == 0) & ~table.genotype.isin(excluded)):
        indices = rng.permutation(table.index[mask])
        ntrain = int(config['sampling']['training_fraction'] * len(indices))
        table.loc[indices[:ntrain], 'split'] = 'training'
        table.loc[indices[ntrain:], 'split'] = 'holdout'
    stop = table[(table.stop_count > 0) & (table.split == 'training')]
    donor = table[(table.stop_count == 0) & (table.split == 'training')]
    if set(donor.genotype) & set(excluded):
        raise ValueError('Exposed genotype leaked into noise estimation')
    depth_breaks = np.quantile(stop.input_count, config['observation']['count_quantiles'])
    cutoff = np.mean(config['measurement']['activity_cutoffs'])
    median_active = float(donor.loc[donor.published_mean > cutoff, 'published_mean'].median())
    fitness_breaks = np.array([0., cutoff, median_active])
    if not np.all(np.diff(fitness_breaks) > 0):
        raise ValueError('Degenerate training fitness strata')
    table['depth_bin'] = np.searchsorted(depth_breaks, table.input_count, side='right')
    table['fitness_bin'] = np.searchsorted(fitness_breaks, table.published_mean, side='right')
    table['difference'] = table.replicate_1 - table.replicate_2
    table.to_csv(output / 'donor_assignments.csv.gz', index=False)
    train = table[(table.stop_count == 0) & (table.split == 'training')]
    controls, common, q = [], [], []
    for depth in range(4):
        for split in ('training', 'holdout'):
            s = table[(table.stop_count > 0) & (table.split == split) & (table.depth_bin == depth)]
            a, b = s.replicate_1.to_numpy(), s.replicate_2.to_numpy()
            covariance = float(cov(a, b))
            controls.append(dict(depth_bin=depth, split=split, n=len(s),
                                 mean_1=float(a.mean()), mean_2=float(b.mean()),
                                 variance_1=float(a.var()), variance_2=float(b.var()),
                                 covariance=covariance, difference_variance=float((a-b).var())))
            if split == 'training':
                means = (a+b)/2
                means -= means.mean()
                variance = means.var()
                common.append(means * math.sqrt(max(0., covariance) / variance) if variance > 0 else means * 0)
                q.append(covariance)
    donors, depths, audits = {}, {}, []
    minimum = config['observation']['minimum_donors_per_cell']
    for f in range(4):
        fs = train[train.fitness_bin == f]
        # The same documented fallback applies if a whole fitness stratum is sparse.
        depth_pool = fs if len(fs) >= minimum else train
        depths[f] = depth_pool.input_count.to_numpy()
        for depth in range(4):
            cell = fs[fs.depth_bin == depth]
            selected, fallback = (cell, 'cell') if len(cell) >= minimum else (
                (fs, 'fitness') if len(fs) >= minimum else (train, 'all_training'))
            values = selected.difference.to_numpy()
            donors[f, depth] = values - values.mean()
            holdout = table[(table.split == 'holdout') & (table.stop_count == 0) &
                            (table.fitness_bin == f) & (table.depth_bin == depth)]
            audits.append(dict(fitness_bin=f, depth_bin=depth, cell_n=len(cell), donor_n=len(selected),
                               fallback=fallback, depth_pool_n=len(depth_pool),
                               depth_fallback='fitness' if len(fs) >= minimum else 'all_training',
                               training_difference_mean=float(values.mean()),
                               training_difference_variance=float(values.var()), holdout_n=len(holdout),
                               holdout_difference_mean=float(holdout.difference.mean()) if len(holdout) else None,
                               holdout_difference_variance=float(holdout.difference.var(ddof=0)) if len(holdout) else None))
    pd.DataFrame(controls).to_csv(output / 'control_diagnostics.csv', index=False)
    pd.DataFrame(audits).to_csv(output / 'noise_calibration.csv', index=False)
    info = dict(depth_breaks=depth_breaks.tolist(), fitness_breaks=fitness_breaks.tolist(),
                shared_covariance_raw=q, cell_audit=audits,
                membership_sha256=sha(output / 'donor_assignments.csv.gz'),
                excluded_genotypes=len(excluded), non_stop_training=len(train),
                stop_training=len(stop), conditioning='Observed donor means proxy latent fitness; not an identified error model.')
    dump(output / 'noise_model.json', info)
    return dict(donors=donors, depths=depths, common=common, q=np.array(q),
                depth_breaks=depth_breaks, fitness_breaks=fitness_breaks), table


def observe(latent_raw, noise, uniforms, shared=False):
    """Common uniforms pair arms; each replicate draws independent residuals."""
    shape = uniforms.shape[:-1]
    latent = np.broadcast_to(latent_raw, shape)
    fitness = np.searchsorted(noise['fitness_breaks'], latent, side='right')
    depth_value = np.empty(shape)
    for f in range(4):
        mask = fitness == f
        depth_value[mask] = sample_array(noise['depths'][f], uniforms[..., 0][mask])
    depth = np.searchsorted(noise['depth_breaks'], depth_value, side='right')
    error1, error2 = np.empty(shape), np.empty(shape)
    for f in range(4):
        for d in range(4):
            mask = (fitness == f) & (depth == d)
            pool = noise['donors'][f, d]
            error1[mask] = sample_array(pool, uniforms[..., 1][mask]) * np.where(uniforms[..., 2][mask] < .5, -1., 1.) / np.sqrt(2)
            error2[mask] = sample_array(pool, uniforms[..., 3][mask]) * np.where(uniforms[..., 4][mask] < .5, -1., 1.) / np.sqrt(2)
    common = np.zeros(shape)
    if shared:
        for d in range(4):
            mask = depth == d
            common[mask] = sample_array(noise['common'][d], uniforms[..., 5][mask])
    return latent + error1 + common, latent + error2 + common, depth_value, noise['q'][depth]


def aggregate_table(stats, group_cols):
    rows = []
    for keys, group in stats.groupby(group_cols, dropna=False, sort=False):
        if len(group_cols) == 1:
            keys = (keys,) if not isinstance(keys, tuple) else keys
        base = dict(zip(group_cols, keys))
        res = float(group.resolution.iloc[0])
        threshold = max(res * res, 1e-18)
        for stratum in STRATA:
            g = group if stratum == 'all' else group[group.activity_stratum == stratum]
            for method in METHODS:
                total = g[f'{method}_total_power'].sum()
                high = sum(g[f'{method}_power_{k}'].sum() for k in (3, 4))
                rows.append({**base, 'stratum': stratum, 'method': method, 'n_cubes': len(g),
                             'mean_total_power': float(total / len(g)) if len(g) else np.nan,
                             'aggregate_ho_fraction': float(high / total) if len(g) and total > threshold * len(g) else np.nan,
                             'mean_cube_ho_fraction': g[f'{method}_ho_fraction'].mean(),
                             'aggregate_epistatic_fraction': float((high + g[f'{method}_power_2'].sum()) / total) if len(g) and total > threshold * len(g) else np.nan,
                             'mean_cube_epistatic_fraction': g[f'{method}_epistatic_fraction'].mean(),
                             'power_available_n': int(g[f'{method}_power_available'].sum()),
                             'mean_slope': g[f'{method}_slope'].mean(),
                             'slope_mutations_n': int(g[f'{method}_slope_n'].sum()),
                             'mean_beneficial_slope': g[f'{method}_beneficial_slope'].mean(),
                             'beneficial_mutations_n': int(g[f'{method}_beneficial_slope_n'].sum()),
                             'mean_sd': g.sd_fitness.mean(),
                             'mean_neutral_edge_fraction': g.neutral_edge_fraction.mean()})
    return pd.DataFrame(rows)


def bootstrap_empirical(stats, draws, seed, output):
    rng = np.random.default_rng(seed)
    rows = []
    for res, group in stats[stats.split == 'confirmation'].groupby('resolution', sort=False):
        # Identical resampling across estimators, and within each activity stratum.
        indices = rng.integers(len(group), size=(draws, len(group)))
        for stratum in STRATA:
            membership = np.ones(len(group), bool) if stratum == 'all' else (group.activity_stratum == stratum).to_numpy()
            counts = membership[indices].sum(axis=1)
            for method in METHODS:
                total = group[f'{method}_total_power'].to_numpy()
                high = sum(group[f'{method}_power_{k}'].to_numpy() for k in (3, 4))
                ts = np.where(membership[indices], total[indices], 0.).sum(axis=1)
                hs = np.where(membership[indices], high[indices], 0.).sum(axis=1)
                fraction = divide(hs, ts, np.maximum(counts, 1) * max(res*res, 1e-18))
                for metric, values in [('aggregate_ho_fraction', fraction)]+[
                    (metric, finite_mean(np.where(membership[indices], group[f'{method}_{col}'].to_numpy()[indices], np.nan)))
                    for metric, col in [('mean_slope', 'slope'), ('mean_beneficial_slope', 'beneficial_slope')]]:
                    finite = values[np.isfinite(values)]
                    rows.append(dict(resolution=res, stratum=stratum, method=method, metric=metric,
                                     n_cubes=int(membership.sum()), valid_draws=len(finite), draws=draws,
                                     lower=np.quantile(finite,.025) if len(finite) else np.nan,
                                     upper=np.quantile(finite,.975) if len(finite) else np.nan))
    result = pd.DataFrame(rows)
    result.to_csv(output / 'empirical_bootstrap.csv', index=False)
    return result


def save_components(path, comp, rep1, rep2, depth, q, resolutions):
    # Cube row order is documented in the matching index table. Mutation dimension
    # is fixed site 0..3, methods average/cross/stop_adjusted; no duplicated arrays.
    np.savez_compressed(path, replicate_1_raw=rep1, replicate_2_raw=rep2,
                        input_depth=depth, shared_covariance_raw=q,
                        power=comp['power'], slope_numerator=comp['numerator'],
                        slope_denominator=comp['denominator'],
                        slope_correction=comp['slope_correction'], mean_effect=comp['mean_effect'],
                        eligible=np.stack([comp['denominator'] > max(r*r,1e-18) for r in resolutions]),
                        beneficial=np.stack([comp['mean_effect'] > r for r in resolutions]),
                        resolutions=np.array(resolutions), methods=np.array(METHODS))


def figures(output, empirical, synthetic, controls, resolutions):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    target = output / 'figures'
    target.mkdir(exist_ok=True)
    def finish(fig, name):
        fig.tight_layout()
        fig.savefig(target / (name+'.png'), dpi=160)
        fig.savefig(target / (name+'.pdf'))
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11,4))
    for split, g in controls.groupby('split'):
        axes[0].plot(g.depth_bin, g.covariance, 'o-', label=split)
        axes[1].plot(g.depth_bin, g.difference_variance/2, 'o-', label=split)
    axes[0].set(ylabel='Stop-control cross covariance (raw units²)', xlabel='Input-depth quartile')
    axes[1].set(ylabel='Stop replicate difference variance / 2', xlabel='Input-depth quartile')
    axes[0].legend(); axes[1].legend()
    fig.suptitle('Control covariance is a shared-error sensitivity, not an identified noise variance')
    finish(fig, '01_control_covariance')
    primary = resolutions[1]
    e = empirical[(empirical.resolution == primary) & (empirical.split == 'confirmation')]
    fig, axes = plt.subplots(1, 3, figsize=(14,4))
    colors = ['#1f77b4','#e69f00','#009e73']
    for j, method in enumerate(METHODS):
        g = e[e.method == method].set_index('stratum').reindex(STRATA)
        x = np.arange(4)+(j-1)*.23
        axes[0].bar(x, g.aggregate_ho_fraction, .22, label=method, color=colors[j])
        axes[1].bar(x, g.mean_beneficial_slope, .22, color=colors[j])
        axes[2].bar(x, g.beneficial_mutations_n, .22, color=colors[j])
    for ax in axes:
        ax.set_xticks(range(4), ['Full','No active','Some active','All active'], rotation=20)
        ax.axhline(0,color='.7',lw=.7)
    axes[0].set_ylabel('Aggregate signed WHT fraction, orders ≥3')
    axes[1].set_ylabel('Mean resolved beneficial slope')
    axes[2].set_ylabel('Resolved beneficial mutations (count)')
    axes[0].legend(fontsize=8)
    fig.suptitle('Same 24 empirical confirmation cubes; empty strata remain empty')
    finish(fig, '02_empirical_replicate_statistics')
    s = synthetic[(synthetic.resolution == primary) & (synthetic.method == 'average')]
    fig, axes = plt.subplots(2, 3, figsize=(14,8))
    metrics = ['aggregate_ho_fraction','mean_beneficial_slope','mean_sd']
    for row, stratum in enumerate(['all','no_active']):
        for (mode, readout), g in s[s.stratum==stratum].groupby(['mode','readout']):
            g = g.groupby('ratio',as_index=False)[metrics].mean()
            for ax,col in zip(axes[row],metrics):
                ax.plot(g.ratio,g[col], 'o-' if readout=='sigmoid' else 's--',label=f'{readout}: {mode}',ms=3)
        emp = e[(e.stratum==stratum) & (e.method=='average')].iloc[0]
        for ax,col in zip(axes[row],metrics):
            ax.axhline(emp[col],color='black',lw=1,label='Empirical average')
            ax.set(xlabel='Frozen interaction/additive ratio',ylabel=col.replace('_',' '),title=stratum.replace('_',' '))
    axes[0,0].legend(fontsize=6)
    fig.suptitle('Frozen landscapes: observation draws share the same 64 cubes per arm')
    finish(fig, '03_frozen_observation_comparison')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--quick', action='store_true')
    args = parser.parse_args()
    start = time.monotonic()
    config = yaml.safe_load(args.config.read_text())
    root = args.data_root or Path(config['data']['root'])
    if not (root/'input_manifest.json').is_file() and args.data_root is None:
        root = Path(config['data']['shared_root'])
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    inp = json.loads((root/'input_manifest.json').read_text())
    for name, expected in inp['input_sha256'].items():
        if sha(root/name) != expected:
            raise ValueError(f'Input hash mismatch: {name}')
    for member in json.loads((root/'source_manifest.json').read_text())['files']:
        # Exact downloaded archive members are preserved and checked, never executed.
        if sha(root/member['file']) != member['sha256']:
            raise ValueError(f'Source hash mismatch: {member}')
    for name in ('input_manifest.json','source_manifest.json','excluded_donor_genotypes.json'):
        shutil.copyfile(root/name, output/name)
    shutil.copyfile(args.config, output/'config.yaml')
    shutil.copyfile(root/'prior/selected_parameters.json', output/'frozen_parameters.json')
    resolutions = config['measurement']['sensitivity_resolutions']
    cutoffs = config['measurement']['activity_cutoffs']
    draws = config['sampling']['quick_observation_draws' if args.quick else 'observation_draws']
    seed = config['sampling']['quick_observation_seed' if args.quick else 'observation_seed']
    print(f"{'SMOKE_ONLY' if args.quick else 'PRODUCTION'}: frozen RUN-166, {draws} observation draws", flush=True)
    table = pd.read_csv(root/'replicates.csv')
    excluded = json.loads((root/'excluded_donor_genotypes.json').read_text())
    noise, assignments = calibrate_noise(table, excluded, config, output)
    empirical = pd.read_csv(root/'prior/empirical_cubes.csv')
    if len(empirical) != 48 or empirical.raw_scale.nunique()!=1 or empirical.raw_offset.nunique()!=1:
        raise ValueError('Unexpected empirical panel or inconsistent affine units')
    offset, scale = float(empirical.raw_offset.iloc[0]), float(empirical.raw_scale.iloc[0])
    if scale <= 0:
        raise ValueError('Nonpositive prior range')
    lookup = table.set_index('genotype')
    nodes, r1, r2, depths, means_old = [], [], [], [], []
    for _, row in empirical.iterrows():
        keys = [''.join(k) for k in json.loads(row.node_keys)]
        selected = lookup.loc[keys].reset_index()
        selected['cube_id'], selected['split'], selected['node'] = row.cube_id, row.split, np.arange(16)
        selected['prior_imported_mean'] = json.loads(row.raw_values)
        nodes.append(selected)
        r1.append(selected.replicate_1.to_numpy()); r2.append(selected.replicate_2.to_numpy())
        depths.append(selected.input_count.to_numpy()); means_old.append(selected.prior_imported_mean.to_numpy())
    r1, r2, depths = np.array(r1), np.array(r2), np.array(depths)
    if not np.allclose((r1+r2)/2, means_old, atol=1e-7, rtol=0):
        raise ValueError('Public measurements disagree with RUN-166')
    pd.concat(nodes).to_csv(output/'empirical_nodes.csv',index=False)
    meta = empirical[['cube_id','split','seed','cube']].copy()
    meta['source'], meta['mode'], meta['readout'], meta['ratio'] = 'empirical', 'published', 'empirical', np.nan
    meta['active_nodes'], meta['activity_stratum'] = activity(r1,r2,cutoffs)
    q = noise['q'][np.searchsorted(noise['depth_breaks'],depths,side='right')]
    comp = pair_components((r1-offset)/scale,(r2-offset)/scale,q/scale**2)
    estats, mutations = statistics(comp,meta,resolutions,mutations=True)
    estats.to_csv(output/'empirical_statistics.csv',index=False)
    mutations.to_csv(output/'empirical_mutations.csv',index=False)
    save_components(output/'empirical_components.npz',comp,r1,r2,depths,q,resolutions)
    esummary = aggregate_table(estats,['split','resolution'])
    esummary.to_csv(output/'empirical_summary.csv',index=False)
    bootstrap_empirical(estats,config['sampling']['quick_bootstrap_draws' if args.quick else 'bootstrap_draws'],
                        config['sampling']['bootstrap_seed'],output)
    prior = pd.read_csv(root/'prior/cube_values.csv')
    frozen = prior[(prior.phase=='active_four_site') & (prior.split=='confirmation')].copy()
    params = json.loads((root/'prior/selected_parameters.json').read_text())['within_ratio_pairs']
    if len(frozen)!=512 or frozen.groupby(['readout','ratio']).size().nunique()!=1:
        raise ValueError('Expected eight selected arms and exactly 64 cubes per arm')
    groups = list(frozen.groupby(['readout','ratio'],sort=True))
    if {(k[0],k[1]) for k,g in groups}!={(p['readout'],p['ratio']) for p in params}:
        raise ValueError('Frozen arm identity mismatch')
    arrays = output/'arrays'; arrays.mkdir(exist_ok=True)
    rng = np.random.default_rng(seed)
    uniforms = rng.random((draws,64,16,6))
    # Expose identifiers separately; npz axis 0 is flattened observation_draw,cube.
    summaries, seed_summaries, array_manifest = [], [], []
    stat_path = output/'synthetic_statistics.csv.gz'
    if stat_path.exists():
        raise FileExistsError(f'Refusing to append over an earlier run: {stat_path}')
    first = True
    reference_pairing = None
    for (readout, ratio), frame in groups:
        frame = frame.sort_values(['seed','cube']).reset_index(drop=True)
        pairing = frame[['seed','cube']].to_numpy()
        if reference_pairing is None:
            reference_pairing = pairing
        elif not np.array_equal(pairing,reference_pairing):
            raise ValueError('Frozen arms are not paired by seed and cube')
        expected = next(p for p in params if p['readout']==readout and p['ratio']==ratio)
        for key in ('offset','scale','beta','p0'):
            if expected[key] is not None and not np.allclose(frame[key],expected[key],rtol=1e-12,atol=1e-12):
                raise ValueError(f'Frozen {key} mismatch')
        latent = np.array([json.loads(s) for s in frame['values']])
        latent_raw = offset + scale*latent
        for mode in ['intrinsic']+config['observation']['modes']:
            if mode=='intrinsic':
                a,b = latent_raw.copy(), latent_raw.copy()
                dep,qraw = np.full_like(a,np.nan),np.zeros_like(a)
                index = frame[['cube_id','seed','cube','readout','ratio','beta','p0','offset','scale']].copy()
                index['observation_draw']=-1
            else:
                a,b,dep,qraw = observe(latent_raw,noise,uniforms,shared=mode=='shared_stop_sensitivity')
                a,b,dep,qraw = [v.reshape(-1,16) for v in (a,b,dep,qraw)]
                index = pd.concat([frame[['cube_id','seed','cube','readout','ratio','beta','p0','offset','scale']]]*draws,ignore_index=True)
                index['observation_draw']=np.repeat(np.arange(draws),64)
            index['source'],index['mode'],index['split']='synthetic',mode,'confirmation'
            index['active_nodes'], index['activity_stratum'] = activity(a,b,cutoffs)
            comp = pair_components((a-offset)/scale,(b-offset)/scale,qraw/scale**2)
            stats,_ = statistics(comp,index,resolutions)
            stats.to_csv(stat_path,index=False,mode='w' if first else 'a',header=first,compression='gzip')
            first=False
            stem=f'{readout}_r{ratio:g}_{mode}'
            index.to_csv(arrays/(stem+'_index.csv.gz'),index=False)
            save_components(arrays/(stem+'.npz'),comp,a,b,dep,qraw,resolutions)
            array_manifest.append(dict(stem=stem,rows=len(index),order='observation_draw then seed,cube; vertices use RUN-166 bit order',
                                       raw_units='replicates and shared_covariance_raw',
                                       transformed_units='powers, slope components and mean_effect in original RUN-166 affine units'))
            summaries.append(aggregate_table(stats,['readout','ratio','mode','observation_draw','resolution']))
            seed_summaries.append(aggregate_table(stats,['readout','ratio','mode','seed','resolution']))
        print(f'Completed frozen arm {readout}, ratio {ratio:g}',flush=True)
    summary=pd.concat(summaries,ignore_index=True)
    summary.to_csv(output/'synthetic_summary.csv',index=False)
    pd.concat(seed_summaries,ignore_index=True).to_csv(output/'seed_summary.csv',index=False)
    dump(output/'array_manifest.json',array_manifest)
    figures(output,esummary,summary,pd.read_csv(output/'control_diagnostics.csv'),resolutions)
    import scipy,matplotlib
    manifest=dict(protocol=config['protocol'],mode='SMOKE_ONLY' if args.quick else 'PRODUCTION',prior_run='RUN-166',
                  observation_draws=draws,observation_seed=seed,empirical_cubes=len(empirical),frozen_model_cubes=len(frozen),
                  input_manifest_sha256=sha(root/'input_manifest.json'),config_sha256=sha(args.config),
                  runner_sha256=sha(__file__),raw_offset=offset,raw_scale=scale,
                  versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
                                scipy=scipy.__version__,matplotlib=matplotlib.__version__,pyyaml=yaml.__version__),
                  elapsed_seconds=time.monotonic()-start,
                  limitations=['Shared error is not identifiable from replicate differences.',
                               'Stop covariance can include control biology and normalization.',
                               'Observed donor means proxy latent fitness; conditional errors need not transfer.',
                               '128 observation draws reuse 64 cubes, not independent landscapes.',
                               'Bootstrap conditions on noise calibration and omits shared-site dependence.'])
    manifest['output_sha256']={str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    dump(output/'manifest.json',manifest)
    print(f'Completed {manifest["mode"]} in {manifest["elapsed_seconds"]:.1f} seconds',flush=True)


if __name__ == '__main__':
    main()
