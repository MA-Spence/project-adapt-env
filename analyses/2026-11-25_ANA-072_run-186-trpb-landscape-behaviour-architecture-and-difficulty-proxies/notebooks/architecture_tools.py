"""Deterministic architecture descriptors from RUN-186's immutable saved cubes.

No fitting, new landscape realization, or empirical coefficient inference occurs.
The small WHT and square-face implementation is independent of the run's code.
"""
from pathlib import Path
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ANA_NAME = '2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies'
MODES = ('none', 'independent_difference', 'shared_stop_sensitivity')
LABELS = {'none': 'No added noise · r = 4',
          'independent_difference': 'Independent difference · r = 2',
          'shared_stop_sensitivity': 'Shared-stop sensitivity · r = 0.5'}
SHORT = {'none': 'No added noise\nr = 4',
         'independent_difference': 'Independent\ndifference\nr = 2',
         'shared_stop_sensitivity': 'Shared-stop\nsensitivity\nr = 0.5'}
COLORS = {'none': '#386cb0', 'independent_difference': '#db8a20',
          'shared_stop_sensitivity': '#148578'}
STAGES = {'latent': ('#8a8b90', 'o', 'Latent additive + pairwise'),
          'intrinsic': ('#386cb0', 's', 'Intrinsic sigmoid fitness'),
          'observed': ('#c87427', '^', 'Observation augmented')}
WIDTH = 16
ORDERS = np.array([i.bit_count() for i in range(WIDTH)])
H = np.array([[(-1.) ** ((i & j).bit_count()) for j in range(WIDTH)]
              for i in range(WIDTH)])


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_geometry(codes):
    """Each bit flip must change one and only one consistent categorical site."""
    codes = np.asarray(codes)
    assert codes.ndim == 3 and codes.shape[1:] == (16, 4)
    for cube in codes:
        assert len(np.unique(cube, axis=0)) == 16
        sites = []
        for bit in range(4):
            difference = cube != cube[np.arange(16) ^ (1 << bit)]
            assert np.all(difference.sum(axis=1) == 1)
            site = np.argmax(difference, axis=1)
            assert len(np.unique(site)) == 1
            sites.append(site[0])
        assert len(set(sites)) == 4


def descriptors(values, epsilon):
    """WHT variance and face sign changes with the historical resolution rule."""
    y = np.asarray(values, dtype=float)
    assert y.ndim == 2 and y.shape[1] == WIDTH and np.isfinite(y).all()
    coefficients = y @ H.T / WIDTH
    powers = np.stack([np.sum(coefficients[:, ORDERS == k] ** 2, axis=1)
                       for k in range(1, 5)], axis=1)
    total = powers.sum(axis=1)
    np.testing.assert_allclose(total, y.var(axis=1), rtol=1e-11, atol=1e-24)
    available = total > max(epsilon ** 2, 1e-24)
    fractions = np.divide(powers, total[:, None], out=np.full_like(powers, np.nan),
                          where=available[:, None])
    equal = np.nanmean(fractions, axis=0) if available.any() else np.full(4, np.nan)
    pooled = powers.sum(axis=0) / total.sum() if total.sum() > len(y) * max(epsilon ** 2, 1e-24) else np.full(4, np.nan)
    sign, reciprocal = [], []
    for first in range(4):
        for second in range(first + 1, 4):
            start = np.array([i for i in range(16) if not i & ((1 << first) | (1 << second))])
            a = y[:, start ^ (1 << first)] - y[:, start]
            b = y[:, start ^ (1 << first) ^ (1 << second)] - y[:, start ^ (1 << second)]
            c = y[:, start ^ (1 << second)] - y[:, start]
            d = y[:, start ^ (1 << first) ^ (1 << second)] - y[:, start ^ (1 << first)]
            switch_first = ((a > epsilon) & (b < -epsilon)) | ((b > epsilon) & (a < -epsilon))
            switch_second = ((c > epsilon) & (d < -epsilon)) | ((d > epsilon) & (c < -epsilon))
            sign.append(switch_first | switch_second)
            reciprocal.append(switch_first & switch_second)
    cube_sign = np.concatenate(sign, axis=1).mean(axis=1)
    cube_reciprocal = np.concatenate(reciprocal, axis=1).mean(axis=1)
    summary = dict(n_cubes=len(y), power_availability=float(available.mean()),
                   epistatic_fraction=float(equal[1:].sum()), ho_fraction=float(equal[2:].sum()),
                   pooled_epistatic_fraction=float(pooled[1:].sum()), pooled_ho_fraction=float(pooled[2:].sum()),
                   sign_epistasis_fraction=float(cube_sign.mean()),
                   reciprocal_sign_epistasis_fraction=float(cube_reciprocal.mean()))
    summary.update({f'equal_power_{k}': float(equal[k-1]) for k in range(1, 5)})
    summary.update({f'pooled_power_{k}': float(pooled[k-1]) for k in range(1, 5)})
    detail = pd.DataFrame({'cube': np.arange(len(y)), 'total_power': total, 'power_available': available,
                           'sign_epistasis_fraction': cube_sign, 'reciprocal_sign_epistasis_fraction': cube_reciprocal})
    for k in range(1, 5):
        detail[f'power_{k}'] = powers[:, k-1]
        detail[f'fraction_{k}'] = fractions[:, k-1]
    return summary, detail


def self_check():
    rng = np.random.default_rng(186)
    coefficient = rng.normal(size=16)
    pair_coefficient = np.where(ORDERS <= 2, coefficient, 0.)
    y = pair_coefficient @ H
    summary, _ = descriptors(y[None], 0.)
    assert summary['pooled_ho_fraction'] < 1e-28
    recovered = y @ H.T / 16
    np.testing.assert_allclose(recovered, pair_coefficient, atol=1e-14)
    purely_additive = np.where(ORDERS == 1, coefficient, 0.) @ H
    additive, _ = descriptors(purely_additive[None], 0.)
    assert additive['epistatic_fraction'] < 1e-28
    assert additive['sign_epistasis_fraction'] == 0
    # A double-peak two-site landscape, duplicated over the other two axes.
    reciprocal_example = np.array([1., 0., 0., 1.] * 4)
    sign, _ = descriptors(reciprocal_example[None], 0.)
    assert sign['sign_epistasis_fraction'] == 1 / 6
    assert sign['reciprocal_sign_epistasis_fraction'] == 1 / 6
    return {'WHT_reconstruction': 'pass', 'pairwise_high_order_zero': 'pass',
            'additive_sign_epistasis_zero': 'pass', 'reciprocal_face_count': 'pass'}


def build(root):
    root = Path(root)
    ana = root / 'analyses' / ANA_NAME
    tables = ana / 'tables'
    tables.mkdir(exist_ok=True)
    output = root / 'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-186/outputs'
    manifest = json.loads((output / 'manifest.json').read_text())
    used = ['TrpB_empirical_cubes.npz', 'selected_parameters.json', 'seed_audit.csv', 'coverage.csv', 'cube_summary.csv']
    used += [str(p.relative_to(output)) for folder in ['cube_values', 'latent_fields']
             for p in sorted((output / folder).glob('*.npz'))]
    for rel in used:
        assert _sha(output / rel) == manifest['files'][rel], rel
    settings = json.loads((output / 'selected_parameters.json').read_text())
    audit = pd.read_csv(output / 'seed_audit.csv')
    scale = float(pd.read_csv(output / 'coverage.csv').iloc[0].scale)
    rows, cube_rows = [], []
    def append(values, label, relative, native_scale):
        summary, detail = descriptors(values, relative * native_scale)
        label = dict(label, relative_resolution=relative, epsilon=relative * native_scale)
        rows.append(dict(label, **summary))
        for key, value in label.items():
            detail[key] = value
        cube_rows.append(detail)
    with np.load(output / 'TrpB_empirical_cubes.npz') as z:
        empirical_codes, empirical_y = z['genotypes'], z['fitness']
    validate_geometry(empirical_codes)
    for resolution in [0., .001, .01]:
        append(empirical_y, dict(seed=-1, cohort='empirical', observation='measured', ratio=-1., layer='empirical'), resolution, scale)
    latent_max_high_order = 0.
    latent_intrinsic_sign_differences = []
    for row in audit.itertuples():
        seed, cohort = int(row.seed), row.cohort
        with np.load(output / f'latent_fields/TrpB_{seed}.npz') as z:
            a, g = z['a'], z['g']
        with np.load(output / f'cube_values/TrpB_{seed}.npz') as z:
            values = {key: z[key] for key in z.files}
        genotypes = values.pop('genotypes')
        validate_geometry(genotypes)
        indices = np.ravel_multi_index(genotypes.reshape(-1, 4).T, (20,) * 4).reshape(genotypes.shape[:2])
        for setting in settings:
            mode, ratio = setting['observation'], setting['ratio']
            base = dict(seed=seed, cohort=cohort, observation=mode, ratio=ratio)
            latent = (a + ratio*g)[indices]
            append(latent, dict(base, layer='latent'), 0., 1.)
            latent_summary, _ = descriptors(latent, 0.)
            latent_max_high_order = max(latent_max_high_order, latent_summary['pooled_ho_fraction'])
            for layer in ['intrinsic'] + ([] if mode == 'none' else ['observed']):
                cube_y = values[f'{mode}_{layer}_r{ratio:g}']
                for resolution in [0., .001, .01]:
                    append(cube_y, dict(base, layer=layer), resolution, scale)
                if layer == 'intrinsic':
                    s, _ = descriptors(cube_y, 0.)
                    latent_intrinsic_sign_differences.append(max(abs(s[k] - latent_summary[k]) for k in
                        ['sign_epistasis_fraction', 'reciprocal_sign_epistasis_fraction']))
    frame = pd.DataFrame(rows)
    archived = pd.read_csv(output / 'cube_summary.csv')
    keys = ['seed', 'cohort', 'ratio', 'observation', 'layer', 'relative_resolution']
    common = ['epistatic_fraction', 'ho_fraction', 'pooled_epistatic_fraction', 'pooled_ho_fraction',
              'sign_epistasis_fraction', 'reciprocal_sign_epistasis_fraction', 'power_availability']
    common += [f'pooled_power_{k}' for k in range(1, 5)]
    replay = archived.merge(frame, on=keys, suffixes=('_archive', '_new'), validate='one_to_one')
    assert len(replay) == len(archived) == 363
    differences = {}
    for metric in common:
        original, new = replay[f'{metric}_archive'], replay[f'{metric}_new']
        np.testing.assert_allclose(original, new, rtol=1e-11, atol=1e-13, equal_nan=True)
        differences[metric] = float(np.nanmax(np.abs(original - new)))
    assert latent_max_high_order < 1e-26
    assert max(latent_intrinsic_sign_differences) == 0.
    frame.to_csv(tables / 'architecture_seed_metrics.csv', index=False)
    pd.concat(cube_rows, ignore_index=True).to_csv(tables / 'architecture_cube_metrics.csv.gz', index=False)
    metrics = ['epistatic_fraction', 'ho_fraction', 'pooled_epistatic_fraction', 'pooled_ho_fraction',
               'sign_epistasis_fraction', 'reciprocal_sign_epistasis_fraction', 'power_availability']
    metrics += [f'{weight}_power_{k}' for weight in ['equal', 'pooled'] for k in range(1, 5)]
    summaries = []
    for group, part in frame.loc[frame.cohort.ne('empirical')].groupby(['cohort', 'observation', 'layer', 'relative_resolution']):
        cohort, mode, layer, rel = group
        for metric in metrics:
            v = part[metric].dropna().to_numpy()
            summaries.append(dict(cohort=cohort, observation=mode, layer=layer, relative_resolution=rel,
                metric=metric, n_available=len(v), n_seeds=len(part), median=float(np.median(v)) if len(v) else np.nan,
                q10=float(np.quantile(v, .1)) if len(v) else np.nan, q90=float(np.quantile(v, .9)) if len(v) else np.nan))
    pd.DataFrame(summaries).to_csv(tables / 'architecture_cohort_summary.csv', index=False)
    evidence = dict(self_checks=self_check(), source_files={r: manifest['files'][r] for r in used},
        source_manifest_sha256=_sha(output / 'manifest.json'), archived_rows_replayed=len(replay),
        maximum_absolute_replay_difference=differences, latent_maximum_order_ge3_fraction=latent_max_high_order,
        zero_threshold_latent_intrinsic_maximum_sign_difference=max(latent_intrinsic_sign_differences),
        confirmation_cubes=32, synthetic_cubes_per_seed=8,
        seeds={'fresh': int(audit.cohort.eq('fresh').sum()), 'reused': int(audit.cohort.eq('reused').sum())},
        interpretation='Synthetic seed ranges include changing assigned cube geometries. They are not confidence intervals. All 24 seeds retained. No biological coefficients inferred.')
    _json(tables / 'architecture_validation.json', evidence)
    return frame, evidence


def style():
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.labelsize': 11, 'axes.titlesize': 12, 'figure.dpi': 140,
                         'savefig.dpi': 200, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})


def save(fig, ana, stem):
    folder = Path(ana) / 'figures'
    folder.mkdir(exist_ok=True)
    for extension in ['png', 'pdf']:
        fig.savefig(folder / f'{stem}.{extension}', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return folder / f'{stem}.png'


def spectrum_figure(frame, ana):
    style()
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.2), sharey=True)
    rel = .001
    empirical = frame.loc[frame.cohort.eq('empirical') & frame.relative_resolution.eq(rel)].iloc[0]
    fresh = frame.loc[frame.cohort.eq('fresh') & frame.relative_resolution.eq(rel)]
    x = np.arange(1, 5)
    for ax, weight, title in zip(axes, ['equal', 'pooled'], ['A  Each resolvable cube has equal weight', 'B  Cubes are weighted by fitness variance']):
        cols = [f'{weight}_power_{k}' for k in x]
        ax.plot(x, empirical[cols].to_numpy(float), 'o-', color='#222222', lw=2.4, ms=6, label='Measured TrpB', zorder=5)
        for index, mode in enumerate(MODES):
            layer = 'intrinsic' if mode == 'none' else 'observed'
            part = fresh.loc[fresh.observation.eq(mode) & fresh.layer.eq(layer), cols]
            median, lo, hi = part.median().to_numpy(), part.quantile(.1).to_numpy(), part.quantile(.9).to_numpy()
            offset = (index - 1) * .065
            ax.errorbar(x + offset, median, yerr=[median-lo, hi-median], fmt='o--', color=COLORS[mode],
                        lw=1.4, markersize=4, capsize=3, label=LABELS[mode], alpha=.95)
        ax.set(xlim=(.72, 4.28), ylim=(-.02, .95), xticks=x, xlabel='Interaction order in local binary cube', title=title)
        ax.grid(axis='y', alpha=.15)
    axes[0].set_ylabel('Fraction of nonconstant fitness variance')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, -.015), ncol=2, frameon=False, fontsize=10)
    fig.suptitle('Local epistasis resembles TrpB conditionally; weighting reveals a different picture', fontsize=14, y=1.02)
    fig.subplots_adjust(bottom=.25, wspace=.13)
    return save(fig, ana, '03_local_genetic_architecture')


def stage_figure(frame, ana):
    style()
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.9))
    descriptors_to_show = [('pooled_ho_fraction', 'A  Order 3–4 variance fraction', 'Pooled higher-order variance fraction'),
                           ('sign_epistasis_fraction', 'B  Sign epistasis', 'Fraction of two-mutation faces'),
                           ('reciprocal_sign_epistasis_fraction', 'C  Reciprocal sign epistasis', 'Fraction of two-mutation faces')]
    zero = frame.loc[frame.relative_resolution.eq(0.)]
    empirical = zero.loc[zero.cohort.eq('empirical')].iloc[0]
    fresh = zero.loc[zero.cohort.eq('fresh')]
    for ax, (metric, title, ylabel) in zip(axes, descriptors_to_show):
        ax.axhline(empirical[metric], color='#222222', lw=1.5, ls='--', label='Measured TrpB')
        for stage_index, (stage, (color, marker, label)) in enumerate(STAGES.items()):
            for mode_index, mode in enumerate(MODES):
                part = fresh.loc[fresh.observation.eq(mode) & fresh.layer.eq(stage), metric].dropna().to_numpy()
                if not len(part):
                    continue
                x = mode_index + (stage_index - 1) * .18
                median, lo, hi = np.median(part), np.quantile(part, .1), np.quantile(part, .9)
                ax.plot(np.full(len(part), x) + np.linspace(-.035, .035, len(part)), part,
                        linestyle='', marker=marker, ms=3.5, color=color, alpha=.25)
                ax.errorbar(x, median, yerr=[[median-lo], [hi-median]], fmt=marker, color=color,
                            lw=1.8, capsize=3, ms=6, label=label if mode_index == 0 or (stage == 'observed' and mode_index == 1) else None)
        ax.set(title=title, ylabel=ylabel, xticks=np.arange(3), xticklabels=[SHORT[m] for m in MODES], xlim=(-.42, 2.42))
        ax.tick_params(axis='x', labelsize=8.3)
        ax.set_ylim(-.035, .83 if metric == 'sign_epistasis_fraction' else .62)
        ax.grid(axis='y', alpha=.15)
    handles, labels = axes[0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(unique.values(), unique.keys(), loc='lower center', bbox_to_anchor=(.5, -.045), ncol=4, frameon=False, fontsize=9.5)
    fig.suptitle('The sigmoid adds higher-order power; observation noise also changes mutation signs', fontsize=14, y=1.03)
    fig.text(.5, .04, 'Zero threshold. Faint points: 8 fresh seeds; large markers and bars: median and 10–90% seed range.', ha='center', fontsize=9)
    fig.subplots_adjust(wspace=.32, bottom=.25)
    return save(fig, ana, '04_latent_readout_observation')
