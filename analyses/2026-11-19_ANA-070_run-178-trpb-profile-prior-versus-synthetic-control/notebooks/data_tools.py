"""Deterministic comparison tables for the two frozen TrpB runs; no refitting."""
from pathlib import Path
import hashlib
import json
import platform
import importlib.metadata
import numpy as np
import pandas as pd

ANA = Path(__file__).resolve().parents[1]
ROOT = ANA.parents[1]
TABLES = ANA / 'tables'
FIGURES = ANA / 'figures'
BUNDLES = {name: ROOT / 'data/processed/aggregate-thermodynamic-gnk-calibration' / run
           for name, run in [('Synthetic', 'RUN-175'), ('Profile', 'RUN-178')]}
KEYS = ['protein', 'seed', 'cohort', 'ratio', 'observation', 'layer', 'relative_resolution']
GROUPS = ['prior', 'cohort', 'ratio', 'observation', 'layer', 'relative_resolution']
METRICS = {
    'uniform_steps': 'Uphill walk length',
    'uniform_hit_top': 'Top-1% success',
    'greedy_steps': 'Greedy walk length',
    'greedy_hit_top': 'Greedy top-1% success',
    'uniform_hit_best': 'Exact-best success',
    'accessible_top': 'Top-1% accessible',
    'epistatic_fraction': 'Equal-cube epistatic power fraction',
    'ho_fraction': 'Equal-cube high-order power fraction',
    'pooled_epistatic_fraction': 'Pooled epistatic power fraction',
    'pooled_ho_fraction': 'Pooled high-order power fraction',
    'sign_epistasis_fraction': 'Sign epistasis',
    'reciprocal_sign_epistasis_fraction': 'Reciprocal sign epistasis',
    'slope_median': 'Median diminishing-returns slope',
    'slope': 'Mean diminishing-returns slope',
    'sd_fitness': 'Confirmation local fitness SD',
    'near_floor_fraction': 'Confirmation near-floor fraction',
    'power_availability': 'Resolvable WHT fraction',
    'slope_availability': 'Estimable slope fraction',
    'neutral_edge_fraction': 'Neutral-edge fraction',
    'above_floor_largest_fraction': 'Largest above-floor component',
    'all_largest_fraction': 'Largest neutral component',
}
FOCUS = ['uniform_steps', 'uniform_hit_top', 'epistatic_fraction', 'ho_fraction',
         'pooled_epistatic_fraction', 'pooled_ho_fraction', 'sign_epistasis_fraction',
         'slope_median', 'sd_fitness']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def load():
    TABLES.mkdir(exist_ok=True); FIGURES.mkdir(exist_ok=True)
    frames, fit_frames, latent_frames = [], [], []
    for prior, bundle in BUNDLES.items():
        output = bundle / 'outputs'
        graph = pd.read_csv(output / 'graph_summary.csv').query("protein == 'TrpB'")
        cubes = pd.read_csv(output / 'cube_summary.csv').query("protein == 'TrpB'")
        merged = graph.merge(cubes, on=KEYS, validate='one_to_one')
        assert len(merged) == 2163
        merged['prior'] = prior
        frames.append(merged)
        fit = pd.json_normalize(json.loads((output / 'selected_parameters.json').read_text()))
        fit.columns = [column.replace('.', '_') for column in fit.columns]
        fit = fit.loc[fit.protein.eq('TrpB')].copy()
        assert len(fit) == 18
        fit['prior'] = prior; fit_frames.append(fit)
        latent = pd.read_csv(output / 'latent_graph_summary.csv').query("protein == 'TrpB'").copy()
        latent['prior'] = prior; latent_frames.append(latent)
    data = pd.concat(frames, ignore_index=True)
    fits = pd.concat(fit_frames, ignore_index=True)
    latent = pd.concat(latent_frames, ignore_index=True)
    assert not data.duplicated(['prior'] + KEYS).any()
    return data, fits, latent


def empirical(data, metric, epsilon=.001):
    selected = data.loc[data.prior.eq('Profile') & data.cohort.eq('empirical')
                        & data.relative_resolution.eq(epsilon), metric]
    assert len(selected) == 1
    return float(selected.iloc[0])


def observed_rows(data, mode='shared_stop_sensitivity', cohort='fresh', epsilon=.001):
    layer = 'intrinsic' if mode == 'none' else 'observed'
    return data.loc[data.observation.eq(mode) & data.layer.eq(layer)
                    & data.cohort.eq(cohort) & data.relative_resolution.eq(epsilon)].copy()


def make_tables(data, fits, latent):
    # Preserve all settings, both cohorts, all three resolutions, and the intrinsic
    # maps associated with noisy calibrations, even when not highlighted in plots.
    data.to_csv(TABLES / 'per_seed_results.csv.gz', index=False,
                compression={'method': 'gzip', 'mtime': 0})
    fits.to_csv(TABLES / 'selected_parameters.csv', index=False)
    latent.to_csv(TABLES / 'latent_graph_results.csv', index=False)
    aggregate = []
    for key, group in data.groupby(GROUPS, sort=True, dropna=False):
        base = dict(zip(GROUPS, key))
        for metric, label in METRICS.items():
            values = group[metric].dropna()
            aggregate.append(dict(**base, metric=metric, label=label, n_total=len(group),
                n_available=len(values), median=values.median(),
                q10=values.quantile(.1), q90=values.quantile(.9)))
    aggregate = pd.DataFrame(aggregate)
    aggregate.to_csv(TABLES / 'aggregate_estimates.csv', index=False)

    model = data.loc[~data.cohort.eq('empirical')]
    joined = model.loc[model.prior.eq('Synthetic'), KEYS + list(METRICS)].merge(
        model.loc[model.prior.eq('Profile'), KEYS + list(METRICS)],
        on=KEYS, suffixes=('_synthetic', '_profile'), validate='one_to_one')
    records = []
    for metric in METRICS:
        frame = joined[KEYS].copy()
        frame['metric'] = metric
        frame['synthetic'] = joined[metric + '_synthetic']
        frame['profile'] = joined[metric + '_profile']
        targets = {epsilon: empirical(data, metric, epsilon) for epsilon in [0., .001, .01]}
        frame['empirical'] = frame.relative_resolution.map(targets)
        frame['profile_minus_synthetic'] = frame.profile - frame.synthetic
        frame['synthetic_absolute_error'] = abs(frame.synthetic - frame.empirical)
        frame['profile_absolute_error'] = abs(frame.profile - frame.empirical)
        frame['absolute_error_change'] = frame.profile_absolute_error - frame.synthetic_absolute_error
        frame['paired_available'] = frame[['synthetic', 'profile', 'empirical']].notna().all(axis=1)
        frame['closer'] = frame.paired_available & frame.absolute_error_change.lt(-1e-10)
        frame['farther'] = frame.paired_available & frame.absolute_error_change.gt(1e-10)
        frame['roundoff_tie'] = frame.paired_available & frame.absolute_error_change.abs().le(1e-10)
        records.append(frame)
    paired = pd.concat(records, ignore_index=True)
    paired.to_csv(TABLES / 'paired_seed_comparisons.csv.gz', index=False,
                  compression={'method': 'gzip', 'mtime': 0})
    groups = ['cohort', 'ratio', 'observation', 'layer', 'relative_resolution', 'metric']
    summaries = []
    for key, group in paired.groupby(groups, sort=True):
        available = group.loc[group.paired_available]
        summaries.append(dict(**dict(zip(groups, key)), n_pairs=len(available),
            closer=int(available.closer.sum()), farther=int(available.farther.sum()),
            roundoff_ties=int(available.roundoff_tie.sum()),
            median_value_change=available.profile_minus_synthetic.median(),
            median_absolute_error_change=available.absolute_error_change.median(),
            error_change_q10=available.absolute_error_change.quantile(.1),
            error_change_q90=available.absolute_error_change.quantile(.9)))
    paired_summary = pd.DataFrame(summaries)
    paired_summary.to_csv(TABLES / 'paired_summary.csv', index=False)
    chosen = aggregate.loc[aggregate.relative_resolution.eq(.001)
        & aggregate.ratio.eq(.5) & aggregate.observation.eq('shared_stop_sensitivity')
        & aggregate.layer.eq('observed') & aggregate.metric.isin(FOCUS)].copy()
    chosen['empirical'] = chosen.metric.map(lambda metric: empirical(data, metric))
    match = paired_summary.loc[paired_summary.relative_resolution.eq(.001)
        & paired_summary.ratio.eq(.5) & paired_summary.observation.eq('shared_stop_sensitivity')
        & paired_summary.layer.eq('observed')]
    chosen = chosen.merge(match[['cohort', 'metric', 'n_pairs', 'closer', 'farther',
        'roundoff_ties', 'median_absolute_error_change']], on=['cohort', 'metric'], validate='many_to_one')
    chosen.to_csv(TABLES / 'inherited_regime_comparison.csv', index=False)
    versions = {name: importlib.metadata.version(name) for name in
        ['numpy', 'pandas', 'scipy', 'pyarrow', 'PyYAML', 'matplotlib', 'nbformat', 'nbclient', 'ipykernel']}
    dump(TABLES / 'analysis_environment.json', dict(python=platform.python_version(),
         platform=platform.platform(), packages=versions))
    dump(TABLES / 'table_manifest.json', {
        'run_inputs': {prior: {name: sha(bundle / 'outputs' / name) for name in
                      ['manifest.json', 'graph_summary.csv', 'cube_summary.csv', 'selected_parameters.json']}
                       for prior, bundle in BUNDLES.items()},
        'primary_display': {'cohort': 'fresh', 'seeds': 8, 'relative_resolution': .001,
                            'inherited_ratio': .5, 'inherited_observation': 'shared_stop_sensitivity'},
        'pairing_keys': KEYS, 'paired_rows': len(paired), 'aggregate_rows': len(aggregate),
        'range_interpretation': '10th–90th percentiles of generator-seed outcomes, not confidence intervals',
        'closer_definition': 'Absolute error to the empirical descriptor decreases by more than 1e-10; this numerical tie tolerance is not a scientific acceptance threshold.',
        'no_aggregate_weighted_score': True, 'scientific_acceptance_gate': False,
        'no_new_fit_or_simulation': True,
    })
    return aggregate, paired_summary, chosen
