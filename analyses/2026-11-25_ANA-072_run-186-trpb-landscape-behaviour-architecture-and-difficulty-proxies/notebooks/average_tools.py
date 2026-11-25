"""Read-only reconstruction and descriptive summaries for ANA-072 / RUN-186.

The archived latent fields already use each seed's computational reference.
Do not re-center them. The observation function is imported from frozen source;
the empirical donor pools are restored from the archived assignments, with no
parameter fitting, new donor selection, or new random seed.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import math

import numpy as np
import pandas as pd
from scipy.special import expit

ANA_NAME = '2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies'
MODE_LABELS = {'none': 'No observation error',
               'independent_difference': 'Independent replicate error',
               'shared_stop_sensitivity': 'Shared-stop sensitivity'}
MODE_COLORS = {'none': '#2878a5', 'independent_difference': '#df8f27',
               'shared_stop_sensitivity': '#9b4d96'}
MODES = list(MODE_LABELS)


class Context:
    def __init__(self, project):
        self.project = Path(project)
        self.analysis = self.project / 'analyses' / ANA_NAME
        self.run = self.project / 'data/processed/aggregate-thermodynamic-gnk-calibration/RUN-186'
        self.output = self.run / 'outputs'
        self.source = self.run / 'frozen_runtime/source'
        self.settings = json.loads((self.output / 'selected_parameters.json').read_text())
        self.audit = pd.read_csv(self.output / 'seed_audit.csv')
        self.cubes = pd.read_csv(self.output / 'cube_summary.csv')
        with np.load(self.output / 'TrpB_empirical.npz') as z:
            self.empirical = z['fitness'].copy()
            self.valid = z['valid'].copy()
            self.genotypes = z['genotypes'].copy()
        with np.load(self.output / 'TrpB_empirical_cubes.npz') as z:
            self.empirical_cubes = z['fitness'].copy()
        self.floor = float(self.settings[0]['near_floor_raw'])
        self.scale = float(self.settings[0]['calibration_scale'])
        self.epsilon = .001 * self.scale
        self.reconstruction_checks = []
        self.noise = restore_noise(self.output)
        script = self.source / 'scripts/replicate_observation_exp065.py'
        spec = importlib.util.spec_from_file_location('_ana072_frozen_observation', script)
        self.observation_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.observation_module)
        self.observation_source_sha256 = hashlib.sha256(script.read_bytes()).hexdigest()


def restore_noise(output):
    """Restore fixed pools in the same row order used by the frozen sampler."""
    info = json.loads((output / 'noise/noise_model.json').read_text())
    table = pd.read_csv(output / 'noise/donor_assignments.csv.gz')
    train = table[(table.stop_count == 0) & (table.split == 'training')]
    donors, depths, common = {}, {}, []
    for entry in info['cell_audit']:
        f, d = entry['fitness_bin'], entry['depth_bin']
        fs = train[train.fitness_bin == f]
        cell = fs[fs.depth_bin == d]
        selected = {'cell': cell, 'fitness': fs, 'all_training': train}[entry['fallback']]
        values = selected.difference.to_numpy()
        donors[f, d] = values - values.mean()
        depth_pool = fs if entry['depth_fallback'] == 'fitness' else train
        depths[f] = depth_pool.input_count.to_numpy()
        assert len(values) == entry['donor_n']
        assert len(depths[f]) == entry['depth_pool_n']
    q = np.asarray(info['shared_covariance_raw'])
    for d in range(4):
        stop = table[(table.stop_count > 0) & (table.split == 'training') & (table.depth_bin == d)]
        means = (stop.replicate_1.to_numpy() + stop.replicate_2.to_numpy()) / 2
        means -= means.mean()
        var = means.var()
        common.append(means * math.sqrt(max(0., q[d]) / var) if var > 0 else means * 0)
    return dict(donors=donors, depths=depths, common=common, q=q,
                depth_breaks=np.asarray(info['depth_breaks']),
                fitness_breaks=np.asarray(info['fitness_breaks']))


def iter_maps(ctx, cohorts=None):
    """Yield (metadata, 160000 native-unit fitness values), verifying saved cubes.

    `cohorts` can restrict seeds; None includes all 24. Each call deterministically
    reconstructs the same static replicate-mean observation maps used in RUN-186.
    The noise here is frozen by genotype, not resampled on repeated queries.
    """
    for row in ctx.audit.itertuples(index=False):
        if cohorts is not None and row.cohort not in cohorts:
            continue
        with np.load(ctx.output / 'latent_fields' / f'TrpB_{row.seed}.npz') as z:
            a, g = z['a'], z['g']
        with np.load(ctx.output / 'cube_values' / f'TrpB_{row.seed}.npz') as z:
            cube_genotypes = z['genotypes']
            indices = np.ravel_multi_index(cube_genotypes.reshape(-1, 4).T, (20,) * 4)
            saved = {k: z[k].copy() for k in z.files if k != 'genotypes'}
        uniforms = np.random.default_rng(int(row.noise_seed)).random((len(a), 6))
        for setting in ctx.settings:
            ratio, mode = setting['ratio'], setting['observation']
            p = setting['parameters']
            raw = setting['calibration_origin'] + setting['calibration_scale'] * (
                p[0] + np.exp(p[1]) * expit(p[2] + np.exp(p[3]) * (a + ratio * g)))
            maps = [('intrinsic', raw)]
            if mode != 'none':
                r1, r2, _, _ = ctx.observation_module.observe(
                    raw, ctx.noise, uniforms, shared=mode == 'shared_stop_sensitivity')
                maps.append(('observed', (r1 + r2) / 2))
            for layer, values in maps:
                key = f'{mode}_{layer}_r{ratio:g}'
                rebuilt = values[indices].reshape(cube_genotypes.shape[:2])
                maxdiff = float(np.max(np.abs(rebuilt - saved[key])))
                np.testing.assert_allclose(rebuilt, saved[key], rtol=1e-11, atol=1e-12)
                meta = dict(seed=int(row.seed), cohort=row.cohort, observation=mode,
                            ratio=float(ratio), layer=layer)
                ctx.reconstruction_checks.append(dict(**meta, nodes=len(values),
                    saved_cube_values=rebuilt.size, maximum_absolute_difference=maxdiff))
                yield meta, values


def fixed_direction_effects(values):
    y = np.asarray(values)
    low = np.array([[i for i in range(16) if not i & (1 << j)] for j in range(4)])
    high = low | (1 << np.arange(4))[:, None]
    return (y[:, high] - y[:, low]).ravel()


def quantile_summary(frame, group_cols, value='value'):
    """Equal weight per generator seed; ranges are variability, not confidence."""
    records = []
    for keys, group in frame.groupby(group_cols, sort=False, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        v = group[value].dropna().to_numpy()
        records.append(dict(zip(group_cols, keys)) | dict(n=len(v),
            q10=float(np.quantile(v, .1)), median=float(np.median(v)),
            q90=float(np.quantile(v, .9))))
    return pd.DataFrame(records)


def save_figure(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    for suffix in ('png', 'pdf'):
        fig.savefig(directory / f'{name}.{suffix}', dpi=180, bbox_inches='tight')
