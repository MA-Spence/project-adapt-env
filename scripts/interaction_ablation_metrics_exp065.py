"""Replicate-aware aggregate diagnostics for the EXP-065 r=0 ablation.

No empirical coefficient identities are fitted. Cross-power retains shared
measurement error; the stop correction is explicitly a sensitivity analysis.
"""
from itertools import combinations
import numpy as np
import pandas as pd

ORDER = np.array([i.bit_count() for i in range(16)])
H = np.array([[(-1.) ** ((i & j).bit_count()) for j in range(16)] for i in range(16)])
FACES = np.array([[s, s | (1 << i), s | (1 << j), s | (1 << i) | (1 << j)]
                  for i, j in combinations(range(4), 2)
                  for s in range(16) if not s & ((1 << i) | (1 << j))])


def face_switches(y, epsilon):
    """Return the two directed sign reversals per square, with strict resolution."""
    y = np.asarray(y, float)
    z = y[:, FACES]
    edges = np.stack([z[..., 1] - z[..., 0], z[..., 3] - z[..., 2],
                      z[..., 2] - z[..., 0], z[..., 3] - z[..., 1]], axis=-1)
    signs = (edges > epsilon).astype(int) - (edges < -epsilon).astype(int)
    return signs, np.stack([signs[..., 0] * signs[..., 1] == -1,
                           signs[..., 2] * signs[..., 3] == -1], axis=-1)


def replicate_metrics(rep1, rep2, shared_covariance, cutoffs, resolutions):
    """Return per-cube absolute powers and per-face reproducibility descriptors.

    A reproducible sign reversal must involve the SAME mutation with the SAME
    two edge signs in both replicates. Reciprocal requires all four edge signs
    to agree. Activity selection is observation-conditioned, not noise-free.
    """
    a, b = np.asarray(rep1, float), np.asarray(rep2, float)
    if a.ndim != 2 or a.shape != b.shape or a.shape[1] != 16 or not np.isfinite([a, b]).all():
        raise ValueError('Expected finite paired four-site cubes')
    q = np.broadcast_to(shared_covariance, a.shape)
    if not np.isfinite(q).all():
        raise ValueError('Nonfinite shared covariance')
    c1, c2 = a @ H.T / 16, b @ H.T / 16
    mean = (a + b) / 2
    cm = (c1 + c2) / 2
    active = (a > cutoffs[0]) & (b > cutoffs[1])
    active_count = active.sum(axis=1)
    cube = pd.DataFrame(dict(cube=np.arange(len(a)), active_nodes=active_count,
        stratum=np.where(active_count == 0, 'no_active', np.where(active_count == 16, 'all_active', 'some_active')),
        mean_fitness=mean.mean(axis=1), sd_fitness=mean.std(axis=1),
        near_floor_fraction=(mean <= np.mean(cutoffs)).mean(axis=1),
        replicate_difference_variance=(a-b).var(axis=1)))
    for k in range(1, 5):
        mask = ORDER == k
        ordinary = np.sum(cm[:, mask] ** 2, axis=1)
        cross = np.sum(c1[:, mask] * c2[:, mask], axis=1)
        correction = mask.sum() * q.sum(axis=1) / 256
        cube[f'mean_power_{k}'] = ordinary
        cube[f'cross_power_{k}'] = cross
        cube[f'stop_adjusted_power_{k}'] = cross - correction
    np.testing.assert_allclose(cube[[f'mean_power_{k}' for k in range(1, 5)]].sum(axis=1),
                               mean.var(axis=1), rtol=1e-11, atol=1e-24)
    face_active = active[:, FACES].sum(axis=-1)
    faces = []
    for epsilon in resolutions:
        s1, w1 = face_switches(a, epsilon)
        s2, w2 = face_switches(b, epsilon)
        _, wm = face_switches(mean, epsilon)
        agree_first = np.all(s1[..., :2] == s2[..., :2], axis=-1)
        agree_second = np.all(s1[..., 2:] == s2[..., 2:], axis=-1)
        reproducible = w1 & w2 & np.stack([agree_first, agree_second], axis=-1)
        faces.append(pd.DataFrame(dict(cube=np.repeat(np.arange(len(a)), 24),
            face=np.tile(np.arange(24), len(a)), epsilon=float(epsilon),
            active_vertices=face_active.ravel(),
            stratum=np.where(face_active == 0, 'no_active', np.where(face_active == 4, 'all_active', 'some_active')).ravel(),
            mean_sign=wm.any(axis=-1).ravel(), mean_reciprocal=wm.all(axis=-1).ravel(),
            replicate1_sign=w1.any(axis=-1).ravel(), replicate2_sign=w2.any(axis=-1).ravel(),
            reproducible_sign=reproducible.any(axis=-1).ravel(),
            reproducible_reciprocal=reproducible.all(axis=-1).ravel(),
            all_edge_signs_agree=np.all(s1 == s2, axis=-1).ravel())))
    return cube, pd.concat(faces, ignore_index=True)


def summarize_replicates(cubes, faces):
    """Keep empty strata and signed powers; never silently select available ratios."""
    rows = []
    for stratum in ['all', 'no_active', 'some_active', 'all_active']:
        c = cubes if stratum == 'all' else cubes[cubes.stratum == stratum]
        row = dict(stratum=stratum, n_cubes=len(c))
        for method in ['mean', 'cross', 'stop_adjusted']:
            powers = c[[f'{method}_power_{k}' for k in range(1, 5)]].to_numpy()
            totals = powers.sum(axis=1)
            pooled = powers.sum(axis=0)
            for k in range(1, 5):
                row[f'{method}_power_{k}'] = float(powers[:, k-1].mean()) if len(c) else np.nan
            row[f'{method}_positive_total_cubes'] = int((totals > 1e-18).sum())
            row[f'{method}_pooled_ho_fraction'] = float(pooled[2:].sum()/pooled.sum()) if len(c) and pooled.sum() > len(c)*1e-18 else np.nan
            row[f'{method}_pooled_epistatic_fraction'] = float(pooled[1:].sum()/pooled.sum()) if len(c) and pooled.sum() > len(c)*1e-18 else np.nan
        for epsilon, f in faces.groupby('epsilon', sort=False):
            f = f if stratum == 'all' else f[f.stratum == stratum]
            item = dict(row, epsilon=float(epsilon), n_faces=len(f))
            for name in ['mean_sign', 'mean_reciprocal', 'replicate1_sign', 'replicate2_sign',
                         'reproducible_sign', 'reproducible_reciprocal', 'all_edge_signs_agree']:
                item[name] = float(f[name].mean()) if len(f) else np.nan
            rows.append(item)
    return pd.DataFrame(rows)


def static_graph(values, valid, epsilon, q=20, dim=4):
    """Exact local traps and edge moments; no walk or optimiser is executed."""
    values, valid = np.asarray(values, float), np.asarray(valid, bool)
    if values.shape != valid.shape or values.size != q**dim or not valid.any() or not np.isfinite(values[valid]).all():
        raise ValueError('Invalid common-mask landscape')
    f = np.where(valid, values, -np.inf).reshape((q,)*dim)
    mask = valid.reshape((q,)*dim)
    centered = np.where(valid, values-values[valid].mean(), 0.).reshape((q,)*dim)
    best_gain = np.full(f.shape, -np.inf)
    edges = 0
    squared_difference = endpoint_sum = endpoint_square = 0.
    for axis in range(dim):
        # Including self guarantees best gain >=0, without creating a beneficial self edge.
        gain = np.zeros_like(f)
        np.subtract(np.max(f, axis=axis, keepdims=True), f, out=gain, where=mask)
        np.maximum(best_gain, gain, out=best_gain)
        n, s, t = mask.sum(axis=axis), centered.sum(axis=axis), (centered**2).sum(axis=axis)
        edges += int(np.sum(n*(n-1)//2))
        squared_difference += float(np.sum(n*t-s*s))
        endpoint_sum += float(np.sum((n-1)*s))
        endpoint_square += float(np.sum((n-1)*t))
    if edges == 0:
        raise ValueError('Common-mask landscape has no edges')
    variance = endpoint_square/(2*edges)-(endpoint_sum/(2*edges))**2
    correlation = 1-squared_difference/edges/(2*variance) if variance > 0 else np.nan
    return dict(nodes=int(valid.sum()), edges=edges, epsilon=float(epsilon),
        no_beneficial_fraction=float((best_gain[mask] <= epsilon).mean()),
        neighbor_correlation=float(correlation),
        edge_rms_native=float(np.sqrt(max(0., squared_difference/edges))),
        mean_fitness=float(values[valid].mean()), sd_fitness=float(values[valid].std()))
