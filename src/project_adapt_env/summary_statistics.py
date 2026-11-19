"""Versioned aggregate statistics on measured complete binary subspaces.

Mutation direction and the background measure are explicit. These are statistics
of the supplied subspace, never estimates of an unobserved whole-protein WHT.
"""
from __future__ import annotations

import math
import numpy as np

METRIC_VERSION = "aggregate_v2"


def cube_values(values, dim):
    array = np.asarray(values, dtype=float)
    if int(dim) != dim or dim < 1 or array.shape != (1 << int(dim),) or not np.isfinite(array).all():
        raise ValueError("a complete finite binary cube with dimension >= 1 is required")
    return array


def plateau_components(values, dim, resolution=1e-12):
    arr = cube_values(values, dim)
    if not np.isfinite(resolution) or resolution < 0:
        raise ValueError("resolution must be finite and nonnegative")
    parent = np.arange(len(arr))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return int(i)
    neutral = 0
    for i in range(len(arr)):
        for bit in range(dim):
            j = i ^ (1 << bit)
            if i < j and abs(arr[i] - arr[j]) <= resolution:
                parent[root(j)] = root(i)
                neutral += 1
    labels = np.array([root(i) for i in range(len(arr))])
    groups = [np.flatnonzero(labels == label) for label in np.unique(labels)]
    return groups, labels, neutral


def extrema_metrics(values, *, dim, resolution=1e-12):
    arr = cube_values(values, dim)
    groups, labels, neutral = plateau_components(arr, dim, resolution)
    maxima, minima, weak_maxima = [], [], []
    for i in range(len(arr)):
        delta = arr[i] - arr[[i ^ (1 << bit) for bit in range(dim)]]
        if np.all(delta > resolution):
            maxima.append(i)
        if np.all(delta < -resolution):
            minima.append(i)
        if np.all(delta >= -resolution):
            weak_maxima.append(i)
    terminal = 0
    for group in groups:
        improving_exit = any(arr[i ^ (1 << bit)] > arr[i] + resolution
                             for i in group for bit in range(dim)
                             if labels[i ^ (1 << bit)] != labels[i])
        terminal += int(not improving_exit)
    return {
        "n_local_maxima": float(len(maxima)), "n_local_minima": float(len(minima)),
        "n_weak_local_maxima": float(len(weak_maxima)),
        "local_maxima": ";".join(map(str, maxima)), "local_minima": ";".join(map(str, minima)),
        "local_maxima_density": len(maxima) / len(arr),
        "neutral_edge_fraction": neutral / (dim * len(arr) / 2),
        "n_plateau_components": float(len(groups)), "n_terminal_plateaus": float(terminal),
        "largest_plateau_fraction": max(map(len, groups)) / len(arr),
        "flat_or_unresolved": bool(np.ptp(arr) <= resolution),
        "resolution": float(resolution),
    }


def mutation_summary(values, *, dim, resolution=1e-12, reference_node=0, replicates=None):
    """DFE at reference; same 0->1 mutations over all other-bit backgrounds.

    Without independent replicates, diminishing-returns slopes are descriptive
    and may contain shared-error bias. With >=2 independent replicates, average
    ordered cross-replicate covariances (including the background denominator)
    to remove independent measurement-error covariance. Negative/noisy background
    variance estimates are unavailable, not forced into a valid slope.
    Beneficial/deleterious strata classify each fixed mutation by its mean effect
    across backgrounds, never by selecting its effect separately at each background.
    """
    arr = cube_values(values, dim)
    if not np.isfinite(resolution) or resolution < 0:
        raise ValueError("resolution must be finite and nonnegative")
    if reference_node < 0 or reference_node >= len(arr):
        raise ValueError("invalid reference node")
    dfe = arr[[reference_node ^ (1 << bit) for bit in range(dim)]] - arr[reference_node]
    reps = None if replicates is None else np.asarray(replicates, dtype=float)
    if reps is not None and (reps.ndim != 2 or reps.shape[0] < 2 or reps.shape[1] != len(arr) or not np.isfinite(reps).all()):
        raise ValueError("replicates must contain >=2 independent finite measurements of the complete cube")
    slopes, r2s, beneficial, deleterious = [], [], [], []
    for bit in range(dim):
        low = np.array([n for n in range(len(arr)) if not n & (1 << bit)])
        high = low | (1 << bit)
        x, y = arr[low], arr[high] - arr[low]
        xc, yc = x - x.mean(), y - y.mean()
        variance = float(np.mean(xc * xc)); covariance = float(np.mean(xc * yc))
        if reps is not None:
            xs = reps[:, low]; ys = reps[:, high] - reps[:, low]
            xs = xs - xs.mean(axis=1, keepdims=True)
            ys = ys - ys.mean(axis=1, keepdims=True)
            pairs = [(a, b) for a in range(len(reps)) for b in range(len(reps)) if a != b]
            variance = float(np.mean([np.mean(xs[a] * xs[b]) for a, b in pairs]))
            covariance = float(np.mean([np.mean(xs[a] * ys[b]) for a, b in pairs]))
        if variance > max(resolution**2, 1e-18):
            slope = covariance / variance
            slopes.append(slope)
            if y.mean() > resolution:
                beneficial.append(slope)
            elif y.mean() < -resolution:
                deleterious.append(slope)
            if np.var(y) > 1e-18 and reps is None:
                r2s.append(float(np.corrcoef(x, y)[0, 1] ** 2))
    return {
        "reference_dfe_mean": float(dfe.mean()), "reference_dfe_sd": float(dfe.std()),
        "reference_beneficial_fraction": float(np.mean(dfe > resolution)),
        "reference_neutral_fraction": float(np.mean(abs(dfe) <= resolution)),
        "reference_deleterious_fraction": float(np.mean(dfe < -resolution)),
        "diminishing_returns_slope": float(np.mean(slopes)) if slopes else math.nan,
        "diminishing_returns_slope_sd": float(np.std(slopes)) if slopes else math.nan,
        "diminishing_returns_n_mutations": len(slopes),
        "beneficial_diminishing_returns_slope": float(np.mean(beneficial)) if beneficial else math.nan,
        "beneficial_slope_n_mutations": len(beneficial),
        "deleterious_background_slope": float(np.mean(deleterious)) if deleterious else math.nan,
        "deleterious_slope_n_mutations": len(deleterious),
        "diminishing_returns_mean_r2": float(np.mean(r2s)) if r2s else math.nan,
        "diminishing_returns_error_corrected": reps is not None,
    }


def accessible_fraction(values, *, dim, resolution=1e-12):
    """Fraction of strictly increasing shortest paths between opposite corners.

    Orient from lower-fitness corner to higher; no claim that an endpoint is
    a global peak. Dynamic programming avoids factorial path enumeration.
    """
    arr = cube_values(values, dim)
    if arr[-1] < arr[0]:
        arr = arr[::-1]
    counts = [0] * len(arr); counts[0] = 1
    for node in range(1, len(arr)):
        for bit in range(dim):
            if node & (1 << bit):
                prev = node ^ (1 << bit)
                if arr[node] > arr[prev] + resolution:
                    counts[node] += counts[prev]
    return counts[-1] / math.factorial(dim)


def population_edge_effects(states, values, *, background_weights=None, max_edges=None):
    """Oriented Hamming-1 effects with explicit source probabilities.

    Uniform backgrounds give a symmetric distribution on a complete regular
    graph. A reference DFE uses a point mass at that reference. Return effects
    and their weights; never silently label the symmetric graph DFE as WT DFE.
    """
    if not states:
        return np.array([]), np.array([])
    if len(states) != len(values) or not np.isfinite(values).all():
        raise ValueError("states and finite values must align")
    if len(set(states)) != len(states):
        raise ValueError("duplicate genotype states must be aggregated first")
    score = dict(zip(states, map(float, values)))
    mass = np.ones(len(states)) if background_weights is None else np.asarray(background_weights, dtype=float)
    if mass.shape != (len(states),) or not np.isfinite(mass).all() or np.any(mass < 0) or mass.sum() <= 0:
        raise ValueError("invalid background probabilities")
    if not states:
        return np.array([]), np.array([])
    alleles = [sorted({s[i] for s in states}) for i in range(len(states[0]))]
    effects, weights = [], []
    for state, p in zip(states, mass / mass.sum()):
        if p == 0:
            continue
        neighbors = []
        for site, alternatives in enumerate(alleles):
            for allele in alternatives:
                other = state[:site] + (allele,) + state[site + 1:]
                if allele != state[site] and other in score:
                    neighbors.append(other)
        for other in neighbors:
            effects.append(score[other] - score[state]); weights.append(p / len(neighbors))
        if max_edges is not None and len(effects) >= max_edges:
            break
    w = np.asarray(weights)
    return np.asarray(effects), w / w.sum() if w.size else w
