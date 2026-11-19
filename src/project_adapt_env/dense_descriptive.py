"""Exact descriptive summaries on complete binary genotype subspaces.

The functions here intentionally do not run ABC or any stochastic simulation.
They evaluate fixed empirical/model values on complete Hamming subcubes and
return deterministic summary failures that can later define an ABC target.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .summary_statistics import METRIC_VERSION, accessible_fraction, cube_values, extrema_metrics, plateau_components


MutationFeature = tuple[int, str]
MutationMask = frozenset[MutationFeature]


@dataclass(frozen=True)
class CompleteGraphlet:
    graphlet_id: str
    dimension: int
    features: tuple[MutationFeature, ...]

    @property
    def mutation_features(self) -> str:
        return ";".join(f"{pos}:{aa}" for pos, aa in self.features)


def mutation_set(wildtype: str, sequence: str) -> MutationMask:
    return frozenset(
        (idx, aa)
        for idx, (wt_aa, aa) in enumerate(zip(str(wildtype), str(sequence)))
        if aa != wt_aa
    )


def node_masks(features: Sequence[MutationFeature]) -> list[MutationMask]:
    """Return masks in binary node order for a hypercube over ``features``."""

    feature_tuple = tuple(features)
    masks: list[MutationMask] = []
    for node in range(1 << len(feature_tuple)):
        masks.append(
            frozenset(feature for bit, feature in enumerate(feature_tuple) if (node >> bit) & 1)
        )
    return masks


def build_score_map(
    frame: pd.DataFrame,
    *,
    wildtype: str,
    score_columns: Sequence[str],
    wt_scores: Mapping[str, float] | None = None,
    sequence_column: str = "mutated_sequence",
) -> dict[MutationMask, dict[str, float]]:
    grouped: dict[MutationMask, dict[str, list[float]]] = {}
    for row in frame.itertuples(index=False):
        sequence = str(getattr(row, sequence_column))
        if len(sequence) != len(wildtype):
            continue
        key = mutation_set(wildtype, sequence)
        by_col = grouped.setdefault(key, {})
        for col in score_columns:
            if not hasattr(row, col):
                continue
            try:
                value = float(getattr(row, col))
            except (TypeError, ValueError):
                continue
            if np.isfinite(value):
                by_col.setdefault(str(col), []).append(value)
    empty = frozenset()
    if wt_scores:
        by_col = grouped.setdefault(empty, {})
        for col, value in wt_scores.items():
            if col not in score_columns:
                continue
            try:
                numeric = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(numeric):
                by_col.setdefault(str(col), []).append(numeric)

    out: dict[MutationMask, dict[str, float]] = {}
    for key, values_by_col in grouped.items():
        out[key] = {
            col: float(np.mean(values))
            for col, values in values_by_col.items()
            if values
        }
    return out


def feature_frequency(score_map: Mapping[MutationMask, Mapping[str, float]]) -> Counter[MutationFeature]:
    counts: Counter[MutationFeature] = Counter()
    for mask in score_map:
        for feature in mask:
            counts[feature] += 1
    return counts


def enumerate_complete_graphlets(
    score_map: Mapping[MutationMask, Mapping[str, float]],
    *,
    assay_id: str,
    score_columns: Sequence[str],
    dimensions: Sequence[int],
    max_graphlets_per_dimension: int,
    top_feature_pool: int,
    max_combinations_per_dimension: int | None = None,
) -> list[CompleteGraphlet]:
    counts = feature_frequency(score_map)
    ranked_features = sorted(
        counts,
        key=lambda feature: (-int(counts[feature]), int(feature[0]), str(feature[1])),
    )
    requested_columns = tuple(str(col) for col in score_columns)
    graphlets: list[CompleteGraphlet] = []
    used: set[tuple[MutationFeature, ...]] = set()
    for dim in sorted({int(x) for x in dimensions}):
        pool = tuple(ranked_features[: max(int(top_feature_pool), dim)])
        found = 0
        checked = 0
        for combo in itertools.combinations(pool, dim):
            checked += 1
            if max_combinations_per_dimension is not None and checked > int(max_combinations_per_dimension):
                break
            features = tuple(sorted(combo, key=lambda feature: (int(feature[0]), str(feature[1]))))
            if features in used:
                continue
            complete = True
            for mask in node_masks(features):
                values = score_map.get(mask)
                if values is None or not all(col in values for col in requested_columns):
                    complete = False
                    break
            if not complete:
                continue
            graphlet_id = f"{assay_id}__exact_g{len(graphlets):04d}__d{dim}"
            graphlets.append(CompleteGraphlet(graphlet_id=graphlet_id, dimension=dim, features=features))
            used.add(features)
            found += 1
            if found >= int(max_graphlets_per_dimension):
                break
    return graphlets


def values_for_graphlet(
    score_map: Mapping[MutationMask, Mapping[str, float]],
    graphlet: CompleteGraphlet,
    *,
    score_column: str,
) -> np.ndarray:
    return np.asarray(
        [float(score_map[mask][score_column]) for mask in node_masks(graphlet.features)],
        dtype=np.float64,
    )


def hypercube_edges(dim: int) -> list[tuple[int, int, int]]:
    edges: list[tuple[int, int, int]] = []
    for source in range(1 << int(dim)):
        for bit in range(int(dim)):
            target = source ^ (1 << bit)
            if source < target:
                edges.append((source, target, bit))
    return edges


def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return math.nan
    if float(np.std(x)) <= 1e-12 or float(np.std(y)) <= 1e-12:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def exact_random_walk_metrics(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    source: list[float] = []
    target: list[float] = []
    deltas: list[float] = []
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            nxt = node ^ (1 << bit)
            source.append(float(arr[node]))
            target.append(float(arr[nxt]))
            deltas.append(float(arr[nxt] - arr[node]))
    x = np.asarray(source, dtype=np.float64)
    y = np.asarray(target, dtype=np.float64)
    lag1 = _safe_corr(x, y)
    corr_len = float(-1.0 / math.log(lag1)) if np.isfinite(lag1) and 0.0 < lag1 < 1.0 else math.nan
    return {
        "random_walk_lag1_autocorrelation": lag1,
        "random_walk_autocorrelation_length": corr_len,
        "random_walk_step_sd": float(np.std(deltas)) if deltas else math.nan,
    }


def sign_epistasis_metrics(values: np.ndarray, *, dim: int, resolution: float = 1e-12) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    total = 0
    sign_count = 0
    reciprocal_count = 0
    pair_rows: list[dict[str, Any]] = []
    for left, right in itertools.combinations(range(int(dim)), 2):
        other_bits = [bit for bit in range(int(dim)) if bit not in {left, right}]
        for other_state in range(1 << len(other_bits)):
            base = 0
            for offset, bit in enumerate(other_bits):
                if (other_state >> offset) & 1:
                    base |= 1 << bit
            n00 = base
            n10 = base | (1 << left)
            n01 = base | (1 << right)
            n11 = base | (1 << left) | (1 << right)
            left_low = float(arr[n10] - arr[n00])
            left_high = float(arr[n11] - arr[n01])
            right_low = float(arr[n01] - arr[n00])
            right_high = float(arr[n11] - arr[n10])
            left_flip = (left_low > resolution and left_high < -resolution) or (left_low < -resolution and left_high > resolution)
            right_flip = (right_low > resolution and right_high < -resolution) or (right_low < -resolution and right_high > resolution)
            sign = bool(left_flip or right_flip)
            reciprocal = bool(left_flip and right_flip)
            sign_count += int(sign)
            reciprocal_count += int(reciprocal)
            total += 1
            pair_rows.append(
                {
                    "left_bit": int(left),
                    "right_bit": int(right),
                    "background_node": int(base),
                    "sign_epistasis": sign,
                    "reciprocal_sign_epistasis": reciprocal,
                    "plaquette_effect": float(arr[n11] - arr[n10] - arr[n01] + arr[n00]),
                }
            )
    return {
        "total_pair_backgrounds": float(total),
        "sign_epistasis_count": float(sign_count),
        "reciprocal_sign_epistasis_count": float(reciprocal_count),
        "sign_epistasis_fraction": float(sign_count / total) if total else math.nan,
        "reciprocal_sign_epistasis_fraction": float(reciprocal_count / total) if total else math.nan,
        "pair_rows": pair_rows,
    }


def local_extrema_metrics(values: np.ndarray, *, dim: int, resolution: float = 1e-12) -> dict[str, Any]:
    """Strict peaks and resolution-connected plateaus (aggregate_v2)."""
    return extrema_metrics(values, dim=dim, resolution=resolution)


def accessible_path_metrics(values: np.ndarray, *, dim: int, max_dim: int = 8, resolution: float = 1e-12) -> dict[str, float]:
    """Exact shortest-path fraction by dynamic programming; max_dim is retained for compatibility."""
    arr = cube_values(values, dim)
    start, end = (0, len(arr) - 1) if arr[-1] >= arr[0] else (len(arr) - 1, 0)
    total = math.factorial(dim)
    fraction = accessible_fraction(arr, dim=dim, resolution=resolution)
    return {"accessible_path_start_node": float(start), "accessible_path_end_node": float(end),
            "accessible_path_total_shortest_paths": float(total), "accessible_path_count": round(fraction * total),
            "accessible_path_fraction": fraction, "accessible_path_endpoint_delta": float(arr[end]-arr[start])}


def basin_metrics(values: np.ndarray, *, dim: int, resolution: float = 1e-12) -> dict[str, Any]:
    """Greedy basins after contracting resolution-connected neutral components."""
    arr = cube_values(values, dim)
    groups, labels, _ = plateau_components(arr, dim, resolution)
    successors = {}
    for group in groups:
        label = int(labels[group[0]])
        ceiling = float(np.max(arr[group]))
        exits = [int(i ^ (1 << bit)) for i in group for bit in range(dim)
                 if labels[i ^ (1 << bit)] != label and arr[i ^ (1 << bit)] > ceiling + resolution]
        successors[label] = int(labels[max(exits, key=lambda n: (arr[n], -n))]) if exits else label
    counts = Counter()
    for group in groups:
        label = int(labels[group[0]])
        while successors[label] != label:
            label = successors[label]
        counts[label] += len(group)
    probs = np.array(list(counts.values()), dtype=float) / len(arr)
    return {"basin_count": float(len(counts)), "basin_entropy": -float(np.sum(probs * np.log2(probs))),
            "basin_largest_fraction": float(probs.max()), "basin_attractors": ";".join(map(str, sorted(counts)))}


def walsh_decomposition(values: np.ndarray, *, dim: int) -> tuple[list[dict[str, float]], dict[int, float]]:
    coeff = np.asarray(values, dtype=np.float64).copy()
    n = 1 << int(dim)
    if coeff.size != n:
        raise ValueError(f"expected {n} values for dim={dim}, got {coeff.size}")
    step = 1
    while step < n:
        for start in range(0, n, step * 2):
            left = coeff[start : start + step].copy()
            right = coeff[start + step : start + step * 2].copy()
            coeff[start : start + step] = left + right
            coeff[start + step : start + step * 2] = left - right
        step *= 2
    coeff = coeff / float(n)
    rows: list[dict[str, float]] = []
    energy_by_order: Counter[int] = Counter()
    for mask, value in enumerate(coeff):
        order = int(mask.bit_count())
        energy = float(value * value)
        rows.append(
            {
                "mask": float(mask),
                "order": float(order),
                "coefficient": float(value),
                "energy": energy,
            }
        )
        if order > 0:
            energy_by_order[order] += energy
    total = float(sum(energy_by_order.values()))
    fractions = {
        order: float(energy / total) if total > 1e-18 else 0.0
        for order, energy in sorted(energy_by_order.items())
    }
    return rows, fractions


def spectral_metrics(values: np.ndarray, *, dim: int, variance_floor: float = 1e-18) -> dict[str, float]:
    cube_values(values, dim)
    rows, fractions = walsh_decomposition(values, dim=dim)
    total_energy = float(sum(float(row["energy"]) for row in rows if int(row["order"]) > 0))
    resolved = total_energy > variance_floor
    fractions = {order: (sum(row["energy"] for row in rows if int(row["order"]) == order) / total_energy
                         if resolved else math.nan) for order in range(1, dim + 1)}
    high = float(sum(value for order, value in fractions.items() if int(order) >= 3)) if resolved else math.nan
    band_energy = {"low": 0.0, "mid": 0.0, "high": 0.0}
    for row in rows:
        order = int(row["order"])
        if order == 0:
            continue
        q = order / max(int(dim), 1)
        if q <= 0.33:
            band = "low"
        elif q <= 0.67:
            band = "mid"
        else:
            band = "high"
        band_energy[band] += float(row["energy"])
    return {
        "total_variance": total_energy,
        "spectral_resolved": resolved,
        "walsh_epistatic_fraction": sum(v for k, v in fractions.items() if k >= 2) if resolved else math.nan,
        "walsh_high_order_fraction": high,
        **{f"walsh_fraction_order_{k}": v for k, v in fractions.items()},
        **{f"walsh_energy_order_{k}": sum(row["energy"] for row in rows if int(row["order"]) == k) for k in range(1, dim + 1)},
        "walsh_high_order_present": float(high > 1e-6) if resolved else math.nan,
        **{f"walsh_variance_order_{order}": float(value) for order, value in fractions.items()},
        "gft_low_energy_fraction": float(band_energy["low"] / total_energy) if resolved else math.nan,
        "gft_mid_energy_fraction": float(band_energy["mid"] / total_energy) if resolved else math.nan,
        "gft_high_energy_fraction": float(band_energy["high"] / total_energy) if resolved else math.nan,
    }


def _concentration(values: Iterable[float], *, prefix: str) -> dict[str, float]:
    arr = np.asarray([max(float(value), 0.0) for value in values if np.isfinite(float(value))], dtype=np.float64)
    n = int(arr.size)
    total = float(np.sum(arr))
    if n == 0 or total <= 1e-18:
        return {
            f"{prefix}_count": float(n),
            f"{prefix}_total": 0.0,
            f"{prefix}_effective_fraction": math.nan,
            f"{prefix}_top1_fraction": math.nan,
            f"{prefix}_top10_fraction": math.nan,
        }
    probs = np.sort(arr / total)[::-1]
    hhi = float(np.sum(probs * probs))

    def top_fraction(frac: float) -> float:
        count = max(1, int(math.ceil(frac * n)))
        return float(np.sum(probs[:count]))

    return {
        f"{prefix}_count": float(n),
        f"{prefix}_total": total,
        f"{prefix}_effective_fraction": float((1.0 / hhi) / n) if hhi > 1e-18 else math.nan,
        f"{prefix}_top1_fraction": float(probs[0]),
        f"{prefix}_top10_fraction": top_fraction(0.10),
    }


def local_energy_metrics(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    edges = hypercube_edges(dim)
    edge_energy = np.asarray([(float(arr[v]) - float(arr[u])) ** 2 for u, v, _bit in edges], dtype=np.float64)
    node_incident = np.zeros(1 << int(dim), dtype=np.float64)
    site_energy = np.zeros(int(dim), dtype=np.float64)
    for (u, v, bit), energy in zip(edges, edge_energy):
        node_incident[u] += float(energy)
        node_incident[v] += float(energy)
        site_energy[bit] += float(energy)

    incident: dict[int, list[int]] = {}
    for idx, (u, v, _bit) in enumerate(edges):
        incident.setdefault(u, []).append(idx)
        incident.setdefault(v, []).append(idx)
    adjacent_pairs = sorted(
        {
            tuple(sorted((left, right)))
            for edge_ids in incident.values()
            for left_index, left in enumerate(edge_ids)
            for right in edge_ids[left_index + 1 :]
        }
    )
    if adjacent_pairs:
        left = np.asarray([edge_energy[i] for i, _j in adjacent_pairs], dtype=np.float64)
        right = np.asarray([edge_energy[j] for _i, j in adjacent_pairs], dtype=np.float64)
        pearson = _safe_corr(left, right)
        mean = float(np.mean(edge_energy))
        denom = float(np.sum((edge_energy - mean) ** 2))
        moran = (
            float(len(edge_energy) / len(adjacent_pairs) * np.sum((left - mean) * (right - mean)) / denom)
            if denom > 1e-18
            else math.nan
        )
    else:
        pearson = math.nan
        moran = math.nan

    plaquette_energy: list[float] = []
    for left, right in itertools.combinations(range(int(dim)), 2):
        other_bits = [bit for bit in range(int(dim)) if bit not in {left, right}]
        for other_state in range(1 << len(other_bits)):
            base = 0
            for offset, bit in enumerate(other_bits):
                if (other_state >> offset) & 1:
                    base |= 1 << bit
            n00 = base
            n10 = base | (1 << left)
            n01 = base | (1 << right)
            n11 = base | (1 << left) | (1 << right)
            plaquette_energy.append(float(arr[n11] - arr[n10] - arr[n01] + arr[n00]) ** 2)

    return {
        "edge_energy_adjacent_pearson": pearson,
        "edge_energy_adjacent_moran_i": moran,
        **_concentration(edge_energy, prefix="edge_energy"),
        **_concentration(node_incident, prefix="node_incident"),
        **_concentration(site_energy, prefix="site_energy"),
        **_concentration(plaquette_energy, prefix="plaquette_energy"),
    }


def exact_summary(values: np.ndarray, *, dim: int, accessibility_max_dim: int = 8, resolution: float = 1e-12) -> dict[str, Any]:
    arr = cube_values(values, dim)
    return {
        "metric_version": METRIC_VERSION,
        "fitness_min": float(np.min(arr)),
        "fitness_mean": float(np.mean(arr)),
        "fitness_max": float(np.max(arr)),
        "fitness_sd": float(np.std(arr)),
        **spectral_metrics(arr, dim=dim, variance_floor=max(1e-18, resolution**2)),
        **{key: value for key, value in sign_epistasis_metrics(arr, dim=dim, resolution=resolution).items() if key != "pair_rows"},
        **local_extrema_metrics(arr, dim=dim, resolution=resolution),
        **accessible_path_metrics(arr, dim=dim, max_dim=accessibility_max_dim, resolution=resolution),
        **basin_metrics(arr, dim=dim, resolution=resolution),
        **exact_random_walk_metrics(arr, dim=dim),
        **local_energy_metrics(arr, dim=dim),
    }


def component_failure_table(
    metric_rows: pd.DataFrame,
    *,
    target_candidate_model: str,
    metric_columns: Sequence[str],
    group_columns: Sequence[str] = ("candidate_model", "score_column"),
    empirical_scale_floor_fraction: float = 0.10,
    absolute_scale_floor: float = 1e-6,
) -> pd.DataFrame:
    if metric_rows.empty:
        return pd.DataFrame()
    group_columns = tuple(group_columns)
    target = metric_rows[metric_rows["candidate_model"] == target_candidate_model]
    rows: list[dict[str, Any]] = []
    for keys, model_sub in metric_rows.groupby(list(group_columns), dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        group = {col: value for col, value in zip(group_columns, keys)}
        candidate = str(group.get("candidate_model", ""))
        if candidate == target_candidate_model:
            continue
        score_column = str(group.get("score_column", ""))
        target_sub = target
        if "score_column" in group_columns and "score_column" in target:
            target_sub = target[target["score_column"] == score_column]
        for metric in metric_columns:
            if metric not in model_sub or metric not in target_sub:
                continue
            model_values = pd.to_numeric(model_sub[metric], errors="coerce").dropna()
            target_values = pd.to_numeric(target_sub[metric], errors="coerce").dropna()
            if model_values.empty or target_values.empty:
                continue
            model_mean = float(model_values.mean())
            target_mean = float(target_values.mean())
            target_sd = float(target_values.std(ddof=0)) if len(target_values) >= 2 else 0.0
            scale = max(target_sd, abs(target_mean) * float(empirical_scale_floor_fraction), float(absolute_scale_floor))
            signed_error = float((model_mean - target_mean) / scale)
            rows.append(
                {
                    **group,
                    "metric": str(metric),
                    "n_model_graphlets": int(len(model_values)),
                    "n_empirical_graphlets": int(len(target_values)),
                    "model_mean": model_mean,
                    "empirical_mean": target_mean,
                    "empirical_sd": target_sd,
                    "scale": scale,
                    "signed_scaled_error": signed_error,
                    "abs_scaled_error": abs(signed_error),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["candidate_model", "score_column", "abs_scaled_error"], ascending=[True, True, False])
    return out


def model_distance_table(failure_rows: pd.DataFrame) -> pd.DataFrame:
    if failure_rows.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for keys, sub in failure_rows.groupby(["candidate_model", "score_column"], dropna=False):
        values = pd.to_numeric(sub["abs_scaled_error"], errors="coerce").dropna()
        if values.empty:
            continue
        severe = sub[sub["abs_scaled_error"] >= 2.0]
        top = sub.sort_values("abs_scaled_error", ascending=False).head(5)
        rows.append(
            {
                "candidate_model": keys[0],
                "score_column": keys[1],
                "mean_abs_scaled_error": float(values.mean()),
                "median_abs_scaled_error": float(values.median()),
                "max_abs_scaled_error": float(values.max()),
                "n_components": int(len(values)),
                "n_components_ge_2sd": int(len(severe)),
                "top_failed_metrics": ";".join(str(x) for x in top["metric"].tolist()),
            }
        )
    return pd.DataFrame(rows).sort_values("mean_abs_scaled_error", ascending=True)
