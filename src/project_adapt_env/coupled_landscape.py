"""Utilities for deterministic coupled-landscape graphlet scans.

The functions here deliberately operate on dense binary graphlet value arrays.
They are meant for inexpensive machinery tests before fitting a full empirical
coupled-landscape model.
"""

from __future__ import annotations

import itertools
import math
import hashlib
from collections import defaultdict
from typing import Any

import numpy as np
import pandas as pd


Feature = tuple[int, str]
Module = dict[str, Any]


def sigmoid(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-x))


def centered(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    return arr - float(np.mean(arr))


def zscore(values: np.ndarray) -> np.ndarray:
    arr = centered(np.asarray(values, dtype=np.float64))
    sd = float(np.std(arr))
    if sd <= 1e-12:
        return np.zeros_like(arr, dtype=np.float64)
    return arr / sd


def safe_corr(xs: np.ndarray, ys: np.ndarray) -> float:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(ok)) < 3:
        return math.nan
    x = x[ok]
    y = y[ok]
    if float(np.std(x)) <= 1e-12 or float(np.std(y)) <= 1e-12:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def walsh_synthesize(coefficients: dict[int, float], *, dim: int) -> np.ndarray:
    values = np.zeros(1 << int(dim), dtype=np.float64)
    for mask, coeff in coefficients.items():
        mask = int(mask)
        coeff = float(coeff)
        for node in range(values.size):
            parity = (node & mask).bit_count() % 2
            values[node] += coeff * (-1.0 if parity else 1.0)
    return values


def low_order_auxiliary_core(
    *,
    dim: int,
    rng: np.random.RandomState,
    base_norm: np.ndarray,
) -> np.ndarray:
    coeffs = {
        1 << bit: float(rng.normal())
        for bit in range(int(dim))
    }
    raw = zscore(walsh_synthesize(coeffs, dim=dim))
    base = zscore(base_norm)
    corr = safe_corr(raw, base)
    if np.isfinite(corr):
        raw = raw - corr * base
    return zscore(raw)


def node_edge_gradients(values: np.ndarray, *, dim: int) -> np.ndarray:
    out: list[float] = []
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            neighbor = node ^ (1 << bit)
            if node < neighbor:
                out.append(float(values[neighbor] - values[node]))
    return np.asarray(out, dtype=np.float64)


def additive_walsh_coefficients(values: np.ndarray, *, dim: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    coeffs: list[float] = []
    n = values.size
    for bit in range(int(dim)):
        mask = 1 << bit
        total = 0.0
        for node in range(n):
            parity = (node & mask).bit_count() % 2
            total += float(values[node]) * (-1.0 if parity else 1.0)
        coeffs.append(total / float(n))
    return np.asarray(coeffs, dtype=np.float64)


def plaquette_energy_by_pair(values: np.ndarray, *, dim: int) -> dict[tuple[int, int], float]:
    values = np.asarray(values, dtype=np.float64)
    out: dict[tuple[int, int], float] = {}
    for i, j in itertools.combinations(range(int(dim)), 2):
        total = 0.0
        count = 0
        for node in range(values.size):
            if ((node >> i) & 1) or ((node >> j) & 1):
                continue
            diff = (
                float(values[node ^ (1 << i) ^ (1 << j)])
                - float(values[node ^ (1 << i)])
                - float(values[node ^ (1 << j)])
                + float(values[node])
            )
            total += diff * diff
            count += 1
        out[(i, j)] = float(total / max(count, 1))
    return out


def edge_energy_by_bit(values: np.ndarray, *, dim: int) -> dict[int, float]:
    values = np.asarray(values, dtype=np.float64)
    out: dict[int, float] = {}
    for bit in range(int(dim)):
        total = 0.0
        count = 0
        for node in range(values.size):
            neighbor = node ^ (1 << bit)
            if node < neighbor:
                delta = float(values[neighbor] - values[node])
                total += delta * delta
                count += 1
        out[bit] = float(total / max(count, 1))
    return out


def concentration_summary(energy_by_module: dict[str, float], *, prefix: str) -> dict[str, float]:
    values = np.asarray([max(float(v), 0.0) for v in energy_by_module.values()], dtype=np.float64)
    total = float(np.sum(values))
    if total <= 1e-18:
        return {
            f"{prefix}_total_energy": 0.0,
            f"{prefix}_effective_module_count": math.nan,
            f"{prefix}_top10_fraction": math.nan,
            f"{prefix}_top25_fraction": math.nan,
            f"{prefix}_hhi": math.nan,
            f"{prefix}_normalized_entropy": math.nan,
        }
    probs = values / total
    sorted_probs = np.sort(probs)[::-1]
    n10 = max(1, int(math.ceil(0.10 * len(sorted_probs))))
    n25 = max(1, int(math.ceil(0.25 * len(sorted_probs))))
    hhi = float(np.sum(probs * probs))
    entropy = -float(np.sum(probs * np.log(np.maximum(probs, 1e-18))))
    return {
        f"{prefix}_total_energy": total,
        f"{prefix}_effective_module_count": float(1.0 / hhi) if hhi > 1e-18 else math.nan,
        f"{prefix}_top10_fraction": float(np.sum(sorted_probs[:n10])),
        f"{prefix}_top25_fraction": float(np.sum(sorted_probs[:n25])),
        f"{prefix}_hhi": hhi,
        f"{prefix}_normalized_entropy": (
            float(entropy / math.log(len(probs))) if len(probs) > 1 else 0.0
        ),
    }


def select_empirical_modules(
    graphlets: list[dict[str, Any]],
    innovation_by_graphlet: dict[str, np.ndarray],
    *,
    max_modules: int,
    module_size: int,
) -> list[Module]:
    """Select calibration modules from empirical innovation plaquette energy."""
    pair_scores: defaultdict[tuple[Feature, Feature], float] = defaultdict(float)
    pair_contexts: defaultdict[tuple[Feature, Feature], defaultdict[Feature, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    feature_counts: defaultdict[Feature, int] = defaultdict(int)
    for graphlet in graphlets:
        if graphlet.get("graphlet_split") != "calibration":
            continue
        features = tuple(graphlet["features"])
        values = innovation_by_graphlet.get(str(graphlet["graphlet_id"]))
        if values is None:
            continue
        for feature in features:
            feature_counts[feature] += 1
        for (i, j), energy in plaquette_energy_by_pair(values, dim=int(graphlet["graphlet_dimension"])).items():
            pair = tuple(sorted((features[i], features[j])))  # type: ignore[arg-type]
            pair_scores[pair] += float(energy)
            for feature in features:
                if feature not in pair:
                    pair_contexts[pair][feature] += float(energy)

    modules: list[Module] = []
    used: set[tuple[Feature, ...]] = set()
    for rank, (pair, score) in enumerate(
        sorted(pair_scores.items(), key=lambda item: (-item[1], str(item[0])))
    ):
        selected = list(pair)
        if int(module_size) > len(selected):
            contexts = pair_contexts[pair]
            extras = sorted(
                contexts,
                key=lambda feature: (-float(contexts[feature]), -feature_counts[feature], str(feature)),
            )
            for feature in extras:
                if feature not in selected:
                    selected.append(feature)
                if len(selected) >= int(module_size):
                    break
        key = tuple(sorted(selected))
        if key in used or len(key) < 2:
            continue
        used.add(key)
        modules.append(
            {
                "module_id": f"empirical_calibration_plaq_{len(modules):03d}",
                "module_source": "empirical_calibration_plaquette",
                "features": list(key),
                "calibration_pair_score": float(score),
                "rank": int(rank),
            }
        )
        if len(modules) >= int(max_modules):
            break
    return modules


def control_modules(
    modules: list[Module],
    feature_pool: list[Feature],
    *,
    source: str,
    rng: np.random.RandomState,
) -> list[Module]:
    if not modules:
        return []
    pool = list(feature_pool)
    if source == "shuffled_module":
        shuffled = pool.copy()
        rng.shuffle(shuffled)
        mapping = {feature: shuffled[idx % len(shuffled)] for idx, feature in enumerate(pool)}
        out = []
        for module in modules:
            out.append(
                {
                    **module,
                    "module_id": str(module["module_id"]).replace("empirical_calibration_plaq", "shuffled"),
                    "module_source": source,
                    "features": [mapping.get(feature, feature) for feature in module["features"]],
                }
            )
        return out
    if source == "random_module":
        out = []
        size = max(2, len(modules[0]["features"]))
        for idx, module in enumerate(modules):
            chosen = list(rng.choice(len(pool), size=min(size, len(pool)), replace=False))
            out.append(
                {
                    **module,
                    "module_id": f"random_module_{idx:03d}",
                    "module_source": source,
                    "features": [pool[item] for item in chosen],
                }
            )
        return out
    raise ValueError(f"Unknown control module source: {source}")


def module_delta_values(
    *,
    graphlet_features: tuple[Feature, ...],
    modules: list[Module],
    dim: int,
    module_residual_scale: float,
    residual_density: float,
    orders: list[int],
    rng: np.random.RandomState,
) -> tuple[np.ndarray, dict[str, float], dict[str, float]]:
    coeffs: dict[int, float] = {}
    active_modules: list[Module] = []
    module_coeff_energy: defaultdict[str, float] = defaultdict(float)
    index_by_feature = {feature: idx for idx, feature in enumerate(graphlet_features)}
    for module in modules:
        local = [index_by_feature[f] for f in module["features"] if f in index_by_feature]
        if len(local) < 2:
            continue
        active_modules.append(module)
        for order in orders:
            order = int(order)
            if order < 2 or order > len(local):
                continue
            for combo in itertools.combinations(local, order):
                if rng.rand() > float(residual_density):
                    continue
                mask = sum(1 << bit for bit in combo)
                coeff = float(rng.normal(scale=float(module_residual_scale) / math.sqrt(max(order, 1))))
                coeffs[mask] = coeffs.get(mask, 0.0) + coeff
                module_coeff_energy[str(module["module_id"])] += coeff * coeff
    values = walsh_synthesize(coeffs, dim=dim)
    values = centered(values)
    if values.size and float(np.std(values)) > 1e-12:
        values = values / float(np.std(values)) * float(module_residual_scale)
    active = {
        "n_active_modules": float(len(active_modules)),
        "n_module_coefficients": float(len(coeffs)),
    }
    return values, active, dict(module_coeff_energy)


def _stable_uint32(*parts: object) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:4], byteorder="big", signed=False)


def _feature_key(feature: Feature) -> str:
    return f"{int(feature[0])}:{str(feature[1])}"


def module_residual_coefficients(
    modules: list[Module],
    *,
    orders: list[int],
    residual_density: float,
    coefficient_seed: int,
    coefficient_scale: float = 1.0,
) -> dict[tuple[Feature, ...], float]:
    """Create deterministic coefficients for module-local residual features.

    Unlike ``module_delta_values``, this coefficient map is stable across assay
    variants and graphlets, so a module-local ``d(x)=b(x)-a(x)`` field is
    evaluated consistently on the full assay and on dense graphlets.
    """
    density = max(0.0, min(float(residual_density), 1.0))
    coeffs: dict[tuple[Feature, ...], float] = {}
    for module in modules:
        features = tuple(sorted(tuple(module["features"]), key=lambda item: (int(item[0]), str(item[1]))))
        module_id = str(module.get("module_id", "module"))
        for order in orders:
            order = int(order)
            if order < 1 or order > len(features):
                continue
            for combo in itertools.combinations(features, order):
                combo_key = ";".join(_feature_key(feature) for feature in combo)
                keep_seed = _stable_uint32("module-residual-keep", coefficient_seed, module_id, order, combo_key)
                if keep_seed / float(2**32 - 1) > density:
                    continue
                coeff_seed = _stable_uint32("module-residual-coeff", coefficient_seed, module_id, order, combo_key)
                rng = np.random.RandomState(coeff_seed)
                coeff = float(rng.normal(scale=float(coefficient_scale) / math.sqrt(max(order, 1))))
                coeffs[tuple(combo)] = coeffs.get(tuple(combo), 0.0) + coeff
    return coeffs


def module_residual_values_for_feature_sets(
    feature_sets: list[tuple[Feature, ...]],
    coefficients: dict[tuple[Feature, ...], float],
) -> np.ndarray:
    """Evaluate module-local residual coefficients on arbitrary feature sets."""
    values = np.zeros(len(feature_sets), dtype=np.float64)
    if not coefficients:
        return values
    normalized_sets = [set(tuple(item) for item in features) for features in feature_sets]
    for idx, active_features in enumerate(normalized_sets):
        total = 0.0
        for combo, coeff in coefficients.items():
            if all(feature in active_features for feature in combo):
                total += float(coeff)
        values[idx] = total
    return values


def module_tied_difference_values(
    *,
    base_values: np.ndarray,
    feature_sets: list[tuple[Feature, ...]],
    modules: list[Module],
    parameter_set: dict[str, Any],
    coefficient_seed: int,
) -> tuple[np.ndarray, dict[str, float], dict[tuple[Feature, ...], float]]:
    """Construct a module-tied ``d(x)=b(x)-a(x)`` field in base-value units."""
    base = np.asarray(base_values, dtype=np.float64)
    if len(feature_sets) != int(base.size):
        raise ValueError("feature_sets length must match base_values length")
    base_sd = max(float(np.std(base)), 1e-8)
    base_norm = zscore(base)
    use_modules = bool(parameter_set.get("use_modules", True))
    coeffs = module_residual_coefficients(
        modules if use_modules else [],
        orders=[int(item) for item in parameter_set.get("module_residual_orders", [2])],
        residual_density=float(parameter_set.get("module_residual_density", 1.0)),
        coefficient_seed=int(coefficient_seed),
        coefficient_scale=1.0,
    )
    module_raw = module_residual_values_for_feature_sets(feature_sets, coeffs)
    module_norm = zscore(module_raw)
    bias = float(parameter_set.get("difference_bias_sd", 0.0))
    base_slope = float(parameter_set.get("difference_base_slope", 0.0))
    module_scale = float(parameter_set.get("difference_module_scale", 0.0))
    raw_difference = bias + base_slope * base_norm + module_scale * module_norm
    difference = base_sd * raw_difference
    info = {
        "difference_bias_sd": bias,
        "difference_base_slope": base_slope,
        "difference_module_scale": module_scale,
        "difference_sd": float(np.std(difference)),
        "difference_to_base_sd_ratio": float(np.std(difference) / base_sd),
        "difference_to_base_correlation": safe_corr(difference, base),
        "module_raw_sd": float(np.std(module_raw)),
        "module_difference_component_sd": float(np.std(base_sd * module_scale * module_norm)),
        "n_module_coefficients": float(len(coeffs)),
        "n_modules": float(len(modules if use_modules else [])),
    }
    return difference, info, coeffs


def softmin_from_difference(
    *,
    base_values: np.ndarray,
    difference_values: np.ndarray,
    tau: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate a soft-min state operator using explicit ``d=b-a`` values."""
    base = np.asarray(base_values, dtype=np.float64)
    difference = np.asarray(difference_values, dtype=np.float64)
    auxiliary = base + difference
    values, innovation, auxiliary_weight = soft_min_values(
        base_values=base,
        auxiliary_values=auxiliary,
        tau=tau,
    )
    return values, innovation, auxiliary_weight, auxiliary


def edge_difference_values(values: np.ndarray, *, dim: int) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    out: list[float] = []
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            neighbor = node ^ (1 << bit)
            if node < neighbor:
                out.append(float(arr[neighbor] - arr[node]))
    return np.asarray(out, dtype=np.float64)


def finite_difference_budget_metrics(
    *,
    base_values: np.ndarray,
    innovation: np.ndarray,
    auxiliary_weight: np.ndarray,
    difference_values: np.ndarray,
    tau: float,
    dim: int,
    empirical_residual: np.ndarray | None = None,
) -> dict[str, float]:
    """Measure whether a soft state operator is creating plausible local energy."""
    base = np.asarray(base_values, dtype=np.float64)
    innov = np.asarray(innovation, dtype=np.float64)
    weight = np.asarray(auxiliary_weight, dtype=np.float64)
    difference = np.asarray(difference_values, dtype=np.float64)
    base_var = max(float(np.var(base)), 1e-12)
    innov_edge = float(sum(edge_energy_by_bit(innov, dim=dim).values()))
    innov_plaq = float(sum(plaquette_energy_by_pair(innov, dim=dim).values()))
    base_edge = float(sum(edge_energy_by_bit(base, dim=dim).values()))
    base_plaq = float(sum(plaquette_energy_by_pair(base, dim=dim).values()))
    weight_edge = float(sum(edge_energy_by_bit(weight, dim=dim).values()))
    weight_plaq = float(sum(plaquette_energy_by_pair(weight, dim=dim).values()))
    diff_edges = np.abs(edge_difference_values(difference, dim=dim))
    temp = max(float(tau), 1e-12)
    metrics = {
        "innovation_var_ratio_to_base": float(np.var(innov) / base_var),
        "innovation_edge_energy": innov_edge,
        "innovation_plaquette_energy": innov_plaq,
        "innovation_edge_ratio_to_base": float(innov_edge / max(base_edge, 1e-12)),
        "innovation_plaquette_ratio_to_base": float(innov_plaq / max(base_plaq, 1e-12)),
        "state_weight_edge_energy": weight_edge,
        "state_weight_plaquette_energy": weight_plaq,
        "difference_edge_median_abs": float(np.median(diff_edges)) if diff_edges.size else math.nan,
        "difference_edge_p90_abs": float(np.quantile(diff_edges, 0.9)) if diff_edges.size else math.nan,
        "difference_edge_median_abs_over_tau": (
            float(np.median(diff_edges) / temp) if diff_edges.size else math.nan
        ),
        "difference_edge_p90_abs_over_tau": (
            float(np.quantile(diff_edges, 0.9) / temp) if diff_edges.size else math.nan
        ),
    }
    if empirical_residual is not None:
        resid = np.asarray(empirical_residual, dtype=np.float64)
        resid_edge = float(sum(edge_energy_by_bit(resid, dim=dim).values()))
        resid_plaq = float(sum(plaquette_energy_by_pair(resid, dim=dim).values()))
        metrics.update(
            {
                "empirical_residual_edge_energy": resid_edge,
                "empirical_residual_plaquette_energy": resid_plaq,
                "innovation_edge_ratio_to_empirical_residual": float(
                    innov_edge / max(resid_edge, 1e-12)
                ),
                "innovation_plaquette_ratio_to_empirical_residual": float(
                    innov_plaq / max(resid_plaq, 1e-12)
                ),
            }
        )
    return metrics


def budget_penalty_terms(metrics: dict[str, float], config: dict[str, Any]) -> dict[str, float]:
    """Quadratic penalties for implausibly sharp or high-energy innovation."""

    def above(name: str, max_key: str, lambda_key: str) -> float:
        value = float(metrics.get(name, math.nan))
        maximum = float(config.get(max_key, math.inf))
        weight = float(config.get(lambda_key, 0.0))
        if not np.isfinite(value) or not np.isfinite(maximum):
            return 0.0
        return float(weight * max(0.0, value - maximum) ** 2)

    return {
        "penalty_high_innovation_var_ratio": above(
            "innovation_var_ratio_to_base",
            "innovation_var_ratio_max",
            "lambda_innovation_var",
        ),
        "penalty_high_edge_energy_ratio": above(
            "innovation_edge_ratio_to_empirical_residual",
            "edge_energy_ratio_to_empirical_max",
            "lambda_edge_energy",
        ),
        "penalty_high_plaquette_energy_ratio": above(
            "innovation_plaquette_ratio_to_empirical_residual",
            "plaquette_energy_ratio_to_empirical_max",
            "lambda_plaquette_energy",
        ),
        "penalty_high_state_weight_edge_energy": above(
            "state_weight_edge_energy",
            "state_weight_edge_energy_max",
            "lambda_state_boundary",
        ),
        "penalty_high_state_weight_plaquette_energy": above(
            "state_weight_plaquette_energy",
            "state_weight_plaquette_energy_max",
            "lambda_state_boundary",
        ),
        "penalty_low_temperature_relative_to_d": above(
            "difference_edge_median_abs_over_tau",
            "difference_edge_median_abs_over_tau_max",
            "lambda_relative_temperature",
        ),
    }


def coupled_values(
    *,
    base_values: np.ndarray,
    graphlet_features: tuple[Feature, ...],
    modules: list[Module],
    parameter_set: dict[str, Any],
    seed: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, float], dict[str, float]]:
    dim = len(graphlet_features)
    rng = np.random.RandomState(int(seed) % (2**32))
    base = np.asarray(base_values, dtype=np.float64)
    base_sd = max(float(np.std(base)), 1e-8)
    base_norm = zscore(base)
    rho = float(parameter_set.get("auxiliary_correlation", 0.75))
    rho = max(0.0, min(0.999, rho))
    if bool(parameter_set.get("uncoupled_auxiliary", False)):
        rho = 0.0
    smooth_core = low_order_auxiliary_core(dim=dim, rng=rng, base_norm=base_norm)
    delta, active, coeff_energy = module_delta_values(
        graphlet_features=graphlet_features,
        modules=modules if bool(parameter_set.get("use_modules", True)) else [],
        dim=dim,
        module_residual_scale=float(parameter_set.get("module_residual_scale", 0.25)),
        residual_density=float(parameter_set.get("module_residual_density", 1.0)),
        orders=[int(x) for x in parameter_set.get("module_residual_orders", [2, 3])],
        rng=rng,
    )
    aux = rho * base_norm + math.sqrt(max(1.0 - rho * rho, 0.0)) * smooth_core + delta
    gate_midpoint = float(parameter_set.get("operator_gate_midpoint", 0.0))
    gate_steepness = max(float(parameter_set.get("operator_gate_steepness", 1.0)), 1e-6)
    gate_direction = str(parameter_set.get("operator_gate_direction", "high"))
    gate_x = (base_norm - gate_midpoint) / gate_steepness
    if gate_direction == "low":
        gate_x = -gate_x
    gate = sigmoid(gate_x)
    innovation = float(parameter_set.get("operator_lambda", 0.25)) * base_sd * gate * aux
    innovation = centered(innovation)
    values = base + innovation
    active.update(
        {
            "innovation_sd": float(np.std(innovation)),
            "innovation_to_base_sd_ratio": float(np.std(innovation) / base_sd),
            "mean_gate": float(np.mean(gate)),
            "auxiliary_phenotype_correlation_to_base": safe_corr(aux, base_norm),
            "innovation_phenotype_correlation_to_base": safe_corr(innovation, base),
            "auxiliary_smooth_sd": float(np.std(smooth_core)),
            "module_delta_sd": float(np.std(delta)),
        }
    )
    return values, innovation, active, coeff_energy


def hard_min_values(
    *,
    base_values: np.ndarray,
    auxiliary_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combine two landscapes through a strict weakest-link minimum."""
    base = np.asarray(base_values, dtype=np.float64)
    aux = np.asarray(auxiliary_values, dtype=np.float64)
    values = np.minimum(base, aux)
    innovation = values - base
    auxiliary_limited = aux < base
    return values, innovation, auxiliary_limited.astype(bool)


def soft_min_values(
    *,
    base_values: np.ndarray,
    auxiliary_values: np.ndarray,
    tau: float,
    auxiliary_bias: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combine two landscapes with a temperature-controlled soft minimum."""
    base = np.asarray(base_values, dtype=np.float64)
    aux = np.asarray(auxiliary_values, dtype=np.float64) + float(auxiliary_bias)
    temp = max(float(tau), 1e-6)
    stacked = np.vstack([base, aux])
    m = np.min(stacked, axis=0)
    values = m - temp * np.log(np.sum(np.exp(-(stacked - m) / temp), axis=0))
    innovation = values - base
    weights = np.exp(-(stacked - m) / temp)
    probs = weights / np.maximum(np.sum(weights, axis=0), 1e-300)
    auxiliary_weight = probs[1]
    return values, innovation, auxiliary_weight


def boltzmann_two_state_values(
    *,
    base_values: np.ndarray,
    auxiliary_values: np.ndarray,
    tau: float,
    auxiliary_bias: float = 0.0,
    mode: str = "min",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Two-state Boltzmann operator; ``mode='min'`` is the soft-min free energy."""
    if str(mode).strip().lower() == "min":
        return soft_min_values(
            base_values=base_values,
            auxiliary_values=auxiliary_values,
            tau=tau,
            auxiliary_bias=auxiliary_bias,
        )
    if str(mode).strip().lower() != "max":
        raise ValueError(f"Unsupported Boltzmann mode: {mode!r}")
    base = np.asarray(base_values, dtype=np.float64)
    aux = np.asarray(auxiliary_values, dtype=np.float64) + float(auxiliary_bias)
    temp = max(float(tau), 1e-6)
    stacked = np.vstack([base, aux])
    m = np.max(stacked, axis=0)
    values = m + temp * np.log(np.sum(np.exp((stacked - m) / temp), axis=0))
    innovation = values - base
    weights = np.exp((stacked - m) / temp)
    probs = weights / np.maximum(np.sum(weights, axis=0), 1e-300)
    auxiliary_weight = probs[1]
    return values, innovation, auxiliary_weight


def soft_state_metrics(
    auxiliary_weight: np.ndarray,
    innovation: np.ndarray,
    *,
    dim: int,
) -> dict[str, float]:
    """Summarize soft state competition and transition-boundary energy."""
    weight = np.asarray(auxiliary_weight, dtype=np.float64)
    innov = np.asarray(innovation, dtype=np.float64)
    n_edges = 0
    edge_delta_total = 0.0
    soft_switch_edges = 0
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            neighbor = node ^ (1 << bit)
            if node >= neighbor:
                continue
            delta = abs(float(weight[neighbor]) - float(weight[node]))
            edge_delta_total += delta
            soft_switch_edges += int(delta >= 0.25)
            n_edges += 1

    n_plaquettes = 0
    soft_switch_plaquettes = 0
    switch_plaq_energy = 0.0
    total_plaq_energy = 0.0
    for i, j in itertools.combinations(range(int(dim)), 2):
        for node in range(1 << int(dim)):
            if ((node >> i) & 1) or ((node >> j) & 1):
                continue
            nodes = [node, node ^ (1 << i), node ^ (1 << j), node ^ (1 << i) ^ (1 << j)]
            weights = [float(weight[item]) for item in nodes]
            mixed = (max(weights) - min(weights)) >= 0.25
            diff = (
                float(innov[nodes[3]])
                - float(innov[nodes[1]])
                - float(innov[nodes[2]])
                + float(innov[nodes[0]])
            )
            energy = diff * diff
            n_plaquettes += 1
            soft_switch_plaquettes += int(mixed)
            total_plaq_energy += energy
            if mixed:
                switch_plaq_energy += energy

    return {
        "auxiliary_weight_mean": float(np.mean(weight)),
        "auxiliary_weight_sd": float(np.std(weight)),
        "auxiliary_weight_min": float(np.min(weight)) if weight.size else math.nan,
        "auxiliary_weight_max": float(np.max(weight)) if weight.size else math.nan,
        "soft_switch_edge_fraction": float(soft_switch_edges / n_edges) if n_edges else math.nan,
        "mean_auxiliary_weight_edge_delta": float(edge_delta_total / n_edges) if n_edges else math.nan,
        "soft_switch_plaquette_fraction": (
            float(soft_switch_plaquettes / n_plaquettes) if n_plaquettes else math.nan
        ),
        "soft_switch_boundary_plaquette_energy_fraction": (
            float(switch_plaq_energy / total_plaq_energy) if total_plaq_energy > 1e-18 else math.nan
        ),
    }


def limiter_switch_metrics(
    limiter_state: np.ndarray,
    innovation: np.ndarray,
    *,
    dim: int,
) -> dict[str, float]:
    """Summarize where hard-min state boundaries fall on graphlet edges/faces."""
    state = np.asarray(limiter_state, dtype=bool)
    innov = np.asarray(innovation, dtype=np.float64)
    n_edges = 0
    n_switch_edges = 0
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            neighbor = node ^ (1 << bit)
            if node >= neighbor:
                continue
            n_edges += 1
            n_switch_edges += int(bool(state[node]) != bool(state[neighbor]))

    n_plaquettes = 0
    n_switch_plaquettes = 0
    switch_plaq_energy = 0.0
    total_plaq_energy = 0.0
    for i, j in itertools.combinations(range(int(dim)), 2):
        for node in range(1 << int(dim)):
            if ((node >> i) & 1) or ((node >> j) & 1):
                continue
            nodes = [node, node ^ (1 << i), node ^ (1 << j), node ^ (1 << i) ^ (1 << j)]
            states = [bool(state[item]) for item in nodes]
            mixed = any(states) and not all(states)
            diff = (
                float(innov[nodes[3]])
                - float(innov[nodes[1]])
                - float(innov[nodes[2]])
                + float(innov[nodes[0]])
            )
            energy = diff * diff
            n_plaquettes += 1
            n_switch_plaquettes += int(mixed)
            total_plaq_energy += energy
            if mixed:
                switch_plaq_energy += energy

    return {
        "auxiliary_limited_fraction": float(np.mean(state)) if state.size else math.nan,
        "switch_edge_fraction": float(n_switch_edges / n_edges) if n_edges else math.nan,
        "switch_plaquette_fraction": (
            float(n_switch_plaquettes / n_plaquettes) if n_plaquettes else math.nan
        ),
        "switch_boundary_plaquette_energy_fraction": (
            float(switch_plaq_energy / total_plaq_energy) if total_plaq_energy > 1e-18 else math.nan
        ),
    }


def energy_by_graphlet_feature(
    *,
    graphlet_id: str,
    graphlet_split: str,
    graphlet_features: tuple[Feature, ...],
    innovation: np.ndarray,
    candidate_model: str,
    parameter_set_id: str,
) -> list[dict[str, Any]]:
    """Assign edge and plaquette innovation energy to graphlet features/pairs."""
    dim = len(graphlet_features)
    edge_by_bit = edge_energy_by_bit(innovation, dim=dim)
    plaq_by_pair = plaquette_energy_by_pair(innovation, dim=dim)
    rows: list[dict[str, Any]] = []
    for bit, feature in enumerate(graphlet_features):
        rows.append(
            {
                "graphlet_id": graphlet_id,
                "graphlet_split": graphlet_split,
                "candidate_model": candidate_model,
                "parameter_set_id": parameter_set_id,
                "energy_kind": "edge",
                "feature_key": f"{feature[0]}:{feature[1]}",
                "energy": float(edge_by_bit.get(bit, 0.0)),
            }
        )
    for (i, j), energy in plaq_by_pair.items():
        f_i = graphlet_features[int(i)]
        f_j = graphlet_features[int(j)]
        rows.append(
            {
                "graphlet_id": graphlet_id,
                "graphlet_split": graphlet_split,
                "candidate_model": candidate_model,
                "parameter_set_id": parameter_set_id,
                "energy_kind": "plaquette",
                "feature_key": f"{f_i[0]}:{f_i[1]}|{f_j[0]}:{f_j[1]}",
                "energy": float(energy),
            }
        )
    return rows


def energy_locality_for_modules(
    *,
    graphlet_id: str,
    graphlet_split: str,
    graphlet_features: tuple[Feature, ...],
    innovation: np.ndarray,
    modules: list[Module],
    candidate_model: str,
    parameter_set_id: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    dim = len(graphlet_features)
    edge_by_bit = edge_energy_by_bit(innovation, dim=dim)
    plaq_by_pair = plaquette_energy_by_pair(innovation, dim=dim)
    index_by_feature = {feature: idx for idx, feature in enumerate(graphlet_features)}
    edge_energy_by_module: defaultdict[str, float] = defaultdict(float)
    plaq_energy_by_module: defaultdict[str, float] = defaultdict(float)
    rows: list[dict[str, Any]] = []
    for module in modules:
        module_id = str(module["module_id"])
        local = [index_by_feature[f] for f in module["features"] if f in index_by_feature]
        edge_energy = float(sum(edge_by_bit.get(bit, 0.0) for bit in local))
        plaq_energy = 0.0
        for pair in itertools.combinations(local, 2):
            plaq_energy += float(plaq_by_pair.get(tuple(sorted(pair)), 0.0))
        edge_energy_by_module[module_id] += edge_energy
        plaq_energy_by_module[module_id] += plaq_energy
        rows.append(
            {
                "graphlet_id": graphlet_id,
                "graphlet_split": graphlet_split,
                "candidate_model": candidate_model,
                "parameter_set_id": parameter_set_id,
                "module_id": module_id,
                "module_source": module.get("module_source"),
                "module_features": ";".join(f"{pos}:{aa}" for pos, aa in module["features"]),
                "module_overlap_count": int(len(local)),
                "edge_energy": edge_energy,
                "plaquette_energy": plaq_energy,
            }
        )
    summary = {
        "graphlet_id": graphlet_id,
        "graphlet_split": graphlet_split,
        "candidate_model": candidate_model,
        "parameter_set_id": parameter_set_id,
        "edge_innovation_total_energy": float(sum(edge_by_bit.values())),
        "plaquette_innovation_total_energy": float(sum(plaq_by_pair.values())),
        **concentration_summary(dict(edge_energy_by_module), prefix="edge_module"),
        **concentration_summary(dict(plaq_energy_by_module), prefix="plaquette_module"),
    }
    return rows, summary


def fit_polynomial_readout(
    base_values: list[float],
    empirical_values: list[float],
    *,
    degree: int,
) -> np.ndarray:
    x = np.asarray(base_values, dtype=np.float64)
    y = np.asarray(empirical_values, dtype=np.float64)
    if float(np.std(x)) <= 1e-12:
        return np.asarray([float(np.mean(y))], dtype=np.float64)
    x0 = (x - float(np.mean(x))) / max(float(np.std(x)), 1e-12)
    columns = [np.ones_like(x0)]
    for power in range(1, int(degree) + 1):
        columns.append(x0**power)
    design = np.vstack(columns).T
    beta, *_rest = np.linalg.lstsq(design, y, rcond=None)
    return beta.astype(np.float64)


def apply_polynomial_readout(base_values: np.ndarray, beta: np.ndarray) -> np.ndarray:
    x = np.asarray(base_values, dtype=np.float64)
    if beta.size <= 1 or float(np.std(x)) <= 1e-12:
        return np.full_like(x, float(beta[0]) if beta.size else float(np.mean(x)), dtype=np.float64)
    x0 = (x - float(np.mean(x))) / max(float(np.std(x)), 1e-12)
    columns = [np.ones_like(x0)]
    for power in range(1, int(beta.size)):
        columns.append(x0**power)
    return np.vstack(columns).T @ beta


def aggregate_module_energy(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(["candidate_model", "parameter_set_id", "graphlet_split"], dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        edge = sub.groupby("module_id")["edge_energy"].sum().to_dict()
        plaq = sub.groupby("module_id")["plaquette_energy"].sum().to_dict()
        out.append(
            {
                "candidate_model": keys[0],
                "parameter_set_id": keys[1],
                "graphlet_split": keys[2],
                "n_graphlet_module_rows": int(len(sub)),
                **concentration_summary(edge, prefix="edge_module"),
                **concentration_summary(plaq, prefix="plaquette_module"),
            }
        )
    return out
