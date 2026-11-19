#!/usr/bin/env python3
"""Fast PHOT fit for the bounded trait-space epistasis model.

This runner tests the minimal model patch in AdaptEnv: epistatic terms may be
bounded before they perturb latent traits. It rebuilds PHOT landscapes with a
small numerical search instead of running broad SMC-ABC, then compares graphlet
summaries against empirical PHOT and the previously computed full-ABC PHOT
comparator.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_ROOT, ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from phot_coupled_family_discriminating_abc import (  # noqa: E402
    ASSAY_ID,
    build_graphlets,
    graphlet_catalog_rows,
    guardrail_penalty,
    module_catalog_rows,
    record_metric_values,
)
from phot_coupled_landscape_grid_search import (  # noqa: E402
    load_config,
    summarize_metric_rows,
    write_csv,
    write_json,
)
from phot_nested_fit_protocol_benchmark import (  # noqa: E402
    apply_affine_readout,
    build_equivalence_rows,
    epistasis_property_metrics,
    fit_positive_affine_readout,
    identity_readout,
    observation_calibration_config,
    observation_calibration_enabled,
    prefixed_metrics,
    target_metric_rows,
)
from project_adapt_env.coupled_landscape import (  # noqa: E402
    additive_walsh_coefficients,
    low_order_auxiliary_core,
    soft_min_values,
    walsh_synthesize,
    zscore,
)
from project_adapt_env.dense_descriptive import (  # noqa: E402
    accessible_path_metrics,
    local_extrema_metrics,
)
from project_adapt_env.smc_abc import ParameterSpec  # noqa: E402
from proteingym_system_specific_current_graphlet_refit import (  # noqa: E402
    build_landscape,
    evaluate_landscape_on_graphlet,
    metric_distance,
)
import tem1_aggregate_graphlet_residual_ruggedness as graphlet_runner  # noqa: E402
from adaptenv.calibration import (  # noqa: E402
    _apply_observation_model,
    _fit_monotone_piecewise_observation_model,
)


EQUIVALENCE_COLUMNS = [
    "fit_method",
    "candidate_family",
    "candidate_model",
    "parameter_set_id",
    "reference",
    "reference_candidate_family",
    "reference_parameter_set_id",
    "candidate_test_distance",
    "reference_test_distance",
    "delta_test_distance",
    "relative_delta_test_distance",
    "max_component_distance",
    "within_absolute_margin",
    "within_relative_margin",
    "no_component_distance_failure",
    "operationally_equivalent",
    "fix_stage",
    "fixes_applied",
]

SOFTMIN_RAW_KEYS = {
    "softmin_tau",
    "softmin_n_traits",
    "softmin_auxiliary_scale",
    "softmin_auxiliary_bias",
    "softmin_auxiliary_correlation",
    "softmin_auxiliary_order",
    "softmin_output_scale",
    "softmin_structured_residual_scale",
}

SPARSE_WALSH_RAW_KEYS = {
    "sparse_walsh_epi_weight",
    "sparse_walsh_terms",
    "sparse_walsh_order2_weight",
    "sparse_walsh_order3_weight",
    "sparse_walsh_order4_weight",
    "sparse_walsh_order5_weight",
    "sparse_walsh_order6_weight",
    "sparse_walsh_order7_weight",
    "sparse_walsh_order8_weight",
}

NON_METRIC_COLUMNS = {
    "graphlet_id",
    "assay_id",
    "graphlet_split",
    "split",
    "graphlet_dimension",
    "mutation_features",
    "score_column",
    "candidate_model",
    "candidate_family",
    "parameter_set_id",
    "local_maxima",
    "local_minima",
}


def latent_generator_label(config: dict[str, Any]) -> str:
    if sparse_walsh_enabled(config):
        return "sparse_walsh_potts_field"
    if softmin_latent_enabled(config):
        return "softmin_two_additive_traits"
    return "current_trait_space"


def candidate_family_label(config: dict[str, Any]) -> str:
    return "sparse_walsh_field" if sparse_walsh_enabled(config) else "bounded_trait_space"


def candidate_model_label(config: dict[str, Any]) -> str:
    return (
        "additive_core_sparse_walsh_potts_field"
        if sparse_walsh_enabled(config)
        else "bounded_trait_space_current_model"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def parameter_specs(config: dict[str, Any]) -> list[ParameterSpec]:
    raw_specs = list(config["bounded_trait_space_fast_fit"]["parameter_specs"])
    specs: list[ParameterSpec] = []
    for raw in raw_specs:
        specs.append(
            ParameterSpec(
                name=str(raw["name"]),
                kind=str(raw.get("kind", "float")),
                low=float(raw["low"]),
                high=float(raw["high"]),
                initial_scale_fraction=float(raw.get("initial_scale_fraction", 0.25)),
                min_scale_fraction=float(raw.get("min_scale_fraction", 0.02)),
            )
        )
    return specs


def raw_key(raw: dict[str, float | int]) -> str:
    return json.dumps(raw, sort_keys=True, separators=(",", ":"))


def write_csv_with_columns(
    path: Path,
    rows: list[dict[str, Any]],
    *,
    columns: list[str],
) -> None:
    if rows:
        write_csv(path, rows)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=columns).to_csv(path, index=False)


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if key in out and isinstance(out[key], dict) and isinstance(value, dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def stage_configs(config: dict[str, Any]) -> list[dict[str, Any]]:
    stages = config.get("cumulative_fix_stages")
    if not stages:
        scoped = copy.deepcopy(config)
        scoped["_fix_stage"] = {
            "label": "baseline",
            "fixes": "",
            "description": "single-stage bounded trait-space fast fit",
            "index": 1,
        }
        return [scoped]

    scoped_stages: list[dict[str, Any]] = []
    current = copy.deepcopy(config)
    current.pop("cumulative_fix_stages", None)
    for index, raw_stage in enumerate(stages, start=1):
        label = str(raw_stage["label"])
        overrides = {
            key: value
            for key, value in dict(raw_stage).items()
            if key not in {"label", "fixes", "description", "overrides"}
        }
        overrides = deep_merge(overrides, dict(raw_stage.get("overrides", {})))
        scoped = deep_merge(current, overrides)
        scoped["_fix_stage"] = {
            "label": label,
            "fixes": str(raw_stage.get("fixes", label)),
            "description": str(raw_stage.get("description", "")),
            "index": index,
        }
        scoped_stages.append(scoped)
        current = copy.deepcopy(scoped)
        current.pop("_fix_stage", None)
    return scoped_stages


def fix_stage(config: dict[str, Any]) -> dict[str, Any]:
    return dict(config.get("_fix_stage", {}))


def stage_label(config: dict[str, Any]) -> str:
    return str(fix_stage(config).get("label", "baseline"))


def stage_prefix(config: dict[str, Any], fit_method: str) -> str:
    label = stage_label(config)
    if label == "baseline":
        return str(fit_method)
    return f"{label}__{fit_method}"


def readout_fit_arrays(
    *,
    values_by_graphlet: dict[str, np.ndarray],
    graphlets: list[dict[str, Any]],
    primary_score: str,
) -> tuple[np.ndarray, np.ndarray]:
    model_values: list[np.ndarray] = []
    empirical_values: list[np.ndarray] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        model_values.append(np.asarray(values_by_graphlet[str(graphlet["graphlet_id"])], dtype=np.float64))
        empirical_values.append(np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64))
    if not model_values or not empirical_values:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)
    return (
        np.concatenate([values.ravel() for values in model_values]),
        np.concatenate([values.ravel() for values in empirical_values]),
    )


def readout_fit_metrics(
    fitted: np.ndarray,
    observed: np.ndarray,
) -> dict[str, float]:
    y_hat = np.asarray(fitted, dtype=np.float64)
    y = np.asarray(observed, dtype=np.float64)
    mask = np.isfinite(y_hat) & np.isfinite(y)
    y_hat = y_hat[mask]
    y = y[mask]
    if y.size == 0:
        return {
            "readout_fit_n": 0,
            "readout_fit_rmse": math.nan,
            "readout_fit_correlation": math.nan,
        }
    rmse = float(np.sqrt(np.mean((y_hat - y) ** 2)))
    corr = (
        float(np.corrcoef(y_hat, y)[0, 1])
        if y.size >= 2 and np.std(y_hat) > 0 and np.std(y) > 0
        else math.nan
    )
    return {
        "readout_fit_n": int(y.size),
        "readout_fit_rmse": rmse,
        "readout_fit_correlation": corr,
    }


def sigmoid(values: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(values, dtype=np.float64), -60.0, 60.0)))


def fit_monotone_sigmoid_readout(
    x: np.ndarray,
    y: np.ndarray,
    *,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Fit a two-parameter monotone sigmoid with empirical-range output fixed."""
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    ok = np.isfinite(x_arr) & np.isfinite(y_arr)
    x_arr = x_arr[ok]
    y_arr = y_arr[ok]
    if x_arr.size == 0 or y_arr.size == 0:
        return {
            "readout_method": "monotone_sigmoid",
            "readout_intercept": 0.0,
            "readout_slope": 1.0,
            "readout_midpoint": 0.0,
            "readout_sigmoid_slope": 1.0,
            "readout_y_min": 0.0,
            "readout_y_max": 1.0,
            "readout_fit_n": 0,
            "readout_fit_rmse": math.nan,
            "readout_fit_correlation": math.nan,
        }

    q_low = float(config.get("y_quantile_low", 0.02))
    q_high = float(config.get("y_quantile_high", 0.98))
    y_min = float(np.quantile(y_arr, q_low))
    y_max = float(np.quantile(y_arr, q_high))
    if y_max <= y_min + 1e-12:
        y_min = float(np.min(y_arr))
        y_max = float(np.max(y_arr))
    if y_max <= y_min + 1e-12:
        fitted = np.full_like(y_arr, float(np.mean(y_arr)), dtype=np.float64)
        return {
            "readout_method": "monotone_sigmoid",
            "readout_intercept": float(np.mean(y_arr)),
            "readout_slope": 0.0,
            "readout_midpoint": float(np.mean(x_arr)),
            "readout_sigmoid_slope": 1.0,
            "readout_y_min": float(np.mean(y_arr)),
            "readout_y_max": float(np.mean(y_arr)),
            **readout_fit_metrics(fitted, y_arr),
        }

    x_sd = max(float(np.std(x_arr)), 1e-8)
    n_midpoints = int(config.get("n_midpoint_grid", 17))
    midpoint_quantiles = np.linspace(0.05, 0.95, max(n_midpoints, 3))
    midpoints = np.unique(
        np.concatenate(
            [
                np.quantile(x_arr, midpoint_quantiles),
                np.asarray([float(np.mean(x_arr)), float(np.median(x_arr))], dtype=np.float64),
            ]
        )
    )
    slope_multipliers = config.get(
        "slope_multipliers",
        [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0],
    )
    slopes = np.asarray([float(value) / x_sd for value in slope_multipliers], dtype=np.float64)
    slopes = slopes[np.isfinite(slopes) & (slopes > 0.0)]
    if slopes.size == 0:
        slopes = np.asarray([1.0 / x_sd], dtype=np.float64)

    best_rmse = math.inf
    best_midpoint = float(np.median(x_arr))
    best_slope = float(slopes[0])
    best_fitted = y_min + (y_max - y_min) * sigmoid(best_slope * (x_arr - best_midpoint))
    for midpoint in midpoints:
        for slope in slopes:
            fitted = y_min + (y_max - y_min) * sigmoid(float(slope) * (x_arr - float(midpoint)))
            rmse = float(np.sqrt(np.mean((fitted - y_arr) ** 2)))
            if rmse < best_rmse:
                best_rmse = rmse
                best_midpoint = float(midpoint)
                best_slope = float(slope)
                best_fitted = fitted

    return {
        "readout_method": "monotone_sigmoid",
        "readout_intercept": float(y_min),
        "readout_slope": float(y_max - y_min),
        "readout_midpoint": best_midpoint,
        "readout_sigmoid_slope": best_slope,
        "readout_y_min": float(y_min),
        "readout_y_max": float(y_max),
        **readout_fit_metrics(best_fitted, y_arr),
    }


def fit_observation_readout(
    *,
    values_by_graphlet: dict[str, np.ndarray],
    graphlets: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    if not observation_calibration_enabled(config):
        return identity_readout()
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    calibration_cfg = observation_calibration_config(config)
    method = str(calibration_cfg.get("method", "positive_affine")).strip().lower()
    x, y = readout_fit_arrays(
        values_by_graphlet=values_by_graphlet,
        graphlets=graphlets,
        primary_score=primary_score,
    )
    if method in {"identity", "none"}:
        return identity_readout()
    if method in {"positive_affine", "affine"}:
        min_slope = float(calibration_cfg.get("min_slope", 1e-12))
        return fit_positive_affine_readout([x], [y], min_slope=min_slope)
    if method == "monotone_piecewise":
        model = _fit_monotone_piecewise_observation_model(
            x,
            y,
            max_knots=int(calibration_cfg.get("max_knots", 16)),
            min_bin_size=int(calibration_cfg.get("min_bin_size", 8)),
        )
        fitted = _apply_observation_model(x, model)
        return {
            "readout_method": "monotone_piecewise",
            "readout_intercept": 0.0,
            "readout_slope": 1.0,
            "readout_x_knots": json.dumps(model.get("x_knots", []), separators=(",", ":")),
            "readout_y_knots": json.dumps(model.get("y_knots", []), separators=(",", ":")),
            **readout_fit_metrics(fitted, y),
        }
    if method in {"monotone_sigmoid", "sigmoid", "sigmoid_2param"}:
        return fit_monotone_sigmoid_readout(x, y, config=calibration_cfg)
    raise ValueError(f"Unsupported observation_calibration.method: {method}")


def apply_observation_readout(values: np.ndarray, readout: dict[str, Any]) -> np.ndarray:
    method = str(readout.get("readout_method", "identity"))
    if method == "monotone_piecewise":
        return _apply_observation_model(
            values,
            {
                "kind": "monotone_piecewise",
                "x_knots": json.loads(str(readout.get("readout_x_knots", "[]"))),
                "y_knots": json.loads(str(readout.get("readout_y_knots", "[]"))),
            },
        )
    if method == "monotone_sigmoid":
        return float(readout.get("readout_y_min", readout.get("readout_intercept", 0.0))) + float(
            readout.get("readout_slope", 1.0)
        ) * sigmoid(
            float(readout.get("readout_sigmoid_slope", 1.0))
            * (np.asarray(values, dtype=np.float64) - float(readout.get("readout_midpoint", 0.0)))
        )
    if method == "identity":
        return np.asarray(values, dtype=np.float64)
    return apply_affine_readout(values, readout)


def prefixed_numeric_metrics(metrics: dict[str, Any], prefix: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in metrics.items():
        if key in NON_METRIC_COLUMNS:
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        out[f"{prefix}{key}"] = numeric
    return out


def raw_diagnostic_metrics(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    metrics: dict[str, Any] = {}
    metrics.update(epistasis_property_metrics(arr, dim=int(dim)))
    metrics.update(local_extrema_metrics(arr, dim=int(dim)))
    metrics.update(accessible_path_metrics(arr, dim=int(dim)))
    return prefixed_numeric_metrics(metrics, "")


def sample_raw(specs: list[ParameterSpec], rng: np.random.RandomState) -> dict[str, float | int]:
    return {spec.name: spec.sample_prior(rng) for spec in specs}


def clamp_raw(raw: dict[str, float | int], specs: list[ParameterSpec]) -> dict[str, float | int]:
    out: dict[str, float | int] = {}
    by_name = {spec.name: spec for spec in specs}
    for name, value in raw.items():
        spec = by_name.get(name)
        if spec is None:
            continue
        clipped = float(np.clip(float(value), float(spec.low), float(spec.high)))
        out[name] = int(round(clipped)) if spec.kind == "int" else clipped
    return out


def center_raw_from_updates(updates: dict[str, Any], specs: list[ParameterSpec]) -> dict[str, float | int]:
    raw: dict[str, float | int] = {}
    for spec in specs:
        if spec.name in updates:
            value = float(updates[spec.name])
        elif spec.name == "trait_space_epistasis_bound":
            value = float((float(spec.low) * float(spec.high)) ** 0.5)
        else:
            value = 0.5 * (float(spec.low) + float(spec.high))
        raw[spec.name] = int(round(value)) if spec.kind == "int" else value
    return clamp_raw(raw, specs)


def jitter_raw(
    center: dict[str, float | int],
    specs: list[ParameterSpec],
    *,
    width_fraction: float,
    rng: np.random.RandomState,
) -> dict[str, float | int]:
    raw: dict[str, float | int] = {}
    for spec in specs:
        support = float(spec.high - spec.low)
        width = max(float(width_fraction) * support, 1.0 if spec.kind == "int" else support * 0.01)
        value = float(center[spec.name]) + rng.normal(0.0, width * 0.5)
        raw[spec.name] = int(round(value)) if spec.kind == "int" else value
    return clamp_raw(raw, specs)


def unique_raws(raws: list[dict[str, float | int]]) -> list[dict[str, float | int]]:
    seen: set[str] = set()
    out: list[dict[str, float | int]] = []
    for raw in raws:
        key = raw_key(raw)
        if key in seen:
            continue
        seen.add(key)
        out.append(raw)
    return out


def softmin_latent_config(config: dict[str, Any]) -> dict[str, Any]:
    return dict(config["bounded_trait_space_fast_fit"].get("softmin_latent", {}))


def softmin_latent_enabled(config: dict[str, Any]) -> bool:
    return bool(softmin_latent_config(config).get("enabled", False))


def sparse_walsh_config(config: dict[str, Any]) -> dict[str, Any]:
    return dict(config["bounded_trait_space_fast_fit"].get("sparse_walsh_field", {}))


def sparse_walsh_enabled(config: dict[str, Any]) -> bool:
    return bool(sparse_walsh_config(config).get("enabled", False))


def sparse_walsh_orders(config: dict[str, Any]) -> list[int]:
    cfg = sparse_walsh_config(config)
    min_order = max(int(cfg.get("min_order", 2)), 1)
    max_order = max(int(cfg.get("max_order", min_order)), min_order)
    return list(range(min_order, max_order + 1))


def sparse_walsh_order_weights(
    *,
    raw: dict[str, float | int],
    config: dict[str, Any],
) -> list[float]:
    cfg = sparse_walsh_config(config)
    orders = sparse_walsh_orders(config)
    configured = list(cfg.get("order_weights", []))
    weights: list[float] = []
    for index, order in enumerate(orders):
        key = f"sparse_walsh_order{order}_weight"
        if key in raw:
            value = float(raw[key])
        elif index < len(configured):
            value = float(configured[index])
        else:
            value = 1.0
        weights.append(max(value, 0.0))
    if sum(weights) <= 0.0:
        weights = [1.0 for _order in orders]
    total = float(sum(weights))
    return [float(value / total) for value in weights]


def sparse_walsh_cluster_positions(
    *,
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    config: dict[str, Any],
) -> list[list[int]]:
    cfg = sparse_walsh_config(config)
    if cfg.get("cluster_positions"):
        return [
            sorted({int(pos) for pos in cluster})
            for cluster in cfg.get("cluster_positions", [])
            if len({int(pos) for pos in cluster}) >= 2
        ]
    source = str(cfg.get("support_source", "empirical_calibration_plaquette"))
    clusters: list[list[int]] = []
    seen: set[tuple[int, ...]] = set()
    for module in modules_by_source.get(source, []):
        positions = tuple(sorted({int(feature[0]) for feature in module.get("features", [])}))
        if len(positions) < 2 or positions in seen:
            continue
        clusters.append(list(positions))
        seen.add(positions)
    if bool(cfg.get("include_graphlet_feature_clusters", True)):
        splits = {str(value) for value in cfg.get("support_graphlet_splits", ["calibration", "test"])}
        for graphlet in graphlets:
            if str(graphlet.get("graphlet_split", "")) not in splits:
                continue
            positions = tuple(sorted({int(feature[0]) for feature in graphlet["features"]}))
            if len(positions) < 2 or positions in seen:
                continue
            clusters.append(list(positions))
            seen.add(positions)
    return clusters


def inject_sparse_walsh_cluster_positions(
    *,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    if not sparse_walsh_enabled(config):
        return config
    scoped = copy.deepcopy(config)
    cfg = dict(scoped["bounded_trait_space_fast_fit"].get("sparse_walsh_field", {}))
    cfg["cluster_positions"] = sparse_walsh_cluster_positions(
        graphlets=graphlets,
        modules_by_source=modules_by_source,
        config=scoped,
    )
    scoped["bounded_trait_space_fast_fit"]["sparse_walsh_field"] = cfg
    return scoped


def additive_projection(values: np.ndarray, *, dim: int) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    coeffs = {
        1 << bit: float(coeff)
        for bit, coeff in enumerate(additive_walsh_coefficients(arr, dim=int(dim)))
    }
    return float(np.mean(arr)) + walsh_synthesize(coeffs, dim=int(dim))


def binary_design_matrix(dim: int) -> np.ndarray:
    rows: list[list[float]] = []
    for node in range(1 << int(dim)):
        rows.append([1.0] + [1.0 if node & (1 << bit) else -1.0 for bit in range(int(dim))])
    return np.asarray(rows, dtype=np.float64)


def random_walsh_core(
    *,
    dim: int,
    rng: np.random.RandomState,
    base_norm: np.ndarray,
    min_order: int = 1,
    max_order: int = 1,
    n_terms: int | None = None,
) -> np.ndarray:
    masks: list[int] = []
    for order in range(max(int(min_order), 1), min(int(max_order), int(dim)) + 1):
        for combo in itertools.combinations(range(int(dim)), order):
            mask = 0
            for bit in combo:
                mask |= 1 << bit
            masks.append(mask)
    if not masks:
        return low_order_auxiliary_core(dim=int(dim), rng=rng, base_norm=base_norm)
    if n_terms is not None and int(n_terms) > 0 and int(n_terms) < len(masks):
        masks = list(rng.choice(np.asarray(masks, dtype=np.int64), size=int(n_terms), replace=False))
    coeffs = {int(mask): float(rng.normal()) for mask in masks}
    raw = zscore(walsh_synthesize(coeffs, dim=int(dim)))
    base = zscore(base_norm)
    corr = float(np.corrcoef(raw, base)[0, 1]) if np.std(raw) > 0.0 and np.std(base) > 0.0 else math.nan
    if np.isfinite(corr):
        raw = raw - corr * base
    return zscore(raw)


def soft_min_many_values(traits: list[np.ndarray], *, tau: float) -> np.ndarray:
    stacked = np.vstack([np.asarray(trait, dtype=np.float64) for trait in traits])
    temp = max(float(tau), 1e-6)
    m = np.min(stacked, axis=0)
    return m - temp * np.log(np.sum(np.exp(-(stacked - m) / temp), axis=0))


def softmin_latent_values(
    *,
    base_values: np.ndarray,
    graphlet: dict[str, Any],
    raw: dict[str, float | int],
    config: dict[str, Any],
) -> np.ndarray:
    cfg = softmin_latent_config(config)
    dim = int(graphlet["graphlet_dimension"])
    base_arr = np.asarray(base_values, dtype=np.float64)
    base_trait = additive_projection(base_arr, dim=dim) if bool(cfg.get("project_base_additive", True)) else base_arr
    base_sd = max(float(np.std(base_trait)), 1e-8)
    base_norm = zscore(base_trait)

    n_traits = max(int(raw.get("softmin_n_traits", cfg.get("n_traits", 2))), 2)
    correlation = float(raw.get("softmin_auxiliary_correlation", cfg.get("auxiliary_correlation", 0.25)))
    correlation = float(np.clip(correlation, -0.99, 0.99))
    auxiliary_scale = float(raw.get("softmin_auxiliary_scale", cfg.get("auxiliary_scale", 1.0)))
    auxiliary_bias = float(raw.get("softmin_auxiliary_bias", cfg.get("auxiliary_bias", 0.0)))
    auxiliary_order = int(raw.get("softmin_auxiliary_order", cfg.get("auxiliary_order", 1)))
    auxiliary_order = int(np.clip(auxiliary_order, 1, max(dim, 1)))

    traits = [base_trait]
    for trait_index in range(1, n_traits):
        seed = graphlet_runner.stable_seed(
            "path-forward-softmin-auxiliary",
            str(graphlet["graphlet_id"]),
            str(cfg.get("seed", 20260617)),
            str(trait_index),
            str(n_traits),
        )
        rng = np.random.RandomState(seed)
        if auxiliary_order <= 1:
            auxiliary_core = low_order_auxiliary_core(dim=dim, rng=rng, base_norm=base_norm)
        else:
            auxiliary_core = random_walsh_core(
                dim=dim,
                rng=rng,
                base_norm=base_norm,
                min_order=1,
                max_order=auxiliary_order,
                n_terms=cfg.get("auxiliary_n_terms"),
            )
        auxiliary_norm = zscore(
            correlation * base_norm + math.sqrt(max(1.0 - correlation * correlation, 0.0)) * auxiliary_core
        )
        trait_bias = auxiliary_bias + float(cfg.get("trait_bias_jitter", 0.0)) * float(rng.normal())
        trait_scale = auxiliary_scale * math.exp(float(cfg.get("trait_scale_jitter", 0.0)) * float(rng.normal()))
        auxiliary_trait = float(np.mean(base_trait)) + base_sd * (trait_bias + trait_scale * auxiliary_norm)
        traits.append(auxiliary_trait)

    tau_raw = float(raw.get("softmin_tau", cfg.get("tau", 0.5)))
    tau_units = str(cfg.get("tau_units", "base_sd")).strip().lower()
    tau = tau_raw * base_sd if tau_units == "base_sd" else tau_raw
    if n_traits == 2:
        values, _innovation, _weight = soft_min_values(
            base_values=traits[0],
            auxiliary_values=traits[1],
            tau=tau,
        )
    else:
        values = soft_min_many_values(traits, tau=tau)
    texture = zscore(values) if bool(cfg.get("normalize_output", False)) else values - float(np.mean(values))
    residual_cfg = dict(cfg.get("structured_residual", {}))
    if bool(residual_cfg.get("enabled", False)):
        residual_seed = graphlet_runner.stable_seed(
            "path-forward-softmin-structured-residual",
            str(graphlet["graphlet_id"]),
            str(cfg.get("seed", 20260617)),
        )
        residual_rng = np.random.RandomState(residual_seed)
        residual = random_walsh_core(
            dim=dim,
            rng=residual_rng,
            base_norm=texture,
            min_order=int(residual_cfg.get("min_order", 2)),
            max_order=int(residual_cfg.get("max_order", min(4, dim))),
            n_terms=residual_cfg.get("n_terms", max(4, dim)),
        )
        residual_scale = float(
            raw.get(
                "softmin_structured_residual_scale",
                residual_cfg.get("scale", 0.0),
            )
        )
        texture = texture + residual_scale * residual
        if bool(residual_cfg.get("renormalize_after_residual", True)):
            texture = zscore(texture)
    output_scale = float(raw.get("softmin_output_scale", cfg.get("output_scale", 1.0)))
    output_center = str(cfg.get("output_center", "base_mean")).strip().lower()
    if output_center == "zero":
        center = 0.0
    elif output_center == "base_median":
        center = float(np.median(base_trait))
    else:
        center = float(np.mean(base_trait))
    return center + output_scale * texture


def candidate_updates(
    *,
    base_updates: dict[str, Any],
    raw: dict[str, float | int],
    config: dict[str, Any],
) -> dict[str, Any]:
    raw_updates = dict(raw)
    global_midpoint = raw_updates.pop("global_readout_midpoint", None)
    global_slope = raw_updates.pop("global_readout_slope", None)
    global_weight = raw_updates.pop("global_readout_weight", None)
    ratio = raw_updates.pop("specific_epistasis_variance_ratio", None)
    for key in SOFTMIN_RAW_KEYS:
        raw_updates.pop(key, None)
    for key in SPARSE_WALSH_RAW_KEYS:
        raw_updates.pop(key, None)
    updates = dict(base_updates)
    updates.update(raw_updates)
    updates["trait_space_epistasis_bound_mode"] = str(
        config["bounded_trait_space_fast_fit"].get("bound_mode", "tanh")
    )
    updates["generic_epistasis_target"] = "stability"
    global_cfg = dict(config["bounded_trait_space_fast_fit"].get("global_readout", {}))
    if bool(global_cfg.get("enabled", False)):
        term = {
            "source": str(global_cfg.get("source", "stability_margin")),
            "transform": str(global_cfg.get("transform", "sigmoid")),
            "midpoint": float(
                global_midpoint
                if global_midpoint is not None
                else global_cfg.get("midpoint", 0.0)
            ),
            "slope": float(
                global_slope
                if global_slope is not None
                else global_cfg.get("slope", 1.0)
            ),
            "weight": float(
                global_weight
                if global_weight is not None
                else global_cfg.get("weight", 1.0)
            ),
            "exponent": float(global_cfg.get("exponent", 1.0)),
        }
        updates["observed_fitness_combine_mode"] = str(global_cfg.get("combine_mode", "sum"))
        updates["observed_fitness_terms"] = [term]
        if "fitness_output_scale" in global_cfg:
            updates["fitness_output_scale"] = float(global_cfg["fitness_output_scale"])
    ratio_cfg = dict(
        config["bounded_trait_space_fast_fit"].get("specific_epistasis_variance_ratio", {})
    )
    if ratio is not None:
        scale = float(ratio_cfg.get("strength_scale", 0.04))
        stability_scale = abs(float(updates.get("stability_scale", base_updates.get("stability_scale", 1.0))))
        term_count = max(float(updates.get("n_epistatic_terms", base_updates.get("n_epistatic_terms", 1))), 1.0)
        reference_terms = max(float(ratio_cfg.get("reference_n_terms", 220.0)), 1.0)
        strength = float(ratio) * scale * stability_scale * math.sqrt(reference_terms / term_count)
        max_strength = float(ratio_cfg.get("max_strength", 0.08))
        strength = float(np.clip(strength, 0.0, max_strength))
        updates["epistasis_strength"] = strength
        updates["empirical_pairwise_strength"] = strength * float(
            ratio_cfg.get("empirical_pairwise_multiplier", 1.0)
        )
    latent_cfg = softmin_latent_config(config)
    if bool(latent_cfg.get("enabled", False)) and bool(latent_cfg.get("disable_specific_epistasis", True)):
        updates["n_epistatic_terms"] = 0
        updates["epistasis_strength"] = 0.0
        updates["empirical_pairwise_strength"] = 0.0
    walsh_cfg = sparse_walsh_config(config)
    if bool(walsh_cfg.get("enabled", False)):
        updates["n_epistatic_terms"] = 0
        updates["n_higher_order_epistatic_terms"] = 0
        updates["empirical_pairwise_strength"] = 0.0
        updates["generic_epistasis_target"] = "stability"
        updates["epistasis_strength"] = float(
            raw.get("sparse_walsh_epi_weight", walsh_cfg.get("epi_weight", updates.get("epistasis_strength", 0.025)))
        )
        updates["sparse_walsh_terms"] = int(
            raw.get("sparse_walsh_terms", walsh_cfg.get("terms", 1024))
        )
        updates["sparse_walsh_min_order"] = int(walsh_cfg.get("min_order", 2))
        updates["sparse_walsh_max_order"] = int(walsh_cfg.get("max_order", 5))
        updates["sparse_walsh_order_weights"] = sparse_walsh_order_weights(
            raw=raw,
            config=config,
        )
        updates["sparse_walsh_cluster_positions"] = list(walsh_cfg.get("cluster_positions", []))
        updates["sparse_walsh_seed"] = int(walsh_cfg.get("seed", 20260617))
    return updates


def graphlet_values_for_updates(
    *,
    updates: dict[str, Any],
    raw: dict[str, float | int],
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    seed_offset: int = 0,
) -> dict[str, np.ndarray]:
    base_config = dict(config["adapt_env"]["base_config"])
    seed = int(base_config.get("seed", 42)) + int(seed_offset)
    landscape = build_landscape(
        base_config=base_config,
        updates=updates,
        wildtype=str(context["wildtype"]),
        alignment_path=Path(str(context["alignment_path"])),
        seed=seed,
    )
    values_by_graphlet: dict[str, np.ndarray] = {}
    for graphlet in graphlets:
        base_values = evaluate_landscape_on_graphlet(
            landscape,
            wildtype=str(context["wildtype"]),
            features=tuple(graphlet["features"]),
        )
        if softmin_latent_enabled(config):
            base_values = softmin_latent_values(
                base_values=base_values,
                graphlet=graphlet,
                raw=raw,
                config=config,
            )
        values_by_graphlet[str(graphlet["graphlet_id"])] = base_values
    return values_by_graphlet


def fit_readout_for_values(
    *,
    values_by_graphlet: dict[str, np.ndarray],
    graphlets: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    return fit_observation_readout(
        values_by_graphlet=values_by_graphlet,
        graphlets=graphlets,
        config=config,
    )


def finite_column_mean(rows: pd.DataFrame, column: str) -> float:
    if column not in rows:
        return math.nan
    values = pd.to_numeric(rows[column], errors="coerce").dropna()
    return float(values.mean()) if not values.empty else math.nan


def normalized_target_distance(
    *,
    model_value: float,
    target_rows: pd.DataFrame,
    target_metric: str,
    scale_floor_fraction: float = 0.10,
    absolute_scale_floor: float = 1e-6,
) -> float:
    target_value = finite_column_mean(target_rows, target_metric)
    target_values = pd.to_numeric(target_rows.get(target_metric), errors="coerce").dropna()
    if not np.isfinite(model_value) or not np.isfinite(target_value):
        return math.nan
    scale = float(target_values.std(ddof=0)) if len(target_values) >= 2 else 0.0
    scale = max(scale, abs(target_value) * float(scale_floor_fraction), float(absolute_scale_floor))
    return float(abs(float(model_value) - target_value) / scale)


def raw_latent_objective(
    *,
    active_rows: list[dict[str, Any]],
    target_rows: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[float, dict[str, Any]]:
    cfg = dict(config["bounded_trait_space_fast_fit"].get("raw_latent_objective", {}))
    if not bool(cfg.get("enabled", False)) or not active_rows:
        return 0.0, {
            "raw_latent_objective_distance": 0.0,
            "raw_latent_objective_penalty": 0.0,
            "raw_latent_objective_weight": 0.0,
        }
    terms = [str(term) for term in cfg.get("terms", [])]
    metric_weights = {str(key): float(value) for key, value in dict(cfg.get("metric_weights", {})).items()}
    active_df = pd.DataFrame(active_rows)
    pieces: list[float] = []
    weights: list[float] = []
    extras: dict[str, Any] = {}
    for term in terms:
        model_metric = f"raw_{term}"
        model_value = finite_column_mean(active_df, model_metric)
        distance = normalized_target_distance(
            model_value=model_value,
            target_rows=target_rows,
            target_metric=term,
            scale_floor_fraction=float(cfg.get("scale_floor_fraction", 0.10)),
            absolute_scale_floor=float(cfg.get("absolute_scale_floor", 1e-6)),
        )
        extras[f"raw_latent_distance_component__{term}"] = distance
        if not np.isfinite(distance):
            continue
        weight = max(float(metric_weights.get(term, 1.0)), 0.0)
        pieces.append(weight * float(distance))
        weights.append(weight)
    raw_distance = float(np.sum(pieces) / max(np.sum(weights), 1e-12)) if pieces else math.nan
    objective_weight = float(cfg.get("objective_weight", 1.0))
    penalty = float(objective_weight * raw_distance) if np.isfinite(raw_distance) else 0.0
    extras.update(
        {
            "raw_latent_objective_distance": raw_distance,
            "raw_latent_objective_penalty": penalty,
            "raw_latent_objective_weight": objective_weight,
        }
    )
    return penalty, extras


def record_candidate_metrics(
    *,
    values_by_graphlet: dict[str, np.ndarray],
    graphlets: list[dict[str, Any]],
    config: dict[str, Any],
    readout: dict[str, Any],
    fit_method: str,
    parameter_set_id: str,
    raw: dict[str, float | int],
    target_rows: pd.DataFrame,
) -> tuple[float, dict[str, Any]]:
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    metric_cols = [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]
    calibration_cfg = observation_calibration_config(config)
    apply_to_objective = bool(
        calibration_cfg.get("apply_to_objective", observation_calibration_enabled(config))
    )
    random_walk_cfg = {"analysis": config.get("analysis", {})}
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    active_rows: list[dict[str, Any]] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        raw_values = np.asarray(values_by_graphlet[str(graphlet["graphlet_id"])], dtype=np.float64)
        values = apply_observation_readout(raw_values, readout) if apply_to_objective else raw_values
        stage = fix_stage(config)
        extra = {
            **readout,
            "fit_method": fit_method,
            "candidate_family": candidate_family_label(config),
            "source_parameter_set_id": parameter_set_id,
            "fix_stage": stage.get("label", "baseline"),
            "fix_stage_index": stage.get("index", 1),
            "fixes_applied": stage.get("fixes", ""),
            "latent_generator": latent_generator_label(config),
            "readout_applied_to_objective": bool(apply_to_objective),
            **{f"raw__{key}": value for key, value in raw.items()},
        }
        metric = record_metric_values(
            graphlet=graphlet,
            candidate_model=candidate_model_label(config),
            candidate_family=candidate_family_label(config),
            parameter_set_id=parameter_set_id,
            values=values,
            score_column=primary_score,
            metric_rows=metric_rows,
            walsh_rows=walsh_rows,
            gft_rows=gft_rows,
            random_walk_cfg=random_walk_cfg,
            extra=extra,
        )
        epi = epistasis_property_metrics(values, dim=int(graphlet["graphlet_dimension"]))
        raw_epi = epistasis_property_metrics(raw_values, dim=int(graphlet["graphlet_dimension"]))
        raw_diag = raw_diagnostic_metrics(raw_values, dim=int(graphlet["graphlet_dimension"]))
        metric.update(epi)
        metric.update(prefixed_metrics(raw_epi, "raw_"))
        metric.update(prefixed_numeric_metrics(raw_diag, "raw_"))
        active_rows.append({**epi, **prefixed_metrics(raw_epi, "raw_"), **prefixed_numeric_metrics(raw_diag, "raw_")})

    model_rows = pd.DataFrame(metric_rows)
    distance, by_metric = metric_distance(model_rows, target_rows, metric_cols)
    penalty, guardrail = guardrail_penalty(by_metric, config)
    raw_penalty, raw_objective = raw_latent_objective(
        active_rows=active_rows,
        target_rows=target_rows,
        config=config,
    )
    if not np.isfinite(distance):
        distance = float("inf")
    objective = float(distance + penalty + raw_penalty)
    extras: dict[str, Any] = {
        **readout,
        "calibration_graphlet_distance": float(distance),
        "objective_distance": objective,
        **guardrail,
        **raw_objective,
        **{f"distance_component__{key}": value for key, value in by_metric.items()},
    }
    if active_rows:
        active_df = pd.DataFrame(active_rows)
        for col in [
            "magnitude_epistasis_mean",
            "raw_magnitude_epistasis_mean",
            "signed_pairwise_epistasis_mean",
            "raw_signed_pairwise_epistasis_mean",
            "raw_single_mutant_effect_sd",
            "raw_n_local_maxima",
            "raw_accessible_path_fraction",
            "sign_epistasis_fraction_local",
            "reciprocal_sign_epistasis_fraction_local",
            "diminishing_returns_epistasis_fraction",
            "single_mutant_effect_sd",
            "magnitude_epistasis_to_single_mutant_sd_ratio",
            "signed_pairwise_epistasis_sd",
            "signed_pairwise_epistasis_q10",
            "signed_pairwise_epistasis_q90",
        ]:
            if col in active_df:
                values = pd.to_numeric(active_df[col], errors="coerce").dropna()
                if not values.empty:
                    extras[f"{col}_mean"] = float(values.mean())
    return objective, extras


def evaluate_raw_candidate(
    *,
    raw: dict[str, float | int],
    fit_method: str,
    evaluation_index: int,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    base_updates: dict[str, Any],
    target_rows: pd.DataFrame,
) -> dict[str, Any]:
    stage = fix_stage(config)
    effective_fit_method = stage_prefix(config, fit_method)
    parameter_set_id = (
        f"{candidate_family_label(config)}__{stage_label(config)}__{fit_method}__"
        f"{graphlet_runner.stable_seed(stage_label(config), fit_method, raw_key(raw)):010d}"
    )
    values_by_graphlet = graphlet_values_for_updates(
        updates=candidate_updates(base_updates=base_updates, raw=raw, config=config),
        raw=raw,
        config=config,
        graphlets=graphlets,
        context=context,
        seed_offset=1000 * int(stage.get("index", 1)),
    )
    readout = fit_readout_for_values(
        values_by_graphlet=values_by_graphlet,
        graphlets=graphlets,
        config=config,
    )
    objective, extras = record_candidate_metrics(
        values_by_graphlet=values_by_graphlet,
        graphlets=graphlets,
        config=config,
        readout=readout,
        fit_method=effective_fit_method,
        parameter_set_id=parameter_set_id,
        raw=raw,
        target_rows=target_rows,
    )
    return {
        "fit_method": effective_fit_method,
        "candidate_family": candidate_family_label(config),
        "candidate_model": candidate_model_label(config),
        "parameter_set_id": parameter_set_id,
        "fix_stage": stage.get("label", "baseline"),
        "fix_stage_index": stage.get("index", 1),
        "fixes_applied": stage.get("fixes", ""),
        "fix_stage_description": stage.get("description", ""),
        "latent_generator": latent_generator_label(config),
        "evaluation_index": int(evaluation_index),
        "objective_distance": float(objective),
        **{f"raw__{key}": value for key, value in raw.items()},
        **extras,
    }


def run_fast_fit(
    *,
    config: dict[str, Any],
    output_dir: Path,
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    target_rows: pd.DataFrame,
) -> list[dict[str, Any]]:
    cfg = dict(config["bounded_trait_space_fast_fit"])
    specs = parameter_specs(config)
    base_updates = dict(context["fit_summary"]["branches"]["smc_abc_best_raw"]["fitted_parameters"])
    rng = np.random.RandomState(int(cfg.get("seed", 20260616)))
    center = center_raw_from_updates(base_updates, specs)
    candidates: list[tuple[str, dict[str, float | int]]] = []
    candidates.append(("current_best_center", center))
    for bound in cfg.get("initial_bounds", [0.05, 0.1, 0.25, 0.5, 1.0]):
        bounded = dict(center)
        bounded["trait_space_epistasis_bound"] = float(bound)
        candidates.append((f"center_bound_{float(bound):.3g}", clamp_raw(bounded, specs)))
    for index, raw_override in enumerate(cfg.get("initial_candidates", []), start=1):
        seeded = dict(center)
        seeded.update({str(key): value for key, value in dict(raw_override).items()})
        candidates.append((f"seeded_candidate_{index:02d}", clamp_raw(seeded, specs)))
    for _ in range(int(cfg.get("random_candidates", 64))):
        candidates.append(("random_prior", sample_raw(specs, rng)))

    rows: list[dict[str, Any]] = []
    for index, raw in enumerate(unique_raws([raw for _fit_method, raw in candidates]), start=1):
        matching_method = next(method for method, candidate in candidates if raw_key(candidate) == raw_key(raw))
        rows.append(
            evaluate_raw_candidate(
                raw=raw,
                fit_method=matching_method,
                evaluation_index=index,
                config=config,
                graphlets=graphlets,
                context=context,
                base_updates=base_updates,
                target_rows=target_rows,
            )
        )
        write_csv(output_dir / "bounded_gate_candidates.csv", rows)

    top_n = max(int(cfg.get("local_top_n", 4)), 0)
    local_rounds = max(int(cfg.get("local_rounds", 1)), 0)
    local_candidates = max(int(cfg.get("local_candidates_per_round", 24)), 0)
    width = float(cfg.get("local_width_fraction", 0.20))
    for local_round in range(local_rounds):
        finite_rows = sorted(
            [row for row in rows if np.isfinite(float(row["objective_distance"]))],
            key=lambda row: float(row["objective_distance"]),
        )
        if not finite_rows or top_n <= 0 or local_candidates <= 0:
            break
        seeds = finite_rows[:top_n]
        new_raws: list[dict[str, float | int]] = []
        for seed_row in seeds:
            center_raw = {
                key.replace("raw__", "", 1): seed_row[key]
                for key in seed_row
                if str(key).startswith("raw__")
            }
            for _ in range(local_candidates):
                new_raws.append(
                    jitter_raw(
                        center_raw,
                        specs,
                        width_fraction=width / float(local_round + 1),
                        rng=rng,
                    )
                )
        existing = {raw_key({key.replace("raw__", "", 1): row[key] for key in row if str(key).startswith("raw__")}) for row in rows}
        for raw in unique_raws(new_raws):
            if raw_key(raw) in existing:
                continue
            rows.append(
                evaluate_raw_candidate(
                    raw=raw,
                    fit_method=f"local_refinement_round_{local_round + 1}",
                    evaluation_index=len(rows) + 1,
                    config=config,
                    graphlets=graphlets,
                    context=context,
                    base_updates=base_updates,
                    target_rows=target_rows,
                )
            )
            existing.add(raw_key(raw))
            write_csv(output_dir / "bounded_gate_candidates.csv", rows)

    return rows


def select_best_rows(rows: list[dict[str, Any]], config: dict[str, Any]) -> list[dict[str, Any]]:
    limit = int(config["bounded_trait_space_fast_fit"].get("select_top_n", 12))
    finite = sorted(
        [row for row in rows if np.isfinite(float(row["objective_distance"]))],
        key=lambda row: float(row["objective_distance"]),
    )
    return finite[:limit]


def evaluate_selected_candidates(
    *,
    selected_rows: list[dict[str, Any]],
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    distance_metric_cols = [
        str(x)
        for x in config["equivalence"].get(
            "distance_metric_columns",
            config["equivalence"]["metric_columns"],
        )
    ]
    report_metric_cols = [
        str(x)
        for x in config["equivalence"].get(
            "report_metric_columns",
            config["equivalence"]["metric_columns"],
        )
    ]
    random_walk_cfg = {"analysis": config.get("analysis", {})}
    calibration_cfg = observation_calibration_config(config)
    apply_to_reporting = bool(
        calibration_cfg.get("apply_to_reporting", observation_calibration_enabled(config))
    )
    base_updates = dict(context["fit_summary"]["branches"]["smc_abc_best_raw"]["fitted_parameters"])
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []

    for graphlet in graphlets:
        empirical_values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        stage = fix_stage(config)
        metric = record_metric_values(
            graphlet=graphlet,
            candidate_model=f"empirical_observed__{primary_score}",
            candidate_family="empirical_target",
            parameter_set_id=primary_score,
            values=empirical_values,
            score_column=primary_score,
            metric_rows=metric_rows,
            walsh_rows=walsh_rows,
            gft_rows=gft_rows,
            random_walk_cfg=random_walk_cfg,
            extra={
                "fit_method": stage_prefix(config, "empirical_target"),
                "fix_stage": stage.get("label", "baseline"),
                "fix_stage_index": stage.get("index", 1),
                "fixes_applied": stage.get("fixes", ""),
                "latent_generator": latent_generator_label(config),
                **identity_readout(),
                "readout_applied_to_reporting": False,
            },
        )
        empirical_epi = epistasis_property_metrics(empirical_values, dim=int(graphlet["graphlet_dimension"]))
        empirical_raw_diag = raw_diagnostic_metrics(empirical_values, dim=int(graphlet["graphlet_dimension"]))
        metric.update(empirical_epi)
        metric.update(prefixed_metrics(empirical_epi, "raw_"))
        metric.update(prefixed_numeric_metrics(empirical_raw_diag, "raw_"))

    cache: dict[str, dict[str, Any]] = {}
    for row in selected_rows:
        raw = {
            key.replace("raw__", "", 1): row[key]
            for key in row
            if str(key).startswith("raw__")
        }
        values_by_graphlet = graphlet_values_for_updates(
            updates=candidate_updates(base_updates=base_updates, raw=raw, config=config),
            raw=raw,
            config=config,
            graphlets=graphlets,
            context=context,
            seed_offset=1000 * int(row.get("fix_stage_index", fix_stage(config).get("index", 1))),
        )
        readout = fit_readout_for_values(
            values_by_graphlet=values_by_graphlet,
            graphlets=graphlets,
            config=config,
        )
        cache[str(row["parameter_set_id"])] = {
            "raw": raw,
            "values_by_graphlet": values_by_graphlet,
            "readout": readout,
        }

    for graphlet in graphlets:
        graphlet_id = str(graphlet["graphlet_id"])
        for row in selected_rows:
            cached = cache[str(row["parameter_set_id"])]
            raw_values = np.asarray(cached["values_by_graphlet"][graphlet_id], dtype=np.float64)
            readout = dict(cached["readout"])
            values = apply_observation_readout(raw_values, readout) if apply_to_reporting else raw_values
            stage = fix_stage(config)
            extra = {
                **readout,
                "fit_method": row["fit_method"],
                "fix_stage": row.get("fix_stage", stage.get("label", "baseline")),
                "fix_stage_index": row.get("fix_stage_index", stage.get("index", 1)),
                "fixes_applied": row.get("fixes_applied", stage.get("fixes", "")),
                "latent_generator": row.get(
                    "latent_generator",
                    latent_generator_label(config),
                ),
                "selection_objective_distance": row["objective_distance"],
                "readout_applied_to_reporting": bool(apply_to_reporting),
                **{f"raw__{key}": value for key, value in cached["raw"].items()},
            }
            metric = record_metric_values(
                graphlet=graphlet,
                candidate_model=candidate_model_label(config),
                candidate_family=candidate_family_label(config),
                parameter_set_id=str(row["parameter_set_id"]),
                values=values,
                score_column=primary_score,
                metric_rows=metric_rows,
                walsh_rows=walsh_rows,
                gft_rows=gft_rows,
                random_walk_cfg=random_walk_cfg,
                extra=extra,
            )
            epi = epistasis_property_metrics(values, dim=int(graphlet["graphlet_dimension"]))
            raw_epi = epistasis_property_metrics(raw_values, dim=int(graphlet["graphlet_dimension"]))
            raw_diag = raw_diagnostic_metrics(raw_values, dim=int(graphlet["graphlet_dimension"]))
            metric.update(epi)
            metric.update(prefixed_metrics(raw_epi, "raw_"))
            metric.update(prefixed_numeric_metrics(raw_diag, "raw_"))

    metric_df = pd.DataFrame(metric_rows)
    distances: list[dict[str, Any]] = []
    for split in ["calibration", "test", "all"]:
        split_df = metric_df if split == "all" else metric_df[metric_df["graphlet_split"] == split]
        target_rows = split_df[
            (split_df["candidate_model"] == f"empirical_observed__{primary_score}")
            & (split_df["score_column"] == primary_score)
        ]
        if target_rows.empty:
            continue
        grouped = split_df[split_df["score_column"] == primary_score].groupby(
            ["fit_method", "candidate_family", "candidate_model", "parameter_set_id"],
            dropna=False,
        )
        for keys, model_rows in grouped:
            fit_method, family, candidate_model, parameter_set_id = keys
            if str(candidate_model).startswith("empirical_observed__"):
                continue
            distance, by_metric = metric_distance(model_rows, target_rows, distance_metric_cols)
            penalty, guardrail = guardrail_penalty(by_metric, config)
            first = model_rows.iloc[0].to_dict()
            distances.append(
                {
                    "assay_id": ASSAY_ID,
                    "graphlet_split": split,
                    "fit_method": fit_method,
                    "candidate_family": family,
                    "candidate_model": candidate_model,
                    "parameter_set_id": parameter_set_id,
                    "fix_stage": first.get("fix_stage", stage_label(config)),
                    "fix_stage_index": first.get("fix_stage_index", fix_stage(config).get("index", 1)),
                    "fixes_applied": first.get("fixes_applied", fix_stage(config).get("fixes", "")),
                    "distance_to_empirical_primary": distance,
                    "guardrail_penalty": penalty,
                    **guardrail,
                    "selection_objective_distance": first.get("selection_objective_distance", math.nan),
                    **{f"distance_component__{key}": value for key, value in by_metric.items()},
                }
            )

    summary_rows = summarize_metric_rows(
        metric_rows,
        group_cols=[
            "graphlet_split",
            "fix_stage",
            "fix_stage_index",
            "fixes_applied",
            "fit_method",
            "candidate_family",
            "candidate_model",
            "parameter_set_id",
            "score_column",
        ],
        metric_cols=report_metric_cols,
    )
    return metric_rows, walsh_rows, gft_rows, distances, summary_rows


def comparator_output_dir(config: dict[str, Any], project_root: Path) -> Path:
    path = Path(str(config["abc_comparator"]["output_dir"]))
    if not path.is_absolute():
        path = project_root / path
    if not path.exists():
        raise FileNotFoundError(f"ABC comparator output is missing: {path}")
    return path


def bounded_gate_preflight(
    *,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    modules_by_source: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    return {
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_result": config["source_run"].get("result_id", "RES-030"),
        "source_output_dir": context["source_output_dir"],
        "mode": str(config.get("experiment", {}).get("mode", "phot_bounded_trait_space_current_model_fast_fit")),
        "model_patch": "trait_space_epistasis_bound",
        "primary_score_column": config["graphlet_evaluation"]["primary_score_column"],
        "score_columns": context["assay_score_columns"],
        "objective_metric_columns": [str(x) for x in config["graphlet_evaluation"]["metric_columns"]],
        "n_graphlets": int(len(graphlets)),
        "n_calibration_graphlets": int(sum(1 for graphlet in graphlets if graphlet["graphlet_split"] == "calibration")),
        "n_test_graphlets": int(sum(1 for graphlet in graphlets if graphlet["graphlet_split"] == "test")),
        "graphlet_dimensions": sorted({int(graphlet["graphlet_dimension"]) for graphlet in graphlets}),
        "module_counts": {source: int(len(modules)) for source, modules in modules_by_source.items()},
        "readout_base_anchor": config["empirical_innovation"]["base_anchor"],
        "run077_phot_best_objective": context["source_phot_best_objective"],
        "observation_calibration": observation_calibration_config(config),
        "latent_generator": (
            latent_generator_label(config)
        ),
        "softmin_latent": softmin_latent_config(config),
        "sparse_walsh_field": sparse_walsh_config(config),
        "abc_comparator_run": config["abc_comparator"]["run_id"],
    }


def affine_fit_values(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, float, float, float]:
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    ok = np.isfinite(x_arr) & np.isfinite(y_arr)
    x_arr = x_arr[ok]
    y_arr = y_arr[ok]
    if x_arr.size == 0:
        return np.asarray([], dtype=np.float64), math.nan, math.nan, math.nan, math.nan
    design = np.column_stack([np.ones_like(x_arr), x_arr])
    intercept, slope = np.linalg.lstsq(design, y_arr, rcond=None)[0]
    fitted = intercept + slope * x_arr
    rmse = float(np.sqrt(np.mean((fitted - y_arr) ** 2)))
    corr = (
        float(np.corrcoef(fitted, y_arr)[0, 1])
        if y_arr.size >= 2 and np.std(fitted) > 0.0 and np.std(y_arr) > 0.0
        else math.nan
    )
    return fitted, float(intercept), float(slope), rmse, corr


def graphlet_metric_bundle(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    walsh_rows, variance_by_order = graphlet_runner.walsh_decomposition(arr, dim=int(dim))
    high_order_fraction = float(
        sum(float(value) for order, value in variance_by_order.items() if int(order) >= 3)
    )
    return {
        "total_variance": float(np.var(arr)),
        "walsh_high_order_fraction": high_order_fraction,
        **{f"walsh_variance_order_{order}": float(value) for order, value in variance_by_order.items()},
        **epistasis_property_metrics(arr, dim=int(dim)),
        **raw_diagnostic_metrics(arr, dim=int(dim)),
        "n_walsh_terms": int(len(walsh_rows)),
    }


def sparse_walsh_spectrum_calibration(
    *,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    clusters: list[list[int]],
) -> tuple[list[dict[str, Any]], list[float]]:
    cfg = sparse_walsh_config(config)
    calib_cfg = dict(cfg.get("spectrum_calibration", {}))
    if not bool(calib_cfg.get("enabled", False)):
        weights = sparse_walsh_order_weights(raw={}, config=config)
        return [], weights
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    splits = {str(value) for value in calib_cfg.get("graphlet_splits", ["calibration"])}
    orders = sparse_walsh_orders(config)
    graphlet_positions = [
        set(int(feature[0]) for feature in graphlet["features"])
        for graphlet in graphlets
        if str(graphlet.get("graphlet_split", "")) in splits
    ]
    rows_by_order: dict[int, list[float]] = {order: [] for order in orders}
    fraction_by_order: dict[int, list[float]] = {order: [] for order in orders}
    for graphlet in graphlets:
        if str(graphlet.get("graphlet_split", "")) not in splits:
            continue
        dim = int(graphlet["graphlet_dimension"])
        values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        _walsh_rows, fractions = graphlet_runner.walsh_decomposition(values, dim=dim)
        total_variance = float(np.var(values))
        for order in orders:
            fraction = float(fractions.get(str(order), fractions.get(order, 0.0)))
            fraction_by_order[order].append(fraction)
            rows_by_order[order].append(fraction * total_variance)

    rng = np.random.RandomState(int(calib_cfg.get("seed", cfg.get("seed", 20260617))))
    bootstrap_iterations = max(int(calib_cfg.get("bootstrap_iterations", 200)), 1)
    floor = float(calib_cfg.get("noise_floor_fraction", 0.01))
    mc_samples = max(int(calib_cfg.get("support_monte_carlo_samples", 2000)), 1)
    support_probs: dict[int, float] = {}
    cluster_pool = [tuple(sorted({int(pos) for pos in cluster})) for cluster in clusters if len(set(cluster)) >= 2]
    all_positions = sorted({pos for positions in graphlet_positions for pos in positions})
    for order in orders:
        hits = []
        for _ in range(mc_samples):
            if cluster_pool:
                cluster = list(cluster_pool[int(rng.randint(len(cluster_pool)))])
            else:
                cluster = list(all_positions)
            if len(cluster) < order:
                support = set(rng.choice(np.asarray(all_positions, dtype=np.int64), size=min(order, len(all_positions)), replace=False))
            else:
                support = set(rng.choice(np.asarray(cluster, dtype=np.int64), size=order, replace=False))
            if not support or not graphlet_positions:
                hits.append(0.0)
            else:
                hits.append(float(np.mean([support.issubset(positions) for positions in graphlet_positions])))
        support_probs[order] = float(np.mean(hits)) if hits else 0.0

    raw_global: dict[int, float] = {}
    out_rows: list[dict[str, Any]] = []
    for order in orders:
        local_values = np.asarray(rows_by_order[order], dtype=np.float64)
        fraction_values = np.asarray(fraction_by_order[order], dtype=np.float64)
        if local_values.size:
            boot_means = []
            for _ in range(bootstrap_iterations):
                sample = rng.choice(local_values, size=local_values.size, replace=True)
                boot_means.append(float(np.mean(sample)))
            lower = float(np.quantile(boot_means, 0.10))
            upper = float(np.quantile(boot_means, 0.90))
            local_mean = float(np.mean(local_values))
            fraction_mean = float(np.mean(fraction_values))
        else:
            lower = upper = local_mean = fraction_mean = 0.0
        passes = bool(fraction_mean > floor and upper > floor)
        hit_prob = max(float(support_probs.get(order, 0.0)), float(calib_cfg.get("min_support_probability", 1e-4)))
        global_sigma2 = float(local_mean / hit_prob) if passes else 0.0
        raw_global[order] = global_sigma2
        out_rows.append(
            {
                "assay_id": ASSAY_ID,
                "calibration_stage": "G",
                "order": int(order),
                "n_graphlets": int(local_values.size),
                "local_variance_mean": local_mean,
                "local_variance_bootstrap_q10": lower,
                "local_variance_bootstrap_q90": upper,
                "local_fraction_mean": fraction_mean,
                "noise_floor_fraction": floor,
                "passes_noise_floor": passes,
                "support_hit_probability": float(support_probs.get(order, 0.0)),
                "global_sigma2_unscaled": global_sigma2,
            }
        )
    total = float(sum(raw_global.values()))
    if total <= 0.0:
        weights = sparse_walsh_order_weights(raw={}, config=config)
    else:
        weights = [float(raw_global[order] / total) for order in orders]
    for row, weight in zip(out_rows, weights):
        row["global_order_weight_normalized"] = float(weight)
    return out_rows, weights


def sparse_walsh_with_order_weights(config: dict[str, Any], weights: list[float]) -> dict[str, Any]:
    scoped = copy.deepcopy(config)
    cfg = dict(scoped["bounded_trait_space_fast_fit"].get("sparse_walsh_field", {}))
    cfg["order_weights"] = [float(value) for value in weights]
    scoped["bounded_trait_space_fast_fit"]["sparse_walsh_field"] = cfg
    return scoped


def scaled_sparse_walsh_weights(
    *,
    config: dict[str, Any],
    base_weights: list[float],
    high_order_multiplier: float,
) -> list[float]:
    orders = sparse_walsh_orders(config)
    weights = []
    for order, weight in zip(orders, base_weights):
        scale = float(high_order_multiplier) if int(order) >= 3 else 1.0
        weights.append(max(float(weight) * scale, 0.0))
    total = float(sum(weights))
    if total <= 0.0:
        return sparse_walsh_order_weights(raw={}, config=config)
    return [float(value / total) for value in weights]


def run_sparse_walsh_forward_probe(
    *,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    base_weights: list[float],
) -> list[dict[str, Any]]:
    cfg = sparse_walsh_config(config)
    probe_cfg = dict(cfg.get("forward_probe", {}))
    if not bool(probe_cfg.get("enabled", False)):
        return []
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    splits = {str(value) for value in probe_cfg.get("graphlet_splits", ["calibration"])}
    target_metrics = [
        str(value)
        for value in probe_cfg.get(
            "target_metrics",
            [
                "n_local_maxima",
                "accessible_path_fraction",
                "magnitude_epistasis_to_single_mutant_sd_ratio",
                "sign_epistasis_fraction_local",
                "reciprocal_sign_epistasis_fraction_local",
                "diminishing_returns_epistasis_fraction",
                "walsh_high_order_fraction",
            ],
        )
    ]
    selected_graphlets = [
        graphlet for graphlet in graphlets if str(graphlet.get("graphlet_split", "")) in splits
    ]
    empirical_by_metric: dict[str, list[float]] = {metric: [] for metric in target_metrics}
    for graphlet in selected_graphlets:
        bundle = graphlet_metric_bundle(
            np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64),
            dim=int(graphlet["graphlet_dimension"]),
        )
        for metric in target_metrics:
            value = bundle.get(metric)
            if value is not None and np.isfinite(float(value)):
                empirical_by_metric[metric].append(float(value))

    base_updates = dict(context["fit_summary"]["branches"]["smc_abc_best_raw"]["fitted_parameters"])
    epi_grid = [float(value) for value in probe_cfg.get("epi_weight_grid", [0.01, 0.025, 0.05])]
    high_order_grid = [float(value) for value in probe_cfg.get("high_order_weight_multipliers", [0.75, 1.0, 1.25])]
    n_realizations = max(int(probe_cfg.get("n_realizations", 8)), 1)
    central_low = float(probe_cfg.get("central_interval_low", 0.10))
    central_high = float(probe_cfg.get("central_interval_high", 0.90))
    rows: list[dict[str, Any]] = []
    stage = fix_stage(config)
    candidate_index = 0
    for high_multiplier in high_order_grid:
        weights = scaled_sparse_walsh_weights(
            config=config,
            base_weights=base_weights,
            high_order_multiplier=high_multiplier,
        )
        scoped = sparse_walsh_with_order_weights(config, weights)
        for epi_weight in epi_grid:
            candidate_index += 1
            model_by_metric: dict[str, list[float]] = {metric: [] for metric in target_metrics}
            for realization in range(n_realizations):
                scoped_realization = copy.deepcopy(scoped)
                field_cfg = dict(scoped_realization["bounded_trait_space_fast_fit"]["sparse_walsh_field"])
                field_cfg["seed"] = int(field_cfg.get("seed", 20260617)) + int(realization)
                scoped_realization["bounded_trait_space_fast_fit"]["sparse_walsh_field"] = field_cfg
                raw = {
                    "sparse_walsh_epi_weight": epi_weight,
                    "sparse_walsh_terms": int(field_cfg.get("terms", 1024)),
                    **{
                        f"sparse_walsh_order{order}_weight": weight
                        for order, weight in zip(sparse_walsh_orders(scoped_realization), weights)
                    },
                }
                values_by_graphlet = graphlet_values_for_updates(
                    updates=candidate_updates(base_updates=base_updates, raw=raw, config=scoped_realization),
                    raw=raw,
                    config=scoped_realization,
                    graphlets=selected_graphlets,
                    context=context,
                    seed_offset=1000 * int(stage.get("index", 1)) + int(realization),
                )
                for graphlet in selected_graphlets:
                    bundle = graphlet_metric_bundle(
                        np.asarray(values_by_graphlet[str(graphlet["graphlet_id"])], dtype=np.float64),
                        dim=int(graphlet["graphlet_dimension"]),
                    )
                    for metric in target_metrics:
                        value = bundle.get(metric)
                        if value is not None and np.isfinite(float(value)):
                            model_by_metric[metric].append(float(value))
            metric_passes: dict[str, bool] = {}
            for metric in target_metrics:
                empirical_values = empirical_by_metric.get(metric, [])
                model_values = model_by_metric.get(metric, [])
                empirical_mean = float(np.mean(empirical_values)) if empirical_values else math.nan
                if model_values:
                    model_q_low = float(np.quantile(model_values, central_low))
                    model_q_high = float(np.quantile(model_values, central_high))
                    model_mean = float(np.mean(model_values))
                    passes = bool(
                        np.isfinite(empirical_mean)
                        and model_q_low <= empirical_mean <= model_q_high
                    )
                else:
                    model_q_low = model_q_high = model_mean = math.nan
                    passes = False
                metric_passes[metric] = passes
                rows.append(
                    {
                        "assay_id": ASSAY_ID,
                        "probe_stage": "H1",
                        "candidate_form": "sparse_walsh_field",
                        "fix_stage": stage.get("label", "baseline"),
                        "fix_stage_index": stage.get("index", 1),
                        "fixes_applied": stage.get("fixes", ""),
                        "probe_candidate_index": int(candidate_index),
                        "epi_weight": float(epi_weight),
                        "high_order_weight_multiplier": float(high_multiplier),
                        "order_weights": json.dumps(weights, separators=(",", ":")),
                        "metric": metric,
                        "empirical_mean": empirical_mean,
                        "model_mean": model_mean,
                        "model_central_q_low": model_q_low,
                        "model_central_q_high": model_q_high,
                        "central_interval_low": central_low,
                        "central_interval_high": central_high,
                        "metric_pass": passes,
                    }
                )
            joint_pass = bool(metric_passes and all(metric_passes.values()))
            for row in rows:
                if int(row.get("probe_candidate_index", -1)) == candidate_index:
                    row["joint_pass"] = joint_pass
    return rows


def run_softmin_capability_probe(
    *,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    cfg = dict(config["bounded_trait_space_fast_fit"].get("softmin_capability_probe", {}))
    if not bool(cfg.get("enabled", False)):
        return []
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    k_values = [int(value) for value in cfg.get("k_values", [2, 3, 5, 8])]
    tau_values = [float(value) for value in cfg.get("tau_values", [0.05, 0.1, 0.25, 0.5, 1.0])]
    n_candidates = max(int(cfg.get("random_candidates", 96)), 1)
    coefficient_noise = float(cfg.get("coefficient_noise", 1.0))
    intercept_spread = float(cfg.get("intercept_spread", 1.0))
    splits = {str(value) for value in cfg.get("graphlet_splits", ["calibration"])}
    random_walk_cfg = {"analysis": config.get("analysis", {})}
    stage = fix_stage(config)
    rows: list[dict[str, Any]] = []
    for graphlet in graphlets:
        if str(graphlet.get("graphlet_split", "")) not in splits:
            continue
        dim = int(graphlet["graphlet_dimension"])
        y = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        design = binary_design_matrix(dim)
        base_beta = np.linalg.lstsq(design, y, rcond=None)[0]
        y_sd = max(float(np.std(y)), 1e-8)
        for k in k_values:
            best: dict[str, Any] | None = None
            seed = graphlet_runner.stable_seed(
                "model-free-softmin-capability",
                str(graphlet["graphlet_id"]),
                str(stage.get("label", "baseline")),
                str(k),
                str(cfg.get("seed", 20260617)),
            )
            rng = np.random.RandomState(seed)
            for tau in tau_values:
                for candidate_index in range(n_candidates):
                    traits: list[np.ndarray] = []
                    for _trait_index in range(max(k, 2)):
                        beta = np.asarray(base_beta, dtype=np.float64).copy()
                        beta[0] += float(rng.normal(0.0, intercept_spread * y_sd))
                        beta[1:] += rng.normal(0.0, coefficient_noise * y_sd, size=dim)
                        traits.append(design @ beta)
                    raw_values = soft_min_many_values(traits, tau=float(tau) * y_sd)
                    fitted, intercept, slope, rmse, corr = affine_fit_values(zscore(raw_values), y)
                    if fitted.size != y.size:
                        continue
                    if best is None or rmse < float(best["probe_fit_rmse"]):
                        best = {
                            "probe_k_traits": int(k),
                            "probe_tau": float(tau),
                            "probe_candidate_index": int(candidate_index),
                            "probe_fit_intercept": intercept,
                            "probe_fit_slope": slope,
                            "probe_fit_rmse": rmse,
                            "probe_fit_correlation": corr,
                            "values": fitted,
                        }
            if best is None:
                continue
            metric_rows: list[dict[str, Any]] = []
            walsh_rows: list[dict[str, Any]] = []
            gft_rows: list[dict[str, Any]] = []
            values = np.asarray(best.pop("values"), dtype=np.float64)
            metric = record_metric_values(
                graphlet=graphlet,
                candidate_model="model_free_k_trait_softmin",
                candidate_family="model_free_capability_probe",
                parameter_set_id=f"model_free_softmin__k{k}__{graphlet['graphlet_id']}",
                values=values,
                score_column=primary_score,
                metric_rows=metric_rows,
                walsh_rows=walsh_rows,
                gft_rows=gft_rows,
                random_walk_cfg=random_walk_cfg,
                extra={
                    "fit_method": stage_prefix(config, f"model_free_k{k}"),
                    "fix_stage": stage.get("label", "baseline"),
                    "fix_stage_index": stage.get("index", 1),
                    "fixes_applied": stage.get("fixes", ""),
                    "latent_generator": "model_free_additive_k_trait_softmin",
                    **best,
                },
            )
            epi = epistasis_property_metrics(values, dim=dim)
            diag = raw_diagnostic_metrics(values, dim=dim)
            target_epi = epistasis_property_metrics(y, dim=dim)
            target_diag = raw_diagnostic_metrics(y, dim=dim)
            row = {
                **metric,
                **epi,
                **prefixed_numeric_metrics(diag, ""),
                **{f"target_{key}": value for key, value in {**target_epi, **target_diag}.items()},
                **best,
                "assay_id": ASSAY_ID,
                "graphlet_id": graphlet["graphlet_id"],
                "graphlet_split": graphlet["graphlet_split"],
                "graphlet_dimension": dim,
                "fix_stage": stage.get("label", "baseline"),
                "fix_stage_index": stage.get("index", 1),
                "fixes_applied": stage.get("fixes", ""),
            }
            rows.append(row)
    return rows


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=args.quick)
    random_walk_cfg = {"analysis": config.get("analysis", {})}

    graphlets, context, modules_by_source, _readout = build_graphlets(
        config=config,
        project_root=project_root,
        output_dir=output_dir,
    )
    config = inject_sparse_walsh_cluster_positions(
        config=config,
        graphlets=graphlets,
        modules_by_source=modules_by_source,
    )
    preflight = bounded_gate_preflight(
        config=config,
        graphlets=graphlets,
        context=context,
        modules_by_source=modules_by_source,
    )
    preflight.update(
        {
            "model_patch": "trait_space_epistasis_bound",
            "fit_protocol": "deterministic_random_search_with_local_refinement",
            "observation_calibration": observation_calibration_config(config),
            "abc_comparator_run": config["abc_comparator"]["run_id"],
        }
    )
    write_json(output_dir / "preflight.json", preflight)

    comparator_dir = comparator_output_dir(config, project_root)
    comparator_selected = pd.read_csv(comparator_dir / "selected_family_distances.csv")
    comparator_selected.to_csv(output_dir / "abc_comparator_selected_family_distances.csv", index=False)

    all_candidate_rows: list[dict[str, Any]] = []
    all_selected_rows: list[dict[str, Any]] = []
    all_metric_rows: list[dict[str, Any]] = []
    all_walsh_rows: list[dict[str, Any]] = []
    all_gft_rows: list[dict[str, Any]] = []
    all_distances: list[dict[str, Any]] = []
    all_summary_rows: list[dict[str, Any]] = []
    all_probe_rows: list[dict[str, Any]] = []
    all_sparse_spectrum_rows: list[dict[str, Any]] = []
    all_sparse_forward_probe_rows: list[dict[str, Any]] = []
    stage_summaries: list[dict[str, Any]] = []
    capability_probe_ran = False
    sparse_walsh_calibrated_weights: list[float] | None = None
    sparse_walsh_probe_passed = False
    for scoped_config in stage_configs(config):
        if sparse_walsh_enabled(scoped_config) and sparse_walsh_calibrated_weights is not None:
            scoped_config = sparse_walsh_with_order_weights(
                scoped_config,
                sparse_walsh_calibrated_weights,
            )
        stage = fix_stage(scoped_config)
        walsh_cfg = sparse_walsh_config(scoped_config)
        spectrum_rows, spectrum_weights = sparse_walsh_spectrum_calibration(
            config=scoped_config,
            graphlets=graphlets,
            clusters=list(walsh_cfg.get("cluster_positions", [])),
        )
        if spectrum_rows:
            all_sparse_spectrum_rows.extend(spectrum_rows)
            sparse_walsh_calibrated_weights = spectrum_weights
            scoped_config = sparse_walsh_with_order_weights(scoped_config, spectrum_weights)
            walsh_cfg = sparse_walsh_config(scoped_config)

        forward_rows = run_sparse_walsh_forward_probe(
            config=scoped_config,
            graphlets=graphlets,
            context=context,
            base_weights=sparse_walsh_calibrated_weights
            if sparse_walsh_calibrated_weights is not None
            else sparse_walsh_order_weights(raw={}, config=scoped_config),
        )
        if forward_rows:
            all_sparse_forward_probe_rows.extend(forward_rows)
            sparse_walsh_probe_passed = sparse_walsh_probe_passed or any(
                bool(row.get("joint_pass", False)) for row in forward_rows
            )

        probe_cfg = dict(scoped_config["bounded_trait_space_fast_fit"].get("softmin_capability_probe", {}))
        if bool(probe_cfg.get("enabled", False)) and (
            not capability_probe_ran or not bool(probe_cfg.get("run_once", True))
        ):
            probe_rows = run_softmin_capability_probe(
                config=scoped_config,
                graphlets=graphlets,
            )
            all_probe_rows.extend(probe_rows)
            capability_probe_ran = True
        fast_fit_enabled = bool(scoped_config["bounded_trait_space_fast_fit"].get("enabled", True))
        if bool(walsh_cfg.get("requires_forward_probe_pass", False)) and not sparse_walsh_probe_passed:
            fast_fit_enabled = False
        if not fast_fit_enabled:
            stage_summaries.append(
                {
                    "fix_stage": stage.get("label", "baseline"),
                    "fix_stage_index": stage.get("index", 1),
                    "fixes_applied": stage.get("fixes", ""),
                    "fast_fit_enabled": False,
                    "sparse_walsh_forward_probe_passed": bool(sparse_walsh_probe_passed),
                    "n_candidate_evaluations": 0,
                    "n_selected_candidates": 0,
                    "best_objective_distance": math.nan,
                }
            )
            continue
        target_rows = target_metric_rows(
            graphlets=graphlets,
            config=scoped_config,
            random_walk_cfg=random_walk_cfg,
        )
        candidate_rows = run_fast_fit(
            config=scoped_config,
            output_dir=output_dir,
            graphlets=graphlets,
            context=context,
            target_rows=target_rows,
        )
        selected_rows = select_best_rows(candidate_rows, scoped_config)
        metric_rows, walsh_rows, gft_rows, distances, summary_rows = evaluate_selected_candidates(
            selected_rows=selected_rows,
            config=scoped_config,
            graphlets=graphlets,
            context=context,
        )
        all_candidate_rows.extend(candidate_rows)
        all_selected_rows.extend(selected_rows)
        all_metric_rows.extend(metric_rows)
        all_walsh_rows.extend(walsh_rows)
        all_gft_rows.extend(gft_rows)
        all_distances.extend(distances)
        all_summary_rows.extend(summary_rows)
        finite = [float(row["objective_distance"]) for row in candidate_rows if np.isfinite(float(row["objective_distance"]))]
        stage_summaries.append(
            {
                "fix_stage": stage.get("label", "baseline"),
                "fix_stage_index": stage.get("index", 1),
                "fixes_applied": stage.get("fixes", ""),
                "fast_fit_enabled": True,
                "n_candidate_evaluations": int(len(candidate_rows)),
                "n_selected_candidates": int(len(selected_rows)),
                "best_objective_distance": float(min(finite)) if finite else math.nan,
            }
        )
    selected_rows = all_selected_rows
    candidate_rows = all_candidate_rows
    metric_rows = all_metric_rows
    walsh_rows = all_walsh_rows
    gft_rows = all_gft_rows
    distances = all_distances
    summary_rows = all_summary_rows
    write_csv(output_dir / "selected_bounded_gate_candidates.csv", selected_rows)

    equivalence_rows = build_equivalence_rows(
        [row for row in distances if row.get("fit_method") != "empirical_target"],
        comparator_selected,
        config,
    )
    distance_by_parameter_set = {str(row["parameter_set_id"]): row for row in distances}
    for row in equivalence_rows:
        source = distance_by_parameter_set.get(str(row.get("parameter_set_id")), {})
        row["fix_stage"] = source.get("fix_stage", "")
        row["fixes_applied"] = source.get("fixes_applied", "")

    write_csv(output_dir / "graphlet_catalog.csv", graphlet_catalog_rows(graphlets))
    write_csv(output_dir / "module_catalog.csv", module_catalog_rows(modules_by_source))
    write_csv(output_dir / "bounded_gate_candidates.csv", candidate_rows)
    write_csv(output_dir / "graphlet_metrics.csv", metric_rows)
    write_csv(output_dir / "walsh_terms.csv", walsh_rows)
    write_csv(output_dir / "gft_coefficients.csv", gft_rows)
    write_csv(output_dir / "graphlet_summary.csv", summary_rows)
    write_csv(output_dir / "bounded_gate_distances.csv", distances)
    write_csv(output_dir / "softmin_capability_probe.csv", all_probe_rows)
    write_csv(output_dir / "sparse_walsh_spectrum_calibration.csv", all_sparse_spectrum_rows)
    write_csv(output_dir / "sparse_walsh_forward_probe.csv", all_sparse_forward_probe_rows)
    write_csv_with_columns(
        output_dir / "bounded_gate_equivalence.csv",
        equivalence_rows,
        columns=EQUIVALENCE_COLUMNS,
    )
    write_json(
        output_dir / "summary.json",
        {
            "experiment_id": config["experiment"]["id"],
            "run_mode": config["experiment"]["mode"],
            "abc_comparator": {
                "run_id": config["abc_comparator"]["run_id"],
                "output_dir": str(comparator_dir),
            },
            "n_graphlets": int(len(graphlets)),
            "n_candidate_evaluations": int(len(candidate_rows)),
            "n_selected_candidates": int(len(selected_rows)),
            "n_softmin_capability_probe_rows": int(len(all_probe_rows)),
            "n_sparse_walsh_spectrum_rows": int(len(all_sparse_spectrum_rows)),
            "n_sparse_walsh_forward_probe_rows": int(len(all_sparse_forward_probe_rows)),
            "sparse_walsh_forward_probe_passed": bool(sparse_walsh_probe_passed),
            "n_equivalence_rows": int(len(equivalence_rows)),
            "n_operationally_equivalent": int(
                sum(bool(row["operationally_equivalent"]) for row in equivalence_rows)
            ),
            "fix_stage_summaries": stage_summaries,
            "best_objective_distance": (
                float(min(row["objective_distance"] for row in candidate_rows))
                if candidate_rows
                else math.nan
            ),
            "output_files": {
                "preflight": str(output_dir / "preflight.json"),
                "bounded_gate_candidates": str(output_dir / "bounded_gate_candidates.csv"),
                "selected_bounded_gate_candidates": str(output_dir / "selected_bounded_gate_candidates.csv"),
                "bounded_gate_distances": str(output_dir / "bounded_gate_distances.csv"),
                "bounded_gate_equivalence": str(output_dir / "bounded_gate_equivalence.csv"),
                "graphlet_metrics": str(output_dir / "graphlet_metrics.csv"),
                "graphlet_summary": str(output_dir / "graphlet_summary.csv"),
                "softmin_capability_probe": str(output_dir / "softmin_capability_probe.csv"),
                "sparse_walsh_spectrum_calibration": str(
                    output_dir / "sparse_walsh_spectrum_calibration.csv"
                ),
                "sparse_walsh_forward_probe": str(output_dir / "sparse_walsh_forward_probe.csv"),
                "walsh_terms": str(output_dir / "walsh_terms.csv"),
                "gft_coefficients": str(output_dir / "gft_coefficients.csv"),
                "abc_comparator_selected_family_distances": str(
                    output_dir / "abc_comparator_selected_family_distances.csv"
                ),
            },
        },
    )


if __name__ == "__main__":
    main()
