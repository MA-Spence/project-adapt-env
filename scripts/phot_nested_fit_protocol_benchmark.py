#!/usr/bin/env python3
"""Compare cheaper PHOT nested-min fitting protocols against full ABC.

This runner keeps the RUN-098 nested soft-min/min residual model family fixed
and asks whether cheaper search protocols select candidates that are
operationally equivalent to the previously computed full ABC fit.
"""

from __future__ import annotations

import argparse
import copy
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
    graphlet_catalog_rows,
    graphlet_runner,
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
from phot_nested_softmin_ensemble_abc import (  # noqa: E402
    build_graphlets,
    decode_second_parameters,
    evaluate_nested_graphlet,
    fit_locked_first_stage_transform,
    specs_for_family,
    write_preflight,
)
from project_adapt_env.coupled_landscape import (  # noqa: E402
    aggregate_module_energy,
    energy_locality_for_modules,
)
from project_adapt_env.smc_abc import (  # noqa: E402
    ParameterSpec,
    SMCABCBackendConfig,
    SMCABCCheckpointConfig,
    SMCABCConfig,
    SimulationResult,
    run_smc_abc,
)
from proteingym_system_specific_current_graphlet_refit import metric_distance  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def config_with_metric_columns(config: dict[str, Any], metric_columns: list[str]) -> dict[str, Any]:
    scoped = copy.deepcopy(config)
    scoped["graphlet_evaluation"]["metric_columns"] = list(metric_columns)
    return scoped


def raw_key(raw: dict[str, float | int]) -> str:
    return json.dumps(raw, sort_keys=True, separators=(",", ":"))


def sample_specs(specs: list[ParameterSpec], rng: np.random.RandomState) -> dict[str, float | int]:
    return {spec.name: spec.sample_prior(rng) for spec in specs}


def raw_to_unit_vector(raw: dict[str, float | int], specs: list[ParameterSpec]) -> np.ndarray:
    values = []
    for spec in specs:
        width = max(float(spec.high - spec.low), 1e-12)
        values.append((float(raw[spec.name]) - float(spec.low)) / width)
    return np.asarray(values, dtype=np.float64)


def unit_vector_to_raw(vector: np.ndarray, specs: list[ParameterSpec]) -> dict[str, float | int]:
    clipped = np.clip(np.asarray(vector, dtype=np.float64), 0.0, 1.0)
    raw: dict[str, float | int] = {}
    for value, spec in zip(clipped, specs):
        scaled = float(spec.low) + float(value) * float(spec.high - spec.low)
        if spec.kind == "int":
            raw[spec.name] = int(np.clip(round(scaled), int(spec.low), int(spec.high)))
        else:
            raw[spec.name] = float(np.clip(scaled, float(spec.low), float(spec.high)))
    return raw


def polynomial_feature_matrix(x: np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=np.float64)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    columns = [np.ones(arr.shape[0], dtype=np.float64), arr]
    columns.append(arr * arr)
    pairwise = []
    for left in range(arr.shape[1]):
        for right in range(left + 1, arr.shape[1]):
            pairwise.append(arr[:, left] * arr[:, right])
    if pairwise:
        columns.append(np.vstack(pairwise).T)
    return np.hstack([col if col.ndim == 2 else col.reshape(-1, 1) for col in columns])


def fit_ridge_predictor(
    raw_rows: list[dict[str, float | int]],
    distances: list[float],
    specs: list[ParameterSpec],
    *,
    ridge_alpha: float,
) -> np.ndarray:
    x = np.vstack([raw_to_unit_vector(raw, specs) for raw in raw_rows])
    y = np.log1p(np.asarray(distances, dtype=np.float64))
    design = polynomial_feature_matrix(x)
    penalty = float(ridge_alpha) * np.eye(design.shape[1], dtype=np.float64)
    penalty[0, 0] = 0.0
    return np.linalg.solve(design.T @ design + penalty, design.T @ y)


def predict_ridge(raw_rows: list[dict[str, float | int]], specs: list[ParameterSpec], beta: np.ndarray) -> np.ndarray:
    x = np.vstack([raw_to_unit_vector(raw, specs) for raw in raw_rows])
    pred = polynomial_feature_matrix(x) @ np.asarray(beta, dtype=np.float64)
    return np.expm1(pred)


def narrow_specs_around(
    specs: list[ParameterSpec],
    center: dict[str, float | int],
    *,
    width_fraction: float,
) -> list[ParameterSpec]:
    narrowed: list[ParameterSpec] = []
    for spec in specs:
        support = float(spec.high - spec.low)
        half_width = max(support * float(width_fraction) * 0.5, 1.0 if spec.kind == "int" else support * 0.01)
        midpoint = float(center[spec.name])
        low = max(float(spec.low), midpoint - half_width)
        high = min(float(spec.high), midpoint + half_width)
        if spec.kind == "int":
            low = math.floor(low)
            high = math.ceil(high)
            if low > high:
                low = high = round(midpoint)
        elif high <= low:
            low = max(float(spec.low), midpoint - support * 0.01)
            high = min(float(spec.high), midpoint + support * 0.01)
        narrowed.append(
            ParameterSpec(
                name=spec.name,
                kind=spec.kind,
                low=float(low),
                high=float(high),
                initial_scale_fraction=spec.initial_scale_fraction,
                min_scale_fraction=spec.min_scale_fraction,
            )
        )
    return narrowed


def summarize_numeric(values: list[float], prefix: str) -> dict[str, float]:
    arr = np.asarray([float(v) for v in values if np.isfinite(float(v))], dtype=np.float64)
    if arr.size == 0:
        return {
            f"{prefix}_mean": math.nan,
            f"{prefix}_median": math.nan,
            f"{prefix}_sd": math.nan,
            f"{prefix}_q10": math.nan,
            f"{prefix}_q90": math.nan,
            f"{prefix}_max": math.nan,
        }
    return {
        f"{prefix}_mean": float(np.mean(arr)),
        f"{prefix}_median": float(np.median(arr)),
        f"{prefix}_sd": float(np.std(arr)),
        f"{prefix}_q10": float(np.quantile(arr, 0.10)),
        f"{prefix}_q90": float(np.quantile(arr, 0.90)),
        f"{prefix}_max": float(np.max(arr)),
    }


def epistasis_property_metrics(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    magnitudes: list[float] = []
    signed: list[float] = []
    single_effects: list[float] = []
    sign_events = 0
    reciprocal_events = 0
    diminishing_events = 0
    count = 0
    for axis in range(int(dim)):
        for base in range(1 << int(dim)):
            if base & (1 << axis):
                continue
            single_effects.append(float(arr[base | (1 << axis)] - arr[base]))
    for left in range(int(dim)):
        for right in range(left + 1, int(dim)):
            others = [idx for idx in range(int(dim)) if idx not in {left, right}]
            for bg_index in range(1 << len(others)):
                base = 0
                for offset, bit in enumerate(others):
                    if bg_index & (1 << offset):
                        base |= 1 << bit
                w00 = float(arr[base])
                w10 = float(arr[base | (1 << left)])
                w01 = float(arr[base | (1 << right)])
                w11 = float(arr[base | (1 << left) | (1 << right)])
                left_on_bg = w10 - w00
                left_on_right = w11 - w01
                right_on_bg = w01 - w00
                right_on_left = w11 - w10
                epsilon = w11 - w10 - w01 + w00
                left_flip = left_on_bg * left_on_right < 0.0
                right_flip = right_on_bg * right_on_left < 0.0
                sign_events += int(left_flip or right_flip)
                reciprocal_events += int(left_flip and right_flip)
                diminishing_events += int(
                    epsilon < 0.0 and (left_on_bg > 0.0 or right_on_bg > 0.0 or left_on_right > 0.0 or right_on_left > 0.0)
                )
                signed.append(epsilon)
                magnitudes.append(abs(epsilon))
                count += 1
    single_effect_arr = np.asarray(single_effects, dtype=np.float64)
    single_effect_sd = (
        float(np.std(single_effect_arr)) if single_effect_arr.size else math.nan
    )
    magnitude_mean = (
        float(np.mean(np.asarray(magnitudes, dtype=np.float64)))
        if magnitudes
        else math.nan
    )
    return {
        **summarize_numeric(magnitudes, "magnitude_epistasis"),
        **summarize_numeric(signed, "signed_pairwise_epistasis"),
        "single_mutant_effect_sd": single_effect_sd,
        "magnitude_epistasis_to_single_mutant_sd_ratio": (
            float(magnitude_mean / max(single_effect_sd, 1e-12))
            if np.isfinite(magnitude_mean) and np.isfinite(single_effect_sd)
            else math.nan
        ),
        "sign_epistasis_fraction_local": float(sign_events / count) if count else math.nan,
        "reciprocal_sign_epistasis_fraction_local": float(reciprocal_events / count) if count else math.nan,
        "diminishing_returns_epistasis_fraction": float(diminishing_events / count) if count else math.nan,
    }


def prefixed_metrics(metrics: dict[str, float], prefix: str) -> dict[str, float]:
    return {f"{prefix}{key}": float(value) for key, value in metrics.items()}


def observation_calibration_config(config: dict[str, Any]) -> dict[str, Any]:
    return dict(config.get("observation_calibration", {}))


def observation_calibration_enabled(config: dict[str, Any]) -> bool:
    return bool(observation_calibration_config(config).get("enabled", False))


def fit_positive_affine_readout(
    model_values: list[np.ndarray],
    empirical_values: list[np.ndarray],
    *,
    min_slope: float = 1e-12,
) -> dict[str, float | str]:
    x = np.concatenate([np.asarray(values, dtype=np.float64).ravel() for values in model_values])
    y = np.concatenate([np.asarray(values, dtype=np.float64).ravel() for values in empirical_values])
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size == 0 or y.size == 0:
        return {
            "readout_method": "positive_affine",
            "readout_intercept": 0.0,
            "readout_slope": 1.0,
            "readout_fit_n": 0,
            "readout_fit_rmse": math.nan,
            "readout_fit_correlation": math.nan,
        }
    x_mean = float(np.mean(x))
    y_mean = float(np.mean(y))
    x_var = float(np.mean((x - x_mean) ** 2))
    if x_var <= 1e-24:
        slope = 0.0
    else:
        slope = float(np.mean((x - x_mean) * (y - y_mean)) / x_var)
    slope = max(float(slope), float(min_slope))
    intercept = y_mean - slope * x_mean
    fitted = intercept + slope * x
    rmse = float(np.sqrt(np.mean((fitted - y) ** 2))) if x.size else math.nan
    corr = float(np.corrcoef(fitted, y)[0, 1]) if x.size >= 2 and np.std(fitted) > 0 and np.std(y) > 0 else math.nan
    return {
        "readout_method": "positive_affine",
        "readout_intercept": float(intercept),
        "readout_slope": float(slope),
        "readout_fit_n": int(x.size),
        "readout_fit_rmse": rmse,
        "readout_fit_correlation": corr,
    }


def apply_affine_readout(values: np.ndarray, readout: dict[str, float | str]) -> np.ndarray:
    return float(readout.get("readout_intercept", 0.0)) + float(readout.get("readout_slope", 1.0)) * np.asarray(
        values,
        dtype=np.float64,
    )


def identity_readout() -> dict[str, float | str]:
    return {
        "readout_method": "identity",
        "readout_intercept": 0.0,
        "readout_slope": 1.0,
        "readout_fit_n": 0,
        "readout_fit_rmse": 0.0,
        "readout_fit_correlation": 1.0,
    }


def target_metric_rows(
    *,
    graphlets: list[dict[str, Any]],
    config: dict[str, Any],
    random_walk_cfg: dict[str, Any],
) -> pd.DataFrame:
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    rows: list[dict[str, Any]] = []
    walsh: list[dict[str, Any]] = []
    gft: list[dict[str, Any]] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        metric = record_metric_values(
            graphlet=graphlet,
            candidate_model=f"empirical_observed__{primary_score}",
            candidate_family="empirical_target",
            parameter_set_id=primary_score,
            values=values,
            score_column=primary_score,
            metric_rows=rows,
            walsh_rows=walsh,
            gft_rows=gft,
            random_walk_cfg=random_walk_cfg,
            extra={"base_anchor": "", "module_source": "", **identity_readout()},
        )
        epi = epistasis_property_metrics(values, dim=int(graphlet["graphlet_dimension"]))
        metric.update(epi)
        metric.update(prefixed_metrics(epi, "raw_"))
    return pd.DataFrame(rows)


def simulate_family_particle(
    *,
    raw: dict[str, float | int],
    family: str,
    parameter_set_id: str,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    target_rows: pd.DataFrame,
    random_walk_cfg: dict[str, Any],
) -> SimulationResult:
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    metric_cols = [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]
    params = decode_second_parameters(raw, config=config, family=family, parameter_set_id=parameter_set_id)
    calibration_cfg = observation_calibration_config(config)
    min_slope = float(calibration_cfg.get("min_slope", 1e-12))
    apply_to_objective = bool(calibration_cfg.get("apply_to_objective", observation_calibration_enabled(config)))

    prepared: list[dict[str, Any]] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        nested = evaluate_nested_graphlet(
            graphlet=graphlet,
            config=config,
            modules_by_source=modules_by_source,
            second_params=params,
        )
        prepared.append(
            {
                "graphlet": graphlet,
                "nested": nested,
                "raw_values": np.asarray(nested["values"], dtype=np.float64),
                "empirical_values": np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64),
            }
        )

    if observation_calibration_enabled(config):
        readout = fit_positive_affine_readout(
            [item["raw_values"] for item in prepared],
            [item["empirical_values"] for item in prepared],
            min_slope=min_slope,
        )
    else:
        readout = identity_readout()

    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    active_rows: list[dict[str, Any]] = []
    for item in prepared:
        graphlet = item["graphlet"]
        raw_values = np.asarray(item["raw_values"], dtype=np.float64)
        values = apply_affine_readout(raw_values, readout) if apply_to_objective else raw_values
        extra = {
            **dict(item["nested"]["extra"]),
            **readout,
            "readout_applied_to_objective": bool(apply_to_objective),
        }
        metric = record_metric_values(
            graphlet=graphlet,
            candidate_model=str(params["candidate_model"]),
            candidate_family=family,
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
        metric.update(epi)
        metric.update(prefixed_metrics(raw_epi, "raw_"))
        active_rows.append({**extra, **epi, **prefixed_metrics(raw_epi, "raw_")})

    model_rows = pd.DataFrame(metric_rows)
    distance, by_metric = metric_distance(model_rows, target_rows, metric_cols)
    penalty, guardrail = guardrail_penalty(by_metric, config)
    if not np.isfinite(distance):
        distance = float("inf")
    objective = float(distance + penalty)
    summary = [float(distance), penalty] + [float(by_metric.get(metric, math.nan)) for metric in metric_cols]
    extras = {
        **params,
        **readout,
        "readout_applied_to_objective": bool(apply_to_objective),
        "calibration_graphlet_distance": float(distance),
        "objective_distance": objective,
        **guardrail,
        **{f"distance_component__{key}": value for key, value in by_metric.items()},
    }
    if active_rows:
        active_df = pd.DataFrame(active_rows)
        for col in [
            "values_to_c_correlation",
            "values_to_a_correlation",
            "second_auxiliary_to_c_correlation",
            "second_innovation_to_c_sd_ratio",
            "total_innovation_to_a_sd_ratio",
            "second_auxiliary_weight_sd",
            "second_raw_sd",
            "second_module_raw_sd",
            "second_n_module_coefficients",
            "magnitude_epistasis_mean",
            "raw_magnitude_epistasis_mean",
            "signed_pairwise_epistasis_mean",
            "raw_signed_pairwise_epistasis_mean",
        ]:
            if col in active_df:
                values = pd.to_numeric(active_df[col], errors="coerce").dropna()
                if not values.empty:
                    extras[f"{col}_mean"] = float(values.mean())
    return SimulationResult(distance=objective, summary_vector=np.asarray(summary, dtype=np.float64), extras=extras)


def evaluate_raw_candidate(
    *,
    raw: dict[str, float | int],
    family: str,
    fit_method: str,
    objective_name: str,
    parameter_set_id: str,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    target_rows: pd.DataFrame,
    random_walk_cfg: dict[str, Any],
    evaluation_index: int,
) -> dict[str, Any]:
    result = simulate_family_particle(
        raw=raw,
        family=family,
        parameter_set_id=parameter_set_id,
        config=config,
        graphlets=graphlets,
        modules_by_source=modules_by_source,
        target_rows=target_rows,
        random_walk_cfg=random_walk_cfg,
    )
    return {
        "fit_method": fit_method,
        "objective_name": objective_name,
        "candidate_family": family,
        "parameter_set_id": parameter_set_id,
        "evaluation_index": int(evaluation_index),
        "objective_distance": float(result.distance),
        "summary_vector": ";".join(
            "nan" if not np.isfinite(float(x)) else f"{float(x):.12g}"
            for x in result.summary_vector
        ),
        **{f"raw__{key}": value for key, value in raw.items()},
        **result.extras,
    }


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


def deterministic_summary_loss_search(
    *,
    family: str,
    specs: list[ParameterSpec],
    config: dict[str, Any],
    method_cfg: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    target_rows: pd.DataFrame,
    random_walk_cfg: dict[str, Any],
    seed: int,
    output_dir: Path,
) -> list[dict[str, Any]]:
    rng = np.random.RandomState(seed)
    rows: list[dict[str, Any]] = []
    candidates = [sample_specs(specs, rng) for _ in range(int(method_cfg["n_initial"]))]
    best_raw: dict[str, float | int] | None = None
    best_distance = math.inf
    eval_index = 0
    for phase in ["initial"] + [f"local_round_{idx + 1}" for idx in range(int(method_cfg["n_local_rounds"]))]:
        for raw in unique_raws(candidates):
            eval_index += 1
            row = evaluate_raw_candidate(
                raw=raw,
                family=family,
                fit_method="deterministic_summary_loss",
                objective_name="full_summary",
                parameter_set_id=f"{family}__deterministic_{eval_index:04d}",
                config=config,
                graphlets=graphlets,
                modules_by_source=modules_by_source,
                target_rows=target_rows,
                random_walk_cfg=random_walk_cfg,
                evaluation_index=eval_index,
            )
            row["search_phase"] = phase
            rows.append(row)
            if float(row["objective_distance"]) < best_distance:
                best_distance = float(row["objective_distance"])
                best_raw = raw
        write_csv(output_dir / "fit_protocol_candidates.csv", rows)
        if best_raw is None or not phase.startswith("local_round_") and phase != "initial":
            continue
        center = raw_to_unit_vector(best_raw, specs)
        scale = float(method_cfg["local_scale_fraction"]) * (0.5 ** max(0, len([r for r in rows if str(r.get("search_phase", "")).startswith("local_round_")])))
        candidates = [
            unit_vector_to_raw(center + rng.normal(scale=scale, size=center.shape), specs)
            for _ in range(int(method_cfg["n_local_candidates"]))
        ]
    return rows


def surrogate_assisted_search(
    *,
    family: str,
    specs: list[ParameterSpec],
    config: dict[str, Any],
    method_cfg: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    target_rows: pd.DataFrame,
    random_walk_cfg: dict[str, Any],
    seed: int,
    output_dir: Path,
    prior_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rng = np.random.RandomState(seed)
    rows: list[dict[str, Any]] = []
    train_raws = [sample_specs(specs, rng) for _ in range(int(method_cfg["n_train"]))]
    eval_index = 0
    for raw in unique_raws(train_raws):
        eval_index += 1
        row = evaluate_raw_candidate(
            raw=raw,
            family=family,
            fit_method="surrogate_assisted",
            objective_name="full_summary",
            parameter_set_id=f"{family}__surrogate_train_{eval_index:04d}",
            config=config,
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=target_rows,
            random_walk_cfg=random_walk_cfg,
            evaluation_index=eval_index,
        )
        row["search_phase"] = "surrogate_train"
        rows.append(row)
    if prior_rows:
        for row in prior_rows:
            raw = {key.replace("raw__", "", 1): row[key] for key in row if str(key).startswith("raw__")}
            if all(spec.name in raw for spec in specs):
                eval_index += 1
                copied = dict(row)
                copied["fit_method"] = "surrogate_assisted"
                copied["parameter_set_id"] = f"{family}__surrogate_reused_{eval_index:04d}"
                copied["search_phase"] = "surrogate_reused"
                copied["evaluation_index"] = eval_index
                rows.append(copied)
    beta = fit_ridge_predictor(
        [{key.replace("raw__", "", 1): row[key] for key in row if str(key).startswith("raw__")} for row in rows],
        [float(row["objective_distance"]) for row in rows],
        specs,
        ridge_alpha=float(method_cfg.get("ridge_alpha", 1e-3)),
    )
    pool = unique_raws([sample_specs(specs, rng) for _ in range(int(method_cfg["n_candidate_pool"]))])
    predicted = predict_ridge(pool, specs, beta)
    order = np.argsort(predicted)[: int(method_cfg["n_verify"])]
    for rank, pool_index in enumerate(order, start=1):
        raw = pool[int(pool_index)]
        eval_index += 1
        row = evaluate_raw_candidate(
            raw=raw,
            family=family,
            fit_method="surrogate_assisted",
            objective_name="full_summary",
            parameter_set_id=f"{family}__surrogate_verify_{rank:04d}",
            config=config,
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=target_rows,
            random_walk_cfg=random_walk_cfg,
            evaluation_index=eval_index,
        )
        row["search_phase"] = "surrogate_verify"
        row["surrogate_predicted_distance"] = float(predicted[int(pool_index)])
        rows.append(row)
    write_csv(output_dir / "fit_protocol_candidates.csv", rows)
    return rows


def exact_gate_search(
    *,
    family: str,
    specs: list[ParameterSpec],
    full_config: dict[str, Any],
    gate_config: dict[str, Any],
    method_cfg: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    gate_target_rows: pd.DataFrame,
    full_target_rows: pd.DataFrame,
    random_walk_cfg: dict[str, Any],
    seed: int,
    output_dir: Path,
) -> list[dict[str, Any]]:
    gate_rows = deterministic_summary_loss_search(
        family=family,
        specs=specs,
        config=gate_config,
        method_cfg={
            "n_initial": int(method_cfg["n_initial"]),
            "n_local_rounds": int(method_cfg["n_local_rounds"]),
            "n_local_candidates": int(method_cfg["n_local_candidates"]),
            "local_scale_fraction": float(method_cfg["local_scale_fraction"]),
        },
        graphlets=graphlets,
        modules_by_source=modules_by_source,
        target_rows=gate_target_rows,
        random_walk_cfg=random_walk_cfg,
        seed=seed,
        output_dir=output_dir,
    )
    rows: list[dict[str, Any]] = []
    best = sorted(gate_rows, key=lambda row: float(row["objective_distance"]))[: int(method_cfg["n_promote"])]
    for index, gate_row in enumerate(best, start=1):
        raw = {key.replace("raw__", "", 1): gate_row[key] for key in gate_row if str(key).startswith("raw__")}
        row = evaluate_raw_candidate(
            raw=raw,
            family=family,
            fit_method="exact_descriptive_gated",
            objective_name="full_summary_after_gate_selection",
            parameter_set_id=f"{family}__exact_gate_promoted_{index:04d}",
            config=full_config,
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=full_target_rows,
            random_walk_cfg=random_walk_cfg,
            evaluation_index=index,
        )
        row["search_phase"] = "gate_promoted_full_eval"
        row["gate_objective_distance"] = float(gate_row["objective_distance"])
        rows.append(row)
    return rows


def local_abc_after_cheap_fit(
    *,
    family: str,
    specs: list[ParameterSpec],
    config: dict[str, Any],
    method_cfg: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    target_rows: pd.DataFrame,
    random_walk_cfg: dict[str, Any],
    seed_row: dict[str, Any],
    output_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    center = {key.replace("raw__", "", 1): seed_row[key] for key in seed_row if str(key).startswith("raw__")}
    narrowed = narrow_specs_around(specs, center, width_fraction=float(method_cfg["prior_width_fraction"]))
    smc = SMCABCConfig(
        n_particles=int(method_cfg["n_particles"]),
        n_rounds=int(method_cfg["n_rounds"]),
        initial_pool_multiplier=int(method_cfg["initial_pool_multiplier"]),
        proposal_pool_multiplier=int(method_cfg["proposal_pool_multiplier"]),
        epsilon_quantile=float(method_cfg["epsilon_quantile"]),
        seed=int(method_cfg["seed"]) + int(graphlet_runner.stable_seed("local-abc", family) % 100000),
    )
    backend = SMCABCBackendConfig(
        kind=str(method_cfg.get("backend", {}).get("kind", "serial")),
        max_workers=int(method_cfg.get("backend", {}).get("max_workers", 1)),
        batch_size=int(method_cfg.get("backend", {}).get("batch_size", 1)),
    )
    checkpoint = SMCABCCheckpointConfig(
        path=str(output_dir / "checkpoints" / f"local_abc_after_cheap_fit_{family}.json"),
        resume=bool(method_cfg.get("resume_checkpoints", True)),
    )

    def simulate(raw: dict[str, float | int]) -> SimulationResult:
        return simulate_family_particle(
            raw=raw,
            family=family,
            parameter_set_id=f"{family}__local_abc_trial",
            config=config,
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=target_rows,
            random_walk_cfg=random_walk_cfg,
        )

    result = run_smc_abc(
        specs=narrowed,
        config=smc,
        simulate=simulate,
        backend=backend,
        checkpoint=checkpoint,
        progress_callback=lambda event: write_json(output_dir / "progress.json", {"local_abc": event}),
        run_label=f"local_abc_after_cheap_fit_{family}",
    )
    rows: list[dict[str, Any]] = []
    for rank, particle in enumerate(sorted(result.particles, key=lambda item: item.distance), start=1):
        extras = dict(particle.extras)
        for key in ["fit_method", "objective_name", "candidate_family", "parameter_set_id", "objective_distance"]:
            extras.pop(key, None)
        rows.append(
            {
                **extras,
                "fit_method": "local_abc_after_cheap_fit",
                "objective_name": "full_summary_local_abc",
                "candidate_family": family,
                "posterior_rank": int(rank),
                "weight": float(particle.weight),
                "objective_distance": float(particle.distance),
                "parameter_set_id": f"{family}__local_abc_after_cheap_fit_rank{rank:03d}",
                "seed_parameter_set_id": seed_row.get("parameter_set_id", ""),
                "seed_fit_method": seed_row.get("fit_method", ""),
                "summary_vector": ";".join(
                    "nan" if not np.isfinite(float(x)) else f"{float(x):.12g}"
                    for x in particle.summary_vector
                ),
                **{f"raw__{key}": value for key, value in particle.parameters.items()},
            }
        )
    round_rows = [
        {
            "fit_method": "local_abc_after_cheap_fit",
            "candidate_family": family,
            "round_index": int(row.round_index),
            "epsilon": float(row.epsilon),
            "attempts": int(row.attempts),
            "best_distance": float(row.best_distance),
            "median_distance": float(row.median_distance),
        }
        for row in result.round_summaries
    ]
    spec_rows = [
        {
            "fit_method": "local_abc_after_cheap_fit",
            "candidate_family": family,
            "name": spec.name,
            "kind": spec.kind,
            "low": float(spec.low),
            "high": float(spec.high),
            "center": center.get(spec.name, math.nan),
        }
        for spec in narrowed
    ]
    return rows, round_rows, spec_rows


def candidate_raw(row: dict[str, Any]) -> dict[str, float | int]:
    return {key.replace("raw__", "", 1): row[key] for key in row if str(key).startswith("raw__")}


def selected_full_abc_particle_rows(
    *,
    config: dict[str, Any],
    comparator_dir: Path,
    comparator_selected: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Return RUN-098 selected posterior particles for same-path reporting.

    The RUN-098 selected-family distance table is the primary comparator for
    predeclared equivalence. These rows let us re-evaluate the selected full
    ABC particles through any newer descriptive readouts added here without
    treating those readouts as part of the original equivalence distance.
    """
    posterior_path = comparator_dir / "posterior_particles.csv"
    if not posterior_path.exists():
        return []
    posterior = pd.read_csv(posterior_path)
    selected = comparator_selected[comparator_selected["graphlet_split"] == "calibration"].copy()
    rows: list[dict[str, Any]] = []
    for family in [str(x) for x in config["abc"]["fit_families"]]:
        family_selected = selected[selected["candidate_family"].astype(str) == family]
        if family_selected.empty:
            continue
        rank_value = family_selected.iloc[0].get("posterior_rank", math.nan)
        if not np.isfinite(float(rank_value)):
            continue
        rank = int(round(float(rank_value)))
        family_particles = posterior[
            (posterior["candidate_family"].astype(str) == family)
            & (pd.to_numeric(posterior["posterior_rank"], errors="coerce") == rank)
        ]
        if family_particles.empty:
            continue
        row = family_particles.iloc[0].to_dict()
        row["fit_method"] = "full_abc_comparator"
        row["objective_name"] = "RUN-098_selected_full_abc"
        row["candidate_family"] = family
        row["source_run_id"] = str(config["abc_comparator"]["run_id"])
        row["source_result_id"] = str(config["abc_comparator"].get("result_id", ""))
        row["source_parameter_set_id"] = row.get("parameter_set_id", "")
        row["parameter_set_id"] = f"{family}__RUN-098_posterior_rank{rank:03d}"
        row["objective_distance"] = float(row.get("objective_distance", math.nan))
        rows.append(row)
    return rows


def select_best_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    df = pd.DataFrame(rows)
    if df.empty:
        return selected
    df["objective_distance"] = pd.to_numeric(df["objective_distance"], errors="coerce")
    for (method, family), sub in df.groupby(["fit_method", "candidate_family"], dropna=False):
        sub = sub[np.isfinite(sub["objective_distance"])]
        if sub.empty:
            continue
        selected.append(sub.sort_values("objective_distance").iloc[0].to_dict())
    return selected


def selected_candidate_cache(
    *,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    selected_rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    calibration_cfg = observation_calibration_config(config)
    min_slope = float(calibration_cfg.get("min_slope", 1e-12))
    cache: dict[str, dict[str, Any]] = {}
    for row in selected_rows:
        family = str(row["candidate_family"])
        parameter_set_id = str(row["parameter_set_id"])
        decoded = decode_second_parameters(
            candidate_raw(row),
            config=config,
            family=family,
            parameter_set_id=parameter_set_id,
        )
        entries: dict[str, dict[str, Any]] = {}
        for graphlet in graphlets:
            nested = evaluate_nested_graphlet(
                graphlet=graphlet,
                config=config,
                modules_by_source=modules_by_source,
                second_params=decoded,
            )
            entries[str(graphlet["graphlet_id"])] = {
                "nested": nested,
                "raw_values": np.asarray(nested["values"], dtype=np.float64),
                "empirical_values": np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64),
            }
        calibration_entries = [
            entries[str(graphlet["graphlet_id"])]
            for graphlet in graphlets
            if graphlet["graphlet_split"] == "calibration"
        ]
        if observation_calibration_enabled(config):
            readout = fit_positive_affine_readout(
                [item["raw_values"] for item in calibration_entries],
                [item["empirical_values"] for item in calibration_entries],
                min_slope=min_slope,
            )
        else:
            readout = identity_readout()
        cache[parameter_set_id] = {
            "decoded": decoded,
            "entries": entries,
            "readout": readout,
        }
    return cache


def evaluate_selected_candidates(
    *,
    config: dict[str, Any],
    output_dir: Path,
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    selected_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
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
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    locality_rows: list[dict[str, Any]] = []
    locality_summary_rows: list[dict[str, Any]] = []
    candidate_cache = selected_candidate_cache(
        config=config,
        graphlets=graphlets,
        modules_by_source=modules_by_source,
        selected_rows=selected_rows,
    )
    calibration_cfg = observation_calibration_config(config)
    apply_to_reporting = bool(calibration_cfg.get("apply_to_reporting", observation_calibration_enabled(config)))

    for graphlet in graphlets:
        empirical_values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
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
            extra={"fit_method": "empirical_target", **identity_readout(), "readout_applied_to_reporting": False},
        )
        empirical_epi = epistasis_property_metrics(empirical_values, dim=int(graphlet["graphlet_dimension"]))
        metric.update(empirical_epi)
        metric.update(prefixed_metrics(empirical_epi, "raw_"))

        for row in selected_rows:
            family = str(row["candidate_family"])
            cached = candidate_cache[str(row["parameter_set_id"])]
            nested = cached["entries"][str(graphlet["graphlet_id"])]["nested"]
            raw_values = np.asarray(cached["entries"][str(graphlet["graphlet_id"])]["raw_values"], dtype=np.float64)
            readout = dict(cached["readout"])
            values = apply_affine_readout(raw_values, readout) if apply_to_reporting else raw_values
            candidate_model = f"{row['fit_method']}__nested_softmin_{family}"
            extra = {
                **dict(nested["extra"]),
                **readout,
                "fit_method": row["fit_method"],
                "objective_name": row.get("objective_name", ""),
                "selection_objective_distance": row.get("objective_distance", math.nan),
                "source_parameter_set_id": row.get("parameter_set_id", ""),
                "readout_applied_to_reporting": bool(apply_to_reporting),
            }
            metric = record_metric_values(
                graphlet=graphlet,
                candidate_model=candidate_model,
                candidate_family=family,
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
            metric.update(epi)
            metric.update(prefixed_metrics(raw_epi, "raw_"))
            loc_rows, loc_summary = energy_locality_for_modules(
                graphlet_id=str(graphlet["graphlet_id"]),
                graphlet_split=str(graphlet["graphlet_split"]),
                graphlet_features=tuple(graphlet["features"]),
                innovation=np.asarray(nested["second_innovation"], dtype=np.float64),
                modules=nested["second_modules"],
                candidate_model=candidate_model,
                parameter_set_id=str(row["parameter_set_id"]),
            )
            for loc_row in loc_rows:
                loc_row.update(extra)
            loc_summary.update(extra)
            locality_rows.extend(loc_rows)
            locality_summary_rows.append(loc_summary)

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
        for keys, model_rows in split_df[split_df["score_column"] == primary_score].groupby(
            ["fit_method", "candidate_family", "candidate_model", "parameter_set_id"],
            dropna=False,
        ):
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
                    "distance_to_empirical_primary": distance,
                    "guardrail_penalty": penalty,
                    **guardrail,
                    "selection_objective_distance": first.get("selection_objective_distance", math.nan),
                    **{f"distance_component__{key}": value for key, value in by_metric.items()},
                }
            )

    summary_rows = summarize_metric_rows(
        metric_rows,
        group_cols=["graphlet_split", "fit_method", "candidate_family", "candidate_model", "parameter_set_id", "score_column"],
        metric_cols=report_metric_cols,
    )
    module_summary = aggregate_module_energy(locality_rows)
    return metric_rows, walsh_rows, gft_rows, distances, summary_rows, module_summary


def comparator_output_dir(config: dict[str, Any], project_root: Path) -> Path:
    path = Path(str(config["abc_comparator"]["output_dir"]))
    if not path.is_absolute():
        path = project_root / path
    if not path.exists():
        raise FileNotFoundError(f"ABC comparator output is missing: {path}")
    return path


def build_equivalence_rows(
    protocol_distances: list[dict[str, Any]],
    comparator_selected: pd.DataFrame,
    config: dict[str, Any],
    recomputed_comparator_distances: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    margins = dict(config["equivalence"])
    abs_margin = float(margins["absolute_distance_margin"])
    rel_margin = float(margins["relative_distance_margin"])
    max_component_margin = float(margins["max_component_distance_margin"])
    test_protocol = [row for row in protocol_distances if row.get("graphlet_split") == "test"]
    if bool(margins.get("use_recomputed_full_abc_reference", False)) and recomputed_comparator_distances is not None:
        comparator_test = pd.DataFrame(
            [
                row
                for row in recomputed_comparator_distances
                if row.get("graphlet_split") == "test" and row.get("fit_method") == "full_abc_comparator"
            ]
        )
    else:
        comparator_test = comparator_selected[comparator_selected["graphlet_split"] == "test"].copy()
    if comparator_test.empty:
        return []
    comparator_test["distance_to_empirical_primary"] = pd.to_numeric(
        comparator_test["distance_to_empirical_primary"],
        errors="coerce",
    )
    best_any = comparator_test.sort_values("distance_to_empirical_primary").head(1)
    rows: list[dict[str, Any]] = []
    for row in test_protocol:
        family = str(row["candidate_family"])
        same = comparator_test[comparator_test["candidate_family"].astype(str) == family]
        references = []
        if not same.empty:
            references.append(("same_family_full_abc", same.sort_values("distance_to_empirical_primary").iloc[0]))
        if not best_any.empty:
            references.append(("best_any_full_abc", best_any.iloc[0]))
        for reference_name, reference in references:
            candidate_distance = float(row["distance_to_empirical_primary"])
            reference_distance = float(reference["distance_to_empirical_primary"])
            delta = candidate_distance - reference_distance
            component_values = []
            for key, value in row.items():
                if not str(key).startswith("distance_component__"):
                    continue
                try:
                    parsed = float(value)
                except (TypeError, ValueError):
                    continue
                if np.isfinite(parsed):
                    component_values.append(parsed)
            max_component = max(component_values) if component_values else math.nan
            within_absolute = delta <= abs_margin
            within_relative = delta <= rel_margin * max(abs(reference_distance), 1e-12)
            no_component_failure = (not np.isfinite(max_component)) or max_component <= max_component_margin
            rows.append(
                {
                    "fit_method": row["fit_method"],
                    "candidate_family": family,
                    "candidate_model": row["candidate_model"],
                    "parameter_set_id": row["parameter_set_id"],
                    "reference": reference_name,
                    "reference_candidate_family": reference["candidate_family"],
                    "reference_parameter_set_id": reference["parameter_set_id"],
                    "candidate_test_distance": candidate_distance,
                    "reference_test_distance": reference_distance,
                    "delta_test_distance": delta,
                    "relative_delta_test_distance": delta / max(abs(reference_distance), 1e-12),
                    "max_component_distance": max_component,
                    "within_absolute_margin": bool(within_absolute),
                    "within_relative_margin": bool(within_relative),
                    "no_component_distance_failure": bool(no_component_failure),
                    "operationally_equivalent": bool((within_absolute or within_relative) and no_component_failure),
                }
            )
    return rows


def run_all_protocols(
    *,
    config: dict[str, Any],
    gate_config: dict[str, Any],
    output_dir: Path,
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    full_target_rows: pd.DataFrame,
    gate_target_rows: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    methods = dict(config["fitting_methods"])
    random_walk_cfg = {"analysis": config.get("analysis", {})}
    all_rows: list[dict[str, Any]] = []
    abc_rounds: list[dict[str, Any]] = []
    abc_specs: list[dict[str, Any]] = []
    for family in [str(x) for x in config["abc"]["fit_families"]]:
        specs = specs_for_family(config, family)
        base_seed = int(methods.get("seed", 20260616)) + int(graphlet_runner.stable_seed("fit-protocol", family) % 100000)
        deterministic = deterministic_summary_loss_search(
            family=family,
            specs=specs,
            config=config,
            method_cfg=dict(methods["deterministic_summary_loss"]),
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=full_target_rows,
            random_walk_cfg=random_walk_cfg,
            seed=base_seed,
            output_dir=output_dir,
        )
        all_rows.extend(deterministic)
        surrogate = surrogate_assisted_search(
            family=family,
            specs=specs,
            config=config,
            method_cfg=dict(methods["surrogate_assisted"]),
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=full_target_rows,
            random_walk_cfg=random_walk_cfg,
            seed=base_seed + 17,
            output_dir=output_dir,
            prior_rows=deterministic,
        )
        all_rows.extend(surrogate)
        gated = exact_gate_search(
            family=family,
            specs=specs,
            full_config=config,
            gate_config=gate_config,
            method_cfg=dict(methods["exact_descriptive_gated"]),
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            gate_target_rows=gate_target_rows,
            full_target_rows=full_target_rows,
            random_walk_cfg=random_walk_cfg,
            seed=base_seed + 31,
            output_dir=output_dir,
        )
        all_rows.extend(gated)
        seed_rows = sorted(
            [row for row in deterministic + surrogate + gated if np.isfinite(float(row["objective_distance"]))],
            key=lambda row: float(row["objective_distance"]),
        )
        if seed_rows:
            local_rows, rounds, narrowed_specs = local_abc_after_cheap_fit(
                family=family,
                specs=specs,
                config=config,
                method_cfg=dict(methods["local_abc_after_cheap_fit"]),
                graphlets=graphlets,
                modules_by_source=modules_by_source,
                target_rows=full_target_rows,
                random_walk_cfg=random_walk_cfg,
                seed_row=seed_rows[0],
                output_dir=output_dir,
            )
            all_rows.extend(local_rows)
            abc_rounds.extend(rounds)
            abc_specs.extend(narrowed_specs)
        write_csv(output_dir / "fit_protocol_candidates.csv", all_rows)
        write_csv(output_dir / "local_abc_rounds.csv", abc_rounds)
        write_csv(output_dir / "local_abc_parameter_specs.csv", abc_specs)
    return all_rows, abc_rounds, abc_specs


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=args.quick)
    random_walk_cfg = {"analysis": config.get("analysis", {})}
    gate_metrics = [str(x) for x in config["fitting_methods"]["exact_descriptive_gated"]["metric_columns"]]
    gate_config = config_with_metric_columns(config, gate_metrics)

    graphlets, context, modules_by_source, _readout = build_graphlets(
        config=config,
        project_root=project_root,
        output_dir=output_dir,
    )
    config["_runtime_first_stage_transform"] = fit_locked_first_stage_transform(
        config=config,
        graphlets=graphlets,
        output_dir=output_dir,
    )
    gate_config["_runtime_first_stage_transform"] = config["_runtime_first_stage_transform"]
    preflight = write_preflight(
        config=config,
        output_dir=output_dir,
        graphlets=graphlets,
        context=context,
        modules_by_source=modules_by_source,
    )
    preflight.update(
        {
            "protocols_tested": list(config["fitting_methods"].keys()),
            "abc_comparator_run": config["abc_comparator"]["run_id"],
            "equivalence_metric_columns": config["equivalence"].get(
                "distance_metric_columns",
                config["equivalence"]["metric_columns"],
            ),
            "report_metric_columns": config["equivalence"].get(
                "report_metric_columns",
                config["equivalence"]["metric_columns"],
            ),
            "exact_gate_metric_columns": gate_metrics,
            "observation_calibration": observation_calibration_config(config),
        }
    )
    write_json(output_dir / "preflight.json", preflight)
    full_target_rows = target_metric_rows(
        graphlets=graphlets,
        config=config,
        random_walk_cfg=random_walk_cfg,
    )
    gate_target_rows = target_metric_rows(
        graphlets=graphlets,
        config=gate_config,
        random_walk_cfg=random_walk_cfg,
    )

    candidate_rows, abc_rounds, abc_specs = run_all_protocols(
        config=config,
        gate_config=gate_config,
        output_dir=output_dir,
        graphlets=graphlets,
        modules_by_source=modules_by_source,
        full_target_rows=full_target_rows,
        gate_target_rows=gate_target_rows,
    )
    selected_rows = select_best_rows(candidate_rows)
    write_csv(output_dir / "selected_fit_protocol_candidates.csv", selected_rows)

    comparator_dir = comparator_output_dir(config, project_root)
    comparator_selected = pd.read_csv(comparator_dir / "selected_family_distances.csv")
    comparator_selected.to_csv(output_dir / "abc_comparator_selected_family_distances.csv", index=False)
    comparator_particle_rows = selected_full_abc_particle_rows(
        config=config,
        comparator_dir=comparator_dir,
        comparator_selected=comparator_selected,
    )
    write_csv(output_dir / "full_abc_comparator_selected_particles.csv", comparator_particle_rows)
    evaluation_rows = selected_rows + comparator_particle_rows

    metric_rows, walsh_rows, gft_rows, distances, summary_rows, module_summary = evaluate_selected_candidates(
        config=config,
        output_dir=output_dir,
        graphlets=graphlets,
        modules_by_source=modules_by_source,
        selected_rows=evaluation_rows,
    )
    equivalence_rows = build_equivalence_rows(
        [row for row in distances if row.get("fit_method") != "full_abc_comparator"],
        comparator_selected,
        config,
        recomputed_comparator_distances=distances,
    )

    write_csv(output_dir / "graphlet_catalog.csv", graphlet_catalog_rows(graphlets))
    write_csv(output_dir / "module_catalog.csv", module_catalog_rows(modules_by_source))
    write_csv(output_dir / "fit_protocol_candidates.csv", candidate_rows)
    write_csv(output_dir / "local_abc_rounds.csv", abc_rounds)
    write_csv(output_dir / "local_abc_parameter_specs.csv", abc_specs)
    write_csv(output_dir / "graphlet_metrics.csv", metric_rows)
    write_csv(output_dir / "walsh_terms.csv", walsh_rows)
    write_csv(output_dir / "gft_coefficients.csv", gft_rows)
    write_csv(output_dir / "graphlet_summary.csv", summary_rows)
    write_csv(output_dir / "fit_protocol_distances.csv", distances)
    write_csv(output_dir / "fit_protocol_equivalence.csv", equivalence_rows)
    write_csv(output_dir / "module_energy_summary.csv", module_summary)
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
            "n_equivalence_rows": int(len(equivalence_rows)),
            "n_operationally_equivalent": int(sum(bool(row["operationally_equivalent"]) for row in equivalence_rows)),
            "protocols": sorted({str(row["fit_method"]) for row in candidate_rows}),
            "output_files": {
                "preflight": str(output_dir / "preflight.json"),
                "fit_protocol_candidates": str(output_dir / "fit_protocol_candidates.csv"),
                "selected_fit_protocol_candidates": str(output_dir / "selected_fit_protocol_candidates.csv"),
                "fit_protocol_distances": str(output_dir / "fit_protocol_distances.csv"),
                "fit_protocol_equivalence": str(output_dir / "fit_protocol_equivalence.csv"),
                "graphlet_metrics": str(output_dir / "graphlet_metrics.csv"),
                "graphlet_summary": str(output_dir / "graphlet_summary.csv"),
                "walsh_terms": str(output_dir / "walsh_terms.csv"),
                "gft_coefficients": str(output_dir / "gft_coefficients.csv"),
                "module_energy_summary": str(output_dir / "module_energy_summary.csv"),
                "abc_comparator_selected_family_distances": str(output_dir / "abc_comparator_selected_family_distances.csv"),
                "full_abc_comparator_selected_particles": str(output_dir / "full_abc_comparator_selected_particles.csv"),
            },
        },
    )


if __name__ == "__main__":
    main()
