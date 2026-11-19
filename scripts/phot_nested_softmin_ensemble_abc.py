#!/usr/bin/env python3
"""PHOT nested soft-min ensemble ABC capacity diagnostic.

This runner tests whether adding a second independently parameterised latent
landscape to an already coupled PHOT model improves held-out graphlet and local
ruggedness statistics. It is deliberately a capacity diagnostic, not evidence
for residue-specific epistatic modules.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
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

from empirical_graphlet_falsification_suite import subset_masks  # noqa: E402
from phot_coupled_family_discriminating_abc import (  # noqa: E402
    ASSAY_ID,
    BASELINE_FAMILIES,
    FAMILY_MODULE_SOURCES,
    LOCAL_RUGGEDNESS_METRIC_COLUMNS,
    build_graphlets,
    calibration_target_rows,
    graphlet_catalog_rows,
    graphlet_runner,
    guardrail_penalty,
    module_catalog_rows,
    parameter_specs_from_rows,
    progress_writer,
    record_metric_values,
    top_distance_rows,
)
from phot_coupled_landscape_grid_search import (  # noqa: E402
    correlation_row,
    load_config,
    resolve_source_output_dir,
    summarize_metric_rows,
    write_csv,
    write_json,
)
from phot_hard_min_coupled_landscape_grid_search import (  # noqa: E402
    apply_auxiliary_transform,
    auxiliary_raw_for_features,
    build_auxiliary_components,
    evaluate_sequences,
    fit_auxiliary_transform,
    parse_feature_lists,
)
from project_adapt_env.coupled_landscape import (  # noqa: E402
    aggregate_module_energy,
    boltzmann_two_state_values,
    edge_difference_values,
    energy_locality_for_modules,
    module_residual_coefficients,
    module_residual_values_for_feature_sets,
    node_edge_gradients,
    safe_corr,
    soft_state_metrics,
    zscore,
)
from project_adapt_env.smc_abc import (  # noqa: E402
    ParameterSpec,
    SMCABCBackendConfig,
    SMCABCCheckpointConfig,
    SMCABCConfig,
    SimulationResult,
    run_smc_abc,
)
from proteingym_system_specific_current_graphlet_refit import (  # noqa: E402
    build_landscape,
    gnk_values,
    metric_distance,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def specs_for_family(config: dict[str, Any], family: str) -> list[ParameterSpec]:
    priors = dict(config["abc"]["parameter_priors"])
    names = list(config["abc"]["family_parameter_names"].get(family, []))
    rows = [{**priors[name], "name": name} for name in names]
    return parameter_specs_from_rows(rows)


def density_choice(*, raw: dict[str, float | int], config: dict[str, Any], prefix: str) -> float:
    choices = [float(x) for x in config["nested_operator"]["second_stage"]["auxiliary_pairwise_density_choices"]]
    idx = int(raw.get(f"{prefix}_auxiliary_pairwise_density_index", 0))
    idx = max(0, min(idx, len(choices) - 1))
    return float(choices[idx])


def tau_from_log(value: float | int) -> float:
    return float(math.exp(float(value)))


def prefixed(values: dict[str, Any], prefix: str) -> dict[str, Any]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def graphlet_node_feature_lists(features: tuple[tuple[int, str], ...]) -> list[tuple[tuple[int, str], ...]]:
    return [tuple(sorted(mask, key=lambda item: (int(item[0]), str(item[1])))) for mask in subset_masks(features)]


def first_stage_parameters(config: dict[str, Any]) -> dict[str, Any]:
    stage = dict(config["nested_operator"]["first_stage"])
    module_source = str(stage.get("module_source", "none"))
    return {
        "candidate_model": "nested_locked_first_stage_c",
        "candidate_family": "locked_first_stage",
        "parameter_set_id": str(stage.get("parameter_set_id", "locked_first_stage")),
        "a_anchor": str(stage["a_anchor"]),
        "module_source": module_source,
        "use_modules": module_source != "none",
        "auxiliary_correlation": float(stage["auxiliary_correlation"]),
        "auxiliary_scale": float(stage["auxiliary_scale"]),
        "auxiliary_offset_sd": float(stage["auxiliary_offset_sd"]),
        "auxiliary_pairwise_scale": float(stage["auxiliary_pairwise_scale"]),
        "auxiliary_pairwise_density": float(stage["auxiliary_pairwise_density"]),
        "tau": float(stage["tau"]),
        "mode": str(stage.get("mode", "min")),
        "coefficient_seed": int(stage.get("coefficient_seed", 20260517)),
        "module_residual_scale": float(stage.get("module_residual_scale", 0.0)),
        "module_residual_density": float(stage.get("module_residual_density", 0.0)),
        "module_residual_orders": [int(x) for x in stage.get("module_residual_orders", [2, 3])],
        "source_run": str(stage.get("source_run", "")),
        "source_parameter_rank": stage.get("source_parameter_rank", ""),
    }


def fit_locked_first_stage_transform(
    *,
    config: dict[str, Any],
    graphlets: list[dict[str, Any]],
    output_dir: Path,
) -> dict[str, Any]:
    """Fit one locked first-stage transform over the experiment graphlet pool."""
    params = first_stage_parameters(config)
    density = float(params["auxiliary_pairwise_density"])
    all_a: list[float] = []
    all_raw: list[float] = []
    for graphlet in graphlets:
        a_values = np.asarray(graphlet[params["a_anchor"]], dtype=np.float64)
        raw_pairs = [
            auxiliary_raw_for_features(
                node_features,
                seed=int(params["coefficient_seed"]),
                pair_density=density,
            )
            for node_features in graphlet_node_feature_lists(tuple(graphlet["features"]))
        ]
        raw_values = np.asarray(
            [add + float(params["auxiliary_pairwise_scale"]) * pair for add, pair in raw_pairs],
            dtype=np.float64,
        )
        all_a.extend(a_values.tolist())
        all_raw.extend(raw_values.tolist())
    transform = fit_auxiliary_transform(
        a_values=np.asarray(all_a, dtype=np.float64),
        raw_values=np.asarray(all_raw, dtype=np.float64),
    )
    payload = {
        "source_run": params["source_run"],
        "source_parameter_rank": params["source_parameter_rank"],
        "a_anchor": params["a_anchor"],
        "transform_basis": "experiment_graphlet_node_pool",
        "n_graphlets": int(len(graphlets)),
        "n_graphlet_nodes": int(len(all_a)),
        "auxiliary_pairwise_density": density,
        "auxiliary_pairwise_scale": float(params["auxiliary_pairwise_scale"]),
        **transform,
    }
    write_json(output_dir / "locked_first_stage_transform.json", payload)
    return payload


def decode_second_parameters(
    raw: dict[str, float | int],
    *,
    config: dict[str, Any],
    family: str,
    parameter_set_id: str,
) -> dict[str, Any]:
    module_source = FAMILY_MODULE_SOURCES[family]
    use_modules = module_source != "none"
    uncoupled = family == "uncoupled_module"
    density = density_choice(raw=raw, config=config, prefix="second")
    return {
        "candidate_model": f"nested_softmin_{family}",
        "candidate_family": family,
        "parameter_set_id": parameter_set_id,
        "second_module_source": module_source,
        "second_use_modules": use_modules,
        "second_uncoupled_auxiliary": uncoupled,
        "second_auxiliary_correlation": 0.0
        if uncoupled
        else float(raw.get("second_auxiliary_correlation", 0.75)),
        "second_auxiliary_scale": float(raw.get("second_auxiliary_scale", 0.75)),
        "second_auxiliary_offset_sd": float(raw.get("second_auxiliary_offset_sd", 0.0)),
        "second_auxiliary_pairwise_scale": float(raw.get("second_auxiliary_pairwise_scale", 0.05)),
        "second_auxiliary_pairwise_density": density,
        "second_tau": tau_from_log(raw.get("second_log_tau", 0.0)),
        "second_log_tau": float(raw.get("second_log_tau", 0.0)),
        "second_mode": str(config["nested_operator"]["second_stage"].get("mode", "min")),
        "second_module_residual_scale": 0.0
        if not use_modules
        else float(raw.get("second_module_residual_scale", 0.25)),
        "second_module_residual_density": 0.0
        if not use_modules
        else float(raw.get("second_module_residual_density", 1.0)),
        "second_module_residual_orders": [int(x) for x in config["nested_operator"]["second_stage"]["module_residual_orders"]],
        "second_coefficient_seed": int(config["nested_operator"]["second_stage"].get("coefficient_seed", 20260518)),
        **{f"raw__{key}": value for key, value in raw.items()},
    }


def auxiliary_values_for_graphlet(
    *,
    base_values: np.ndarray,
    features: tuple[tuple[int, str], ...],
    modules: list[dict[str, Any]],
    params: dict[str, Any],
    prefix: str,
    graphlet_id: str,
    parameter_set_id: str,
    transform_override: dict[str, Any] | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    base = np.asarray(base_values, dtype=np.float64)
    pair_density = float(params[f"{prefix}_auxiliary_pairwise_density"])
    pair_scale = float(params[f"{prefix}_auxiliary_pairwise_scale"])
    coeff_seed = int(params.get(f"{prefix}_coefficient_seed", 20260517))
    raw_pairs = [
        auxiliary_raw_for_features(
            node_features,
            seed=coeff_seed,
            pair_density=pair_density,
        )
        for node_features in graphlet_node_feature_lists(features)
    ]
    raw_values = np.asarray([add + pair_scale * pair for add, pair in raw_pairs], dtype=np.float64)

    module_coefficients: dict[tuple[tuple[int, str], ...], float] = {}
    module_raw = np.zeros_like(raw_values, dtype=np.float64)
    module_scale = float(params.get(f"{prefix}_module_residual_scale", 0.0))
    if modules and module_scale > 0.0:
        module_seed = int(
            graphlet_runner.stable_seed(
                "nested-softmin-module",
                prefix,
                parameter_set_id,
                graphlet_id,
                coeff_seed,
            )
        )
        module_coefficients = module_residual_coefficients(
            modules,
            orders=[int(x) for x in params.get(f"{prefix}_module_residual_orders", [2, 3])],
            residual_density=float(params.get(f"{prefix}_module_residual_density", 1.0)),
            coefficient_seed=module_seed,
            coefficient_scale=1.0,
        )
        module_raw = module_residual_values_for_feature_sets(
            graphlet_node_feature_lists(features),
            module_coefficients,
        )
        raw_values = raw_values + module_scale * zscore(module_raw)

    transform = (
        transform_override
        if transform_override is not None
        else fit_auxiliary_transform(a_values=base, raw_values=raw_values)
    )
    auxiliary = apply_auxiliary_transform(
        a_values=base,
        raw_values=raw_values,
        transform=transform,
        rho=float(params[f"{prefix}_auxiliary_correlation"]),
        auxiliary_scale=float(params[f"{prefix}_auxiliary_scale"]),
        offset_sd=float(params[f"{prefix}_auxiliary_offset_sd"]),
    )
    info = {
        f"{prefix}_raw_sd": float(np.std(raw_values)),
        f"{prefix}_raw_to_base_correlation": safe_corr(raw_values, base),
        f"{prefix}_module_raw_sd": float(np.std(module_raw)),
        f"{prefix}_n_module_coefficients": int(len(module_coefficients)),
    }
    return auxiliary, info


def first_stage_for_graphlet(
    *,
    graphlet: dict[str, Any],
    config: dict[str, Any],
    modules_by_source: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    params = first_stage_parameters(config)
    a_values = np.asarray(graphlet[params["a_anchor"]], dtype=np.float64)
    stage_params = {
        "first_auxiliary_correlation": params["auxiliary_correlation"],
        "first_auxiliary_scale": params["auxiliary_scale"],
        "first_auxiliary_offset_sd": params["auxiliary_offset_sd"],
        "first_auxiliary_pairwise_scale": params["auxiliary_pairwise_scale"],
        "first_auxiliary_pairwise_density": params["auxiliary_pairwise_density"],
        "first_module_residual_scale": params["module_residual_scale"],
        "first_module_residual_density": params["module_residual_density"],
        "first_module_residual_orders": params["module_residual_orders"],
        "first_coefficient_seed": params["coefficient_seed"],
    }
    b_values, aux_info = auxiliary_values_for_graphlet(
        base_values=a_values,
        features=tuple(graphlet["features"]),
        modules=modules_by_source[str(params["module_source"])],
        params=stage_params,
        prefix="first",
        graphlet_id=str(graphlet["graphlet_id"]),
        parameter_set_id=str(params["parameter_set_id"]),
        transform_override=config.get("_runtime_first_stage_transform"),
    )
    c_values, first_innovation, first_weight = boltzmann_two_state_values(
        base_values=a_values,
        auxiliary_values=b_values,
        tau=float(params["tau"]),
        mode=str(params["mode"]),
    )
    dim = int(graphlet["graphlet_dimension"])
    extra = {
        "first_stage_source_run": params["source_run"],
        "first_stage_source_parameter_rank": params["source_parameter_rank"],
        "a_anchor": params["a_anchor"],
        "first_module_source": params["module_source"],
        "first_tau": params["tau"],
        "first_mode": params["mode"],
        "first_auxiliary_to_a_correlation": safe_corr(b_values, a_values),
        "c_to_a_correlation": safe_corr(c_values, a_values),
        "first_edge_gradient_correlation_to_a": safe_corr(
            node_edge_gradients(c_values, dim=dim),
            node_edge_gradients(a_values, dim=dim),
        ),
        "first_innovation_to_a_sd_ratio": float(np.std(first_innovation) / max(float(np.std(a_values)), 1e-12)),
        **prefixed(soft_state_metrics(first_weight, first_innovation, dim=dim), "first"),
        **aux_info,
    }
    return {
        "params": params,
        "a_values": a_values,
        "b_values": b_values,
        "c_values": c_values,
        "first_innovation": first_innovation,
        "first_weight": first_weight,
        "extra": extra,
    }


def evaluate_nested_graphlet(
    *,
    graphlet: dict[str, Any],
    config: dict[str, Any],
    modules_by_source: dict[str, list[dict[str, Any]]],
    second_params: dict[str, Any],
) -> dict[str, Any]:
    first = first_stage_for_graphlet(graphlet=graphlet, config=config, modules_by_source=modules_by_source)
    c_values = np.asarray(first["c_values"], dtype=np.float64)
    modules = modules_by_source[str(second_params["second_module_source"])]
    d_values, aux_info = auxiliary_values_for_graphlet(
        base_values=c_values,
        features=tuple(graphlet["features"]),
        modules=modules,
        params=second_params,
        prefix="second",
        graphlet_id=str(graphlet["graphlet_id"]),
        parameter_set_id=str(second_params["parameter_set_id"]),
    )
    values, second_innovation, second_weight = boltzmann_two_state_values(
        base_values=c_values,
        auxiliary_values=d_values,
        tau=float(second_params["second_tau"]),
        mode=str(second_params["second_mode"]),
    )
    a_values = np.asarray(first["a_values"], dtype=np.float64)
    dim = int(graphlet["graphlet_dimension"])
    total_innovation = values - a_values
    extra = {
        **{key: value for key, value in second_params.items() if key != "second_module_residual_orders"},
        "second_module_residual_orders": ";".join(str(x) for x in second_params["second_module_residual_orders"]),
        **first["extra"],
        "second_auxiliary_to_c_correlation": safe_corr(d_values, c_values),
        "values_to_c_correlation": safe_corr(values, c_values),
        "values_to_a_correlation": safe_corr(values, a_values),
        "second_edge_gradient_correlation_to_c": safe_corr(
            node_edge_gradients(values, dim=dim),
            node_edge_gradients(c_values, dim=dim),
        ),
        "second_innovation_to_c_sd_ratio": float(np.std(second_innovation) / max(float(np.std(c_values)), 1e-12)),
        "total_innovation_to_a_sd_ratio": float(np.std(total_innovation) / max(float(np.std(a_values)), 1e-12)),
        "second_innovation_edge_sd": float(np.std(edge_difference_values(second_innovation, dim=dim))),
        **prefixed(soft_state_metrics(second_weight, second_innovation, dim=dim), "second"),
        **aux_info,
    }
    return {
        **first,
        "d_values": d_values,
        "values": values,
        "second_innovation": second_innovation,
        "second_weight": second_weight,
        "total_innovation": total_innovation,
        "second_modules": modules,
        "extra": extra,
    }


def regularization_penalty(rows: list[dict[str, Any]], config: dict[str, Any]) -> tuple[float, dict[str, Any]]:
    cfg = dict(config["abc"].get("regularization", {}))
    if not rows:
        return 0.0, {"nested_regularization_penalty": 0.0}
    df = pd.DataFrame(rows)

    def mean_col(name: str) -> float:
        if name not in df:
            return math.nan
        values = pd.to_numeric(df[name], errors="coerce").dropna()
        return float(values.mean()) if not values.empty else math.nan

    mean_corr = mean_col("values_to_c_correlation")
    mean_edge_corr = mean_col("second_edge_gradient_correlation_to_c")
    mean_innov_ratio = mean_col("second_innovation_to_c_sd_ratio")
    mean_weight_sd = mean_col("second_auxiliary_weight_sd")

    penalties = {
        "penalty_low_values_to_c_correlation": float(
            cfg.get("lambda_pheno", 0.0)
            * max(0.0, float(cfg.get("values_to_c_correlation_min", -1.0)) - mean_corr) ** 2
            if np.isfinite(mean_corr)
            else 0.0
        ),
        "penalty_low_edge_gradient_correlation_to_c": float(
            cfg.get("lambda_edge", 0.0)
            * max(0.0, float(cfg.get("edge_gradient_correlation_to_c_min", -1.0)) - mean_edge_corr) ** 2
            if np.isfinite(mean_edge_corr)
            else 0.0
        ),
        "penalty_high_second_innovation_ratio": float(
            cfg.get("lambda_energy", 0.0)
            * max(0.0, mean_innov_ratio - float(cfg.get("second_innovation_to_c_sd_max", math.inf))) ** 2
            if np.isfinite(mean_innov_ratio)
            else 0.0
        ),
        "penalty_low_second_state_use": float(
            cfg.get("lambda_state_use", 0.0)
            * max(0.0, float(cfg.get("second_auxiliary_weight_sd_min", 0.0)) - mean_weight_sd) ** 2
            if np.isfinite(mean_weight_sd)
            else 0.0
        ),
    }
    total = float(np.nansum(list(penalties.values())))
    return total, {
        "nested_regularization_penalty": total,
        "mean_values_to_c_correlation": mean_corr,
        "mean_second_edge_gradient_correlation_to_c": mean_edge_corr,
        "mean_second_innovation_to_c_sd_ratio": mean_innov_ratio,
        "mean_second_auxiliary_weight_sd": mean_weight_sd,
        **penalties,
    }


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
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    active_rows: list[dict[str, Any]] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        nested = evaluate_nested_graphlet(
            graphlet=graphlet,
            config=config,
            modules_by_source=modules_by_source,
            second_params=params,
        )
        extra = dict(nested["extra"])
        record_metric_values(
            graphlet=graphlet,
            candidate_model=str(params["candidate_model"]),
            candidate_family=family,
            parameter_set_id=parameter_set_id,
            values=np.asarray(nested["values"], dtype=np.float64),
            score_column=primary_score,
            metric_rows=metric_rows,
            walsh_rows=walsh_rows,
            gft_rows=gft_rows,
            random_walk_cfg=random_walk_cfg,
            extra=extra,
        )
        active_rows.append(extra)

    model_rows = pd.DataFrame(metric_rows)
    distance, by_metric = metric_distance(model_rows, target_rows, metric_cols)
    penalty, guardrail = guardrail_penalty(by_metric, config)
    reg_penalty, reg_details = regularization_penalty(active_rows, config)
    if not np.isfinite(distance):
        distance = float("inf")
    objective = float(distance + penalty + reg_penalty)
    summary = [float(distance), penalty, reg_penalty] + [float(by_metric.get(metric, math.nan)) for metric in metric_cols]
    extras = {
        **params,
        "calibration_graphlet_distance": float(distance),
        "objective_distance": objective,
        **guardrail,
        **reg_details,
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
        ]:
            if col in active_df:
                values = pd.to_numeric(active_df[col], errors="coerce").dropna()
                if not values.empty:
                    extras[f"{col}_mean"] = float(values.mean())
    return SimulationResult(distance=objective, summary_vector=np.asarray(summary, dtype=np.float64), extras=extras)


def run_family_abc(
    *,
    family: str,
    config: dict[str, Any],
    output_dir: Path,
    graphlets: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
    target_rows: pd.DataFrame,
    random_walk_cfg: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    smc_cfg = dict(config["abc"]["smc"])
    backend_cfg = dict(config["abc"].get("backend", {}))
    specs = specs_for_family(config, family)
    seed = int(smc_cfg.get("seed", 0)) + int(graphlet_runner.stable_seed("nested-family", family) % 100000)
    smc = SMCABCConfig(
        n_particles=int(smc_cfg["n_particles"]),
        n_rounds=int(smc_cfg["n_rounds"]),
        initial_pool_multiplier=int(smc_cfg["initial_pool_multiplier"]),
        proposal_pool_multiplier=int(smc_cfg["proposal_pool_multiplier"]),
        epsilon_quantile=float(smc_cfg["epsilon_quantile"]),
        seed=seed,
    )
    backend = SMCABCBackendConfig(
        kind=str(backend_cfg.get("kind", "serial")),
        max_workers=int(backend_cfg.get("max_workers", 1)),
        batch_size=int(backend_cfg.get("batch_size", 1)),
    )
    checkpoint = SMCABCCheckpointConfig(
        path=str(output_dir / "checkpoints" / f"smc_nested_{family}.json"),
        resume=bool(config["abc"].get("resume_checkpoints", True)),
    )

    def simulate(raw: dict[str, float | int]) -> SimulationResult:
        raw_payload = json.dumps(raw, sort_keys=True)
        parameter_set_id = f"{family}__nested_abc_{graphlet_runner.stable_seed(family, raw_payload):010d}"
        return simulate_family_particle(
            raw=raw,
            family=family,
            parameter_set_id=parameter_set_id,
            config=config,
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=target_rows,
            random_walk_cfg=random_walk_cfg,
        )

    result = run_smc_abc(
        specs=specs,
        config=smc,
        simulate=simulate,
        backend=backend,
        checkpoint=checkpoint,
        progress_callback=progress_writer(output_dir),
        run_label=f"nested_abc_{family}",
    )
    particle_rows: list[dict[str, Any]] = []
    for rank, particle in enumerate(sorted(result.particles, key=lambda item: item.distance), start=1):
        particle_rows.append(
            {
                "candidate_family": family,
                "posterior_rank": int(rank),
                "weight": float(particle.weight),
                "objective_distance": float(particle.distance),
                "summary_vector": ";".join(
                    "nan" if not np.isfinite(float(x)) else f"{float(x):.12g}"
                    for x in particle.summary_vector
                ),
                **{f"raw__{key}": value for key, value in particle.parameters.items()},
                **particle.extras,
            }
        )
    round_rows = [
        {
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
            "candidate_family": family,
            "name": spec.name,
            "kind": spec.kind,
            "low": float(spec.low),
            "high": float(spec.high),
            "initial_scale_fraction": float(spec.initial_scale_fraction),
            "min_scale_fraction": float(spec.min_scale_fraction),
        }
        for spec in specs
    ]
    return particle_rows, round_rows, spec_rows


def selected_family_summaries(
    *,
    distances: list[dict[str, Any]],
    metric_cols: list[str],
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    if not distances:
        return [], [], []
    df = pd.DataFrame(distances)
    df["distance_to_empirical_primary"] = pd.to_numeric(df["distance_to_empirical_primary"], errors="coerce")
    selected: list[dict[str, Any]] = []
    families = [str(x) for x in config["abc"]["fit_families"]]
    static_families = [
        "locked_first_stage",
        "fitted_current",
        "smooth_core",
        "fitted_current_posterior_mean",
        "residual_only_diagnostic",
        "gnk_local_null",
        "empirical_readout_residual",
    ]
    for family in families + static_families:
        calib = df[(df["candidate_family"] == family) & (df["graphlet_split"] == "calibration")].copy()
        if calib.empty:
            continue
        calib = calib[np.isfinite(calib["distance_to_empirical_primary"])]
        if calib.empty:
            continue
        best = calib.sort_values("distance_to_empirical_primary").iloc[0]
        parameter_set_id = str(best["parameter_set_id"])
        for split in ["calibration", "test", "all"]:
            sub = df[
                (df["candidate_family"] == family)
                & (df["parameter_set_id"].astype(str) == parameter_set_id)
                & (df["graphlet_split"] == split)
            ]
            if sub.empty:
                continue
            row = sub.iloc[0].to_dict()
            row["selection_split"] = "calibration"
            selected.append(row)

    test_rows = pd.DataFrame([row for row in selected if row.get("graphlet_split") == "test"])
    pairwise: list[dict[str, Any]] = []
    if not test_rows.empty:
        values = {
            str(row["candidate_family"]): float(row["distance_to_empirical_primary"])
            for row in test_rows.to_dict(orient="records")
            if np.isfinite(float(row["distance_to_empirical_primary"]))
        }
        references = [
            "locked_first_stage",
            "fitted_current",
            "smooth_core",
            "no_module",
            "empirical_module",
            "random_module",
            "shuffled_module",
            "uncoupled_module",
        ]
        for family, value in values.items():
            for reference in references:
                if reference == family or reference not in values:
                    continue
                pairwise.append(
                    {
                        "graphlet_split": "test",
                        "candidate_family": family,
                        "reference_family": reference,
                        "candidate_distance": value,
                        "reference_distance": values[reference],
                        "delta_distance": value - values[reference],
                        "candidate_better": bool(value < values[reference]),
                    }
                )

    guardrail_rows: list[dict[str, Any]] = []
    for row in selected:
        if row.get("graphlet_split") != "test":
            continue
        by_metric = {metric: float(row.get(f"distance_component__{metric}", math.nan)) for metric in metric_cols}
        penalty, details = guardrail_penalty(by_metric, config)
        guardrail_rows.append(
            {
                "graphlet_split": "test",
                "candidate_family": row.get("candidate_family"),
                "candidate_model": row.get("candidate_model"),
                "parameter_set_id": row.get("parameter_set_id"),
                "distance_to_empirical_primary": row.get("distance_to_empirical_primary"),
                "guardrail_penalty": penalty,
                **details,
                **{f"distance_component__{metric}": by_metric[metric] for metric in metric_cols},
            }
        )
    return selected, pairwise, guardrail_rows


def record_all_outputs(
    *,
    config: dict[str, Any],
    output_dir: Path,
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    modules_by_source: dict[str, list[dict[str, Any]]],
    posterior_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    metric_cols = [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]
    random_walk_cfg = {"analysis": config.get("analysis", {})}
    empirical_innovation = context["empirical_innovation"]
    readout_base_anchor = str(config["empirical_innovation"]["base_anchor"])

    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    locality_rows: list[dict[str, Any]] = []
    locality_summary_rows: list[dict[str, Any]] = []
    correlation_rows: list[dict[str, Any]] = []

    candidate_particle_rows = []
    posterior_df = pd.DataFrame(posterior_rows)
    n_predictive = int(config["abc"].get("posterior_predictive_particles", 8))
    if not posterior_df.empty:
        for family, sub in posterior_df.groupby("candidate_family", dropna=False):
            sub = sub.sort_values("objective_distance").head(n_predictive)
            for rank, row in enumerate(sub.to_dict(orient="records"), start=1):
                raw = {
                    key.replace("raw__", "", 1): row[key]
                    for key in row
                    if str(key).startswith("raw__") and not str(key).startswith("raw__raw__")
                }
                parameter_set_id = f"{family}__posterior_rank{rank:03d}"
                decoded = decode_second_parameters(
                    raw,
                    config=config,
                    family=str(family),
                    parameter_set_id=parameter_set_id,
                )
                decoded["posterior_rank"] = int(rank)
                decoded["posterior_objective_distance"] = float(row["objective_distance"])
                candidate_particle_rows.append(decoded)

    first_params = first_stage_parameters(config)
    for graphlet in graphlets:
        empirical_values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        record_metric_values(
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
            extra={"a_anchor": "", "first_module_source": "", "second_module_source": ""},
        )

        base_for_empirical = np.asarray(graphlet[readout_base_anchor], dtype=np.float64)
        empirical_resid = empirical_innovation[str(graphlet["graphlet_id"])]
        loc_rows, loc_summary = energy_locality_for_modules(
            graphlet_id=str(graphlet["graphlet_id"]),
            graphlet_split=str(graphlet["graphlet_split"]),
            graphlet_features=tuple(graphlet["features"]),
            innovation=empirical_resid,
            modules=modules_by_source["empirical_calibration_plaquette"],
            candidate_model="empirical_global_readout_residual",
            parameter_set_id=f"{primary_score}__minus_poly_{readout_base_anchor}",
        )
        locality_rows.extend(loc_rows)
        locality_summary_rows.append(loc_summary)
        correlation_rows.append(
            correlation_row(
                graphlet=graphlet,
                candidate_model="empirical_global_readout_residual",
                parameter_set_id=f"{primary_score}__minus_poly_{readout_base_anchor}",
                base_values=base_for_empirical,
                values=empirical_values,
                innovation=empirical_resid,
            )
        )
        record_metric_values(
            graphlet=graphlet,
            candidate_model="empirical_global_readout_residual",
            candidate_family="empirical_readout_residual",
            parameter_set_id=f"{primary_score}__minus_poly_{readout_base_anchor}",
            values=empirical_resid,
            score_column=primary_score,
            metric_rows=metric_rows,
            walsh_rows=walsh_rows,
            gft_rows=gft_rows,
            random_walk_cfg=random_walk_cfg,
            extra={"a_anchor": readout_base_anchor, "first_module_source": "empirical_calibration_plaquette", "second_module_source": ""},
        )

        for candidate_model, family in BASELINE_FAMILIES.items():
            if candidate_model in {"gnk_local_null", "empirical_global_readout_residual"}:
                continue
            values = np.asarray(graphlet[candidate_model], dtype=np.float64)
            record_metric_values(
                graphlet=graphlet,
                candidate_model=candidate_model,
                candidate_family=family,
                parameter_set_id=f"{ASSAY_ID}__{candidate_model}",
                values=values,
                score_column=primary_score,
                metric_rows=metric_rows,
                walsh_rows=walsh_rows,
                gft_rows=gft_rows,
                random_walk_cfg=random_walk_cfg,
                extra={"a_anchor": candidate_model, "first_module_source": "baseline_control", "second_module_source": ""},
            )

        gnk = gnk_values(
            dim=int(graphlet["graphlet_dimension"]),
            seed=graphlet_runner.stable_seed("phot-nested-softmin-gnk", graphlet["graphlet_id"]),
            target_mean=float(np.mean(empirical_values)),
            target_sd=float(np.std(empirical_values)),
            k=2,
        )
        record_metric_values(
            graphlet=graphlet,
            candidate_model="gnk_local_null",
            candidate_family="gnk_local_null",
            parameter_set_id=f"{ASSAY_ID}__gnk_k2",
            values=gnk,
            score_column=primary_score,
            metric_rows=metric_rows,
            walsh_rows=walsh_rows,
            gft_rows=gft_rows,
            random_walk_cfg=random_walk_cfg,
            extra={"a_anchor": "", "first_module_source": "gnk_control", "second_module_source": ""},
        )

        first = first_stage_for_graphlet(graphlet=graphlet, config=config, modules_by_source=modules_by_source)
        first_extra = {
            **first["extra"],
            "candidate_model": first_params["candidate_model"],
            "candidate_family": first_params["candidate_family"],
            "parameter_set_id": first_params["parameter_set_id"],
            "second_module_source": "",
        }
        record_metric_values(
            graphlet=graphlet,
            candidate_model=str(first_params["candidate_model"]),
            candidate_family=str(first_params["candidate_family"]),
            parameter_set_id=str(first_params["parameter_set_id"]),
            values=np.asarray(first["c_values"], dtype=np.float64),
            score_column=primary_score,
            metric_rows=metric_rows,
            walsh_rows=walsh_rows,
            gft_rows=gft_rows,
            random_walk_cfg=random_walk_cfg,
            extra=first_extra,
        )
        loc_rows, loc_summary = energy_locality_for_modules(
            graphlet_id=str(graphlet["graphlet_id"]),
            graphlet_split=str(graphlet["graphlet_split"]),
            graphlet_features=tuple(graphlet["features"]),
            innovation=np.asarray(first["first_innovation"], dtype=np.float64),
            modules=modules_by_source[str(first_params["module_source"])],
            candidate_model=str(first_params["candidate_model"]),
            parameter_set_id=str(first_params["parameter_set_id"]),
        )
        for row in loc_rows:
            row.update(first_extra)
        loc_summary.update(first_extra)
        locality_rows.extend(loc_rows)
        locality_summary_rows.append(loc_summary)
        corr = correlation_row(
            graphlet=graphlet,
            candidate_model=str(first_params["candidate_model"]),
            parameter_set_id=str(first_params["parameter_set_id"]),
            base_values=np.asarray(first["a_values"], dtype=np.float64),
            values=np.asarray(first["c_values"], dtype=np.float64),
            innovation=np.asarray(first["first_innovation"], dtype=np.float64),
        )
        corr.update(first_extra)
        correlation_rows.append(corr)

        for params in candidate_particle_rows:
            nested = evaluate_nested_graphlet(
                graphlet=graphlet,
                config=config,
                modules_by_source=modules_by_source,
                second_params=params,
            )
            extra = dict(nested["extra"])
            record_metric_values(
                graphlet=graphlet,
                candidate_model=str(params["candidate_model"]),
                candidate_family=str(params["candidate_family"]),
                parameter_set_id=str(params["parameter_set_id"]),
                values=np.asarray(nested["values"], dtype=np.float64),
                score_column=primary_score,
                metric_rows=metric_rows,
                walsh_rows=walsh_rows,
                gft_rows=gft_rows,
                random_walk_cfg=random_walk_cfg,
                extra=extra,
            )
            loc_rows, loc_summary = energy_locality_for_modules(
                graphlet_id=str(graphlet["graphlet_id"]),
                graphlet_split=str(graphlet["graphlet_split"]),
                graphlet_features=tuple(graphlet["features"]),
                innovation=np.asarray(nested["second_innovation"], dtype=np.float64),
                modules=nested["second_modules"],
                candidate_model=str(params["candidate_model"]),
                parameter_set_id=str(params["parameter_set_id"]),
            )
            for row in loc_rows:
                row.update(extra)
            loc_summary.update(extra)
            locality_rows.extend(loc_rows)
            locality_summary_rows.append(loc_summary)
            corr = correlation_row(
                graphlet=graphlet,
                candidate_model=str(params["candidate_model"]),
                parameter_set_id=str(params["parameter_set_id"]),
                base_values=np.asarray(nested["c_values"], dtype=np.float64),
                values=np.asarray(nested["values"], dtype=np.float64),
                innovation=np.asarray(nested["second_innovation"], dtype=np.float64),
            )
            corr.update(extra)
            correlation_rows.append(corr)

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
            ["candidate_family", "candidate_model", "parameter_set_id"],
            dropna=False,
        ):
            candidate_family, candidate_model, parameter_set_id = keys
            if str(candidate_model).startswith("empirical_observed__"):
                continue
            distance, by_metric = metric_distance(model_rows, target_rows, metric_cols)
            first = model_rows.iloc[0].to_dict()
            distances.append(
                {
                    "assay_id": ASSAY_ID,
                    "graphlet_split": split,
                    "candidate_family": candidate_family,
                    "candidate_model": candidate_model,
                    "parameter_set_id": parameter_set_id,
                    "a_anchor": first.get("a_anchor", ""),
                    "first_module_source": first.get("first_module_source", ""),
                    "second_module_source": first.get("second_module_source", ""),
                    "posterior_rank": first.get("posterior_rank", ""),
                    "distance_to_empirical_primary": distance,
                    **{f"distance_component__{key}": value for key, value in by_metric.items()},
                }
            )

    selected_rows, pairwise_rows, guardrail_rows = selected_family_summaries(
        distances=distances,
        metric_cols=metric_cols,
        config=config,
    )
    summary_rows = summarize_metric_rows(
        metric_rows,
        group_cols=["graphlet_split", "candidate_family", "candidate_model", "parameter_set_id", "score_column"],
        metric_cols=metric_cols
        + [
            "second_innovation_to_c_sd_ratio",
            "total_innovation_to_a_sd_ratio",
            "second_auxiliary_weight_sd",
            "plaquette_innovation_total_energy",
            "edge_innovation_total_energy",
        ],
    )
    module_summary = aggregate_module_energy(locality_rows)
    top_rows = top_distance_rows(distances)

    write_csv(output_dir / "graphlet_catalog.csv", graphlet_catalog_rows(graphlets))
    write_csv(output_dir / "module_catalog.csv", module_catalog_rows(modules_by_source))
    write_csv(output_dir / "graphlet_metrics.csv", metric_rows)
    write_csv(output_dir / "walsh_terms.csv", walsh_rows)
    write_csv(output_dir / "gft_coefficients.csv", gft_rows)
    write_csv(output_dir / "graphlet_summary.csv", summary_rows)
    write_csv(output_dir / "model_distances.csv", distances)
    write_csv(output_dir / "top_model_distances.csv", top_rows)
    write_csv(output_dir / "selected_family_distances.csv", selected_rows)
    write_csv(output_dir / "selected_family_pairwise_deltas.csv", pairwise_rows)
    write_csv(output_dir / "component_guardrail_summary.csv", guardrail_rows)
    write_csv(output_dir / "module_energy_rows.csv", locality_rows)
    write_csv(output_dir / "module_energy_by_graphlet.csv", locality_summary_rows)
    write_csv(output_dir / "module_energy_summary.csv", module_summary)
    write_csv(output_dir / "base_correlation.csv", correlation_rows)
    write_csv(output_dir / "evaluated_parameter_sets.csv", candidate_particle_rows)
    write_json(output_dir / "empirical_readout_fit.json", context["readout"])

    return {
        "top_model_distances": top_rows[:20],
        "selected_family_distances": selected_rows,
        "selected_family_pairwise_deltas": pairwise_rows,
    }


def write_preflight(
    *,
    config: dict[str, Any],
    output_dir: Path,
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    modules_by_source: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    first = first_stage_parameters(config)
    preflight = {
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_result": config["source_run"].get("result_id", "RES-030"),
        "source_output_dir": context["source_output_dir"],
        "mode": str(config.get("experiment", {}).get("mode", "nested_softmin_ensemble_capacity_diagnostic")),
        "primary_score_column": config["graphlet_evaluation"]["primary_score_column"],
        "score_columns": context["assay_score_columns"],
        "objective_metric_columns": [str(x) for x in config["graphlet_evaluation"]["metric_columns"]],
        "local_ruggedness_metric_columns": [
            metric
            for metric in LOCAL_RUGGEDNESS_METRIC_COLUMNS
            if metric in [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]
        ],
        "fit_families": [str(x) for x in config["abc"]["fit_families"]],
        "n_graphlets": int(len(graphlets)),
        "n_calibration_graphlets": int(sum(1 for graphlet in graphlets if graphlet["graphlet_split"] == "calibration")),
        "n_test_graphlets": int(sum(1 for graphlet in graphlets if graphlet["graphlet_split"] == "test")),
        "graphlet_dimensions": sorted({int(graphlet["graphlet_dimension"]) for graphlet in graphlets}),
        "module_counts": {source: int(len(modules)) for source, modules in modules_by_source.items()},
        "readout_base_anchor": config["empirical_innovation"]["base_anchor"],
        "run077_phot_best_objective": context["source_phot_best_objective"],
        "first_stage_locked_parameters": first,
        "first_stage_locked_transform": config.get("_runtime_first_stage_transform", {}),
    }
    write_json(output_dir / "preflight.json", preflight)
    return preflight


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
    config["_runtime_first_stage_transform"] = fit_locked_first_stage_transform(
        config=config,
        graphlets=graphlets,
        output_dir=output_dir,
    )
    target_rows = calibration_target_rows(
        graphlets=graphlets,
        config=config,
        random_walk_cfg=random_walk_cfg,
    )
    preflight = write_preflight(
        config=config,
        output_dir=output_dir,
        graphlets=graphlets,
        context=context,
        modules_by_source=modules_by_source,
    )

    posterior_rows: list[dict[str, Any]] = []
    round_rows: list[dict[str, Any]] = []
    spec_rows: list[dict[str, Any]] = []
    for family in [str(x) for x in config["abc"]["fit_families"]]:
        particles, rounds, specs = run_family_abc(
            family=family,
            config=config,
            output_dir=output_dir,
            graphlets=graphlets,
            modules_by_source=modules_by_source,
            target_rows=target_rows,
            random_walk_cfg=random_walk_cfg,
        )
        posterior_rows.extend(particles)
        round_rows.extend(rounds)
        spec_rows.extend(specs)
        write_csv(output_dir / "posterior_particles.csv", posterior_rows)
        write_csv(output_dir / "posterior_rounds.csv", round_rows)
        write_csv(output_dir / "parameter_specs.csv", spec_rows)

    output_summary = record_all_outputs(
        config=config,
        output_dir=output_dir,
        graphlets=graphlets,
        context=context,
        modules_by_source=modules_by_source,
        posterior_rows=posterior_rows,
    )

    write_json(
        output_dir / "summary.json",
        {
            "experiment_id": config["experiment"]["id"],
            "preflight": preflight,
            "output_files": {
                "preflight": str(output_dir / "preflight.json"),
                "progress": str(output_dir / "progress.json"),
                "locked_first_stage_transform": str(output_dir / "locked_first_stage_transform.json"),
                "parameter_specs": str(output_dir / "parameter_specs.csv"),
                "posterior_particles": str(output_dir / "posterior_particles.csv"),
                "posterior_rounds": str(output_dir / "posterior_rounds.csv"),
                "evaluated_parameter_sets": str(output_dir / "evaluated_parameter_sets.csv"),
                "graphlet_catalog": str(output_dir / "graphlet_catalog.csv"),
                "module_catalog": str(output_dir / "module_catalog.csv"),
                "graphlet_metrics": str(output_dir / "graphlet_metrics.csv"),
                "walsh_terms": str(output_dir / "walsh_terms.csv"),
                "gft_coefficients": str(output_dir / "gft_coefficients.csv"),
                "graphlet_summary": str(output_dir / "graphlet_summary.csv"),
                "model_distances": str(output_dir / "model_distances.csv"),
                "top_model_distances": str(output_dir / "top_model_distances.csv"),
                "selected_family_distances": str(output_dir / "selected_family_distances.csv"),
                "selected_family_pairwise_deltas": str(output_dir / "selected_family_pairwise_deltas.csv"),
                "component_guardrail_summary": str(output_dir / "component_guardrail_summary.csv"),
                "module_energy_rows": str(output_dir / "module_energy_rows.csv"),
                "module_energy_by_graphlet": str(output_dir / "module_energy_by_graphlet.csv"),
                "module_energy_summary": str(output_dir / "module_energy_summary.csv"),
                "base_correlation": str(output_dir / "base_correlation.csv"),
                "empirical_readout_fit": str(output_dir / "empirical_readout_fit.json"),
            },
            **output_summary,
        },
    )


if __name__ == "__main__":
    main()
