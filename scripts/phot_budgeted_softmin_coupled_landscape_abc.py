#!/usr/bin/env python3
"""PHOT budgeted module-tied soft-min coupled-landscape SMC-ABC test.

This HYP-008 experiment keeps the RUN-077 PHOT base landscape fixed and fits a
small two-state operator family against calibration graphlets. The fitted object
is the state-difference field ``d(x)=b(x)-a(x)``, with module-local coefficients
selected from calibration empirical residual structure and explicit
finite-difference energy/state-boundary budgets.
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
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_ROOT, ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import tem1_aggregate_graphlet_residual_ruggedness as graphlet_runner  # noqa: E402
from empirical_graphlet_falsification_suite import (  # noqa: E402
    find_complete_graphlets,
    metric_rows_for_values,
    score_map_for_assay,
    subset_masks,
)
from phot_hard_min_coupled_landscape_grid_search import (  # noqa: E402
    ASSAY_ID,
    affine_fit,
    apply_affine,
    apply_readout_model,
    assay_prediction_metrics,
    deterministic_mask,
    feature_label,
    fit_readout_model,
    graphlet_energy_concentration,
    load_config,
    locality_covariance_summary,
    nrmse,
    parse_feature_lists,
    rank_spearman,
    resolve_source_output_dir,
    stable_unit,
    summarize_metric_rows,
    write_csv,
    write_json,
)
from project_adapt_env.coupled_landscape import (  # noqa: E402
    budget_penalty_terms,
    control_modules,
    edge_energy_by_bit,
    energy_by_graphlet_feature,
    finite_difference_budget_metrics,
    module_tied_difference_values,
    node_edge_gradients,
    plaquette_energy_by_pair,
    safe_corr,
    select_empirical_modules,
    softmin_from_difference,
    soft_state_metrics,
)
from project_adapt_env.proteingym_panel import merge_dict  # noqa: E402
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
    evaluate_landscape_on_graphlet,
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


def parameter_specs(config: dict[str, Any]) -> list[ParameterSpec]:
    specs: list[ParameterSpec] = []
    for item in config["abc"]["parameters"]:
        specs.append(
            ParameterSpec(
                name=str(item["name"]),
                kind=str(item["kind"]),
                low=float(item["low"]),
                high=float(item["high"]),
                initial_scale_fraction=float(item.get("initial_scale_fraction", 0.15)),
                min_scale_fraction=float(item.get("min_scale_fraction", 0.05)),
            )
        )
    return specs


def tau_from_parameters(parameters: dict[str, float | int]) -> float:
    return float(math.exp(float(parameters["log_tau"])))


def parameter_payload(
    parameters: dict[str, float | int],
    *,
    config: dict[str, Any],
    label: str,
) -> dict[str, Any]:
    return {
        "candidate_model": "budgeted_module_tied_softmin",
        "parameter_set_id": label,
        "operator": "budgeted_module_tied_softmin",
        "a_anchor": str(config["operator"]["a_anchor"]),
        "tau": tau_from_parameters(parameters),
        "module_residual_orders": ";".join(str(x) for x in config["module_difference"].get("module_residual_orders", [2])),
        **{key: value for key, value in parameters.items()},
    }


def posterior_mean_parameters(
    *,
    specs: list[ParameterSpec],
    particle_rows: list[dict[str, Any]],
) -> dict[str, float | int]:
    weights = np.asarray([float(row["weight"]) for row in particle_rows], dtype=np.float64)
    total = float(np.sum(weights))
    weights = weights / total if total > 0 and np.isfinite(total) else np.full(len(particle_rows), 1.0 / len(particle_rows))
    out: dict[str, float | int] = {}
    for spec in specs:
        values = np.asarray([float(row[spec.name]) for row in particle_rows], dtype=np.float64)
        mean = float(np.sum(weights * values))
        if spec.kind == "int":
            out[spec.name] = int(np.clip(round(mean), int(spec.low), int(spec.high)))
        else:
            out[spec.name] = float(np.clip(mean, spec.low, spec.high))
    return out


def ks_statistic(xs: np.ndarray, ys: np.ndarray) -> float:
    x = np.sort(np.asarray(xs, dtype=np.float64)[np.isfinite(xs)])
    y = np.sort(np.asarray(ys, dtype=np.float64)[np.isfinite(ys)])
    if x.size == 0 or y.size == 0:
        return math.nan
    values = np.sort(np.unique(np.concatenate([x, y])))
    cdf_x = np.searchsorted(x, values, side="right") / float(x.size)
    cdf_y = np.searchsorted(y, values, side="right") / float(y.size)
    return float(np.max(np.abs(cdf_x - cdf_y)))


def assay_guardrail_distance(
    *,
    observed: np.ndarray,
    predicted: np.ndarray,
    holdout: np.ndarray,
) -> float:
    train = ~holdout
    intercept, slope = affine_fit(predicted[train], observed[train])
    mapped = apply_affine(predicted, intercept, slope)
    return float(
        np.nanmean(
            [
                1.0 - rank_spearman(observed[holdout], predicted[holdout]),
                nrmse(observed[holdout], mapped[holdout]),
                ks_statistic(observed, mapped),
            ]
        )
    )


def graphlet_node_feature_lists(features: tuple[tuple[int, str], ...]) -> list[tuple[tuple[int, str], ...]]:
    return [tuple(sorted(mask, key=lambda item: (item[0], item[1]))) for mask in subset_masks(features)]


def softmin_prediction(
    *,
    a_values: np.ndarray,
    feature_sets: list[tuple[tuple[int, str], ...]],
    modules: list[dict[str, Any]],
    parameters: dict[str, float | int],
    config: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    difference_params = {
        **parameters,
        "use_modules": bool(config["module_difference"].get("use_modules", True)),
        "module_residual_orders": [int(x) for x in config["module_difference"].get("module_residual_orders", [2])],
    }
    difference_values, difference_info, _coeffs = module_tied_difference_values(
        base_values=np.asarray(a_values, dtype=np.float64),
        feature_sets=feature_sets,
        modules=modules,
        parameter_set=difference_params,
        coefficient_seed=int(config["operator"].get("coefficient_seed", 20260517)),
    )
    values, innovation, aux_weight, b_values = softmin_from_difference(
        base_values=np.asarray(a_values, dtype=np.float64),
        difference_values=difference_values,
        tau=tau_from_parameters(parameters),
    )
    return values, innovation, aux_weight, b_values, difference_values, difference_info


def graphlet_extra_metrics(
    *,
    graphlet: dict[str, Any],
    values: np.ndarray,
    innovation: np.ndarray,
    aux_weight: np.ndarray,
    b_values: np.ndarray,
    difference_values: np.ndarray,
    difference_info: dict[str, float],
    budget_metrics: dict[str, float] | None,
    params_payload: dict[str, Any],
) -> dict[str, Any]:
    dim = int(graphlet["graphlet_dimension"])
    a_values = np.asarray(graphlet[params_payload["a_anchor"]], dtype=np.float64)
    extra = dict(params_payload)
    extra.update(
        {
            "auxiliary_to_a_correlation": safe_corr(b_values, a_values),
            "difference_to_a_correlation": safe_corr(difference_values, a_values),
            "values_to_a_correlation": safe_corr(values, a_values),
            "edge_gradient_correlation_to_a": safe_corr(
                node_edge_gradients(values, dim=dim),
                node_edge_gradients(a_values, dim=dim),
            ),
            "innovation_to_a_sd_ratio": float(np.std(innovation) / max(float(np.std(a_values)), 1e-12)),
            **difference_info,
            **(budget_metrics or {}),
            **soft_state_metrics(aux_weight, innovation, dim=dim),
        }
    )
    return extra


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=args.quick)
    source_output = resolve_source_output_dir(config, project_root)
    system_dir = source_output / "systems" / ASSAY_ID
    fit_dir = system_dir / "fit"

    fit_summary = json.loads((fit_dir / "summary.json").read_text(encoding="utf-8"))
    source_summary = json.loads((source_output / "summary.json").read_text(encoding="utf-8"))
    panel = pd.read_csv(fit_dir / "selected_panel.csv")
    wildtype = str(panel.loc[panel["DMS_id"] == ASSAY_ID, "target_seq"].iloc[0])
    assay_path = fit_dir / "checkpoints" / "panel_preparation" / ASSAY_ID / "prepared_assay.csv"
    metrics_path = fit_dir / "checkpoints" / "panel_preparation" / ASSAY_ID / "mavenn_metrics.json"
    full_assay_df = pd.read_csv(assay_path)
    assay_df = full_assay_df.copy()
    if int(config["assay_fit"].get("max_variants", 0)) > 0:
        limit = int(config["assay_fit"]["max_variants"])
        assay_df = assay_df.assign(_sample_key=assay_df["mutant"].map(lambda value: stable_unit("quick-assay", value)))
        assay_df = assay_df.sort_values("_sample_key").head(limit).drop(columns=["_sample_key"]).reset_index(drop=True)
    mavenn_metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    score_columns = [str(x) for x in config["graphlet_evaluation"]["score_columns"]]
    wt_scores = {
        "latent_phi": float(mavenn_metrics.get("wildtype_phi", math.nan)),
        "mavenn_yhat": float(mavenn_metrics.get("wildtype_yhat", math.nan)),
    }
    assay_score_columns = [col for col in score_columns if col in full_assay_df.columns or col in wt_scores]
    if primary_score not in assay_score_columns:
        raise ValueError(f"Primary score column {primary_score!r} is not available in PHOT assay.")

    alignment_path = system_dir / f"{ASSAY_ID}.alignment.fasta"
    if not alignment_path.exists():
        candidate = Path(str(fit_summary.get("panel", [{}])[0].get("alignment_path", "")))
        if candidate.exists():
            alignment_path = candidate
    if not alignment_path.exists():
        raise FileNotFoundError(f"Alignment FASTA not found for PHOT: {alignment_path}")
    shutil.copy2(alignment_path, output_dir / f"{ASSAY_ID}.alignment.fasta")

    base_config = dict(config["adapt_env"]["base_config"])
    seed = int(base_config.get("seed", 42))
    best_updates = dict(fit_summary["branches"]["smc_abc_best_raw"]["fitted_parameters"])
    mean_updates = dict(fit_summary["branches"]["smc_abc_posterior_mean_raw"]["fitted_parameters"])
    fitted_current = build_landscape(
        base_config=base_config,
        updates=best_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed,
    )
    fitted_mean = build_landscape(
        base_config=base_config,
        updates=mean_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed + 1,
    )
    smooth_core = build_landscape(
        base_config=base_config,
        updates=best_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed,
        smooth_core_only=True,
    )

    feature_lists = parse_feature_lists(assay_df, wildtype)
    observed = pd.to_numeric(assay_df[primary_score], errors="coerce").to_numpy(dtype=np.float64)
    sequences = assay_df["mutated_sequence"].astype(str).tolist()
    chunk_size = int(config["assay_fit"].get("evaluation_chunk_size", 4096))
    smooth_pred = evaluate_sequences_cached(smooth_core, sequences, chunk_size=chunk_size)
    current_pred = evaluate_sequences_cached(fitted_current, sequences, chunk_size=chunk_size)
    mean_pred = evaluate_sequences_cached(fitted_mean, sequences, chunk_size=chunk_size)
    wt_smooth = float(evaluate_sequences_cached(smooth_core, [wildtype], chunk_size=1)[0])
    wt_current = float(evaluate_sequences_cached(fitted_current, [wildtype], chunk_size=1)[0])
    wt_mean = float(evaluate_sequences_cached(fitted_mean, [wildtype], chunk_size=1)[0])
    wt_observed = float(wt_scores[primary_score])

    a_anchor = str(config["operator"]["a_anchor"])
    a_assay = {
        "fitted_smooth_core_only": smooth_pred,
        "fitted_current_best": current_pred,
    }[a_anchor]
    wt_a = {
        "fitted_smooth_core_only": wt_smooth,
        "fitted_current_best": wt_current,
    }[a_anchor]

    coeff_seed = int(config["operator"].get("coefficient_seed", 20260517))

    score_map = score_map_for_assay(
        full_assay_df,
        wildtype=wildtype,
        wt_scores=wt_scores,
        score_columns=assay_score_columns,
    )
    complete = []
    graphlet_cfg = config["phot_graphlets"]
    for dim in [int(x) for x in graphlet_cfg["dimensions"]]:
        complete.extend(
            find_complete_graphlets(
                score_map,
                score_columns=assay_score_columns,
                dimensions=[dim],
                max_per_dimension=int(graphlet_cfg["max_graphlets_per_dimension"]),
                top_feature_pool=int(graphlet_cfg["top_feature_pool"]),
            )
        )

    graphlets: list[dict[str, Any]] = []
    for graphlet_index, (dim, features) in enumerate(complete):
        split = "calibration" if graphlet_index % 2 == 0 else "test"
        masks = subset_masks(features)
        empirical = {
            score_col: np.asarray([score_map[mask][score_col] for mask in masks], dtype=np.float64)
            for score_col in assay_score_columns
        }
        smooth_values = evaluate_landscape_on_graphlet(smooth_core, wildtype=wildtype, features=features)
        current_values = evaluate_landscape_on_graphlet(fitted_current, wildtype=wildtype, features=features)
        mean_values = evaluate_landscape_on_graphlet(fitted_mean, wildtype=wildtype, features=features)
        graphlets.append(
            {
                "graphlet_id": f"{ASSAY_ID}__g{graphlet_index:04d}__d{dim}",
                "assay_id": ASSAY_ID,
                "graphlet_split": split,
                "split": split,
                "graphlet_dimension": int(dim),
                "features": tuple(features),
                "mutation_features": ";".join(feature_label(feature) for feature in features),
                "empirical_values": empirical,
                "fitted_smooth_core_only": smooth_values,
                "fitted_current_best": current_values,
                "fitted_current_posterior_mean": mean_values,
                "fitted_current_vs_smooth_innovation": current_values - smooth_values,
            }
        )

    calib_base: list[float] = []
    calib_empirical: list[float] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        calib_base.extend(np.asarray(graphlet[a_anchor], dtype=np.float64).tolist())
        calib_empirical.extend(np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64).tolist())
    readout = fit_readout_model(
        calib_base,
        calib_empirical,
        degree=int(config["empirical_local_signal"].get("global_readout_polynomial_degree", 3)),
    )
    empirical_residual_by_graphlet: dict[str, np.ndarray] = {}
    for graphlet in graphlets:
        a_values = np.asarray(graphlet[a_anchor], dtype=np.float64)
        empirical_values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        empirical_residual_by_graphlet[str(graphlet["graphlet_id"])] = empirical_values - apply_readout_model(
            a_values,
            readout,
        )
    module_cfg = config["module_selection"]
    modules = select_empirical_modules(
        graphlets,
        empirical_residual_by_graphlet,
        max_modules=int(module_cfg["max_modules"]),
        module_size=int(module_cfg["module_size"]),
    )
    feature_pool = sorted(
        {feature for graphlet in graphlets for feature in tuple(graphlet["features"])},
        key=lambda item: (int(item[0]), str(item[1])),
    )
    shuffled_modules = control_modules(
        modules,
        feature_pool,
        source="shuffled_module",
        rng=np.random.RandomState(graphlet_runner.stable_seed("budgeted-softmin-shuffled-modules", ASSAY_ID)),
    )
    random_modules = control_modules(
        modules,
        feature_pool,
        source="random_module",
        rng=np.random.RandomState(graphlet_runner.stable_seed("budgeted-softmin-random-modules", ASSAY_ID)),
    )
    modules_by_source = {
        "empirical_calibration_plaquette": modules,
        "shuffled_module": shuffled_modules,
        "random_module": random_modules,
        "none": [],
    }

    random_walk_cfg = {"analysis": config.get("analysis", {})}
    metric_cols = [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]

    def metric_for_graphlet(
        *,
        graphlet: dict[str, Any],
        candidate_model: str,
        parameter_set_id: str,
        values: np.ndarray,
        extra: dict[str, Any] | None = None,
    ) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
        record_graphlet = {
            "graphlet_id": graphlet["graphlet_id"],
            "assay_id": ASSAY_ID,
            "graphlet_split": graphlet["graphlet_split"],
            "split": graphlet["split"],
            "graphlet_dimension": graphlet["graphlet_dimension"],
            "mutation_features": graphlet["mutation_features"],
            "score_column": primary_score,
        }
        metric, walsh, gft, _access = metric_rows_for_values(
            graphlet=record_graphlet,
            values=np.asarray(values, dtype=np.float64),
            candidate_model=candidate_model,
            parameter_set_id=parameter_set_id,
            cfg=random_walk_cfg,
        )
        if extra:
            metric.update(extra)
            for row in walsh:
                row.update(extra)
            for row in gft:
                row.update(extra)
        return metric, walsh, gft

    target_metric_rows: list[dict[str, Any]] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        metric, _walsh, _gft = metric_for_graphlet(
            graphlet=graphlet,
            candidate_model=f"empirical_observed__{primary_score}",
            parameter_set_id=primary_score,
            values=np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64),
            extra={"operator": "empirical", "a_anchor": ""},
        )
        target_metric_rows.append(metric)
    target_metric_df = pd.DataFrame(target_metric_rows)

    guard_limit = int(config["abc"].get("assay_guardrail_max_variants", 8000))
    guard_df = assay_df.reset_index(names="_assay_index")
    if guard_limit > 0 and len(assay_df) > guard_limit:
        guard_df = guard_df.assign(_guard_key=guard_df["mutant"].map(lambda value: stable_unit("abc-guardrail", value)))
        guard_df = guard_df.sort_values("_guard_key").head(guard_limit).drop(columns=["_guard_key"]).reset_index(drop=True)
    guard_index = guard_df["_assay_index"].to_numpy(dtype=np.int64)
    guard_observed = pd.to_numeric(guard_df[primary_score], errors="coerce").to_numpy(dtype=np.float64)
    guard_holdout = deterministic_mask(
        guard_df["mutant"].astype(str).tolist(),
        holdout_fraction=float(config["assay_fit"].get("holdout_fraction", 0.2)),
        salt="softmin-abc-guardrail",
    )

    cache: dict[tuple[tuple[str, float | int], ...], SimulationResult] = {}
    specs = parameter_specs(config)

    def simulate(parameters: dict[str, float | int]) -> SimulationResult:
        key = tuple(
            (spec.name, round(float(parameters[spec.name]), 8) if spec.kind == "float" else int(parameters[spec.name]))
            for spec in specs
        )
        if key in cache:
            return cache[key]
        model_rows: list[dict[str, Any]] = []
        corr_values: list[float] = []
        edge_corr_values: list[float] = []
        innovation_ratios: list[float] = []
        weight_sd_values: list[float] = []
        budget_accumulator: dict[str, list[float]] = {}
        payload = parameter_payload(parameters, config=config, label="abc_candidate")
        for graphlet in graphlets:
            if graphlet["graphlet_split"] != "calibration":
                continue
            a_values = np.asarray(graphlet[a_anchor], dtype=np.float64)
            feature_sets = graphlet_node_feature_lists(tuple(graphlet["features"]))
            values, innovation, aux_weight, b_values, difference_values, difference_info = softmin_prediction(
                a_values=a_values,
                feature_sets=feature_sets,
                modules=modules,
                parameters=parameters,
                config=config,
            )
            budget_metrics = finite_difference_budget_metrics(
                base_values=a_values,
                innovation=innovation,
                auxiliary_weight=aux_weight,
                difference_values=difference_values,
                tau=tau_from_parameters(parameters),
                dim=int(graphlet["graphlet_dimension"]),
                empirical_residual=empirical_residual_by_graphlet[str(graphlet["graphlet_id"])],
            )
            extra = graphlet_extra_metrics(
                graphlet=graphlet,
                values=values,
                innovation=innovation,
                aux_weight=aux_weight,
                b_values=b_values,
                difference_values=difference_values,
                difference_info=difference_info,
                budget_metrics=budget_metrics,
                params_payload=payload,
            )
            metric, _walsh, _gft = metric_for_graphlet(
                graphlet=graphlet,
                candidate_model="budgeted_module_tied_softmin",
                parameter_set_id="abc_candidate",
                values=values,
                extra=extra,
            )
            model_rows.append(metric)
            corr_values.append(float(extra["values_to_a_correlation"]))
            edge_corr_values.append(float(extra["edge_gradient_correlation_to_a"]))
            innovation_ratios.append(float(extra["innovation_to_a_sd_ratio"]))
            weight_sd_values.append(float(extra["auxiliary_weight_sd"]))
            for metric_name, metric_value in budget_metrics.items():
                budget_accumulator.setdefault(metric_name, []).append(float(metric_value))

        graphlet_distance, by_metric = metric_distance(pd.DataFrame(model_rows), target_metric_df, metric_cols)
        guard_feature_sets = [
            tuple(sorted(feature_lists[int(idx)], key=lambda item: (int(item[0]), str(item[1]))))
            for idx in guard_index
        ]
        guard_pred, _guard_innov, _guard_weight, _guard_b, _guard_d, _guard_info = softmin_prediction(
            a_values=a_assay[guard_index],
            feature_sets=guard_feature_sets,
            modules=modules,
            parameters=parameters,
            config=config,
        )
        guard_distance = assay_guardrail_distance(
            observed=guard_observed,
            predicted=guard_pred,
            holdout=guard_holdout,
        )
        reg_cfg = config["abc"].get("regularization", {})
        mean_corr = float(np.nanmean(corr_values))
        mean_edge_corr = float(np.nanmean(edge_corr_values))
        mean_innovation_ratio = float(np.nanmean(innovation_ratios))
        mean_weight_sd = float(np.nanmean(weight_sd_values))
        mean_budget_metrics = {
            f"mean_{metric_name}": float(np.nanmean(values))
            for metric_name, values in budget_accumulator.items()
        }
        budget_inputs = {
            metric_name: float(np.nanmean(values))
            for metric_name, values in budget_accumulator.items()
        }
        penalties = {
            "penalty_low_phenotype_correlation": float(
                reg_cfg.get("lambda_pheno", 0.0)
                * max(0.0, float(reg_cfg.get("phenotype_correlation_min", 0.0)) - mean_corr) ** 2
            ),
            "penalty_low_edge_correlation": float(
                reg_cfg.get("lambda_edge", 0.0)
                * max(0.0, float(reg_cfg.get("edge_gradient_correlation_min", 0.0)) - mean_edge_corr) ** 2
            ),
            "penalty_high_innovation_ratio": float(
                reg_cfg.get("lambda_energy", 0.0)
                * max(0.0, mean_innovation_ratio - float(reg_cfg.get("innovation_to_base_sd_max", 1.0))) ** 2
            ),
            "penalty_low_state_transition": float(
                reg_cfg.get("lambda_state_use", 0.0)
                * max(0.0, float(reg_cfg.get("auxiliary_weight_sd_min", 0.0)) - mean_weight_sd) ** 2
            ),
            **budget_penalty_terms(budget_inputs, reg_cfg),
        }
        total_penalty = float(np.nansum(list(penalties.values())))
        assay_weight = float(config["abc"].get("assay_guardrail_weight", 0.15))
        distance = float(graphlet_distance + assay_weight * guard_distance + total_penalty)
        summary_vector = np.asarray(
            [
                graphlet_distance,
                guard_distance,
                mean_corr,
                mean_edge_corr,
                mean_innovation_ratio,
                mean_weight_sd,
                total_penalty,
                *mean_budget_metrics.values(),
                *[by_metric.get(metric, math.nan) for metric in metric_cols],
            ],
            dtype=np.float64,
        )
        result = SimulationResult(
            distance=distance,
            summary_vector=summary_vector,
            extras={
                "graphlet_distance": graphlet_distance,
                "assay_guardrail_distance": guard_distance,
                "mean_values_to_a_correlation": mean_corr,
                "mean_edge_gradient_correlation_to_a": mean_edge_corr,
                "mean_innovation_to_a_sd_ratio": mean_innovation_ratio,
                "mean_auxiliary_weight_sd": mean_weight_sd,
                "total_penalty": total_penalty,
                **penalties,
                **mean_budget_metrics,
                **{f"distance_component__{key}": value for key, value in by_metric.items()},
            },
        )
        cache[key] = result
        return result

    progress_events: list[dict[str, Any]] = []
    smc_result = run_smc_abc(
        specs=specs,
        config=SMCABCConfig(**config["abc"]["smc"]),
        simulate=simulate,
        backend=SMCABCBackendConfig(**config["abc"].get("backend", {"kind": "serial", "max_workers": 1})),
        checkpoint=SMCABCCheckpointConfig(
            path=str((output_dir / "checkpoints" / "smc_budgeted_softmin_abc.json").resolve()),
            resume=True,
        ),
        progress_callback=lambda event: progress_events.append(dict(event)),
        run_label="phot_budgeted_softmin_coupled_landscape_abc",
    )

    particle_rows: list[dict[str, Any]] = []
    for rank, particle in enumerate(sorted(smc_result.particles, key=lambda item: item.distance), start=1):
        particle_rows.append(
            {
                "rank": rank,
                "distance": float(particle.distance),
                "weight": float(particle.weight),
                **{name: particle.parameters[name] for name in sorted(particle.parameters)},
                **particle.extras,
            }
        )
    round_rows = [
        {
            "round_index": item.round_index,
            "epsilon": item.epsilon,
            "attempts": item.attempts,
            "best_distance": item.best_distance,
            "median_distance": item.median_distance,
        }
        for item in smc_result.round_summaries
    ]
    mean_params = posterior_mean_parameters(specs=specs, particle_rows=particle_rows)

    evaluated_parameters: list[tuple[str, dict[str, float | int]]] = [
        ("abc_best_particle", dict(smc_result.best_particle.parameters)),
        ("abc_posterior_mean", mean_params),
    ]
    max_eval = int(config["abc"].get("posterior_predictive_particles", len(particle_rows)))
    for idx, row in enumerate(particle_rows[:max_eval], start=1):
        evaluated_parameters.append(
            (
                f"abc_particle_rank_{idx:03d}",
                {spec.name: int(row[spec.name]) if spec.kind == "int" else float(row[spec.name]) for spec in specs},
            )
        )

    neutral_threshold = float(config["assay_fit"].get("neutral_effect_threshold", 0.05))
    holdout_fraction = float(config["assay_fit"].get("holdout_fraction", 0.2))
    assay_fit_rows: list[dict[str, Any]] = []
    dfe_rows: list[dict[str, Any]] = []
    epistasis_rows: list[dict[str, Any]] = []

    baseline_predictions = [
        ("fitted_smooth_core_only", f"{ASSAY_ID}__smooth_core_A", smooth_pred, wt_smooth, {"operator": "none", "a_anchor": "fitted_smooth_core_only"}),
        ("fitted_current_best", f"{ASSAY_ID}__current_best", current_pred, wt_current, {"operator": "none", "a_anchor": "fitted_current_best"}),
        ("fitted_current_posterior_mean", f"{ASSAY_ID}__posterior_mean", mean_pred, wt_mean, {"operator": "none", "a_anchor": "fitted_current_posterior_mean"}),
    ]
    for candidate_model, parameter_set_id, pred, wt_pred, params in baseline_predictions:
        fit_row, dfe_row, epi_row = assay_prediction_metrics(
            assay_df=assay_df,
            feature_lists=feature_lists,
            observed=observed,
            pred=pred,
            wt_observed=wt_observed,
            wt_pred=wt_pred,
            candidate_model=candidate_model,
            parameter_set_id=parameter_set_id,
            params=params,
            neutral_threshold=neutral_threshold,
            holdout_fraction=holdout_fraction,
        )
        assay_fit_rows.append(fit_row)
        dfe_rows.append(dfe_row)
        epistasis_rows.append(epi_row)

    for label, parameters in evaluated_parameters:
        assay_feature_sets = list(feature_lists)
        values, innovation, aux_weight, b_pred, d_pred, d_info = softmin_prediction(
            a_values=a_assay,
            feature_sets=assay_feature_sets,
            modules=modules,
            parameters=parameters,
            config=config,
        )
        wt_value, wt_innovation, wt_weight, wt_b, wt_d, wt_info = softmin_prediction(
            a_values=np.asarray([wt_a], dtype=np.float64),
            feature_sets=[tuple()],
            modules=modules,
            parameters=parameters,
            config=config,
        )
        payload = parameter_payload(parameters, config=config, label=label)
        payload.update(
            {
                "assay_auxiliary_weight_mean": float(np.mean(aux_weight)),
                "assay_auxiliary_weight_sd": float(np.std(aux_weight)),
                "assay_auxiliary_to_a_correlation": safe_corr(b_pred, a_assay),
                "assay_difference_to_a_correlation": safe_corr(d_pred, a_assay),
                "assay_innovation_to_a_correlation": safe_corr(innovation, a_assay),
                "assay_innovation_to_a_sd_ratio": float(np.std(innovation) / max(float(np.std(a_assay)), 1e-12)),
                **{f"assay_{key}": value for key, value in d_info.items()},
            }
        )
        fit_row, dfe_row, epi_row = assay_prediction_metrics(
            assay_df=assay_df,
            feature_lists=feature_lists,
            observed=observed,
            pred=values,
            wt_observed=wt_observed,
            wt_pred=float(wt_value[0]),
            candidate_model="budgeted_module_tied_softmin",
            parameter_set_id=label,
            params=payload,
            neutral_threshold=neutral_threshold,
            holdout_fraction=holdout_fraction,
        )
        assay_fit_rows.append(fit_row)
        dfe_rows.append(dfe_row)
        epistasis_rows.append(epi_row)

    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    local_signal_rows: list[dict[str, Any]] = []
    feature_energy_rows: list[dict[str, Any]] = []

    def record_values(
        *,
        graphlet: dict[str, Any],
        candidate_model: str,
        parameter_set_id: str,
        values: np.ndarray,
        extra: dict[str, Any] | None = None,
    ) -> None:
        metric, walsh, gft = metric_for_graphlet(
            graphlet=graphlet,
            candidate_model=candidate_model,
            parameter_set_id=parameter_set_id,
            values=values,
            extra=extra,
        )
        metric_rows.append(metric)
        walsh_rows.extend(walsh)
        gft_rows.extend(gft)

    for graphlet in graphlets:
        empirical_values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        a_values = np.asarray(graphlet[a_anchor], dtype=np.float64)
        current_values = np.asarray(graphlet["fitted_current_best"], dtype=np.float64)
        empirical_pred = apply_readout_model(a_values, readout)
        empirical_resid = empirical_values - empirical_pred
        record_values(
            graphlet=graphlet,
            candidate_model=f"empirical_observed__{primary_score}",
            parameter_set_id=primary_score,
            values=empirical_values,
            extra={"operator": "empirical", "a_anchor": ""},
        )
        record_values(
            graphlet=graphlet,
            candidate_model="empirical_global_readout_residual",
            parameter_set_id=f"{primary_score}__minus_poly_{a_anchor}",
            values=empirical_resid,
            extra={"operator": "empirical_residual", "a_anchor": a_anchor},
        )
        for candidate_model in [
            "fitted_smooth_core_only",
            "fitted_current_best",
            "fitted_current_posterior_mean",
            "fitted_current_vs_smooth_innovation",
        ]:
            values = np.asarray(graphlet[candidate_model], dtype=np.float64)
            record_values(
                graphlet=graphlet,
                candidate_model=candidate_model,
                parameter_set_id=f"{ASSAY_ID}__{candidate_model}",
                values=values,
                extra={"operator": "baseline_control", "a_anchor": candidate_model},
            )
        gnk = gnk_values(
            dim=int(graphlet["graphlet_dimension"]),
            seed=graphlet_runner.stable_seed("phot-softmin-abc-gnk", graphlet["graphlet_id"]),
            target_mean=float(np.mean(empirical_values)),
            target_sd=float(np.std(empirical_values)),
            k=2,
        )
        record_values(
            graphlet=graphlet,
            candidate_model="gnk_local_null",
            parameter_set_id=f"{ASSAY_ID}__gnk_k2",
            values=gnk,
            extra={"operator": "gnk_control", "a_anchor": ""},
        )
        feature_energy_rows.extend(
            energy_by_graphlet_feature(
                graphlet_id=str(graphlet["graphlet_id"]),
                graphlet_split=str(graphlet["graphlet_split"]),
                graphlet_features=tuple(graphlet["features"]),
                innovation=empirical_resid,
                candidate_model="empirical_global_readout_residual",
                parameter_set_id=f"{primary_score}__minus_poly_{a_anchor}",
            )
        )
        local_signal_rows.append(
            {
                "graphlet_id": graphlet["graphlet_id"],
                "graphlet_split": graphlet["graphlet_split"],
                "mutation_features": graphlet["mutation_features"],
                "candidate_model": "empirical_global_readout_residual",
                "parameter_set_id": f"{primary_score}__minus_poly_{a_anchor}",
                "operator": "empirical_residual",
                "edge_innovation_total_energy": float(sum(edge_energy_by_bit(empirical_resid, dim=int(graphlet["graphlet_dimension"])).values())),
                "plaquette_innovation_total_energy": float(sum(plaquette_energy_by_pair(empirical_resid, dim=int(graphlet["graphlet_dimension"])).values())),
                "soft_switch_edge_fraction": math.nan,
                "soft_switch_plaquette_fraction": math.nan,
                "soft_switch_boundary_plaquette_energy_fraction": math.nan,
                "auxiliary_weight_mean": math.nan,
            }
        )
        current_innovation = current_values - np.asarray(graphlet["fitted_smooth_core_only"], dtype=np.float64)
        feature_energy_rows.extend(
            energy_by_graphlet_feature(
                graphlet_id=str(graphlet["graphlet_id"]),
                graphlet_split=str(graphlet["graphlet_split"]),
                graphlet_features=tuple(graphlet["features"]),
                innovation=current_innovation,
                candidate_model="fitted_current_vs_smooth_innovation",
                parameter_set_id=f"{ASSAY_ID}__current_minus_smooth",
            )
        )
        local_signal_rows.append(
            {
                "graphlet_id": graphlet["graphlet_id"],
                "graphlet_split": graphlet["graphlet_split"],
                "mutation_features": graphlet["mutation_features"],
                "candidate_model": "fitted_current_vs_smooth_innovation",
                "parameter_set_id": f"{ASSAY_ID}__current_minus_smooth",
                "operator": "baseline_current_minus_smooth",
                "edge_innovation_total_energy": float(sum(edge_energy_by_bit(current_innovation, dim=int(graphlet["graphlet_dimension"])).values())),
                "plaquette_innovation_total_energy": float(sum(plaquette_energy_by_pair(current_innovation, dim=int(graphlet["graphlet_dimension"])).values())),
                "soft_switch_edge_fraction": math.nan,
                "soft_switch_plaquette_fraction": math.nan,
                "soft_switch_boundary_plaquette_energy_fraction": math.nan,
                "auxiliary_weight_mean": math.nan,
            }
        )

        for label, parameters in evaluated_parameters:
            module_evaluations: list[tuple[str, str, list[dict[str, Any]], str]] = [
                ("budgeted_module_tied_softmin", label, modules, "empirical_calibration_plaquette")
            ]
            if label in {"abc_best_particle", "abc_posterior_mean"}:
                module_evaluations.extend(
                    [
                        ("budgeted_softmin_no_module_control", f"{label}__no_module", modules_by_source["none"], "none"),
                        ("budgeted_softmin_shuffled_module_control", f"{label}__shuffled_module", shuffled_modules, "shuffled_module"),
                        ("budgeted_softmin_random_module_control", f"{label}__random_module", random_modules, "random_module"),
                    ]
                )
            for candidate_model, parameter_set_id, active_modules, module_source in module_evaluations:
                payload = parameter_payload(parameters, config=config, label=parameter_set_id)
                payload["module_source"] = module_source
                feature_sets = graphlet_node_feature_lists(tuple(graphlet["features"]))
                values, innovation, aux_weight, b_values, difference_values, difference_info = softmin_prediction(
                    a_values=a_values,
                    feature_sets=feature_sets,
                    modules=active_modules,
                    parameters=parameters,
                    config=config,
                )
                budget_metrics = finite_difference_budget_metrics(
                    base_values=a_values,
                    innovation=innovation,
                    auxiliary_weight=aux_weight,
                    difference_values=difference_values,
                    tau=tau_from_parameters(parameters),
                    dim=int(graphlet["graphlet_dimension"]),
                    empirical_residual=empirical_resid,
                )
                extra = graphlet_extra_metrics(
                    graphlet=graphlet,
                    values=values,
                    innovation=innovation,
                    aux_weight=aux_weight,
                    b_values=b_values,
                    difference_values=difference_values,
                    difference_info=difference_info,
                    budget_metrics=budget_metrics,
                    params_payload=payload,
                )
                record_values(
                    graphlet=graphlet,
                    candidate_model=candidate_model,
                    parameter_set_id=parameter_set_id,
                    values=values,
                    extra=extra,
                )
                feature_energy_rows.extend(
                    energy_by_graphlet_feature(
                        graphlet_id=str(graphlet["graphlet_id"]),
                        graphlet_split=str(graphlet["graphlet_split"]),
                        graphlet_features=tuple(graphlet["features"]),
                        innovation=innovation,
                        candidate_model=candidate_model,
                        parameter_set_id=parameter_set_id,
                    )
                )
                local_signal_rows.append(
                    {
                        "graphlet_id": graphlet["graphlet_id"],
                        "graphlet_split": graphlet["graphlet_split"],
                        "mutation_features": graphlet["mutation_features"],
                        "candidate_model": candidate_model,
                        "parameter_set_id": parameter_set_id,
                        **extra,
                        "edge_innovation_total_energy": float(sum(edge_energy_by_bit(innovation, dim=int(graphlet["graphlet_dimension"])).values())),
                        "plaquette_innovation_total_energy": float(sum(plaquette_energy_by_pair(innovation, dim=int(graphlet["graphlet_dimension"])).values())),
                    }
                )

    metric_df = pd.DataFrame(metric_rows)
    distances: list[dict[str, Any]] = []
    for split in ["calibration", "test", "all"]:
        split_df = metric_df if split == "all" else metric_df[metric_df["graphlet_split"] == split]
        target_rows = split_df[
            (split_df["candidate_model"] == f"empirical_observed__{primary_score}")
        ]
        if target_rows.empty:
            continue
        for keys, model_rows in split_df.groupby(["candidate_model", "parameter_set_id"], dropna=False):
            candidate_model, parameter_set_id = keys
            if str(candidate_model).startswith("empirical_observed__"):
                continue
            distance, by_metric = metric_distance(model_rows, target_rows, metric_cols)
            first = model_rows.iloc[0].to_dict()
            distances.append(
                {
                    "assay_id": ASSAY_ID,
                    "graphlet_split": split,
                    "candidate_model": candidate_model,
                    "parameter_set_id": parameter_set_id,
                    "distance_to_empirical_primary": distance,
                    **{
                        key: first.get(key, "")
                        for key in [
                            "operator",
                            "a_anchor",
                            "module_source",
                            "tau",
                            "log_tau",
                            "difference_bias_sd",
                            "difference_base_slope",
                            "difference_module_scale",
                            "module_residual_density",
                            "module_residual_orders",
                            "innovation_var_ratio_to_base",
                            "innovation_edge_ratio_to_empirical_residual",
                            "innovation_plaquette_ratio_to_empirical_residual",
                            "difference_edge_median_abs_over_tau",
                        ]
                    },
                    **{f"distance_component__{key}": value for key, value in by_metric.items()},
                }
            )

    summary_rows = summarize_metric_rows(
        metric_rows,
        group_cols=["graphlet_split", "candidate_model", "parameter_set_id", "score_column"],
        metric_cols=metric_cols
        + [
            "soft_switch_edge_fraction",
            "soft_switch_plaquette_fraction",
            "soft_switch_boundary_plaquette_energy_fraction",
            "auxiliary_weight_mean",
            "auxiliary_weight_sd",
            "innovation_to_a_sd_ratio",
            "values_to_a_correlation",
            "edge_gradient_correlation_to_a",
            "difference_to_base_sd_ratio",
            "difference_to_a_correlation",
            "innovation_var_ratio_to_base",
            "innovation_edge_ratio_to_empirical_residual",
            "innovation_plaquette_ratio_to_empirical_residual",
            "state_weight_edge_energy",
            "state_weight_plaquette_energy",
            "difference_edge_median_abs_over_tau",
        ],
    )
    catalog_rows = [
        {
            "graphlet_id": graphlet["graphlet_id"],
            "assay_id": ASSAY_ID,
            "graphlet_split": graphlet["graphlet_split"],
            "graphlet_dimension": graphlet["graphlet_dimension"],
            "mutation_features": graphlet["mutation_features"],
        }
        for graphlet in graphlets
    ]
    feature_energy_summary = aggregate_feature_energy_soft(feature_energy_rows)
    graphlet_energy_summary = graphlet_energy_concentration(local_signal_rows)
    covariance_summary = locality_covariance_summary_soft(local_signal_rows)

    top_rows: list[dict[str, Any]] = []
    distance_df = pd.DataFrame(distances)
    if not distance_df.empty:
        for split in ["calibration", "test", "all"]:
            sub = distance_df[distance_df["graphlet_split"] == split].copy()
            sub["distance_to_empirical_primary"] = pd.to_numeric(sub["distance_to_empirical_primary"], errors="coerce")
            sub = sub[np.isfinite(sub["distance_to_empirical_primary"])]
            if sub.empty:
                continue
            for rank, row in enumerate(sub.sort_values("distance_to_empirical_primary").head(20).to_dict(orient="records"), start=1):
                top_rows.append({"rank": rank, **row})

    assay_rank_rows: list[dict[str, Any]] = []
    assay_fit_df = pd.DataFrame(assay_fit_rows)
    if not assay_fit_df.empty:
        assay_fit_df["assay_fit_objective_distance"] = pd.to_numeric(
            assay_fit_df["assay_fit_objective_distance"],
            errors="coerce",
        )
        for rank, row in enumerate(
            assay_fit_df.sort_values("assay_fit_objective_distance").head(30).to_dict(orient="records"),
            start=1,
        ):
            assay_rank_rows.append({"rank": rank, **row})

    module_rows = [
        {
            "module_id": str(module["module_id"]),
            "module_source": str(module.get("module_source", "empirical_calibration_plaquette")),
            "rank": module.get("rank", math.nan),
            "calibration_pair_score": module.get("calibration_pair_score", math.nan),
            "module_features": ";".join(feature_label(feature) for feature in module["features"]),
        }
        for module in modules
    ]

    spec_rows = [spec.__dict__ for spec in specs]
    write_csv(output_dir / "parameter_specs.csv", spec_rows)
    write_csv(output_dir / "module_catalog.csv", module_rows)
    write_csv(output_dir / "posterior_particles.csv", particle_rows)
    write_csv(output_dir / "posterior_rounds.csv", round_rows)
    write_csv(output_dir / "evaluated_parameter_sets.csv", [parameter_payload(params, config=config, label=label) for label, params in evaluated_parameters])
    write_csv(output_dir / "assay_fit_metrics.csv", assay_fit_rows)
    write_csv(output_dir / "assay_dfe_metrics.csv", dfe_rows)
    write_csv(output_dir / "assay_epistasis_metrics.csv", epistasis_rows)
    write_csv(output_dir / "assay_fit_rankings.csv", assay_rank_rows)
    write_csv(output_dir / "graphlet_catalog.csv", catalog_rows)
    write_csv(output_dir / "graphlet_metrics.csv", metric_rows)
    write_csv(output_dir / "walsh_terms.csv", walsh_rows)
    write_csv(output_dir / "gft_coefficients.csv", gft_rows)
    write_csv(output_dir / "graphlet_summary.csv", summary_rows)
    write_csv(output_dir / "model_distances.csv", distances)
    write_csv(output_dir / "top_model_distances.csv", top_rows)
    write_csv(output_dir / "local_signal_by_graphlet.csv", local_signal_rows)
    write_csv(output_dir / "local_signal_energy_by_feature.csv", feature_energy_rows)
    write_csv(output_dir / "local_signal_feature_energy_summary.csv", feature_energy_summary)
    write_csv(output_dir / "local_signal_graphlet_energy_summary.csv", graphlet_energy_summary)
    write_csv(output_dir / "locality_covariance_summary.csv", covariance_summary)
    write_json(output_dir / "empirical_readout_fit.json", readout)
    write_json(output_dir / "progress.json", progress_events)

    source_phot_best = math.nan
    for system_summary in source_summary.get("systems", []):
        if str(system_summary.get("assay_id")) == ASSAY_ID:
            source_phot_best = float(system_summary["fit_summary"]["best_distance"])
            break
    preflight = {
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_output_dir": str(source_output),
        "source_results": ["RES-030", "RES-031", "RES-032", "RES-033"],
        "mode": "frozen_base_budgeted_module_tied_softmin_smc_abc",
        "operator": "budgeted module-tied soft-min over d(x)=b(x)-a(x)",
        "a_anchor": a_anchor,
        "n_modules": int(len(modules)),
        "module_selection": dict(module_cfg),
        "regularization": dict(config["abc"].get("regularization", {})),
        "primary_score_column": primary_score,
        "score_columns": assay_score_columns,
        "n_particles": int(config["abc"]["smc"]["n_particles"]),
        "n_rounds": int(config["abc"]["smc"]["n_rounds"]),
        "n_posterior_evaluated": int(len(evaluated_parameters)),
        "n_assay_variants_evaluated": int(len(assay_df)),
        "n_assay_guardrail_variants": int(len(guard_df)),
        "n_graphlets": int(len(graphlets)),
        "n_calibration_graphlets": int(sum(g["graphlet_split"] == "calibration" for g in graphlets)),
        "n_test_graphlets": int(sum(g["graphlet_split"] == "test" for g in graphlets)),
        "graphlet_dimensions": sorted({int(row["graphlet_dimension"]) for row in catalog_rows}),
        "run077_phot_best_objective": source_phot_best,
    }
    write_json(output_dir / "preflight.json", preflight)
    write_json(
        output_dir / "summary.json",
        {
            "experiment_id": config["experiment"]["id"],
            "preflight": preflight,
            "smc_best_particle": particle_rows[0] if particle_rows else {},
            "smc_rounds": round_rows,
            "output_files": {
                "parameter_specs": str(output_dir / "parameter_specs.csv"),
                "module_catalog": str(output_dir / "module_catalog.csv"),
                "posterior_particles": str(output_dir / "posterior_particles.csv"),
                "posterior_rounds": str(output_dir / "posterior_rounds.csv"),
                "evaluated_parameter_sets": str(output_dir / "evaluated_parameter_sets.csv"),
                "assay_fit_metrics": str(output_dir / "assay_fit_metrics.csv"),
                "graphlet_metrics": str(output_dir / "graphlet_metrics.csv"),
                "model_distances": str(output_dir / "model_distances.csv"),
                "local_signal_by_graphlet": str(output_dir / "local_signal_by_graphlet.csv"),
                "locality_covariance_summary": str(output_dir / "locality_covariance_summary.csv"),
            },
            "top_assay_fit_rankings": assay_rank_rows[:20],
            "top_graphlet_distances": top_rows[:20],
        },
    )


def evaluate_sequences_cached(landscape: Any, sequences: list[str], *, chunk_size: int) -> np.ndarray:
    from phot_hard_min_coupled_landscape_grid_search import evaluate_sequences

    return evaluate_sequences(landscape, sequences, chunk_size=chunk_size)


def aggregate_feature_energy_soft(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    from phot_hard_min_coupled_landscape_grid_search import aggregate_feature_energy

    return aggregate_feature_energy(rows)


def locality_covariance_summary_soft(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    metrics = [
        "edge_innovation_total_energy",
        "plaquette_innovation_total_energy",
        "soft_switch_edge_fraction",
        "soft_switch_plaquette_fraction",
    ]
    for keys, sub in df.groupby(["candidate_model", "parameter_set_id", "graphlet_split"], dropna=False):
        if len(sub) < 3:
            continue
        if not isinstance(keys, tuple):
            keys = (keys,)
        feature_sets = [
            set(str(item).split(";")) if str(item) else set()
            for item in sub["mutation_features"].tolist()
        ]
        shared_counts: list[float] = []
        out_row: dict[str, Any] = {
            "candidate_model": keys[0],
            "parameter_set_id": keys[1],
            "graphlet_split": keys[2],
            "n_graphlets": int(len(sub)),
        }
        metric_values = {metric: pd.to_numeric(sub.get(metric), errors="coerce").to_numpy(dtype=np.float64) for metric in metrics}
        diffs: dict[str, list[float]] = {metric: [] for metric in metrics}
        for i in range(len(sub)):
            for j in range(i + 1, len(sub)):
                shared_counts.append(float(len(feature_sets[i] & feature_sets[j])))
                for metric in metrics:
                    vi = metric_values[metric][i]
                    vj = metric_values[metric][j]
                    diffs[metric].append(-abs(float(vi - vj)) if np.isfinite(vi) and np.isfinite(vj) else math.nan)
        shared = np.asarray(shared_counts, dtype=np.float64)
        for metric, values in diffs.items():
            arr = np.asarray(values, dtype=np.float64)
            ok = np.isfinite(shared) & np.isfinite(arr)
            out_row[f"shared_feature_count_vs_{metric}_similarity_spearman"] = (
                rank_spearman(shared[ok], arr[ok]) if int(np.sum(ok)) >= 3 else math.nan
            )
        out.append(out_row)
    return out


if __name__ == "__main__":
    main()
