#!/usr/bin/env python3
"""PHOT coupled-family ABC diagnostic with held-out graphlet controls."""

from __future__ import annotations

import argparse
import json
import math
import os
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
from phot_coupled_landscape_grid_search import (  # noqa: E402
    apply_readout_model,
    correlation_row,
    feature_label,
    fit_readout_model,
    graphlet_split,
    load_config,
    resolve_source_output_dir,
    summarize_metric_rows,
    write_csv,
    write_json,
)
from project_adapt_env.coupled_landscape import (  # noqa: E402
    aggregate_module_energy,
    control_modules,
    coupled_values,
    energy_locality_for_modules,
    plaquette_energy_by_pair,
    safe_corr,
    select_empirical_modules,
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
    evaluate_landscape_on_graphlet,
    gnk_values,
    metric_distance,
)


ASSAY_ID = "PHOT_CHLRE_Chen_2023"

FAMILY_MODULE_SOURCES = {
    "no_module": "none",
    "empirical_module": "empirical_calibration_plaquette",
    "random_module": "random_module",
    "shuffled_module": "shuffled_module",
    "uncoupled_module": "uncoupled_empirical_module",
}

BASELINE_FAMILIES = {
    "fitted_current_best": "fitted_current",
    "fitted_smooth_core_only": "smooth_core",
    "fitted_current_posterior_mean": "fitted_current_posterior_mean",
    "fitted_residual_only_diagnostic": "residual_only_diagnostic",
    "gnk_local_null": "gnk_local_null",
    "empirical_global_readout_residual": "empirical_readout_residual",
}

LOCAL_RUGGEDNESS_METRIC_COLUMNS = [
    "edge_energy_adjacent_moran_i",
    "edge_energy_adjacent_pearson",
    "edge_energy_top10_fraction",
    "edge_energy_effective_fraction",
    "node_incident_top10_fraction",
    "node_incident_effective_fraction",
    "site_energy_top1_fraction",
    "site_energy_effective_fraction",
    "plaquette_energy_top1_fraction",
    "plaquette_energy_effective_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def parameter_specs_from_rows(rows: list[dict[str, Any]]) -> list[ParameterSpec]:
    specs: list[ParameterSpec] = []
    for row in rows:
        specs.append(
            ParameterSpec(
                name=str(row["name"]),
                kind=str(row["kind"]),
                low=float(row["low"]),
                high=float(row["high"]),
                initial_scale_fraction=float(row.get("initial_scale_fraction", 0.15)),
                min_scale_fraction=float(row.get("min_scale_fraction", 0.05)),
            )
        )
    return specs


def specs_for_family(config: dict[str, Any], family: str) -> list[ParameterSpec]:
    priors = dict(config["abc"]["parameter_priors"])
    names = list(config["abc"]["family_parameter_names"].get(family, []))
    rows = [{**priors[name], "name": name} for name in names]
    return parameter_specs_from_rows(rows)


def hypercube_edges(dim: int) -> list[dict[str, int]]:
    edges: list[dict[str, int]] = []
    for source in range(1 << int(dim)):
        for bit in range(int(dim)):
            target = source ^ (1 << bit)
            if source < target:
                edges.append(
                    {
                        "source_node": int(source),
                        "target_node": int(target),
                        "bit_index": int(bit),
                    }
                )
    return edges


def adjacent_edge_pairs(edges: list[dict[str, int]]) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    incident: dict[int, list[int]] = {}
    for idx, edge in enumerate(edges):
        incident.setdefault(edge["source_node"], []).append(idx)
        incident.setdefault(edge["target_node"], []).append(idx)
    seen: set[tuple[int, int]] = set()
    for edge_ids in incident.values():
        for left_index, left in enumerate(edge_ids):
            for right in edge_ids[left_index + 1 :]:
                pair = tuple(sorted((int(left), int(right))))
                if pair not in seen:
                    pairs.append(pair)
                    seen.add(pair)
    return pairs


def concentration_metrics(values: np.ndarray, *, prefix: str) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    arr = np.where(np.isfinite(arr), np.maximum(arr, 0.0), 0.0)
    total = float(np.sum(arr))
    n = int(arr.size)
    if n == 0 or total <= 1e-18:
        return {
            f"{prefix}_count": float(n),
            f"{prefix}_total": 0.0,
            f"{prefix}_mean": math.nan,
            f"{prefix}_sd": math.nan,
            f"{prefix}_hhi": math.nan,
            f"{prefix}_effective_count": math.nan,
            f"{prefix}_effective_fraction": math.nan,
            f"{prefix}_top1_fraction": math.nan,
            f"{prefix}_top5_fraction": math.nan,
            f"{prefix}_top10_fraction": math.nan,
            f"{prefix}_top25_fraction": math.nan,
        }
    probs = arr / total
    sorted_probs = np.sort(probs)[::-1]
    hhi = float(np.sum(probs * probs))

    def top_fraction(frac: float) -> float:
        count = max(1, int(math.ceil(frac * n)))
        return float(np.sum(sorted_probs[:count]))

    return {
        f"{prefix}_count": float(n),
        f"{prefix}_total": total,
        f"{prefix}_mean": float(np.mean(arr)),
        f"{prefix}_sd": float(np.std(arr)),
        f"{prefix}_hhi": hhi,
        f"{prefix}_effective_count": float(1.0 / hhi) if hhi > 1e-18 else math.nan,
        f"{prefix}_effective_fraction": float((1.0 / hhi) / n) if hhi > 1e-18 else math.nan,
        f"{prefix}_top1_fraction": float(sorted_probs[0]) if sorted_probs.size else math.nan,
        f"{prefix}_top5_fraction": top_fraction(0.05),
        f"{prefix}_top10_fraction": top_fraction(0.10),
        f"{prefix}_top25_fraction": top_fraction(0.25),
    }


def line_graph_autocorrelation(edge_energy: np.ndarray, edges: list[dict[str, int]]) -> dict[str, float]:
    pairs = adjacent_edge_pairs(edges)
    if not pairs:
        return {
            "edge_adjacency_pair_count": 0.0,
            "edge_energy_adjacent_pearson": math.nan,
            "edge_energy_adjacent_moran_i": math.nan,
        }
    arr = np.asarray(edge_energy, dtype=np.float64)
    left = np.asarray([arr[i] for i, _j in pairs], dtype=np.float64)
    right = np.asarray([arr[j] for _i, j in pairs], dtype=np.float64)
    mean = float(np.mean(arr))
    denom = float(np.sum((arr - mean) ** 2))
    moran = (
        float(len(arr) / len(pairs) * np.sum((left - mean) * (right - mean)) / denom)
        if denom > 1e-18
        else math.nan
    )
    return {
        "edge_adjacency_pair_count": float(len(pairs)),
        "edge_energy_adjacent_pearson": safe_corr(left, right),
        "edge_energy_adjacent_moran_i": moran,
    }


def local_ruggedness_metrics_for_values(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    edges = hypercube_edges(dim)
    edge_energy = np.asarray(
        [
            (float(arr[edge["target_node"]]) - float(arr[edge["source_node"]])) ** 2
            for edge in edges
        ],
        dtype=np.float64,
    )
    node_incident = np.zeros(1 << int(dim), dtype=np.float64)
    site_energy = np.zeros(int(dim), dtype=np.float64)
    for edge, energy in zip(edges, edge_energy):
        node_incident[edge["source_node"]] += float(energy)
        node_incident[edge["target_node"]] += float(energy)
        site_energy[edge["bit_index"]] += float(energy)
    plaquettes = plaquette_energy_by_pair(arr, dim=int(dim))
    pair_energy = np.asarray([float(value) for value in plaquettes.values()], dtype=np.float64)
    return {
        **line_graph_autocorrelation(edge_energy, edges),
        **concentration_metrics(edge_energy, prefix="edge_energy"),
        **concentration_metrics(node_incident, prefix="node_incident"),
        **concentration_metrics(site_energy, prefix="site_energy"),
        **concentration_metrics(pair_energy, prefix="plaquette_energy"),
    }


def decode_parameters(
    raw: dict[str, float | int],
    *,
    config: dict[str, Any],
    family: str,
    parameter_set_id: str,
) -> dict[str, Any]:
    base_anchors = [str(x) for x in config["abc"]["base_anchors"]]
    directions = [str(x) for x in config["abc"]["operator_gate_directions"]]
    base_index = int(raw.get("base_anchor_index", 0))
    direction_index = int(raw.get("operator_gate_direction_index", 0))
    base_index = max(0, min(base_index, len(base_anchors) - 1))
    direction_index = max(0, min(direction_index, len(directions) - 1))
    module_source = FAMILY_MODULE_SOURCES[family]
    use_modules = module_source != "none"
    uncoupled = family == "uncoupled_module"
    return {
        "candidate_model": f"abc_{family}",
        "candidate_family": family,
        "parameter_set_id": parameter_set_id,
        "base_anchor": base_anchors[base_index],
        "module_source": module_source,
        "use_modules": use_modules,
        "uncoupled_auxiliary": uncoupled,
        "auxiliary_correlation": 0.0 if uncoupled else float(raw.get("auxiliary_correlation", 0.75)),
        "operator_lambda": float(raw.get("operator_lambda", 0.25)),
        "operator_gate_midpoint": float(raw.get("operator_gate_midpoint", 0.0)),
        "operator_gate_steepness": float(raw.get("operator_gate_steepness", 1.0)),
        "operator_gate_direction": directions[direction_index],
        "module_residual_scale": 0.0 if not use_modules else float(raw.get("module_residual_scale", 0.25)),
        "module_residual_density": 0.0 if not use_modules else float(raw.get("module_residual_density", 1.0)),
        "module_residual_orders": [int(x) for x in config["coupled_family"]["module_residual_orders"]],
        **{f"raw__{key}": value for key, value in raw.items()},
    }


def build_graphlets(
    *,
    config: dict[str, Any],
    project_root: Path,
    output_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, list[dict[str, Any]]], dict[str, Any]]:
    source_output = resolve_source_output_dir(config, project_root)
    system_dir = source_output / "systems" / ASSAY_ID
    fit_dir = system_dir / "fit"

    fit_summary = json.loads((fit_dir / "summary.json").read_text(encoding="utf-8"))
    source_summary = json.loads((source_output / "summary.json").read_text(encoding="utf-8"))
    panel = pd.read_csv(fit_dir / "selected_panel.csv")
    wildtype = str(panel.loc[panel["DMS_id"] == ASSAY_ID, "target_seq"].iloc[0])
    assay_path = fit_dir / "checkpoints" / "panel_preparation" / ASSAY_ID / "prepared_assay.csv"
    metrics_path = fit_dir / "checkpoints" / "panel_preparation" / ASSAY_ID / "mavenn_metrics.json"
    assay_df = pd.read_csv(assay_path)
    mavenn_metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    score_columns = [str(x) for x in config["graphlet_evaluation"]["score_columns"]]
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    wt_scores = {
        "latent_phi": float(mavenn_metrics.get("wildtype_phi", math.nan)),
        "mavenn_yhat": float(mavenn_metrics.get("wildtype_yhat", math.nan)),
    }
    assay_score_columns = [col for col in score_columns if col in assay_df.columns or col in wt_scores]
    score_map = score_map_for_assay(
        assay_df,
        wildtype=wildtype,
        wt_scores=wt_scores,
        score_columns=assay_score_columns,
    )

    graphlet_cfg = config["phot_graphlets"]
    complete = []
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
    best_landscape = build_landscape(
        base_config=base_config,
        updates=best_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed,
    )
    mean_landscape = build_landscape(
        base_config=base_config,
        updates=mean_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed + 1,
    )
    smooth_landscape = build_landscape(
        base_config=base_config,
        updates=best_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed,
        smooth_core_only=True,
    )

    graphlets: list[dict[str, Any]] = []
    for graphlet_index, (dim, features) in enumerate(complete):
        split = graphlet_split(graphlet_index)
        masks = subset_masks(features)
        empirical = {
            score_col: np.asarray([score_map[mask][score_col] for mask in masks], dtype=np.float64)
            for score_col in assay_score_columns
        }
        best_values = evaluate_landscape_on_graphlet(best_landscape, wildtype=wildtype, features=features)
        mean_values = evaluate_landscape_on_graphlet(mean_landscape, wildtype=wildtype, features=features)
        smooth_values = evaluate_landscape_on_graphlet(smooth_landscape, wildtype=wildtype, features=features)
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
                "fitted_current_best": best_values,
                "fitted_current_posterior_mean": mean_values,
                "fitted_smooth_core_only": smooth_values,
                "fitted_residual_only_diagnostic": best_values - smooth_values,
            }
        )

    readout_base_anchor = str(config["empirical_innovation"]["base_anchor"])
    calib_base: list[float] = []
    calib_empirical: list[float] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        calib_base.extend(np.asarray(graphlet[readout_base_anchor], dtype=np.float64).tolist())
        calib_empirical.extend(np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64).tolist())
    readout = fit_readout_model(
        calib_base,
        calib_empirical,
        degree=int(config["empirical_innovation"].get("global_readout_polynomial_degree", 3)),
    )

    empirical_innovation: dict[str, np.ndarray] = {}
    for graphlet in graphlets:
        pred = apply_readout_model(np.asarray(graphlet[readout_base_anchor], dtype=np.float64), readout)
        empirical = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        empirical_innovation[str(graphlet["graphlet_id"])] = empirical - pred

    modules = select_empirical_modules(
        graphlets,
        empirical_innovation,
        max_modules=int(config["module_selection"]["max_modules"]),
        module_size=int(config["module_selection"]["module_size"]),
    )
    feature_pool = sorted({feature for graphlet in graphlets for feature in graphlet["features"]}, key=str)
    shuffled_modules = control_modules(
        modules,
        feature_pool,
        source="shuffled_module",
        rng=np.random.RandomState(graphlet_runner.stable_seed("discriminating-abc-shuffled", ASSAY_ID)),
    )
    random_modules = control_modules(
        modules,
        feature_pool,
        source="random_module",
        rng=np.random.RandomState(graphlet_runner.stable_seed("discriminating-abc-random", ASSAY_ID)),
    )
    modules_by_source = {
        "empirical_calibration_plaquette": modules,
        "shuffled_module": shuffled_modules,
        "random_module": random_modules,
        "uncoupled_empirical_module": modules,
        "none": [],
    }

    source_phot_best = math.nan
    for system_summary in source_summary.get("systems", []):
        if str(system_summary.get("assay_id")) == ASSAY_ID:
            source_phot_best = float(system_summary["fit_summary"]["best_distance"])
            break

    context = {
        "source_output_dir": str(source_output),
        "source_phot_best_objective": source_phot_best,
        "fit_summary": fit_summary,
        "readout": readout,
        "empirical_innovation": empirical_innovation,
        "assay_score_columns": assay_score_columns,
        "wildtype": wildtype,
        "alignment_path": str(alignment_path),
    }
    return graphlets, context, modules_by_source, readout


def record_metric_values(
    *,
    graphlet: dict[str, Any],
    candidate_model: str,
    candidate_family: str,
    parameter_set_id: str,
    values: np.ndarray,
    score_column: str,
    metric_rows: list[dict[str, Any]],
    walsh_rows: list[dict[str, Any]],
    gft_rows: list[dict[str, Any]],
    random_walk_cfg: dict[str, Any],
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    record_graphlet = {
        "graphlet_id": graphlet["graphlet_id"],
        "assay_id": ASSAY_ID,
        "graphlet_split": graphlet["graphlet_split"],
        "split": graphlet["split"],
        "graphlet_dimension": graphlet["graphlet_dimension"],
        "mutation_features": graphlet["mutation_features"],
        "score_column": score_column,
    }
    metric, walsh, gft, _access = metric_rows_for_values(
        graphlet=record_graphlet,
        values=np.asarray(values, dtype=np.float64),
        candidate_model=candidate_model,
        parameter_set_id=parameter_set_id,
        cfg=random_walk_cfg,
    )
    metric.update(
        local_ruggedness_metrics_for_values(
            np.asarray(values, dtype=np.float64),
            dim=int(graphlet["graphlet_dimension"]),
        )
    )
    metric["candidate_family"] = candidate_family
    for row in walsh:
        row["candidate_family"] = candidate_family
    for row in gft:
        row["candidate_family"] = candidate_family
    if extra:
        metric.update(extra)
        for row in walsh:
            row.update(extra)
        for row in gft:
            row.update(extra)
    metric_rows.append(metric)
    walsh_rows.extend(walsh)
    gft_rows.extend(gft)
    return metric


def calibration_target_rows(
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
        record_metric_values(
            graphlet=graphlet,
            candidate_model=f"empirical_observed__{primary_score}",
            candidate_family="empirical_target",
            parameter_set_id=primary_score,
            values=np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64),
            score_column=primary_score,
            metric_rows=rows,
            walsh_rows=walsh,
            gft_rows=gft,
            random_walk_cfg=random_walk_cfg,
            extra={"base_anchor": "", "module_source": ""},
        )
    return pd.DataFrame(rows)


def guardrail_penalty(by_metric: dict[str, float], config: dict[str, Any]) -> tuple[float, dict[str, Any]]:
    guardrails = dict(config["abc"].get("component_guardrails", {}))
    default_max = float(guardrails.get("default_max_component_distance", math.inf))
    penalty_weight = float(guardrails.get("penalty_weight", 1.0))
    metric_max = {str(k): float(v) for k, v in dict(guardrails.get("metric_max_component_distance", {})).items()}
    penalty = 0.0
    violations: list[str] = []
    worst_margin = 0.0
    for metric, value in by_metric.items():
        limit = metric_max.get(metric, default_max)
        if not np.isfinite(value) or not np.isfinite(limit) or limit <= 0.0:
            continue
        if float(value) > limit:
            margin = (float(value) - limit) / limit
            penalty += penalty_weight * margin
            worst_margin = max(worst_margin, margin)
            violations.append(metric)
    return float(penalty), {
        "guardrail_penalty": float(penalty),
        "n_guardrail_violations": int(len(violations)),
        "guardrail_violations": ";".join(violations),
        "worst_guardrail_fractional_margin": float(worst_margin),
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
    params = decode_parameters(raw, config=config, family=family, parameter_set_id=parameter_set_id)
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    active_rows: list[dict[str, Any]] = []
    for graphlet in graphlets:
        if graphlet["graphlet_split"] != "calibration":
            continue
        base_values = np.asarray(graphlet[params["base_anchor"]], dtype=np.float64)
        values, innovation, active, coeff_energy = coupled_values(
            base_values=base_values,
            graphlet_features=tuple(graphlet["features"]),
            modules=modules_by_source[str(params["module_source"])],
            parameter_set=params,
            seed=graphlet_runner.stable_seed(parameter_set_id, graphlet["graphlet_id"]),
        )
        extra = {
            key: value
            for key, value in params.items()
            if key not in {"module_residual_orders"}
        }
        extra["module_residual_orders"] = ";".join(str(x) for x in params["module_residual_orders"])
        extra.update(active)
        extra["module_coefficient_energy_total"] = float(sum(coeff_energy.values()))
        extra["innovation_sd_observed"] = float(np.std(innovation))
        record_metric_values(
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
        active_rows.append(extra)

    model_rows = pd.DataFrame(metric_rows)
    distance, by_metric = metric_distance(model_rows, target_rows, metric_cols)
    penalty, guardrail = guardrail_penalty(by_metric, config)
    if not np.isfinite(distance):
        distance = float("inf")
    objective = float(distance + penalty)
    summary = [float(distance), penalty] + [float(by_metric.get(metric, math.nan)) for metric in metric_cols]
    extras = {
        **params,
        "calibration_graphlet_distance": float(distance),
        "objective_distance": objective,
        **guardrail,
        **{f"distance_component__{key}": value for key, value in by_metric.items()},
    }
    if active_rows:
        active_df = pd.DataFrame(active_rows)
        for col in [
            "innovation_to_base_sd_ratio",
            "auxiliary_phenotype_correlation_to_base",
            "innovation_phenotype_correlation_to_base",
            "mean_gate",
            "module_delta_sd",
            "module_coefficient_energy_total",
        ]:
            if col in active_df:
                values = pd.to_numeric(active_df[col], errors="coerce").dropna()
                if not values.empty:
                    extras[f"{col}_mean"] = float(values.mean())
    return SimulationResult(distance=objective, summary_vector=np.asarray(summary, dtype=np.float64), extras=extras)


def progress_writer(output_dir: Path) -> Any:
    progress_path = output_dir / "progress.json"

    def _write(event: dict[str, Any]) -> None:
        write_json(progress_path, event)

    return _write


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
    seed = int(smc_cfg.get("seed", 0)) + int(graphlet_runner.stable_seed("family", family) % 100000)
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
        path=str(output_dir / "checkpoints" / f"smc_{family}.json"),
        resume=bool(config["abc"].get("resume_checkpoints", True)),
    )

    def simulate(raw: dict[str, float | int]) -> SimulationResult:
        raw_payload = json.dumps(raw, sort_keys=True)
        parameter_set_id = f"{family}__abc_{graphlet_runner.stable_seed(family, raw_payload):010d}"
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
        run_label=f"abc_{family}",
    )
    particle_rows: list[dict[str, Any]] = []
    for rank, particle in enumerate(sorted(result.particles, key=lambda item: item.distance), start=1):
        row = {
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
        particle_rows.append(row)
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
    for family, sub in posterior_df.groupby("candidate_family", dropna=False):
        sub = sub.sort_values("objective_distance").head(n_predictive)
        for rank, row in enumerate(sub.to_dict(orient="records"), start=1):
            raw = {
                key.replace("raw__", "", 1): row[key]
                for key in row
                if str(key).startswith("raw__") and not str(key).startswith("raw__raw__")
            }
            parameter_set_id = f"{family}__posterior_rank{rank:03d}"
            decoded = decode_parameters(
                raw,
                config=config,
                family=str(family),
                parameter_set_id=parameter_set_id,
            )
            decoded["posterior_rank"] = int(rank)
            decoded["posterior_objective_distance"] = float(row["objective_distance"])
            candidate_particle_rows.append(decoded)

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
            extra={"base_anchor": "", "module_source": ""},
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
            extra={"base_anchor": readout_base_anchor, "module_source": "empirical_calibration_plaquette"},
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
                extra={"base_anchor": candidate_model, "module_source": "baseline_control"},
            )

        gnk = gnk_values(
            dim=int(graphlet["graphlet_dimension"]),
            seed=graphlet_runner.stable_seed("phot-discriminating-abc-gnk", graphlet["graphlet_id"]),
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
            extra={"base_anchor": "", "module_source": "gnk_control"},
        )

        for params in candidate_particle_rows:
            base_values = np.asarray(graphlet[params["base_anchor"]], dtype=np.float64)
            values, innovation, active, coeff_energy = coupled_values(
                base_values=base_values,
                graphlet_features=tuple(graphlet["features"]),
                modules=modules_by_source[str(params["module_source"])],
                parameter_set=params,
                seed=graphlet_runner.stable_seed(params["parameter_set_id"], graphlet["graphlet_id"]),
            )
            extra = {key: value for key, value in params.items() if key != "module_residual_orders"}
            extra["module_residual_orders"] = ";".join(str(x) for x in params["module_residual_orders"])
            extra.update(active)
            extra["module_coefficient_energy_total"] = float(sum(coeff_energy.values()))
            record_metric_values(
                graphlet=graphlet,
                candidate_model=str(params["candidate_model"]),
                candidate_family=str(params["candidate_family"]),
                parameter_set_id=str(params["parameter_set_id"]),
                values=values,
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
                innovation=innovation,
                modules=modules_by_source[str(params["module_source"])],
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
                base_values=base_values,
                values=values,
                innovation=innovation,
            )
            corr.update(extra)
            corr["module_coefficient_energy_total"] = float(sum(coeff_energy.values()))
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
                    "base_anchor": first.get("base_anchor", ""),
                    "module_source": first.get("module_source", ""),
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
            "innovation_to_base_sd_ratio",
            "n_active_modules",
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
    static_families = ["fitted_current", "smooth_core", "fitted_current_posterior_mean", "residual_only_diagnostic", "gnk_local_null", "empirical_readout_residual"]
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
        references = ["fitted_current", "smooth_core", "no_module", "empirical_module", "random_module", "shuffled_module", "uncoupled_module"]
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
        by_metric = {
            metric: float(row.get(f"distance_component__{metric}", math.nan))
            for metric in metric_cols
        }
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


def top_distance_rows(distances: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not distances:
        return []
    df = pd.DataFrame(distances)
    df["distance_to_empirical_primary"] = pd.to_numeric(df["distance_to_empirical_primary"], errors="coerce")
    rows: list[dict[str, Any]] = []
    for split in ["calibration", "test", "all"]:
        sub = df[df["graphlet_split"] == split].copy()
        sub = sub[np.isfinite(sub["distance_to_empirical_primary"])]
        for rank, row in enumerate(sub.sort_values("distance_to_empirical_primary").head(30).to_dict(orient="records"), start=1):
            rows.append({"rank": rank, **row})
    return rows


def graphlet_catalog_rows(graphlets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "graphlet_id": graphlet["graphlet_id"],
            "assay_id": ASSAY_ID,
            "graphlet_split": graphlet["graphlet_split"],
            "graphlet_dimension": graphlet["graphlet_dimension"],
            "mutation_features": graphlet["mutation_features"],
        }
        for graphlet in graphlets
    ]


def module_catalog_rows(modules_by_source: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source, modules in modules_by_source.items():
        for module in modules:
            rows.append(
                {
                    "module_id": module["module_id"],
                    "module_source": source,
                    "rank": module.get("rank", math.nan),
                    "calibration_pair_score": module.get("calibration_pair_score", math.nan),
                    "module_features": ";".join(feature_label(feature) for feature in module["features"]),
                }
            )
    return rows


def write_preflight(
    *,
    config: dict[str, Any],
    output_dir: Path,
    graphlets: list[dict[str, Any]],
    context: dict[str, Any],
    modules_by_source: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    preflight = {
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_result": config["source_run"].get("result_id", "RES-030"),
        "source_output_dir": context["source_output_dir"],
        "mode": str(config.get("experiment", {}).get("mode", "discriminating_abc_over_run079_coupled_family_controls")),
        "primary_score_column": config["graphlet_evaluation"]["primary_score_column"],
        "score_columns": context["assay_score_columns"],
        "objective_metric_columns": [str(x) for x in config["graphlet_evaluation"]["metric_columns"]],
        "local_ruggedness_metric_columns": [
            metric
            for metric in LOCAL_RUGGEDNESS_METRIC_COLUMNS
            if metric in [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]
        ],
        "fit_families": [str(x) for x in config["abc"]["fit_families"]],
        "base_anchors": [str(x) for x in config["abc"]["base_anchors"]],
        "n_graphlets": int(len(graphlets)),
        "n_calibration_graphlets": int(sum(1 for graphlet in graphlets if graphlet["graphlet_split"] == "calibration")),
        "n_test_graphlets": int(sum(1 for graphlet in graphlets if graphlet["graphlet_split"] == "test")),
        "graphlet_dimensions": sorted({int(graphlet["graphlet_dimension"]) for graphlet in graphlets}),
        "module_counts": {source: int(len(modules)) for source, modules in modules_by_source.items()},
        "readout_base_anchor": config["empirical_innovation"]["base_anchor"],
        "run077_phot_best_objective": context["source_phot_best_objective"],
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
                "module_energy_summary": str(output_dir / "module_energy_summary.csv"),
                "base_correlation": str(output_dir / "base_correlation.csv"),
                "empirical_readout_fit": str(output_dir / "empirical_readout_fit.json"),
            },
            **output_summary,
        },
    )


if __name__ == "__main__":
    main()
