#!/usr/bin/env python3
"""Deterministic PHOT coupled-landscape grid search.

This is the first HYP-008 machinery test. It does not run ABC. It freezes the
current PHOT model from RUN-077/RES-030, constructs coupled auxiliary landscape
candidates over the same dense PHOT graphlets, and records both graphlet fit and
innovation-locality diagnostics.
"""

from __future__ import annotations

import argparse
import itertools
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
from project_adapt_env.coupled_landscape import (  # noqa: E402
    aggregate_module_energy,
    apply_polynomial_readout,
    control_modules,
    coupled_values,
    edge_energy_by_bit,
    energy_locality_for_modules,
    fit_polynomial_readout,
    node_edge_gradients,
    plaquette_energy_by_pair,
    safe_corr,
    select_empirical_modules,
)
from project_adapt_env.proteingym_panel import merge_dict  # noqa: E402
from proteingym_system_specific_current_graphlet_refit import (  # noqa: E402
    build_landscape,
    evaluate_landscape_on_graphlet,
    gnk_values,
    metric_distance,
)


ASSAY_ID = "PHOT_CHLRE_Chen_2023"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=graphlet_runner.json_default)
        + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if quick:
        merged = merge_dict(raw, raw.get("quick", {}))
        merged.pop("quick", None)
        return merged
    raw.pop("quick", None)
    return raw


def resolve_source_output_dir(config: dict[str, Any], project_root: Path) -> Path:
    source = Path(str(config["source_run"]["output_dir"]))
    if not source.is_absolute():
        source = project_root / source
    if not source.exists():
        raise FileNotFoundError(
            f"RUN-077 source output directory is missing: {source}. "
            "Restore it with DVC before running this deterministic scan."
        )
    return source


def feature_label(feature: tuple[int, str]) -> str:
    return f"{int(feature[0])}:{feature[1]}"


def graphlet_split(index: int) -> str:
    return "calibration" if int(index) % 2 == 0 else "test"


def fit_readout_model(base_values: list[float], empirical_values: list[float], *, degree: int) -> dict[str, Any]:
    x = np.asarray(base_values, dtype=np.float64)
    y = np.asarray(empirical_values, dtype=np.float64)
    x_mean = float(np.mean(x))
    x_sd = max(float(np.std(x)), 1e-12)
    x0 = (x - x_mean) / x_sd
    columns = [np.ones_like(x0)]
    for power in range(1, int(degree) + 1):
        columns.append(x0**power)
    design = np.vstack(columns).T
    beta, *_rest = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ beta
    sst = float(np.sum((y - float(np.mean(y))) ** 2))
    ssr = float(np.sum((y - pred) ** 2))
    return {
        "degree": int(degree),
        "x_mean": x_mean,
        "x_sd": x_sd,
        "beta": [float(v) for v in beta],
        "calibration_r2": float(1.0 - ssr / sst) if sst > 1e-12 else 1.0,
    }


def apply_readout_model(base_values: np.ndarray, model: dict[str, Any]) -> np.ndarray:
    x = np.asarray(base_values, dtype=np.float64)
    x0 = (x - float(model["x_mean"])) / max(float(model["x_sd"]), 1e-12)
    beta = np.asarray(model["beta"], dtype=np.float64)
    columns = [np.ones_like(x0)]
    for power in range(1, beta.size):
        columns.append(x0**power)
    return np.vstack(columns).T @ beta


def parameter_grid(config: dict[str, Any]) -> list[dict[str, Any]]:
    grid = config["coupled_grid"]
    base_anchors = [str(x) for x in grid["base_anchors"]]
    module_sources = [str(x) for x in grid["module_sources"]]
    rows: list[dict[str, Any]] = []
    idx = 0
    for base_anchor, module_source in itertools.product(base_anchors, module_sources):
        if module_source == "none":
            local_scales = [0.0]
            densities = [0.0]
        else:
            local_scales = [float(x) for x in grid["module_residual_scales"]]
            densities = [float(x) for x in grid["module_residual_densities"]]
        for rho, lam, scale, density, midpoint, steepness, direction in itertools.product(
            [float(x) for x in grid["auxiliary_correlations"]],
            [float(x) for x in grid["operator_lambdas"]],
            local_scales,
            densities,
            [float(x) for x in grid["operator_gate_midpoints"]],
            [float(x) for x in grid["operator_gate_steepnesses"]],
            [str(x) for x in grid["operator_gate_directions"]],
        ):
            idx += 1
            candidate = (
                "coupled_no_module"
                if module_source == "none"
                else "uncoupled_local_module"
                if module_source == "uncoupled_empirical_module"
                else f"coupled_{module_source}"
            )
            rows.append(
                {
                    "candidate_model": candidate,
                    "parameter_set_id": f"{candidate}__{base_anchor}__grid{idx:04d}",
                    "base_anchor": base_anchor,
                    "module_source": module_source,
                    "use_modules": module_source != "none",
                    "uncoupled_auxiliary": module_source == "uncoupled_empirical_module",
                    "auxiliary_correlation": rho,
                    "operator_lambda": lam,
                    "operator_gate_midpoint": midpoint,
                    "operator_gate_steepness": steepness,
                    "operator_gate_direction": direction,
                    "module_residual_scale": scale,
                    "module_residual_density": density,
                    "module_residual_orders": [int(x) for x in grid["module_residual_orders"]],
                }
            )
    return rows


def summarize_metric_rows(
    rows: list[dict[str, Any]],
    *,
    group_cols: list[str],
    metric_cols: list[str],
) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {col: value for col, value in zip(group_cols, keys)}
        row["n_graphlets"] = int(len(sub))
        for col in metric_cols:
            if col not in sub:
                continue
            values = pd.to_numeric(sub[col], errors="coerce").dropna()
            if values.empty:
                continue
            row[f"{col}_mean"] = float(values.mean())
            row[f"{col}_sd"] = float(values.std(ddof=0))
        out.append(row)
    return out


def correlation_row(
    *,
    graphlet: dict[str, Any],
    candidate_model: str,
    parameter_set_id: str,
    base_values: np.ndarray,
    values: np.ndarray,
    innovation: np.ndarray,
) -> dict[str, Any]:
    dim = int(graphlet["graphlet_dimension"])
    return {
        "graphlet_id": graphlet["graphlet_id"],
        "graphlet_split": graphlet["graphlet_split"],
        "candidate_model": candidate_model,
        "parameter_set_id": parameter_set_id,
        "base_phenotype_correlation": safe_corr(values, base_values),
        "innovation_phenotype_correlation": safe_corr(innovation, base_values),
        "edge_gradient_correlation": safe_corr(
            node_edge_gradients(values, dim=dim),
            node_edge_gradients(base_values, dim=dim),
        ),
        "innovation_edge_gradient_correlation": safe_corr(
            node_edge_gradients(innovation, dim=dim),
            node_edge_gradients(base_values, dim=dim),
        ),
        "edge_innovation_energy": float(sum(edge_energy_by_bit(innovation, dim=dim).values())),
        "plaquette_innovation_energy": float(sum(plaquette_energy_by_pair(innovation, dim=dim).values())),
        "innovation_sd": float(np.std(innovation)),
        "base_sd": float(np.std(base_values)),
        "innovation_to_base_sd_ratio": float(np.std(innovation) / max(float(np.std(base_values)), 1e-12)),
    }


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
        rng=np.random.RandomState(graphlet_runner.stable_seed("shuffled-modules", ASSAY_ID)),
    )
    random_modules = control_modules(
        modules,
        feature_pool,
        source="random_module",
        rng=np.random.RandomState(graphlet_runner.stable_seed("random-modules", ASSAY_ID)),
    )
    modules_by_source = {
        "empirical_calibration_plaquette": modules,
        "shuffled_module": shuffled_modules,
        "random_module": random_modules,
        "uncoupled_empirical_module": modules,
        "none": [],
    }

    random_walk_cfg = {"analysis": config.get("analysis", {})}
    metric_cols = [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    locality_rows: list[dict[str, Any]] = []
    locality_summary_rows: list[dict[str, Any]] = []
    correlation_rows: list[dict[str, Any]] = []
    parameter_rows = parameter_grid(config)

    def record_values(
        *,
        graphlet: dict[str, Any],
        candidate_model: str,
        parameter_set_id: str,
        values: np.ndarray,
        score_column: str,
        extra: dict[str, Any] | None = None,
    ) -> None:
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
        if extra:
            metric.update(extra)
            for row in walsh:
                row.update(extra)
            for row in gft:
                row.update(extra)
        metric_rows.append(metric)
        walsh_rows.extend(walsh)
        gft_rows.extend(gft)

    for graphlet in graphlets:
        base_for_empirical = np.asarray(graphlet[readout_base_anchor], dtype=np.float64)
        empirical_resid = empirical_innovation[str(graphlet["graphlet_id"])]
        emp_loc_rows, emp_loc_summary = energy_locality_for_modules(
            graphlet_id=str(graphlet["graphlet_id"]),
            graphlet_split=str(graphlet["graphlet_split"]),
            graphlet_features=tuple(graphlet["features"]),
            innovation=empirical_resid,
            modules=modules,
            candidate_model="empirical_global_readout_residual",
            parameter_set_id=f"{primary_score}__minus_poly_{readout_base_anchor}",
        )
        locality_rows.extend(emp_loc_rows)
        locality_summary_rows.append(emp_loc_summary)
        correlation_rows.append(
            correlation_row(
                graphlet=graphlet,
                candidate_model="empirical_global_readout_residual",
                parameter_set_id=f"{primary_score}__minus_poly_{readout_base_anchor}",
                base_values=base_for_empirical,
                values=np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64),
                innovation=empirical_resid,
            )
        )
        record_values(
            graphlet=graphlet,
            candidate_model="empirical_global_readout_residual",
            parameter_set_id=f"{primary_score}__minus_poly_{readout_base_anchor}",
            values=empirical_resid,
            score_column=primary_score,
            extra={"base_anchor": readout_base_anchor, "module_source": "empirical_calibration_plaquette"},
        )

        for score_col in assay_score_columns:
            record_values(
                graphlet=graphlet,
                candidate_model=f"empirical_observed__{score_col}",
                parameter_set_id=score_col,
                values=np.asarray(graphlet["empirical_values"][score_col], dtype=np.float64),
                score_column=score_col,
                extra={"base_anchor": "", "module_source": ""},
            )

        for candidate_model in [
            "fitted_current_best",
            "fitted_current_posterior_mean",
            "fitted_smooth_core_only",
            "fitted_residual_only_diagnostic",
        ]:
            values = np.asarray(graphlet[candidate_model], dtype=np.float64)
            record_values(
                graphlet=graphlet,
                candidate_model=candidate_model,
                parameter_set_id=f"{ASSAY_ID}__{candidate_model}",
                values=values,
                score_column=primary_score,
                extra={"base_anchor": candidate_model, "module_source": "baseline_control"},
            )

        gnk = gnk_values(
            dim=int(graphlet["graphlet_dimension"]),
            seed=graphlet_runner.stable_seed("phot-coupled-gnk", graphlet["graphlet_id"]),
            target_mean=float(np.mean(graphlet["empirical_values"][primary_score])),
            target_sd=float(np.std(graphlet["empirical_values"][primary_score])),
            k=2,
        )
        record_values(
            graphlet=graphlet,
            candidate_model="gnk_local_null",
            parameter_set_id=f"{ASSAY_ID}__gnk_k2",
            values=gnk,
            score_column=primary_score,
            extra={"base_anchor": "", "module_source": "gnk_control"},
        )

        for params in parameter_rows:
            base_values = np.asarray(graphlet[params["base_anchor"]], dtype=np.float64)
            values, innovation, active, coeff_energy = coupled_values(
                base_values=base_values,
                graphlet_features=tuple(graphlet["features"]),
                modules=modules_by_source[str(params["module_source"])],
                parameter_set=params,
                seed=graphlet_runner.stable_seed(params["parameter_set_id"], graphlet["graphlet_id"]),
            )
            extra = {
                key: value
                for key, value in params.items()
                if key != "module_residual_orders"
            }
            extra["module_residual_orders"] = ";".join(str(x) for x in params["module_residual_orders"])
            extra.update(active)
            record_values(
                graphlet=graphlet,
                candidate_model=str(params["candidate_model"]),
                parameter_set_id=str(params["parameter_set_id"]),
                values=values,
                score_column=primary_score,
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
            ["candidate_model", "parameter_set_id"],
            dropna=False,
        ):
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
                    "base_anchor": first.get("base_anchor", ""),
                    "module_source": first.get("module_source", ""),
                    "distance_to_empirical_primary": distance,
                    **{f"distance_component__{key}": value for key, value in by_metric.items()},
                }
            )

    summary_rows = summarize_metric_rows(
        metric_rows,
        group_cols=["graphlet_split", "candidate_model", "parameter_set_id", "score_column"],
        metric_cols=metric_cols
        + [
            "innovation_to_base_sd_ratio",
            "n_active_modules",
            "plaquette_innovation_total_energy",
            "edge_innovation_total_energy",
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
    module_rows = [
        {
            "module_id": module["module_id"],
            "module_source": module["module_source"],
            "rank": module.get("rank", math.nan),
            "calibration_pair_score": module.get("calibration_pair_score", math.nan),
            "module_features": ";".join(feature_label(feature) for feature in module["features"]),
        }
        for module in modules
    ]
    module_summary = aggregate_module_energy(locality_rows)

    write_csv(output_dir / "graphlet_catalog.csv", catalog_rows)
    write_csv(output_dir / "module_catalog.csv", module_rows)
    write_csv(output_dir / "parameter_grid.csv", parameter_rows)
    write_csv(output_dir / "graphlet_metrics.csv", metric_rows)
    write_csv(output_dir / "walsh_terms.csv", walsh_rows)
    write_csv(output_dir / "gft_coefficients.csv", gft_rows)
    write_csv(output_dir / "model_distances.csv", distances)
    write_csv(output_dir / "graphlet_summary.csv", summary_rows)
    write_csv(output_dir / "module_energy_rows.csv", locality_rows)
    write_csv(output_dir / "module_energy_by_graphlet.csv", locality_summary_rows)
    write_csv(output_dir / "module_energy_summary.csv", module_summary)
    write_csv(output_dir / "base_correlation.csv", correlation_rows)
    write_json(output_dir / "empirical_readout_fit.json", readout)

    distance_df = pd.DataFrame(distances)
    top_rows: list[dict[str, Any]] = []
    if not distance_df.empty:
        for split in ["calibration", "test", "all"]:
            sub = distance_df[distance_df["graphlet_split"] == split].copy()
            sub["distance_to_empirical_primary"] = pd.to_numeric(sub["distance_to_empirical_primary"], errors="coerce")
            sub = sub[np.isfinite(sub["distance_to_empirical_primary"])]
            if sub.empty:
                continue
            for rank, row in enumerate(sub.sort_values("distance_to_empirical_primary").head(20).to_dict(orient="records"), start=1):
                top_rows.append({"rank": rank, **row})
    write_csv(output_dir / "top_model_distances.csv", top_rows)

    source_phot_best = math.nan
    for system_summary in source_summary.get("systems", []):
        if str(system_summary.get("assay_id")) == ASSAY_ID:
            source_phot_best = float(system_summary["fit_summary"]["best_distance"])
            break
    preflight = {
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_output_dir": str(source_output),
        "source_result": "RES-030",
        "mode": "deterministic_frozen_base_grid_no_abc",
        "primary_score_column": primary_score,
        "score_columns": assay_score_columns,
        "base_anchors": sorted({str(row["base_anchor"]) for row in parameter_rows}),
        "n_parameter_sets": int(len(parameter_rows)),
        "n_modules": int(len(modules)),
        "n_graphlets": int(len(graphlets)),
        "graphlet_dimensions": sorted({int(row["graphlet_dimension"]) for row in catalog_rows}),
        "readout_base_anchor": readout_base_anchor,
        "run077_phot_best_objective": source_phot_best,
    }
    write_json(output_dir / "preflight.json", preflight)
    write_json(
        output_dir / "summary.json",
        {
            "experiment_id": config["experiment"]["id"],
            "preflight": preflight,
            "empirical_readout_fit": readout,
            "output_files": {
                "graphlet_catalog": str(output_dir / "graphlet_catalog.csv"),
                "module_catalog": str(output_dir / "module_catalog.csv"),
                "parameter_grid": str(output_dir / "parameter_grid.csv"),
                "graphlet_metrics": str(output_dir / "graphlet_metrics.csv"),
                "walsh_terms": str(output_dir / "walsh_terms.csv"),
                "gft_coefficients": str(output_dir / "gft_coefficients.csv"),
                "model_distances": str(output_dir / "model_distances.csv"),
                "module_energy_rows": str(output_dir / "module_energy_rows.csv"),
                "module_energy_summary": str(output_dir / "module_energy_summary.csv"),
                "base_correlation": str(output_dir / "base_correlation.csv"),
                "top_model_distances": str(output_dir / "top_model_distances.csv"),
            },
            "top_model_distances": top_rows[:20],
        },
    )


if __name__ == "__main__":
    main()
