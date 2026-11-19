#!/usr/bin/env python3
"""Exact PHOT descriptive failure analysis before any new ABC fit."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
for path in (SRC_ROOT, ADAPT_ENV_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from adaptenv import AlignmentProfile, FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from adaptenv.amino_acids import AA_TO_IDX  # noqa: E402
from adaptenv.calibration import CalibrationOptions, _synthetic_readout_from_components  # noqa: E402

from project_adapt_env.dense_descriptive import (  # noqa: E402
    CompleteGraphlet,
    build_score_map,
    component_failure_table,
    enumerate_complete_graphlets,
    exact_summary,
    model_distance_table,
    node_masks,
    values_for_graphlet,
)
from project_adapt_env.utils import atomic_write_dataframe_csv, atomic_write_json, to_builtin  # noqa: E402


ASSAY_ID = "PHOT_CHLRE_Chen_2023"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def merge_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge_dict(dict(out[key]), value)
        else:
            out[key] = value
    return out


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if quick:
        raw = merge_dict(raw, raw.get("quick", {}))
    raw.pop("quick", None)
    return raw


def resolve_path(project_root: Path, raw: str | os.PathLike[str]) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else project_root / path


def emit_progress(output_dir: Path, payload: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    enriched = {
        "run_id": os.environ.get("LABPROJ_RUN_ID"),
        "experiment_id": os.environ.get("LABPROJ_EXPERIMENT_ID"),
        **payload,
    }
    atomic_write_json(output_dir / "progress.json", enriched)
    with (output_dir / "progress.log").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(to_builtin(enriched), sort_keys=True))
        handle.write("\n")
    print(json.dumps(to_builtin(enriched), sort_keys=True), flush=True)


def sequence_indices_for_graphlet(wildtype: str, graphlet: CompleteGraphlet) -> np.ndarray:
    rows: list[list[int]] = []
    for mask in node_masks(graphlet.features):
        residues = list(wildtype)
        for pos, aa in mask:
            residues[int(pos)] = str(aa)
        rows.append([AA_TO_IDX[str(aa)] for aa in residues])
    return np.asarray(rows, dtype=np.int64)


def build_landscape(
    *,
    base_config: dict[str, Any],
    updates: dict[str, Any],
    wildtype: str,
    alignment_path: Path,
    seed_offset: int,
) -> FitnessLandscapeEnv:
    payload = dict(base_config)
    payload.update(dict(updates))
    payload.update(
        {
            "L": len(wildtype),
            "seed": int(payload.get("seed", 42)) + int(seed_offset),
        }
    )
    profile = AlignmentProfile.from_alignment(alignment_path, max_sequences=512)
    return FitnessLandscapeEnv(
        LandscapeConfig(**payload),
        alignment_profile=profile,
        reference_sequence=wildtype,
    )


def evaluate_landscape(
    landscape: FitnessLandscapeEnv,
    seqs: np.ndarray,
    *,
    readout_mode: str,
) -> np.ndarray:
    mode = str(readout_mode or "fitness").strip()
    if mode == "fitness":
        return np.asarray(landscape.evaluate_batch(seqs), dtype=np.float64)
    components = landscape.evaluate_batch_components(seqs)
    return np.asarray(
        _synthetic_readout_from_components(
            components,
            options=CalibrationOptions(synthetic_readout_mode=mode),
        ),
        dtype=np.float64,
    )


def graphlet_catalog_rows(graphlets: list[CompleteGraphlet], *, assay_id: str) -> list[dict[str, Any]]:
    return [
        {
            "graphlet_id": graphlet.graphlet_id,
            "assay_id": assay_id,
            "graphlet_dimension": int(graphlet.dimension),
            "mutation_features": graphlet.mutation_features,
        }
        for graphlet in graphlets
    ]


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config_path, quick=bool(args.quick))

    phot_cfg = dict(cfg["phot"])
    assay_id = str(phot_cfg.get("assay_id", ASSAY_ID))
    score_columns = [str(col) for col in cfg["score_columns"]]
    primary_score = str(cfg["primary_score_column"])
    metric_columns = [str(col) for col in cfg["metric_columns"]]

    selected_panel_path = resolve_path(project_root, phot_cfg["selected_panel_path"])
    prepared_assay_path = resolve_path(project_root, phot_cfg["prepared_assay_path"])
    mavenn_metrics_path = resolve_path(project_root, phot_cfg["mavenn_metrics_path"])
    alignment_path = resolve_path(project_root, phot_cfg["alignment_path"])
    source_config_path = resolve_path(project_root, cfg["source_full_abc"]["config_path"])
    source_summary_path = resolve_path(project_root, cfg["source_full_abc"]["summary_path"])

    emit_progress(output_dir, {"event": "loading_inputs", "prepared_assay_path": str(prepared_assay_path)})
    selected_panel = pd.read_csv(selected_panel_path)
    wildtype = str(selected_panel.loc[selected_panel["DMS_id"] == assay_id, "target_seq"].iloc[0])
    assay_df = pd.read_csv(prepared_assay_path)
    mavenn_metrics = json.loads(mavenn_metrics_path.read_text(encoding="utf-8"))
    source_config = yaml.safe_load(source_config_path.read_text(encoding="utf-8"))
    source_summary = json.loads(source_summary_path.read_text(encoding="utf-8"))

    wt_scores = {
        "latent_phi": float(mavenn_metrics.get("wildtype_phi", math.nan)),
        "mavenn_yhat": float(mavenn_metrics.get("wildtype_yhat", math.nan)),
    }
    available_score_columns = [col for col in score_columns if col in assay_df.columns or col in wt_scores]
    if primary_score not in available_score_columns:
        raise ValueError(f"primary_score_column {primary_score!r} is not available in prepared PHOT assay")
    score_map = build_score_map(
        assay_df,
        wildtype=wildtype,
        score_columns=available_score_columns,
        wt_scores=wt_scores,
    )

    graphlet_cfg = dict(cfg["graphlets"])
    graphlets = enumerate_complete_graphlets(
        score_map,
        assay_id=assay_id,
        score_columns=available_score_columns,
        dimensions=[int(x) for x in graphlet_cfg["dimensions"]],
        max_graphlets_per_dimension=int(graphlet_cfg["max_graphlets_per_dimension"]),
        top_feature_pool=int(graphlet_cfg["top_feature_pool"]),
        max_combinations_per_dimension=(
            int(graphlet_cfg["max_combinations_per_dimension"])
            if graphlet_cfg.get("max_combinations_per_dimension") is not None
            else None
        ),
    )
    if not graphlets:
        raise RuntimeError("no complete PHOT graphlets were found for exact descriptive analysis")
    emit_progress(
        output_dir,
        {
            "event": "graphlets_enumerated",
            "n_graphlets": len(graphlets),
            "dimensions": sorted({graphlet.dimension for graphlet in graphlets}),
        },
    )

    base_config = dict(source_config["adapt_env"]["base_config"])
    candidate_cfgs = list(cfg["model_candidates"])
    landscapes: dict[str, tuple[FitnessLandscapeEnv, str, dict[str, Any]]] = {}
    for index, candidate in enumerate(candidate_cfgs):
        branch = str(candidate["source_branch"])
        branch_payload = dict(source_summary["branches"][branch])
        landscapes[str(candidate["candidate_model"])] = (
            build_landscape(
                base_config=base_config,
                updates=dict(branch_payload["fitted_parameters"]),
                wildtype=wildtype,
                alignment_path=alignment_path,
                seed_offset=int(candidate.get("seed_offset", index)),
            ),
            str(candidate.get("readout_mode", cfg.get("model_readout_mode", "stability_margin"))),
            {
                "source_branch": branch,
                "source_fit_mode": branch_payload.get("fit_mode"),
                "source_objective_distance": branch_payload.get("objective_terms", {}).get("distance"),
            },
        )

    catalog_rows = graphlet_catalog_rows(graphlets, assay_id=assay_id)
    metric_rows: list[dict[str, Any]] = []
    value_rows: list[dict[str, Any]] = []
    for graphlet_index, graphlet in enumerate(graphlets):
        for score_col in available_score_columns:
            empirical_values = values_for_graphlet(score_map, graphlet, score_column=score_col)
            summary = exact_summary(
                empirical_values,
                dim=int(graphlet.dimension),
                accessibility_max_dim=int(cfg.get("accessibility_max_dim", 8)),
            )
            metric_rows.append(
                {
                    "assay_id": assay_id,
                    "graphlet_id": graphlet.graphlet_id,
                    "graphlet_index": int(graphlet_index),
                    "graphlet_dimension": int(graphlet.dimension),
                    "mutation_features": graphlet.mutation_features,
                    "score_column": score_col,
                    "candidate_model": f"empirical_observed__{score_col}",
                    "parameter_set_id": score_col,
                    **summary,
                }
            )
            for node, value in enumerate(empirical_values):
                value_rows.append(
                    {
                        "graphlet_id": graphlet.graphlet_id,
                        "node": int(node),
                        "score_column": score_col,
                        "candidate_model": f"empirical_observed__{score_col}",
                        "value": float(value),
                    }
                )

        seqs = sequence_indices_for_graphlet(wildtype, graphlet)
        for candidate_model, (landscape, readout_mode, candidate_meta) in landscapes.items():
            values = evaluate_landscape(landscape, seqs, readout_mode=readout_mode)
            summary = exact_summary(
                values,
                dim=int(graphlet.dimension),
                accessibility_max_dim=int(cfg.get("accessibility_max_dim", 8)),
            )
            metric_rows.append(
                {
                    "assay_id": assay_id,
                    "graphlet_id": graphlet.graphlet_id,
                    "graphlet_index": int(graphlet_index),
                    "graphlet_dimension": int(graphlet.dimension),
                    "mutation_features": graphlet.mutation_features,
                    "score_column": primary_score,
                    "candidate_model": candidate_model,
                    "parameter_set_id": candidate_model,
                    "readout_mode": readout_mode,
                    **candidate_meta,
                    **summary,
                }
            )
            for node, value in enumerate(values):
                value_rows.append(
                    {
                        "graphlet_id": graphlet.graphlet_id,
                        "node": int(node),
                        "score_column": primary_score,
                        "candidate_model": candidate_model,
                        "value": float(value),
                    }
                )

    metric_df = pd.DataFrame(metric_rows)
    failure_df = component_failure_table(
        metric_df[metric_df["score_column"] == primary_score].copy(),
        target_candidate_model=f"empirical_observed__{primary_score}",
        metric_columns=metric_columns,
    )
    distance_df = model_distance_table(failure_df)
    severe_threshold = float(cfg["failure_gate"].get("severe_abs_scaled_error", 2.0))
    recommended_rows: list[dict[str, Any]] = []
    for candidate_model, sub in failure_df.groupby("candidate_model", dropna=False):
        severe = sub[sub["abs_scaled_error"] >= severe_threshold].sort_values("abs_scaled_error", ascending=False)
        recommended_rows.append(
            {
                "candidate_model": candidate_model,
                "severe_abs_scaled_error_threshold": severe_threshold,
                "n_failed_summary_components": int(len(severe)),
                "abc_allowed_by_descriptive_gate": bool(len(severe) > 0),
                "recommended_abc_summary_columns": ";".join(str(x) for x in severe["metric"].tolist()),
                "top_failure_signed_scaled_errors": ";".join(
                    f"{row.metric}:{float(row.signed_scaled_error):.3g}"
                    for row in severe.itertuples(index=False)
                ),
            }
        )
    recommended_df = pd.DataFrame(recommended_rows)

    atomic_write_dataframe_csv(output_dir / "graphlet_catalog.csv", pd.DataFrame(catalog_rows))
    atomic_write_dataframe_csv(output_dir / "exact_graphlet_values.csv", pd.DataFrame(value_rows))
    atomic_write_dataframe_csv(output_dir / "exact_graphlet_metrics.csv", metric_df)
    atomic_write_dataframe_csv(output_dir / "exact_model_failures.csv", failure_df)
    atomic_write_dataframe_csv(output_dir / "exact_model_distances.csv", distance_df)
    atomic_write_dataframe_csv(output_dir / "abc_descriptive_gate.csv", recommended_df)

    summary = {
        "experiment_id": cfg["experiment"]["id"],
        "assay_id": assay_id,
        "source_full_abc_run": cfg["source_full_abc"].get("run_id"),
        "source_full_abc_summary_path": str(source_summary_path),
        "source_full_abc_config_path": str(source_config_path),
        "n_prepared_rows": int(len(assay_df)),
        "n_mutation_sets": int(len(score_map)),
        "score_columns": available_score_columns,
        "primary_score_column": primary_score,
        "n_complete_graphlets": int(len(graphlets)),
        "graphlet_dimensions_observed": sorted({int(graphlet.dimension) for graphlet in graphlets}),
        "metric_columns": metric_columns,
        "descriptive_gate": recommended_rows,
        "best_distance_rows": distance_df.to_dict(orient="records") if not distance_df.empty else [],
        "output_files": {
            "graphlet_catalog": str(output_dir / "graphlet_catalog.csv"),
            "exact_graphlet_values": str(output_dir / "exact_graphlet_values.csv"),
            "exact_graphlet_metrics": str(output_dir / "exact_graphlet_metrics.csv"),
            "exact_model_failures": str(output_dir / "exact_model_failures.csv"),
            "exact_model_distances": str(output_dir / "exact_model_distances.csv"),
            "abc_descriptive_gate": str(output_dir / "abc_descriptive_gate.csv"),
        },
    }
    atomic_write_json(output_dir / "summary.json", summary)
    emit_progress(
        output_dir,
        {
            "event": "complete",
            "n_complete_graphlets": len(graphlets),
            "n_metric_rows": len(metric_rows),
            "n_failure_rows": int(len(failure_df)),
        },
    )
    print(json.dumps(to_builtin(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
