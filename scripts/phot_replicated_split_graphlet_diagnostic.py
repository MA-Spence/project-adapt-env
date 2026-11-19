#!/usr/bin/env python3
"""Replicated PHOT held-out graphlet split diagnostic.

This diagnostic reuses the fixed RUN-079 candidate graphlet metrics and
repeatedly resamples calibration/test graphlet splits. It asks whether the
RUN-079 no-module coupled-landscape advantage is stable under calibration
selection, without introducing a new operator or refitting the PHOT base.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if quick and "quick" in config:
        config = merge_dict(config, config["quick"])
    return config


def merge_dict(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in updates.items():
        if key == "quick":
            continue
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_dict(dict(merged[key]), value)
        else:
            merged[key] = value
    return merged


def finite_mean(rows: pd.DataFrame, metric: str) -> float:
    if metric not in rows:
        return math.nan
    values = pd.to_numeric(rows[metric], errors="coerce").dropna()
    return float(values.mean()) if not values.empty else math.nan


def metric_distance(
    model_rows: pd.DataFrame,
    target_rows: pd.DataFrame,
    metrics: list[str],
) -> tuple[float, dict[str, float], dict[str, float], dict[str, float]]:
    pieces: list[float] = []
    by_metric: dict[str, float] = {}
    model_means: dict[str, float] = {}
    target_means: dict[str, float] = {}
    for metric in metrics:
        model_value = finite_mean(model_rows, metric)
        target_value = finite_mean(target_rows, metric)
        target_values = pd.to_numeric(target_rows.get(metric), errors="coerce").dropna()
        if not np.isfinite(model_value) or not np.isfinite(target_value):
            continue
        scale = float(target_values.std(ddof=0)) if len(target_values) >= 2 else 0.0
        scale = max(scale, abs(target_value) * 0.10, 1e-6)
        value = abs(model_value - target_value) / scale
        by_metric[metric] = float(value)
        model_means[metric] = float(model_value)
        target_means[metric] = float(target_value)
        pieces.append(float(value))
    return (
        float(np.mean(pieces)) if pieces else math.nan,
        by_metric,
        model_means,
        target_means,
    )


def candidate_family(candidate_model: str, config: dict[str, Any]) -> str | None:
    mapping = dict(config["candidate_families"])
    return mapping.get(str(candidate_model))


def split_rows(
    graphlet_catalog: pd.DataFrame,
    *,
    rng: np.random.Generator,
    split_id: str,
    split_kind: str,
    calibration_fraction: float,
    stratify_column: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if split_kind == "original":
        for item in graphlet_catalog.to_dict(orient="records"):
            rows.append(
                {
                    "split_id": split_id,
                    "split_kind": split_kind,
                    "graphlet_id": item["graphlet_id"],
                    "graphlet_dimension": item.get("graphlet_dimension", ""),
                    "assigned_split": item["graphlet_split"],
                }
            )
        return rows

    for _key, sub in graphlet_catalog.groupby(stratify_column, dropna=False):
        ids = [str(item) for item in sub["graphlet_id"].tolist()]
        rng.shuffle(ids)
        n_calibration = int(round(len(ids) * calibration_fraction))
        n_calibration = min(max(n_calibration, 1), len(ids) - 1)
        calibration_ids = set(ids[:n_calibration])
        for item in sub.to_dict(orient="records"):
            graphlet_id = str(item["graphlet_id"])
            rows.append(
                {
                    "split_id": split_id,
                    "split_kind": split_kind,
                    "graphlet_id": graphlet_id,
                    "graphlet_dimension": item.get("graphlet_dimension", ""),
                    "assigned_split": "calibration" if graphlet_id in calibration_ids else "test",
                }
            )
    return rows


def distances_for_split(
    metric_df: pd.DataFrame,
    assignment: pd.DataFrame,
    *,
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    metric_cols = [str(item) for item in config["graphlet_evaluation"]["metric_columns"]]
    family_order = [str(item) for item in config["diagnostic"]["family_order"]]
    split_id = str(assignment["split_id"].iloc[0])
    split_kind = str(assignment["split_kind"].iloc[0])
    assigned = assignment[["graphlet_id", "assigned_split"]]
    split_metric_df = metric_df.merge(assigned, on="graphlet_id", how="inner")
    selected_rows: list[dict[str, Any]] = []

    for family in family_order:
        family_models = [
            model
            for model, mapped_family in config["candidate_families"].items()
            if str(mapped_family) == family
        ]
        family_df = split_metric_df[
            (split_metric_df["score_column"] == primary_score)
            & (split_metric_df["candidate_model"].isin(family_models))
        ]
        if family_df.empty:
            continue
        calibration_distances: list[dict[str, Any]] = []
        for keys, model_rows in family_df.groupby(["candidate_model", "parameter_set_id"], dropna=False):
            target_rows = split_metric_df[
                (split_metric_df["score_column"] == primary_score)
                & (split_metric_df["candidate_model"] == f"empirical_observed__{primary_score}")
                & (split_metric_df["assigned_split"] == "calibration")
            ]
            candidate_rows = model_rows[model_rows["assigned_split"] == "calibration"]
            distance, by_metric, _model_means, _target_means = metric_distance(candidate_rows, target_rows, metric_cols)
            first = model_rows.iloc[0].to_dict()
            calibration_distances.append(
                {
                    "candidate_model": keys[0],
                    "parameter_set_id": keys[1],
                    "calibration_distance": distance,
                    "base_anchor": first.get("base_anchor", ""),
                    "module_source": first.get("module_source", ""),
                    **{f"calibration_component__{key}": value for key, value in by_metric.items()},
                }
            )
        valid = [row for row in calibration_distances if np.isfinite(float(row["calibration_distance"]))]
        if not valid:
            continue
        selected = sorted(valid, key=lambda row: float(row["calibration_distance"]))[0]
        selected_metric_rows = family_df[
            (family_df["candidate_model"] == selected["candidate_model"])
            & (family_df["parameter_set_id"] == selected["parameter_set_id"])
        ]
        output_row = {
            "split_id": split_id,
            "split_kind": split_kind,
            "candidate_family": family,
            **selected,
        }
        for assigned_split in ["test", "all"]:
            target_rows = split_metric_df[
                (split_metric_df["score_column"] == primary_score)
                & (split_metric_df["candidate_model"] == f"empirical_observed__{primary_score}")
                & ((split_metric_df["assigned_split"] == assigned_split) if assigned_split != "all" else True)
            ]
            candidate_rows = selected_metric_rows[
                (selected_metric_rows["assigned_split"] == assigned_split) if assigned_split != "all" else slice(None)
            ]
            distance, by_metric, model_means, target_means = metric_distance(candidate_rows, target_rows, metric_cols)
            output_row[f"{assigned_split}_distance"] = distance
            for key, value in by_metric.items():
                output_row[f"{assigned_split}_component__{key}"] = value
            for key, value in model_means.items():
                output_row[f"{assigned_split}_model_mean__{key}"] = value
            for key, value in target_means.items():
                output_row[f"{assigned_split}_empirical_mean__{key}"] = value
                output_row[f"{assigned_split}_raw_delta__{key}"] = model_means.get(key, math.nan) - value
        selected_rows.append(output_row)
    return selected_rows


def summarize_family(rows: pd.DataFrame, *, primary_families: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    random_rows = rows[rows["split_kind"] == "random"].copy()
    if random_rows.empty:
        return out
    pivot = random_rows.pivot(index="split_id", columns="candidate_family", values="test_distance")
    for family, sub in random_rows.groupby("candidate_family", dropna=False):
        values = pd.to_numeric(sub["test_distance"], errors="coerce").dropna()
        if values.empty:
            continue
        row: dict[str, Any] = {
            "candidate_family": family,
            "is_primary_family": str(family) in primary_families,
            "n_random_splits": int(values.shape[0]),
            "mean_test_distance": float(values.mean()),
            "median_test_distance": float(values.median()),
            "sd_test_distance": float(values.std(ddof=0)),
            "q10_test_distance": float(values.quantile(0.10)),
            "q90_test_distance": float(values.quantile(0.90)),
        }
        if "fitted_current" in pivot and family in pivot:
            delta = pivot[family] - pivot["fitted_current"]
            row["mean_delta_vs_fitted_current"] = float(delta.mean())
            row["median_delta_vs_fitted_current"] = float(delta.median())
            row["win_rate_vs_fitted_current"] = float((delta < 0).mean())
        if "smooth_core" in pivot and family in pivot:
            delta = pivot[family] - pivot["smooth_core"]
            row["mean_delta_vs_smooth_core"] = float(delta.mean())
            row["median_delta_vs_smooth_core"] = float(delta.median())
            row["win_rate_vs_smooth_core"] = float((delta < 0).mean())
        if family in pivot:
            row["best_family_rate"] = float((pivot.idxmin(axis=1) == family).mean())
        out.append(row)
    return out


def pairwise_deltas(rows: pd.DataFrame, *, comparisons: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    pivot = rows.pivot(index=["split_id", "split_kind"], columns="candidate_family", values="test_distance")
    for family in pivot.columns:
        for comparison in comparisons:
            if family == comparison or comparison not in pivot.columns:
                continue
            delta = pivot[family] - pivot[comparison]
            for key, value in delta.dropna().items():
                out.append(
                    {
                        "split_id": key[0],
                        "split_kind": key[1],
                        "candidate_family": family,
                        "comparison_family": comparison,
                        "delta_test_distance": float(value),
                        "beats_comparison": bool(value < 0),
                    }
                )
    return out


def summarize_components(rows: pd.DataFrame, *, metrics: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    random_rows = rows[rows["split_kind"] == "random"].copy()
    for family, sub in random_rows.groupby("candidate_family", dropna=False):
        for metric in metrics:
            component = pd.to_numeric(sub.get(f"test_component__{metric}"), errors="coerce").dropna()
            raw_delta = pd.to_numeric(sub.get(f"test_raw_delta__{metric}"), errors="coerce").dropna()
            model_mean = pd.to_numeric(sub.get(f"test_model_mean__{metric}"), errors="coerce").dropna()
            empirical_mean = pd.to_numeric(sub.get(f"test_empirical_mean__{metric}"), errors="coerce").dropna()
            if component.empty:
                continue
            out.append(
                {
                    "candidate_family": family,
                    "metric": metric,
                    "median_component_distance": float(component.median()),
                    "mean_component_distance": float(component.mean()),
                    "median_raw_delta": float(raw_delta.median()) if not raw_delta.empty else math.nan,
                    "mean_raw_delta": float(raw_delta.mean()) if not raw_delta.empty else math.nan,
                    "median_model_mean": float(model_mean.median()) if not model_mean.empty else math.nan,
                    "median_empirical_mean": float(empirical_mean.median()) if not empirical_mean.empty else math.nan,
                }
            )
    return out


def selection_counts(rows: pd.DataFrame) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    random_rows = rows[rows["split_kind"] == "random"].copy()
    for keys, sub in random_rows.groupby(["candidate_family", "candidate_model", "parameter_set_id"], dropna=False):
        family = str(keys[0])
        family_total = int((random_rows["candidate_family"] == family).sum())
        out.append(
            {
                "candidate_family": family,
                "candidate_model": keys[1],
                "parameter_set_id": keys[2],
                "n_selected_random_splits": int(len(sub)),
                "selection_rate_within_family": float(len(sub) / family_total) if family_total else math.nan,
                "median_test_distance": float(pd.to_numeric(sub["test_distance"], errors="coerce").median()),
            }
        )
    return out


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config = load_config(Path(args.config), quick=bool(args.quick))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    source_dir = project_root / str(config["source_run"]["output_dir"])
    graphlet_metrics_path = source_dir / "graphlet_metrics.csv"
    graphlet_catalog_path = source_dir / "graphlet_catalog.csv"
    parameter_grid_path = source_dir / "parameter_grid.csv"
    if not graphlet_metrics_path.exists():
        raise FileNotFoundError(f"Missing RUN-079 graphlet metrics: {graphlet_metrics_path}")
    if not graphlet_catalog_path.exists():
        raise FileNotFoundError(f"Missing RUN-079 graphlet catalog: {graphlet_catalog_path}")

    metric_df = pd.read_csv(graphlet_metrics_path)
    graphlet_catalog = pd.read_csv(graphlet_catalog_path)
    if parameter_grid_path.exists():
        parameter_grid = pd.read_csv(parameter_grid_path)
    else:
        parameter_grid = pd.DataFrame()

    metric_df["graphlet_id"] = metric_df["graphlet_id"].astype(str)
    graphlet_catalog["graphlet_id"] = graphlet_catalog["graphlet_id"].astype(str)
    diagnostic_cfg = config["diagnostic"]
    rng = np.random.default_rng(int(diagnostic_cfg["seed"]))
    n_splits = int(diagnostic_cfg["n_random_splits"])
    calibration_fraction = float(diagnostic_cfg["calibration_fraction"])
    stratify_column = str(diagnostic_cfg["stratify_column"])

    assignments: list[dict[str, Any]] = []
    selected_distances: list[dict[str, Any]] = []
    if bool(diagnostic_cfg.get("include_original_split", True)):
        original = split_rows(
            graphlet_catalog,
            rng=rng,
            split_id="original",
            split_kind="original",
            calibration_fraction=calibration_fraction,
            stratify_column=stratify_column,
        )
        assignments.extend(original)
        selected_distances.extend(distances_for_split(metric_df, pd.DataFrame(original), config=config))

    for idx in range(n_splits):
        split_id = f"split_{idx:04d}"
        rows = split_rows(
            graphlet_catalog,
            rng=rng,
            split_id=split_id,
            split_kind="random",
            calibration_fraction=calibration_fraction,
            stratify_column=stratify_column,
        )
        assignments.extend(rows)
        selected_distances.extend(distances_for_split(metric_df, pd.DataFrame(rows), config=config))

    selected_df = pd.DataFrame(selected_distances)
    primary_families = [str(item) for item in diagnostic_cfg["primary_families"]]
    comparison_families = [str(item) for item in diagnostic_cfg["comparison_families"]]
    metric_cols = [str(item) for item in config["graphlet_evaluation"]["metric_columns"]]

    family_summary = summarize_family(selected_df, primary_families=primary_families)
    pairwise_rows = pairwise_deltas(selected_df, comparisons=comparison_families)
    component_rows = summarize_components(selected_df, metrics=metric_cols)
    counts = selection_counts(selected_df)

    write_csv(output_dir / "split_assignments.csv", assignments)
    write_csv(output_dir / "split_selected_distances.csv", selected_distances)
    write_csv(output_dir / "split_family_summary.csv", family_summary)
    write_csv(output_dir / "split_family_pairwise_deltas.csv", pairwise_rows)
    write_csv(output_dir / "component_summary.csv", component_rows)
    write_csv(output_dir / "selection_counts.csv", counts)

    family_summary_df = pd.DataFrame(family_summary)
    no_module_row = family_summary_df[family_summary_df["candidate_family"] == "no_module"]
    preflight = {
        "experiment_id": config["experiment_id"],
        "source_run_id": config["source_run"]["run_id"],
        "source_output_dir": str(source_dir),
        "n_graphlets": int(graphlet_catalog["graphlet_id"].nunique()),
        "n_graphlet_metric_rows": int(len(metric_df)),
        "n_candidate_parameter_rows": int(len(parameter_grid)),
        "n_random_splits": n_splits,
        "include_original_split": bool(diagnostic_cfg.get("include_original_split", True)),
        "calibration_fraction": calibration_fraction,
        "stratify_column": stratify_column,
        "candidate_families": config["candidate_families"],
        "metric_columns": metric_cols,
    }
    summary = {
        **preflight,
        "outputs": {
            "preflight": str(output_dir / "preflight.json"),
            "split_assignments": str(output_dir / "split_assignments.csv"),
            "split_selected_distances": str(output_dir / "split_selected_distances.csv"),
            "split_family_summary": str(output_dir / "split_family_summary.csv"),
            "split_family_pairwise_deltas": str(output_dir / "split_family_pairwise_deltas.csv"),
            "component_summary": str(output_dir / "component_summary.csv"),
            "selection_counts": str(output_dir / "selection_counts.csv"),
        },
        "primary_diagnostic": (
            no_module_row.to_dict(orient="records")[0] if not no_module_row.empty else {}
        ),
        "family_summary": family_summary,
    }
    write_json(output_dir / "preflight.json", preflight)
    write_json(output_dir / "summary.json", summary)


if __name__ == "__main__":
    main()
