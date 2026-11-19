"""Run empirical graphlet falsification checks for the residual-ruggedness model.

This wraps the TEM-1 aggregate graphlet runner used in RUN-067/RUN-069 and adds
an empirical ProteinGym graphlet panel. The goal is not to fit a new model, but
to test the explicit falsifiers recorded in scratch/model_literature_falsification_review.md.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

import tem1_aggregate_graphlet_residual_ruggedness as graphlet_runner


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(Path.cwd()))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=graphlet_runner.json_default),
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def mutation_set(wildtype: str, sequence: str) -> frozenset[tuple[int, str]]:
    return frozenset(
        (idx, aa)
        for idx, (wt_aa, aa) in enumerate(zip(wildtype, str(sequence)))
        if aa != wt_aa
    )


def subset_masks(features: tuple[tuple[int, str], ...]) -> list[frozenset[tuple[int, str]]]:
    masks: list[frozenset[tuple[int, str]]] = []
    for size in range(len(features) + 1):
        for sub in itertools.combinations(features, size):
            masks.append(frozenset(sub))
    return masks


def score_map_for_assay(
    assay_df: pd.DataFrame,
    *,
    wildtype: str,
    wt_scores: dict[str, float],
    score_columns: list[str],
) -> dict[frozenset[tuple[int, str]], dict[str, float]]:
    grouped: dict[frozenset[tuple[int, str]], dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in assay_df.itertuples(index=False):
        seq = str(getattr(row, "mutated_sequence"))
        if len(seq) != len(wildtype):
            continue
        key = mutation_set(wildtype, seq)
        for col in score_columns:
            value = getattr(row, col, math.nan)
            if np.isfinite(float(value)):
                grouped[key][col].append(float(value))
    empty = frozenset()
    for col, value in wt_scores.items():
        if col in score_columns and np.isfinite(float(value)):
            grouped[empty][col].append(float(value))
    out: dict[frozenset[tuple[int, str]], dict[str, float]] = {}
    for key, values_by_col in grouped.items():
        out[key] = {
            col: float(np.mean(values))
            for col, values in values_by_col.items()
            if values
        }
    return out


def feature_frequency(score_map: dict[frozenset[tuple[int, str]], dict[str, float]]) -> Counter:
    counts: Counter = Counter()
    for feature_set in score_map:
        for feature in feature_set:
            counts[feature] += 1
    return counts


def find_complete_graphlets(
    score_map: dict[frozenset[tuple[int, str]], dict[str, float]],
    *,
    score_columns: list[str],
    dimensions: list[int],
    max_per_dimension: int,
    top_feature_pool: int,
) -> list[tuple[int, tuple[tuple[int, str], ...]]]:
    counts = feature_frequency(score_map)
    top_features = [feature for feature, _count in counts.most_common(top_feature_pool)]
    found: list[tuple[int, tuple[tuple[int, str], ...]]] = []
    used: set[tuple[tuple[int, str], ...]] = set()
    for dim in sorted(set(int(x) for x in dimensions)):
        dim_found = 0
        pool = top_features[: max(top_feature_pool, dim)]
        for combo in itertools.combinations(pool, dim):
            combo = tuple(combo)
            if combo in used:
                continue
            ok = True
            for sub in subset_masks(combo):
                values = score_map.get(sub)
                if values is None or not all(col in values for col in score_columns):
                    ok = False
                    break
            if not ok:
                continue
            found.append((dim, combo))
            used.add(combo)
            dim_found += 1
            if dim_found >= max_per_dimension:
                break
    return found


def walsh_high_fraction(values: np.ndarray, dim: int) -> float:
    _rows, fractions = graphlet_runner.walsh_decomposition(values, dim=dim)
    return float(sum(value for order, value in fractions.items() if int(order) >= 3))


def additive_projection(values: np.ndarray, dim: int) -> np.ndarray:
    rows, _fractions = graphlet_runner.walsh_decomposition(values, dim=dim)
    coeffs = {int(row["mask"]): float(row["coefficient"]) for row in rows if int(row["order"]) <= 1}
    out = np.zeros_like(values, dtype=np.float64)
    for mask, coeff in coeffs.items():
        for node in range(values.size):
            parity = (node & mask).bit_count() % 2
            out[node] += coeff * (-1.0 if parity else 1.0)
    return out


def fit_global_readout(values: np.ndarray, dim: int, degree: int) -> tuple[np.ndarray, np.ndarray, float]:
    additive = additive_projection(values, dim)
    x = additive.astype(np.float64)
    if float(np.std(x)) < 1e-12:
        pred = np.full_like(values, float(np.mean(values)), dtype=np.float64)
    else:
        x0 = (x - float(np.mean(x))) / max(float(np.std(x)), 1e-12)
        columns = [np.ones_like(x0)]
        for power in range(1, int(degree) + 1):
            columns.append(x0**power)
        design = np.vstack(columns).T
        beta, *_rest = np.linalg.lstsq(design, values.astype(np.float64), rcond=None)
        pred = design @ beta
    resid = values.astype(np.float64) - pred
    sst = float(np.sum((values - float(np.mean(values))) ** 2))
    ssr = float(np.sum((values - pred) ** 2))
    r2 = 1.0 - ssr / sst if sst > 1e-12 else 1.0
    return pred.astype(np.float64), resid.astype(np.float64), float(r2)


def metric_rows_for_values(
    *,
    graphlet: dict[str, Any],
    values: np.ndarray,
    candidate_model: str,
    parameter_set_id: str,
    cfg: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    dim = int(graphlet["graphlet_dimension"])
    walsh_rows, variance_by_order = graphlet_runner.walsh_decomposition(values, dim=dim)
    total_variance = float(np.var(values))
    sign_counts = graphlet_runner.sign_epistasis_counts(values, dim=dim)
    extrema = graphlet_runner.local_extrema_counts(values, dim=dim)
    access = graphlet_runner.accessible_path_metrics(values, dim=dim)
    basins = graphlet_runner.basin_metrics(values, dim=dim)
    random_walk = graphlet_runner.graph_random_walk_metrics(
        values,
        dim=dim,
        cfg=cfg,
        rng=np.random.RandomState(graphlet_runner.stable_seed("empirical-rw", graphlet["graphlet_id"], candidate_model)),
    )
    parameter_set = {"candidate_model": candidate_model, "parameter_set_id": parameter_set_id}
    gft_rows, gft_summary = graphlet_runner.gft_rows_from_walsh(
        walsh_rows,
        graphlet=graphlet,
        parameter_set=parameter_set,
        total_variance=total_variance,
    )
    high_order_fraction = float(
        sum(float(value) for order, value in variance_by_order.items() if int(order) >= 3)
    )
    metric = {
        **graphlet,
        **parameter_set,
        "split": "empirical",
        "fitness_min": float(np.min(values)),
        "fitness_mean": float(np.mean(values)),
        "fitness_max": float(np.max(values)),
        "fitness_sd": float(np.std(values)),
        "total_variance": total_variance,
        "walsh_high_order_fraction": high_order_fraction,
        "walsh_high_order_present": bool(high_order_fraction > 1e-6),
        **{f"walsh_variance_order_{order}": float(value) for order, value in variance_by_order.items()},
        **{key: value for key, value in sign_counts.items() if key != "pair_rows"},
        **{key: value for key, value in extrema.items() if not isinstance(value, list)},
        "local_maxima": ";".join(str(x) for x in extrema["local_maxima"]),
        "local_minima": ";".join(str(x) for x in extrema["local_minima"]),
        **access,
        **basins,
        **random_walk,
        **gft_summary,
    }
    for row in walsh_rows:
        row.update(
            {
                "graphlet_id": graphlet["graphlet_id"],
                "split": "empirical",
                "candidate_model": candidate_model,
                "parameter_set_id": parameter_set_id,
                "assay_id": graphlet["assay_id"],
                "score_column": graphlet["score_column"],
            }
        )
    return metric, walsh_rows, gft_rows, access


def summarize_rows(rows: list[dict[str, Any]], group_cols: list[str], metric_cols: list[str]) -> list[dict[str, Any]]:
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
            row[f"{col}_median"] = float(values.median())
            row[f"{col}_q10"] = float(values.quantile(0.10))
            row[f"{col}_q90"] = float(values.quantile(0.90))
            row[f"{col}_sd"] = float(values.std(ddof=0))
        out.append(row)
    return out


def spearman(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2 or len(xs) != len(ys):
        return math.nan
    xr = pd.Series(xs).rank(method="average").to_numpy(dtype=np.float64)
    yr = pd.Series(ys).rank(method="average").to_numpy(dtype=np.float64)
    if float(np.std(xr)) < 1e-12 or float(np.std(yr)) < 1e-12:
        return math.nan
    return float(np.corrcoef(xr, yr)[0, 1])


def empirical_graphlet_analysis(cfg: dict[str, Any], *, project_root: Path, output_dir: Path) -> dict[str, Any]:
    empirical_cfg = cfg["empirical_graphlets"]
    metric_cols = list(empirical_cfg["metric_columns"])
    score_columns = list(empirical_cfg["score_columns"])
    random_walk_cfg = {"analysis": cfg.get("analysis", {})}

    graphlet_catalog: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    access_rows: list[dict[str, Any]] = []
    readout_rows: list[dict[str, Any]] = []
    assay_rows: list[dict[str, Any]] = []

    for assay_cfg in empirical_cfg["assays"]:
        assay_id = str(assay_cfg["assay_id"])
        panel_path = project_root / str(assay_cfg["selected_panel_path"])
        assay_path = project_root / str(assay_cfg["prepared_assay_path"])
        metrics_path = project_root / str(assay_cfg["mavenn_metrics_path"])
        panel = pd.read_csv(panel_path)
        wildtype = str(panel.loc[panel["DMS_id"] == assay_id, "target_seq"].iloc[0])
        assay_df = pd.read_csv(assay_path)
        mavenn_metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
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
        dims = [int(x) for x in assay_cfg.get("dimensions", empirical_cfg["dimensions"])]
        complete = find_complete_graphlets(
            score_map,
            score_columns=assay_score_columns,
            dimensions=dims,
            max_per_dimension=int(assay_cfg.get("max_graphlets_per_dimension", empirical_cfg["max_graphlets_per_assay_per_dimension"])),
            top_feature_pool=int(assay_cfg.get("top_feature_pool", empirical_cfg["top_feature_pool"])),
        )
        assay_rows.append(
            {
                "assay_id": assay_id,
                "prepared_assay_path": str(assay_path),
                "n_rows": int(len(assay_df)),
                "n_mutation_sets": int(len(score_map)),
                "n_complete_graphlets": int(len(complete)),
                "dimensions_requested": ";".join(str(x) for x in dims),
            }
        )
        for graphlet_index, (dim, features) in enumerate(complete):
            base_graphlet = {
                "graphlet_id": f"{assay_id}__g{graphlet_index:04d}__d{dim}",
                "assay_id": assay_id,
                "split": "empirical",
                "graphlet_dimension": int(dim),
                "mutation_features": ";".join(f"{pos}:{aa}" for pos, aa in features),
            }
            graphlet_catalog.append(base_graphlet)
            masks = subset_masks(features)
            for score_col in assay_score_columns:
                values = np.asarray([score_map[mask][score_col] for mask in masks], dtype=np.float64)
                graphlet = {**base_graphlet, "score_column": score_col}
                metric, walsh, gft, access = metric_rows_for_values(
                    graphlet=graphlet,
                    values=values,
                    candidate_model=f"empirical_observed__{score_col}",
                    parameter_set_id=score_col,
                    cfg=random_walk_cfg,
                )
                metric_rows.append(metric)
                walsh_rows.extend(walsh)
                gft_rows.extend(gft)
                access_rows.append({"graphlet_id": graphlet["graphlet_id"], "score_column": score_col, **access})

                pred, resid, r2 = fit_global_readout(
                    values,
                    dim=dim,
                    degree=int(empirical_cfg.get("global_readout_polynomial_degree", 3)),
                )
                fit_metric, fit_walsh, fit_gft, _fit_access = metric_rows_for_values(
                    graphlet=graphlet,
                    values=pred,
                    candidate_model=f"empirical_additive_global_readout_fit__{score_col}",
                    parameter_set_id=f"{score_col}__poly",
                    cfg=random_walk_cfg,
                )
                residual_metric, residual_walsh, residual_gft, _resid_access = metric_rows_for_values(
                    graphlet=graphlet,
                    values=resid,
                    candidate_model=f"empirical_global_readout_residual__{score_col}",
                    parameter_set_id=f"{score_col}__residual",
                    cfg=random_walk_cfg,
                )
                metric_rows.extend([fit_metric, residual_metric])
                walsh_rows.extend(fit_walsh)
                walsh_rows.extend(residual_walsh)
                gft_rows.extend(fit_gft)
                gft_rows.extend(residual_gft)
                observed_high = walsh_high_fraction(values, dim)
                residual_high = walsh_high_fraction(resid, dim)
                readout_rows.append(
                    {
                        **graphlet,
                        "global_readout_r2": r2,
                        "observed_high_order_fraction": observed_high,
                        "post_global_readout_residual_high_order_fraction": residual_high,
                        "high_order_fraction_retained_after_readout_fit": (
                            residual_high / observed_high if observed_high > 1e-12 else math.nan
                        ),
                    }
                )

    summary_rows = summarize_rows(
        metric_rows,
        group_cols=["candidate_model", "parameter_set_id", "score_column"],
        metric_cols=metric_cols,
    )
    readout_summary = summarize_rows(
        readout_rows,
        group_cols=["score_column"],
        metric_cols=[
            "global_readout_r2",
            "observed_high_order_fraction",
            "post_global_readout_residual_high_order_fraction",
            "high_order_fraction_retained_after_readout_fit",
        ],
    )

    write_csv(output_dir / "empirical_assay_audit.csv", assay_rows)
    write_csv(output_dir / "empirical_graphlet_catalog.csv", graphlet_catalog)
    write_csv(output_dir / "empirical_graphlet_metrics.csv", metric_rows)
    write_csv(output_dir / "empirical_walsh_terms.csv", walsh_rows)
    write_csv(output_dir / "empirical_gft_coefficients.csv", gft_rows)
    write_csv(output_dir / "empirical_graphlet_accessibility.csv", access_rows)
    write_csv(output_dir / "empirical_global_readout_confounding.csv", readout_rows)
    write_csv(output_dir / "empirical_graphlet_summary.csv", summary_rows)
    write_csv(output_dir / "empirical_global_readout_summary.csv", readout_summary)

    return {
        "n_assays": int(len(empirical_cfg["assays"])),
        "n_empirical_graphlets": int(len(graphlet_catalog)),
        "n_empirical_metric_rows": int(len(metric_rows)),
        "assay_audit": assay_rows,
        "summary_rows": summary_rows,
        "readout_summary": readout_summary,
    }


def metric_distance(row: pd.Series, empirical_target: dict[str, float], empirical_scale: dict[str, float], metrics: list[str]) -> float:
    pieces: list[float] = []
    for metric in metrics:
        model_col = f"{metric}_mean"
        if model_col not in row or metric not in empirical_target:
            continue
        model_value = float(row[model_col])
        target_value = float(empirical_target[metric])
        scale = max(float(empirical_scale.get(metric, 0.0)), abs(target_value) * 0.10, 1e-6)
        pieces.append(abs(model_value - target_value) / scale)
    return float(np.mean(pieces)) if pieces else math.nan


def compare_synthetic_to_empirical(
    cfg: dict[str, Any],
    *,
    synthetic_dir: Path,
    empirical_dir: Path,
    output_dir: Path,
) -> dict[str, Any]:
    metrics = list(cfg["falsification_tests"]["comparison_metrics"])
    primary_score = str(cfg["empirical_graphlets"]["primary_score_column"])
    synthetic = pd.read_csv(synthetic_dir / "module_parameter_sweep.csv")
    synthetic_metrics = pd.read_csv(synthetic_dir / "test_graphlet_metrics.csv")
    empirical_metrics = pd.read_csv(empirical_dir / "empirical_graphlet_metrics.csv")
    empirical_obs = empirical_metrics[
        (empirical_metrics["candidate_model"] == f"empirical_observed__{primary_score}")
        & (empirical_metrics["score_column"] == primary_score)
        & (empirical_metrics["graphlet_dimension"] >= int(cfg["empirical_graphlets"].get("primary_min_dimension", 3)))
    ].copy()
    if empirical_obs.empty:
        empirical_obs = empirical_metrics[
            (empirical_metrics["candidate_model"] == f"empirical_observed__{primary_score}")
            & (empirical_metrics["score_column"] == primary_score)
        ].copy()
    empirical_target = {
        metric: float(pd.to_numeric(empirical_obs[metric], errors="coerce").mean())
        for metric in metrics
        if metric in empirical_obs
    }
    empirical_scale = {
        metric: float(pd.to_numeric(empirical_obs[metric], errors="coerce").std(ddof=0))
        for metric in metrics
        if metric in empirical_obs
    }
    rows: list[dict[str, Any]] = []
    test_rows = synthetic[synthetic["split"] == "test"].copy()
    for _idx, row in test_rows.iterrows():
        rows.append(
            {
                "split": row["split"],
                "candidate_model": row["candidate_model"],
                "parameter_set_id": row["parameter_set_id"],
                "distance_to_empirical_primary": metric_distance(row, empirical_target, empirical_scale, metrics),
                **{f"model_{metric}_mean": float(row.get(f"{metric}_mean", math.nan)) for metric in metrics},
                **{f"empirical_{metric}_mean": empirical_target.get(metric, math.nan) for metric in metrics},
            }
        )
    write_csv(output_dir / "synthetic_empirical_model_distances.csv", rows)

    distance_df = pd.DataFrame(rows)
    best_by_model = (
        distance_df.sort_values("distance_to_empirical_primary")
        .groupby("candidate_model", as_index=False)
        .first()
        if not distance_df.empty
        else pd.DataFrame()
    )
    best_rows = best_by_model.to_dict(orient="records") if not best_by_model.empty else []
    write_csv(output_dir / "synthetic_empirical_best_by_model.csv", best_rows)
    return {
        "empirical_primary_score_column": primary_score,
        "empirical_target_means": empirical_target,
        "empirical_target_scales": empirical_scale,
        "best_by_model": best_rows,
    }


def evaluate_falsifiers(
    cfg: dict[str, Any],
    *,
    synthetic_dir: Path,
    empirical_dir: Path,
    comparison: dict[str, Any],
) -> dict[str, Any]:
    tests_cfg = cfg["falsification_tests"]
    synthetic = pd.read_csv(synthetic_dir / "module_parameter_sweep.csv")
    empirical_summary = pd.read_csv(empirical_dir / "empirical_graphlet_summary.csv")
    readout_summary = pd.read_csv(empirical_dir / "empirical_global_readout_summary.csv")
    synthetic_metrics = pd.read_csv(synthetic_dir / "test_graphlet_metrics.csv")
    primary_score = str(cfg["empirical_graphlets"]["primary_score_column"])

    out: dict[str, Any] = {}

    scale_rows = synthetic[
        (synthetic["split"] == "test")
        & (synthetic["candidate_model"] == "contact_residual_ensemble")
    ].copy()
    scale_lookup = (
        synthetic_metrics[synthetic_metrics["candidate_model"] == "contact_residual_ensemble"]
        .groupby("parameter_set_id")["residual_scale_multiplier"]
        .mean()
        .to_dict()
        if "residual_scale_multiplier" in synthetic_metrics.columns
        else {}
    )
    scale_values: list[float] = []
    for value in scale_rows["parameter_set_id"]:
        if value in scale_lookup:
            scale_values.append(float(scale_lookup[value]))
        else:
            # RUN-069 style pset ids do not embed the scale; fall back to row order.
            scale_values.append(float(len(scale_values) + 1))
    scale_effects = {}
    direction_cfg = tests_cfg.get("scale_effect_expected_direction", {})
    directional_passes = 0
    for metric in tests_cfg["scale_effect_metrics"]:
        col = f"{metric}_mean"
        if col in scale_rows:
            rho = spearman(scale_values, [float(x) for x in scale_rows[col]])
            expected_direction = float(direction_cfg.get(metric, 1.0))
            scale_effects[metric] = {
                "spearman": rho,
                "expected_direction": expected_direction,
                "directional_score": (
                    float(rho * expected_direction) if np.isfinite(rho) else math.nan
                ),
            }
            if np.isfinite(rho) and rho * expected_direction > 0.5:
                directional_passes += 1
    out["no_residual_scale_effect"] = {
        "falsifier": "Residual scale is not directional across the fixed-scale grid.",
        "scale_values_by_parameter_set": {
            str(parameter_set_id): float(scale)
            for parameter_set_id, scale in zip(scale_rows["parameter_set_id"], scale_values)
        },
        "spearman_by_metric": scale_effects,
        "passes_directional_screen": bool(directional_passes >= int(tests_cfg.get("min_directional_scale_metrics", 4))),
    }

    primary_emp = empirical_summary[
        (empirical_summary["candidate_model"] == f"empirical_observed__{primary_score}")
        & (empirical_summary["score_column"] == primary_score)
    ]
    spectral = {}
    if not primary_emp.empty:
        row = primary_emp.iloc[0]
        spectral = {
            "gft_low_energy_fraction_mean": float(row.get("gft_low_energy_fraction_mean", math.nan)),
            "gft_high_energy_fraction_mean": float(row.get("gft_high_energy_fraction_mean", math.nan)),
            "walsh_high_order_fraction_mean": float(row.get("walsh_high_order_fraction_mean", math.nan)),
        }
    out["wrong_spectral_envelope"] = {
        "falsifier": "Empirical graphlets lack low-order dominance with a measurable high-order tail.",
        "empirical_primary_spectrum": spectral,
        "passes_envelope_screen": bool(
            spectral
            and spectral["gft_low_energy_fraction_mean"] >= float(tests_cfg["spectral_envelope"]["min_low_fraction"])
            and spectral["walsh_high_order_fraction_mean"] >= float(tests_cfg["spectral_envelope"]["min_high_order_fraction"])
            and spectral["walsh_high_order_fraction_mean"] <= float(tests_cfg["spectral_envelope"]["max_high_order_fraction"])
        ),
    }

    best_by_model = pd.DataFrame(comparison.get("best_by_model", []))
    contact_distance = math.nan
    if not best_by_model.empty and "contact_residual_ensemble" in set(best_by_model["candidate_model"]):
        contact_distance = float(
            best_by_model.loc[
                best_by_model["candidate_model"] == "contact_residual_ensemble",
                "distance_to_empirical_primary",
            ].iloc[0]
        )
    control_models = ["random_position_residual", "shuffled_contact_residual", "gnk_local_null"]
    control_equivalent = []
    for model in control_models:
        if not best_by_model.empty and model in set(best_by_model["candidate_model"]):
            distance = float(
                best_by_model.loc[
                    best_by_model["candidate_model"] == model,
                    "distance_to_empirical_primary",
                ].iloc[0]
            )
            if np.isfinite(contact_distance) and distance <= contact_distance * float(tests_cfg["control_equivalence_margin"]):
                control_equivalent.append({"candidate_model": model, "distance": distance})
    out["control_equivalence"] = {
        "falsifier": "Matched controls match empirical graphlets as well as contact residuals.",
        "contact_distance": contact_distance,
        "equivalent_or_better_controls": control_equivalent,
        "passes_specificity_screen": bool(len(control_equivalent) == 0 and np.isfinite(contact_distance)),
    }

    readout_row = readout_summary[readout_summary["score_column"] == primary_score]
    retained = math.nan
    if not readout_row.empty:
        retained = float(readout_row["high_order_fraction_retained_after_readout_fit_mean"].iloc[0])
    out["scale_confounding"] = {
        "falsifier": "Observed high-order WHT signal disappears after nonlinear global-readout correction.",
        "mean_high_order_fraction_retained_after_readout_fit": retained,
        "passes_not_scale_artifact_screen": bool(
            np.isfinite(retained)
            and retained >= float(tests_cfg["scale_confounding"]["min_retained_high_order_fraction"])
        ),
    }

    smooth_distance = math.nan
    if not best_by_model.empty and "smooth_core_only" in set(best_by_model["candidate_model"]):
        smooth_distance = float(
            best_by_model.loc[
                best_by_model["candidate_model"] == "smooth_core_only",
                "distance_to_empirical_primary",
            ].iloc[0]
        )
    out["smooth_core_only_empirical_fit"] = {
        "falsifier": "Smooth core or additive global-readout model matches empirical topology and spectra.",
        "smooth_core_distance": smooth_distance,
        "contact_distance": contact_distance,
        "passes_smooth_core_rejection_screen": bool(
            np.isfinite(smooth_distance)
            and np.isfinite(contact_distance)
            and smooth_distance > contact_distance * float(tests_cfg["smooth_core_margin"])
        ),
    }

    topology_metrics = ["accessible_path_fraction", "n_local_maxima", "random_walk_lag1_autocorrelation"]
    spectral_metrics = ["walsh_high_order_fraction", "gft_high_energy_fraction"]
    best_spectral_model = None
    if not best_by_model.empty:
        best_spectral_model = best_by_model.sort_values("distance_to_empirical_primary").iloc[0].to_dict()
    out["wrong_evolutionary_dynamics"] = {
        "falsifier": "Closest spectral models fail topology/accessibility/random-walk dynamics.",
        "comparison_metrics_include_topology": topology_metrics,
        "comparison_metrics_include_spectra": spectral_metrics,
        "best_overall_model": best_spectral_model,
        "passes_joint_metric_screen": bool(best_spectral_model is not None),
    }
    return out


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = graphlet_runner.load_config(config_path, quick=bool(args.quick))
    synthetic_dir = output_dir / "synthetic_graphlet_suite"
    empirical_dir = output_dir / "empirical_graphlet_suite"
    synthetic_summary = graphlet_runner.run_experiment(
        cfg,
        config_path=config_path,
        project_root=project_root,
        output_dir=synthetic_dir,
    )
    empirical_summary = empirical_graphlet_analysis(cfg, project_root=project_root, output_dir=empirical_dir)
    comparison = compare_synthetic_to_empirical(
        cfg,
        synthetic_dir=synthetic_dir,
        empirical_dir=empirical_dir,
        output_dir=output_dir,
    )
    falsifiers = evaluate_falsifiers(
        cfg,
        synthetic_dir=synthetic_dir,
        empirical_dir=empirical_dir,
        comparison=comparison,
    )
    summary = {
        "experiment_id": cfg.get("experiment_id"),
        "synthetic_graphlet_suite": synthetic_summary,
        "empirical_graphlet_suite": empirical_summary,
        "synthetic_empirical_comparison": comparison,
        "falsification_tests": falsifiers,
        "output_layout": {
            "synthetic_graphlet_suite": str(synthetic_dir),
            "empirical_graphlet_suite": str(empirical_dir),
            "synthetic_empirical_model_distances": str(output_dir / "synthetic_empirical_model_distances.csv"),
            "falsification_tests": str(output_dir / "falsification_tests.json"),
        },
    }
    write_json(output_dir / "falsification_tests.json", falsifiers)
    write_json(output_dir / "summary.json", summary)
    console_summary = {
        "experiment_id": cfg.get("experiment_id"),
        "synthetic_graphlet_rows": synthetic_summary.get("n_metric_rows"),
        "empirical_graphlet_rows": empirical_summary.get("n_empirical_metric_rows"),
        "empirical_primary_score": cfg["empirical_graphlets"].get("primary_score_column"),
        "falsification_tests": falsifiers,
        "summary_json": str(output_dir / "summary.json"),
    }
    print(json.dumps(console_summary, indent=2, sort_keys=True, default=graphlet_runner.json_default))


if __name__ == "__main__":
    main()
