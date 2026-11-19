#!/usr/bin/env python3
"""System-specific current-model ProteinGym refits plus graphlet falsifiers.

This runner fixes the RUN-073 pooling error. It fits one current AdaptEnv
landscape per ProteinGym assay, then evaluates empirical and fitted-model dense
Hamming graphlets for the RES-028 breakpoint tests.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import shutil
import subprocess
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

from adaptenv import AlignmentProfile, FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from adaptenv.amino_acids import AA_TO_IDX  # noqa: E402
from empirical_graphlet_falsification_suite import (  # noqa: E402
    find_complete_graphlets,
    fit_global_readout,
    metric_rows_for_values,
    score_map_for_assay,
    subset_masks,
)
import tem1_aggregate_graphlet_residual_ruggedness as graphlet_runner  # noqa: E402
from project_adapt_env.proteingym_panel import merge_dict  # noqa: E402


EXPECTED_ASSAYS = {
    "PHOT_CHLRE_Chen_2023",
    "SPTN1_CHICK_Tsuboyama_2023_1TUD",
}

PROHIBITED_PARAMETER_NAMES = {
    "n_functional_dims",
    "functional_sigma_base",
    "functional_sigma_anisotropy",
    "functional_weight",
    "peak_distance_from_consensus",
    "noise_amplitude",
    "noise_heterogeneity",
    "noise_taper_mode",
    "noise_taper_exponent",
    "strict_peak_lock",
    "peak_lock_margin",
    "n_higher_order_epistatic_terms",
    "epistasis_order_max",
    "epistasis_order_decay",
    "readout_stability_midpoint",
    "readout_stability_slope",
    "readout_function_exponent",
    "readout_function_midpoint",
    "readout_function_slope",
}

LEGACY_BASE_DEFAULTS: dict[str, Any] = {
    "n_functional_dims": 4,
    "functional_sigma_base": 20.0,
    "functional_sigma_anisotropy": 0.5,
    "functional_weight": 1.0,
    "peak_distance_from_consensus": 0,
    "noise_amplitude": 0.0,
    "noise_taper_mode": "uniform",
    "noise_taper_exponent": 1.0,
    "strict_peak_lock": False,
    "n_higher_order_epistatic_terms": 0,
}


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


def read_csv_if_nonempty(path: Path) -> list[dict[str, Any]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    try:
        return pd.read_csv(path).to_dict(orient="records")
    except pd.errors.EmptyDataError:
        return []


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if quick:
        merged = merge_dict(raw, raw.get("quick", {}))
        merged.pop("quick", None)
        return merged
    raw.pop("quick", None)
    return raw


def as_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def assert_current_semantics(config: dict[str, Any]) -> dict[str, Any]:
    assay_ids = {str(system["assay_id"]) for system in config.get("systems", [])}
    if assay_ids != EXPECTED_ASSAYS:
        raise ValueError(
            "EXP-030 must fit exactly PHOT and SPTN1 as separate systems; "
            f"got {sorted(assay_ids)}"
        )
    if int(config.get("panel", {}).get("panel_size", 0)) != 1:
        raise ValueError("EXP-030 per-system child configs must use panel_size=1")

    base = dict(config.get("adapt_env", {}).get("base_config", {}))
    combine_mode = str(base.get("observed_fitness_combine_mode", "")).strip().lower()
    if combine_mode != "stability_margin":
        raise ValueError("Current-model refit must use observed_fitness_combine_mode=stability_margin")
    target = str(base.get("generic_epistasis_target", "")).strip().lower()
    if target != "stability":
        raise ValueError("Current-model epistasis must be routed into stability")

    for key, default in LEGACY_BASE_DEFAULTS.items():
        if key not in base:
            continue
        value = base[key]
        if isinstance(default, bool):
            if bool(value) != default:
                raise ValueError(f"Legacy base_config key {key} is active: {value!r}")
        elif isinstance(default, str):
            if str(value).strip().lower() != default:
                raise ValueError(f"Legacy base_config key {key} is active: {value!r}")
        else:
            numeric = as_float(value)
            if numeric is None or abs(numeric - float(default)) > 1e-12:
                raise ValueError(f"Legacy base_config key {key} is active: {value!r}")

    if base.get("latent_trait_blocks"):
        raise ValueError("Current-model refit must not configure latent_trait_blocks")
    if base.get("observed_fitness_terms"):
        raise ValueError("Current-model refit must not configure observed_fitness_terms")

    fitted_names = {
        str(item.get("name"))
        for item in config.get("abc", {}).get("parameters", [])
        if isinstance(item, dict)
    }
    prohibited_fit = sorted(fitted_names & PROHIBITED_PARAMETER_NAMES)
    if prohibited_fit:
        raise ValueError(f"Legacy parameters are prohibited in EXP-030 ABC fit: {prohibited_fit}")

    calibration = dict(config.get("abc", {}).get("calibration", {}))
    if str(calibration.get("synthetic_readout_mode", "")).strip().lower() != "stability_margin":
        raise ValueError("abc.calibration.synthetic_readout_mode must be stability_margin")
    if str(calibration.get("empirical_pairwise_target", "")).strip().lower() != "stability":
        raise ValueError("abc.calibration.empirical_pairwise_target must be stability")
    if config.get("deterministic_branches"):
        raise ValueError("EXP-030 must use the ABC pathway only; deterministic branches stay disabled")

    return {
        "assay_ids": sorted(assay_ids),
        "fitted_parameter_names": sorted(fitted_names),
        "observed_fitness_combine_mode": combine_mode,
        "generic_epistasis_target": target,
        "fit_mode": "two_independent_system_specific_landscapes",
        "legacy_parameter_guard": "passed",
    }


def child_config_for_system(config: dict[str, Any], system: dict[str, Any]) -> dict[str, Any]:
    child = copy.deepcopy(config)
    assay_id = str(system["assay_id"])
    child["experiment"] = dict(child.get("experiment", {}))
    child["experiment"]["id"] = f"{config['experiment']['id']}::{assay_id}"
    child["panel"] = dict(child["panel"])
    child["panel"]["panel_size"] = 1
    child["panel"]["assay_ids"] = [assay_id]
    child["panel"]["assay_roles"] = {assay_id: system.get("role", "system_specific_fit")}
    child["panel"]["stale_cache_use"] = {assay_id: "comparison_only_not_parameter_source"}
    child["abc"] = copy.deepcopy(child["abc"])
    child["abc"]["calibration"] = dict(child["abc"]["calibration"])
    child["abc"]["calibration"]["max_landscapes_per_fit"] = 1
    child["systems"] = [copy.deepcopy(system)]
    child["graphlet_evaluation"] = copy.deepcopy(config["graphlet_evaluation"])
    child.pop("resume", None)
    return child


def resume_source_output_dir(config: dict[str, Any], output_dir: Path) -> Path | None:
    resume_cfg = dict(config.get("resume") or {})
    if not bool(resume_cfg.get("reuse_completed_system_fits", False)):
        return None
    explicit = str(resume_cfg.get("source_output_dir") or "").strip()
    if explicit:
        path = Path(explicit).expanduser()
        return path.resolve() if path.exists() else path
    source_run_id = str(resume_cfg.get("source_run_id") or "").strip()
    if not source_run_id:
        return None
    current_run_id = str(os.environ.get("LABPROJ_RUN_ID") or output_dir.parent.name).strip()
    if current_run_id == source_run_id:
        return None
    if output_dir.name == "outputs":
        return output_dir.parent.parent / source_run_id / "outputs"
    if output_dir.parent.name == "outputs":
        return output_dir.parent / source_run_id
    return None


def copy_completed_fit_from_resume_source(
    *,
    config: dict[str, Any],
    assay_id: str,
    output_dir: Path,
    fit_dir: Path,
) -> dict[str, Any] | None:
    if (fit_dir / "summary.json").exists():
        return json.loads((fit_dir / "summary.json").read_text(encoding="utf-8"))
    source_output = resume_source_output_dir(config, output_dir)
    if source_output is None:
        return None
    source_fit_dir = source_output / "systems" / assay_id / "fit"
    source_summary = source_fit_dir / "summary.json"
    source_checkpoint = source_fit_dir / "checkpoints" / "smc_empirical.json"
    if not source_summary.exists() or not source_checkpoint.exists():
        return None
    fit_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source_fit_dir, fit_dir, dirs_exist_ok=True)
    print(
        json.dumps(
            {
                "event": "reused_completed_system_fit",
                "assay_id": assay_id,
                "source_fit_dir": str(source_fit_dir),
                "target_fit_dir": str(fit_dir),
            },
            sort_keys=True,
        ),
        file=sys.stderr,
        flush=True,
    )
    return json.loads((fit_dir / "summary.json").read_text(encoding="utf-8"))


def run_single_system_fit(
    *,
    config: dict[str, Any],
    system: dict[str, Any],
    project_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    assay_id = str(system["assay_id"])
    system_dir = output_dir / "systems" / assay_id
    fit_dir = system_dir / "fit"
    config_dir = output_dir / "system_configs"
    fit_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)
    reused = copy_completed_fit_from_resume_source(
        config=config,
        assay_id=assay_id,
        output_dir=output_dir,
        fit_dir=fit_dir,
    )
    if reused is not None:
        return reused
    child_config = child_config_for_system(config, system)
    child_config_path = config_dir / f"{assay_id}.yaml"
    child_config_path.write_text(yaml.safe_dump(child_config, sort_keys=False), encoding="utf-8")
    command = [
        sys.executable,
        "-u",
        str(project_root / "scripts" / "proteingym_bayesian_summary_calibration.py"),
        "--config",
        str(child_config_path),
        "--output-dir",
        str(fit_dir),
        "--project-root",
        str(project_root),
    ]
    subprocess.run(command, check=True)
    return json.loads((fit_dir / "summary.json").read_text(encoding="utf-8"))


def sequence_indices_for_masks(
    *,
    wildtype: str,
    masks: list[frozenset[tuple[int, str]]],
) -> np.ndarray:
    rows: list[list[int]] = []
    for mask in masks:
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
    seed: int,
    smooth_core_only: bool = False,
) -> FitnessLandscapeEnv:
    payload = dict(base_config)
    payload.update(updates)
    payload.update({"L": len(wildtype), "seed": int(seed)})
    if smooth_core_only:
        payload.update(
            {
                "n_epistatic_terms": 0,
                "epistasis_strength": 0.0,
                "epistasis_covariation_bias": 0.0,
                "empirical_pairwise_strength": 0.0,
            }
        )
    config = LandscapeConfig(**payload)
    profile = AlignmentProfile.from_alignment(alignment_path, max_sequences=512)
    return FitnessLandscapeEnv(
        config,
        alignment_profile=profile,
        reference_sequence=wildtype,
    )


def evaluate_landscape_on_graphlet(
    landscape: FitnessLandscapeEnv,
    *,
    wildtype: str,
    features: tuple[tuple[int, str], ...],
) -> np.ndarray:
    masks = subset_masks(features)
    seqs = sequence_indices_for_masks(wildtype=wildtype, masks=masks)
    return np.asarray(landscape.evaluate_batch(seqs), dtype=np.float64)


def gnk_values(
    *,
    dim: int,
    seed: int,
    target_mean: float,
    target_sd: float,
    k: int = 2,
) -> np.ndarray:
    rng = np.random.RandomState(int(seed) % (2**32))
    k = max(0, min(int(k), max(int(dim) - 1, 0)))
    deps: list[list[int]] = []
    tables: list[np.ndarray] = []
    for locus in range(dim):
        dep = [locus] + [int((locus + offset + 1) % dim) for offset in range(k)]
        deps.append(dep)
        tables.append(rng.normal(size=2 ** len(dep)))
    values = np.zeros(2**dim, dtype=np.float64)
    for node in range(2**dim):
        total = 0.0
        for dep, table in zip(deps, tables):
            idx = 0
            for bit_index, bit in enumerate(dep):
                idx |= ((node >> bit) & 1) << bit_index
            total += float(table[idx])
        values[node] = total / max(dim, 1)
    sd = float(np.std(values))
    if sd > 1e-12:
        values = (values - float(np.mean(values))) / sd
    else:
        values = values - float(np.mean(values))
    return float(target_mean) + values * max(float(target_sd), 1e-8)


def finite_mean(rows: pd.DataFrame, metric: str) -> float:
    if metric not in rows:
        return math.nan
    values = pd.to_numeric(rows[metric], errors="coerce").dropna()
    return float(values.mean()) if not values.empty else math.nan


def metric_distance(
    model_rows: pd.DataFrame,
    target_rows: pd.DataFrame,
    metrics: list[str],
) -> tuple[float, dict[str, float]]:
    pieces: list[float] = []
    by_metric: dict[str, float] = {}
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
        pieces.append(float(value))
    return (float(np.mean(pieces)) if pieces else math.nan), by_metric


def graphlet_split(index: int) -> str:
    return "calibration" if int(index) % 2 == 0 else "test"


def graphlet_analysis_for_system(
    *,
    config: dict[str, Any],
    system: dict[str, Any],
    fit_summary: dict[str, Any],
    project_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    assay_id = str(system["assay_id"])
    system_dir = output_dir / "systems" / assay_id
    fit_dir = system_dir / "fit"
    graphlet_cfg = dict(system["graphlets"])
    eval_cfg = dict(config["graphlet_evaluation"])
    metric_cols = [str(item) for item in eval_cfg["metric_columns"]]
    primary_score = str(eval_cfg["primary_score_column"])
    score_columns = [str(item) for item in eval_cfg["score_columns"]]
    random_walk_cfg = {"analysis": config.get("analysis", {})}

    panel = pd.read_csv(fit_dir / "selected_panel.csv")
    wildtype = str(panel.loc[panel["DMS_id"] == assay_id, "target_seq"].iloc[0])
    assay_path = fit_dir / "checkpoints" / "panel_preparation" / assay_id / "prepared_assay.csv"
    metrics_path = fit_dir / "checkpoints" / "panel_preparation" / assay_id / "mavenn_metrics.json"
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
    complete: list[tuple[int, tuple[tuple[int, str], ...]]] = []
    for dim in [int(item) for item in graphlet_cfg["dimensions"]]:
        complete.extend(
            find_complete_graphlets(
                score_map,
                score_columns=assay_score_columns,
                dimensions=[dim],
                max_per_dimension=int(graphlet_cfg["max_graphlets_per_dimension"]),
                top_feature_pool=int(graphlet_cfg["top_feature_pool"]),
            )
        )

    alignment_path = project_root / config["mmseqs"]["alignment_dir"] / f"{assay_id}.fasta"
    if not alignment_path.exists() and fit_summary.get("panel"):
        candidate = Path(str(fit_summary["panel"][0].get("alignment_path", "")))
        if candidate.exists():
            alignment_path = candidate
    if not alignment_path.exists():
        raise FileNotFoundError(f"Alignment FASTA not found for {assay_id}: {alignment_path}")
    shutil.copy2(alignment_path, system_dir / f"{assay_id}.alignment.fasta")

    branches = fit_summary["branches"]
    best_updates = dict(branches["smc_abc_best_raw"]["fitted_parameters"])
    mean_updates = dict(branches["smc_abc_posterior_mean_raw"]["fitted_parameters"])
    seed = int(config["adapt_env"]["base_config"].get("seed", 42))
    base_config = dict(config["adapt_env"]["base_config"])
    fitted_best = build_landscape(
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

    catalog_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    readout_rows: list[dict[str, Any]] = []

    for graphlet_index, (dim, features) in enumerate(complete):
        split = graphlet_split(graphlet_index)
        base_graphlet = {
            "graphlet_id": f"{assay_id}__g{graphlet_index:04d}__d{dim}",
            "assay_id": assay_id,
            "graphlet_split": split,
            "graphlet_dimension": int(dim),
            "mutation_features": ";".join(f"{pos}:{aa}" for pos, aa in features),
        }
        catalog_rows.append(base_graphlet)
        masks = subset_masks(features)
        best_values = evaluate_landscape_on_graphlet(
            fitted_best,
            wildtype=wildtype,
            features=features,
        )
        mean_values = evaluate_landscape_on_graphlet(
            fitted_mean,
            wildtype=wildtype,
            features=features,
        )
        smooth_values = evaluate_landscape_on_graphlet(
            smooth_core,
            wildtype=wildtype,
            features=features,
        )
        residual_values = best_values - smooth_values

        for score_col in assay_score_columns:
            empirical_values = np.asarray([score_map[mask][score_col] for mask in masks], dtype=np.float64)
            graphlet = {**base_graphlet, "split": split, "score_column": score_col}
            candidates = [
                (f"empirical_observed__{score_col}", score_col, empirical_values),
                ("fitted_current_best", f"{assay_id}__best", best_values),
                ("fitted_current_posterior_mean", f"{assay_id}__posterior_mean", mean_values),
                ("fitted_smooth_core_only", f"{assay_id}__smooth_core", smooth_values),
                ("fitted_residual_only_diagnostic", f"{assay_id}__residual_only", residual_values),
                (
                    "gnk_local_null",
                    f"{assay_id}__gnk_k2",
                    gnk_values(
                        dim=int(dim),
                        seed=graphlet_runner.stable_seed("gnk", assay_id, graphlet_index, score_col),
                        target_mean=float(np.mean(empirical_values)),
                        target_sd=float(np.std(empirical_values)),
                        k=2,
                    ),
                ),
            ]
            pred, resid, r2 = fit_global_readout(
                empirical_values,
                dim=int(dim),
                degree=int(eval_cfg.get("global_readout_polynomial_degree", 3)),
            )
            candidates.extend(
                [
                    (f"empirical_additive_global_readout_fit__{score_col}", f"{score_col}__poly", pred),
                    (f"empirical_global_readout_residual__{score_col}", f"{score_col}__residual", resid),
                ]
            )
            observed_high = graphlet_runner.walsh_decomposition(empirical_values, dim=int(dim))[1]
            residual_high = graphlet_runner.walsh_decomposition(resid, dim=int(dim))[1]
            observed_high_fraction = float(
                sum(value for order, value in observed_high.items() if int(order) >= 3)
            )
            residual_high_fraction = float(
                sum(value for order, value in residual_high.items() if int(order) >= 3)
            )
            readout_rows.append(
                {
                    **graphlet,
                    "global_readout_r2": float(r2),
                    "observed_high_order_fraction": observed_high_fraction,
                    "post_global_readout_residual_high_order_fraction": residual_high_fraction,
                    "high_order_fraction_retained_after_readout_fit": (
                        residual_high_fraction / observed_high_fraction
                        if observed_high_fraction > 1e-12
                        else math.nan
                    ),
                }
            )
            for candidate_model, parameter_set_id, values in candidates:
                metric, walsh, gft, _access = metric_rows_for_values(
                    graphlet=graphlet,
                    values=np.asarray(values, dtype=np.float64),
                    candidate_model=candidate_model,
                    parameter_set_id=parameter_set_id,
                    cfg=random_walk_cfg,
                )
                metric["fit_mode"] = "system_specific"
                metric_rows.append(metric)
                for row in walsh:
                    row["graphlet_split"] = split
                    row["fit_mode"] = "system_specific"
                for row in gft:
                    row["graphlet_split"] = split
                    row["fit_mode"] = "system_specific"
                walsh_rows.extend(walsh)
                gft_rows.extend(gft)

    metric_df = pd.DataFrame(metric_rows)
    summary_rows: list[dict[str, Any]] = []
    if not metric_df.empty:
        for keys, sub in metric_df.groupby(
            ["assay_id", "graphlet_split", "candidate_model", "parameter_set_id", "score_column"],
            dropna=False,
        ):
            row = {
                "assay_id": keys[0],
                "graphlet_split": keys[1],
                "candidate_model": keys[2],
                "parameter_set_id": keys[3],
                "score_column": keys[4],
                "n_graphlets": int(len(sub)),
            }
            for col in metric_cols:
                if col in sub:
                    values = pd.to_numeric(sub[col], errors="coerce").dropna()
                    if not values.empty:
                        row[f"{col}_mean"] = float(values.mean())
                        row[f"{col}_sd"] = float(values.std(ddof=0))
            summary_rows.append(row)

    distances: list[dict[str, Any]] = []
    for split in ["calibration", "test", "all"]:
        if split == "all":
            split_df = metric_df
        else:
            split_df = metric_df[metric_df["graphlet_split"] == split]
        if split_df.empty:
            continue
        target_rows = split_df[
            (split_df["candidate_model"] == f"empirical_observed__{primary_score}")
            & (split_df["score_column"] == primary_score)
        ]
        if target_rows.empty:
            continue
        for candidate_model in [
            "fitted_current_best",
            "fitted_current_posterior_mean",
            "fitted_smooth_core_only",
            "fitted_residual_only_diagnostic",
            "gnk_local_null",
        ]:
            model_rows = split_df[
                (split_df["candidate_model"] == candidate_model)
                & (split_df["score_column"] == primary_score)
            ]
            if model_rows.empty:
                continue
            distance, by_metric = metric_distance(model_rows, target_rows, metric_cols)
            distances.append(
                {
                    "assay_id": assay_id,
                    "graphlet_split": split,
                    "candidate_model": candidate_model,
                    "distance_to_empirical_primary": distance,
                    **{f"distance_component__{key}": value for key, value in by_metric.items()},
                }
            )

    write_csv(system_dir / "graphlet_catalog.csv", catalog_rows)
    write_csv(system_dir / "graphlet_metrics.csv", metric_rows)
    write_csv(system_dir / "walsh_terms.csv", walsh_rows)
    write_csv(system_dir / "gft_coefficients.csv", gft_rows)
    write_csv(system_dir / "global_readout_confounding.csv", readout_rows)
    write_csv(system_dir / "graphlet_summary.csv", summary_rows)
    write_csv(system_dir / "model_distances.csv", distances)

    return {
        "assay_id": assay_id,
        "n_complete_graphlets": int(len(catalog_rows)),
        "graphlet_dimensions_observed": sorted({int(row["graphlet_dimension"]) for row in catalog_rows}),
        "n_metric_rows": int(len(metric_rows)),
        "n_walsh_rows": int(len(walsh_rows)),
        "n_gft_rows": int(len(gft_rows)),
        "catalog_rows": catalog_rows,
        "summary_rows": summary_rows,
        "distance_rows": distances,
        "readout_rows": readout_rows,
        "fit_summary": {
            "best_distance": fit_summary["smc_abc"]["best_particle"]["distance"],
            "round_summaries": fit_summary["smc_abc"]["round_summaries"],
            "branch_validations": fit_summary["branches"],
        },
    }


def mean_from_rows(rows: list[dict[str, Any]], key: str) -> float:
    values = [as_float(row.get(key)) for row in rows]
    values = [value for value in values if value is not None]
    return float(np.mean(values)) if values else math.nan


def evaluate_breakpoints(config: dict[str, Any], system_summaries: list[dict[str, Any]]) -> dict[str, Any]:
    eval_cfg = config["graphlet_evaluation"]
    primary_score = str(eval_cfg["primary_score_column"])
    spectral_cfg = eval_cfg["spectral_envelope"]
    tests: dict[str, Any] = {}
    for system in system_summaries:
        assay_id = system["assay_id"]
        metrics = pd.DataFrame(system["summary_rows"])
        distances = pd.DataFrame(system["distance_rows"])
        readout = pd.DataFrame(system["readout_rows"])
        test_metrics = metrics[
            (metrics["graphlet_split"] == "test")
            & (metrics["candidate_model"] == f"empirical_observed__{primary_score}")
            & (metrics["score_column"] == primary_score)
        ]
        if test_metrics.empty:
            test_metrics = metrics[
                (metrics["candidate_model"] == f"empirical_observed__{primary_score}")
                & (metrics["score_column"] == primary_score)
            ]
        spectral = {}
        if not test_metrics.empty:
            row = test_metrics.iloc[0]
            spectral = {
                "gft_low_energy_fraction_mean": as_float(row.get("gft_low_energy_fraction_mean")),
                "gft_high_energy_fraction_mean": as_float(row.get("gft_high_energy_fraction_mean")),
                "walsh_high_order_fraction_mean": as_float(row.get("walsh_high_order_fraction_mean")),
            }
        max_dim = max(system["graphlet_dimensions_observed"], default=0)
        spectral_testable = bool(max_dim >= 3 and spectral)
        passes_spectral = bool(
            spectral_testable
            and spectral["gft_low_energy_fraction_mean"] is not None
            and spectral["walsh_high_order_fraction_mean"] is not None
            and spectral["gft_low_energy_fraction_mean"] >= float(spectral_cfg["min_low_fraction"])
            and spectral["walsh_high_order_fraction_mean"] >= float(spectral_cfg["min_high_order_fraction"])
            and spectral["walsh_high_order_fraction_mean"] <= float(spectral_cfg["max_high_order_fraction"])
        )

        distance_split = distances[distances["graphlet_split"] == "test"]
        if distance_split.empty:
            distance_split = distances[distances["graphlet_split"] == "all"]
        best_row = None
        fitted_distance = math.nan
        if not distance_split.empty:
            finite = distance_split[np.isfinite(pd.to_numeric(distance_split["distance_to_empirical_primary"], errors="coerce"))]
            if not finite.empty:
                best_row = finite.sort_values("distance_to_empirical_primary").iloc[0].to_dict()
            fitted = distance_split[distance_split["candidate_model"] == "fitted_current_best"]
            if not fitted.empty:
                fitted_distance = float(fitted["distance_to_empirical_primary"].iloc[0])
        best_distance = float(best_row["distance_to_empirical_primary"]) if best_row else math.nan
        wrong_dynamics_fires = bool(
            not np.isfinite(fitted_distance)
            or (
                np.isfinite(best_distance)
                and fitted_distance > best_distance * float(eval_cfg["control_equivalence_margin"])
            )
        )

        smooth_distance = math.nan
        residual_distance = math.nan
        gnk_distance = math.nan
        for model_name, target_name in [
            ("fitted_smooth_core_only", "smooth_distance"),
            ("fitted_residual_only_diagnostic", "residual_only_distance"),
            ("gnk_local_null", "gnk_distance"),
        ]:
            rows = distance_split[distance_split["candidate_model"] == model_name]
            value = float(rows["distance_to_empirical_primary"].iloc[0]) if not rows.empty else math.nan
            if target_name == "smooth_distance":
                smooth_distance = value
            elif target_name == "residual_only_distance":
                residual_distance = value
            else:
                gnk_distance = value

        retained = math.nan
        if not readout.empty:
            sub = readout[readout["score_column"] == primary_score]
            retained = mean_from_rows(sub.to_dict(orient="records"), "high_order_fraction_retained_after_readout_fit")

        tests[assay_id] = {
            "wrong_spectral_envelope": {
                "falsifier": "Empirical graphlets do not show low-order dominance with a measurable but non-dominant high-order WHT/GFT tail.",
                "testable": spectral_testable,
                "max_graphlet_dimension": int(max_dim),
                "empirical_primary_spectrum": spectral,
                "passes_expected_envelope": passes_spectral,
                "fires": bool(spectral_testable and not passes_spectral),
            },
            "wrong_evolutionary_dynamics": {
                "falsifier": "The fitted embedded current model fails the joint graphlet package: spectra, SE/RSE, local maxima, accessibility and random-walk autocorrelation.",
                "distance_split_used": str(distance_split["graphlet_split"].iloc[0]) if not distance_split.empty else None,
                "best_overall_model": best_row,
                "fitted_current_best_distance": fitted_distance,
                "smooth_core_distance": smooth_distance,
                "residual_only_distance": residual_distance,
                "gnk_distance": gnk_distance,
                "fires": wrong_dynamics_fires,
            },
            "scale_confounding_control": {
                "falsifier": "High-order WHT signal disappears after nonlinear global-readout correction.",
                "mean_high_order_fraction_retained_after_readout_fit": retained,
                "passes_not_scale_artifact_screen": bool(
                    np.isfinite(retained)
                    and retained >= float(eval_cfg["scale_confounding"]["min_retained_high_order_fraction"])
                ),
                "fires": bool(
                    np.isfinite(retained)
                    and retained < float(eval_cfg["scale_confounding"]["min_retained_high_order_fraction"])
                ),
            },
            "smooth_core_control": {
                "falsifier": "Smooth-core-only fitted model explains the empirical graphlet statistics as well as the residual current model.",
                "smooth_core_distance": smooth_distance,
                "fitted_current_best_distance": fitted_distance,
                "passes_smooth_core_rejection_screen": bool(
                    np.isfinite(smooth_distance)
                    and np.isfinite(fitted_distance)
                    and smooth_distance > fitted_distance * float(eval_cfg["smooth_core_margin"])
                ),
            },
        }
    return tests


def collect_fit_rows(system_summaries: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    fit_rows: list[dict[str, Any]] = []
    branch_rows: list[dict[str, Any]] = []
    for system in system_summaries:
        assay_id = system["assay_id"]
        fit = system["fit_summary"]
        rounds = fit["round_summaries"]
        fit_rows.append(
            {
                "assay_id": assay_id,
                "best_distance": fit["best_distance"],
                "round_count": len(rounds),
                "round0_best_distance": rounds[0]["best_distance"] if rounds else math.nan,
                "final_round_best_distance": rounds[-1]["best_distance"] if rounds else math.nan,
                "min_round_best_distance": min(float(row["best_distance"]) for row in rounds) if rounds else math.nan,
            }
        )
        for branch, payload in fit["branch_validations"].items():
            row = {
                "assay_id": assay_id,
                "branch": branch,
                "fit_mode": "system_specific",
                "source": payload.get("source", "raw"),
            }
            for prefix, mapping in [
                ("fitted", payload.get("fitted_parameters", {})),
                ("validation", payload.get("validation", {})),
                ("objective", payload.get("objective_terms", {})),
            ]:
                for key, value in dict(mapping).items():
                    row[f"{prefix}__{key}"] = value
            branch_rows.append(row)
    return fit_rows, branch_rows


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(config_path, quick=bool(args.quick))
    preflight = assert_current_semantics(config)
    preflight["sampling_budget"] = config["sampling_budget"]
    write_json(output_dir / "current_model_preflight.json", preflight)

    system_summaries: list[dict[str, Any]] = []
    for system in config["systems"]:
        fit_summary = run_single_system_fit(
            config=config,
            system=system,
            project_root=project_root,
            output_dir=output_dir,
        )
        graphlet_summary = graphlet_analysis_for_system(
            config=config,
            system=system,
            fit_summary=fit_summary,
            project_root=project_root,
            output_dir=output_dir,
        )
        system_summaries.append(graphlet_summary)

    fit_rows, branch_rows = collect_fit_rows(system_summaries)
    write_csv(output_dir / "system_fit_summary.csv", fit_rows)
    write_csv(output_dir / "system_branch_validations.csv", branch_rows)

    all_metric_rows: list[dict[str, Any]] = []
    all_walsh_rows: list[dict[str, Any]] = []
    all_gft_rows: list[dict[str, Any]] = []
    all_distance_rows: list[dict[str, Any]] = []
    all_catalog_rows: list[dict[str, Any]] = []
    for system in system_summaries:
        assay_id = system["assay_id"]
        system_dir = output_dir / "systems" / assay_id
        all_metric_rows.extend(read_csv_if_nonempty(system_dir / "graphlet_metrics.csv"))
        all_walsh_rows.extend(read_csv_if_nonempty(system_dir / "walsh_terms.csv"))
        all_gft_rows.extend(read_csv_if_nonempty(system_dir / "gft_coefficients.csv"))
        all_distance_rows.extend(read_csv_if_nonempty(system_dir / "model_distances.csv"))
        all_catalog_rows.extend(read_csv_if_nonempty(system_dir / "graphlet_catalog.csv"))
    write_csv(output_dir / "system_graphlet_catalog.csv", all_catalog_rows)
    write_csv(output_dir / "system_graphlet_metrics.csv", all_metric_rows)
    write_csv(output_dir / "system_walsh_terms.csv", all_walsh_rows)
    write_csv(output_dir / "system_gft_coefficients.csv", all_gft_rows)
    write_csv(output_dir / "system_model_distances.csv", all_distance_rows)

    falsifiers = evaluate_breakpoints(config, system_summaries)
    write_json(output_dir / "falsification_tests.json", falsifiers)
    summary = {
        "experiment_id": config["experiment"]["id"],
        "preflight": preflight,
        "systems": [
            {
                "assay_id": system["assay_id"],
                "n_complete_graphlets": system["n_complete_graphlets"],
                "graphlet_dimensions_observed": system["graphlet_dimensions_observed"],
                "n_metric_rows": system["n_metric_rows"],
                "n_walsh_rows": system["n_walsh_rows"],
                "n_gft_rows": system["n_gft_rows"],
                "fit_summary": system["fit_summary"],
            }
            for system in system_summaries
        ],
        "falsification_tests": falsifiers,
        "output_files": {
            "system_fit_summary": str(output_dir / "system_fit_summary.csv"),
            "system_branch_validations": str(output_dir / "system_branch_validations.csv"),
            "system_graphlet_metrics": str(output_dir / "system_graphlet_metrics.csv"),
            "system_walsh_terms": str(output_dir / "system_walsh_terms.csv"),
            "system_gft_coefficients": str(output_dir / "system_gft_coefficients.csv"),
            "system_model_distances": str(output_dir / "system_model_distances.csv"),
            "falsification_tests": str(output_dir / "falsification_tests.json"),
        },
    }
    write_json(output_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True, default=graphlet_runner.json_default))


if __name__ == "__main__":
    main()
