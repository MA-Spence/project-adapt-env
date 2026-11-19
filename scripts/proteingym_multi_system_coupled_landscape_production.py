#!/usr/bin/env python3
"""Production ProteinGym coupled-landscape calibration panel.

This runner intentionally generalizes the PHOT-only graphlet/coupled-landscape
diagnostics to a fixed, preregistered multi-system ProteinGym panel.  It keeps
the base thermodynamic residual-epistasis fit on the same project SMC-ABC path
used by earlier ProteinGym work, then fits lightweight alignment-routed two-way
and three-way auxiliary ruggedness layers against calibration graphlets.

The design goal is descriptive calibration of bulk landscape behaviour, not
residue-level interaction recovery.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
import os
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for _path in (SRC_ROOT, ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from adaptenv import AlignmentProfile  # noqa: E402
from adaptenv.amino_acids import AA_TO_IDX  # noqa: E402
from empirical_graphlet_falsification_suite import (  # noqa: E402
    fit_global_readout,
    metric_rows_for_values,
    mutation_set,
    score_map_for_assay,
    subset_masks,
)
import tem1_aggregate_graphlet_residual_ruggedness as graphlet_runner  # noqa: E402
from phot_coupled_family_discriminating_abc import (  # noqa: E402
    local_ruggedness_metrics_for_values,
)
from project_adapt_env.proteingym_panel import (  # noqa: E402
    canonical_sequence,
    download_file,
    merge_dict,
    normalize_includes_multiple,
)
from project_adapt_env.smc_abc import (  # noqa: E402
    ParameterSpec,
    SMCABCBackendConfig,
    SMCABCCheckpointConfig,
    SMCABCConfig,
    SimulationResult,
    run_smc_abc,
)
from project_adapt_env.utils import hash_mod, to_builtin  # noqa: E402
from proteingym_system_specific_current_graphlet_refit import (  # noqa: E402
    build_landscape,
    evaluate_landscape_on_graphlet,
)


PRIMARY_SCORE_DEFAULT = "latent_phi"


@dataclass(frozen=True)
class VariantSample:
    frame: pd.DataFrame
    base_values: np.ndarray
    feature_sets: list[tuple[tuple[int, str], ...]]


@dataclass(frozen=True)
class AlignmentCoefficients:
    pair_coefficients: dict[tuple[int, int], float]
    triplet_coefficients: dict[tuple[int, int, int], float]
    top_pairs: list[dict[str, float]]


@dataclass(frozen=True)
class CoupledStandardization:
    base_center: float
    base_sd: float
    raw2_center: float
    raw2_sd: float
    raw3_center: float
    raw3_sd: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    tmp_path.write_text(
        json.dumps(to_builtin(payload), indent=2, sort_keys=True, default=graphlet_runner.json_default)
        + "\n",
        encoding="utf-8",
    )
    tmp_path.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    pd.DataFrame(rows).to_csv(tmp_path, index=False)
    tmp_path.replace(path)


def read_csv_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    try:
        return pd.read_csv(path).to_dict(orient="records")
    except pd.errors.EmptyDataError:
        return []


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if quick:
        raw = merge_dict(raw, raw.get("quick", {}))
    raw.pop("quick", None)
    return raw


def emit_progress(output_dir: Path, payload: dict[str, Any]) -> None:
    progress_path = output_dir / "progress.json"
    log_path = output_dir / "progress.log"
    payload = {
        "run_id": os.environ.get("LABPROJ_RUN_ID"),
        "experiment_id": os.environ.get("LABPROJ_EXPERIMENT_ID"),
        **payload,
    }
    write_json(progress_path, payload)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(to_builtin(payload), sort_keys=True, default=graphlet_runner.json_default))
        handle.write("\n")
    print(json.dumps(to_builtin(payload), sort_keys=True, default=graphlet_runner.json_default), flush=True)


def seed_resume_outputs(output_dir: Path, config: dict[str, Any]) -> None:
    resume_cfg = dict(config.get("resume", {}))
    source_run_id = str(resume_cfg.get("source_run_id") or "").strip()
    if not source_run_id:
        return
    source_output_dir_raw = str(resume_cfg.get("source_output_dir") or "").strip()
    if source_output_dir_raw:
        source_output_dir = Path(source_output_dir_raw)
    else:
        source_output_dir = output_dir.parent.parent / source_run_id / "outputs"
    if not source_output_dir.exists():
        emit_progress(
            output_dir,
            {
                "event": "resume_seed_missing",
                "source_run_id": source_run_id,
                "source_output_dir": str(source_output_dir),
            },
        )
        return
    copied: list[str] = []
    for dirname in ["systems", "system_configs"]:
        source = source_output_dir / dirname
        destination = output_dir / dirname
        if source.exists() and not destination.exists():
            shutil.copytree(source, destination)
            copied.append(dirname)
    for filename in ["selected_systems.csv", "target_selection_audit.json"]:
        source = source_output_dir / filename
        destination = output_dir / filename
        if source.exists() and not destination.exists():
            shutil.copy2(source, destination)
            copied.append(filename)
    emit_progress(
        output_dir,
        {
            "event": "resume_seed_complete",
            "source_run_id": source_run_id,
            "source_output_dir": str(source_output_dir),
            "copied": copied,
        },
    )


def finite_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def rank_spearman(observed: np.ndarray, predicted: np.ndarray) -> float:
    obs = np.asarray(observed, dtype=np.float64)
    pred = np.asarray(predicted, dtype=np.float64)
    mask = np.isfinite(obs) & np.isfinite(pred)
    obs = obs[mask]
    pred = pred[mask]
    if obs.size < 3 or float(np.std(obs)) <= 1e-12 or float(np.std(pred)) <= 1e-12:
        return math.nan
    obs_rank = pd.Series(obs).rank(method="average").to_numpy(dtype=np.float64)
    pred_rank = pd.Series(pred).rank(method="average").to_numpy(dtype=np.float64)
    return float(np.corrcoef(obs_rank, pred_rank)[0, 1])


def nrmse(observed: np.ndarray, predicted: np.ndarray) -> float:
    obs = np.asarray(observed, dtype=np.float64)
    pred = np.asarray(predicted, dtype=np.float64)
    mask = np.isfinite(obs) & np.isfinite(pred)
    obs = obs[mask]
    pred = pred[mask]
    if obs.size == 0:
        return math.nan
    scale = max(float(np.std(obs)), 1e-8)
    return float(np.sqrt(np.mean((pred - obs) ** 2)) / scale)


def summarize_numeric(values: Iterable[float], prefix: str) -> dict[str, float]:
    arr = np.asarray([float(v) for v in values if np.isfinite(float(v))], dtype=np.float64)
    if arr.size == 0:
        return {
            f"{prefix}_n": 0.0,
            f"{prefix}_mean": math.nan,
            f"{prefix}_sd": math.nan,
            f"{prefix}_q10": math.nan,
            f"{prefix}_q50": math.nan,
            f"{prefix}_q90": math.nan,
        }
    return {
        f"{prefix}_n": float(arr.size),
        f"{prefix}_mean": float(np.mean(arr)),
        f"{prefix}_sd": float(np.std(arr)),
        f"{prefix}_q10": float(np.quantile(arr, 0.10)),
        f"{prefix}_q50": float(np.quantile(arr, 0.50)),
        f"{prefix}_q90": float(np.quantile(arr, 0.90)),
    }


def actual_mutation_counts(reference_df: pd.DataFrame, assay_df: pd.DataFrame) -> pd.DataFrame:
    wt_by_id = {
        str(row.DMS_id): canonical_sequence(str(row.target_seq))
        for row in reference_df.itertuples(index=False)
    }
    rows: list[dict[str, Any]] = []
    for dms_id, sub in assay_df.groupby("DMS_id"):
        wildtype = wt_by_id[str(dms_id)]
        counts = []
        for sequence in sub["mutated_sequence"].astype(str):
            if len(sequence) != len(wildtype):
                counts.append(-1)
            else:
                counts.append(sum(a != b for a, b in zip(sequence, wildtype)))
        series = pd.Series(counts, dtype=np.int64)
        rows.append(
            {
                "DMS_id": str(dms_id),
                "n_loaded": int(len(series)),
                "n_single_actual": int((series == 1).sum()),
                "n_double_actual": int((series == 2).sum()),
                "n_triple_actual": int((series == 3).sum()),
                "n_fourplus_actual": int((series >= 4).sum()),
                "max_mutation_count_actual": int(series.max()) if len(series) else 0,
            }
        )
    return pd.DataFrame(rows)


def selection_preflight(
    *,
    project_root: Path,
    config: dict[str, Any],
    output_dir: Path,
) -> pd.DataFrame:
    proteingym_cfg = dict(config["proteingym"])
    reference_path = project_root / proteingym_cfg["reference_csv_path"]
    substitutions_path = project_root / proteingym_cfg["substitutions_parquet_path"]
    download_file(str(proteingym_cfg["reference_csv_url"]), reference_path)
    download_file(str(proteingym_cfg["substitutions_parquet_url"]), substitutions_path)
    reference = pd.read_csv(reference_path)
    for col in [
        "seq_len",
        "DMS_number_single_mutants",
        "DMS_number_multiple_mutants",
        "DMS_total_number_mutants",
        "year",
    ]:
        reference[col] = pd.to_numeric(reference[col], errors="coerce")

    configured = pd.DataFrame(config["systems"])
    assay_ids = [str(item) for item in configured["assay_id"]]
    if len(assay_ids) != len(set(assay_ids)):
        raise ValueError("Duplicate assay_id values in systems config")
    if len(assay_ids) != int(config["panel_design"]["expected_system_count"]):
        raise ValueError(
            f"Expected {config['panel_design']['expected_system_count']} systems, got {len(assay_ids)}"
        )

    selected = reference[reference["DMS_id"].isin(assay_ids)].copy()
    missing = sorted(set(assay_ids) - set(selected["DMS_id"].astype(str)))
    if missing:
        raise ValueError(f"Configured assay IDs not found in ProteinGym reference: {missing}")
    selected["assay_id"] = selected["DMS_id"].astype(str)
    selected = configured.merge(selected, on="assay_id", how="left", suffixes=("", "_reference"))

    if bool(config["panel_design"].get("require_multiple_mutants", True)):
        bad = selected[
            ~selected["includes_multiple_mutants"].map(normalize_includes_multiple)
            | (selected["DMS_number_multiple_mutants"].fillna(0).astype(int) <= 0)
        ]
        if not bad.empty:
            raise ValueError(f"Configured systems without multiple-mutant coverage: {bad['assay_id'].tolist()}")

    max_tsub = int(config["panel_design"]["max_tsuboyama_systems"])
    n_tsub = int(selected["assay_id"].str.contains("Tsuboyama").sum())
    if n_tsub > max_tsub:
        raise ValueError(f"Too many Tsuboyama systems: {n_tsub} > {max_tsub}")
    required = set(str(item) for item in config["panel_design"].get("required_assays", []))
    absent = sorted(required - set(assay_ids))
    if absent:
        raise ValueError(f"Required assays absent from selected panel: {absent}")

    assay_df = pd.read_parquet(
        substitutions_path,
        columns=["DMS_id", "mutated_sequence"],
        filters=[("DMS_id", "in", assay_ids)],
    )
    reference_for_counts = selected[["assay_id", "target_seq"]].rename(
        columns={"assay_id": "DMS_id"}
    )
    actual = actual_mutation_counts(reference_for_counts, assay_df)
    selected = selected.merge(actual, left_on="assay_id", right_on="DMS_id", how="left", suffixes=("", "_actual"))
    selected = selected.drop(columns=[col for col in ["DMS_id_actual"] if col in selected])
    selected["panel_order"] = np.arange(len(selected), dtype=int)
    selected["is_tsuboyama"] = selected["assay_id"].str.contains("Tsuboyama")
    selected["selection_source"] = "ProteinGym reference DMS_substitutions.csv plus local parquet mutation counts"

    write_csv(output_dir / "selected_systems.csv", selected.to_dict(orient="records"))
    write_json(
        output_dir / "target_selection_audit.json",
        {
            "protein_gym_reference_csv": str(reference_path),
            "protein_gym_substitutions_parquet": str(substitutions_path),
            "n_systems": int(len(selected)),
            "n_tsuboyama": n_tsub,
            "n_non_tsuboyama": int(len(selected) - n_tsub),
            "by_coarse_selection_type": selected["coarse_selection_type"].value_counts(dropna=False).to_dict(),
            "by_taxon": selected["taxon"].value_counts(dropna=False).to_dict(),
            "selection_criteria": config["panel_design"],
        },
    )
    return selected


def base_child_config(config: dict[str, Any], assay_id: str) -> dict[str, Any]:
    child = copy.deepcopy(config)
    child["experiment"] = dict(child.get("experiment", {}))
    child["experiment"]["id"] = f"{config['experiment']['id']}::{assay_id}::base"
    child["panel"] = dict(child["panel"])
    child["panel"]["panel_size"] = 1
    child["panel"]["assay_ids"] = [assay_id]
    child["panel"]["allowed_taxa"] = []
    child["panel"]["assays_per_taxon"] = 0
    child["deterministic_branches"] = {}
    child["abc"] = copy.deepcopy(config["base_abc"])
    child["synthetic_truth_recovery"] = {"enabled": False}
    child["resume"] = {}
    for key in [
        "systems",
        "panel_design",
        "graphlets",
        "variant_validation",
        "coupled_fit",
        "walks",
        "outputs",
        "base_abc",
    ]:
        child.pop(key, None)
    return child


def run_base_fit(
    *,
    config: dict[str, Any],
    assay_id: str,
    project_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    system_dir = output_dir / "systems" / assay_id
    fit_dir = system_dir / "base_fit"
    summary_path = fit_dir / "summary.json"
    if summary_path.is_file():
        return json.loads(summary_path.read_text(encoding="utf-8"))

    fit_dir.mkdir(parents=True, exist_ok=True)
    config_dir = output_dir / "system_configs"
    config_dir.mkdir(parents=True, exist_ok=True)
    child = base_child_config(config, assay_id)
    child_config_path = config_dir / f"{assay_id}.base.yaml"
    child_config_path.write_text(yaml.safe_dump(child, sort_keys=False), encoding="utf-8")
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
    return json.loads(summary_path.read_text(encoding="utf-8"))


def sequence_indices_for_strings(sequences: list[str]) -> np.ndarray:
    return np.asarray([[AA_TO_IDX[str(aa)] for aa in seq] for seq in sequences], dtype=np.int64)


def evaluate_landscape_sequences(landscape: Any, sequences: list[str], *, batch_size: int) -> np.ndarray:
    values: list[np.ndarray] = []
    for start in range(0, len(sequences), max(int(batch_size), 1)):
        batch = sequences[start : start + max(int(batch_size), 1)]
        values.append(np.asarray(landscape.evaluate_batch(sequence_indices_for_strings(batch)), dtype=np.float64))
    return np.concatenate(values) if values else np.empty(0, dtype=np.float64)


def variant_feature_sets(wildtype: str, sequences: Iterable[str]) -> list[tuple[tuple[int, str], ...]]:
    out: list[tuple[tuple[int, str], ...]] = []
    for sequence in sequences:
        out.append(tuple(sorted(mutation_set(wildtype, str(sequence)), key=lambda item: (int(item[0]), str(item[1])))))
    return out


def load_prepared_assay(fit_dir: Path, assay_id: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    panel = pd.read_csv(fit_dir / "selected_panel.csv")
    assay_path = fit_dir / "checkpoints" / "panel_preparation" / assay_id / "prepared_assay.csv"
    metrics_path = fit_dir / "checkpoints" / "panel_preparation" / assay_id / "mavenn_metrics.json"
    assay_df = pd.read_csv(assay_path)
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    return panel, assay_df, metrics


def bounded_complete_graphlets(
    score_map: dict[frozenset[tuple[int, str]], dict[str, float]],
    *,
    score_columns: list[str],
    dimensions: list[int],
    max_per_dimension: int,
    top_feature_pool: int,
    max_combinations_per_dimension: int,
) -> list[tuple[int, tuple[tuple[int, str], ...], int]]:
    counts: Counter = Counter()
    for feature_set in score_map:
        for feature in feature_set:
            counts[feature] += 1
    top_features = [feature for feature, _count in counts.most_common(int(top_feature_pool))]
    found: list[tuple[int, tuple[tuple[int, str], ...], int]] = []
    used: set[tuple[tuple[int, str], ...]] = set()
    for dim in sorted(set(int(item) for item in dimensions)):
        checked = 0
        dim_found = 0
        pool = top_features[: max(int(top_feature_pool), dim)]
        for combo in itertools.combinations(pool, dim):
            checked += 1
            if checked > int(max_combinations_per_dimension):
                break
            combo = tuple(combo)
            if combo in used:
                continue
            ok = True
            for mask in subset_masks(combo):
                values = score_map.get(mask)
                if values is None or not all(col in values for col in score_columns):
                    ok = False
                    break
            if not ok:
                continue
            found.append((dim, combo, checked))
            used.add(combo)
            dim_found += 1
            if dim_found >= int(max_per_dimension):
                break
    return found


def build_graphlet_catalog(
    *,
    assay_id: str,
    wildtype: str,
    assay_df: pd.DataFrame,
    mavenn_metrics: dict[str, Any],
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[frozenset[tuple[int, str]], dict[str, float]]]:
    graphlet_cfg = dict(config["graphlets"])
    score_columns = [str(item) for item in graphlet_cfg["score_columns"]]
    wt_scores = {
        "latent_phi": float(mavenn_metrics.get("wildtype_phi", math.nan)),
        "mavenn_yhat": float(mavenn_metrics.get("wildtype_yhat", math.nan)),
        "DMS_score": float(mavenn_metrics.get("wildtype_yhat", math.nan)),
    }
    assay_score_columns = [col for col in score_columns if col in assay_df.columns or col in wt_scores]
    score_map = score_map_for_assay(
        assay_df,
        wildtype=wildtype,
        wt_scores=wt_scores,
        score_columns=assay_score_columns,
    )
    complete = bounded_complete_graphlets(
        score_map,
        score_columns=[str(graphlet_cfg["primary_score_column"])],
        dimensions=[int(item) for item in graphlet_cfg["dimensions"]],
        max_per_dimension=int(graphlet_cfg["max_graphlets_per_dimension"]),
        top_feature_pool=int(graphlet_cfg["top_feature_pool"]),
        max_combinations_per_dimension=int(graphlet_cfg["max_combinations_per_dimension"]),
    )
    rows = []
    for index, (dim, features, checked) in enumerate(complete):
        rows.append(
            {
                "graphlet_id": f"{assay_id}__g{index:04d}__d{dim}",
                "assay_id": assay_id,
                "graphlet_index": int(index),
                "graphlet_split": "calibration" if index % 2 == 0 else "test",
                "graphlet_dimension": int(dim),
                "mutation_features": ";".join(f"{pos}:{aa}" for pos, aa in features),
                "features": features,
                "combinations_checked_for_dimension": int(checked),
            }
        )
    return rows, score_map


def graphlet_feature_sets(features: tuple[tuple[int, str], ...]) -> list[tuple[tuple[int, str], ...]]:
    return [tuple(sorted(mask, key=lambda item: (int(item[0]), str(item[1])))) for mask in subset_masks(features)]


def edge_effects(values: np.ndarray, *, dim: int) -> list[float]:
    arr = np.asarray(values, dtype=np.float64)
    effects = []
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            if node & (1 << bit):
                continue
            effects.append(float(arr[node ^ (1 << bit)] - arr[node]))
    return effects


def dfe_summary(values: np.ndarray, *, dim: int, threshold: float) -> dict[str, float]:
    effects = np.asarray(edge_effects(values, dim=dim), dtype=np.float64)
    if effects.size == 0:
        return {}
    return {
        **summarize_numeric(effects, "dfe_edge_effect"),
        "dfe_fraction_beneficial": float(np.mean(effects > float(threshold))),
        "dfe_fraction_deleterious": float(np.mean(effects < -float(threshold))),
        "dfe_fraction_neutral": float(np.mean(np.abs(effects) <= float(threshold))),
    }


def magnitude_epistasis_summary(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    magnitudes: list[float] = []
    signed: list[float] = []
    for var_i, var_j in itertools.combinations(range(int(dim)), 2):
        other_vars = [idx for idx in range(int(dim)) if idx not in {var_i, var_j}]
        for bg_bits in itertools.product([0, 1], repeat=len(other_vars)):
            base = 0
            for idx, bit in zip(other_vars, bg_bits):
                if bit:
                    base |= 1 << idx
            epsilon = (
                float(arr[base | (1 << var_i) | (1 << var_j)])
                - float(arr[base | (1 << var_i)])
                - float(arr[base | (1 << var_j)])
                + float(arr[base])
            )
            signed.append(epsilon)
            magnitudes.append(abs(epsilon))
    return {
        **summarize_numeric(magnitudes, "magnitude_epistasis"),
        **summarize_numeric(signed, "signed_pairwise_epistasis"),
    }


def graphlet_metric_rows(
    *,
    graphlet: dict[str, Any],
    values: np.ndarray,
    candidate_model: str,
    parameter_set_id: str,
    config: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    cfg = {"analysis": config.get("analysis", {})}
    metric, walsh, gft, _access = metric_rows_for_values(
        graphlet=graphlet,
        values=np.asarray(values, dtype=np.float64),
        candidate_model=candidate_model,
        parameter_set_id=parameter_set_id,
        cfg=cfg,
    )
    dim = int(graphlet["graphlet_dimension"])
    metric.update(local_ruggedness_metrics_for_values(np.asarray(values, dtype=np.float64), dim=dim))
    metric.update(dfe_summary(values, dim=dim, threshold=float(config["graphlets"]["dfe_neutral_threshold"])))
    metric.update(magnitude_epistasis_summary(values, dim=dim))
    return metric, walsh, gft


def top_alignment_coefficients(
    *,
    alignment_path: Path,
    max_pairs: int,
    max_triplets: int,
    seed: int,
) -> AlignmentCoefficients:
    profile = AlignmentProfile.from_alignment(alignment_path, max_sequences=512)
    top_pairs_raw = profile.get_top_couplings(n_pairs=max(int(max_pairs), 1))
    max_score = max((float(score) for _i, _j, score in top_pairs_raw), default=1.0)
    rng = np.random.RandomState(int(seed) % (2**32))
    pair_coefficients: dict[tuple[int, int], float] = {}
    pair_scores: dict[tuple[int, int], float] = {}
    top_pairs: list[dict[str, float]] = []
    for i, j, score in top_pairs_raw:
        pair = tuple(sorted((int(i), int(j))))
        normalized = float(score) / max(max_score, 1e-12)
        coeff = float(rng.normal(scale=math.sqrt(max(normalized, 1e-8))))
        pair_coefficients[pair] = coeff
        pair_scores[pair] = normalized
        top_pairs.append({"position_i": float(pair[0]), "position_j": float(pair[1]), "score": float(score), "normalized_score": normalized, "coefficient": coeff})

    neighbor_map: dict[int, set[int]] = {}
    for i, j in pair_scores:
        neighbor_map.setdefault(i, set()).add(j)
        neighbor_map.setdefault(j, set()).add(i)
    triplet_candidates: list[tuple[tuple[int, int, int], float]] = []
    positions = sorted(neighbor_map)
    for triplet in itertools.combinations(positions, 3):
        pairs = [tuple(sorted(pair)) for pair in itertools.combinations(triplet, 2)]
        score = sum(pair_scores.get(pair, 0.0) for pair in pairs)
        if score <= 0:
            continue
        if sum(pair in pair_scores for pair in pairs) < 2:
            continue
        triplet_candidates.append((tuple(int(x) for x in triplet), float(score)))
    triplet_candidates.sort(key=lambda item: (-item[1], item[0]))
    triplet_coefficients: dict[tuple[int, int, int], float] = {}
    max_triplet_score = max((score for _triplet, score in triplet_candidates[:max_triplets]), default=1.0)
    for triplet, score in triplet_candidates[: int(max_triplets)]:
        triplet_coefficients[triplet] = float(rng.normal(scale=math.sqrt(score / max(max_triplet_score, 1e-12))))

    return AlignmentCoefficients(
        pair_coefficients=pair_coefficients,
        triplet_coefficients=triplet_coefficients,
        top_pairs=top_pairs,
    )


def raw_position_module_values(
    feature_sets: list[tuple[tuple[int, str], ...]],
    coefficients: dict[tuple[int, ...], float],
) -> np.ndarray:
    values = np.zeros(len(feature_sets), dtype=np.float64)
    if not coefficients:
        return values
    for index, features in enumerate(feature_sets):
        active_positions = {int(pos) for pos, _aa in features}
        total = 0.0
        for positions, coeff in coefficients.items():
            if all(int(pos) in active_positions for pos in positions):
                total += float(coeff)
        values[index] = total
    return values


def standardize(raw: np.ndarray, center: float, scale: float) -> np.ndarray:
    return (np.asarray(raw, dtype=np.float64) - float(center)) / max(float(scale), 1e-8)


def coupled_standardization(
    *,
    base_values: np.ndarray,
    feature_sets: list[tuple[tuple[int, str], ...]],
    coeffs: AlignmentCoefficients,
) -> CoupledStandardization:
    raw2 = raw_position_module_values(feature_sets, coeffs.pair_coefficients)
    raw3 = raw_position_module_values(feature_sets, coeffs.triplet_coefficients)
    return CoupledStandardization(
        base_center=float(np.mean(base_values)) if base_values.size else 0.0,
        base_sd=max(float(np.std(base_values)), 1e-8) if base_values.size else 1.0,
        raw2_center=float(np.mean(raw2)) if raw2.size else 0.0,
        raw2_sd=max(float(np.std(raw2)), 1e-8) if raw2.size else 1.0,
        raw3_center=float(np.mean(raw3)) if raw3.size else 0.0,
        raw3_sd=max(float(np.std(raw3)), 1e-8) if raw3.size else 1.0,
    )


def apply_coupled_transform(
    *,
    base_values: np.ndarray,
    feature_sets: list[tuple[tuple[int, str], ...]],
    coeffs: AlignmentCoefficients,
    stats: CoupledStandardization,
    parameters: dict[str, Any],
    allow_three_way: bool,
) -> np.ndarray:
    base = np.asarray(base_values, dtype=np.float64)
    raw2 = raw_position_module_values(feature_sets, coeffs.pair_coefficients)
    z2 = standardize(raw2, stats.raw2_center, stats.raw2_sd)
    scale2 = float(parameters.get("alignment_pair_scale", 0.0))
    base_rescale = float(parameters.get("base_rescale", 1.0))
    out = stats.base_center + base_rescale * (base - stats.base_center) + stats.base_sd * scale2 * z2
    if allow_three_way:
        raw3 = raw_position_module_values(feature_sets, coeffs.triplet_coefficients)
        z3 = standardize(raw3, stats.raw3_center, stats.raw3_sd)
        scale3 = float(parameters.get("alignment_triplet_scale", 0.0))
        out = out + stats.base_sd * scale3 * z3
    return np.asarray(out, dtype=np.float64)


def sample_variants(
    *,
    assay_df: pd.DataFrame,
    wildtype: str,
    score_column: str,
    config: dict[str, Any],
    rng: np.random.RandomState,
) -> pd.DataFrame:
    cfg = dict(config["variant_validation"])
    df = assay_df.copy()
    df = df[df["mutated_sequence"].map(lambda seq: len(str(seq)) == len(wildtype))].copy()
    if "mutation_count" not in df.columns:
        df["mutation_count"] = [
            len(mutation_set(wildtype, str(seq))) for seq in df["mutated_sequence"].astype(str)
        ]
    df = df[np.isfinite(pd.to_numeric(df[score_column], errors="coerce"))].copy()
    frames = []
    for mutation_count, max_count in [
        (1, int(cfg["max_single_variants"])),
        (2, int(cfg["max_double_variants"])),
        (3, int(cfg["max_triple_variants"])),
    ]:
        sub = df[df["mutation_count"] == mutation_count].copy()
        if max_count > 0 and len(sub) > max_count:
            sub = sub.sample(n=max_count, random_state=int(rng.randint(0, 2**31 - 1)))
        frames.append(sub)
    multi = df[df["mutation_count"] >= 4].copy()
    max_multi = int(cfg.get("max_high_order_variants", 0))
    if max_multi > 0 and len(multi) > max_multi:
        multi = multi.sample(n=max_multi, random_state=int(rng.randint(0, 2**31 - 1)))
    frames.append(multi)
    out = pd.concat(frames, ignore_index=True) if frames else df.iloc[0:0].copy()
    holdout_modulo = max(int(cfg.get("holdout_modulo", 5)), 2)
    holdout_remainder = int(cfg.get("holdout_remainder", 0)) % holdout_modulo
    out["validation_split"] = [
        "holdout" if hash_mod(str(seq), holdout_modulo) == holdout_remainder else "calibration"
        for seq in out["mutated_sequence"].astype(str)
    ]
    return out.reset_index(drop=True)


def validation_rows_for_predictions(
    *,
    assay_id: str,
    sample: pd.DataFrame,
    score_column: str,
    predictions: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for model_name, predicted in predictions.items():
        for split in ["calibration", "holdout", "all"]:
            split_df = sample if split == "all" else sample[sample["validation_split"] == split]
            if split_df.empty:
                continue
            split_pred = np.asarray(predicted[split_df.index.to_numpy(dtype=int)], dtype=np.float64)
            observed = split_df[score_column].to_numpy(dtype=np.float64)
            for label, sub_mask in [
                ("all", np.ones(len(split_df), dtype=bool)),
                ("single", split_df["mutation_count"].to_numpy(dtype=int) == 1),
                ("double", split_df["mutation_count"].to_numpy(dtype=int) == 2),
                ("triple", split_df["mutation_count"].to_numpy(dtype=int) == 3),
                ("high_order", split_df["mutation_count"].to_numpy(dtype=int) >= 4),
            ]:
                if int(np.sum(sub_mask)) == 0:
                    continue
                rows.append(
                    {
                        "assay_id": assay_id,
                        "candidate_model": model_name,
                        "validation_split": split,
                        "mutation_class": label,
                        "n_variants": int(np.sum(sub_mask)),
                        "spearman": rank_spearman(observed[sub_mask], split_pred[sub_mask]),
                        "nrmse": nrmse(observed[sub_mask], split_pred[sub_mask]),
                        "empirical_score_column": score_column,
                    }
                )
    return rows


def metric_distance_from_rows(
    *,
    model_rows: list[dict[str, Any]],
    target_rows: list[dict[str, Any]],
    metric_columns: list[str],
) -> tuple[float, dict[str, float]]:
    model_df = pd.DataFrame(model_rows)
    target_df = pd.DataFrame(target_rows)
    pieces: list[float] = []
    by_metric: dict[str, float] = {}
    for metric in metric_columns:
        if metric not in model_df or metric not in target_df:
            continue
        model_values = pd.to_numeric(model_df[metric], errors="coerce").dropna()
        target_values = pd.to_numeric(target_df[metric], errors="coerce").dropna()
        if model_values.empty or target_values.empty:
            continue
        model_mean = float(model_values.mean())
        target_mean = float(target_values.mean())
        scale = max(float(target_values.std(ddof=0)), abs(target_mean) * 0.10, 1e-6)
        value = abs(model_mean - target_mean) / scale
        by_metric[metric] = float(value)
        pieces.append(float(value))
    return (float(np.mean(pieces)) if pieces else math.nan), by_metric


def parameter_specs(rows: list[dict[str, Any]]) -> list[ParameterSpec]:
    return [
        ParameterSpec(
            name=str(row["name"]),
            kind=str(row["kind"]),
            low=float(row["low"]),
            high=float(row["high"]),
            initial_scale_fraction=float(row.get("initial_scale_fraction", 0.15)),
            min_scale_fraction=float(row.get("min_scale_fraction", 0.05)),
        )
        for row in rows
    ]


def s_mc_cfg(cfg: dict[str, Any]) -> SMCABCConfig:
    return SMCABCConfig(
        n_particles=int(cfg["n_particles"]),
        n_rounds=int(cfg["n_rounds"]),
        initial_pool_multiplier=int(cfg["initial_pool_multiplier"]),
        proposal_pool_multiplier=int(cfg["proposal_pool_multiplier"]),
        epsilon_quantile=float(cfg["epsilon_quantile"]),
        seed=int(cfg["seed"]),
    )


def backend_cfg(cfg: dict[str, Any]) -> SMCABCBackendConfig:
    return SMCABCBackendConfig(
        kind=str(cfg.get("kind", "serial")),
        max_workers=int(cfg.get("max_workers", 1)),
        batch_size=int(cfg.get("batch_size", 1)),
    )


def fit_coupled_layer(
    *,
    assay_id: str,
    layer_name: str,
    allow_three_way: bool,
    graphlets: list[dict[str, Any]],
    empirical_values_by_graphlet: dict[str, np.ndarray],
    base_values_by_graphlet: dict[str, np.ndarray],
    feature_sets_by_graphlet: dict[str, list[tuple[tuple[int, str], ...]]],
    variant_sample: VariantSample,
    coeffs: AlignmentCoefficients,
    stats: CoupledStandardization,
    config: dict[str, Any],
    checkpoint_path: Path,
    progress_callback: Any,
) -> tuple[dict[str, Any], dict[str, Any]]:
    fit_cfg = dict(config["coupled_fit"][layer_name])
    specs = parameter_specs(list(fit_cfg["parameters"]))
    metric_columns = [str(item) for item in config["graphlets"]["metric_columns"]]
    primary_score = str(config["graphlets"]["primary_score_column"])
    calibration_graphlets = [g for g in graphlets if str(g["graphlet_split"]) == "calibration"]
    if not calibration_graphlets:
        calibration_graphlets = list(graphlets)

    target_metric_rows: list[dict[str, Any]] = []
    for graphlet in calibration_graphlets:
        graphlet_payload = {**graphlet, "split": graphlet["graphlet_split"], "score_column": primary_score}
        metric, _walsh, _gft = graphlet_metric_rows(
            graphlet=graphlet_payload,
            values=empirical_values_by_graphlet[str(graphlet["graphlet_id"])],
            candidate_model=f"empirical_observed__{primary_score}",
            parameter_set_id=primary_score,
            config=config,
        )
        target_metric_rows.append(metric)

    sample = variant_sample.frame
    calibration_mask = sample["validation_split"].to_numpy(dtype=object) == "calibration"
    empirical_variant = sample[primary_score].to_numpy(dtype=np.float64)
    base_variant = variant_sample.base_values
    feature_variant = variant_sample.feature_sets
    weight_graphlets = float(fit_cfg.get("graphlet_distance_weight", 1.0))
    weight_variant = float(fit_cfg.get("variant_rank_weight", 0.25))
    corr_min = float(fit_cfg.get("base_correlation_min", 0.35))
    corr_penalty_weight = float(fit_cfg.get("base_correlation_penalty_weight", 2.0))

    def simulate(parameters: dict[str, float | int]) -> SimulationResult:
        candidate_metric_rows: list[dict[str, Any]] = []
        for graphlet in calibration_graphlets:
            graphlet_id = str(graphlet["graphlet_id"])
            values = apply_coupled_transform(
                base_values=base_values_by_graphlet[graphlet_id],
                feature_sets=feature_sets_by_graphlet[graphlet_id],
                coeffs=coeffs,
                stats=stats,
                parameters=parameters,
                allow_three_way=allow_three_way,
            )
            graphlet_payload = {**graphlet, "split": graphlet["graphlet_split"], "score_column": primary_score}
            metric, _walsh, _gft = graphlet_metric_rows(
                graphlet=graphlet_payload,
                values=values,
                candidate_model=layer_name,
                parameter_set_id="candidate",
                config=config,
            )
            candidate_metric_rows.append(metric)
        graphlet_distance, by_metric = metric_distance_from_rows(
            model_rows=candidate_metric_rows,
            target_rows=target_metric_rows,
            metric_columns=metric_columns,
        )
        variant_pred = apply_coupled_transform(
            base_values=base_variant,
            feature_sets=feature_variant,
            coeffs=coeffs,
            stats=stats,
            parameters=parameters,
            allow_three_way=allow_three_way,
        )
        variant_rho = rank_spearman(empirical_variant[calibration_mask], variant_pred[calibration_mask])
        variant_distance = 1.0 if not np.isfinite(variant_rho) else float((1.0 - variant_rho) / 2.0)
        base_corr = rank_spearman(base_variant[calibration_mask], variant_pred[calibration_mask])
        corr_penalty = 0.0
        if np.isfinite(base_corr) and base_corr < corr_min:
            corr_penalty = corr_penalty_weight * float((corr_min - base_corr) ** 2)
        elif not np.isfinite(base_corr):
            corr_penalty = corr_penalty_weight
        distance = weight_graphlets * graphlet_distance + weight_variant * variant_distance + corr_penalty
        vector = np.asarray(
            [
                graphlet_distance,
                variant_distance,
                base_corr if np.isfinite(base_corr) else -1.0,
                corr_penalty,
            ],
            dtype=np.float64,
        )
        return SimulationResult(
            distance=float(distance),
            summary_vector=vector,
            extras={
                "graphlet_distance": graphlet_distance,
                "variant_distance": variant_distance,
                "variant_calibration_spearman": variant_rho,
                "base_spearman": base_corr,
                "correlation_penalty": corr_penalty,
                "component_distances": by_metric,
            },
        )

    result = run_smc_abc(
        specs=specs,
        config=s_mc_cfg(dict(fit_cfg["smc"])),
        simulate=simulate,
        backend=backend_cfg(dict(fit_cfg.get("backend", {"kind": "serial"}))),
        checkpoint=SMCABCCheckpointConfig(path=str(checkpoint_path), resume=bool(fit_cfg.get("resume_checkpoints", True))),
        progress_callback=progress_callback,
        run_label=f"{assay_id}::{layer_name}",
    )
    best = result.best_particle
    mean_params: dict[str, float | int] = {}
    weights = np.asarray([particle.weight for particle in result.particles], dtype=np.float64)
    weights = weights / max(float(np.sum(weights)), 1e-12)
    for spec in specs:
        vals = np.asarray([particle.parameters[spec.name] for particle in result.particles], dtype=np.float64)
        mean = float(np.sum(weights * vals))
        mean_params[spec.name] = int(round(mean)) if spec.kind == "int" else mean
    summary = {
        "assay_id": assay_id,
        "layer_name": layer_name,
        "allow_three_way": allow_three_way,
        "best_particle": {
            "distance": best.distance,
            "parameters": dict(best.parameters),
            "extras": dict(best.extras),
        },
        "posterior_mean_parameters": mean_params,
        "round_summaries": [
            {
                "round_index": item.round_index,
                "epsilon": item.epsilon,
                "attempts": item.attempts,
                "best_distance": item.best_distance,
                "median_distance": item.median_distance,
            }
            for item in result.round_summaries
        ],
    }
    particles = {
        "rows": [
            {
                "assay_id": assay_id,
                "layer_name": layer_name,
                "distance": particle.distance,
                "weight": particle.weight,
                **{name: particle.parameters[name] for name in sorted(particle.parameters)},
                **{f"extra__{key}": value for key, value in particle.extras.items() if not isinstance(value, dict)},
            }
            for particle in result.particles
        ]
    }
    return summary, particles


def local_node_epistasis(values: np.ndarray, *, dim: int, node: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    single_effects = np.asarray([float(arr[node ^ (1 << bit)] - arr[node]) for bit in range(dim)], dtype=np.float64)
    magnitudes = []
    sign_events = 0
    reciprocal_events = 0
    pair_count = 0
    for i, j in itertools.combinations(range(dim), 2):
        fij = float(arr[node ^ (1 << i) ^ (1 << j)])
        fi = float(arr[node ^ (1 << i)])
        fj = float(arr[node ^ (1 << j)])
        f0 = float(arr[node])
        eps = fij - fi - fj + f0
        magnitudes.append(abs(eps))
        effect_i = fi - f0
        effect_i_after_j = fij - fj
        effect_j = fj - f0
        effect_j_after_i = fij - fi
        sign_i = np.sign(effect_i) != 0 and np.sign(effect_i) != np.sign(effect_i_after_j)
        sign_j = np.sign(effect_j) != 0 and np.sign(effect_j) != np.sign(effect_j_after_i)
        sign_events += int(sign_i or sign_j)
        reciprocal_events += int(sign_i and sign_j)
        pair_count += 1
    beneficial = single_effects[single_effects > 0]
    return {
        "node_fitness": float(arr[node]),
        "node_rank_fraction": float(np.mean(arr <= arr[node])) if arr.size else math.nan,
        "beneficial_neighbor_fraction": float(np.mean(single_effects > 0)) if single_effects.size else math.nan,
        "mean_beneficial_effect": float(np.mean(beneficial)) if beneficial.size else 0.0,
        "mean_single_effect": float(np.mean(single_effects)) if single_effects.size else math.nan,
        "sign_epistasis_fraction_local": float(sign_events / pair_count) if pair_count else math.nan,
        "reciprocal_sign_epistasis_fraction_local": float(reciprocal_events / pair_count) if pair_count else math.nan,
        "magnitude_epistasis_mean_local": float(np.mean(magnitudes)) if magnitudes else math.nan,
    }


def walk_rows_for_values(
    *,
    assay_id: str,
    graphlet: dict[str, Any],
    values: np.ndarray,
    candidate_model: str,
    parameter_set_id: str,
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    cfg = dict(config["walks"])
    dim = int(graphlet["graphlet_dimension"])
    if dim <= 0:
        return []
    arr = np.asarray(values, dtype=np.float64)
    rows: list[dict[str, Any]] = []
    walsh_rows, fractions = graphlet_runner.walsh_decomposition(arr, dim=dim)
    gft_rows, gft_summary = graphlet_runner.gft_rows_from_walsh(
        walsh_rows,
        graphlet=graphlet,
        parameter_set={"candidate_model": candidate_model, "parameter_set_id": parameter_set_id},
        total_variance=float(np.var(arr)),
    )
    del gft_rows
    high_order = float(sum(value for order, value in fractions.items() if int(order) >= 3))
    peak_node = int(np.argmax(arr))
    rng = np.random.RandomState(graphlet_runner.stable_seed("walks", assay_id, graphlet["graphlet_id"], candidate_model))
    max_steps = int(cfg["steps"])
    for walk_type, n_walks in [
        ("random", int(cfg["n_random_walks_per_graphlet"])),
        ("adaptive_greedy", int(cfg["n_adaptive_walks_per_graphlet"])),
    ]:
        for walk_index in range(n_walks):
            node = int(rng.randint(0, 1 << dim))
            for step in range(max_steps + 1):
                row = {
                    "assay_id": assay_id,
                    "graphlet_id": graphlet["graphlet_id"],
                    "graphlet_split": graphlet["graphlet_split"],
                    "graphlet_dimension": dim,
                    "candidate_model": candidate_model,
                    "parameter_set_id": parameter_set_id,
                    "walk_type": walk_type,
                    "walk_index": walk_index,
                    "step_index": step,
                    "node": node,
                    "hamming_distance_to_graphlet_peak": int((node ^ peak_node).bit_count()),
                    "walsh_high_order_fraction_graphlet": high_order,
                    **{f"{key}_graphlet": value for key, value in gft_summary.items()},
                    **local_node_epistasis(arr, dim=dim, node=node),
                }
                rows.append(row)
                if step == max_steps:
                    break
                neighbors = [node ^ (1 << bit) for bit in range(dim)]
                if walk_type == "random":
                    node = int(neighbors[int(rng.randint(0, len(neighbors)))])
                else:
                    best = max(neighbors, key=lambda item: float(arr[item]))
                    if float(arr[best]) <= float(arr[node]):
                        break
                    node = int(best)
    return rows


def process_system(
    *,
    system: dict[str, Any],
    config: dict[str, Any],
    project_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    assay_id = str(system["assay_id"])
    system_dir = output_dir / "systems" / assay_id
    system_dir.mkdir(parents=True, exist_ok=True)
    system_summary_path = system_dir / "system_summary.json"
    if system_summary_path.is_file() and bool(config.get("resume", {}).get("reuse_completed_systems", True)):
        return json.loads(system_summary_path.read_text(encoding="utf-8"))

    emit_progress(output_dir, {"event": "system_start", "assay_id": assay_id})
    fit_summary = run_base_fit(
        config=config,
        assay_id=assay_id,
        project_root=project_root,
        output_dir=output_dir,
    )
    fit_dir = system_dir / "base_fit"
    panel, assay_df, mavenn_metrics = load_prepared_assay(fit_dir, assay_id)
    wildtype = canonical_sequence(str(panel.loc[panel["DMS_id"] == assay_id, "target_seq"].iloc[0]))
    alignment_path = project_root / config["mmseqs"]["alignment_dir"] / f"{assay_id}.fasta"
    if not alignment_path.exists() and fit_summary.get("panel"):
        candidate = Path(str(fit_summary["panel"][0].get("alignment_path", "")))
        if candidate.exists():
            alignment_path = candidate
    if not alignment_path.exists():
        raise FileNotFoundError(f"Missing alignment for {assay_id}: {alignment_path}")
    shutil.copy2(alignment_path, system_dir / f"{assay_id}.alignment.fasta")

    branches = fit_summary["branches"]
    base_updates = dict(branches["smc_abc_best_raw"]["fitted_parameters"])
    posterior_mean_updates = dict(branches["smc_abc_posterior_mean_raw"]["fitted_parameters"])
    seed = int(config["adapt_env"]["base_config"].get("seed", 42))
    base_config = dict(config["adapt_env"]["base_config"])
    base_landscape = build_landscape(
        base_config=base_config,
        updates=base_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed,
    )
    posterior_mean_landscape = build_landscape(
        base_config=base_config,
        updates=posterior_mean_updates,
        wildtype=wildtype,
        alignment_path=alignment_path,
        seed=seed + 1,
    )

    primary_score = str(config["graphlets"].get("primary_score_column", PRIMARY_SCORE_DEFAULT))
    rng = np.random.RandomState(graphlet_runner.stable_seed("variant-sample", assay_id))
    variant_df = sample_variants(
        assay_df=assay_df,
        wildtype=wildtype,
        score_column=primary_score,
        config=config,
        rng=rng,
    )
    variant_sequences = variant_df["mutated_sequence"].astype(str).tolist()
    base_variant_values = evaluate_landscape_sequences(
        base_landscape,
        variant_sequences,
        batch_size=int(config["variant_validation"]["evaluation_batch_size"]),
    )
    posterior_mean_variant_values = evaluate_landscape_sequences(
        posterior_mean_landscape,
        variant_sequences,
        batch_size=int(config["variant_validation"]["evaluation_batch_size"]),
    )
    variant_features = variant_feature_sets(wildtype, variant_sequences)
    variant_sample = VariantSample(
        frame=variant_df,
        base_values=base_variant_values,
        feature_sets=variant_features,
    )

    graphlets, score_map = build_graphlet_catalog(
        assay_id=assay_id,
        wildtype=wildtype,
        assay_df=assay_df,
        mavenn_metrics=mavenn_metrics,
        config=config,
    )
    write_csv(
        system_dir / "graphlet_catalog.csv",
        [{k: v for k, v in row.items() if k != "features"} for row in graphlets],
    )
    if not graphlets:
        raise ValueError(f"No complete graphlets found for {assay_id}; target is not appropriate")

    empirical_values_by_graphlet: dict[str, np.ndarray] = {}
    base_values_by_graphlet: dict[str, np.ndarray] = {}
    mean_values_by_graphlet: dict[str, np.ndarray] = {}
    feature_sets_by_graphlet: dict[str, list[tuple[tuple[int, str], ...]]] = {}
    for graphlet in graphlets:
        graphlet_id = str(graphlet["graphlet_id"])
        features = tuple(graphlet["features"])
        masks = subset_masks(features)
        empirical_values_by_graphlet[graphlet_id] = np.asarray(
            [score_map[mask][primary_score] for mask in masks],
            dtype=np.float64,
        )
        base_values_by_graphlet[graphlet_id] = evaluate_landscape_on_graphlet(
            base_landscape,
            wildtype=wildtype,
            features=features,
        )
        mean_values_by_graphlet[graphlet_id] = evaluate_landscape_on_graphlet(
            posterior_mean_landscape,
            wildtype=wildtype,
            features=features,
        )
        feature_sets_by_graphlet[graphlet_id] = graphlet_feature_sets(features)

    coeffs = top_alignment_coefficients(
        alignment_path=alignment_path,
        max_pairs=int(config["coupled_fit"]["alignment_modules"]["max_pairs"]),
        max_triplets=int(config["coupled_fit"]["alignment_modules"]["max_triplets"]),
        seed=graphlet_runner.stable_seed("alignment-coefficients", assay_id),
    )
    write_csv(system_dir / "alignment_top_pairs.csv", coeffs.top_pairs)
    calibration_variant_mask = variant_df["validation_split"].to_numpy(dtype=object) == "calibration"
    stats = coupled_standardization(
        base_values=base_variant_values[calibration_variant_mask],
        feature_sets=[features for features, keep in zip(variant_features, calibration_variant_mask) if keep],
        coeffs=coeffs,
    )

    layer_summaries: dict[str, dict[str, Any]] = {}
    particle_rows: list[dict[str, Any]] = []
    for layer_name, allow_three_way in [("coupled_2way", False), ("coupled_3way", True)]:
        checkpoint_path = system_dir / "checkpoints" / f"smc_{layer_name}.json"

        def progress(payload: dict[str, Any], *, layer: str = layer_name) -> None:
            emit_progress(output_dir, {"event": "coupled_fit_progress", "assay_id": assay_id, "layer": layer, **payload})

        summary, particles = fit_coupled_layer(
            assay_id=assay_id,
            layer_name=layer_name,
            allow_three_way=allow_three_way,
            graphlets=graphlets,
            empirical_values_by_graphlet=empirical_values_by_graphlet,
            base_values_by_graphlet=base_values_by_graphlet,
            feature_sets_by_graphlet=feature_sets_by_graphlet,
            variant_sample=variant_sample,
            coeffs=coeffs,
            stats=stats,
            config=config,
            checkpoint_path=checkpoint_path,
            progress_callback=progress,
        )
        layer_summaries[layer_name] = summary
        particle_rows.extend(particles["rows"])
        write_json(system_dir / f"{layer_name}_summary.json", summary)
    write_csv(system_dir / "coupled_particles.csv", particle_rows)

    model_parameter_sets: list[tuple[str, str, np.ndarray | None, dict[str, Any] | None, bool]] = [
        ("base_abc_best", "smc_abc_best_raw", None, None, False),
        ("base_abc_posterior_mean", "smc_abc_posterior_mean_raw", None, None, False),
        (
            "coupled_2way_best",
            "coupled_2way_best",
            None,
            dict(layer_summaries["coupled_2way"]["best_particle"]["parameters"]),
            False,
        ),
        (
            "coupled_2way_posterior_mean",
            "coupled_2way_posterior_mean",
            None,
            dict(layer_summaries["coupled_2way"]["posterior_mean_parameters"]),
            False,
        ),
        (
            "coupled_3way_best",
            "coupled_3way_best",
            None,
            dict(layer_summaries["coupled_3way"]["best_particle"]["parameters"]),
            True,
        ),
        (
            "coupled_3way_posterior_mean",
            "coupled_3way_posterior_mean",
            None,
            dict(layer_summaries["coupled_3way"]["posterior_mean_parameters"]),
            True,
        ),
    ]

    graphlet_metric_rows_all: list[dict[str, Any]] = []
    walsh_rows_all: list[dict[str, Any]] = []
    gft_rows_all: list[dict[str, Any]] = []
    walk_rows_all: list[dict[str, Any]] = []
    max_walk_graphlets = int(config["walks"]["max_graphlets_per_system"])
    for graphlet_index, graphlet in enumerate(graphlets):
        graphlet_id = str(graphlet["graphlet_id"])
        graphlet_payload = {**graphlet, "split": graphlet["graphlet_split"], "score_column": primary_score}
        empirical_values = empirical_values_by_graphlet[graphlet_id]
        candidates: list[tuple[str, str, np.ndarray]] = [
            (f"empirical_observed__{primary_score}", primary_score, empirical_values),
            ("base_abc_best", "smc_abc_best_raw", base_values_by_graphlet[graphlet_id]),
            ("base_abc_posterior_mean", "smc_abc_posterior_mean_raw", mean_values_by_graphlet[graphlet_id]),
        ]
        for model_name, parameter_id, _unused_values, parameters, allow_three_way in model_parameter_sets[2:]:
            values = apply_coupled_transform(
                base_values=base_values_by_graphlet[graphlet_id],
                feature_sets=feature_sets_by_graphlet[graphlet_id],
                coeffs=coeffs,
                stats=stats,
                parameters=dict(parameters or {}),
                allow_three_way=allow_three_way,
            )
            candidates.append((model_name, parameter_id, values))
        for model_name, parameter_id, values in candidates:
            metric, walsh, gft = graphlet_metric_rows(
                graphlet=graphlet_payload,
                values=values,
                candidate_model=model_name,
                parameter_set_id=parameter_id,
                config=config,
            )
            metric["assay_role"] = system.get("role", "")
            graphlet_metric_rows_all.append(metric)
            walsh_rows_all.extend(walsh)
            gft_rows_all.extend(gft)
            if graphlet_index < max_walk_graphlets:
                walk_rows_all.extend(
                    walk_rows_for_values(
                        assay_id=assay_id,
                        graphlet=graphlet_payload,
                        values=values,
                        candidate_model=model_name,
                        parameter_set_id=parameter_id,
                        config=config,
                    )
                )

    write_csv(system_dir / "graphlet_metrics.csv", graphlet_metric_rows_all)
    write_csv(system_dir / "walsh_terms.csv", walsh_rows_all)
    write_csv(system_dir / "gft_coefficients.csv", gft_rows_all)
    write_csv(system_dir / "walk_step_metrics.csv", walk_rows_all)

    predictions = {
        "base_abc_best": base_variant_values,
        "base_abc_posterior_mean": posterior_mean_variant_values,
    }
    coupling_rows: list[dict[str, Any]] = []
    for model_name, _parameter_id, _unused_values, parameters, allow_three_way in model_parameter_sets[2:]:
        pred = apply_coupled_transform(
            base_values=base_variant_values,
            feature_sets=variant_features,
            coeffs=coeffs,
            stats=stats,
            parameters=dict(parameters or {}),
            allow_three_way=allow_three_way,
        )
        predictions[model_name] = pred
        for split in ["calibration", "holdout", "all"]:
            mask = np.ones(len(variant_df), dtype=bool) if split == "all" else variant_df["validation_split"].to_numpy(dtype=object) == split
            if int(np.sum(mask)) < 3:
                continue
            coupling_rows.append(
                {
                    "assay_id": assay_id,
                    "candidate_model": model_name,
                    "validation_split": split,
                    "base_spearman_to_coupled": rank_spearman(base_variant_values[mask], pred[mask]),
                    "base_pearson_to_coupled": float(np.corrcoef(base_variant_values[mask], pred[mask])[0, 1])
                    if float(np.std(base_variant_values[mask])) > 1e-12 and float(np.std(pred[mask])) > 1e-12
                    else math.nan,
                    "innovation_sd_ratio_to_base": float(np.std(pred[mask] - base_variant_values[mask]) / max(float(np.std(base_variant_values[mask])), 1e-8)),
                }
            )
    variant_validation_rows = validation_rows_for_predictions(
        assay_id=assay_id,
        sample=variant_df,
        score_column=primary_score,
        predictions=predictions,
    )
    write_csv(system_dir / "variant_validation.csv", variant_validation_rows)
    write_csv(system_dir / "coupling_function_stats.csv", coupling_rows)

    distance_rows: list[dict[str, Any]] = []
    metric_columns = [str(item) for item in config["graphlets"]["metric_columns"]]
    for split in ["calibration", "test", "all"]:
        split_rows = graphlet_metric_rows_all if split == "all" else [
            row for row in graphlet_metric_rows_all if str(row.get("graphlet_split")) == split
        ]
        target_rows = [
            row
            for row in split_rows
            if str(row.get("candidate_model")) == f"empirical_observed__{primary_score}"
        ]
        if not target_rows:
            continue
        for model_name in predictions:
            model_rows = [row for row in split_rows if str(row.get("candidate_model")) == model_name]
            if not model_rows:
                continue
            distance, by_metric = metric_distance_from_rows(
                model_rows=model_rows,
                target_rows=target_rows,
                metric_columns=metric_columns,
            )
            distance_rows.append(
                {
                    "assay_id": assay_id,
                    "graphlet_split": split,
                    "candidate_model": model_name,
                    "distance_to_empirical_primary": distance,
                    **{f"distance_component__{key}": value for key, value in by_metric.items()},
                }
            )
    write_csv(system_dir / "model_distances.csv", distance_rows)

    base_validation = branches.get("smc_abc_best_raw", {}).get("validation", {})
    system_summary = {
        "assay_id": assay_id,
        "role": system.get("role"),
        "n_graphlets": int(len(graphlets)),
        "graphlet_dimensions_observed": sorted({int(row["graphlet_dimension"]) for row in graphlets}),
        "n_variant_sample": int(len(variant_df)),
        "n_alignment_pairs": int(len(coeffs.pair_coefficients)),
        "n_alignment_triplets": int(len(coeffs.triplet_coefficients)),
        "base_fit": {
            "best_distance": fit_summary.get("smc_abc", {}).get("best_particle", {}).get("distance"),
            "validation": base_validation,
            "fitted_parameters": base_updates,
        },
        "coupled": layer_summaries,
        "outputs": {
            "graphlet_metrics": str(system_dir / "graphlet_metrics.csv"),
            "variant_validation": str(system_dir / "variant_validation.csv"),
            "walk_step_metrics": str(system_dir / "walk_step_metrics.csv"),
            "model_distances": str(system_dir / "model_distances.csv"),
            "coupling_function_stats": str(system_dir / "coupling_function_stats.csv"),
        },
    }
    write_json(system_summary_path, system_summary)
    emit_progress(output_dir, {"event": "system_complete", "assay_id": assay_id, "n_graphlets": len(graphlets)})
    return system_summary


def aggregate_outputs(output_dir: Path, system_summaries: list[dict[str, Any]]) -> None:
    aggregate_specs = [
        ("system_graphlet_metrics.csv", "graphlet_metrics.csv"),
        ("system_walsh_terms.csv", "walsh_terms.csv"),
        ("system_gft_coefficients.csv", "gft_coefficients.csv"),
        ("system_model_distances.csv", "model_distances.csv"),
        ("system_variant_validation.csv", "variant_validation.csv"),
        ("system_walk_step_metrics.csv", "walk_step_metrics.csv"),
        ("coupling_function_stats.csv", "coupling_function_stats.csv"),
        ("coupled_particles.csv", "coupled_particles.csv"),
    ]
    for aggregate_name, child_name in aggregate_specs:
        rows: list[dict[str, Any]] = []
        for system in system_summaries:
            rows.extend(read_csv_rows(output_dir / "systems" / str(system["assay_id"]) / child_name))
        write_csv(output_dir / aggregate_name, rows)

    fit_rows = []
    coupled_rows = []
    for system in system_summaries:
        assay_id = str(system["assay_id"])
        fit = dict(system.get("base_fit", {}))
        validation = dict(fit.get("validation", {}))
        fit_rows.append(
            {
                "assay_id": assay_id,
                "role": system.get("role"),
                "n_graphlets": system.get("n_graphlets"),
                "base_best_distance": fit.get("best_distance"),
                **{f"base_validation__{key}": value for key, value in validation.items()},
            }
        )
        for layer_name, layer in dict(system.get("coupled", {})).items():
            best = dict(layer.get("best_particle", {}))
            params = dict(best.get("parameters", {}))
            extras = dict(best.get("extras", {}))
            coupled_rows.append(
                {
                    "assay_id": assay_id,
                    "layer_name": layer_name,
                    "best_distance": best.get("distance"),
                    **{f"param__{key}": value for key, value in params.items()},
                    **{f"extra__{key}": value for key, value in extras.items() if not isinstance(value, dict)},
                }
            )
    write_csv(output_dir / "system_fit_summary.csv", fit_rows)
    write_csv(output_dir / "coupled_parameter_summary.csv", coupled_rows)


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(config_path, quick=bool(args.quick))
    write_json(output_dir / "resolved_config.json", config)
    selected = selection_preflight(project_root=project_root, config=config, output_dir=output_dir)
    emit_progress(output_dir, {"event": "preflight_complete", "n_systems": int(len(selected))})
    seed_resume_outputs(output_dir, config)

    system_summaries: list[dict[str, Any]] = []
    for system in config["systems"]:
        summary = process_system(
            system=system,
            config=config,
            project_root=project_root,
            output_dir=output_dir,
        )
        system_summaries.append(summary)
        aggregate_outputs(output_dir, system_summaries)

    aggregate_outputs(output_dir, system_summaries)
    summary = {
        "experiment_id": config["experiment"]["id"],
        "n_systems": int(len(system_summaries)),
        "systems": system_summaries,
        "target_selection": {
            "selected_systems_csv": str(output_dir / "selected_systems.csv"),
            "audit_json": str(output_dir / "target_selection_audit.json"),
        },
        "output_files": {
            name: str(output_dir / name)
            for name in [
                "system_fit_summary.csv",
                "system_variant_validation.csv",
                "system_graphlet_metrics.csv",
                "system_model_distances.csv",
                "system_walk_step_metrics.csv",
                "coupled_parameter_summary.csv",
                "coupling_function_stats.csv",
            ]
        },
    }
    write_json(output_dir / "summary.json", summary)
    emit_progress(output_dir, {"event": "run_complete", "n_systems": len(system_summaries)})
    print(json.dumps(to_builtin(summary), indent=2, sort_keys=True, default=graphlet_runner.json_default))


if __name__ == "__main__":
    main()
