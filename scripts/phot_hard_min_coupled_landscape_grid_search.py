#!/usr/bin/env python3
"""Deterministic PHOT hard-min coupled-landscape grid search.

This HYP-008 test freezes the current PHOT system-specific fit from
RUN-077/RES-030. The primary base ``A`` is the fitted smooth/additive PHOT
landscape. The grid constructs correlated, offset auxiliary landscapes ``B``
and evaluates the weakest-link operator ``F(x) = min(A(x), B(x))`` on both the
full prepared PHOT assay and complete dense empirical graphlets.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import shutil
import sys
from collections import defaultdict
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
from adaptenv.amino_acids import AA_TO_IDX  # noqa: E402
from empirical_graphlet_falsification_suite import (  # noqa: E402
    find_complete_graphlets,
    metric_rows_for_values,
    score_map_for_assay,
    subset_masks,
)
from project_adapt_env.coupled_landscape import (  # noqa: E402
    concentration_summary,
)
from project_adapt_env.coupled_landscape import (  # noqa: E402
    edge_energy_by_bit,
    energy_by_graphlet_feature,
    hard_min_values,
    limiter_switch_metrics,
    node_edge_gradients,
    plaquette_energy_by_pair,
    safe_corr,
    zscore,
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


def feature_set_label(features: tuple[tuple[int, str], ...] | frozenset[tuple[int, str]]) -> str:
    ordered = sorted(tuple(features), key=lambda item: (int(item[0]), str(item[1])))
    return ";".join(feature_label(feature) for feature in ordered)


def mutation_set(wildtype: str, sequence: str) -> frozenset[tuple[int, str]]:
    return frozenset(
        (idx, aa)
        for idx, (wt_aa, aa) in enumerate(zip(wildtype, str(sequence)))
        if aa != wt_aa
    )


def sequences_to_indices(sequences: list[str]) -> np.ndarray:
    return np.asarray([[AA_TO_IDX[str(aa)] for aa in seq] for seq in sequences], dtype=np.int64)


def evaluate_sequences(landscape: Any, sequences: list[str], *, chunk_size: int) -> np.ndarray:
    chunks: list[np.ndarray] = []
    for start in range(0, len(sequences), int(chunk_size)):
        seq_chunk = sequences[start : start + int(chunk_size)]
        chunks.append(np.asarray(landscape.evaluate_batch(sequences_to_indices(seq_chunk)), dtype=np.float64))
    return np.concatenate(chunks) if chunks else np.asarray([], dtype=np.float64)


def stable_unit(*parts: Any) -> float:
    seed = graphlet_runner.stable_seed(*parts)
    return float((int(seed) % 10_000_000) / 10_000_000.0)


def stable_normal(*parts: Any) -> float:
    rng = np.random.RandomState(graphlet_runner.stable_seed(*parts))
    return float(rng.normal())


def rank_spearman(xs: np.ndarray, ys: np.ndarray) -> float:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(ok)) < 3:
        return math.nan
    xr = pd.Series(x[ok]).rank(method="average").to_numpy(dtype=np.float64)
    yr = pd.Series(y[ok]).rank(method="average").to_numpy(dtype=np.float64)
    return safe_corr(xr, yr)


def pearson(xs: np.ndarray, ys: np.ndarray) -> float:
    return safe_corr(np.asarray(xs, dtype=np.float64), np.asarray(ys, dtype=np.float64))


def ks_statistic(xs: np.ndarray, ys: np.ndarray) -> float:
    x = np.sort(np.asarray(xs, dtype=np.float64)[np.isfinite(xs)])
    y = np.sort(np.asarray(ys, dtype=np.float64)[np.isfinite(ys)])
    if x.size == 0 or y.size == 0:
        return math.nan
    values = np.sort(np.unique(np.concatenate([x, y])))
    cdf_x = np.searchsorted(x, values, side="right") / float(x.size)
    cdf_y = np.searchsorted(y, values, side="right") / float(y.size)
    return float(np.max(np.abs(cdf_x - cdf_y)))


def affine_fit(train_pred: np.ndarray, train_obs: np.ndarray) -> tuple[float, float]:
    x = np.asarray(train_pred, dtype=np.float64)
    y = np.asarray(train_obs, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(ok)) < 3 or float(np.std(x[ok])) <= 1e-12:
        return float(np.nanmean(y[ok])) if int(np.sum(ok)) else 0.0, 0.0
    design = np.vstack([np.ones(np.sum(ok)), x[ok]]).T
    beta, *_ = np.linalg.lstsq(design, y[ok], rcond=None)
    return float(beta[0]), float(beta[1])


def apply_affine(pred: np.ndarray, intercept: float, slope: float) -> np.ndarray:
    return float(intercept) + float(slope) * np.asarray(pred, dtype=np.float64)


def nrmse(obs: np.ndarray, pred: np.ndarray) -> float:
    x = np.asarray(obs, dtype=np.float64)
    y = np.asarray(pred, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(ok)) < 1:
        return math.nan
    rmse = float(np.sqrt(np.mean((x[ok] - y[ok]) ** 2)))
    scale = max(float(np.std(x[ok])), 1e-12)
    return float(rmse / scale)


def deterministic_mask(keys: list[str], *, holdout_fraction: float, salt: str) -> np.ndarray:
    values = [stable_unit("holdout", salt, key) < float(holdout_fraction) for key in keys]
    mask = np.asarray(values, dtype=bool)
    if mask.size and (not np.any(mask) or np.all(mask)):
        mask[np.arange(mask.size) % max(int(round(1.0 / max(holdout_fraction, 1e-6))), 2) == 0] = True
        mask[0] = False
    return mask


def effect_summary(values: np.ndarray, *, neutral_threshold: float) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            "mean": math.nan,
            "sd": math.nan,
            "beneficial_fraction": math.nan,
            "neutral_fraction": math.nan,
            "deleterious_fraction": math.nan,
            "skewness": math.nan,
        }
    mean = float(np.mean(arr))
    sd = float(np.std(arr))
    centered = arr - mean
    skew = float(np.mean(centered**3) / (sd**3)) if sd > 1e-12 else 0.0
    threshold = abs(float(neutral_threshold))
    return {
        "mean": mean,
        "sd": sd,
        "beneficial_fraction": float(np.mean(arr > threshold)),
        "neutral_fraction": float(np.mean(np.abs(arr) <= threshold)),
        "deleterious_fraction": float(np.mean(arr < -threshold)),
        "skewness": skew,
    }


def parse_feature_lists(assay_df: pd.DataFrame, wildtype: str) -> list[tuple[tuple[int, str], ...]]:
    features: list[tuple[tuple[int, str], ...]] = []
    for sequence in assay_df["mutated_sequence"].astype(str).tolist():
        features.append(tuple(sorted(mutation_set(wildtype, sequence), key=lambda item: (item[0], item[1]))))
    return features


def build_auxiliary_components(
    feature_lists: list[tuple[tuple[int, str], ...]],
    *,
    seed: int,
    pair_densities: list[float],
) -> tuple[np.ndarray, dict[float, np.ndarray]]:
    unique_features = sorted({feature_label(feature) for features in feature_lists for feature in features})
    feature_coeff = {
        label: stable_normal("hard-min-aux-additive", seed, label)
        for label in unique_features
    }
    additive = np.zeros(len(feature_lists), dtype=np.float64)
    pair_values = {float(density): np.zeros(len(feature_lists), dtype=np.float64) for density in pair_densities}
    for idx, features in enumerate(feature_lists):
        labels = [feature_label(feature) for feature in features]
        additive[idx] = float(sum(feature_coeff[label] for label in labels))
        for left, right in itertools.combinations(sorted(labels), 2):
            pair_label = f"{left}|{right}"
            coeff = stable_normal("hard-min-aux-pair", seed, pair_label)
            unit = stable_unit("hard-min-aux-pair-density", seed, pair_label)
            for density in pair_values:
                if unit <= float(density):
                    pair_values[density][idx] += coeff
    return additive, pair_values


def auxiliary_raw_for_features(
    features: tuple[tuple[int, str], ...],
    *,
    seed: int,
    pair_density: float,
) -> tuple[float, float]:
    labels = [feature_label(feature) for feature in features]
    additive = float(sum(stable_normal("hard-min-aux-additive", seed, label) for label in labels))
    pair = 0.0
    for left, right in itertools.combinations(sorted(labels), 2):
        pair_label = f"{left}|{right}"
        if stable_unit("hard-min-aux-pair-density", seed, pair_label) <= float(pair_density):
            pair += stable_normal("hard-min-aux-pair", seed, pair_label)
    return additive, float(pair)


def fit_auxiliary_transform(
    *,
    a_values: np.ndarray,
    raw_values: np.ndarray,
) -> dict[str, float]:
    a = np.asarray(a_values, dtype=np.float64)
    raw = np.asarray(raw_values, dtype=np.float64)
    a_mean = float(np.mean(a))
    a_sd = max(float(np.std(a)), 1e-12)
    raw_mean = float(np.mean(raw))
    raw_sd = max(float(np.std(raw)), 1e-12)
    za = (a - a_mean) / a_sd
    zraw = (raw - raw_mean) / raw_sd
    corr = safe_corr(za, zraw)
    if not np.isfinite(corr):
        corr = 0.0
    orth = zraw - float(corr) * za
    orth_sd = max(float(np.std(orth)), 1e-12)
    return {
        "a_mean": a_mean,
        "a_sd": a_sd,
        "raw_mean": raw_mean,
        "raw_sd": raw_sd,
        "raw_corr_to_a": float(corr),
        "orth_sd": orth_sd,
    }


def apply_auxiliary_transform(
    *,
    a_values: np.ndarray,
    raw_values: np.ndarray,
    transform: dict[str, float],
    rho: float,
    auxiliary_scale: float,
    offset_sd: float,
) -> np.ndarray:
    a = np.asarray(a_values, dtype=np.float64)
    raw = np.asarray(raw_values, dtype=np.float64)
    za = (a - float(transform["a_mean"])) / max(float(transform["a_sd"]), 1e-12)
    zraw = (raw - float(transform["raw_mean"])) / max(float(transform["raw_sd"]), 1e-12)
    orth = (zraw - float(transform["raw_corr_to_a"]) * za) / max(float(transform["orth_sd"]), 1e-12)
    rho = max(-0.999, min(0.999, float(rho)))
    aux_unit = rho * za + math.sqrt(max(1.0 - rho * rho, 0.0)) * orth
    return (
        float(transform["a_mean"])
        + float(offset_sd) * float(transform["a_sd"])
        + float(auxiliary_scale) * float(transform["a_sd"]) * aux_unit
    )


def parameter_grid(config: dict[str, Any]) -> list[dict[str, Any]]:
    grid = config["hard_min_grid"]
    rows: list[dict[str, Any]] = []
    idx = 0
    for rho, scale, offset, pair_scale, pair_density, seed in itertools.product(
        [float(x) for x in grid["auxiliary_correlations"]],
        [float(x) for x in grid["auxiliary_scales"]],
        [float(x) for x in grid["auxiliary_offset_sds"]],
        [float(x) for x in grid["auxiliary_pairwise_scales"]],
        [float(x) for x in grid["auxiliary_pairwise_densities"]],
        [int(x) for x in grid["auxiliary_seeds"]],
    ):
        idx += 1
        rows.append(
            {
                "candidate_model": "hard_min_correlated_auxiliary",
                "parameter_set_id": f"hard_min__smooth_core_A__grid{idx:04d}",
                "a_anchor": "fitted_smooth_core_only",
                "operator": "hard_min",
                "auxiliary_correlation": rho,
                "auxiliary_scale": scale,
                "auxiliary_offset_sd": offset,
                "auxiliary_pairwise_scale": pair_scale,
                "auxiliary_pairwise_density": pair_density,
                "auxiliary_seed": seed,
            }
        )
    return rows


def assay_prediction_metrics(
    *,
    assay_df: pd.DataFrame,
    feature_lists: list[tuple[tuple[int, str], ...]],
    observed: np.ndarray,
    pred: np.ndarray,
    wt_observed: float,
    wt_pred: float,
    candidate_model: str,
    parameter_set_id: str,
    params: dict[str, Any],
    neutral_threshold: float,
    holdout_fraction: float,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    keys = assay_df["mutant"].astype(str).tolist()
    holdout = deterministic_mask(keys, holdout_fraction=holdout_fraction, salt="phot-assay")
    train = ~holdout
    intercept, slope = affine_fit(pred[train], observed[train])
    mapped = apply_affine(pred, intercept, slope)
    mapped_wt = float(apply_affine(np.asarray([wt_pred], dtype=np.float64), intercept, slope)[0])
    mutation_count = pd.to_numeric(assay_df["mutation_count"], errors="coerce").to_numpy()
    single = mutation_count == 1
    double = mutation_count == 2
    fit_row: dict[str, Any] = {
        "assay_id": ASSAY_ID,
        "candidate_model": candidate_model,
        "parameter_set_id": parameter_set_id,
        **params,
        "n_variants": int(len(assay_df)),
        "n_train": int(np.sum(train)),
        "n_holdout": int(np.sum(holdout)),
        "affine_intercept_train": float(intercept),
        "affine_slope_train": float(slope),
        "validation__functional_ks": ks_statistic(observed, mapped),
        "validation__functional_train_nrmse": nrmse(observed[train], mapped[train]),
        "validation__functional_holdout_nrmse": nrmse(observed[holdout], mapped[holdout]),
        "validation__functional_holdout_spearman": rank_spearman(observed[holdout], pred[holdout]),
        "validation__functional_all_spearman": rank_spearman(observed, pred),
        "validation__functional_all_pearson": pearson(observed, mapped),
        "validation__functional_single_spearman": rank_spearman(observed[single], pred[single]),
        "validation__functional_double_spearman": rank_spearman(observed[double], pred[double]),
        "validation__functional_double_nrmse": nrmse(observed[double], mapped[double]),
        "validation__functional_double_holdout_spearman": rank_spearman(
            observed[double & holdout],
            pred[double & holdout],
        ),
        "validation__functional_double_holdout_nrmse": nrmse(
            observed[double & holdout],
            mapped[double & holdout],
        ),
    }

    obs_by_feature = {features: float(value) for features, value in zip(feature_lists, observed)}
    pred_by_feature = {features: float(value) for features, value in zip(feature_lists, mapped)}
    single_obs_effects: list[float] = []
    single_pred_effects: list[float] = []
    for features in feature_lists:
        if len(features) != 1:
            continue
        single_obs_effects.append(float(obs_by_feature[features] - wt_observed))
        single_pred_effects.append(float(pred_by_feature[features] - mapped_wt))
    obs_effect = np.asarray(single_obs_effects, dtype=np.float64)
    pred_effect = np.asarray(single_pred_effects, dtype=np.float64)
    dfe_row = {
        "assay_id": ASSAY_ID,
        "candidate_model": candidate_model,
        "parameter_set_id": parameter_set_id,
        **params,
        "n_single_effects": int(obs_effect.size),
        "dfe_effect_spearman": rank_spearman(obs_effect, pred_effect),
        "dfe_effect_ks": ks_statistic(obs_effect, pred_effect),
        **{f"empirical_dfe_{key}": value for key, value in effect_summary(obs_effect, neutral_threshold=neutral_threshold).items()},
        **{f"model_dfe_{key}": value for key, value in effect_summary(pred_effect, neutral_threshold=neutral_threshold).items()},
    }

    obs_epi: list[float] = []
    pred_epi: list[float] = []
    obs_sign = 0
    pred_sign = 0
    obs_rse = 0
    pred_rse = 0
    complete = 0
    for features in feature_lists:
        if len(features) != 2 or features not in obs_by_feature:
            continue
        f1 = (features[0],)
        f2 = (features[1],)
        if f1 not in obs_by_feature or f2 not in obs_by_feature:
            continue
        complete += 1
        obs_wt = float(wt_observed)
        pred_wt = float(mapped_wt)
        obs_i = float(obs_by_feature[f1])
        obs_j = float(obs_by_feature[f2])
        obs_ij = float(obs_by_feature[features])
        pred_i = float(pred_by_feature[f1])
        pred_j = float(pred_by_feature[f2])
        pred_ij = float(pred_by_feature[features])
        obs_epi.append(obs_ij - obs_i - obs_j + obs_wt)
        pred_epi.append(pred_ij - pred_i - pred_j + pred_wt)
        obs_i_base = obs_i - obs_wt
        obs_j_base = obs_j - obs_wt
        obs_i_on_j = obs_ij - obs_j
        obs_j_on_i = obs_ij - obs_i
        pred_i_base = pred_i - pred_wt
        pred_j_base = pred_j - pred_wt
        pred_i_on_j = pred_ij - pred_j
        pred_j_on_i = pred_ij - pred_i
        obs_i_flip = obs_i_base * obs_i_on_j < 0
        obs_j_flip = obs_j_base * obs_j_on_i < 0
        pred_i_flip = pred_i_base * pred_i_on_j < 0
        pred_j_flip = pred_j_base * pred_j_on_i < 0
        obs_sign += int(obs_i_flip or obs_j_flip)
        pred_sign += int(pred_i_flip or pred_j_flip)
        obs_rse += int(obs_i_flip and obs_j_flip)
        pred_rse += int(pred_i_flip and pred_j_flip)
    obs_epi_arr = np.asarray(obs_epi, dtype=np.float64)
    pred_epi_arr = np.asarray(pred_epi, dtype=np.float64)
    epi_row = {
        "assay_id": ASSAY_ID,
        "candidate_model": candidate_model,
        "parameter_set_id": parameter_set_id,
        **params,
        "n_complete_double_epistasis": int(complete),
        "validation__functional_epistasis_prediction_spearman": rank_spearman(obs_epi_arr, pred_epi_arr),
        "validation__functional_epistasis_prediction_ks": ks_statistic(obs_epi_arr, pred_epi_arr),
        "validation__functional_epistasis_prediction_nrmse": nrmse(obs_epi_arr, pred_epi_arr),
        "empirical_double_sign_epistasis_fraction": float(obs_sign / complete) if complete else math.nan,
        "model_double_sign_epistasis_fraction": float(pred_sign / complete) if complete else math.nan,
        "empirical_double_reciprocal_sign_epistasis_fraction": float(obs_rse / complete) if complete else math.nan,
        "model_double_reciprocal_sign_epistasis_fraction": float(pred_rse / complete) if complete else math.nan,
        **{f"empirical_epistasis_{key}": value for key, value in effect_summary(obs_epi_arr, neutral_threshold=neutral_threshold).items()},
        **{f"model_epistasis_{key}": value for key, value in effect_summary(pred_epi_arr, neutral_threshold=neutral_threshold).items()},
    }
    fit_row["validation__functional_epistasis_prediction_spearman"] = epi_row[
        "validation__functional_epistasis_prediction_spearman"
    ]
    fit_row["validation__functional_epistasis_prediction_ks"] = epi_row[
        "validation__functional_epistasis_prediction_ks"
    ]
    fit_row["assay_fit_objective_distance"] = float(
        np.nanmean(
            [
                1.0 - float(fit_row["validation__functional_holdout_spearman"]),
                float(fit_row["validation__functional_holdout_nrmse"]),
                float(fit_row["validation__functional_ks"]),
                1.0 - float(dfe_row["dfe_effect_spearman"]),
                float(dfe_row["dfe_effect_ks"]),
                1.0 - float(epi_row["validation__functional_epistasis_prediction_spearman"]),
                float(epi_row["validation__functional_epistasis_prediction_ks"]),
            ]
        )
    )
    return fit_row, dfe_row, epi_row


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
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
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


def summarize_metric_rows(rows: list[dict[str, Any]], *, group_cols: list[str], metric_cols: list[str]) -> list[dict[str, Any]]:
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


def aggregate_feature_energy(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(["candidate_model", "parameter_set_id", "graphlet_split", "energy_kind"], dropna=False):
        energy = sub.groupby("feature_key")["energy"].sum().to_dict()
        if not isinstance(keys, tuple):
            keys = (keys,)
        out.append(
            {
                "candidate_model": keys[0],
                "parameter_set_id": keys[1],
                "graphlet_split": keys[2],
                "energy_kind": keys[3],
                "n_energy_rows": int(len(sub)),
                "n_features_or_pairs": int(len(energy)),
                **concentration_summary(energy, prefix="feature_energy"),
            }
        )
    return out


def graphlet_energy_concentration(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(["candidate_model", "parameter_set_id", "graphlet_split"], dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        edge = sub.groupby("graphlet_id")["edge_innovation_total_energy"].mean().to_dict()
        plaq = sub.groupby("graphlet_id")["plaquette_innovation_total_energy"].mean().to_dict()
        out.append(
            {
                "candidate_model": keys[0],
                "parameter_set_id": keys[1],
                "graphlet_split": keys[2],
                "n_graphlets": int(len(edge)),
                **concentration_summary(edge, prefix="graphlet_edge_energy"),
                **concentration_summary(plaq, prefix="graphlet_plaquette_energy"),
            }
        )
    return out


def locality_covariance_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    metrics = [
        "edge_innovation_total_energy",
        "plaquette_innovation_total_energy",
        "switch_edge_fraction",
        "switch_plaquette_fraction",
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
        for i, j in itertools.combinations(range(len(sub)), 2):
            shared_counts.append(float(len(feature_sets[i] & feature_sets[j])))
            for metric in metrics:
                vi = metric_values[metric][i]
                vj = metric_values[metric][j]
                if np.isfinite(vi) and np.isfinite(vj):
                    diffs[metric].append(-abs(float(vi - vj)))
                else:
                    diffs[metric].append(math.nan)
        shared = np.asarray(shared_counts, dtype=np.float64)
        for metric, values in diffs.items():
            arr = np.asarray(values, dtype=np.float64)
            ok = np.isfinite(shared) & np.isfinite(arr)
            out_row[f"shared_feature_count_vs_{metric}_similarity_spearman"] = (
                rank_spearman(shared[ok], arr[ok]) if int(np.sum(ok)) >= 3 else math.nan
            )
        out.append(out_row)
    return out


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
    assay_score_columns = [col for col in score_columns if col in assay_df.columns or col in wt_scores]
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
    smooth_pred = evaluate_sequences(smooth_core, sequences, chunk_size=chunk_size)
    current_pred = evaluate_sequences(fitted_current, sequences, chunk_size=chunk_size)
    mean_pred = evaluate_sequences(fitted_mean, sequences, chunk_size=chunk_size)
    wt_smooth = float(evaluate_sequences(smooth_core, [wildtype], chunk_size=1)[0])
    wt_current = float(evaluate_sequences(fitted_current, [wildtype], chunk_size=1)[0])
    wt_mean = float(evaluate_sequences(fitted_mean, [wildtype], chunk_size=1)[0])
    wt_observed = float(wt_scores[primary_score])

    parameter_rows = parameter_grid(config)
    pair_densities = sorted({float(row["auxiliary_pairwise_density"]) for row in parameter_rows})
    aux_add, aux_pairs = build_auxiliary_components(
        feature_lists,
        seed=int(config["hard_min_grid"].get("coefficient_seed", 20260516)),
        pair_densities=pair_densities,
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

    auxiliary_transforms: dict[str, dict[str, float]] = {}
    for params in parameter_rows:
        raw = aux_add + float(params["auxiliary_pairwise_scale"]) * aux_pairs[float(params["auxiliary_pairwise_density"])]
        transform = fit_auxiliary_transform(a_values=smooth_pred, raw_values=raw)
        auxiliary_transforms[str(params["parameter_set_id"])] = transform
        b_pred = apply_auxiliary_transform(
            a_values=smooth_pred,
            raw_values=raw,
            transform=transform,
            rho=float(params["auxiliary_correlation"]),
            auxiliary_scale=float(params["auxiliary_scale"]),
            offset_sd=float(params["auxiliary_offset_sd"]),
        )
        values, innovation, limited = hard_min_values(base_values=smooth_pred, auxiliary_values=b_pred)
        extra = {key: value for key, value in params.items() if key not in {"candidate_model", "parameter_set_id"}}
        extra.update(
            {
                "operator": "hard_min",
                "assay_auxiliary_limited_fraction": float(np.mean(limited)),
                "assay_auxiliary_to_a_correlation": safe_corr(b_pred, smooth_pred),
                "assay_innovation_to_a_correlation": safe_corr(innovation, smooth_pred),
                "assay_innovation_to_a_sd_ratio": float(np.std(innovation) / max(float(np.std(smooth_pred)), 1e-12)),
                "auxiliary_raw_corr_to_a": float(transform["raw_corr_to_a"]),
            }
        )
        fit_row, dfe_row, epi_row = assay_prediction_metrics(
            assay_df=assay_df,
            feature_lists=feature_lists,
            observed=observed,
            pred=values,
            wt_observed=wt_observed,
            wt_pred=min(wt_smooth, float(apply_auxiliary_transform(
                a_values=np.asarray([wt_smooth], dtype=np.float64),
                raw_values=np.asarray([0.0], dtype=np.float64),
                transform=transform,
                rho=float(params["auxiliary_correlation"]),
                auxiliary_scale=float(params["auxiliary_scale"]),
                offset_sd=float(params["auxiliary_offset_sd"]),
            )[0])),
            candidate_model=str(params["candidate_model"]),
            parameter_set_id=str(params["parameter_set_id"]),
            params=extra,
            neutral_threshold=neutral_threshold,
            holdout_fraction=holdout_fraction,
        )
        assay_fit_rows.append(fit_row)
        dfe_rows.append(dfe_row)
        epistasis_rows.append(epi_row)

    score_map = score_map_for_assay(
        pd.read_csv(assay_path),
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
        calib_base.extend(np.asarray(graphlet["fitted_smooth_core_only"], dtype=np.float64).tolist())
        calib_empirical.extend(np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64).tolist())
    readout = fit_readout_model(
        calib_base,
        calib_empirical,
        degree=int(config["empirical_local_signal"].get("global_readout_polynomial_degree", 3)),
    )

    random_walk_cfg = {"analysis": config.get("analysis", {})}
    metric_cols = [str(x) for x in config["graphlet_evaluation"]["metric_columns"]]
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
        score_column: str,
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

    for graphlet in graphlets:
        empirical_values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        smooth_values = np.asarray(graphlet["fitted_smooth_core_only"], dtype=np.float64)
        current_values = np.asarray(graphlet["fitted_current_best"], dtype=np.float64)
        empirical_pred = apply_readout_model(smooth_values, readout)
        empirical_resid = empirical_values - empirical_pred

        record_values(
            graphlet=graphlet,
            candidate_model=f"empirical_observed__{primary_score}",
            parameter_set_id=primary_score,
            values=empirical_values,
            score_column=primary_score,
            extra={"operator": "empirical", "a_anchor": ""},
        )
        record_values(
            graphlet=graphlet,
            candidate_model="empirical_global_readout_residual",
            parameter_set_id=f"{primary_score}__minus_poly_smooth_core_A",
            values=empirical_resid,
            score_column=primary_score,
            extra={"operator": "empirical_residual", "a_anchor": "fitted_smooth_core_only"},
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
                score_column=primary_score,
                extra={"operator": "baseline_control", "a_anchor": candidate_model},
            )
        gnk = gnk_values(
            dim=int(graphlet["graphlet_dimension"]),
            seed=graphlet_runner.stable_seed("phot-hard-min-gnk", graphlet["graphlet_id"]),
            target_mean=float(np.mean(empirical_values)),
            target_sd=float(np.std(empirical_values)),
            k=2,
        )
        record_values(
            graphlet=graphlet,
            candidate_model="gnk_local_null",
            parameter_set_id=f"{ASSAY_ID}__gnk_k2",
            values=gnk,
            score_column=primary_score,
            extra={"operator": "gnk_control", "a_anchor": ""},
        )

        empirical_feature_rows = energy_by_graphlet_feature(
            graphlet_id=str(graphlet["graphlet_id"]),
            graphlet_split=str(graphlet["graphlet_split"]),
            graphlet_features=tuple(graphlet["features"]),
            innovation=empirical_resid,
            candidate_model="empirical_global_readout_residual",
            parameter_set_id=f"{primary_score}__minus_poly_smooth_core_A",
        )
        feature_energy_rows.extend(empirical_feature_rows)
        local_signal_rows.append(
            {
                "graphlet_id": graphlet["graphlet_id"],
                "graphlet_split": graphlet["graphlet_split"],
                "mutation_features": graphlet["mutation_features"],
                "candidate_model": "empirical_global_readout_residual",
                "parameter_set_id": f"{primary_score}__minus_poly_smooth_core_A",
                "operator": "empirical_residual",
                "edge_innovation_total_energy": float(sum(edge_energy_by_bit(empirical_resid, dim=int(graphlet["graphlet_dimension"])).values())),
                "plaquette_innovation_total_energy": float(sum(plaquette_energy_by_pair(empirical_resid, dim=int(graphlet["graphlet_dimension"])).values())),
                "switch_edge_fraction": math.nan,
                "switch_plaquette_fraction": math.nan,
                "switch_boundary_plaquette_energy_fraction": math.nan,
                "auxiliary_limited_fraction": math.nan,
            }
        )
        current_innovation = current_values - smooth_values
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
                "switch_edge_fraction": math.nan,
                "switch_plaquette_fraction": math.nan,
                "switch_boundary_plaquette_energy_fraction": math.nan,
                "auxiliary_limited_fraction": math.nan,
            }
        )

        for params in parameter_rows:
            node_features = [tuple(sorted(mask, key=lambda item: (item[0], item[1]))) for mask in subset_masks(tuple(graphlet["features"]))]
            raw_add_pair = [
                auxiliary_raw_for_features(
                    features,
                    seed=int(config["hard_min_grid"].get("coefficient_seed", 20260516)),
                    pair_density=float(params["auxiliary_pairwise_density"]),
                )
                for features in node_features
            ]
            raw_values = np.asarray(
                [
                    add + float(params["auxiliary_pairwise_scale"]) * pair
                    for add, pair in raw_add_pair
                ],
                dtype=np.float64,
            )
            b_values = apply_auxiliary_transform(
                a_values=smooth_values,
                raw_values=raw_values,
                transform=auxiliary_transforms[str(params["parameter_set_id"])],
                rho=float(params["auxiliary_correlation"]),
                auxiliary_scale=float(params["auxiliary_scale"]),
                offset_sd=float(params["auxiliary_offset_sd"]),
            )
            values, innovation, limited = hard_min_values(base_values=smooth_values, auxiliary_values=b_values)
            dim = int(graphlet["graphlet_dimension"])
            switch = limiter_switch_metrics(limited, innovation, dim=dim)
            extra = {
                key: value
                for key, value in params.items()
                if key not in {"candidate_model", "parameter_set_id"}
            }
            extra.update(
                {
                    "auxiliary_to_a_correlation": safe_corr(b_values, smooth_values),
                    "values_to_a_correlation": safe_corr(values, smooth_values),
                    "edge_gradient_correlation_to_a": safe_corr(
                        node_edge_gradients(values, dim=dim),
                        node_edge_gradients(smooth_values, dim=dim),
                    ),
                    "innovation_to_a_sd_ratio": float(np.std(innovation) / max(float(np.std(smooth_values)), 1e-12)),
                    **switch,
                }
            )
            record_values(
                graphlet=graphlet,
                candidate_model=str(params["candidate_model"]),
                parameter_set_id=str(params["parameter_set_id"]),
                values=values,
                score_column=primary_score,
                extra=extra,
            )
            feature_energy_rows.extend(
                energy_by_graphlet_feature(
                    graphlet_id=str(graphlet["graphlet_id"]),
                    graphlet_split=str(graphlet["graphlet_split"]),
                    graphlet_features=tuple(graphlet["features"]),
                    innovation=innovation,
                    candidate_model=str(params["candidate_model"]),
                    parameter_set_id=str(params["parameter_set_id"]),
                )
            )
            local_signal_rows.append(
                {
                    "graphlet_id": graphlet["graphlet_id"],
                    "graphlet_split": graphlet["graphlet_split"],
                    "mutation_features": graphlet["mutation_features"],
                    "candidate_model": str(params["candidate_model"]),
                    "parameter_set_id": str(params["parameter_set_id"]),
                    **extra,
                    "edge_innovation_total_energy": float(sum(edge_energy_by_bit(innovation, dim=dim).values())),
                    "plaquette_innovation_total_energy": float(sum(plaquette_energy_by_pair(innovation, dim=dim).values())),
                }
            )

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
                    "distance_to_empirical_primary": distance,
                    **{
                        key: first.get(key, "")
                        for key in [
                            "operator",
                            "a_anchor",
                            "auxiliary_correlation",
                            "auxiliary_scale",
                            "auxiliary_offset_sd",
                            "auxiliary_pairwise_scale",
                            "auxiliary_pairwise_density",
                            "auxiliary_seed",
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
            "auxiliary_limited_fraction",
            "switch_edge_fraction",
            "switch_plaquette_fraction",
            "switch_boundary_plaquette_energy_fraction",
            "innovation_to_a_sd_ratio",
            "values_to_a_correlation",
            "edge_gradient_correlation_to_a",
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
    feature_energy_summary = aggregate_feature_energy(feature_energy_rows)
    graphlet_energy_summary = graphlet_energy_concentration(local_signal_rows)
    covariance_summary = locality_covariance_summary(local_signal_rows)

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

    write_csv(output_dir / "parameter_grid.csv", parameter_rows)
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

    source_phot_best = math.nan
    for system_summary in source_summary.get("systems", []):
        if str(system_summary.get("assay_id")) == ASSAY_ID:
            source_phot_best = float(system_summary["fit_summary"]["best_distance"])
            break
    preflight = {
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_output_dir": str(source_output),
        "source_results": ["RES-030", "RES-031"],
        "mode": "deterministic_frozen_smooth_core_A_hard_min_grid_no_abc",
        "operator": "f(x)=min(a(x),b(x))",
        "a_landscape": "fitted_smooth_core_only from RUN-077/RES-030 PHOT best branch",
        "primary_score_column": primary_score,
        "score_columns": assay_score_columns,
        "n_parameter_sets": int(len(parameter_rows)),
        "n_assay_variants_evaluated": int(len(assay_df)),
        "n_graphlets": int(len(graphlets)),
        "graphlet_dimensions": sorted({int(row["graphlet_dimension"]) for row in catalog_rows}),
        "run077_phot_best_objective": source_phot_best,
    }
    write_json(output_dir / "preflight.json", preflight)
    write_json(
        output_dir / "summary.json",
        {
            "experiment_id": config["experiment"]["id"],
            "preflight": preflight,
            "output_files": {
                "parameter_grid": str(output_dir / "parameter_grid.csv"),
                "assay_fit_metrics": str(output_dir / "assay_fit_metrics.csv"),
                "assay_dfe_metrics": str(output_dir / "assay_dfe_metrics.csv"),
                "assay_epistasis_metrics": str(output_dir / "assay_epistasis_metrics.csv"),
                "assay_fit_rankings": str(output_dir / "assay_fit_rankings.csv"),
                "graphlet_catalog": str(output_dir / "graphlet_catalog.csv"),
                "graphlet_metrics": str(output_dir / "graphlet_metrics.csv"),
                "walsh_terms": str(output_dir / "walsh_terms.csv"),
                "gft_coefficients": str(output_dir / "gft_coefficients.csv"),
                "graphlet_summary": str(output_dir / "graphlet_summary.csv"),
                "model_distances": str(output_dir / "model_distances.csv"),
                "top_model_distances": str(output_dir / "top_model_distances.csv"),
                "local_signal_by_graphlet": str(output_dir / "local_signal_by_graphlet.csv"),
                "local_signal_energy_by_feature": str(output_dir / "local_signal_energy_by_feature.csv"),
                "local_signal_feature_energy_summary": str(output_dir / "local_signal_feature_energy_summary.csv"),
                "local_signal_graphlet_energy_summary": str(output_dir / "local_signal_graphlet_energy_summary.csv"),
                "locality_covariance_summary": str(output_dir / "locality_covariance_summary.csv"),
            },
            "top_assay_fit_rankings": assay_rank_rows[:20],
            "top_graphlet_distances": top_rows[:20],
        },
    )


if __name__ == "__main__":
    main()
