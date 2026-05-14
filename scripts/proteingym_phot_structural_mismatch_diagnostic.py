#!/usr/bin/env python3
"""Run PHOT_CHLRE structural-mismatch diagnostics for HYP-007."""

from __future__ import annotations

import argparse
import itertools
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.stats import ks_2samp, spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from project_adapt_env._compat import ensure_external_paths  # noqa: E402
from project_adapt_env.proteingym_panel import (  # noqa: E402
    load_config,
    merge_dict,
    prepare_proteingym_panel,
)
from project_adapt_env.utils import (  # noqa: E402
    ProgressTracker,
    atomic_write_dataframe_csv,
    atomic_write_json,
    build_logger,
    hash_mod,
    to_builtin,
)

ensure_external_paths()

from adaptenv.alignment import AlignmentProfile  # noqa: E402
from adaptenv.amino_acids import AA_TO_IDX  # noqa: E402


@dataclass(frozen=True)
class DiagnosticData:
    frame: pd.DataFrame
    wildtype: str
    wt_score: float
    alignment_profile: AlignmentProfile
    train_mask: np.ndarray
    single_holdout_mask: np.ndarray
    double_holdout_mask: np.ndarray


@dataclass
class OneDimensionalMonotone:
    centers: np.ndarray
    values: np.ndarray

    def predict(self, x: np.ndarray) -> np.ndarray:
        return np.interp(
            np.asarray(x, dtype=np.float64),
            self.centers,
            self.values,
            left=float(self.values[0]),
            right=float(self.values[-1]),
        )


@dataclass
class TwoDimensionalMonotone:
    x_edges: np.ndarray
    y_edges: np.ndarray
    grid: np.ndarray

    def predict(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        x_bin = np.searchsorted(self.x_edges, np.asarray(x), side="right") - 1
        y_bin = np.searchsorted(self.y_edges, np.asarray(y), side="right") - 1
        x_bin = np.clip(x_bin, 0, self.grid.shape[0] - 1)
        y_bin = np.clip(y_bin, 0, self.grid.shape[1] - 1)
        return self.grid[x_bin, y_bin]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def prediction_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, Any]:
    obs = np.asarray(observed, dtype=np.float64)
    pred = np.asarray(predicted, dtype=np.float64)
    mask = np.isfinite(obs) & np.isfinite(pred)
    obs = obs[mask]
    pred = pred[mask]
    if obs.size == 0:
        return {
            "count": 0,
            "rmse": None,
            "nrmse": None,
            "pearson": None,
            "spearman": None,
        }
    rmse = float(np.sqrt(np.mean((pred - obs) ** 2)))
    scale = max(float(np.std(obs)), 1e-8)
    pearson = (
        float(np.corrcoef(obs, pred)[0, 1])
        if obs.size >= 2 and np.std(pred) > 0
        else None
    )
    spearman = float(spearmanr(obs, pred).statistic) if obs.size >= 2 else None
    return {
        "count": int(obs.size),
        "rmse": rmse,
        "nrmse": float(rmse / scale),
        "pearson": pearson,
        "spearman": spearman,
    }


def mutation_key(wildtype: str, sequence: str) -> tuple[tuple[int, str], ...]:
    return tuple(
        (idx, aa)
        for idx, (wt_aa, aa) in enumerate(zip(wildtype, str(sequence)))
        if wt_aa != aa
    )


def single_mutant_sequence(wildtype: str, mutation: tuple[int, str]) -> str:
    residues = list(wildtype)
    residues[int(mutation[0])] = str(mutation[1])
    return "".join(residues)


def build_mutation_vocab(
    sequences: pd.Series,
    wildtype: str,
) -> dict[tuple[int, str], int]:
    vocab: dict[tuple[int, str], int] = {}
    for sequence in sequences.astype(str):
        for mutation in mutation_key(wildtype, sequence):
            if mutation not in vocab:
                vocab[mutation] = len(vocab)
    return vocab


def build_sparse_mutation_matrix(
    sequences: pd.Series | list[str] | np.ndarray,
    *,
    wildtype: str,
    vocab: dict[tuple[int, str], int],
) -> sparse.csr_matrix:
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    for row_idx, sequence in enumerate(pd.Series(sequences).astype(str)):
        for mutation in mutation_key(wildtype, sequence):
            col = vocab.get(mutation)
            if col is not None:
                rows.append(row_idx)
                cols.append(col)
                values.append(1.0)
    return sparse.csr_matrix(
        (values, (rows, cols)),
        shape=(len(pd.Series(sequences)), len(vocab)),
        dtype=np.float64,
    )


def fit_ridge_additive_score(
    *,
    train_sequences: pd.Series,
    train_scores: np.ndarray,
    all_sequences: pd.Series,
    wildtype: str,
    wt_score: float,
    ridge_alpha: float,
) -> tuple[np.ndarray, dict[tuple[int, str], int], np.ndarray]:
    vocab = build_mutation_vocab(train_sequences, wildtype)
    if not vocab:
        return np.full(len(all_sequences), wt_score), vocab, np.zeros(0)
    x_train = build_sparse_mutation_matrix(
        train_sequences,
        wildtype=wildtype,
        vocab=vocab,
    )
    centered_y = np.asarray(train_scores, dtype=np.float64) - float(wt_score)
    xtx = x_train.T @ x_train
    penalty = sparse.eye(xtx.shape[0], dtype=np.float64) * float(ridge_alpha)
    beta = np.asarray(spsolve(xtx + penalty, x_train.T @ centered_y), dtype=np.float64)
    x_all = build_sparse_mutation_matrix(
        all_sequences,
        wildtype=wildtype,
        vocab=vocab,
    )
    return float(wt_score) + np.asarray(x_all @ beta).reshape(-1), vocab, beta


def score_sequences_from_vocab(
    sequences: list[str],
    *,
    wildtype: str,
    wt_score: float,
    vocab: dict[tuple[int, str], int],
    beta: np.ndarray,
) -> np.ndarray:
    if len(vocab) == 0:
        return np.full(len(sequences), wt_score, dtype=np.float64)
    x = build_sparse_mutation_matrix(sequences, wildtype=wildtype, vocab=vocab)
    return float(wt_score) + np.asarray(x @ beta).reshape(-1)


def stability_proxy(
    sequences: pd.Series | list[str] | np.ndarray,
    *,
    profile: AlignmentProfile,
    wildtype: str,
) -> np.ndarray:
    freq = np.asarray(profile.frequencies, dtype=np.float64)
    wt_value = 0.0
    for pos, aa in enumerate(wildtype):
        wt_value += float(np.log(freq[pos, AA_TO_IDX[aa]] + 1e-9))
    values = []
    for sequence in pd.Series(sequences).astype(str):
        score = 0.0
        for pos, aa in enumerate(sequence):
            idx = AA_TO_IDX.get(aa)
            if idx is None:
                score += np.log(1e-9)
            else:
                score += float(np.log(freq[pos, idx] + 1e-9))
        values.append(score - wt_value)
    return np.asarray(values, dtype=np.float64)


def orient_feature(
    all_values: np.ndarray,
    *,
    train_mask: np.ndarray,
    train_scores: np.ndarray,
) -> tuple[np.ndarray, float]:
    train_values = np.asarray(all_values, dtype=np.float64)[train_mask]
    if train_values.size < 3 or np.std(train_values) <= 0:
        return np.asarray(all_values, dtype=np.float64), 1.0
    corr = spearmanr(train_values, train_scores).statistic
    if np.isfinite(corr) and corr < 0:
        return -np.asarray(all_values, dtype=np.float64), -1.0
    return np.asarray(all_values, dtype=np.float64), 1.0


def fit_1d_monotone(
    x: np.ndarray,
    y: np.ndarray,
    *,
    n_bins: int,
) -> OneDimensionalMonotone:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    order = np.argsort(x, kind="mergesort")
    x = x[order]
    y = y[order]
    n_bins = max(3, min(int(n_bins), int(x.size)))
    edges = np.linspace(0, x.size, n_bins + 1, dtype=int)
    centers = []
    values = []
    for start, stop in zip(edges[:-1], edges[1:]):
        if stop <= start:
            continue
        centers.append(float(np.mean(x[start:stop])))
        values.append(float(np.mean(y[start:stop])))
    centers_arr = np.asarray(centers, dtype=np.float64)
    values_arr = np.maximum.accumulate(np.asarray(values, dtype=np.float64))
    uniq, idx = np.unique(centers_arr, return_index=True)
    return OneDimensionalMonotone(centers=uniq, values=values_arr[idx])


def quantile_edges(values: np.ndarray, n_bins: int) -> np.ndarray:
    probs = np.linspace(0.0, 1.0, max(int(n_bins), 2) + 1)
    edges = np.quantile(np.asarray(values, dtype=np.float64), probs)
    edges[0] = -np.inf
    edges[-1] = np.inf
    for idx in range(1, len(edges) - 1):
        if edges[idx] <= edges[idx - 1]:
            edges[idx] = np.nextafter(edges[idx - 1], np.inf)
    return edges


def fit_2d_monotone(
    x: np.ndarray,
    z: np.ndarray,
    y: np.ndarray,
    *,
    n_bins: int,
) -> TwoDimensionalMonotone:
    n_bins = max(int(n_bins), 3)
    x_edges = quantile_edges(x, n_bins)
    z_edges = quantile_edges(z, n_bins)
    x_bin = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, n_bins - 1)
    z_bin = np.clip(np.searchsorted(z_edges, z, side="right") - 1, 0, n_bins - 1)
    grid = np.full((n_bins, n_bins), np.nan, dtype=np.float64)
    global_mean = float(np.mean(y))
    for i in range(n_bins):
        for j in range(n_bins):
            mask = (x_bin == i) & (z_bin == j)
            if np.any(mask):
                grid[i, j] = float(np.mean(y[mask]))
    grid = fill_missing_grid(grid, global_mean)
    for _ in range(4):
        grid = np.maximum.accumulate(grid, axis=0)
        grid = np.maximum.accumulate(grid, axis=1)
    return TwoDimensionalMonotone(x_edges=x_edges, y_edges=z_edges, grid=grid)


def fill_missing_grid(grid: np.ndarray, default: float) -> np.ndarray:
    filled = grid.copy()
    coords = np.argwhere(np.isfinite(filled))
    if coords.size == 0:
        return np.full_like(filled, default)
    for i in range(filled.shape[0]):
        for j in range(filled.shape[1]):
            if np.isfinite(filled[i, j]):
                continue
            distances = np.abs(coords[:, 0] - i) + np.abs(coords[:, 1] - j)
            nearest = coords[int(np.argmin(distances))]
            filled[i, j] = filled[nearest[0], nearest[1]]
    return filled


def fit_pair_residuals(
    frame: pd.DataFrame,
    *,
    wildtype: str,
    train_mask: np.ndarray,
    residuals: np.ndarray,
    shrinkage: float,
) -> dict[tuple[tuple[int, str], tuple[int, str]], float]:
    buckets: dict[tuple[tuple[int, str], tuple[int, str]], list[float]] = {}
    train_frame = frame.loc[train_mask].copy()
    train_residuals = np.asarray(residuals, dtype=np.float64)[train_mask]
    for seq, mut_count, residual in zip(
        train_frame["mutated_sequence"].astype(str),
        train_frame["mutation_count"].astype(int),
        train_residuals,
    ):
        if int(mut_count) != 2:
            continue
        mutations = mutation_key(wildtype, seq)
        if len(mutations) != 2:
            continue
        key = tuple(sorted(mutations))
        buckets.setdefault(key, []).append(float(residual))
    effects = {}
    for key, values in buckets.items():
        n = len(values)
        effects[key] = float(n / (n + float(shrinkage)) * np.mean(values))
    return effects


def predict_pair_residuals(
    sequences: pd.Series | list[str] | np.ndarray,
    *,
    wildtype: str,
    effects: dict[tuple[tuple[int, str], tuple[int, str]], float],
) -> np.ndarray:
    out = []
    for sequence in pd.Series(sequences).astype(str):
        mutations = mutation_key(wildtype, sequence)
        total = 0.0
        for pair in itertools.combinations(mutations, 2):
            total += effects.get(tuple(sorted(pair)), 0.0)
        out.append(total)
    return np.asarray(out, dtype=np.float64)


def epistasis_metrics(
    frame: pd.DataFrame,
    *,
    wildtype: str,
    wt_score: float,
    predictions_by_sequence: dict[str, float],
    eval_mask: np.ndarray,
    max_pairs: int,
) -> dict[str, Any]:
    score_by_sequence = {
        str(row.mutated_sequence): float(row.DMS_score)
        for row in frame.itertuples(index=False)
    }
    observed = []
    predicted = []
    eval_frame = frame.loc[eval_mask].copy()
    if max_pairs > 0 and len(eval_frame) > max_pairs:
        eval_frame = eval_frame.sample(n=max_pairs, random_state=91)
    for row in eval_frame.itertuples(index=False):
        mutations = mutation_key(wildtype, str(row.mutated_sequence))
        if len(mutations) != 2:
            continue
        single_a = single_mutant_sequence(wildtype, mutations[0])
        single_b = single_mutant_sequence(wildtype, mutations[1])
        if single_a not in score_by_sequence or single_b not in score_by_sequence:
            continue
        if (
            str(row.mutated_sequence) not in predictions_by_sequence
            or single_a not in predictions_by_sequence
            or single_b not in predictions_by_sequence
        ):
            continue
        observed.append(
            float(row.DMS_score)
            - score_by_sequence[single_a]
            - score_by_sequence[single_b]
            + float(wt_score)
        )
        predicted.append(
            predictions_by_sequence[str(row.mutated_sequence)]
            - predictions_by_sequence[single_a]
            - predictions_by_sequence[single_b]
            + float(wt_score)
        )
    if not observed:
        return {"count": 0, "spearman": None, "ks": None}
    obs = np.asarray(observed, dtype=np.float64)
    pred = np.asarray(predicted, dtype=np.float64)
    return {
        "count": int(obs.size),
        "spearman": float(spearmanr(obs, pred).statistic) if obs.size >= 2 else None,
        "ks": float(ks_2samp(obs, pred).statistic),
    }


def evaluate_model(
    data: DiagnosticData,
    *,
    model_name: str,
    predictions: np.ndarray,
    max_epistasis_pairs: int,
) -> dict[str, Any]:
    frame = data.frame
    pred_by_sequence = {
        str(sequence): float(pred)
        for sequence, pred in zip(frame["mutated_sequence"].astype(str), predictions)
    }
    single = data.single_holdout_mask
    double = data.double_holdout_mask
    return {
        "model": model_name,
        "single_holdout": prediction_metrics(
            frame.loc[single, "DMS_score"].to_numpy(dtype=np.float64),
            np.asarray(predictions, dtype=np.float64)[single],
        ),
        "double_holdout": prediction_metrics(
            frame.loc[double, "DMS_score"].to_numpy(dtype=np.float64),
            np.asarray(predictions, dtype=np.float64)[double],
        ),
        "epistasis_prediction": epistasis_metrics(
            frame,
            wildtype=data.wildtype,
            wt_score=data.wt_score,
            predictions_by_sequence=pred_by_sequence,
            eval_mask=double,
            max_pairs=max_epistasis_pairs,
        ),
    }


def build_data(
    *,
    config: dict[str, Any],
    project_root: Path,
    output_dir: Path,
    progress: ProgressTracker,
    logger: Any,
) -> DiagnosticData:
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    proteingym_cfg = dict(config["proteingym"])
    panel_cfg = dict(config["panel"])
    mmseqs_cfg = dict(config["mmseqs"])
    mmseqs_cfg["cache_dir"] = str(project_root / mmseqs_cfg["cache_dir"])
    mavenn_cfg = dict(config["mavenn"])

    started = time.time()
    progress.write(
        stage="panel_preparation",
        completed=0,
        total=1,
        message="Starting PHOT panel preparation.",
        stage_started_at=started,
    )
    panel = prepare_proteingym_panel(
        project_root=project_root,
        proteingym_cfg=proteingym_cfg,
        panel_cfg=panel_cfg,
        mmseqs_cfg=mmseqs_cfg,
        mavenn_cfg=mavenn_cfg,
        calibration_max_mutation_count=int(config["calibration_max_mutation_count"]),
        checkpoint_dir=checkpoint_dir,
        progress_callback=lambda event: progress.write(
            stage="panel_preparation",
            completed=int(event["completed"]),
            total=int(event["total"]),
            message=f"Prepared assay {event['dms_id']}.",
            details=event,
            stage_started_at=started,
        ),
    )
    if len(panel.assays) != 1:
        raise ValueError("This diagnostic currently expects exactly one PHOT assay")
    assay = panel.assays[0]
    atomic_write_dataframe_csv(output_dir / "selected_panel.csv", panel.panel_df, index=False)
    atomic_write_dataframe_csv(
        output_dir / "mavenn_assay_metrics.csv",
        panel.mavenn_metrics_frame(),
        index=False,
    )
    cache_frame = checkpoint_dir / "panel_preparation" / assay.dms_id / "prepared_assay.csv"
    frame = pd.read_csv(cache_frame if cache_frame.is_file() else assay.assay_csv_path)
    frame = frame[frame["mutation_count"] >= 1].copy().reset_index(drop=True)

    split_cfg = dict(config["diagnostic"].get("splits", {}))
    single_mod = int(split_cfg.get("single_holdout_modulo", 5))
    single_rem = int(split_cfg.get("single_holdout_remainder", 0))
    double_mod = int(split_cfg.get("double_holdout_modulo", 5))
    double_rem = int(split_cfg.get("double_holdout_remainder", 0))
    seqs = frame["mutated_sequence"].astype(str)
    single_holdout = (
        (frame["mutation_count"].astype(int) == 1).to_numpy()
        & np.asarray([hash_mod(seq, single_mod) == single_rem for seq in seqs], dtype=bool)
    )
    double_holdout = (
        (frame["mutation_count"].astype(int) == 2).to_numpy()
        & np.asarray([hash_mod(seq, double_mod) == double_rem for seq in seqs], dtype=bool)
    )
    max_train_mutation_count = int(config["diagnostic"].get("max_train_mutation_count", 5))
    train_mask = (
        (frame["mutation_count"].astype(int).to_numpy() >= 1)
        & (frame["mutation_count"].astype(int).to_numpy() <= max_train_mutation_count)
        & ~single_holdout
        & ~double_holdout
    )
    max_train_variants = int(config["diagnostic"].get("max_train_variants", 0))
    if max_train_variants > 0 and int(np.sum(train_mask)) > max_train_variants:
        train_indices = np.flatnonzero(train_mask)
        rng = np.random.default_rng(int(config["diagnostic"].get("sample_seed", 17)))
        keep = rng.choice(train_indices, size=max_train_variants, replace=False)
        reduced = np.zeros_like(train_mask, dtype=bool)
        reduced[keep] = True
        train_mask = reduced
        logger.info("Downsampled training variants to %d.", max_train_variants)

    profile = AlignmentProfile(assay.alignment_path)
    return DiagnosticData(
        frame=frame,
        wildtype=assay.wildtype_sequence,
        wt_score=float(assay.mavenn_metrics["wildtype_yhat"]),
        alignment_profile=profile,
        train_mask=train_mask,
        single_holdout_mask=single_holdout,
        double_holdout_mask=double_holdout,
    )


def run_oracle_global_epistasis(data: DiagnosticData, cfg: dict[str, Any]) -> list[dict[str, Any]]:
    frame = data.frame
    y = frame["DMS_score"].to_numpy(dtype=np.float64)
    train_scores = y[data.train_mask]
    phi, _, _ = fit_ridge_additive_score(
        train_sequences=frame.loc[data.train_mask, "mutated_sequence"],
        train_scores=train_scores,
        all_sequences=frame["mutated_sequence"],
        wildtype=data.wildtype,
        wt_score=data.wt_score,
        ridge_alpha=float(cfg.get("ridge_alpha", 2.0)),
    )
    phi, _ = orient_feature(phi, train_mask=data.train_mask, train_scores=train_scores)
    mono = fit_1d_monotone(
        phi[data.train_mask],
        train_scores,
        n_bins=int(cfg.get("n_bins", 32)),
    )
    global_epistasis_pred = mono.predict(phi)
    return [
        evaluate_model(
            data,
            model_name="oracle_additive_linear",
            predictions=phi,
            max_epistasis_pairs=int(cfg.get("max_epistasis_pairs", 20000)),
        ),
        evaluate_model(
            data,
            model_name="oracle_additive_monotone_global_epistasis",
            predictions=global_epistasis_pred,
            max_epistasis_pairs=int(cfg.get("max_epistasis_pairs", 20000)),
        ),
    ]


def run_flexible_two_latent(data: DiagnosticData, cfg: dict[str, Any]) -> list[dict[str, Any]]:
    frame = data.frame
    y = frame["DMS_score"].to_numpy(dtype=np.float64)
    train_scores = y[data.train_mask]
    stability = stability_proxy(
        frame["mutated_sequence"],
        profile=data.alignment_profile,
        wildtype=data.wildtype,
    )
    activity, _, _ = fit_ridge_additive_score(
        train_sequences=frame.loc[data.train_mask, "mutated_sequence"],
        train_scores=train_scores,
        all_sequences=frame["mutated_sequence"],
        wildtype=data.wildtype,
        wt_score=data.wt_score,
        ridge_alpha=float(cfg.get("ridge_alpha", 2.0)),
    )
    stability, _ = orient_feature(
        stability,
        train_mask=data.train_mask,
        train_scores=train_scores,
    )
    activity, _ = orient_feature(
        activity,
        train_mask=data.train_mask,
        train_scores=train_scores,
    )
    surface = fit_2d_monotone(
        stability[data.train_mask],
        activity[data.train_mask],
        train_scores,
        n_bins=int(cfg.get("n_bins_2d", 14)),
    )
    predictions = surface.predict(stability, activity)
    return [
        evaluate_model(
            data,
            model_name="flexible_monotone_stability_activity_surface",
            predictions=predictions,
            max_epistasis_pairs=int(cfg.get("max_epistasis_pairs", 20000)),
        )
    ]


def run_photophysical_sparse_residual(data: DiagnosticData, cfg: dict[str, Any]) -> list[dict[str, Any]]:
    frame = data.frame
    y = frame["DMS_score"].to_numpy(dtype=np.float64)
    train_scores = y[data.train_mask]
    stability = stability_proxy(
        frame["mutated_sequence"],
        profile=data.alignment_profile,
        wildtype=data.wildtype,
    )
    stability, _ = orient_feature(
        stability,
        train_mask=data.train_mask,
        train_scores=train_scores,
    )
    stability_model = fit_1d_monotone(
        stability[data.train_mask],
        train_scores,
        n_bins=int(cfg.get("n_bins", 32)),
    )
    stability_pred = stability_model.predict(stability)
    residual_train = y[data.train_mask] - stability_pred[data.train_mask]
    photophysical_raw, vocab, beta = fit_ridge_additive_score(
        train_sequences=frame.loc[data.train_mask, "mutated_sequence"],
        train_scores=residual_train + data.wt_score,
        all_sequences=frame["mutated_sequence"],
        wildtype=data.wildtype,
        wt_score=data.wt_score,
        ridge_alpha=float(cfg.get("ridge_alpha", 2.0)),
    )
    photophysical = photophysical_raw - data.wt_score
    photophysical, _ = orient_feature(
        photophysical,
        train_mask=data.train_mask,
        train_scores=residual_train,
    )
    surface = fit_2d_monotone(
        stability[data.train_mask],
        photophysical[data.train_mask],
        train_scores,
        n_bins=int(cfg.get("n_bins_2d", 14)),
    )
    two_latent_pred = surface.predict(stability, photophysical)
    residual_after_two = y - two_latent_pred
    pair_effects = fit_pair_residuals(
        frame,
        wildtype=data.wildtype,
        train_mask=data.train_mask,
        residuals=residual_after_two,
        shrinkage=float(cfg.get("pair_residual_shrinkage", 4.0)),
    )
    sparse_pred = two_latent_pred + predict_pair_residuals(
        frame["mutated_sequence"],
        wildtype=data.wildtype,
        effects=pair_effects,
    )
    return [
        evaluate_model(
            data,
            model_name="stability_first_monotone",
            predictions=stability_pred,
            max_epistasis_pairs=int(cfg.get("max_epistasis_pairs", 20000)),
        ),
        evaluate_model(
            data,
            model_name="stability_plus_photophysical_residual_trait",
            predictions=two_latent_pred,
            max_epistasis_pairs=int(cfg.get("max_epistasis_pairs", 20000)),
        ),
        {
            **evaluate_model(
                data,
                model_name="stability_plus_photophysical_trait_sparse_pair_residual",
                predictions=sparse_pred,
                max_epistasis_pairs=int(cfg.get("max_epistasis_pairs", 20000)),
            ),
            "sparse_pair_effect_count": int(len(pair_effects)),
            "photophysical_mutation_feature_count": int(len(vocab)),
            "photophysical_beta_nonzero_count": int(np.sum(np.abs(beta) > 1e-9)),
        },
    ]


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config))
    if args.quick and "quick" in config:
        config = merge_dict(config, config["quick"])
    experiment_id = str(config.get("experiment_id") or config.get("experiment", {}).get("id"))
    logger = build_logger(
        "project_adapt_env.phot_structural_mismatch",
        output_dir / "progress.log",
    )
    progress = ProgressTracker(
        progress_path=output_dir / "progress.json",
        logger=logger,
        run_label=experiment_id,
    )
    data = build_data(
        config=config,
        project_root=project_root,
        output_dir=output_dir,
        progress=progress,
        logger=logger,
    )
    mode = str(config["diagnostic"]["mode"])
    mode_cfg = dict(config["diagnostic"].get("model", {}))
    progress.write(
        stage="diagnostic_fit",
        completed=0,
        total=1,
        message=f"Running diagnostic mode {mode}.",
        stage_started_at=time.time(),
    )
    runners: dict[str, Callable[[DiagnosticData, dict[str, Any]], list[dict[str, Any]]]] = {
        "flexible_monotone_two_latent": run_flexible_two_latent,
        "photophysical_sparse_residual": run_photophysical_sparse_residual,
        "oracle_global_epistasis": run_oracle_global_epistasis,
    }
    if mode not in runners:
        raise ValueError(f"Unsupported diagnostic mode: {mode}")
    model_rows = runners[mode](data, mode_cfg)
    flat_rows = []
    for row in model_rows:
        flat = {"model": row["model"]}
        for key, value in row.items():
            if key == "model":
                continue
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    flat[f"{key}_{sub_key}"] = sub_value
            else:
                flat[key] = value
        flat_rows.append(flat)
    atomic_write_dataframe_csv(
        output_dir / "model_validations.csv",
        pd.DataFrame(flat_rows),
        index=False,
    )
    summary = {
        "experiment_id": experiment_id,
        "diagnostic_mode": mode,
        "assay_id": str(config["panel"]["assay_ids"][0]),
        "train_variant_count": int(np.sum(data.train_mask)),
        "single_holdout_count": int(np.sum(data.single_holdout_mask)),
        "double_holdout_count": int(np.sum(data.double_holdout_mask)),
        "models": to_builtin(model_rows),
        "config": {
            "max_train_mutation_count": int(config["diagnostic"].get("max_train_mutation_count", 5)),
            "max_train_variants": int(config["diagnostic"].get("max_train_variants", 0)),
        },
    }
    atomic_write_json(output_dir / "summary.json", summary)
    progress.write(
        stage="run_complete",
        completed=1,
        total=1,
        message=f"Run completed; summary written to {output_dir / 'summary.json'}.",
    )
    print(output_dir / "summary.json")


if __name__ == "__main__":
    main()
