#!/usr/bin/env python3
"""Run EXP-020: adaptive-walk and peak-proximal epistasis baselines."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
if str(ADAPT_ENV_ROOT) not in sys.path:
    sys.path.insert(0, str(ADAPT_ENV_ROOT))

from adaptenv import AlignmentProfile, FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from adaptenv.amino_acids import AA_ORDER  # noqa: E402


LITERATURE_BASELINES = [
    {
        "id": "declining_adaptability",
        "expectation": (
            "Higher-fitness genetic backgrounds should adapt more slowly, with "
            "smaller available beneficial effects than lower-fitness backgrounds."
        ),
        "sources": [
            "Kryazhimskiy et al. 2014, Science, doi:10.1126/science.1250939",
            "Wiser et al. 2013, Science, doi:10.1126/science.1243357",
        ],
        "metrics": [
            "start_fraction_of_peak versus adaptive-walk total gain",
            "background fraction_of_peak versus beneficial fraction",
            "background fraction_of_peak versus mean beneficial effect",
        ],
    },
    {
        "id": "decelerating_adaptive_trajectories",
        "expectation": (
            "Adaptive trajectories should show concave or decelerating fitness "
            "gain rather than unlimited constant gains near the peak."
        ),
        "sources": [
            "Wiser et al. 2013, Science, doi:10.1126/science.1243357",
            "Kryazhimskiy et al. 2014, Science, doi:10.1126/science.1250939",
        ],
        "metrics": [
            "late adaptive-walk gain divided by early adaptive-walk gain",
            "walk length and stalling probability near high-fitness backgrounds",
        ],
    },
    {
        "id": "beneficial_dfe_near_optimum",
        "expectation": (
            "Near an optimum, the distribution of beneficial effects should be "
            "bounded and depleted relative to more distant backgrounds."
        ),
        "sources": [
            "Martin and Lenormand 2008, Genetics, doi:10.1534/genetics.108.087122",
        ],
        "metrics": [
            "single-mutant beneficial fraction by fraction_of_peak",
            "single-mutant best and mean beneficial effect by fraction_of_peak",
        ],
    },
    {
        "id": "ruggedness_and_sign_epistasis",
        "expectation": (
            "Sign epistasis and ruggedness should constrain accessible paths; "
            "peak-proximal neighborhoods should show stronger or more "
            "constraining pairwise interactions than distant low-fitness "
            "neighborhoods."
        ),
        "sources": [
            "de Visser and Krug 2014, Nature Reviews Genetics, doi:10.1038/nrg3744",
            "Weinreich et al. 2006, Science, doi:10.1126/science.1123539",
        ],
        "metrics": [
            "mean absolute pairwise epistasis by fraction_of_peak",
            "negative, diminishing-returns, and sign-epistasis fractions",
        ],
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if quick:
        merged = merge_dict(raw, raw.get("quick", {}))
        merged.pop("quick", None)
        return merged
    raw.pop("quick", None)
    return raw


def merge_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_dict(dict(merged[key]), value)
        else:
            merged[key] = value
    return merged


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return str(value)


def write_json(path: Path, payload: dict[str, Any] | list[Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def finite_mean(values: list[float] | np.ndarray) -> float | None:
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    return float(np.mean(arr))


def rankdata(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    order = np.argsort(arr, kind="mergesort")
    ranks = np.empty(arr.shape[0], dtype=np.float64)
    sorted_arr = arr[order]
    start = 0
    while start < sorted_arr.size:
        stop = start + 1
        while stop < sorted_arr.size and sorted_arr[stop] == sorted_arr[start]:
            stop += 1
        rank = 0.5 * (start + stop - 1) + 1.0
        ranks[order[start:stop]] = rank
        start = stop
    return ranks


def pearson(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 3:
        return None
    x = x[mask]
    y = y[mask]
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 3:
        return None
    return pearson(rankdata(x[mask]), rankdata(y[mask]))


def sign_with_threshold(value: float, threshold: float) -> int:
    if value > threshold:
        return 1
    if value < -threshold:
        return -1
    return 0


def threshold_for(parent_fitness: float, cfg: dict[str, Any], prefix: str) -> float:
    abs_key = f"{prefix}_abs"
    rel_key = f"{prefix}_fraction_of_parent"
    absolute = float(cfg.get(abs_key, 0.0))
    relative = float(cfg.get(rel_key, 0.0)) * abs(float(parent_fitness))
    return max(absolute, relative)


def sequence_to_string(seq: np.ndarray) -> str:
    return "".join(AA_ORDER[int(aa)] for aa in np.asarray(seq, dtype=np.int64))


def generate_synthetic_alignment(
    seed_sequence: np.ndarray,
    *,
    n_sequences: int,
    seed: int,
    n_coupled_pairs: int | None = None,
) -> list[str]:
    rng = np.random.RandomState(seed)
    seed_idx = np.asarray(seed_sequence, dtype=np.int64)
    length = int(seed_idx.size)
    conservation = rng.beta(2.8, 2.0, size=length)
    if n_coupled_pairs is None:
        n_coupled_pairs = max(4, length // 12)

    candidate_pairs = [(i, j) for i in range(length) for j in range(i + 4, length)]
    rng.shuffle(candidate_pairs)
    coupled_pairs = candidate_pairs[:n_coupled_pairs]
    paired_states: list[list[tuple[int, int]]] = []
    for i, j in coupled_pairs:
        states: list[tuple[int, int]] = []
        for _ in range(int(rng.randint(2, 4))):
            aa_i = int(rng.randint(0, 20))
            aa_j = int(rng.randint(0, 20))
            while aa_i == int(seed_idx[i]):
                aa_i = int(rng.randint(0, 20))
            while aa_j == int(seed_idx[j]):
                aa_j = int(rng.randint(0, 20))
            states.append((aa_i, aa_j))
        paired_states.append(states)

    sequences = [sequence_to_string(seed_idx)]
    independent_rate = 0.45 * (1.0 - conservation)
    for _ in range(max(0, int(n_sequences) - 1)):
        seq = seed_idx.copy()
        independent_mask = rng.random_sample(length) < independent_rate
        for pos in np.flatnonzero(independent_mask):
            aa = int(rng.randint(0, 19))
            if aa >= int(seq[pos]):
                aa += 1
            seq[pos] = aa

        for pair_idx, (i, j) in enumerate(coupled_pairs):
            pair_conservation = 0.5 * (conservation[i] + conservation[j])
            if rng.random_sample() < 0.30 * (1.0 - 0.5 * pair_conservation):
                aa_i, aa_j = paired_states[pair_idx][
                    int(rng.randint(len(paired_states[pair_idx])))
                ]
                seq[i] = aa_i
                seq[j] = aa_j
        sequences.append(sequence_to_string(seq))
    return sequences


def build_landscape(
    cfg: dict[str, Any],
    *,
    seed: int,
    condition: dict[str, Any],
) -> FitnessLandscapeEnv:
    length = int(cfg["landscape"]["sequence_length"])
    base_config = dict(cfg["adapt_env"]["base_config"])
    base_config["L"] = length
    base_config["seed"] = int(seed)
    landscape_config = LandscapeConfig(**base_config)

    alignment_mode = str(condition.get("alignment", "none")).lower()
    if alignment_mode in {"none", "false", "unconditioned"}:
        return FitnessLandscapeEnv(landscape_config)
    if alignment_mode != "synthetic":
        raise ValueError(f"Unsupported alignment mode: {alignment_mode}")

    rng = np.random.RandomState(int(seed) + int(condition.get("alignment_seed_offset", 10_000)))
    reference = rng.randint(0, 20, size=length, dtype=np.int64)
    sequences = generate_synthetic_alignment(
        reference,
        n_sequences=int(condition.get("alignment_size", 320)),
        seed=int(seed) + int(condition.get("alignment_generation_seed_offset", 20_000)),
        n_coupled_pairs=condition.get("n_coupled_pairs"),
    )
    profile = AlignmentProfile.from_sequences(sequences)
    return FitnessLandscapeEnv(landscape_config, alignment_profile=profile)


def enumerate_single_mutants(
    seq: np.ndarray,
    alphabet_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    seq = np.asarray(seq, dtype=np.int64)
    length = int(seq.shape[0])
    mutant_count = length * (alphabet_size - 1)
    mutants = np.repeat(seq[None, :], mutant_count, axis=0)
    positions = np.repeat(np.arange(length, dtype=np.int64), alphabet_size - 1)
    amino_acids = np.empty(mutant_count, dtype=np.int64)
    offset = 0
    for pos in range(length):
        alts = np.delete(np.arange(alphabet_size, dtype=np.int64), seq[pos])
        mutants[offset:offset + alts.size, pos] = alts
        amino_acids[offset:offset + alts.size] = alts
        offset += alts.size
    return mutants, positions, amino_acids


def evaluate_batch(
    landscape: FitnessLandscapeEnv,
    seqs: np.ndarray,
    *,
    batch_size: int,
) -> np.ndarray:
    parts: list[np.ndarray] = []
    for start in range(0, seqs.shape[0], batch_size):
        stop = min(start + batch_size, seqs.shape[0])
        parts.append(np.asarray(landscape.evaluate_batch(seqs[start:stop]), dtype=np.float64))
    return np.concatenate(parts) if parts else np.empty(0, dtype=np.float64)


def single_mutant_scan(
    landscape: FitnessLandscapeEnv,
    seq: np.ndarray,
    *,
    sampling_cfg: dict[str, Any],
) -> dict[str, Any]:
    seq = np.asarray(seq, dtype=np.int64)
    parent_fitness = float(landscape.evaluate(seq))
    threshold = threshold_for(parent_fitness, sampling_cfg, "beneficial_threshold")
    neutral_threshold = threshold_for(parent_fitness, sampling_cfg, "neutral_threshold")
    threshold = max(threshold, neutral_threshold)

    mutants, positions, amino_acids = enumerate_single_mutants(seq, landscape.A)
    fitnesses = evaluate_batch(
        landscape,
        mutants,
        batch_size=int(sampling_cfg["batch_size"]),
    )
    deltas = fitnesses - parent_fitness
    beneficial = deltas > threshold
    neutral = np.abs(deltas) <= neutral_threshold
    deleterious = deltas < -threshold
    return {
        "parent_fitness": parent_fitness,
        "threshold": float(threshold),
        "neutral_threshold": float(neutral_threshold),
        "mutants": mutants,
        "positions": positions,
        "amino_acids": amino_acids,
        "fitnesses": fitnesses,
        "deltas": deltas,
        "beneficial": beneficial,
        "neutral": neutral,
        "deleterious": deleterious,
        "fraction_beneficial": float(np.mean(beneficial)),
        "fraction_neutral": float(np.mean(neutral)),
        "fraction_deleterious": float(np.mean(deleterious)),
        "mean_delta": float(np.mean(deltas)),
        "best_delta": float(np.max(deltas)) if deltas.size else 0.0,
        "mean_beneficial_effect": float(np.mean(deltas[beneficial])) if np.any(beneficial) else 0.0,
    }


def sequence_at_distance_from_peak(
    landscape: FitnessLandscapeEnv,
    distance: int,
    rng: np.random.RandomState,
) -> np.ndarray:
    distance = int(max(0, min(distance, landscape.L)))
    seq = landscape.peak_sequence.copy()
    if distance == 0:
        return seq
    positions = rng.choice(landscape.L, size=distance, replace=False)
    for pos in positions:
        choices = np.delete(np.arange(landscape.A, dtype=np.int64), seq[pos])
        seq[pos] = int(choices[int(rng.randint(choices.size))])
    return seq


def sequence_near_fraction_of_peak(
    landscape: FitnessLandscapeEnv,
    target_fraction: float,
    rng: np.random.RandomState,
    *,
    n_candidates: int,
) -> tuple[np.ndarray, float]:
    peak_fitness = float(landscape.peak_fitness)
    target = float(target_fraction) * peak_fitness
    candidates: list[np.ndarray] = [landscape.reference.copy(), landscape.peak_sequence.copy()]
    max_peak_radius = max(1, int(round((1.05 - min(float(target_fraction), 1.0)) * landscape.L)))
    max_peak_radius = min(landscape.L, max_peak_radius + 4)
    for _ in range(int(n_candidates)):
        if target_fraction < 0.35 and rng.random_sample() < 0.45:
            candidates.append(landscape.random_sequence())
        else:
            distance = int(rng.randint(0, max_peak_radius + 1))
            candidates.append(sequence_at_distance_from_peak(landscape, distance, rng))
    batch = np.asarray(candidates, dtype=np.int64)
    fitnesses = evaluate_batch(landscape, batch, batch_size=4096)
    idx = int(np.argmin(np.abs(fitnesses - target)))
    return batch[idx].copy(), float(fitnesses[idx])


def run_adaptive_walk(
    landscape: FitnessLandscapeEnv,
    start_seq: np.ndarray,
    *,
    mode: str,
    rng: np.random.RandomState,
    sampling_cfg: dict[str, Any],
) -> dict[str, Any]:
    seq = np.asarray(start_seq, dtype=np.int64).copy()
    max_steps = int(sampling_cfg["max_walk_steps"])
    trajectory: list[dict[str, Any]] = []
    previous_fitness: float | None = None
    stop_reason = "max_steps"

    for step in range(max_steps + 1):
        scan = single_mutant_scan(landscape, seq, sampling_cfg=sampling_cfg)
        current_fitness = float(scan["parent_fitness"])
        beneficial_idx = np.flatnonzero(scan["beneficial"])
        row = {
            "step": int(step),
            "fitness": current_fitness,
            "fraction_of_peak": (
                current_fitness / float(landscape.peak_fitness)
                if abs(float(landscape.peak_fitness)) > 1e-12
                else math.nan
            ),
            "distance_to_peak": int(landscape.hamming_distance(seq, landscape.peak_sequence)),
            "step_gain": (
                0.0 if previous_fitness is None else current_fitness - previous_fitness
            ),
            "available_beneficial_fraction": float(scan["fraction_beneficial"]),
            "available_mean_beneficial_effect": float(scan["mean_beneficial_effect"]),
            "best_available_effect": float(scan["best_delta"]),
            "chosen_effect": math.nan,
        }
        trajectory.append(row)
        if step == max_steps:
            break
        if beneficial_idx.size == 0:
            stop_reason = "no_beneficial_single_mutant"
            break
        if mode == "greedy":
            chosen = int(beneficial_idx[np.argmax(scan["deltas"][beneficial_idx])])
        elif mode == "random_beneficial":
            chosen = int(beneficial_idx[int(rng.randint(beneficial_idx.size))])
        else:
            raise ValueError(f"Unsupported adaptive walk mode: {mode}")
        trajectory[-1]["chosen_effect"] = float(scan["deltas"][chosen])
        seq = scan["mutants"][chosen].copy()
        previous_fitness = current_fitness

    fitnesses = np.asarray([row["fitness"] for row in trajectory], dtype=np.float64)
    if fitnesses.size >= 4:
        midpoint = max(1, int((fitnesses.size - 1) // 2))
        early_gain = float(fitnesses[midpoint] - fitnesses[0])
        late_gain = float(fitnesses[-1] - fitnesses[midpoint])
        deceleration_ratio = late_gain / (early_gain + 1e-12) if early_gain > 0 else math.nan
    else:
        early_gain = math.nan
        late_gain = math.nan
        deceleration_ratio = math.nan
    return {
        "trajectory": trajectory,
        "final_sequence": seq,
        "summary": {
            "mode": mode,
            "n_steps": int(max(0, len(trajectory) - 1)),
            "start_fitness": float(fitnesses[0]),
            "final_fitness": float(fitnesses[-1]),
            "total_gain": float(fitnesses[-1] - fitnesses[0]),
            "early_gain": early_gain,
            "late_gain": late_gain,
            "deceleration_ratio": deceleration_ratio,
            "stop_reason": stop_reason,
        },
    }


def sample_epistasis_pairs(
    landscape: FitnessLandscapeEnv,
    background_seq: np.ndarray,
    *,
    background_id: str,
    row_context: dict[str, Any],
    rng: np.random.RandomState,
    sampling_cfg: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n_pairs = int(sampling_cfg["epistasis_pairs_per_background"])
    batch_size = int(sampling_cfg["batch_size"])
    background_seq = np.asarray(background_seq, dtype=np.int64)
    background_fitness = float(landscape.evaluate(background_seq))
    threshold = max(
        threshold_for(background_fitness, sampling_cfg, "beneficial_threshold"),
        threshold_for(background_fitness, sampling_cfg, "neutral_threshold"),
        float(sampling_cfg.get("min_epistasis_effect_threshold", 0.0)),
    )

    pos_i = rng.randint(0, landscape.L, size=n_pairs)
    pos_j = rng.randint(0, landscape.L - 1, size=n_pairs)
    pos_j = pos_j + (pos_j >= pos_i)
    aa_i = rng.randint(0, landscape.A - 1, size=n_pairs)
    aa_j = rng.randint(0, landscape.A - 1, size=n_pairs)
    aa_i = aa_i + (aa_i >= background_seq[pos_i])
    aa_j = aa_j + (aa_j >= background_seq[pos_j])

    single_i = np.repeat(background_seq[None, :], n_pairs, axis=0)
    single_j = np.repeat(background_seq[None, :], n_pairs, axis=0)
    double = np.repeat(background_seq[None, :], n_pairs, axis=0)
    rows_idx = np.arange(n_pairs, dtype=np.int64)
    single_i[rows_idx, pos_i] = aa_i
    single_j[rows_idx, pos_j] = aa_j
    double[rows_idx, pos_i] = aa_i
    double[rows_idx, pos_j] = aa_j

    stacked = np.vstack([single_i, single_j, double])
    fitness = evaluate_batch(landscape, stacked, batch_size=batch_size)
    fit_i = fitness[:n_pairs]
    fit_j = fitness[n_pairs:2 * n_pairs]
    fit_double = fitness[2 * n_pairs:]

    effect_i = fit_i - background_fitness
    effect_j = fit_j - background_fitness
    effect_i_on_j = fit_double - fit_j
    effect_j_on_i = fit_double - fit_i
    epistasis = fit_double - fit_i - fit_j + background_fitness

    rows: list[dict[str, Any]] = []
    sign_flags = []
    reciprocal_flags = []
    diminishing_flags = []
    negative_flags = []
    positive_flags = []
    for idx in range(n_pairs):
        sign_i = sign_with_threshold(float(effect_i[idx]), threshold)
        sign_j = sign_with_threshold(float(effect_j[idx]), threshold)
        sign_i_on_j = sign_with_threshold(float(effect_i_on_j[idx]), threshold)
        sign_j_on_i = sign_with_threshold(float(effect_j_on_i[idx]), threshold)
        flip_i = sign_i != 0 and sign_i_on_j != 0 and sign_i != sign_i_on_j
        flip_j = sign_j != 0 and sign_j_on_i != 0 and sign_j != sign_j_on_i
        sign_epistasis = bool(flip_i or flip_j)
        reciprocal_sign = bool(flip_i and flip_j)
        negative_epistasis = bool(epistasis[idx] < -threshold)
        positive_epistasis = bool(epistasis[idx] > threshold)
        beneficial_context = effect_i[idx] > threshold or effect_j[idx] > threshold
        diminishing_returns = bool(beneficial_context and negative_epistasis)
        sign_flags.append(sign_epistasis)
        reciprocal_flags.append(reciprocal_sign)
        diminishing_flags.append(diminishing_returns)
        negative_flags.append(negative_epistasis)
        positive_flags.append(positive_epistasis)
        rows.append({
            **row_context,
            "background_id": background_id,
            "pair_index": int(idx),
            "pos_i": int(pos_i[idx]),
            "aa_i": int(aa_i[idx]),
            "pos_j": int(pos_j[idx]),
            "aa_j": int(aa_j[idx]),
            "background_fitness": background_fitness,
            "fitness_i": float(fit_i[idx]),
            "fitness_j": float(fit_j[idx]),
            "fitness_double": float(fit_double[idx]),
            "effect_i": float(effect_i[idx]),
            "effect_j": float(effect_j[idx]),
            "effect_i_on_j": float(effect_i_on_j[idx]),
            "effect_j_on_i": float(effect_j_on_i[idx]),
            "epistasis": float(epistasis[idx]),
            "abs_epistasis": float(abs(epistasis[idx])),
            "negative_epistasis": int(negative_epistasis),
            "positive_epistasis": int(positive_epistasis),
            "sign_epistasis": int(sign_epistasis),
            "reciprocal_sign_epistasis": int(reciprocal_sign),
            "diminishing_returns_epistasis": int(diminishing_returns),
        })

    summary = {
        "epistasis_pair_count": int(n_pairs),
        "epistasis_mean": float(np.mean(epistasis)) if epistasis.size else math.nan,
        "epistasis_abs_mean": float(np.mean(np.abs(epistasis))) if epistasis.size else math.nan,
        "epistasis_variance": float(np.var(epistasis)) if epistasis.size else math.nan,
        "negative_epistasis_fraction": float(np.mean(negative_flags)) if negative_flags else math.nan,
        "positive_epistasis_fraction": float(np.mean(positive_flags)) if positive_flags else math.nan,
        "sign_epistasis_fraction": float(np.mean(sign_flags)) if sign_flags else math.nan,
        "reciprocal_sign_epistasis_fraction": float(np.mean(reciprocal_flags)) if reciprocal_flags else math.nan,
        "diminishing_returns_epistasis_fraction": float(np.mean(diminishing_flags)) if diminishing_flags else math.nan,
    }
    return rows, summary


def low_high_mean(rows: list[dict[str, Any]], key: str) -> tuple[float | None, float | None]:
    fractions = np.asarray([row["background_fraction_of_peak"] for row in rows], dtype=np.float64)
    values = np.asarray([row[key] for row in rows], dtype=np.float64)
    mask = np.isfinite(fractions) & np.isfinite(values)
    if np.count_nonzero(mask) < 4:
        return None, None
    low_cut = float(np.quantile(fractions[mask], 1 / 3))
    high_cut = float(np.quantile(fractions[mask], 2 / 3))
    low = values[mask & (fractions <= low_cut)]
    high = values[mask & (fractions >= high_cut)]
    return finite_mean(low), finite_mean(high)


def relation_test(
    rows: list[dict[str, Any]],
    *,
    key: str,
    direction: str,
    min_abs_spearman: float,
) -> dict[str, Any]:
    fractions = np.asarray([row["background_fraction_of_peak"] for row in rows], dtype=np.float64)
    values = np.asarray([row[key] for row in rows], dtype=np.float64)
    rho = spearman(fractions, values)
    low, high = low_high_mean(rows, key)
    if direction == "decrease":
        pass_flag = (
            rho is not None
            and rho <= -float(min_abs_spearman)
            and low is not None
            and high is not None
            and high < low
        )
    elif direction == "increase":
        pass_flag = (
            rho is not None
            and rho >= float(min_abs_spearman)
            and low is not None
            and high is not None
            and high > low
        )
    else:
        raise ValueError(f"Unsupported relation direction: {direction}")
    return {
        "metric": key,
        "direction": direction,
        "spearman_fraction_of_peak": rho,
        "low_fraction_mean": low,
        "high_fraction_mean": high,
        "min_abs_spearman": float(min_abs_spearman),
        "pass": bool(pass_flag),
    }


def walk_relation_test(
    summaries: list[dict[str, Any]],
    *,
    key: str,
    direction: str,
    min_abs_spearman: float,
) -> dict[str, Any]:
    starts = np.asarray([row["start_fraction_of_peak"] for row in summaries], dtype=np.float64)
    values = np.asarray([row[key] for row in summaries], dtype=np.float64)
    rho = spearman(starts, values)
    mask = np.isfinite(starts) & np.isfinite(values)
    low = high = None
    if np.count_nonzero(mask) >= 4:
        low_cut = float(np.quantile(starts[mask], 1 / 3))
        high_cut = float(np.quantile(starts[mask], 2 / 3))
        low = finite_mean(values[mask & (starts <= low_cut)])
        high = finite_mean(values[mask & (starts >= high_cut)])
    if direction == "decrease":
        pass_flag = rho is not None and rho <= -min_abs_spearman and low is not None and high is not None and high < low
    elif direction == "increase":
        pass_flag = rho is not None and rho >= min_abs_spearman and low is not None and high is not None and high > low
    else:
        raise ValueError(f"Unsupported relation direction: {direction}")
    return {
        "metric": key,
        "direction": direction,
        "spearman_start_fraction_of_peak": rho,
        "low_start_mean": low,
        "high_start_mean": high,
        "min_abs_spearman": float(min_abs_spearman),
        "pass": bool(pass_flag),
    }


def summarize_validation(
    background_rows: list[dict[str, Any]],
    walk_summaries: list[dict[str, Any]],
    cfg: dict[str, Any],
) -> dict[str, Any]:
    thresholds = cfg["validation_thresholds"]
    min_abs = float(thresholds["min_abs_spearman"])
    decel_values = np.asarray(
        [
            row["deceleration_ratio"]
            for row in walk_summaries
            if row["mode"] == "greedy" and np.isfinite(row["deceleration_ratio"])
        ],
        dtype=np.float64,
    )
    median_deceleration = finite_mean(decel_values)
    if decel_values.size:
        median_deceleration = float(np.median(decel_values))
    max_deceleration = float(thresholds["max_median_deceleration_ratio"])

    tests = {
        "beneficial_fraction_declines_near_peak": relation_test(
            background_rows,
            key="fraction_beneficial",
            direction="decrease",
            min_abs_spearman=min_abs,
        ),
        "mean_beneficial_effect_declines_near_peak": relation_test(
            background_rows,
            key="mean_beneficial_effect",
            direction="decrease",
            min_abs_spearman=min_abs,
        ),
        "adaptive_total_gain_declines_with_start_fitness": walk_relation_test(
            [row for row in walk_summaries if row["mode"] == "greedy"],
            key="total_gain",
            direction="decrease",
            min_abs_spearman=min_abs,
        ),
        "adaptive_walks_decelerate": {
            "metric": "greedy_deceleration_ratio",
            "median": median_deceleration,
            "expected_max": max_deceleration,
            "pass": bool(
                median_deceleration is not None
                and np.isfinite(median_deceleration)
                and median_deceleration < max_deceleration
            ),
        },
        "mean_abs_epistasis_increases_near_peak": relation_test(
            background_rows,
            key="epistasis_abs_mean",
            direction="increase",
            min_abs_spearman=min_abs,
        ),
        "sign_epistasis_increases_near_peak": relation_test(
            background_rows,
            key="sign_epistasis_fraction",
            direction="increase",
            min_abs_spearman=float(thresholds["min_sign_epistasis_spearman"]),
        ),
        "negative_epistasis_increases_near_peak": relation_test(
            background_rows,
            key="negative_epistasis_fraction",
            direction="increase",
            min_abs_spearman=float(thresholds["min_negative_epistasis_spearman"]),
        ),
    }
    pass_count = sum(1 for test in tests.values() if bool(test.get("pass")))
    return {
        "tests": tests,
        "pass_count": int(pass_count),
        "test_count": int(len(tests)),
        "all_pass": bool(pass_count == len(tests)),
    }


def run_experiment(cfg: dict[str, Any]) -> dict[str, Any]:
    sampling_cfg = cfg["sampling"]
    walk_rows: list[dict[str, Any]] = []
    walk_summaries: list[dict[str, Any]] = []
    background_rows: list[dict[str, Any]] = []
    epistasis_rows: list[dict[str, Any]] = []

    for condition_idx, condition in enumerate(cfg["conditions"]):
        condition_name = str(condition["name"])
        for seed in cfg["sampling"]["seeds"]:
            landscape_seed = int(seed) + condition_idx * 100_000
            landscape = build_landscape(cfg, seed=landscape_seed, condition=condition)
            peak_fitness = float(landscape.peak_fitness)
            for level_idx, level in enumerate(cfg["sampling"]["target_fractions"]):
                level_name = str(level["name"])
                target_fraction = float(level["fraction"])
                for replicate in range(int(sampling_cfg["starts_per_level"])):
                    rng_seed = (
                        int(seed)
                        + condition_idx * 1_000_000
                        + level_idx * 10_000
                        + replicate * 101
                    )
                    rng = np.random.RandomState(rng_seed)
                    start_seq, start_fitness = sequence_near_fraction_of_peak(
                        landscape,
                        target_fraction,
                        rng,
                        n_candidates=int(sampling_cfg["candidate_count"]),
                    )
                    background_id = (
                        f"{condition_name}__seed-{seed}__{level_name}__rep-{replicate}"
                    )
                    start_fraction = start_fitness / peak_fitness if abs(peak_fitness) > 1e-12 else math.nan
                    distance_to_peak = int(landscape.hamming_distance(start_seq, landscape.peak_sequence))
                    row_context = {
                        "condition": condition_name,
                        "seed": int(seed),
                        "landscape_seed": int(landscape_seed),
                        "target_level": level_name,
                        "target_fraction_of_peak": target_fraction,
                        "replicate": int(replicate),
                    }

                    scan = single_mutant_scan(
                        landscape,
                        start_seq,
                        sampling_cfg=sampling_cfg,
                    )
                    pair_rows, pair_summary = sample_epistasis_pairs(
                        landscape,
                        start_seq,
                        background_id=background_id,
                        row_context=row_context,
                        rng=rng,
                        sampling_cfg=sampling_cfg,
                    )
                    epistasis_rows.extend(pair_rows)
                    background_rows.append({
                        **row_context,
                        "background_id": background_id,
                        "background_fitness": float(scan["parent_fitness"]),
                        "peak_fitness": peak_fitness,
                        "background_fraction_of_peak": start_fraction,
                        "distance_to_peak": distance_to_peak,
                        "fraction_beneficial": float(scan["fraction_beneficial"]),
                        "fraction_neutral": float(scan["fraction_neutral"]),
                        "fraction_deleterious": float(scan["fraction_deleterious"]),
                        "mean_delta": float(scan["mean_delta"]),
                        "best_delta": float(scan["best_delta"]),
                        "mean_beneficial_effect": float(scan["mean_beneficial_effect"]),
                        **pair_summary,
                    })

                    for mode in cfg["sampling"]["adaptive_walk_modes"]:
                        walk_rng = np.random.RandomState(rng_seed + (13 if mode == "greedy" else 37))
                        walk = run_adaptive_walk(
                            landscape,
                            start_seq,
                            mode=str(mode),
                            rng=walk_rng,
                            sampling_cfg=sampling_cfg,
                        )
                        walk_summary = {
                            **row_context,
                            "background_id": background_id,
                            "mode": str(mode),
                            "start_fraction_of_peak": start_fraction,
                            "start_distance_to_peak": distance_to_peak,
                            **walk["summary"],
                        }
                        walk_summaries.append(walk_summary)
                        start_fit = float(walk["trajectory"][0]["fitness"])
                        for row in walk["trajectory"]:
                            walk_rows.append({
                                **row_context,
                                "background_id": background_id,
                                "mode": str(mode),
                                "start_fraction_of_peak": start_fraction,
                                "start_distance_to_peak": distance_to_peak,
                                "step": int(row["step"]),
                                "fitness": float(row["fitness"]),
                                "fraction_of_peak": float(row["fraction_of_peak"]),
                                "distance_to_peak": int(row["distance_to_peak"]),
                                "cumulative_gain": float(row["fitness"] - start_fit),
                                "step_gain": float(row["step_gain"]),
                                "available_beneficial_fraction": float(row["available_beneficial_fraction"]),
                                "available_mean_beneficial_effect": float(row["available_mean_beneficial_effect"]),
                                "best_available_effect": float(row["best_available_effect"]),
                                "chosen_effect": float(row["chosen_effect"]),
                            })

    validation = summarize_validation(background_rows, walk_summaries, cfg)
    return {
        "walk_rows": walk_rows,
        "walk_summaries": walk_summaries,
        "background_rows": background_rows,
        "epistasis_rows": epistasis_rows,
        "summary": {
            "experiment_id": cfg["experiment_id"],
            "hypothesis_id": "HYP-004",
            "title": cfg["title"],
            "n_backgrounds": len(background_rows),
            "n_walk_summaries": len(walk_summaries),
            "n_walk_trajectory_rows": len(walk_rows),
            "n_epistasis_pairs": len(epistasis_rows),
            "conditions": [condition["name"] for condition in cfg["conditions"]],
            "literature_baselines": LITERATURE_BASELINES,
            "validation": validation,
        },
    }


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config_path, quick=bool(args.quick))

    payload = run_experiment(cfg)
    write_json(output_dir / "literature_baselines.json", LITERATURE_BASELINES)
    write_json(output_dir / "summary.json", payload["summary"])
    write_json(output_dir / "walk_summaries.json", payload["walk_summaries"])

    write_csv(
        output_dir / "walk_trajectories.csv",
        payload["walk_rows"],
        [
            "condition",
            "seed",
            "landscape_seed",
            "target_level",
            "target_fraction_of_peak",
            "replicate",
            "background_id",
            "mode",
            "start_fraction_of_peak",
            "start_distance_to_peak",
            "step",
            "fitness",
            "fraction_of_peak",
            "distance_to_peak",
            "cumulative_gain",
            "step_gain",
            "available_beneficial_fraction",
            "available_mean_beneficial_effect",
            "best_available_effect",
            "chosen_effect",
        ],
    )
    write_csv(
        output_dir / "background_neighborhoods.csv",
        payload["background_rows"],
        [
            "condition",
            "seed",
            "landscape_seed",
            "target_level",
            "target_fraction_of_peak",
            "replicate",
            "background_id",
            "background_fitness",
            "peak_fitness",
            "background_fraction_of_peak",
            "distance_to_peak",
            "fraction_beneficial",
            "fraction_neutral",
            "fraction_deleterious",
            "mean_delta",
            "best_delta",
            "mean_beneficial_effect",
            "epistasis_pair_count",
            "epistasis_mean",
            "epistasis_abs_mean",
            "epistasis_variance",
            "negative_epistasis_fraction",
            "positive_epistasis_fraction",
            "sign_epistasis_fraction",
            "reciprocal_sign_epistasis_fraction",
            "diminishing_returns_epistasis_fraction",
        ],
    )
    write_csv(
        output_dir / "epistasis_pairs.csv",
        payload["epistasis_rows"],
        [
            "condition",
            "seed",
            "landscape_seed",
            "target_level",
            "target_fraction_of_peak",
            "replicate",
            "background_id",
            "pair_index",
            "pos_i",
            "aa_i",
            "pos_j",
            "aa_j",
            "background_fitness",
            "fitness_i",
            "fitness_j",
            "fitness_double",
            "effect_i",
            "effect_j",
            "effect_i_on_j",
            "effect_j_on_i",
            "epistasis",
            "abs_epistasis",
            "negative_epistasis",
            "positive_epistasis",
            "sign_epistasis",
            "reciprocal_sign_epistasis",
            "diminishing_returns_epistasis",
        ],
    )
    print(json.dumps(payload["summary"]["validation"], indent=2, sort_keys=True, default=json_default))


if __name__ == "__main__":
    main()
