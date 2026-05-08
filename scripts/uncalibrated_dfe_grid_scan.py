#!/usr/bin/env python3
"""Run EXP-002: uncalibrated single- and double-mutant DFE grid scan."""

from __future__ import annotations

import argparse
import csv
import itertools
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

from adaptenv import FitnessLandscapeEnv, LandscapeConfig  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


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
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def setting_slug(parts: dict[str, Any]) -> str:
    formatted = []
    for key, value in parts.items():
        if isinstance(value, float):
            token = f"{value:.3f}".rstrip("0").rstrip(".")
        else:
            token = str(value)
        formatted.append(f"{key}-{token}")
    return "__".join(formatted)


def enumerate_single_mutants(seq: np.ndarray, alphabet_size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    seq = np.asarray(seq, dtype=np.int64)
    L = seq.shape[0]
    mutant_count = L * (alphabet_size - 1)
    mutants = np.repeat(seq[None, :], mutant_count, axis=0)
    positions = np.repeat(np.arange(L, dtype=np.int64), alphabet_size - 1)
    amino_acids = np.empty(mutant_count, dtype=np.int64)

    offset = 0
    for pos in range(L):
        alts = np.delete(np.arange(alphabet_size, dtype=np.int64), seq[pos])
        n_alts = alts.size
        mutants[offset:offset + n_alts, pos] = alts
        amino_acids[offset:offset + n_alts] = alts
        offset += n_alts
    return mutants, positions, amino_acids


def build_pair_index(L: int) -> np.ndarray:
    return np.asarray(list(itertools.combinations(range(L), 2)), dtype=np.int64)


def sample_double_mutants(
    seq: np.ndarray,
    *,
    pair_index: np.ndarray,
    n_samples: int,
    alphabet_size: int,
    rng: np.random.RandomState,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    seq = np.asarray(seq, dtype=np.int64)
    chosen_pairs = pair_index[rng.randint(0, pair_index.shape[0], size=n_samples)]
    left_positions = chosen_pairs[:, 0]
    right_positions = chosen_pairs[:, 1]

    left_aas = rng.randint(0, alphabet_size - 1, size=n_samples)
    right_aas = rng.randint(0, alphabet_size - 1, size=n_samples)
    left_current = seq[left_positions]
    right_current = seq[right_positions]
    left_aas = left_aas + (left_aas >= left_current)
    right_aas = right_aas + (right_aas >= right_current)

    mutants = np.repeat(seq[None, :], n_samples, axis=0)
    rows = np.arange(n_samples, dtype=np.int64)
    mutants[rows, left_positions] = left_aas
    mutants[rows, right_positions] = right_aas

    positions = np.column_stack([left_positions, right_positions]).astype(np.int64)
    amino_acids = np.column_stack([left_aas, right_aas]).astype(np.int64)
    return mutants, positions, amino_acids


def evaluate_components_in_batches(
    landscape: FitnessLandscapeEnv,
    seqs: np.ndarray,
    *,
    batch_size: int,
) -> dict[str, np.ndarray]:
    fitness_parts: list[np.ndarray] = []
    stability_parts: list[np.ndarray] = []
    for start in range(0, seqs.shape[0], batch_size):
        stop = min(start + batch_size, seqs.shape[0])
        batch_components = landscape.evaluate_batch_components(seqs[start:stop])
        fitness_parts.append(np.asarray(batch_components["fitness"], dtype=np.float64))
        stability_parts.append(np.asarray(batch_components["stability"], dtype=np.float64))
    return {
        "fitness": np.concatenate(fitness_parts, axis=0),
        "stability": np.concatenate(stability_parts, axis=0),
    }


def threshold_from_parent(
    parent_fitness: float,
    *,
    minimum: float,
    fraction: float,
) -> float:
    return float(max(minimum, fraction * abs(parent_fitness)))


def quantiles(values: np.ndarray) -> dict[str, float]:
    return {
        "q05": float(np.quantile(values, 0.05)),
        "q25": float(np.quantile(values, 0.25)),
        "q50": float(np.quantile(values, 0.50)),
        "q75": float(np.quantile(values, 0.75)),
        "q95": float(np.quantile(values, 0.95)),
    }


def dfe_stats(
    deltas: np.ndarray,
    stabilities: np.ndarray,
    *,
    neutral_threshold: float,
    stability_threshold: float,
) -> dict[str, float]:
    deltas = np.asarray(deltas, dtype=np.float64)
    stabilities = np.asarray(stabilities, dtype=np.float64)
    beneficial = deltas > neutral_threshold
    deleterious = deltas < -neutral_threshold
    neutral = np.abs(deltas) <= neutral_threshold
    lethal = stabilities < stability_threshold

    centered = deltas - np.mean(deltas)
    variance = float(np.var(deltas))
    std = math.sqrt(max(variance, 0.0))
    if deltas.size > 2 and std > 1e-12:
        skewness = float(np.mean(centered ** 3) / (std ** 3))
    else:
        skewness = 0.0

    payload = {
        "count": int(deltas.size),
        "neutral_threshold": float(neutral_threshold),
        "fraction_beneficial": float(np.mean(beneficial)),
        "fraction_neutral": float(np.mean(neutral)),
        "fraction_deleterious": float(np.mean(deleterious)),
        "fraction_lethal": float(np.mean(lethal)),
        "delta_mean": float(np.mean(deltas)),
        "delta_variance": variance,
        "delta_std": std,
        "delta_skewness": skewness,
        "beneficial_mean": float(np.mean(deltas[beneficial])) if np.any(beneficial) else 0.0,
        "deleterious_mean": float(np.mean(deltas[deleterious])) if np.any(deleterious) else 0.0,
        "delta_min": float(np.min(deltas)),
        "delta_max": float(np.max(deltas)),
    }
    payload.update(quantiles(deltas))
    return payload


def epistasis_stats(epistasis: np.ndarray, *, threshold: float) -> dict[str, float]:
    epistasis = np.asarray(epistasis, dtype=np.float64)
    variance = float(np.var(epistasis))
    return {
        "epistasis_mean": float(np.mean(epistasis)),
        "epistasis_std": math.sqrt(max(variance, 0.0)),
        "epistasis_variance": variance,
        "epistasis_abs_mean": float(np.mean(np.abs(epistasis))),
        "fraction_positive_epistasis": float(np.mean(epistasis > threshold)),
        "fraction_negative_epistasis": float(np.mean(epistasis < -threshold)),
        "fraction_near_additive": float(np.mean(np.abs(epistasis) <= threshold)),
        "epistasis_q05": float(np.quantile(epistasis, 0.05)),
        "epistasis_q50": float(np.quantile(epistasis, 0.50)),
        "epistasis_q95": float(np.quantile(epistasis, 0.95)),
    }


def build_single_lookup(
    positions: np.ndarray,
    amino_acids: np.ndarray,
    deltas: np.ndarray,
    *,
    L: int,
    alphabet_size: int,
) -> np.ndarray:
    lookup = np.full((L, alphabet_size), np.nan, dtype=np.float64)
    lookup[positions, amino_acids] = deltas
    return lookup


def flatten_row(
    setting_id: str,
    parameter_values: dict[str, Any],
    seed: int,
    parent_payload: dict[str, float],
    single_stats_payload: dict[str, float],
    double_stats_payload: dict[str, float],
    epistasis_payload: dict[str, float],
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "setting_id": setting_id,
        "seed": int(seed),
        **parameter_values,
        **parent_payload,
    }
    for key, value in single_stats_payload.items():
        row[f"single_{key}"] = value
    for key, value in double_stats_payload.items():
        row[f"double_{key}"] = value
    for key, value in epistasis_payload.items():
        row[key] = value
    return row


def aggregate_numeric_rows(rows: list[dict[str, Any]], keys: list[str]) -> dict[str, float]:
    payload: dict[str, float] = {}
    for key in keys:
        values = [float(row[key]) for row in rows]
        payload[f"{key}_mean"] = float(np.mean(values))
        payload[f"{key}_std"] = float(np.std(values))
    return payload


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows available for CSV output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_experiment(config: dict[str, Any], *, quick: bool) -> dict[str, Any]:
    effective = merge_dict(config, config.get("quick", {})) if quick else dict(config)
    base_config = dict(effective["adapt_env"]["base_config"])
    landscape_cfg = effective["landscape"]
    if bool(landscape_cfg.get("conditioned", False)):
        raise ValueError("This experiment currently supports only unconditioned landscapes.")

    L = int(landscape_cfg["sequence_length"])
    seeds = [int(seed) for seed in effective["sampling"]["seeds"]]
    n_double_samples = int(effective["sampling"]["double_mutants_per_landscape"])
    batch_size = int(effective["sampling"]["batch_size"])
    neutral_minimum = float(effective["sampling"]["neutral_threshold_min"])
    neutral_fraction = float(effective["sampling"]["neutral_threshold_fraction_of_parent"])
    grid = effective["grid"]
    pair_index = build_pair_index(L)

    grid_names = [
        "stability_margin",
        "functional_sigma_base",
        "n_functional_dims",
        "epistasis_strength",
    ]
    grid_values = [grid[name] for name in grid_names]
    per_landscape_rows: list[dict[str, Any]] = []
    per_setting_payload: list[dict[str, Any]] = []

    for combination in itertools.product(*grid_values):
        parameter_values = {
            name: float(value) if isinstance(value, float) else int(value)
            for name, value in zip(grid_names, combination)
        }
        setting_id = setting_slug(parameter_values)
        replicate_rows: list[dict[str, Any]] = []
        for seed in seeds:
            config_kwargs = dict(base_config)
            config_kwargs.update(parameter_values)
            config_kwargs["L"] = L
            config_kwargs["seed"] = int(seed)

            landscape = FitnessLandscapeEnv(LandscapeConfig(**config_kwargs))
            parent = landscape.evaluate(landscape.reference.copy(), return_components=True)
            parent_fitness = float(parent["fitness"])
            neutral_threshold = threshold_from_parent(
                parent_fitness,
                minimum=neutral_minimum,
                fraction=neutral_fraction,
            )

            single_mutants, single_positions, single_aas = enumerate_single_mutants(
                landscape.reference.copy(),
                landscape.A,
            )
            single_components = evaluate_components_in_batches(
                landscape,
                single_mutants,
                batch_size=batch_size,
            )
            single_deltas = single_components["fitness"] - parent_fitness
            single_lookup = build_single_lookup(
                single_positions,
                single_aas,
                single_deltas,
                L=landscape.L,
                alphabet_size=landscape.A,
            )

            double_mutants, double_positions, double_aas = sample_double_mutants(
                landscape.reference.copy(),
                pair_index=pair_index,
                n_samples=n_double_samples,
                alphabet_size=landscape.A,
                rng=np.random.RandomState(seed + 10_000),
            )
            double_components = evaluate_components_in_batches(
                landscape,
                double_mutants,
                batch_size=batch_size,
            )
            double_deltas = double_components["fitness"] - parent_fitness
            additive_expectation = (
                single_lookup[double_positions[:, 0], double_aas[:, 0]]
                + single_lookup[double_positions[:, 1], double_aas[:, 1]]
            )
            epistasis = double_deltas - additive_expectation

            parent_payload = {
                "sequence_length": int(landscape.L),
                "parent_fitness": parent_fitness,
                "parent_stability": float(parent["stability"]),
                "parent_margin": float(parent["stability_margin_used"]),
                "stability_threshold": float(landscape.stability_threshold),
                "peak_distance_from_reference": int(
                    landscape.hamming_distance(landscape.reference, landscape.peak_sequence)
                ),
            }
            single_stats_payload = dfe_stats(
                single_deltas,
                single_components["stability"],
                neutral_threshold=neutral_threshold,
                stability_threshold=float(landscape.stability_threshold),
            )
            double_stats_payload = dfe_stats(
                double_deltas,
                double_components["stability"],
                neutral_threshold=neutral_threshold,
                stability_threshold=float(landscape.stability_threshold),
            )
            epistasis_payload = epistasis_stats(epistasis, threshold=neutral_threshold)

            row = flatten_row(
                setting_id,
                parameter_values,
                seed,
                parent_payload,
                single_stats_payload,
                double_stats_payload,
                epistasis_payload,
            )
            per_landscape_rows.append(row)
            replicate_rows.append(row)

        metric_keys = [
            "parent_fitness",
            "parent_stability",
            "parent_margin",
            "peak_distance_from_reference",
            "single_fraction_beneficial",
            "single_fraction_neutral",
            "single_fraction_deleterious",
            "single_fraction_lethal",
            "single_delta_mean",
            "single_delta_variance",
            "single_delta_skewness",
            "single_q05",
            "single_q50",
            "single_q95",
            "double_fraction_beneficial",
            "double_fraction_neutral",
            "double_fraction_deleterious",
            "double_fraction_lethal",
            "double_delta_mean",
            "double_delta_variance",
            "double_delta_skewness",
            "double_q05",
            "double_q50",
            "double_q95",
            "epistasis_mean",
            "epistasis_std",
            "epistasis_abs_mean",
            "fraction_positive_epistasis",
            "fraction_negative_epistasis",
            "fraction_near_additive",
        ]
        per_setting_payload.append(
            {
                "setting_id": setting_id,
                **parameter_values,
                "n_replicates": len(replicate_rows),
                **aggregate_numeric_rows(replicate_rows, metric_keys),
            }
        )

    axis_summaries: dict[str, list[dict[str, Any]]] = {}
    for axis_name in grid_names:
        axis_entries: list[dict[str, Any]] = []
        unique_values = sorted({entry[axis_name] for entry in per_landscape_rows})
        for value in unique_values:
            subset = [row for row in per_landscape_rows if row[axis_name] == value]
            axis_entries.append(
                {
                    axis_name: value,
                    "n_landscapes": len(subset),
                    **aggregate_numeric_rows(
                        subset,
                        [
                            "single_fraction_beneficial",
                            "single_fraction_neutral",
                            "single_fraction_lethal",
                            "double_fraction_beneficial",
                            "double_fraction_neutral",
                            "double_fraction_lethal",
                            "epistasis_abs_mean",
                            "fraction_positive_epistasis",
                            "fraction_negative_epistasis",
                        ],
                    ),
                }
            )
        axis_summaries[axis_name] = axis_entries

    summary = {
        "n_settings": len(per_setting_payload),
        "n_landscapes": len(per_landscape_rows),
        "double_more_deleterious_than_single_fraction": float(
            np.mean(
                [
                    row["double_delta_mean"] < row["single_delta_mean"]
                    for row in per_landscape_rows
                ]
            )
        ),
        "double_more_lethal_than_single_fraction": float(
            np.mean(
                [
                    row["double_fraction_lethal"] > row["single_fraction_lethal"]
                    for row in per_landscape_rows
                ]
            )
        ),
        "highest_epistasis_setting": max(
            per_setting_payload,
            key=lambda row: row["epistasis_abs_mean_mean"],
        )["setting_id"],
        "lowest_epistasis_setting": min(
            per_setting_payload,
            key=lambda row: row["epistasis_abs_mean_mean"],
        )["setting_id"],
    }

    return {
        "metadata": {
            "experiment_id": effective["experiment"]["id"],
            "quick": bool(quick),
            "sequence_length": L,
            "seeds": seeds,
            "double_mutants_per_landscape": n_double_samples,
            "batch_size": batch_size,
            "neutral_threshold_min": neutral_minimum,
            "neutral_threshold_fraction_of_parent": neutral_fraction,
            "base_config": base_config,
            "grid": grid,
        },
        "results": {
            "per_landscape_rows": per_landscape_rows,
            "per_setting_rows": per_setting_payload,
            "axis_summaries": axis_summaries,
        },
        "summary": summary,
    }


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config))

    payload = run_experiment(config, quick=bool(args.quick))
    per_landscape_rows = payload["results"]["per_landscape_rows"]
    per_setting_rows = payload["results"]["per_setting_rows"]

    write_csv(output_dir / "single_double_dfe_per_landscape.csv", per_landscape_rows)
    write_csv(output_dir / "single_double_dfe_per_setting.csv", per_setting_rows)
    (output_dir / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    print(output_dir / "summary.json")


if __name__ == "__main__":
    main()
