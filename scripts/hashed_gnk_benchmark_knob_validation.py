#!/usr/bin/env python3
"""Validate hashed GNK benchmark-generator knob maps on enumerable coordinates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from statistics import NormalDist
from typing import Any, Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from phot_coupled_landscape_grid_search import load_config, write_csv, write_json  # noqa: E402
from phot_nested_fit_protocol_benchmark import epistasis_property_metrics  # noqa: E402
from project_adapt_env.dense_descriptive import exact_summary  # noqa: E402


NORMAL = NormalDist()


@dataclass(frozen=True)
class GNKTerm:
    term_id: int
    order: int
    support: tuple[int, ...]
    coefficient: float
    support_mode: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def stable_uint64(*parts: Any, seed: int = 0) -> int:
    payload = json.dumps([seed, *parts], sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little")


def stable_uniform(*parts: Any, seed: int = 0) -> float:
    return (stable_uint64(*parts, seed=seed) + 0.5) / float(2**64)


def stable_gaussian(*parts: Any, seed: int = 0) -> float:
    u = min(max(stable_uniform(*parts, seed=seed), 1e-15), 1.0 - 1e-15)
    return float(NORMAL.inv_cdf(u))


def genotype_matrix(length: int) -> np.ndarray:
    nodes = np.arange(1 << int(length), dtype=np.uint32)
    bits = ((nodes[:, None] >> np.arange(int(length), dtype=np.uint32)[None, :]) & 1).astype(np.float64)
    return bits


def sigmoid(x: np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=np.float64)
    return np.where(arr >= 0, 1.0 / (1.0 + np.exp(-arr)), np.exp(arr) / (1.0 + np.exp(arr)))


def deterministic_permutation(length: int, *, seed: int, label: str) -> list[int]:
    return sorted(range(int(length)), key=lambda pos: stable_uint64(label, pos, seed=seed))


def structured_support_pool(length: int, order: int, *, seed: int, mode: str) -> list[tuple[int, ...]]:
    length = int(length)
    order = int(order)
    if order > length:
        return []
    if mode == "random":
        pool = [tuple(range(length))]
    elif mode == "shuffled_block":
        perm = deterministic_permutation(length, seed=seed, label="shuffled-block")
        pool = []
        window = max(order + 1, min(length, order + 3))
        for start in range(length):
            block = [perm[(start + offset) % length] for offset in range(window)]
            pool.append(tuple(sorted(block)))
    else:
        pool = []
        window = max(order + 1, min(length, order + 3))
        for start in range(length):
            block = [(start + offset) % length for offset in range(window)]
            pool.append(tuple(sorted(block)))
    supports: set[tuple[int, ...]] = set()
    for block in pool:
        candidates = list(block)
        if mode == "random":
            # Add deterministic pseudo-random supports from the whole sequence.
            for term in range(max(length * 4, 16)):
                ranked = sorted(candidates, key=lambda pos: stable_uint64("random-support", order, term, pos, seed=seed))
                supports.add(tuple(sorted(ranked[:order])))
        else:
            for offset in range(len(candidates)):
                ranked = sorted(candidates, key=lambda pos: stable_uint64("support", mode, order, offset, pos, seed=seed))
                supports.add(tuple(sorted(ranked[:order])))
    return sorted(supports)


def allocate_counts(order_weights: dict[int, float], total_terms: int) -> dict[int, int]:
    positive = {int(order): max(float(weight), 0.0) for order, weight in order_weights.items()}
    total = sum(positive.values())
    if total <= 0.0:
        return {}
    expected = {order: float(total_terms) * weight / total for order, weight in positive.items()}
    counts = {order: int(math.floor(value)) for order, value in expected.items()}
    remaining = int(total_terms) - sum(counts.values())
    ranked = sorted(expected, key=lambda order: (-(expected[order] - counts[order]), order))
    for order in ranked[: max(remaining, 0)]:
        counts[order] += 1
    return {order: count for order, count in counts.items() if count > 0}


def build_terms(
    *,
    length: int,
    seed: int,
    support_mode: str,
    order_weights: dict[int, float],
    total_terms: int,
) -> list[GNKTerm]:
    counts = allocate_counts(order_weights, total_terms)
    terms: list[GNKTerm] = []
    term_id = 0
    for order, count in sorted(counts.items()):
        pool = structured_support_pool(length, order, seed=seed, mode=support_mode)
        if not pool:
            continue
        coeff_sd = math.sqrt(max(float(order_weights[order]), 0.0) / max(int(count), 1))
        for local_index in range(int(count)):
            support = pool[stable_uint64("term-support", support_mode, order, local_index, seed=seed) % len(pool)]
            coeff = coeff_sd * stable_gaussian("term-coeff", support_mode, order, local_index, support, seed=seed)
            terms.append(
                GNKTerm(
                    term_id=term_id,
                    order=int(order),
                    support=tuple(int(pos) for pos in support),
                    coefficient=float(coeff),
                    support_mode=str(support_mode),
                )
            )
            term_id += 1
    return terms


def term_rows(parameter_set_id: str, terms: Iterable[GNKTerm]) -> list[dict[str, Any]]:
    return [
        {
            "parameter_set_id": parameter_set_id,
            "term_id": int(term.term_id),
            "support_mode": term.support_mode,
            "order": int(term.order),
            "support_positions": ";".join(str(pos) for pos in term.support),
            "coefficient": float(term.coefficient),
        }
        for term in terms
    ]


def epistatic_field(bits: np.ndarray, terms: list[GNKTerm]) -> np.ndarray:
    if not terms:
        return np.zeros(bits.shape[0], dtype=np.float64)
    total = np.zeros(bits.shape[0], dtype=np.float64)
    contrast = 2.0 * np.asarray(bits, dtype=np.float64) - 1.0
    for term in terms:
        product = np.prod(contrast[:, list(term.support)], axis=1)
        total += float(term.coefficient) * product
    return total


def hamming_layer_rows(
    *,
    parameter: dict[str, Any],
    layer: str,
    values: np.ndarray,
    length: int,
) -> list[dict[str, Any]]:
    arr = np.asarray(values, dtype=np.float64)
    rows: list[dict[str, Any]] = []
    by_layer: dict[int, list[dict[str, float]]] = {hamming: [] for hamming in range(int(length) + 1)}
    for node in range(1 << int(length)):
        hamming = int(node.bit_count())
        effects = np.asarray([float(arr[node ^ (1 << bit)] - arr[node]) for bit in range(int(length))], dtype=np.float64)
        beneficial = effects[effects > 1e-12]
        by_layer[hamming].append(
            {
                "fitness": float(arr[node]),
                "beneficial_fraction": float(len(beneficial) / int(length)),
                "mean_beneficial_effect": float(np.mean(beneficial)) if beneficial.size else 0.0,
                "mean_single_step_effect": float(np.mean(effects)),
            }
        )
    for hamming, records in by_layer.items():
        if not records:
            continue
        rows.append(
            {
                "phase": "WS-B",
                "parameter_set_id": str(parameter["parameter_set_id"]),
                "score_layer": layer,
                "support_mode": str(parameter["support_mode"]),
                "seed": int(parameter["seed"]),
                "epi_weight": float(parameter["epi_weight"]),
                "hamming_distance_from_peak": int(hamming),
                "fraction_of_peak": float(1.0 - (hamming / max(int(length), 1))),
                "n_nodes": int(len(records)),
                "mean_fitness": float(np.mean([record["fitness"] for record in records])),
                "beneficial_fraction": float(np.mean([record["beneficial_fraction"] for record in records])),
                "mean_beneficial_effect": float(np.mean([record["mean_beneficial_effect"] for record in records])),
            }
        )
    return rows


def threshold_band(value: float, *, transition_half_width: float) -> str:
    if value < -float(transition_half_width):
        return "deep_unstable"
    if value > float(transition_half_width):
        return "deep_stable"
    return "fold_transition"


def stratified_pairwise_rows(
    *,
    parameter: dict[str, Any],
    layer: str,
    values: np.ndarray,
    latent_values: np.ndarray,
    length: int,
    transition_half_width: float,
) -> list[dict[str, Any]]:
    arr = np.asarray(values, dtype=np.float64)
    latent = np.asarray(latent_values, dtype=np.float64)
    groups: dict[tuple[int, str], list[dict[str, float]]] = {}
    for left in range(int(length)):
        for right in range(left + 1, int(length)):
            others = [idx for idx in range(int(length)) if idx not in {left, right}]
            for bg_index in range(1 << len(others)):
                base = 0
                for offset, bit in enumerate(others):
                    if bg_index & (1 << offset):
                        base |= 1 << bit
                n00 = base
                n10 = base | (1 << left)
                n01 = base | (1 << right)
                n11 = base | (1 << left) | (1 << right)
                left_low = float(arr[n10] - arr[n00])
                left_high = float(arr[n11] - arr[n01])
                right_low = float(arr[n01] - arr[n00])
                right_high = float(arr[n11] - arr[n10])
                left_flip = left_low * left_high < -1e-24
                right_flip = right_low * right_high < -1e-24
                epsilon = float(arr[n11] - arr[n10] - arr[n01] + arr[n00])
                hamming = int(base.bit_count())
                band = threshold_band(
                    float(np.mean([latent[n00], latent[n10], latent[n01], latent[n11]])),
                    transition_half_width=transition_half_width,
                )
                groups.setdefault((hamming, band), []).append(
                    {
                        "abs_epsilon": abs(epsilon),
                        "epsilon": epsilon,
                        "sign": float(left_flip or right_flip),
                        "reciprocal": float(left_flip and right_flip),
                    }
                )
    rows: list[dict[str, Any]] = []
    for (hamming, band), records in sorted(groups.items()):
        rows.append(
            {
                "phase": "WS-B",
                "parameter_set_id": str(parameter["parameter_set_id"]),
                "score_layer": layer,
                "support_mode": str(parameter["support_mode"]),
                "seed": int(parameter["seed"]),
                "epi_weight": float(parameter["epi_weight"]),
                "background_hamming_distance_from_peak": int(hamming),
                "fraction_of_peak": float(1.0 - (hamming / max(int(length), 1))),
                "threshold_band": band,
                "n_plaquettes": int(len(records)),
                "mean_abs_pairwise_epistasis": float(np.mean([record["abs_epsilon"] for record in records])),
                "mean_signed_pairwise_epistasis": float(np.mean([record["epsilon"] for record in records])),
                "sign_epistasis_fraction": float(np.mean([record["sign"] for record in records])),
                "reciprocal_sign_epistasis_fraction": float(np.mean([record["reciprocal"] for record in records])),
            }
        )
    return rows


def evaluate_parameter_set(
    config: dict[str, Any],
    parameter: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    gen_cfg = config["generator"]
    length = int(gen_cfg["sequence_length"])
    bits = genotype_matrix(length)
    hamming = np.sum(bits, axis=1)
    readout_scale = max(float(gen_cfg.get("readout_scale", 1.0)), 1e-12)
    additive_margin = float(parameter["stability_offset"]) - float(parameter["core_slope"]) * hamming
    order_weights = {int(k): float(v) for k, v in dict(parameter.get("order_weights", {})).items()}
    terms = build_terms(
        length=length,
        seed=int(parameter.get("seed", gen_cfg.get("seed", 20260618))),
        support_mode=str(parameter["support_mode"]),
        order_weights=order_weights,
        total_terms=int(parameter.get("n_terms", gen_cfg.get("n_terms", 96))),
    )
    epi = epistatic_field(bits, terms)
    latent = additive_margin + float(parameter["epi_weight"]) * epi
    observable = sigmoid(latent / readout_scale)
    low_threshold = float(gen_cfg.get("low_plateau_threshold", 0.05))
    high_threshold = float(gen_cfg.get("high_plateau_threshold", 0.95))
    rows: list[dict[str, Any]] = []
    hamming_rows: list[dict[str, Any]] = []
    epistasis_rows: list[dict[str, Any]] = []
    transition_half_width = float(gen_cfg.get("transition_half_width", readout_scale))
    for layer, values in [("latent_margin", latent), ("observable_fitness", observable)]:
        metric = exact_summary(values, dim=length)
        metric.update(epistasis_property_metrics(values, dim=length))
        finite_metric = {
            key: float(value)
            for key, value in metric.items()
            if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(float(value))
        }
        rows.append(
            {
                "phase": "WS-B",
                "parameter_set_id": str(parameter["parameter_set_id"]),
                "score_layer": layer,
                "support_mode": str(parameter["support_mode"]),
                "seed": int(parameter.get("seed", gen_cfg.get("seed", 20260618))),
                "sequence_length": length,
                "n_terms": int(parameter.get("n_terms", gen_cfg.get("n_terms", 96))),
                "epi_weight": float(parameter["epi_weight"]),
                "core_slope": float(parameter["core_slope"]),
                "stability_offset": float(parameter["stability_offset"]),
                "order_weight_2": float(order_weights.get(2, 0.0)),
                "order_weight_3": float(order_weights.get(3, 0.0)),
                "order_weight_4": float(order_weights.get(4, 0.0)),
                "low_plateau_fraction": float(np.mean(observable <= low_threshold)) if layer == "observable_fitness" else math.nan,
                "high_plateau_fraction": float(np.mean(observable >= high_threshold)) if layer == "observable_fitness" else math.nan,
                **finite_metric,
            }
        )
        hamming_rows.extend(
            hamming_layer_rows(
                parameter=parameter,
                layer=layer,
                values=values,
                length=length,
            )
        )
        epistasis_rows.extend(
            stratified_pairwise_rows(
                parameter=parameter,
                layer=layer,
                values=values,
                latent_values=latent,
                length=length,
                transition_half_width=transition_half_width,
            )
        )
    return rows, term_rows(str(parameter["parameter_set_id"]), terms), hamming_rows, epistasis_rows


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    sorted_values = values[order]
    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and sorted_values[end] == sorted_values[start]:
            end += 1
        rank = 0.5 * (start + end - 1)
        ranks[order[start:end]] = rank
        start = end
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    xa = np.asarray(x, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(xa) & np.isfinite(ya)
    xa = xa[mask]
    ya = ya[mask]
    if xa.size < 2 or float(np.std(xa)) <= 1e-12 or float(np.std(ya)) <= 1e-12:
        return math.nan
    return float(np.corrcoef(rankdata(xa), rankdata(ya))[0, 1])


def summarize_knobs(config: dict[str, Any], metric_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    df = pd.DataFrame(metric_rows)
    obs = df[df["score_layer"] == "observable_fitness"].copy()
    key_metrics = [
        "n_local_maxima",
        "accessible_path_fraction",
        "walsh_high_order_fraction",
        "sign_epistasis_fraction",
        "reciprocal_sign_epistasis_fraction",
        "magnitude_epistasis_to_single_mutant_sd_ratio",
        "random_walk_lag1_autocorrelation",
        "high_plateau_fraction",
    ]
    summary_rows: list[dict[str, Any]] = []
    for row in obs.to_dict(orient="records"):
        summary_rows.append(
            {
                "phase": "WS-B",
                "parameter_set_id": row["parameter_set_id"],
                "support_mode": row["support_mode"],
                "epi_weight": row["epi_weight"],
                **{metric: row.get(metric, math.nan) for metric in key_metrics},
            }
        )
    checks: list[dict[str, Any]] = []
    expectations = dict(config["expectations"])
    for support_mode, sub in obs.groupby("support_mode", dropna=False):
        if len(sub) < 2:
            continue
        weights = [float(value) for value in sub["epi_weight"].tolist()]
        for metric, direction in expectations.items():
            if metric not in sub:
                continue
            values = [float(value) for value in sub[metric].tolist()]
            rho = spearman(weights, values)
            expected_sign = 1.0 if direction == "increase" else -1.0
            checks.append(
                {
                    "phase": "WS-B",
                    "support_mode": support_mode,
                    "knob": "epi_weight",
                    "metric": metric,
                    "expected_direction": direction,
                    "spearman_rho": rho,
                    "direction_pass": bool(np.isfinite(rho) and rho * expected_sign > 0.0),
                    "n_parameter_sets": int(len(sub)),
                }
            )
    return summary_rows, checks


def summarize_core_guardrails(
    metric_rows: list[dict[str, Any]],
    hamming_rows: list[dict[str, Any]],
    epistasis_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    metrics = pd.DataFrame(metric_rows)
    hamming = pd.DataFrame(hamming_rows)
    epistasis = pd.DataFrame(epistasis_rows)
    if metrics.empty:
        return []
    rows: list[dict[str, Any]] = []
    zero = metrics[
        (metrics["score_layer"] == "observable_fitness")
        & (pd.to_numeric(metrics["epi_weight"], errors="coerce").abs() <= 1e-12)
    ]
    for parameter_set_id, sub in zero.groupby("parameter_set_id", dropna=False):
        support_mode = str(sub["support_mode"].iloc[0])
        seed = int(sub["seed"].iloc[0])
        hsub = hamming[
            (hamming["parameter_set_id"] == parameter_set_id)
            & (hamming["score_layer"] == "observable_fitness")
        ]
        esub = epistasis[
            (epistasis["parameter_set_id"] == parameter_set_id)
            & (epistasis["score_layer"] == "observable_fitness")
        ]
        beneficial_rho = spearman(
            [float(value) for value in hsub.get("fraction_of_peak", pd.Series(dtype=float)).tolist()],
            [float(value) for value in hsub.get("beneficial_fraction", pd.Series(dtype=float)).tolist()],
        )
        high_plateau = float(pd.to_numeric(sub.get("high_plateau_fraction", pd.Series(dtype=float)), errors="coerce").mean())
        low_plateau = float(pd.to_numeric(sub.get("low_plateau_fraction", pd.Series(dtype=float)), errors="coerce").mean())
        fold = esub[esub.get("threshold_band", pd.Series(dtype=str)) == "fold_transition"]
        deep = esub[esub.get("threshold_band", pd.Series(dtype=str)).isin(["deep_stable", "deep_unstable"])]
        fold_abs = float(pd.to_numeric(fold.get("mean_abs_pairwise_epistasis", pd.Series(dtype=float)), errors="coerce").mean())
        deep_abs = float(pd.to_numeric(deep.get("mean_abs_pairwise_epistasis", pd.Series(dtype=float)), errors="coerce").mean())
        rows.append(
            {
                "phase": "WS-B",
                "parameter_set_id": str(parameter_set_id),
                "support_mode": support_mode,
                "seed": seed,
                "epi_weight": 0.0,
                "n_local_maxima": float(pd.to_numeric(sub.get("n_local_maxima", pd.Series(dtype=float)), errors="coerce").mean()),
                "accessible_path_fraction": float(
                    pd.to_numeric(sub.get("accessible_path_fraction", pd.Series(dtype=float)), errors="coerce").mean()
                ),
                "beneficial_fraction_vs_fraction_of_peak_spearman": beneficial_rho,
                "beneficial_fraction_declines_toward_peak": bool(np.isfinite(beneficial_rho) and beneficial_rho < 0.0),
                "high_plateau_fraction": high_plateau,
                "low_plateau_fraction": low_plateau,
                "has_low_and_high_readout_plateaus": bool(high_plateau > 0.0 and low_plateau > 0.0),
                "fold_transition_mean_abs_pairwise_epistasis": fold_abs,
                "deep_regime_mean_abs_pairwise_epistasis": deep_abs,
                "transition_epistasis_enriched_over_deep_regimes": bool(
                    np.isfinite(fold_abs) and np.isfinite(deep_abs) and fold_abs > deep_abs
                ),
            }
        )
    return rows


def parameter_grid(config: dict[str, Any]) -> list[dict[str, Any]]:
    grid = []
    gen_cfg = config["generator"]
    for support_mode in gen_cfg["support_modes"]:
        for epi_weight in gen_cfg["epi_weight_grid"]:
            for seed in gen_cfg["seeds"]:
                parameter_set_id = f"{support_mode}__ew{float(epi_weight):.3f}__seed{int(seed)}"
                grid.append(
                    {
                        "phase": "WS-B",
                        "parameter_set_id": parameter_set_id,
                        "support_mode": str(support_mode),
                        "seed": int(seed),
                        "sequence_length": int(gen_cfg["sequence_length"]),
                        "n_terms": int(gen_cfg["n_terms"]),
                        "epi_weight": float(epi_weight),
                        "core_slope": float(gen_cfg["core_slope"]),
                        "stability_offset": float(gen_cfg["stability_offset"]),
                        "order_weights": {int(k): float(v) for k, v in dict(gen_cfg["order_weights"]).items()},
                    }
                )
    return grid


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=bool(args.quick))
    grid = parameter_grid(config)
    metric_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    hamming_rows: list[dict[str, Any]] = []
    epistasis_rows: list[dict[str, Any]] = []
    for parameter in grid:
        rows, terms, hamming, epistasis = evaluate_parameter_set(config, parameter)
        metric_rows.extend(rows)
        support_rows.extend(terms)
        hamming_rows.extend(hamming)
        epistasis_rows.extend(epistasis)
    knob_rows, monotonicity_rows = summarize_knobs(config, metric_rows)
    core_guardrail_rows = summarize_core_guardrails(metric_rows, hamming_rows, epistasis_rows)
    write_csv(output_dir / "parameter_grid.csv", grid)
    write_csv(output_dir / "graphlet_metrics.csv", metric_rows)
    write_csv(output_dir / "knob_summary.csv", knob_rows)
    write_csv(output_dir / "monotonicity_checks.csv", monotonicity_rows)
    write_csv(output_dir / "support_catalog.csv", support_rows)
    write_csv(output_dir / "hamming_layer_metrics.csv", hamming_rows)
    write_csv(output_dir / "threshold_epistasis_metrics.csv", epistasis_rows)
    write_csv(output_dir / "core_readout_guardrail_checks.csv", core_guardrail_rows)
    preflight = {
        "experiment_id": config["experiment"]["id"],
        "mode": config["experiment"]["mode"],
        "scratch_source": config["experiment"]["scientific_boundary"]["advisory_source"],
        "sequence_length": int(config["generator"]["sequence_length"]),
        "support_modes": [str(value) for value in config["generator"]["support_modes"]],
        "epi_weight_grid": [float(value) for value in config["generator"]["epi_weight_grid"]],
        "seeds": [int(value) for value in config["generator"]["seeds"]],
        "n_parameter_sets": int(len(grid)),
        "quick": bool(args.quick),
    }
    write_json(output_dir / "preflight.json", preflight)
    write_json(
        output_dir / "summary.json",
        {
            **preflight,
            "n_graphlet_metric_rows": int(len(metric_rows)),
            "n_support_terms": int(len(support_rows)),
            "n_hamming_layer_rows": int(len(hamming_rows)),
            "n_threshold_epistasis_rows": int(len(epistasis_rows)),
            "n_monotonicity_checks": int(len(monotonicity_rows)),
            "n_direction_pass": int(sum(1 for row in monotonicity_rows if row.get("direction_pass"))),
            "n_core_readout_guardrail_rows": int(len(core_guardrail_rows)),
            "n_core_readout_guardrails_with_declining_beneficial_fraction": int(
                sum(1 for row in core_guardrail_rows if row.get("beneficial_fraction_declines_toward_peak"))
            ),
        },
    )


if __name__ == "__main__":
    main()
