#!/usr/bin/env python3
"""EXP-057 design-matched, support-aware Stage 4B probe.

This runner is intentionally separate from the completed Stage 4A-4E runner.
It fixes the measurement defects identified in
scratch/RES-054_assessment_and_path_to_finishing_component1.md section 3:

- sample graphlets on instantiated GNK supports, not positions 0..dim-1;
- match empirical design dimension and alphabet cardinality per Tier-1 assay;
- stratify anchors across unfolded, fold-transition, and plateau regimes;
- compute plateau/noise-robust topography with a shared neutrality tolerance;
- gate on public-composite versus additive-core separation before overlap work.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
EXTERNAL_ADAPTENV = PROJECT_ROOT / "external" / "Adapt-Env"
for path in (SRC_ROOT, SCRIPT_DIR, EXTERNAL_ADAPTENV):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ensemble_landscape_validation_stage4 import build_landscape_config, generator_parameter_grid  # noqa: E402
from hashed_gnk_benchmark_knob_validation import rankdata, spearman  # noqa: E402
from phot_coupled_landscape_grid_search import load_config, write_csv, write_json  # noqa: E402

from adaptenv import FitnessLandscapeEnv, PUBLIC_MODEL_SYSTEM  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def finite(value: Any, default: float = math.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def stable_uint64(*parts: Any, seed: int = 0) -> int:
    payload = json.dumps([seed, *parts], sort_keys=True, separators=(",", ":"), default=str).encode()
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little")


def rng_for(*parts: Any, seed: int = 0) -> np.random.Generator:
    return np.random.default_rng(stable_uint64(*parts, seed=seed) % (2**32))


def state_columns(frame: pd.DataFrame) -> list[str]:
    columns = [col for col in frame.columns if str(col).startswith("pos")]
    return sorted(columns, key=lambda col: int(str(col).removeprefix("pos") or 0))


def resolve_processed_root(data_cfg: dict[str, Any], project_root: Path) -> Path:
    configured = Path(str(data_cfg["processed_root"]))
    processed_root = configured if configured.is_absolute() else project_root / configured
    if processed_root.exists():
        return processed_root
    override = os.environ.get("STAGE4_PROCESSED_ROOT")
    if override and Path(override).exists():
        return Path(override)
    return processed_root


def load_dataset(processed_root: Path, dataset_id: str) -> pd.DataFrame:
    path = processed_root / "datasets" / dataset_id / "variants.parquet"
    if not path.exists():
        raise FileNotFoundError(f"missing normalized dataset: {path}")
    frame = pd.read_parquet(path)
    if "fitness" not in frame:
        raise ValueError(f"{dataset_id} has no fitness column")
    frame = frame[np.isfinite(pd.to_numeric(frame["fitness"], errors="coerce"))].copy()
    frame["fitness"] = pd.to_numeric(frame["fitness"], errors="coerce")
    if "fitness_is_higher_better" in frame and not bool(frame["fitness_is_higher_better"].iloc[0]):
        frame["fitness"] = -frame["fitness"]
    return frame


def dataset_template(frame: pd.DataFrame, dataset_id: str) -> dict[str, Any]:
    pos_cols = state_columns(frame)
    if pos_cols:
        states = [tuple(str(v) for v in row) for row in frame[pos_cols].itertuples(index=False, name=None)]
        per_site = [sorted({state[i] for state in states}) for i in range(len(pos_cols))]
        state_source = "pos_columns"
    elif "genotype" in frame:
        genotypes = frame["genotype"].dropna().astype(str)
        lengths = {len(value) for value in genotypes.head(1000)}
        if len(lengths) != 1:
            raise ValueError(f"{dataset_id} has no rectangular genotype template")
        dim = next(iter(lengths))
        states = [tuple(str(ch) for ch in value) for value in frame["genotype"].astype(str)]
        per_site = [sorted({state[i] for state in states}) for i in range(dim)]
        state_source = "genotype_chars"
    else:
        raise ValueError(f"{dataset_id} has no design coordinates")

    cardinalities = [len(values) for values in per_site]
    n_expected = int(np.prod(cardinalities, dtype=np.int64)) if cardinalities else 0
    unique_states = len(set(states))
    tier = str(frame.get("tier", pd.Series(["unknown"])).iloc[0])
    alphabet = str(frame.get("alphabet", pd.Series(["unknown"])).iloc[0])
    return {
        "dataset_id": dataset_id,
        "tier": tier,
        "landscape_name": str(frame.get("landscape_name", pd.Series([dataset_id])).iloc[0]),
        "readout": str(frame.get("readout", pd.Series(["unknown"])).iloc[0]),
        "state_source": state_source,
        "dimension": int(len(cardinalities)),
        "alphabet_cardinalities": cardinalities,
        "max_alphabet_cardinality": int(max(cardinalities)) if cardinalities else 0,
        "n_expected_nodes": n_expected,
        "n_observed_unique_states": unique_states,
        "observed_state_fraction": float(unique_states / n_expected) if n_expected else math.nan,
        "complete_factorial": bool(unique_states == n_expected),
        "alphabet": alphabet,
    }


def load_templates(config: dict[str, Any], project_root: Path) -> list[dict[str, Any]]:
    data_cfg = config["data"]
    processed_root = resolve_processed_root(data_cfg, project_root)
    templates: list[dict[str, Any]] = []
    for dataset_id in data_cfg["dataset_ids"]:
        frame = load_dataset(processed_root, str(dataset_id))
        template = dataset_template(frame, str(dataset_id))
        if template["dimension"] < 2:
            continue
        templates.append(template)
    return templates


def complete_empirical_values(frame: pd.DataFrame, dataset_id: str) -> tuple[np.ndarray, tuple[int, ...], str, int]:
    """Return empirical values in complete factorial order, or raise on missing states."""
    pos_cols = state_columns(frame)
    if pos_cols:
        states = [tuple(str(v) for v in row) for row in frame[pos_cols].itertuples(index=False, name=None)]
        per_site = [sorted({state[i] for state in states}) for i in range(len(pos_cols))]
        state_source = "pos_columns"
    elif "genotype" in frame:
        genotypes = frame["genotype"].dropna().astype(str)
        lengths = {len(value) for value in genotypes}
        if len(lengths) != 1:
            raise ValueError(f"{dataset_id} has no rectangular genotype template")
        dim = next(iter(lengths))
        states = [tuple(str(ch) for ch in value) for value in frame["genotype"].astype(str)]
        per_site = [sorted({state[i] for state in states}) for i in range(dim)]
        state_source = "genotype_chars"
    else:
        raise ValueError(f"{dataset_id} has no design coordinates")

    shape = tuple(len(values) for values in per_site)
    expected = int(np.prod(shape, dtype=np.int64)) if shape else 0
    values_by_state: dict[tuple[str, ...], list[float]] = defaultdict(list)
    for state, value in zip(states, frame["fitness"].to_numpy(dtype=np.float64)):
        if np.isfinite(value):
            values_by_state[state].append(float(value))
    missing = expected - len(values_by_state)
    if missing != 0:
        raise ValueError(f"{dataset_id} is not complete ({missing} states missing from template)")
    ordered = [
        float(np.mean(values_by_state[tuple(state)]))
        for state in itertools.product(*per_site)
    ]
    return np.asarray(ordered, dtype=np.float64), shape, state_source, expected


def empirical_template_metrics(config: dict[str, Any], project_root: Path) -> list[dict[str, Any]]:
    data_cfg = config["data"]
    processed_root = resolve_processed_root(data_cfg, project_root)
    rows: list[dict[str, Any]] = []
    for dataset_id in data_cfg["dataset_ids"]:
        frame = load_dataset(processed_root, str(dataset_id))
        template = dataset_template(frame, str(dataset_id))
        base = {
            "source": "empirical",
            "model_family": "empirical_observed",
            "dataset_id": str(dataset_id),
            "landscape_name": template["landscape_name"],
            "readout": template["readout"],
            "tier": template["tier"],
            "template_dimension": int(template["dimension"]),
            "template_max_alphabet_cardinality": int(template["max_alphabet_cardinality"]),
            "node_count": int(template["n_observed_unique_states"]),
            "expected_node_count": int(template["n_expected_nodes"]),
            "observed_state_fraction": finite(template["observed_state_fraction"]),
            "complete_factorial": bool(template["complete_factorial"]),
        }
        try:
            values, shape, state_source, expected = complete_empirical_values(frame, str(dataset_id))
        except ValueError as exc:
            rows.append(
                base
                | {
                    "descriptor_status": "skipped_incomplete_template",
                    "state_source": template["state_source"],
                    "skip_reason": str(exc),
                }
            )
            continue
        rows.append(
            base
            | {
                "descriptor_status": "computed",
                "state_source": state_source,
                "node_count": int(expected),
            }
            | robust_metrics(
                values,
                shape,
                config["descriptors"],
                seed=stable_uint64("empirical", dataset_id) % (2**32),
            )
        )
    return rows


def support_catalog(env: FitnessLandscapeEnv, parameter_set_id: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    term_id = 0
    for order, count in sorted(getattr(env, "_gnk_order_counts", {}).items()):
        for _ in range(int(count)):
            support = tuple(int(pos) for pos in env._gnk_support_for_term(term_id, int(order)))
            rows.append(
                {
                    "parameter_set_id": parameter_set_id,
                    "term_id": int(term_id),
                    "order": int(order),
                    "support_positions": ";".join(str(pos) for pos in support),
                    "support_size": int(len(support)),
                }
            )
            term_id += 1
    return rows


def choose_support_sites(
    terms: list[dict[str, Any]],
    *,
    dimension: int,
    graphlet_index: int,
    sequence_length: int,
    seed: int,
) -> tuple[int, ...]:
    if not terms:
        ranked_all = sorted(
            range(sequence_length),
            key=lambda pos: stable_uint64("fill-no-terms", graphlet_index, pos, seed=seed),
        )
        return tuple(sorted(ranked_all[:dimension]))

    ranked_terms = sorted(
        terms,
        key=lambda row: stable_uint64("term-choice", graphlet_index, row["term_id"], row["support_positions"], seed=seed),
    )
    selected: set[int] = set()
    start = ranked_terms[graphlet_index % len(ranked_terms)]
    queue = [start] + ranked_terms
    for row in queue:
        for pos in str(row["support_positions"]).split(";"):
            if pos != "":
                selected.add(int(pos))
        if len(selected) >= dimension:
            break

    if len(selected) < dimension:
        support_positions = sorted(
            {
                int(pos)
                for row in terms
                for pos in str(row["support_positions"]).split(";")
                if pos != ""
            },
            key=lambda pos: stable_uint64("support-fill", graphlet_index, pos, seed=seed),
        )
        for pos in support_positions:
            selected.add(int(pos))
            if len(selected) >= dimension:
                break

    if len(selected) < dimension:
        ranked_all = sorted(
            range(sequence_length),
            key=lambda pos: stable_uint64("global-fill", graphlet_index, pos, seed=seed),
        )
        for pos in ranked_all:
            selected.add(int(pos))
            if len(selected) >= dimension:
                break

    ranked_selected = sorted(
        selected,
        key=lambda pos: stable_uint64("site-order", graphlet_index, pos, seed=seed),
    )[:dimension]
    return tuple(sorted(ranked_selected))


def allele_choices(anchor: np.ndarray, sites: tuple[int, ...], cardinalities: list[int], alphabet_size: int, seed: int) -> list[list[int]]:
    out: list[list[int]] = []
    for offset, site in enumerate(sites):
        k = min(int(cardinalities[offset]), int(alphabet_size))
        anchor_aa = int(anchor[int(site)])
        if k >= alphabet_size:
            choices = list(range(alphabet_size))
        else:
            ranked = sorted(
                [aa for aa in range(alphabet_size) if aa != anchor_aa],
                key=lambda aa: stable_uint64("allele-choice", offset, site, aa, seed=seed),
            )
            choices = [anchor_aa] + ranked[: max(k - 1, 0)]
        out.append([int(x) for x in choices])
    return out


def factorial_sequences(anchor: np.ndarray, sites: tuple[int, ...], choices: list[list[int]]) -> tuple[np.ndarray, tuple[int, ...]]:
    shape = tuple(len(values) for values in choices)
    n_nodes = int(np.prod(shape, dtype=np.int64))
    seqs = np.tile(np.asarray(anchor, dtype=np.int64)[None, :], (n_nodes, 1))
    if n_nodes:
        coordinates = np.indices(shape, dtype=np.int16).reshape(len(shape), -1).T
        for axis, site in enumerate(sites):
            alleles = np.asarray(choices[axis], dtype=np.int64)
            seqs[:, int(site)] = alleles[coordinates[:, axis]]
    return seqs, shape


def evaluate_components(env: FitnessLandscapeEnv, seqs: np.ndarray, *, chunk_size: int) -> dict[str, np.ndarray]:
    chunks: dict[str, list[np.ndarray]] = defaultdict(list)
    for start in range(0, seqs.shape[0], int(chunk_size)):
        comp = env.evaluate_batch_components(seqs[start : start + int(chunk_size)])
        for key in ("fitness", "gnk_epistatic_score", "effective_stability_margin_used"):
            chunks[key].append(np.asarray(comp[key], dtype=np.float64))
    return {key: np.concatenate(parts) for key, parts in chunks.items()}


def neutrality_tolerance(values: np.ndarray, cfg: dict[str, Any]) -> float:
    arr = np.asarray(values, dtype=np.float64)
    sd = float(np.std(arr)) if arr.size else 0.0
    return max(float(cfg.get("abs_tolerance", 1e-8)), float(cfg.get("sd_fraction", 0.02)) * sd)


def single_step_effects(values: np.ndarray, shape: tuple[int, ...], *, n_samples: int, seed: int) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64).reshape(shape)
    rng = rng_for("single-step", seed=seed)
    effects: list[float] = []
    dim = len(shape)
    if dim == 0:
        return np.empty(0, dtype=np.float64)
    total_edges = int(np.prod(shape, dtype=np.int64)) * sum(k - 1 for k in shape)
    n = min(int(n_samples), max(total_edges, 0))
    for _ in range(n):
        idx = [int(rng.integers(0, k)) for k in shape]
        axis_choices = [axis for axis, k in enumerate(shape) if k > 1]
        axis = int(axis_choices[int(rng.integers(0, len(axis_choices)))])
        new = int(rng.integers(0, shape[axis] - 1))
        if new >= idx[axis]:
            new += 1
        idx2 = list(idx)
        idx2[axis] = new
        effects.append(float(arr[tuple(idx2)] - arr[tuple(idx)]))
    return np.asarray(effects, dtype=np.float64)


def high_order_fraction(values: np.ndarray, shape: tuple[int, ...]) -> float:
    tensor = np.asarray(values, dtype=np.float64).reshape(shape)
    centered = tensor - float(np.mean(tensor))
    total = float(np.sum(centered * centered))
    if total <= 1e-18:
        return math.nan
    mean = float(np.mean(tensor))
    main_ss = 0.0
    for axis in range(len(shape)):
        site_mean = np.mean(tensor, axis=tuple(ax for ax in range(len(shape)) if ax != axis))
        effect = site_mean - mean
        repeats = int(np.prod([shape[ax] for ax in range(len(shape)) if ax != axis], dtype=np.int64))
        main_ss += float(repeats * np.sum(effect * effect))
    return float(max(0.0, min(1.0, 1.0 - main_ss / total)))


def strict_local_maxima(values: np.ndarray, shape: tuple[int, ...], tol: float) -> tuple[int, float]:
    tensor = np.asarray(values, dtype=np.float64).reshape(shape)
    is_max = np.ones(shape, dtype=bool)
    for axis, k in enumerate(shape):
        if k <= 1:
            continue
        moved = np.moveaxis(tensor, axis, 0)
        moved_mask = np.moveaxis(is_max, axis, 0)
        if k == 2:
            moved_mask[0] &= moved[0] > moved[1] + float(tol)
            moved_mask[1] &= moved[1] > moved[0] + float(tol)
            continue
        top1_index = np.argmax(moved, axis=0)
        top1 = np.max(moved, axis=0)
        top2 = np.partition(moved, -2, axis=0)[-2]
        for allele in range(k):
            neighbor_max = np.where(top1_index == allele, top2, top1)
            moved_mask[allele] &= moved[allele] > neighbor_max + float(tol)
    count = int(np.sum(is_max))
    return count, float(count / max(tensor.size, 1))


def unravel_index(index: int, shape: tuple[int, ...]) -> list[int]:
    out: list[int] = []
    remainder = int(index)
    for k in shape:
        out.append(remainder % int(k))
        remainder //= int(k)
    return out


def ravel_index(index: Iterable[int], shape: tuple[int, ...]) -> int:
    stride = 1
    out = 0
    for value, k in zip(index, shape):
        out += int(value) * stride
        stride *= int(k)
    return int(out)


def peak_targeted_reachability(values: np.ndarray, shape: tuple[int, ...], tol: float, *, n_starts: int, seed: int) -> float:
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        return math.nan
    best = int(np.argmax(arr))
    rng = rng_for("reachability", seed=seed)
    starts = list(rng.choice(arr.size, size=min(int(n_starts), arr.size), replace=False))
    if best not in starts:
        starts.append(best)
    successes = 0
    for start in starts:
        current = int(start)
        for _ in range(max(4 * len(shape) * max(shape), 1)):
            idx = unravel_index(current, shape)
            best_neighbor = current
            best_value = float(arr[current])
            for axis, k in enumerate(shape):
                for allele in range(k):
                    if allele == idx[axis]:
                        continue
                    idx2 = list(idx)
                    idx2[axis] = allele
                    neighbor = ravel_index(idx2, shape)
                    value = float(arr[neighbor])
                    if value > best_value + float(tol):
                        best_value = value
                        best_neighbor = neighbor
            if best_neighbor == current:
                break
            current = best_neighbor
        successes += int(current == best)
    return float(successes / max(len(starts), 1))


def random_walk_autocorr(values: np.ndarray, shape: tuple[int, ...], *, n_steps: int, seed: int) -> float:
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    if arr.size < 3:
        return math.nan
    rng = rng_for("random-walk", seed=seed)
    current = int(rng.integers(0, arr.size))
    series: list[float] = []
    for _ in range(int(n_steps)):
        series.append(float(arr[current]))
        idx = unravel_index(current, shape)
        axes = [axis for axis, k in enumerate(shape) if k > 1]
        axis = int(axes[int(rng.integers(0, len(axes)))])
        allele = int(rng.integers(0, shape[axis] - 1))
        if allele >= idx[axis]:
            allele += 1
        idx[axis] = allele
        current = ravel_index(idx, shape)
    x = np.asarray(series[:-1], dtype=np.float64)
    y = np.asarray(series[1:], dtype=np.float64)
    if x.size < 2 or float(np.std(x)) <= 1e-12 or float(np.std(y)) <= 1e-12:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def pairwise_epistasis_sample(values: np.ndarray, shape: tuple[int, ...], effects: np.ndarray, *, n_samples: int, seed: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64).reshape(shape)
    rng = rng_for("pairwise-epistasis", seed=seed)
    dim = len(shape)
    if dim < 2:
        return {
            "sign_epistasis_fraction": math.nan,
            "reciprocal_sign_epistasis_fraction": math.nan,
            "mean_abs_pairwise_epistasis": math.nan,
            "magnitude_epistasis_to_single_step_sd_ratio": math.nan,
        }
    epsilons: list[float] = []
    sign_flags: list[float] = []
    reciprocal_flags: list[float] = []
    axes_with_choices = [axis for axis, k in enumerate(shape) if k > 1]
    if len(axes_with_choices) < 2:
        return {
            "sign_epistasis_fraction": math.nan,
            "reciprocal_sign_epistasis_fraction": math.nan,
            "mean_abs_pairwise_epistasis": math.nan,
            "magnitude_epistasis_to_single_step_sd_ratio": math.nan,
        }
    for _ in range(int(n_samples)):
        left, right = rng.choice(axes_with_choices, size=2, replace=False)
        idx = [int(rng.integers(0, k)) for k in shape]
        left_alt = int(rng.integers(0, shape[left] - 1))
        if left_alt >= idx[left]:
            left_alt += 1
        right_alt = int(rng.integers(0, shape[right] - 1))
        if right_alt >= idx[right]:
            right_alt += 1
        idx00 = list(idx)
        idx10 = list(idx)
        idx01 = list(idx)
        idx11 = list(idx)
        idx10[left] = left_alt
        idx01[right] = right_alt
        idx11[left] = left_alt
        idx11[right] = right_alt
        f00 = float(arr[tuple(idx00)])
        f10 = float(arr[tuple(idx10)])
        f01 = float(arr[tuple(idx01)])
        f11 = float(arr[tuple(idx11)])
        left_low = f10 - f00
        left_high = f11 - f01
        right_low = f01 - f00
        right_high = f11 - f10
        eps = f11 - f10 - f01 + f00
        epsilons.append(eps)
        left_flip = left_low * left_high < -1e-24
        right_flip = right_low * right_high < -1e-24
        sign_flags.append(float(left_flip or right_flip))
        reciprocal_flags.append(float(left_flip and right_flip))
    eps_arr = np.asarray(epsilons, dtype=np.float64)
    effect_sd = float(np.std(effects)) if effects.size else math.nan
    return {
        "sign_epistasis_fraction": float(np.mean(sign_flags)) if sign_flags else math.nan,
        "reciprocal_sign_epistasis_fraction": float(np.mean(reciprocal_flags)) if reciprocal_flags else math.nan,
        "mean_abs_pairwise_epistasis": float(np.mean(np.abs(eps_arr))) if eps_arr.size else math.nan,
        "magnitude_epistasis_to_single_step_sd_ratio": (
            float(np.mean(np.abs(eps_arr)) / effect_sd) if np.isfinite(effect_sd) and effect_sd > 1e-12 and eps_arr.size else math.nan
        ),
    }


def robust_metrics(values: np.ndarray, shape: tuple[int, ...], cfg: dict[str, Any], *, seed: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    tol = neutrality_tolerance(arr, cfg.get("neutrality", {}))
    effects = single_step_effects(
        arr,
        shape,
        n_samples=int(cfg.get("single_step_effect_samples", 50000)),
        seed=seed,
    )
    finite_effects = effects[np.isfinite(effects)]
    max_count, max_density = strict_local_maxima(arr, shape, tol)
    epistasis = pairwise_epistasis_sample(
        arr,
        shape,
        finite_effects,
        n_samples=int(cfg.get("pairwise_epistasis_samples", 20000)),
        seed=seed,
    )
    if finite_effects.size:
        beneficial = float(np.mean(finite_effects > tol))
        neutral = float(np.mean(np.abs(finite_effects) <= tol))
        deleterious = float(np.mean(finite_effects < -tol))
        effect_sd = float(np.std(finite_effects))
    else:
        beneficial = neutral = deleterious = effect_sd = math.nan
    high_order = high_order_fraction(arr, shape)
    return {
        "n_nodes": float(arr.size),
        "fitness_mean": float(np.mean(arr)) if arr.size else math.nan,
        "fitness_sd": float(np.std(arr)) if arr.size else math.nan,
        "single_step_effect_sd": effect_sd,
        "beneficial_fraction": beneficial,
        "neutral_fraction": neutral,
        "deleterious_fraction": deleterious,
        "neutrality_tolerance": float(tol),
        "n_local_maxima_strict": float(max_count),
        "local_maxima_density": max_density,
        "peak_targeted_reachability_fraction": peak_targeted_reachability(
            arr,
            shape,
            tol,
            n_starts=int(cfg.get("reachability_starts", 512)),
            seed=seed,
        ),
        "random_walk_lag1_autocorrelation": random_walk_autocorr(
            arr,
            shape,
            n_steps=int(cfg.get("random_walk_steps", 4096)),
            seed=seed,
        ),
        "high_order_fraction": high_order,
        "walsh_high_order_fraction": high_order,
        **epistasis,
    }


def sample_anchor_rows(env: FitnessLandscapeEnv, config: dict[str, Any], parameter: dict[str, Any]) -> list[dict[str, Any]]:
    probe_cfg = config["probe"]
    n_candidates = int(probe_cfg.get("anchor_candidate_samples", 4096))
    seed = int(parameter["seed"]) + int(probe_cfg.get("anchor_seed_offset", 41000))
    rng = rng_for("anchors", parameter["parameter_set_id"], seed=seed)
    candidates = [getattr(env, "_designed_peak").copy(), env.reference.copy()]
    candidates.extend(env.random_sequences(max(n_candidates // 2, 1)))
    peak = getattr(env, "_designed_peak").copy()
    for _ in range(max(n_candidates - len(candidates), 0)):
        seq = peak.copy()
        n_mut = int(rng.integers(1, max(2, min(env.L, 24))))
        positions = rng.choice(env.L, size=min(n_mut, env.L), replace=False)
        for pos in positions:
            aa = int(rng.integers(0, env.A - 1))
            if aa >= seq[int(pos)]:
                aa += 1
            seq[int(pos)] = aa
        candidates.append(seq)
    seqs = np.asarray(candidates, dtype=np.int64)
    margins = evaluate_components(env, seqs, chunk_size=int(probe_cfg.get("eval_chunk_size", 65536)))["effective_stability_margin_used"]
    rows: list[dict[str, Any]] = []
    for regime in probe_cfg["anchor_regimes"]:
        name = str(regime["name"])
        target = float(regime["target_margin"])
        idx = int(np.argmin(np.abs(margins - target)))
        anchor = seqs[idx].copy()
        rows.append(
            {
                "anchor_regime": name,
                "target_effective_stability_margin": target,
                "actual_effective_stability_margin": float(margins[idx]),
                "distance_from_designed_peak": int(np.count_nonzero(anchor != peak)),
                "anchor_sequence": ";".join(str(int(x)) for x in anchor),
                "anchor_array": anchor,
            }
        )
    return rows


def safe_checkpoint_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._=-" else "_" for ch in str(value))


def write_csv_atomic(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    pd.DataFrame(rows).to_csv(tmp, index=False)
    tmp.replace(path)


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    tmp.replace(path)


def read_csv_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return pd.read_csv(path).to_dict(orient="records")


def checkpoint_root(output_dir: Path) -> Path:
    return output_dir / "checkpoints" / "generator_parameters"


def checkpoint_source_roots(output_dir: Path) -> list[Path]:
    roots: list[Path] = []
    for env_name in ("LABPROJ_CHECKPOINT_PATH", "LABPROJ_CHECKPOINT_DIR"):
        value = os.environ.get(env_name)
        if not value:
            continue
        path = Path(value)
        candidates = [
            path,
            path / "generator_parameters",
            path / "checkpoints" / "generator_parameters",
        ]
        if path.is_file():
            candidates.append(path.parent)
        for candidate in candidates:
            if candidate.is_dir() and candidate.name == "generator_parameters" and candidate not in roots:
                roots.append(candidate)
    current = checkpoint_root(output_dir)
    if current.is_dir() and current not in roots:
        roots.append(current)
    return roots


def parameter_checkpoint_dir(root: Path, parameter_set_id: str) -> Path:
    return root / safe_checkpoint_name(parameter_set_id)


def load_parameter_checkpoint(param_dir: Path) -> dict[str, list[dict[str, Any]]] | None:
    done = param_dir / "done.json"
    if not done.exists():
        return None
    return {
        "metric_rows": read_csv_rows(param_dir / "metric_rows.csv"),
        "control_rows": read_csv_rows(param_dir / "control_rows.csv"),
        "support_rows": read_csv_rows(param_dir / "support_rows.csv"),
        "graphlet_rows": read_csv_rows(param_dir / "graphlet_rows.csv"),
        "anchor_audit_rows": read_csv_rows(param_dir / "anchor_audit_rows.csv"),
    }


def copy_resume_checkpoints(output_dir: Path) -> list[str]:
    current = checkpoint_root(output_dir)
    current.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for source_root in checkpoint_source_roots(output_dir):
        if source_root.resolve() == current.resolve():
            continue
        for source_dir in sorted(source_root.iterdir()):
            if not source_dir.is_dir() or not (source_dir / "done.json").exists():
                continue
            destination = current / source_dir.name
            if destination.exists():
                continue
            shutil.copytree(source_dir, destination)
            copied.append(source_dir.name)
    return copied


def existing_parameter_checkpoints(output_dir: Path) -> dict[str, dict[str, list[dict[str, Any]]]]:
    completed: dict[str, dict[str, list[dict[str, Any]]]] = {}
    root = checkpoint_root(output_dir)
    if not root.exists():
        return completed
    for param_dir in sorted(root.iterdir()):
        if not param_dir.is_dir():
            continue
        payload = load_parameter_checkpoint(param_dir)
        if payload is None:
            continue
        metadata_path = param_dir / "metadata.json"
        if metadata_path.exists():
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            parameter_set_id = str(metadata.get("parameter_set_id") or param_dir.name)
        else:
            parameter_set_id = param_dir.name
        completed[parameter_set_id] = payload
    return completed


def write_parameter_checkpoint(
    output_dir: Path,
    parameter: dict[str, Any],
    *,
    metric_rows: list[dict[str, Any]],
    control_rows: list[dict[str, Any]],
    support_rows: list[dict[str, Any]],
    graphlet_rows: list[dict[str, Any]],
    anchor_audit_rows: list[dict[str, Any]],
) -> None:
    param_dir = parameter_checkpoint_dir(checkpoint_root(output_dir), str(parameter["parameter_set_id"]))
    param_dir.mkdir(parents=True, exist_ok=True)
    done = param_dir / "done.json"
    if done.exists():
        done.unlink()
    write_json_atomic(
        param_dir / "metadata.json",
        {
            "parameter_set_id": parameter["parameter_set_id"],
            "parameter": parameter,
            "status": "complete",
        },
    )
    write_csv_atomic(param_dir / "metric_rows.csv", metric_rows)
    write_csv_atomic(param_dir / "control_rows.csv", control_rows)
    write_csv_atomic(param_dir / "support_rows.csv", support_rows)
    write_csv_atomic(param_dir / "graphlet_rows.csv", graphlet_rows)
    write_csv_atomic(param_dir / "anchor_audit_rows.csv", anchor_audit_rows)
    write_json_atomic(
        done,
        {
            "parameter_set_id": parameter["parameter_set_id"],
            "n_metric_rows": len(metric_rows),
            "n_control_rows": len(control_rows),
            "n_support_rows": len(support_rows),
            "n_graphlet_rows": len(graphlet_rows),
            "n_anchor_audit_rows": len(anchor_audit_rows),
        },
    )


def append_checkpoint_payload(
    payload: dict[str, list[dict[str, Any]]],
    *,
    metric_rows: list[dict[str, Any]],
    control_rows: list[dict[str, Any]],
    support_rows: list[dict[str, Any]],
    graphlet_rows: list[dict[str, Any]],
    anchor_audit_rows: list[dict[str, Any]],
) -> None:
    metric_rows.extend(payload["metric_rows"])
    control_rows.extend(payload["control_rows"])
    support_rows.extend(payload["support_rows"])
    graphlet_rows.extend(payload["graphlet_rows"])
    anchor_audit_rows.extend(payload["anchor_audit_rows"])


def generator_metrics(
    config: dict[str, Any],
    *,
    output_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    templates = config["_templates"]
    probe_cfg = config["probe"]
    parameters = generator_parameter_grid(config)
    copied_checkpoints = copy_resume_checkpoints(output_dir)
    completed_checkpoints = existing_parameter_checkpoints(output_dir)
    metric_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    graphlet_rows: list[dict[str, Any]] = []
    anchor_audit_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    resumed_parameter_sets: list[str] = []
    computed_parameter_sets: list[str] = []

    for index, parameter in enumerate(parameters, start=1):
        parameter_set_id = str(parameter["parameter_set_id"])
        if parameter_set_id in completed_checkpoints:
            append_checkpoint_payload(
                completed_checkpoints[parameter_set_id],
                metric_rows=metric_rows,
                control_rows=control_rows,
                support_rows=support_rows,
                graphlet_rows=graphlet_rows,
                anchor_audit_rows=anchor_audit_rows,
            )
            resumed_parameter_sets.append(parameter_set_id)
            print(f"[EXP-057] resume {index}/{len(parameters)} {parameter_set_id}", flush=True)
            continue

        print(f"[EXP-057] compute {index}/{len(parameters)} {parameter_set_id}", flush=True)
        parameter_metric_rows: list[dict[str, Any]] = []
        parameter_control_rows: list[dict[str, Any]] = []
        parameter_support_rows: list[dict[str, Any]] = []
        parameter_graphlet_rows: list[dict[str, Any]] = []
        parameter_anchor_audit_rows: list[dict[str, Any]] = []
        full_env = FitnessLandscapeEnv(build_landscape_config(config, parameter))
        core_env = FitnessLandscapeEnv(build_landscape_config(config, parameter, epi_weight=0.0))
        terms = support_catalog(full_env, str(parameter["parameter_set_id"]))
        parameter_support_rows.extend(terms)
        anchor_rows = sample_anchor_rows(full_env, config, parameter)
        for anchor_row in anchor_rows:
            parameter_anchor_audit_rows.append(
                {
                    key: value
                    for key, value in anchor_row.items()
                    if key != "anchor_array"
                }
                | {
                    "parameter_set_id": parameter["parameter_set_id"],
                    "support_mode": parameter["support_mode"],
                    "stability_margin": parameter["stability_margin"],
                    "epi_weight": parameter["epi_weight"],
                }
            )

        for template in templates:
            dim = int(template["dimension"])
            graphlet_count = int(probe_cfg.get("support_graphlets_per_template", 4))
            if int(template["max_alphabet_cardinality"]) >= 10:
                graphlet_count = int(probe_cfg.get("potts_graphlets_per_template", graphlet_count))
            for anchor_row in anchor_rows:
                anchor = np.asarray(anchor_row["anchor_array"], dtype=np.int64)
                for graphlet_index in range(graphlet_count):
                    sites = choose_support_sites(
                        terms,
                        dimension=dim,
                        graphlet_index=graphlet_index,
                        sequence_length=int(parameter["sequence_length"]),
                        seed=int(parameter["seed"]) + int(probe_cfg.get("site_seed_offset", 51000)),
                    )
                    choices = allele_choices(
                        anchor,
                        sites,
                        [int(x) for x in template["alphabet_cardinalities"]],
                        alphabet_size=int(full_env.cfg.alphabet_size),
                        seed=int(parameter["seed"]) + graphlet_index,
                    )
                    seqs, shape = factorial_sequences(anchor, sites, choices)
                    full_components = evaluate_components(
                        full_env,
                        seqs,
                        chunk_size=int(probe_cfg.get("eval_chunk_size", 65536)),
                    )
                    core_components = evaluate_components(
                        core_env,
                        seqs,
                        chunk_size=int(probe_cfg.get("eval_chunk_size", 65536)),
                    )
                    layers = {
                        "full_public_observable": full_components["fitness"],
                        "additive_core_only_observable": core_components["fitness"],
                        "bare_gnk_texture": full_components["gnk_epistatic_score"],
                        "effective_stability_margin": full_components["effective_stability_margin_used"],
                    }
                    graphlet_id = (
                        f"{parameter['parameter_set_id']}__{template['dataset_id']}"
                        f"__{anchor_row['anchor_regime']}__g{graphlet_index:03d}"
                    )
                    parameter_graphlet_rows.append(
                        {
                            "graphlet_id": graphlet_id,
                            "parameter_set_id": parameter["parameter_set_id"],
                            "dataset_id": template["dataset_id"],
                            "anchor_regime": anchor_row["anchor_regime"],
                            "graphlet_index": int(graphlet_index),
                            "dimension": dim,
                            "alphabet_cardinalities": ";".join(str(k) for k in shape),
                            "n_nodes": int(seqs.shape[0]),
                            "varied_sites": ";".join(str(pos) for pos in sites),
                            "support_aware": True,
                        }
                    )
                    for layer, values in layers.items():
                        metrics = robust_metrics(
                            values,
                            shape,
                            config["descriptors"],
                            seed=stable_uint64(graphlet_id, layer, seed=int(parameter["seed"])) % (2**32),
                        )
                        row = {
                            "source": "generator",
                            "model_family": layer,
                            "graphlet_id": graphlet_id,
                            "dataset_id": template["dataset_id"],
                            "anchor_regime": anchor_row["anchor_regime"],
                            "graphlet_index": int(graphlet_index),
                            "template_dimension": dim,
                            "template_max_alphabet_cardinality": int(template["max_alphabet_cardinality"]),
                            "node_count": int(seqs.shape[0]),
                            "varied_sites": ";".join(str(pos) for pos in sites),
                            **parameter,
                            **metrics,
                        }
                        parameter_metric_rows.append(row)
                        if layer in {"additive_core_only_observable", "bare_gnk_texture"}:
                            parameter_control_rows.append(row)

        write_parameter_checkpoint(
            output_dir,
            parameter,
            metric_rows=parameter_metric_rows,
            control_rows=parameter_control_rows,
            support_rows=parameter_support_rows,
            graphlet_rows=parameter_graphlet_rows,
            anchor_audit_rows=parameter_anchor_audit_rows,
        )
        metric_rows.extend(parameter_metric_rows)
        control_rows.extend(parameter_control_rows)
        support_rows.extend(parameter_support_rows)
        graphlet_rows.extend(parameter_graphlet_rows)
        anchor_audit_rows.extend(parameter_anchor_audit_rows)
        computed_parameter_sets.append(parameter_set_id)

    resume_summary = {
        "checkpoint_root": str(checkpoint_root(output_dir)),
        "copied_checkpoint_dirs": copied_checkpoints,
        "n_resumed_parameter_sets": len(resumed_parameter_sets),
        "n_computed_parameter_sets": len(computed_parameter_sets),
        "resumed_parameter_sets": resumed_parameter_sets,
        "computed_parameter_sets": computed_parameter_sets,
    }
    return parameters, metric_rows, control_rows, support_rows, graphlet_rows + anchor_audit_rows, resume_summary


def monotonicity_checks(config: dict[str, Any], rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(rows)
    sub = frame[frame["model_family"] == "full_public_observable"].copy()
    expectations = dict(config["expectations"])
    rho_floor = float(config["gates"].get("monotonicity_abs_rho_floor", 0.5))
    checks: list[dict[str, Any]] = []
    group_cols = ["dataset_id", "support_mode", "stability_margin", "anchor_regime"]
    for keys, group in sub.groupby(group_cols, dropna=False):
        dataset_id, support_mode, stability_margin, anchor_regime = keys
        collapsed = group.groupby("epi_weight", as_index=False).mean(numeric_only=True)
        for metric, direction in expectations.items():
            if metric not in collapsed:
                continue
            x = [finite(v) for v in collapsed["epi_weight"].tolist()]
            y = [finite(v) for v in collapsed[metric].tolist()]
            rho = spearman(x, y)
            sign = 1.0 if str(direction) == "increase" else -1.0
            checks.append(
                {
                    "dataset_id": dataset_id,
                    "support_mode": support_mode,
                    "stability_margin": finite(stability_margin),
                    "anchor_regime": anchor_regime,
                    "knob": "epi_weight",
                    "metric": metric,
                    "expected_direction": direction,
                    "spearman_rho": rho,
                    "abs_rho_floor": rho_floor,
                    "direction_pass": bool(np.isfinite(rho) and rho * sign >= rho_floor),
                    "n_epi_weight_levels": int(len(collapsed)),
                }
            )
    return checks


def separation_gates(config: dict[str, Any], rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(rows)
    metrics = list(config["gates"]["core_gnk_separation_metrics"])
    thresholds = {
        str(item["metric"]): float(item["min_abs_delta"])
        for item in config["gates"].get("core_gnk_effect_thresholds", [])
    }
    rows_out: list[dict[str, Any]] = []
    key_cols = ["dataset_id", "support_mode", "stability_margin", "anchor_regime", "epi_weight"]
    if frame.empty:
        return rows_out
    for keys, group in frame.groupby(key_cols, dropna=False):
        pass_flags: list[bool] = []
        out = dict(zip(key_cols, keys))
        for metric in metrics:
            full = finite(group.loc[group["model_family"] == "full_public_observable", metric].mean())
            core = finite(group.loc[group["model_family"] == "additive_core_only_observable", metric].mean())
            delta = full - core if np.isfinite(full) and np.isfinite(core) else math.nan
            threshold = float(thresholds.get(metric, 0.0))
            out[f"{metric}_full_minus_core"] = delta
            out[f"{metric}_min_abs_delta"] = threshold
            out[f"{metric}_separation_pass"] = bool(np.isfinite(delta) and abs(delta) >= threshold)
            pass_flags.append(bool(np.isfinite(delta) and abs(delta) >= threshold))
        out["core_gnk_separation_pass"] = bool(any(pass_flags) and finite(out.get("epi_weight"), 0.0) > 0.0)
        rows_out.append(out)
    return rows_out


def summarize(
    config: dict[str, Any],
    parameters: list[dict[str, Any]],
    metric_rows: list[dict[str, Any]],
    empirical_rows: list[dict[str, Any]],
    checks: list[dict[str, Any]],
    gates: list[dict[str, Any]],
    resume_summary: dict[str, Any],
) -> dict[str, Any]:
    check_pass = sum(1 for row in checks if row.get("direction_pass"))
    gate_positive = [row for row in gates if finite(row.get("epi_weight"), 0.0) > 0.0]
    gate_pass = sum(1 for row in gate_positive if row.get("core_gnk_separation_pass"))
    empirical_computed = sum(1 for row in empirical_rows if row.get("descriptor_status") == "computed")
    return {
        "experiment_id": config["experiment"]["id"],
        "stage": config["experiment"]["stage"],
        "mode": config["experiment"]["mode"],
        "public_model_system": PUBLIC_MODEL_SYSTEM,
        "n_templates": len(config["_templates"]),
        "n_empirical_template_descriptor_rows": len(empirical_rows),
        "n_empirical_template_descriptors_computed": int(empirical_computed),
        "n_parameter_sets": len(parameters),
        "n_generator_metric_rows": len(metric_rows),
        "n_monotonicity_checks": len(checks),
        "n_monotonicity_direction_pass": int(check_pass),
        "monotonicity_pass_fraction": float(check_pass / len(checks)) if checks else math.nan,
        "n_positive_core_gnk_gate_rows": len(gate_positive),
        "n_positive_core_gnk_gate_pass": int(gate_pass),
        "core_gnk_gate_pass_fraction": float(gate_pass / len(gate_positive)) if gate_positive else math.nan,
        "checkpoint_root": resume_summary.get("checkpoint_root"),
        "n_resumed_parameter_sets": int(resume_summary.get("n_resumed_parameter_sets", 0)),
        "n_computed_parameter_sets": int(resume_summary.get("n_computed_parameter_sets", 0)),
        "stop_rule": "do_not_run_overlap_until_monotone_knobs_and_core_gnk_separation_pass",
    }


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=bool(args.quick))
    config["_templates"] = load_templates(config, project_root)
    write_json(
        output_dir / "preflight.json",
        {
            "experiment_id": config["experiment"]["id"],
            "stage": config["experiment"]["stage"],
            "mode": config["experiment"]["mode"],
            "quick": bool(args.quick),
            "public_model_system": PUBLIC_MODEL_SYSTEM,
            "project_root": str(project_root),
            "scratch_design_source": "scratch/RES-054_assessment_and_path_to_finishing_component1.md#3",
        },
    )
    write_csv(output_dir / "empirical_design_templates.csv", config["_templates"])
    empirical_rows = empirical_template_metrics(config, project_root)
    write_csv(output_dir / "empirical_template_descriptor_metrics.csv", empirical_rows)
    parameters, metric_rows, controls, support_rows, probe_rows, resume_summary = generator_metrics(
        config,
        output_dir=output_dir,
    )
    checks = monotonicity_checks(config, metric_rows)
    gates = separation_gates(config, metric_rows)
    write_csv(output_dir / "generator_parameter_grid.csv", parameters)
    write_csv(output_dir / "support_catalog.csv", support_rows)
    write_csv(output_dir / "support_aware_graphlet_audit.csv", [row for row in probe_rows if "graphlet_id" in row])
    write_csv(output_dir / "anchor_regime_audit.csv", [row for row in probe_rows if "anchor_sequence" in row])
    write_csv(output_dir / "generator_descriptor_metrics.csv", metric_rows)
    write_csv(output_dir / "generator_knob_monotonicity.csv", checks)
    write_csv(output_dir / "core_texture_ablation_controls.csv", controls)
    write_csv(output_dir / "core_gnk_separation_gates.csv", gates)
    write_json(output_dir / "checkpoint_resume_summary.json", resume_summary)
    write_json(
        output_dir / "summary.json",
        summarize(config, parameters, metric_rows, empirical_rows, checks, gates, resume_summary),
    )


if __name__ == "__main__":
    main()
