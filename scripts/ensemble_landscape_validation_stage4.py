#!/usr/bin/env python3
"""Stage 4 ensemble validation for the additive-core + GNK public protocol."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
from collections import Counter
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

from hashed_gnk_benchmark_knob_validation import rankdata, spearman  # noqa: E402
from phot_coupled_landscape_grid_search import load_config, write_csv, write_json  # noqa: E402
from phot_nested_fit_protocol_benchmark import epistasis_property_metrics  # noqa: E402
from project_adapt_env.dense_descriptive import exact_summary  # noqa: E402

from adaptenv import FitnessLandscapeEnv, LandscapeConfig, PUBLIC_MODEL_SYSTEM  # noqa: E402


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


def summarize_numeric(values: Iterable[float], prefix: str) -> dict[str, float]:
    arr = np.asarray([finite(v) for v in values], dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            f"{prefix}_n": 0.0,
            f"{prefix}_mean": math.nan,
            f"{prefix}_sd": math.nan,
            f"{prefix}_skew": math.nan,
            f"{prefix}_kurtosis": math.nan,
            f"{prefix}_q05": math.nan,
            f"{prefix}_q50": math.nan,
            f"{prefix}_q95": math.nan,
        }
    mean = float(np.mean(arr))
    sd = float(np.std(arr))
    centered = arr - mean
    skew = float(np.mean(centered**3) / (sd**3)) if sd > 1e-12 else math.nan
    kurt = float(np.mean(centered**4) / (sd**4) - 3.0) if sd > 1e-12 else math.nan
    return {
        f"{prefix}_n": float(arr.size),
        f"{prefix}_mean": mean,
        f"{prefix}_sd": sd,
        f"{prefix}_skew": skew,
        f"{prefix}_kurtosis": kurt,
        f"{prefix}_q05": float(np.quantile(arr, 0.05)),
        f"{prefix}_q50": float(np.quantile(arr, 0.50)),
        f"{prefix}_q95": float(np.quantile(arr, 0.95)),
    }


def descriptor_row(
    *,
    stage: str,
    dataset_id: str,
    descriptor_block: str,
    descriptor_name: str,
    value: Any,
    scope: str,
    n_observations: int,
    source: str = "empirical",
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    row = {
        "stage": stage,
        "source": source,
        "dataset_id": dataset_id,
        "descriptor_block": descriptor_block,
        "descriptor_name": descriptor_name,
        "value": finite(value),
        "scope": scope,
        "n_observations": int(n_observations),
    }
    if extra:
        row.update(extra)
    return row


def state_columns(frame: pd.DataFrame) -> list[str]:
    columns = [col for col in frame.columns if str(col).startswith("pos")]
    return sorted(columns, key=lambda col: int(str(col).removeprefix("pos") or 0))


def state_matrix(frame: pd.DataFrame) -> tuple[list[tuple[str, ...]], list[str], str]:
    pos_cols = state_columns(frame)
    if pos_cols:
        states = [tuple(str(v) for v in row) for row in frame[pos_cols].itertuples(index=False, name=None)]
        return states, pos_cols, "pos_columns"
    if "genotype" in frame:
        genotypes = frame["genotype"].dropna().astype(str)
        lengths = {len(value) for value in genotypes.head(1000)}
        if len(lengths) == 1 and next(iter(lengths)) <= 32 and not genotypes.str.contains(":").any():
            states = [tuple(str(ch) for ch in value) for value in frame["genotype"].astype(str)]
            return states, [f"pos{i + 1}" for i in range(next(iter(lengths)))], "genotype_chars"
    return [], [], "none"


def read_dataset_index(processed_root: Path) -> pd.DataFrame:
    index_path = processed_root / "dataset_index.csv"
    if not index_path.exists():
        raise FileNotFoundError(f"missing dataset index: {index_path}")
    return pd.read_csv(index_path)


def resolve_processed_root(data_cfg: dict[str, Any], project_root: Path) -> Path:
    configured = Path(str(data_cfg["processed_root"]))
    processed_root = configured if configured.is_absolute() else project_root / configured
    if processed_root.exists():
        return processed_root
    override = os.environ.get("STAGE4_PROCESSED_ROOT")
    if override:
        override_path = Path(override)
        if override_path.exists():
            return override_path
    return processed_root


def load_dataset(processed_root: Path, dataset_id: str) -> pd.DataFrame:
    path = processed_root / "datasets" / dataset_id / "variants.parquet"
    if not path.exists():
        raise FileNotFoundError(f"missing normalized dataset: {path}")
    frame = pd.read_parquet(path)
    if "fitness" not in frame:
        raise ValueError(f"{dataset_id} has no fitness column")
    frame = frame[np.isfinite(pd.to_numeric(frame["fitness"], errors="coerce"))].copy()
    if "fitness_is_higher_better" in frame and not bool(frame["fitness_is_higher_better"].iloc[0]):
        frame["fitness"] = -pd.to_numeric(frame["fitness"], errors="coerce")
    else:
        frame["fitness"] = pd.to_numeric(frame["fitness"], errors="coerce")
    return frame


def mutation_count_from_genotype(value: Any) -> int:
    text = str(value)
    if not text or text.lower() == "nan":
        return 0
    if ":" in text:
        return sum(1 for part in text.split(":") if part)
    return sum(1 for ch in text if ch not in {"0", ".", "-", "_"})


def symmetric_edge_effects(
    states: list[tuple[str, ...]],
    values: np.ndarray,
    *,
    max_edges: int | None = None,
) -> list[float]:
    if not states:
        return []
    state_to_value: dict[tuple[str, ...], float] = {}
    for state, value in zip(states, values):
        if state not in state_to_value:
            state_to_value[state] = float(value)
    alleles_by_site: list[list[str]] = []
    for site in range(len(states[0])):
        alleles_by_site.append(sorted({state[site] for state in state_to_value}))
    effects: list[float] = []
    seen: set[tuple[tuple[str, ...], tuple[str, ...]]] = set()
    for state, value in state_to_value.items():
        state_list = list(state)
        for site, alleles in enumerate(alleles_by_site):
            original = state_list[site]
            for allele in alleles:
                if allele == original:
                    continue
                state_list[site] = allele
                neighbor = tuple(state_list)
                state_list[site] = original
                if neighbor not in state_to_value:
                    continue
                edge = tuple(sorted((state, neighbor)))
                if edge in seen:
                    continue
                seen.add(edge)
                delta = float(state_to_value[neighbor] - value)
                effects.append(delta)
                effects.append(-delta)
                if max_edges is not None and len(effects) >= int(max_edges):
                    return effects
    return effects


def directed_edge_effects(states, values, *, reference_state, max_edges=None):
    """DFE of observed one-step mutations away from an explicitly supplied reference."""
    if reference_state not in states:
        raise ValueError("reference genotype is not observed")
    lookup = dict(zip(states, map(float, values)))
    effects = [value - lookup[reference_state] for state, value in lookup.items()
               if sum(a != b for a, b in zip(state, reference_state)) == 1]
    return effects if max_edges is None else effects[:int(max_edges)]


def add_order1_descriptors(
    *,
    rows: list[dict[str, Any]],
    stage: str,
    dataset_id: str,
    frame: pd.DataFrame,
    states: list[tuple[str, ...]],
    state_source: str,
    max_edges: int | None,
    reference_state: tuple[str, ...] | None = None,
) -> None:
    fitness = frame["fitness"].to_numpy(dtype=np.float64)
    fitness_stats = summarize_numeric(fitness, "fitness")
    for name, value in fitness_stats.items():
        rows.append(
            descriptor_row(
                stage=stage,
                dataset_id=dataset_id,
                descriptor_block="order1_or_marginal",
                descriptor_name=name,
                value=value,
                scope="fitness_distribution_all_observed_variants",
                n_observations=len(frame),
            )
        )
    if states and reference_state is not None:
        effects = directed_edge_effects(states, fitness, reference_state=reference_state, max_edges=max_edges)
        effect_scope = f"reference_conditioned_hamming1_from_{state_source}"
    elif states:
        effects = []
        effect_scope = "DFE_unavailable_without_explicit_reference"
    elif "genotype" in frame:
        counts = frame["genotype"].map(mutation_count_from_genotype)
        single = frame[counts == 1]
        effects = []
        effect_scope = "DFE_unavailable_without_measured_reference"
    else:
        effects = []
        effect_scope = "not_available"
    effect_stats = summarize_numeric(effects, "single_step_effect")
    neutral_threshold = max(1e-8, 0.05 * finite(effect_stats.get("single_step_effect_sd"), 0.0))
    arr = np.asarray([finite(v) for v in effects], dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    fractions = {
        "beneficial_fraction": float(np.mean(arr > neutral_threshold)) if arr.size else math.nan,
        "neutral_fraction": float(np.mean(np.abs(arr) <= neutral_threshold)) if arr.size else math.nan,
        "deleterious_fraction": float(np.mean(arr < -neutral_threshold)) if arr.size else math.nan,
        "neutral_threshold": neutral_threshold if arr.size else math.nan,
    }
    for name, value in {**effect_stats, **fractions}.items():
        rows.append(
            descriptor_row(
                stage=stage,
                dataset_id=dataset_id,
                descriptor_block="order1_or_marginal",
                descriptor_name=name,
                value=value,
                scope=effect_scope,
                n_observations=int(arr.size),
            )
        )


def accessible_paths_dp(values: np.ndarray, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    start = 0
    end = (1 << int(dim)) - 1
    if float(arr[end]) < float(arr[start]):
        start, end = end, start
    bits = [bit for bit in range(int(dim)) if ((start ^ end) >> bit) & 1]
    bit_count = len(bits)
    path_counts = {0: 1}
    for mask_size in range(bit_count):
        for mask in [m for m in range(1 << bit_count) if int(m).bit_count() == mask_size]:
            count = path_counts.get(mask, 0)
            if count == 0:
                continue
            node = start
            for offset, bit in enumerate(bits):
                if (mask >> offset) & 1:
                    node ^= 1 << bit
            current = float(arr[node])
            for offset, bit in enumerate(bits):
                if (mask >> offset) & 1:
                    continue
                nxt_mask = mask | (1 << offset)
                nxt = node ^ (1 << bit)
                if float(arr[nxt]) > current + 1e-12:
                    path_counts[nxt_mask] = path_counts.get(nxt_mask, 0) + count
    total = math.factorial(bit_count)
    accessible = path_counts.get((1 << bit_count) - 1, 0)
    return {
        "accessible_path_count_dp": float(accessible),
        "accessible_path_total_shortest_paths_dp": float(total),
        "accessible_path_fraction_dp": float(accessible / total) if total else math.nan,
    }


def metrics_from_cube(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    metrics = exact_summary(arr, dim=dim)
    metrics.update(epistasis_property_metrics(arr, dim=dim))
    metrics.update(accessible_paths_dp(arr, dim=dim))
    return {
        key: finite(value)
        for key, value in metrics.items()
        if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(finite(value))
    }


def graphlet_candidates(
    states: list[tuple[str, ...]],
    values: np.ndarray,
    *,
    dataset_id: str,
    max_dim: int,
    max_graphlets: int,
) -> list[dict[str, Any]]:
    if not states:
        return []
    state_to_value: dict[tuple[str, ...], float] = {
        state: float(value) for state, value in zip(states, values)
    }
    n_sites = len(states[0])
    alleles_by_site = [sorted({state[site] for state in states}) for site in range(n_sites)]
    rows: list[dict[str, Any]] = []

    if all(len(alleles) == 2 for alleles in alleles_by_site) and len(state_to_value) >= int(0.9 * (2**n_sites)):
        if len(state_to_value) == 2**n_sites and n_sites <= max(16, max_dim):
            ordered_values = []
            for node in range(1 << n_sites):
                state = tuple(alleles_by_site[site][(node >> site) & 1] for site in range(n_sites))
                ordered_values.append(state_to_value[state])
            metrics = metrics_from_cube(np.asarray(ordered_values), dim=n_sites)
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "graphlet_id": f"{dataset_id}__full_binary_cube",
                    "graphlet_scope": "full_binary_cube",
                    "dimension": n_sites,
                    "sites": ";".join(str(i) for i in range(n_sites)),
                    "alleles": "|".join("/".join(alleles_by_site[i]) for i in range(n_sites)),
                    **metrics,
                }
            )

    dim = min(int(max_dim), n_sites)
    if dim < 2:
        return rows
    attempts = 0
    seen: set[tuple[tuple[int, ...], tuple[tuple[str, str], ...]]] = set()
    site_combos = list(itertools.combinations(range(n_sites), dim))
    ranked_site_combos = sorted(site_combos, key=lambda combo: (sum(combo), combo))
    for sites_tuple in ranked_site_combos:
        allele_pairs_by_site = []
        for site in sites_tuple:
            alleles = alleles_by_site[site]
            if len(alleles) < 2:
                break
            pairs = list(itertools.combinations(alleles[: min(len(alleles), 8)], 2))
            allele_pairs_by_site.append(pairs[:8])
        else:
            for pair_tuple in itertools.product(*allele_pairs_by_site):
                key = (tuple(sites_tuple), tuple(tuple(pair) for pair in pair_tuple))
                if key in seen:
                    continue
                seen.add(key)
                cube_values = []
                complete = True
                for node in range(1 << dim):
                    template = [states[0][site] for site in range(n_sites)]
                    for offset, site in enumerate(sites_tuple):
                        template[site] = pair_tuple[offset][(node >> offset) & 1]
                    state = tuple(template)
                    if state not in state_to_value:
                        complete = False
                        break
                    cube_values.append(state_to_value[state])
                attempts += 1
                if not complete:
                    continue
                metrics = metrics_from_cube(np.asarray(cube_values), dim=dim)
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "graphlet_id": f"{dataset_id}__binary_subcube_{len(rows):04d}",
                        "graphlet_scope": "complete_binary_subcube",
                        "dimension": dim,
                        "sites": ";".join(str(site) for site in sites_tuple),
                        "alleles": "|".join("/".join(pair) for pair in pair_tuple),
                        **metrics,
                    }
                )
                if len(rows) >= int(max_graphlets):
                    return rows
                if attempts >= int(max_graphlets) * 100:
                    return rows
    return rows


def empirical_scope_audit_row(dataset_id: str, frame: pd.DataFrame, states: list[tuple[str, ...]], state_source: str) -> dict[str, Any]:
    tier = str(frame.get("tier", pd.Series(["unknown"])).iloc[0])
    has_state_graph = bool(states)
    mutation_counts = frame["genotype"].map(mutation_count_from_genotype) if "genotype" in frame else pd.Series(dtype=int)
    sparse_only = (not has_state_graph) and (not mutation_counts.empty)
    higher_order_allowed = has_state_graph and tier.startswith("tier1")
    return {
        "dataset_id": dataset_id,
        "tier": tier,
        "n_rows": int(len(frame)),
        "state_source": state_source,
        "n_state_sites": int(len(states[0])) if states else 0,
        "order1_allowed": True,
        "higher_order_allowed": bool(higher_order_allowed),
        "topography_allowed": bool(higher_order_allowed),
        "sparse_proteingym_higher_order_blocked": bool(sparse_only),
        "reason": (
            "tier1_state_graph_operator_free"
            if higher_order_allowed
            else "sparse_or_non_tier1_order1_only"
        ),
    }


def empirical_reference(
    config: dict[str, Any],
    project_root: Path,
    *,
    stage: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    data_cfg = config["data"]
    processed_root = resolve_processed_root(data_cfg, project_root)
    ids = list(data_cfg["dataset_ids"])
    descriptor_rows: list[dict[str, Any]] = []
    graphlet_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    for dataset_id in ids:
        frame = load_dataset(processed_root, str(dataset_id))
        states, _columns, state_source = state_matrix(frame)
        audit = empirical_scope_audit_row(str(dataset_id), frame, states, state_source)
        audit_rows.append(audit)
        add_order1_descriptors(
            rows=descriptor_rows,
            stage=stage,
            dataset_id=str(dataset_id),
            frame=frame,
            states=states,
            state_source=state_source,
            max_edges=data_cfg.get("max_hamming1_directed_effects"),
        )
        if audit["higher_order_allowed"]:
            graphlets = graphlet_candidates(
                states,
                frame["fitness"].to_numpy(dtype=np.float64),
                dataset_id=str(dataset_id),
                max_dim=int(data_cfg.get("binary_subcube_dim", 4)),
                max_graphlets=int(data_cfg.get("max_binary_subcubes_per_dataset", 64)),
            )
            graphlet_rows.extend(graphlets)
    return descriptor_rows, graphlet_rows, audit_rows


def empirical_bands(descriptor_rows: list[dict[str, Any]], graphlet_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    dframe = pd.DataFrame(descriptor_rows)
    if not dframe.empty:
        for (block, name), sub in dframe.groupby(["descriptor_block", "descriptor_name"], dropna=False):
            vals = pd.to_numeric(sub["value"], errors="coerce").dropna()
            if vals.empty:
                continue
            rows.append(
                {
                    "source_table": "descriptor_metrics",
                    "descriptor_block": block,
                    "descriptor_name": name,
                    "n_values": int(len(vals)),
                    "empirical_min": float(vals.min()),
                    "empirical_q05": float(vals.quantile(0.05)),
                    "empirical_median": float(vals.quantile(0.50)),
                    "empirical_q95": float(vals.quantile(0.95)),
                    "empirical_max": float(vals.max()),
                }
            )
    gframe = pd.DataFrame(graphlet_rows)
    metric_cols = [col for col in gframe.columns if col not in {"dataset_id", "graphlet_id", "graphlet_scope", "dimension", "sites", "alleles"}]
    for metric in metric_cols:
        vals = pd.to_numeric(gframe[metric], errors="coerce").dropna()
        if vals.empty:
            continue
        rows.append(
            {
                "source_table": "graphlet_metrics",
                "descriptor_block": "specific_epistasis_or_topography",
                "descriptor_name": metric,
                "n_values": int(len(vals)),
                "empirical_min": float(vals.min()),
                "empirical_q05": float(vals.quantile(0.05)),
                "empirical_median": float(vals.quantile(0.50)),
                "empirical_q95": float(vals.quantile(0.95)),
                "empirical_max": float(vals.max()),
            }
        )
    return rows


def build_landscape_config(config: dict[str, Any], parameter: dict[str, Any], *, epi_weight: float | None = None) -> LandscapeConfig:
    gen = config["generator"]
    order_weights = tuple(float(x) for x in gen.get("gnk_order_weights", [0.35, 0.45, 0.20]))
    return LandscapeConfig(
        model_system=PUBLIC_MODEL_SYSTEM,
        observed_fitness_combine_mode="stability_gate",
        L=int(gen["sequence_length"]),
        alphabet_size=20,
        seed=int(parameter["seed"]),
        stability_margin=float(parameter["stability_margin"]),
        stability_scale=float(gen.get("stability_scale", 1.0)),
        threshold_steepness=float(gen.get("threshold_steepness", 20.0)),
        fitness_output_scale=float(gen.get("fitness_output_scale", 10.0)),
        gnk_epistasis_weight=float(parameter["epi_weight"] if epi_weight is None else epi_weight),
        gnk_terms=int(parameter["gnk_terms"]),
        gnk_min_order=int(gen.get("gnk_min_order", 2)),
        gnk_max_order=int(gen.get("gnk_max_order", 4)),
        gnk_order_weights=order_weights,
        gnk_support_mode=str(parameter["support_mode"]),
        gnk_seed=int(parameter["seed"]) + int(gen.get("gnk_seed_offset", 10000)),
        n_epistatic_terms=0,
        noise_amplitude=0.0,
    )


def generator_parameter_grid(config: dict[str, Any]) -> list[dict[str, Any]]:
    gen = config["generator"]
    rows = []
    for support_mode in gen["support_modes"]:
        for seed in gen["seeds"]:
            for stability_margin in gen["stability_margin_grid"]:
                for epi_weight in gen["epi_weight_grid"]:
                    rows.append(
                        {
                            "parameter_set_id": (
                                f"{support_mode}__sm{float(stability_margin):.3f}"
                                f"__ew{float(epi_weight):.3f}__seed{int(seed)}"
                            ),
                            "support_mode": str(support_mode),
                            "seed": int(seed),
                            "sequence_length": int(gen["sequence_length"]),
                            "subcube_dim": int(gen["subcube_dim"]),
                            "stability_margin": float(stability_margin),
                            "epi_weight": float(epi_weight),
                            "gnk_terms": int(gen["gnk_terms"]),
                        }
                    )
    return rows


def subcube_sequences(length: int, dim: int, *, base: np.ndarray | None = None, alphabet_size: int = 20) -> np.ndarray:
    if base is None:
        base = np.zeros(int(length), dtype=np.int64)
    base = np.asarray(base, dtype=np.int64).copy()
    if base.size != int(length):
        raise ValueError(f"base sequence length {base.size} does not match L={length}")
    seqs = np.tile(base[None, :], (1 << int(dim), 1))
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            if (node >> bit) & 1:
                seqs[node, bit] = (int(base[bit]) + 1) % int(alphabet_size)
    return seqs


def generator_metrics(config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    parameters = generator_parameter_grid(config)
    for parameter in parameters:
        dim = int(parameter["subcube_dim"])
        full_env = FitnessLandscapeEnv(build_landscape_config(config, parameter))
        core_env = FitnessLandscapeEnv(build_landscape_config(config, parameter, epi_weight=0.0))
        base = getattr(full_env, "_designed_peak", None)
        seqs = subcube_sequences(
            int(parameter["sequence_length"]),
            dim,
            base=base,
            alphabet_size=int(getattr(full_env.cfg, "alphabet_size", 20)),
        )
        full_components = full_env.evaluate_batch_components(seqs)
        core_components = core_env.evaluate_batch_components(seqs)
        layers = {
            "full_public_observable": np.asarray(full_components["fitness"], dtype=np.float64),
            "additive_core_only_observable": np.asarray(core_components["fitness"], dtype=np.float64),
            "bare_gnk_texture": np.asarray(full_components["gnk_epistatic_score"], dtype=np.float64),
            "effective_stability_margin": np.asarray(full_components["effective_stability_margin_used"], dtype=np.float64),
        }
        for layer, values in layers.items():
            metrics = metrics_from_cube(values, dim=dim)
            row = {
                "source": "generator",
                "model_family": layer,
                **parameter,
                **metrics,
            }
            rows.append(row)
            if layer in {"additive_core_only_observable", "bare_gnk_texture"}:
                control_rows.append(row)
    return parameters, rows, control_rows


def generator_monotonicity(config: dict[str, Any], rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(rows)
    expectations = dict(config.get("expectations", {}))
    checks: list[dict[str, Any]] = []
    sub = frame[frame["model_family"] == "full_public_observable"]
    for keys, group in sub.groupby(["support_mode", "stability_margin"], dropna=False):
        support_mode, stability_margin = keys
        for metric, direction in expectations.items():
            if metric not in group:
                continue
            x = [finite(v) for v in group["epi_weight"].tolist()]
            y = [finite(v) for v in group[metric].tolist()]
            rho = spearman(x, y)
            sign = 1.0 if direction == "increase" else -1.0
            checks.append(
                {
                    "support_mode": support_mode,
                    "stability_margin": finite(stability_margin),
                    "knob": "epi_weight",
                    "metric": metric,
                    "expected_direction": direction,
                    "spearman_rho": rho,
                    "direction_pass": bool(np.isfinite(rho) and rho * sign > 0.0),
                    "n_parameter_sets": int(len(group)),
                }
            )
    return checks


def target_dataset_metric_rows(descriptor_rows: list[dict[str, Any]], graphlet_rows: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in descriptor_rows:
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "descriptor_name": row["descriptor_name"],
                "value": row["value"],
                "descriptor_block": row["descriptor_block"],
            }
        )
    for row in graphlet_rows:
        dataset_id = row["dataset_id"]
        for key, value in row.items():
            if key in {"dataset_id", "graphlet_id", "graphlet_scope", "dimension", "sites", "alleles"}:
                continue
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "descriptor_name": key,
                    "value": value,
                    "descriptor_block": "specific_epistasis_or_topography",
                }
            )
    return pd.DataFrame(rows)


def z_distance(value: float, target_values: pd.Series) -> float:
    vals = pd.to_numeric(target_values, errors="coerce").dropna().to_numpy(dtype=np.float64)
    if vals.size == 0 or not np.isfinite(value):
        return math.nan
    scale = max(float(np.std(vals)), abs(float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05))) / 3.29, 1e-9)
    center = float(np.median(vals))
    return abs(float(value) - center) / scale


def overlap_tables(
    config: dict[str, Any],
    descriptor_rows: list[dict[str, Any]],
    graphlet_rows: list[dict[str, Any]],
    generator_rows: list[dict[str, Any]],
    *,
    heldout_dataset: str | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    metrics = list(config["overlap"]["metrics"])
    target = target_dataset_metric_rows(descriptor_rows, graphlet_rows)
    if heldout_dataset is not None:
        target = target[target["dataset_id"] != heldout_dataset]
    gen = pd.DataFrame(generator_rows)
    coverage: list[dict[str, Any]] = []
    posterior: list[dict[str, Any]] = []
    for dataset_id, dsub in target.groupby("dataset_id", dropna=False):
        for gen_row in gen.to_dict(orient="records"):
            distances = []
            passes = 0
            considered = 0
            for metric in metrics:
                if metric not in gen_row:
                    continue
                vals = dsub[dsub["descriptor_name"] == metric]["value"]
                if vals.dropna().empty:
                    continue
                considered += 1
                value = finite(gen_row.get(metric))
                q05 = float(pd.to_numeric(vals, errors="coerce").quantile(0.05))
                q95 = float(pd.to_numeric(vals, errors="coerce").quantile(0.95))
                in_band = bool(np.isfinite(value) and q05 <= value <= q95)
                passes += int(in_band)
                distances.append(z_distance(value, vals))
                coverage.append(
                    {
                        "dataset_id": dataset_id,
                        "parameter_set_id": gen_row["parameter_set_id"],
                        "model_family": gen_row["model_family"],
                        "support_mode": gen_row["support_mode"],
                        "stability_margin": gen_row["stability_margin"],
                        "epi_weight": gen_row["epi_weight"],
                        "metric": metric,
                        "generator_value": value,
                        "empirical_q05": q05,
                        "empirical_q95": q95,
                        "within_empirical_band": in_band,
                    }
                )
            finite_dist = [d for d in distances if np.isfinite(d)]
            posterior.append(
                {
                    "dataset_id": dataset_id,
                    "parameter_set_id": gen_row["parameter_set_id"],
                    "model_family": gen_row["model_family"],
                    "support_mode": gen_row["support_mode"],
                    "stability_margin": gen_row["stability_margin"],
                    "epi_weight": gen_row["epi_weight"],
                    "n_metrics_considered": considered,
                    "n_metrics_in_band": passes,
                    "fraction_metrics_in_band": float(passes / considered) if considered else math.nan,
                    "mean_z_distance": float(np.mean(finite_dist)) if finite_dist else math.nan,
                }
            )
    return coverage, posterior


def run_stage_4a(config: dict[str, Any], project_root: Path, output_dir: Path) -> None:
    descriptors, graphlets, audit = empirical_reference(config, project_root, stage="4A")
    bands = empirical_bands(descriptors, graphlets)
    write_csv(output_dir / "descriptor_metrics.csv", descriptors)
    write_csv(output_dir / "graphlet_metrics.csv", graphlets)
    write_csv(output_dir / "empirical_bands.csv", bands)
    write_csv(output_dir / "descriptor_scope_audit.csv", audit)
    write_json(
        output_dir / "summary.json",
        {
            "stage": "4A",
            "n_datasets": len(config["data"]["dataset_ids"]),
            "n_descriptor_rows": len(descriptors),
            "n_graphlet_rows": len(graphlets),
            "n_scope_audit_rows": len(audit),
            "higher_order_dataset_ids": [row["dataset_id"] for row in audit if row["higher_order_allowed"]],
        },
    )


def run_stage_4b(config: dict[str, Any], output_dir: Path) -> None:
    parameters, rows, controls = generator_metrics(config)
    checks = generator_monotonicity(config, rows)
    write_csv(output_dir / "generator_parameter_grid.csv", parameters)
    write_csv(output_dir / "generator_descriptor_metrics.csv", rows)
    write_csv(output_dir / "generator_knob_monotonicity.csv", checks)
    write_csv(output_dir / "core_texture_ablation_controls.csv", controls)
    write_json(
        output_dir / "summary.json",
        {
            "stage": "4B",
            "public_model_system": PUBLIC_MODEL_SYSTEM,
            "n_parameter_sets": len(parameters),
            "n_generator_metric_rows": len(rows),
            "n_monotonicity_checks": len(checks),
            "n_direction_pass": sum(1 for row in checks if row.get("direction_pass")),
        },
    )


def run_stage_4c(config: dict[str, Any], project_root: Path, output_dir: Path) -> None:
    descriptors, graphlets, audit = empirical_reference(config, project_root, stage="4C")
    parameters, generator_rows, controls = generator_metrics(config)
    coverage, posterior = overlap_tables(config, descriptors, graphlets, generator_rows)
    write_csv(output_dir / "empirical_bands.csv", empirical_bands(descriptors, graphlets))
    write_csv(output_dir / "coverage_metrics.csv", coverage)
    write_csv(output_dir / "posterior_knob_regions.csv", posterior)
    write_csv(output_dir / "control_overlap_failures.csv", [row for row in posterior if row["model_family"] != "full_public_observable"])
    write_csv(output_dir / "descriptor_scope_audit.csv", audit)
    write_json(
        output_dir / "summary.json",
        {
            "stage": "4C",
            "n_parameter_sets": len(parameters),
            "n_coverage_rows": len(coverage),
            "n_posterior_rows": len(posterior),
        },
    )


def run_stage_4d(config: dict[str, Any], project_root: Path, output_dir: Path) -> None:
    descriptors, graphlets, audit = empirical_reference(config, project_root, stage="4D")
    _parameters, generator_rows, _controls = generator_metrics(config)
    coverage, posterior = overlap_tables(config, descriptors, graphlets, generator_rows)
    post = pd.DataFrame(posterior)
    min_fraction = float(config["joint_realism"].get("min_fraction_metrics_in_band", 0.6))
    min_assays = int(config["joint_realism"].get("min_assays_hit", 2))
    joint_rows: list[dict[str, Any]] = []
    if not post.empty:
        eligible = post[
            (post["model_family"] == "full_public_observable")
            & (pd.to_numeric(post["fraction_metrics_in_band"], errors="coerce") >= min_fraction)
        ]
        for parameter_set_id, sub in eligible.groupby("parameter_set_id", dropna=False):
            first = sub.iloc[0].to_dict()
            joint_rows.append(
                {
                    "parameter_set_id": parameter_set_id,
                    "support_mode": first["support_mode"],
                    "stability_margin": first["stability_margin"],
                    "epi_weight": first["epi_weight"],
                    "n_assays_hit": int(sub["dataset_id"].nunique()),
                    "joint_realism_pass": bool(sub["dataset_id"].nunique() >= min_assays),
                    "hit_dataset_ids": ";".join(sorted(str(x) for x in sub["dataset_id"].unique())),
                }
            )
    write_csv(output_dir / "joint_realism_coverage.csv", coverage)
    write_csv(output_dir / "joint_realism_candidates.csv", joint_rows)
    write_csv(output_dir / "posterior_knob_regions.csv", posterior)
    write_csv(output_dir / "descriptor_scope_audit.csv", audit)
    write_json(
        output_dir / "summary.json",
        {
            "stage": "4D",
            "n_joint_candidate_rows": len(joint_rows),
            "n_joint_pass": sum(1 for row in joint_rows if row["joint_realism_pass"]),
            "min_assays_hit": min_assays,
            "min_fraction_metrics_in_band": min_fraction,
        },
    )


def run_stage_4e(config: dict[str, Any], project_root: Path, output_dir: Path) -> None:
    descriptors, graphlets, audit = empirical_reference(config, project_root, stage="4E")
    _parameters, generator_rows, _controls = generator_metrics(config)
    target_ids = [str(x) for x in config["data"]["dataset_ids"] if str(x) in {row["dataset_id"] for row in audit if row["higher_order_allowed"]}]
    prediction_rows: list[dict[str, Any]] = []
    for heldout in target_ids:
        _coverage_train, posterior_train = overlap_tables(
            config,
            descriptors,
            graphlets,
            generator_rows,
            heldout_dataset=heldout,
        )
        train = pd.DataFrame(posterior_train)
        train = train[train["model_family"] == "full_public_observable"].copy()
        if train.empty:
            continue
        train["rank_score"] = pd.to_numeric(train["mean_z_distance"], errors="coerce")
        selected = train.sort_values(["rank_score", "fraction_metrics_in_band"], ascending=[True, False]).iloc[0]
        heldout_config = dict(config)
        heldout_config["data"] = dict(config["data"])
        heldout_config["data"]["dataset_ids"] = [heldout]
        h_desc, h_graph, _h_audit = empirical_reference(heldout_config, project_root, stage="4E-heldout")
        coverage, posterior = overlap_tables(config, h_desc, h_graph, generator_rows)
        post = pd.DataFrame(posterior)
        match = post[
            (post["parameter_set_id"] == selected["parameter_set_id"])
            & (post["model_family"] == "full_public_observable")
            & (post["dataset_id"] == heldout)
        ]
        if match.empty:
            continue
        row = match.iloc[0].to_dict()
        prediction_rows.append(
            {
                "heldout_dataset_id": heldout,
                "selected_parameter_set_id": selected["parameter_set_id"],
                "selected_support_mode": selected["support_mode"],
                "selected_stability_margin": selected["stability_margin"],
                "selected_epi_weight": selected["epi_weight"],
                "train_mean_z_distance": selected["mean_z_distance"],
                "heldout_mean_z_distance": row["mean_z_distance"],
                "heldout_fraction_metrics_in_band": row["fraction_metrics_in_band"],
                "heldout_n_metrics_in_band": row["n_metrics_in_band"],
                "heldout_n_metrics_considered": row["n_metrics_considered"],
            }
        )
    write_csv(output_dir / "heldout_predictions.csv", prediction_rows)
    write_csv(output_dir / "descriptor_scope_audit.csv", audit)
    write_json(
        output_dir / "summary.json",
        {
            "stage": "4E",
            "n_heldout_datasets": len(target_ids),
            "n_prediction_rows": len(prediction_rows),
            "median_heldout_fraction_metrics_in_band": (
                float(np.median([row["heldout_fraction_metrics_in_band"] for row in prediction_rows]))
                if prediction_rows
                else math.nan
            ),
        },
    )


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=bool(args.quick))
    stage = str(config["experiment"]["stage"])
    write_json(
        output_dir / "preflight.json",
        {
            "experiment_id": config["experiment"]["id"],
            "stage": stage,
            "mode": config["experiment"]["mode"],
            "quick": bool(args.quick),
            "public_model_system": PUBLIC_MODEL_SYSTEM,
            "project_root": str(project_root),
        },
    )
    if stage == "4A":
        run_stage_4a(config, project_root, output_dir)
    elif stage == "4B":
        run_stage_4b(config, output_dir)
    elif stage == "4C":
        run_stage_4c(config, project_root, output_dir)
    elif stage == "4D":
        run_stage_4d(config, project_root, output_dir)
    elif stage == "4E":
        run_stage_4e(config, project_root, output_dir)
    else:
        raise ValueError(f"unsupported Stage 4 experiment: {stage}")


if __name__ == "__main__":
    main()
