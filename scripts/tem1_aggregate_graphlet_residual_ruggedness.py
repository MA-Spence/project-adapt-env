#!/usr/bin/env python3
"""Run EXP-025 aggregate graphlet residual-ruggedness diagnostics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = Path(__file__).resolve().parent
for path in (ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from adaptenv import AlignmentProfile, FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from adaptenv.amino_acids import AA_ORDER  # noqa: E402

from adaptive_walk_peak_epistasis_baselines import (  # noqa: E402
    enumerate_single_mutants,
    evaluate_batch,
    json_default,
    merge_dict,
    single_mutant_scan,
    write_json,
)
from tem1_dca_clique_haplotype_epistasis import (  # noqa: E402
    dca_graph,
    enumerate_cliques,
    local_extrema_counts,
    marginal_frequency,
    mutation_name,
    sign_epistasis_counts,
    walsh_decomposition,
)
from tem1_homologue_walk_parameter_grid import (  # noqa: E402
    build_alignment_profile,
    components_for_seq,
    consensus_metrics,
    filled_unique_homologues,
    profile_log_odds_bits,
    seq_to_string,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not quick:
        raw.pop("quick", None)
        return raw

    base_models = list(raw.get("candidate_models", []))
    quick_cfg = dict(raw.get("quick", {}))
    merged = merge_dict(raw, quick_cfg)
    include = quick_cfg.get("candidate_models", {}).get("include")
    if include is not None:
        include_set = {str(model_id) for model_id in include}
        merged["candidate_models"] = [
            model for model in base_models if str(model.get("id")) in include_set
        ]
    merged.pop("quick", None)
    return merged


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    preferred_fields: list[str] | None = None,
) -> None:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for key in preferred_fields or []:
        if key not in seen:
            seen.add(key)
            fieldnames.append(key)
    for row in rows:
        for key in row:
            if key.startswith("_"):
                continue
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    if not fieldnames:
        fieldnames = ["empty"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: csv_value(row.get(key)) for key in fieldnames})


def csv_value(value: Any) -> Any:
    if isinstance(value, (list, tuple, np.ndarray)):
        return ";".join(str(csv_value(item)) for item in list(value))
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if value is None:
        return ""
    return value


def stable_seed(*parts: Any) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:8], 16)


def stable_fraction(*parts: Any) -> float:
    return stable_seed(*parts) / float(0xFFFFFFFF)


def choose_evenly(values: list[Any], count: int) -> list[Any]:
    if count <= 0:
        return []
    if len(values) <= count:
        return list(values)
    idx = np.unique(np.rint(np.linspace(0, len(values) - 1, count)).astype(np.int64))
    while idx.size < count:
        for candidate in range(len(values)):
            if candidate not in set(int(x) for x in idx):
                idx = np.append(idx, candidate)
                break
    return [values[int(i)] for i in sorted(idx[:count])]


def rank_quantiles(rows: list[dict[str, Any]], source_key: str, target_key: str) -> None:
    values = np.asarray([float(row.get(source_key, math.nan)) for row in rows], dtype=np.float64)
    valid = np.isfinite(values)
    if np.count_nonzero(valid) == 0:
        for row in rows:
            row[target_key] = math.nan
        return
    order = np.argsort(values[valid], kind="mergesort")
    ranks = np.empty(np.count_nonzero(valid), dtype=np.float64)
    ranks[order] = np.linspace(0.0, 1.0, np.count_nonzero(valid))
    cursor = 0
    for row, ok in zip(rows, valid):
        if ok:
            row[target_key] = float(ranks[cursor])
            cursor += 1
        else:
            row[target_key] = math.nan


def base_landscape_config(cfg: dict[str, Any], *, seed: int) -> LandscapeConfig:
    base = dict(cfg["adapt_env"]["base_config"])
    base.update(
        {
            "L": len(str(cfg["reference"]["sequence"])),
            "seed": int(seed),
            "n_epistatic_terms": 0,
            "epistasis_strength": 0.0,
            "epistasis_covariation_bias": 0.0,
            "n_higher_order_epistatic_terms": 0,
            "noise_amplitude": 0.0,
            "noise_heterogeneity": 0.0,
            "generic_epistasis_target": "stability",
            "empirical_pairwise_strength": 0.0,
        }
    )
    return LandscapeConfig(**base)


def build_base_landscape(
    cfg: dict[str, Any],
    *,
    seed: int,
    profile: AlignmentProfile,
) -> FitnessLandscapeEnv:
    return FitnessLandscapeEnv(
        base_landscape_config(cfg, seed=seed),
        alignment_profile=profile,
        reference_sequence=str(cfg["reference"]["sequence"]),
    )


def scan_config(cfg: dict[str, Any]) -> dict[str, Any]:
    return {
        "batch_size": int(cfg.get("batch_size", 4096)),
        "beneficial_threshold_abs": 1e-9,
        "beneficial_threshold_fraction_of_parent": 0.0,
        "neutral_threshold_abs": 1e-9,
        "neutral_threshold_fraction_of_parent": 0.0,
    }


def local_scan_summary(
    landscape: FitnessLandscapeEnv,
    seq: np.ndarray,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    scan = single_mutant_scan(landscape, seq, sampling_cfg=scan_config(cfg))
    deltas = np.asarray(scan["deltas"], dtype=np.float64)
    return {
        "background_fitness": float(scan["parent_fitness"]),
        "local_beneficial_fraction": float(scan["fraction_beneficial"]),
        "local_neutral_fraction": float(scan["fraction_neutral"]),
        "local_deleterious_fraction": float(scan["fraction_deleterious"]),
        "local_mean_delta": float(np.mean(deltas)) if deltas.size else math.nan,
        "local_sd_delta": float(np.std(deltas)) if deltas.size else math.nan,
        "local_best_delta": float(scan["best_delta"]),
        "local_mean_beneficial_effect": float(scan["mean_beneficial_effect"]),
    }


def start_candidate_rows(
    profile: AlignmentProfile,
    landscape: FitnessLandscapeEnv,
) -> list[dict[str, Any]]:
    seqs, source_indices = filled_unique_homologues(profile)
    reference = np.asarray(landscape.reference, dtype=np.int64)
    consensus = np.asarray(profile.consensus_idx, dtype=np.int64)
    log_odds = profile_log_odds_bits(profile)
    positions = np.arange(profile.L, dtype=np.int64)
    consensus_score = float(np.sum(log_odds[positions, consensus]))

    candidates: list[tuple[str, int | None, np.ndarray]] = [
        ("reference", None, reference),
        ("alignment_consensus", None, consensus),
    ]
    for unique_idx, seq in enumerate(seqs):
        candidates.append(("mmseqs_homologue", int(source_indices[unique_idx]), seq))

    rows: list[dict[str, Any]] = []
    seen: set[tuple[int, ...]] = set()
    for source, source_index, seq in candidates:
        key = tuple(int(x) for x in seq)
        if key in seen:
            continue
        seen.add(key)
        comp = components_for_seq(landscape, seq)
        metrics = consensus_metrics(
            seq,
            consensus=consensus,
            log_odds_bits=log_odds,
            consensus_score=consensus_score,
        )
        rows.append(
            {
                "source": source,
                "source_alignment_index": source_index,
                "fitness": float(landscape.evaluate(seq)),
                "sequence": seq_to_string(seq),
                "_seq": np.asarray(seq, dtype=np.int64).copy(),
                **comp,
                **metrics,
            }
        )
    return rows


def pick_rows(
    rows: list[dict[str, Any]],
    *,
    count: int,
    used_sequences: set[str],
    key,
) -> list[dict[str, Any]]:
    ordered = sorted(rows, key=key)
    out: list[dict[str, Any]] = []
    for row in ordered:
        seq = str(row["sequence"])
        if seq in used_sequences:
            continue
        out.append(row)
        used_sequences.add(seq)
        if len(out) >= count:
            break
    if len(out) < count:
        selected_ids = {id(row) for row in out}
        for row in ordered:
            seq = str(row["sequence"])
            if id(row) in selected_ids:
                continue
            out.append(row)
            used_sequences.add(seq)
            selected_ids.add(id(row))
            if len(out) >= count:
                break
    return out


def select_starts(
    cfg: dict[str, Any],
    *,
    profile: AlignmentProfile,
    landscape: FitnessLandscapeEnv,
    landscape_seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows = start_candidate_rows(profile, landscape)
    starts_per_class = int(cfg["start_selection"]["starts_per_class"])
    used_sequences: set[str] = set()
    selected: list[dict[str, Any]] = []
    fitnesses = np.asarray([float(row["fitness"]) for row in rows], dtype=np.float64)
    distances = np.asarray([float(row["consensus_distance_bits"]) for row in rows], dtype=np.float64)
    close_cut = float(np.quantile(distances, 0.35))
    far_cut = float(np.quantile(distances, 0.75))
    median_fit = float(np.median(fitnesses))

    close_rows = [row for row in rows if float(row["consensus_distance_bits"]) <= close_cut]
    far_rows = [row for row in rows if float(row["consensus_distance_bits"]) >= far_cut]
    viable_rows = [
        row
        for row in rows
        if row["source"] not in {"reference", "alignment_consensus"}
        and float(row["fitness"]) >= median_fit
    ]
    classes: list[tuple[str, list[dict[str, Any]], Any]] = [
        (
            "reference_or_consensus_like",
            rows,
            lambda row: (
                0 if row["source"] in {"reference", "alignment_consensus"} else 1,
                -float(row["consensus_identity"]),
                -float(row["fitness"]),
            ),
        ),
        (
            "high_fitness_close_homologue",
            close_rows or rows,
            lambda row: (-float(row["fitness"]), float(row["consensus_distance_bits"])),
        ),
        (
            "marginal_stability_homologue",
            rows,
            lambda row: (
                abs(float(row["effective_stability_margin_used"]) - 1.0),
                -float(row["consensus_distance_bits"]),
            ),
        ),
        (
            "distant_low_fitness_homologue",
            far_rows or rows,
            lambda row: (float(row["fitness"]), -float(row["consensus_distance_bits"])),
        ),
        (
            "viable_off_consensus",
            viable_rows or rows,
            lambda row: (
                -float(row["fitness"]),
                -float(row["consensus_distance_bits"]),
            ),
        ),
    ]

    for start_class, pool, key in classes:
        chosen = pick_rows(pool, count=starts_per_class, used_sequences=used_sequences, key=key)
        for replicate, row in enumerate(chosen):
            scan = local_scan_summary(landscape, row["_seq"], cfg)
            out = dict(row)
            out.update(scan)
            out["landscape_seed"] = int(landscape_seed)
            out["start_class"] = start_class
            out["start_replicate"] = int(replicate)
            out["start_id"] = f"seed{landscape_seed}__{start_class}__rep{replicate}"
            selected.append(out)

    max_reported = int(cfg["start_selection"].get("max_reported_candidates", 5000))
    reported = []
    for idx, row in enumerate(rows[:max_reported]):
        out = dict(row)
        out.pop("_seq", None)
        out["candidate_index"] = int(idx)
        out["landscape_seed"] = int(landscape_seed)
        reported.append(out)
    return selected, reported


def adaptive_mode(raw_mode: str) -> str:
    mode = str(raw_mode)
    if mode == "greedy_best_step":
        return "greedy"
    if mode == "stochastic_beneficial_step":
        return "random_beneficial"
    return mode


def run_adaptive_walk(
    landscape: FitnessLandscapeEnv,
    start_seq: np.ndarray,
    *,
    mode: str,
    rng: np.random.RandomState,
    max_steps: int,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    seq = np.asarray(start_seq, dtype=np.int64).copy()
    sampling_cfg = scan_config(cfg)
    trajectory: list[dict[str, Any]] = []
    previous_fitness: float | None = None
    stop_reason = "max_steps"

    for step in range(int(max_steps) + 1):
        scan = single_mutant_scan(landscape, seq, sampling_cfg=sampling_cfg)
        current_fitness = float(scan["parent_fitness"])
        beneficial_idx = np.flatnonzero(scan["beneficial"])
        row = {
            "step": int(step),
            "fitness": current_fitness,
            "step_gain": 0.0 if previous_fitness is None else current_fitness - previous_fitness,
            "available_beneficial_fraction": float(scan["fraction_beneficial"]),
            "available_neutral_fraction": float(scan["fraction_neutral"]),
            "available_deleterious_fraction": float(scan["fraction_deleterious"]),
            "available_mean_beneficial_effect": float(scan["mean_beneficial_effect"]),
            "best_available_effect": float(scan["best_delta"]),
            "chosen_effect": math.nan,
            "sequence": seq_to_string(seq),
            "_seq": seq.copy(),
            **components_for_seq(landscape, seq),
        }
        trajectory.append(row)
        if step == int(max_steps):
            break
        if beneficial_idx.size == 0:
            stop_reason = "no_beneficial_single_mutant"
            break
        if mode == "greedy":
            chosen = int(beneficial_idx[np.argmax(scan["deltas"][beneficial_idx])])
        elif mode == "random_beneficial":
            chosen = int(beneficial_idx[int(rng.randint(beneficial_idx.size))])
        else:
            raise ValueError(f"Unsupported adaptive-walk mode: {mode}")
        trajectory[-1]["chosen_effect"] = float(scan["deltas"][chosen])
        previous_fitness = current_fitness
        seq = scan["mutants"][chosen].copy()

    fitnesses = np.asarray([row["fitness"] for row in trajectory], dtype=np.float64)
    return {
        "trajectory": trajectory,
        "final_sequence": seq,
        "summary": {
            "mode": mode,
            "n_steps": int(max(0, len(trajectory) - 1)),
            "start_fitness": float(fitnesses[0]),
            "final_fitness": float(fitnesses[-1]),
            "total_gain": float(fitnesses[-1] - fitnesses[0]),
            "stop_reason": stop_reason,
        },
    }


def graphlet_step_indices(
    trajectory: list[dict[str, Any]],
    cfg: dict[str, Any] | None = None,
) -> list[tuple[int, str, bool]]:
    n_steps = len(trajectory) - 1
    selected: dict[int, tuple[str, bool]] = {}
    step_cfg = (
        cfg.get("adaptive_walks", {}).get("graphlet_step_classes", {})
        if cfg is not None
        else {}
    )
    early_indices = step_cfg.get("early", {}).get("step_indices", [0, 1, 2])
    middle_fracs = step_cfg.get("middle", {}).get(
        "fractional_walk_positions",
        [0.33, 0.66],
    )
    include_final = bool(step_cfg.get("late", {}).get("include_final", True))

    for idx in early_indices:
        idx = int(idx)
        if idx <= n_steps:
            selected[idx] = ("early", idx == n_steps)
    for frac in middle_fracs:
        idx = int(round(frac * n_steps))
        if 0 <= idx <= n_steps:
            selected.setdefault(idx, ("middle", idx == n_steps))
    if include_final:
        selected[n_steps] = ("late", True)
    return [(idx, cls, final) for idx, (cls, final) in sorted(selected.items())]


def build_position_library(
    profile: AlignmentProfile,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    max_dim = max(int(dim) for dim in cfg["graphlets"]["dimensions"])
    try:
        adjacency, edge_scores, vertices = dca_graph(
            profile,
            scan_top_couplings=int(cfg.get("dca_scan_top_couplings", 8000)),
            max_positions=int(min(profile.L, cfg.get("dca_max_positions", 160))),
        )
    except Exception:
        adjacency, edge_scores, vertices = {}, {}, []
    weighted_degree: defaultdict[int, float] = defaultdict(float)
    for (i, j), score in edge_scores.items():
        weighted_degree[int(i)] += float(score)
        weighted_degree[int(j)] += float(score)
    top_vertices = [
        int(pos)
        for pos, _ in sorted(weighted_degree.items(), key=lambda item: (-item[1], item[0]))
    ]
    if not top_vertices:
        top_vertices = list(range(profile.L))
    edges = sorted(edge_scores.items(), key=lambda item: (-float(item[1]), item[0]))
    cliques_by_dim: dict[int, list[dict[str, Any]]] = {}
    for dim in range(2, max_dim + 1):
        try:
            cliques_by_dim[dim] = enumerate_cliques(
                order=dim,
                adjacency=adjacency,
                edge_scores=edge_scores,
                vertices=vertices,
                max_candidates=500,
            )
        except Exception:
            cliques_by_dim[dim] = []
    return {
        "adjacency": adjacency,
        "edge_scores": edge_scores,
        "edges": edges,
        "top_vertices": top_vertices,
        "cliques_by_dim": cliques_by_dim,
        "L": int(profile.L),
    }


def contact_positions(
    library: dict[str, Any],
    *,
    dim: int,
    rank: int,
    rng: np.random.RandomState,
) -> list[int]:
    edges = library["edges"]
    positions: list[int] = []
    if edges:
        (i, j), _score = edges[int(rank) % len(edges)]
        positions = [int(i), int(j)]
        adjacency = library["adjacency"]
        edge_scores = library["edge_scores"]
        candidates: list[tuple[float, int]] = []
        for pos in positions:
            for neighbor in adjacency.get(pos, set()):
                if int(neighbor) in positions:
                    continue
                a, b = sorted((int(pos), int(neighbor)))
                candidates.append((float(edge_scores.get((a, b), 0.0)), int(neighbor)))
        for _score, pos in sorted(candidates, key=lambda item: (-item[0], item[1])):
            if pos not in positions:
                positions.append(pos)
            if len(positions) >= dim:
                break
    for pos in library["top_vertices"]:
        if len(positions) >= dim:
            break
        if int(pos) not in positions:
            positions.append(int(pos))
    while len(positions) < dim:
        pos = int(rng.randint(library["L"]))
        if pos not in positions:
            positions.append(pos)
    return positions[:dim]


def choose_positions(
    library: dict[str, Any],
    *,
    dim: int,
    contact_source: str,
    rank: int,
    rng: np.random.RandomState,
) -> tuple[list[int], dict[str, Any]]:
    if contact_source == "dca_clique_neighborhood":
        cliques = library["cliques_by_dim"].get(int(dim), [])
        if cliques:
            clique = cliques[int(rank) % len(cliques)]
            return list(clique["positions"]), {
                "contact_neighborhood_id": f"clique-{dim}-{int(rank) % len(cliques)}",
                "contact_edge_count": int(clique["edge_count"]),
                "contact_mean_edge_score": float(clique["mean_edge_score"]),
                "contact_min_edge_score": float(clique["min_edge_score"]),
            }
    if contact_source == "dca_contact_neighborhood":
        pos = contact_positions(library, dim=dim, rank=rank, rng=rng)
        return pos, {
            "contact_neighborhood_id": f"contact-{dim}-{rank}",
            "contact_edge_count": sum(
                1
                for i, j in itertools.combinations(pos, 2)
                if tuple(sorted((i, j))) in library["edge_scores"]
            ),
            "contact_mean_edge_score": mean_edge_score(pos, library["edge_scores"]),
            "contact_min_edge_score": min_edge_score(pos, library["edge_scores"]),
        }
    if contact_source == "shuffled_contact_neighborhood":
        pool = library["top_vertices"]
        if len(pool) >= dim:
            pos = [int(x) for x in rng.choice(np.asarray(pool), size=dim, replace=False)]
        else:
            pos = contact_positions(library, dim=dim, rank=rank, rng=rng)
            rng.shuffle(pos)
        return pos, {
            "contact_neighborhood_id": f"shuffled-contact-{dim}-{rank}",
            "contact_edge_count": 0,
            "contact_mean_edge_score": 0.0,
            "contact_min_edge_score": 0.0,
        }
    pos = [int(x) for x in rng.choice(np.arange(library["L"]), size=dim, replace=False)]
    return pos, {
        "contact_neighborhood_id": f"random-{dim}-{rank}",
        "contact_edge_count": 0,
        "contact_mean_edge_score": 0.0,
        "contact_min_edge_score": 0.0,
    }


def mean_edge_score(positions: list[int], edge_scores: dict[tuple[int, int], float]) -> float:
    scores = [
        float(edge_scores[tuple(sorted((i, j)))])
        for i, j in itertools.combinations(positions, 2)
        if tuple(sorted((i, j))) in edge_scores
    ]
    return float(np.mean(scores)) if scores else 0.0


def min_edge_score(positions: list[int], edge_scores: dict[tuple[int, int], float]) -> float:
    scores = [
        float(edge_scores[tuple(sorted((i, j)))])
        for i, j in itertools.combinations(positions, 2)
        if tuple(sorted((i, j))) in edge_scores
    ]
    return float(np.min(scores)) if scores else 0.0


def best_alternate_state(
    cfg: dict[str, Any],
    *,
    profile: AlignmentProfile,
    landscape: FitnessLandscapeEnv,
    background_seq: np.ndarray,
    parent_fitness: float,
    position: int,
    rng: np.random.RandomState,
) -> dict[str, Any]:
    state_cfg = cfg["graphlets"].get("state_selection", {})
    min_frequency = float(state_cfg.get("min_state_frequency", 0.0))
    require_observed = bool(state_cfg.get("require_observed_alignment_state", True))
    target_delta = float(state_cfg.get("target_single_delta", -20.0))
    tolerance = max(float(state_cfg.get("single_delta_tolerance", 120.0)), 1e-9)
    background_aa = int(background_seq[int(position)])
    candidates: list[dict[str, Any]] = []
    for aa in range(len(AA_ORDER)):
        if aa == background_aa:
            continue
        freq = marginal_frequency(profile, int(position), aa)
        if require_observed and freq < min_frequency:
            continue
        mutant = np.asarray(background_seq, dtype=np.int64).copy()
        mutant[int(position)] = aa
        delta = float(landscape.evaluate(mutant) - parent_fitness)
        score = math.log(max(freq, 1e-12)) - abs(delta - target_delta) / tolerance
        candidates.append(
            {
                "amino_acid": int(aa),
                "state_frequency": float(freq),
                "single_delta": delta,
                "state_score": float(score),
            }
        )
    if not candidates:
        aa_choices = np.delete(np.arange(len(AA_ORDER), dtype=np.int64), background_aa)
        aa = int(aa_choices[int(rng.randint(aa_choices.size))])
        mutant = np.asarray(background_seq, dtype=np.int64).copy()
        mutant[int(position)] = aa
        return {
            "amino_acid": aa,
            "state_frequency": marginal_frequency(profile, int(position), aa),
            "single_delta": float(landscape.evaluate(mutant) - parent_fitness),
            "state_score": 0.0,
        }
    candidates.sort(key=lambda row: (-float(row["state_score"]), -float(row["state_frequency"])))
    return candidates[0]


def sequence_from_bits(
    background_seq: np.ndarray,
    positions: list[int],
    amino_acids: list[int],
    mask: int,
) -> np.ndarray:
    seq = np.asarray(background_seq, dtype=np.int64).copy()
    for idx, (pos, aa) in enumerate(zip(positions, amino_acids)):
        if (int(mask) >> idx) & 1:
            seq[int(pos)] = int(aa)
    return seq


def node_sequences(
    background_seq: np.ndarray,
    positions: list[int],
    amino_acids: list[int],
) -> np.ndarray:
    dim = len(positions)
    return np.asarray(
        [sequence_from_bits(background_seq, positions, amino_acids, mask) for mask in range(1 << dim)],
        dtype=np.int64,
    )


def interleaved_balance_keys(keys: list[tuple[Any, ...]]) -> list[tuple[Any, ...]]:
    by_primary: defaultdict[Any, list[tuple[Any, ...]]] = defaultdict(list)
    for key in sorted(keys):
        primary = key[0] if key else ""
        by_primary[primary].append(key)
    for primary in by_primary:
        by_primary[primary].sort(key=lambda item: item[1:])
    ordered: list[tuple[Any, ...]] = []
    primaries = sorted(by_primary)
    while primaries:
        next_primaries: list[Any] = []
        for primary in primaries:
            bucket = by_primary[primary]
            if bucket:
                ordered.append(bucket.pop(0))
            if bucket:
                next_primaries.append(primary)
        primaries = next_primaries
    return ordered


def balanced_background_order(
    backgrounds: list[dict[str, Any]],
    cfg: dict[str, Any],
) -> list[dict[str, Any]]:
    """Order backgrounds so a graphlet cap cannot exhaust one start class first."""
    sampling_cfg = cfg.get("graphlets", {}).get("background_sampling", {})
    if str(sampling_cfg.get("mode", "balanced_round_robin")) == "input_order":
        return list(backgrounds)
    balance_keys = [
        str(key)
        for key in sampling_cfg.get("balance_keys", ["start_class", "walk_step_class"])
    ]
    groups: defaultdict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in backgrounds:
        key = tuple(row.get(field, "") for field in balance_keys)
        groups[key].append(row)
    for group in groups.values():
        group.sort(
            key=lambda row: (
                str(row.get("walk_step_class", "")),
                int(row.get("landscape_seed", 0)),
                str(row.get("walk_mode", "")),
                str(row.get("walk_id", "")),
                int(row.get("step", 0)),
                str(row.get("background_id", "")),
            )
        )

    ordered: list[dict[str, Any]] = []
    keys = interleaved_balance_keys(list(groups))
    while keys:
        next_keys: list[tuple[Any, ...]] = []
        for key in keys:
            group = groups[key]
            if group:
                ordered.append(group.pop(0))
            if group:
                next_keys.append(key)
        keys = next_keys
    return ordered


def balanced_graphlet_task_order(
    backgrounds: list[dict[str, Any]],
    cfg: dict[str, Any],
    *,
    dims: list[int],
    contact_sources: list[str],
    per_dim: int,
) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for bg in backgrounds:
        for dim in dims:
            for contact_source in contact_sources:
                for rep in range(per_dim):
                    tasks.append(
                        {
                            "background": bg,
                            "dim": int(dim),
                            "contact_source": str(contact_source),
                            "rep": int(rep),
                        }
                    )

    sampling_cfg = cfg.get("graphlets", {}).get("background_sampling", {})
    if str(sampling_cfg.get("mode", "balanced_round_robin")) == "input_order":
        return tasks
    balance_keys = [
        str(key)
        for key in sampling_cfg.get("balance_keys", ["start_class", "walk_step_class"])
    ]
    groups: defaultdict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for task in tasks:
        bg = task["background"]
        key = tuple(bg.get(field, "") for field in balance_keys)
        groups[key].append(task)
    for group in groups.values():
        group.sort(
            key=lambda task: (
                int(task["rep"]),
                int(task["dim"]),
                str(task["contact_source"]),
                int(task["background"].get("landscape_seed", 0)),
                str(task["background"].get("walk_mode", "")),
                str(task["background"].get("walk_id", "")),
                int(task["background"].get("step", 0)),
                str(task["background"].get("background_id", "")),
            )
        )

    ordered: list[dict[str, Any]] = []
    keys = interleaved_balance_keys(list(groups))
    while keys:
        next_keys: list[tuple[Any, ...]] = []
        for key in keys:
            group = groups[key]
            if group:
                ordered.append(group.pop(0))
            if group:
                next_keys.append(key)
        keys = next_keys
    return ordered


def annotate_global_background_scales(backgrounds: list[dict[str, Any]]) -> None:
    def finite_positive(key: str) -> np.ndarray:
        values = np.asarray([float(row.get(key, math.nan)) for row in backgrounds], dtype=np.float64)
        return values[np.isfinite(values) & (values > 1e-12)]

    median_local_sd = finite_positive("local_sd_delta")
    median_mean_beneficial = finite_positive("available_mean_beneficial_effect")
    median_best = finite_positive("best_available_effect")
    global_values = {
        "global_median_local_sd_delta": float(np.median(median_local_sd))
        if median_local_sd.size
        else math.nan,
        "global_median_available_mean_beneficial_effect": float(
            np.median(median_mean_beneficial)
        )
        if median_mean_beneficial.size
        else math.nan,
        "global_median_best_available_effect": float(np.median(median_best))
        if median_best.size
        else math.nan,
    }
    for row in backgrounds:
        row.update(global_values)


def sample_graphlets(
    cfg: dict[str, Any],
    *,
    profile: AlignmentProfile,
    landscape: FitnessLandscapeEnv | None = None,
    landscapes_by_seed: dict[int, FitnessLandscapeEnv] | None = None,
    position_library: dict[str, Any],
    backgrounds: list[dict[str, Any]],
    start_index: int = 0,
) -> list[dict[str, Any]]:
    graphlets: list[dict[str, Any]] = []
    dims = [int(dim) for dim in cfg["graphlets"]["dimensions"]]
    contact_sources = [str(source) for source in cfg["graphlets"]["contact_sources"]]
    per_dim = int(cfg["graphlets"].get("graphlets_per_background_per_dimension", 1))
    max_graphlets = int(cfg.get("execution_limits", {}).get("max_graphlets_total", 10**9))
    heldout_classes = {
        str(x)
        for x in cfg["graphlets"].get("calibration_test_split", {}).get("heldout_start_classes", [])
    }
    calibration_fraction = float(
        cfg["graphlets"].get("calibration_test_split", {}).get("calibration_fraction", 0.5)
    )
    tasks = balanced_graphlet_task_order(
        backgrounds,
        cfg,
        dims=dims,
        contact_sources=contact_sources,
        per_dim=per_dim,
    )

    for task in tasks:
        if start_index + len(graphlets) >= max_graphlets:
            break
        bg = task["background"]
        dim = int(task["dim"])
        contact_source = str(task["contact_source"])
        rep = int(task["rep"])
        bg_landscape = landscape
        if bg_landscape is None:
            if landscapes_by_seed is None:
                raise ValueError("sample_graphlets requires a landscape or landscapes_by_seed")
            bg_landscape = landscapes_by_seed[int(bg["landscape_seed"])]
        background_seq = np.asarray(bg["_seq"], dtype=np.int64)
        parent_fitness = float(bg["fitness"])
        graphlet_index = start_index + len(graphlets)
        rng = np.random.RandomState(
            stable_seed("graphlet", bg["background_id"], dim, contact_source, rep)
        )
        positions, contact_meta = choose_positions(
            position_library,
            dim=dim,
            contact_source=contact_source,
            rank=graphlet_index + rep,
            rng=rng,
        )
        state_rows = [
            best_alternate_state(
                cfg,
                profile=profile,
                landscape=bg_landscape,
                background_seq=background_seq,
                parent_fitness=parent_fitness,
                position=pos,
                rng=rng,
            )
            for pos in positions
        ]
        amino_acids = [int(row["amino_acid"]) for row in state_rows]
        mutations = [
            mutation_name(pos, int(background_seq[pos]), aa)
            for pos, aa in zip(positions, amino_acids)
        ]
        split_key = (
            bg["walk_id"],
            contact_meta["contact_neighborhood_id"],
            int(dim),
            int(rep),
        )
        split = (
            "test"
            if str(bg["start_class"]) in heldout_classes
            or stable_fraction("split", *split_key) >= calibration_fraction
            else "calibration"
        )
        graphlets.append(
            {
                "graphlet_id": f"graphlet-{graphlet_index:06d}",
                "split": split,
                "calibration_or_test": split,
                "graphlet_dimension": int(dim),
                "node_count": int(1 << dim),
                "contact_source": contact_source,
                "positions": ";".join(str(pos) for pos in positions),
                "amino_acids": ";".join(str(AA_ORDER[aa]) for aa in amino_acids),
                "mutations": ";".join(mutations),
                "state_frequencies": ";".join(
                    f"{float(row['state_frequency']):.8g}" for row in state_rows
                ),
                "single_deltas": ";".join(
                    f"{float(row['single_delta']):.8g}" for row in state_rows
                ),
                "_positions": positions,
                "_amino_acids": amino_acids,
                "_background_seq": background_seq.copy(),
                **{key: value for key, value in bg.items() if not key.startswith("_")},
                **contact_meta,
            }
        )
    return graphlets


def residual_parameter_sets(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    models = {str(row["id"]): row for row in cfg["candidate_models"]}
    limits = cfg.get("execution_limits", {})
    rows: list[dict[str, Any]] = []

    def add(model_id: str, params: dict[str, Any], role: str | None = None) -> None:
        if model_id not in models:
            return
        pset = dict(params)
        pset["candidate_model"] = model_id
        pset["candidate_role"] = role or str(models[model_id].get("role", ""))
        pset.setdefault("parameter_set_id", f"{model_id}__{len(rows):03d}")
        rows.append(pset)

    add("smooth_core_only", {"parameter_set_id": "smooth_core_only__base"})

    exact_grid = [float(x) for x in cfg["residual_parameter_sweep"].get("exact_match_strength_grid", [0.0])]
    for strength in choose_evenly(exact_grid, int(limits.get("max_exact_match_strengths", len(exact_grid)))):
        add(
            "res024_exact_match_haplotype",
            {
                "parameter_set_id": f"res024_exact_match_haplotype__strength-{strength:g}",
                "exact_match_strength": float(strength),
            },
        )

    pair_grid = [float(x) for x in cfg["residual_parameter_sweep"].get("dca_pair_strength_grid", [0.0])]
    for strength in choose_evenly(pair_grid, int(limits.get("max_pair_strengths", len(pair_grid)))):
        add(
            "dca_pair_override",
            {
                "parameter_set_id": f"dca_pair_override__strength-{strength:g}",
                "pair_strength": float(strength),
            },
        )

    sweep = cfg["residual_parameter_sweep"].get("contact_residual_ensemble", {})
    keys = [
        "residual_scale_relative_to_local_dfe_sd",
        "residual_scale_mode",
        "pairwise_density",
        "higher_order_density",
        "max_interaction_order",
        "coefficient_distribution",
        "sign_balance",
    ]
    combos = [
        {key: value for key, value in zip(keys, values)}
        for values in itertools.product(*(sweep.get(key, [None]) for key in keys))
    ]
    max_residual = int(limits.get("max_residual_parameter_sets_per_candidate", len(combos)))
    combos = choose_evenly(combos, max_residual)
    for model_id in (
        "contact_residual_ensemble",
        "shuffled_contact_residual",
        "random_position_residual",
        "gnk_local_null",
        "hard_certified_graphlet_positive_control",
        "residual_only_diagnostic",
    ):
        for idx, combo in enumerate(combos):
            params = dict(combo)
            params["parameter_set_id"] = f"{model_id}__pset-{idx:03d}"
            add(model_id, params)
    return rows


def walsh_residual_from_coefficients(coefficients: dict[int, float], *, dim: int) -> np.ndarray:
    n = 1 << dim
    values = np.zeros(n, dtype=np.float64)
    for mask, coeff in coefficients.items():
        for node in range(n):
            parity = (node & int(mask)).bit_count() % 2
            values[node] += float(coeff) * (-1.0 if parity else 1.0)
    values -= float(np.mean(values))
    return values


def random_residual_coefficients(
    parameter_set: dict[str, Any],
    *,
    dim: int,
    scale: float,
    rng: np.random.RandomState,
    force_high_order: bool,
) -> dict[int, float]:
    max_order = int(parameter_set.get("max_interaction_order") or 3)
    max_order = min(max_order, dim)
    pair_density = float(parameter_set.get("pairwise_density") or 0.1)
    higher_density = float(parameter_set.get("higher_order_density") or 0.0)
    distribution = str(parameter_set.get("coefficient_distribution") or "laplace")
    sign_balance = str(parameter_set.get("sign_balance") or "symmetric")
    masks_by_order: dict[int, list[int]] = defaultdict(list)
    for mask in range(1, 1 << dim):
        order = int(mask.bit_count())
        if order < 2 or order > max_order:
            continue
        density = pair_density if order == 2 else higher_density
        if rng.random_sample() <= density:
            masks_by_order[order].append(mask)
    if not masks_by_order.get(2):
        pair_masks = [mask for mask in range(1, 1 << dim) if mask.bit_count() == 2]
        if pair_masks:
            masks_by_order[2].append(int(pair_masks[int(rng.randint(len(pair_masks)))]))
    if force_high_order and max_order >= 3 and higher_density > 0.0:
        higher_masks = [
            mask
            for mask in range(1, 1 << dim)
            if 3 <= mask.bit_count() <= max_order
        ]
        if higher_masks and not any(order >= 3 for order in masks_by_order):
            masks_by_order[min(max_order, 3)].append(
                int(higher_masks[int(rng.randint(len(higher_masks)))])
            )
    selected = [mask for masks in masks_by_order.values() for mask in masks]
    if not selected:
        return {}
    normalizer = math.sqrt(max(len(selected), 1))
    coefficients: dict[int, float] = {}
    for mask in selected:
        order = int(mask.bit_count())
        if distribution == "student_t":
            raw = float(rng.standard_t(df=3))
        elif distribution == "laplace":
            raw = float(rng.laplace())
        else:
            raw = float(rng.normal())
        if sign_balance == "dca_enrichment_biased" and order == 2:
            raw = -abs(raw) if rng.random_sample() < 0.65 else abs(raw)
        coeff = raw * float(scale) / (normalizer * (order ** 1.25))
        coefficients[int(mask)] = float(coeff)
    return coefficients


def residual_scale_context(
    parameter_set: dict[str, Any],
    graphlet: dict[str, Any],
    *,
    base_values: np.ndarray,
) -> dict[str, float | str]:
    base_sd = float(np.std(base_values))
    local_dfe_sd = abs(float(graphlet.get("local_sd_delta", 0.0)))
    local_default = max(local_dfe_sd, base_sd, 1e-6)
    mode = str(parameter_set.get("residual_scale_mode") or "local_dfe_sd")

    def positive_graphlet_value(key: str) -> float | None:
        value = float(graphlet.get(key, math.nan))
        return value if np.isfinite(value) and value > 1e-12 else None

    selected: float | None
    if mode == "constant_walk_initial_beneficial_effect":
        selected = positive_graphlet_value("walk_initial_available_mean_beneficial_effect")
    elif mode == "constant_walk_initial_best_effect":
        selected = positive_graphlet_value("walk_initial_best_available_effect")
    elif mode == "constant_walk_initial_dfe_sd":
        selected = positive_graphlet_value("walk_initial_local_sd_delta")
    elif mode == "constant_global_median_beneficial_effect":
        selected = positive_graphlet_value("global_median_available_mean_beneficial_effect")
    elif mode == "constant_global_median_dfe_sd":
        selected = positive_graphlet_value("global_median_local_sd_delta")
    elif mode == "local_available_beneficial_effect":
        selected = positive_graphlet_value("available_mean_beneficial_effect")
    elif mode == "local_best_available_effect":
        selected = positive_graphlet_value("best_available_effect")
    else:
        mode = "local_dfe_sd"
        selected = local_default

    scale = max(float(selected if selected is not None else local_default), 1e-6)
    available_mean = positive_graphlet_value("available_mean_beneficial_effect")
    best_available = positive_graphlet_value("best_available_effect")
    return {
        "residual_scale_mode": mode,
        "residual_reference_scale": scale,
        "local_dfe_scale": local_default,
        "local_base_fitness_sd": max(base_sd, 1e-6),
        "local_available_mean_beneficial_effect_scale": float(available_mean or math.nan),
        "local_best_available_effect_scale": float(best_available or math.nan),
    }


def walsh_variance_summary(values: np.ndarray, *, dim: int) -> dict[str, Any]:
    rows, fractions = walsh_decomposition(values, dim=dim)
    absolute_by_order: defaultdict[int, float] = defaultdict(float)
    for row in rows:
        order = int(row["order"])
        if order <= 0:
            continue
        absolute_by_order[order] += float(row["coefficient"]) ** 2
    high_absolute = sum(
        float(value) for order, value in absolute_by_order.items() if int(order) >= 3
    )
    high_fraction = sum(
        float(value) for order, value in fractions.items() if int(order) >= 3
    )
    return {
        "rows": rows,
        "fraction_by_order": {int(order): float(value) for order, value in fractions.items()},
        "absolute_by_order": dict(absolute_by_order),
        "high_order_fraction": float(high_fraction),
        "high_order_absolute_variance": float(high_absolute),
        "order1_absolute_variance": float(absolute_by_order.get(1, 0.0)),
        "total_variance": float(np.var(values)),
    }


def residual_values(
    candidate_model: str,
    parameter_set: dict[str, Any],
    graphlet: dict[str, Any],
    *,
    base_values: np.ndarray,
    local_scale: float,
    rng: np.random.RandomState,
    cfg: dict[str, Any],
) -> np.ndarray:
    dim = int(graphlet["graphlet_dimension"])
    n = 1 << dim
    if candidate_model == "smooth_core_only":
        return np.zeros(n, dtype=np.float64)
    if candidate_model == "res024_exact_match_haplotype":
        strength = float(parameter_set.get("exact_match_strength", 0.0))
        values = np.zeros(n, dtype=np.float64)
        values[n - 1] = strength
        return values - float(np.mean(values))
    if candidate_model == "dca_pair_override":
        strength = float(parameter_set.get("pair_strength", 0.0)) * local_scale
        values = np.zeros(n, dtype=np.float64)
        for i, j in itertools.combinations(range(dim), 2):
            sign = -1.0 if (i + j) % 2 == 0 else 1.0
            for node in range(n):
                if ((node >> i) & 1) and ((node >> j) & 1):
                    values[node] += sign * strength
        return values - float(np.mean(values))
    if candidate_model == "gnk_local_null":
        raw = rng.normal(size=n)
        raw -= float(np.mean(raw))
        sd = float(np.std(raw))
        return raw * (local_scale / sd) if sd > 1e-12 else raw

    scale = (
        float(parameter_set.get("residual_scale_relative_to_local_dfe_sd") or 0.25)
        * local_scale
    )
    force_high = candidate_model in {"contact_residual_ensemble", "residual_only_diagnostic"}
    coeffs = random_residual_coefficients(
        parameter_set,
        dim=dim,
        scale=scale,
        rng=rng,
        force_high_order=force_high,
    )
    values = walsh_residual_from_coefficients(coeffs, dim=dim)
    if candidate_model == "hard_certified_graphlet_positive_control":
        attempts = int(cfg.get("execution_limits", {}).get("max_hard_certified_attempts", 32))
        best = values
        best_score = -math.inf
        for _attempt in range(attempts):
            coeffs = random_residual_coefficients(
                parameter_set,
                dim=dim,
                scale=scale,
                rng=rng,
                force_high_order=True,
            )
            candidate = walsh_residual_from_coefficients(coeffs, dim=dim)
            fitness = base_values + candidate
            signs = sign_epistasis_counts(fitness, dim=dim)
            extrema = local_extrema_counts(fitness, dim=dim)
            access = accessible_path_metrics(fitness, dim=dim)
            score = (
                float(signs["reciprocal_sign_epistasis_fraction"])
                + 0.05 * max(0, int(extrema["n_local_maxima"]) - 1)
                + (1.0 - float(access["accessible_path_fraction"]))
            )
            if score > best_score:
                best = candidate
                best_score = score
            if (
                float(signs["reciprocal_sign_epistasis_fraction"]) > 0.0
                and int(extrema["n_local_maxima"]) > 1
                and float(access["accessible_path_fraction"]) < 1.0
            ):
                return candidate
        return best
    return values


def autocorrelation(values: np.ndarray, lag: int) -> float | None:
    if lag <= 0 or values.size <= lag + 2:
        return None
    x = values[:-lag]
    y = values[lag:]
    if float(np.std(x)) <= 1e-12 or float(np.std(y)) <= 1e-12:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def graph_random_walk_metrics(
    values: np.ndarray,
    *,
    dim: int,
    cfg: dict[str, Any],
    rng: np.random.RandomState,
) -> dict[str, Any]:
    walk_cfg = cfg["analysis"]["random_walk"]
    n_walks = int(walk_cfg.get("n_walks_per_graphlet", 40))
    steps = int(walk_cfg.get("steps", 64))
    lag_max = int(walk_cfg.get("lag_max", 8))
    lag_values: defaultdict[int, list[float]] = defaultdict(list)
    deltas: list[float] = []
    for _ in range(n_walks):
        node = int(rng.randint(1 << dim))
        path = [float(values[node])]
        for _step in range(steps):
            bit = int(rng.randint(dim))
            node ^= 1 << bit
            path.append(float(values[node]))
        arr = np.asarray(path, dtype=np.float64)
        deltas.extend(np.diff(arr).tolist())
        for lag in range(1, lag_max + 1):
            corr = autocorrelation(arr, lag)
            if corr is not None and np.isfinite(corr):
                lag_values[lag].append(float(corr))
    lag1 = float(np.mean(lag_values[1])) if lag_values[1] else math.nan
    lag5 = float(np.mean(lag_values[5])) if lag_values[5] else math.nan
    corr_len = (
        float(-1.0 / math.log(lag1))
        if np.isfinite(lag1) and 0.0 < lag1 < 1.0
        else math.nan
    )
    out = {
        "random_walk_lag1_autocorrelation": lag1,
        "random_walk_lag5_autocorrelation": lag5,
        "random_walk_autocorrelation_length": corr_len,
        "random_walk_step_sd": float(np.std(deltas)) if deltas else math.nan,
    }
    for lag in range(1, lag_max + 1):
        out[f"random_walk_lag{lag}_autocorrelation"] = (
            float(np.mean(lag_values[lag])) if lag_values[lag] else math.nan
        )
    return out


def accessible_path_metrics(values: np.ndarray, *, dim: int) -> dict[str, Any]:
    start = 0
    end = (1 << dim) - 1
    if float(values[end]) < float(values[start]):
        start, end = end, start
    diff_bits = [idx for idx in range(dim) if ((start ^ end) >> idx) & 1]
    total = math.factorial(len(diff_bits))
    accessible = 0
    for order in itertools.permutations(diff_bits):
        node = start
        previous = float(values[node])
        ok = True
        for bit in order:
            node ^= 1 << bit
            current = float(values[node])
            if current <= previous + 1e-12:
                ok = False
                break
            previous = current
        accessible += int(ok)
    return {
        "accessible_path_start_node": int(start),
        "accessible_path_end_node": int(end),
        "accessible_path_total_shortest_paths": int(total),
        "accessible_path_count": int(accessible),
        "accessible_path_fraction": float(accessible / total) if total else math.nan,
        "accessible_path_endpoint_delta": float(values[end] - values[start]),
    }


def basin_metrics(values: np.ndarray, *, dim: int) -> dict[str, Any]:
    attractors: list[int] = []
    for start in range(1 << dim):
        node = start
        seen: set[int] = set()
        while node not in seen:
            seen.add(node)
            neighbors = [node ^ (1 << bit) for bit in range(dim)]
            best_neighbor = max(neighbors, key=lambda idx: float(values[idx]))
            if float(values[best_neighbor]) <= float(values[node]) + 1e-12:
                break
            node = int(best_neighbor)
        attractors.append(int(node))
    counts = Counter(attractors)
    probs = np.asarray(list(counts.values()), dtype=np.float64) / float(1 << dim)
    entropy = -float(np.sum(probs * np.log2(np.maximum(probs, 1e-12)))) if probs.size else 0.0
    return {
        "basin_count": int(len(counts)),
        "basin_entropy": entropy,
        "basin_largest_fraction": float(np.max(probs)) if probs.size else math.nan,
        "basin_attractors": ";".join(str(node) for node, _count in counts.most_common()),
    }


def gft_rows_from_walsh(
    walsh_rows: list[dict[str, Any]],
    *,
    graphlet: dict[str, Any],
    parameter_set: dict[str, Any],
    total_variance: float,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    dim = int(graphlet["graphlet_dimension"])
    rows: list[dict[str, Any]] = []
    energy_by_band: defaultdict[str, float] = defaultdict(float)
    for row in walsh_rows:
        order = int(row["order"])
        eigenvalue = 2 * order
        if order == 0:
            band = "dc"
        else:
            q = order / max(dim, 1)
            if q <= 0.33:
                band = "low"
            elif q <= 0.67:
                band = "mid"
            else:
                band = "high"
        energy = float(row["coefficient"]) ** 2
        if order > 0:
            energy_by_band[band] += energy
        rows.append(
            {
                "graphlet_id": graphlet["graphlet_id"],
                "split": graphlet["split"],
                "candidate_model": parameter_set["candidate_model"],
                "parameter_set_id": parameter_set["parameter_set_id"],
                "graphlet_dimension": int(dim),
                "gft_basis": "binary_hamming_walsh_equivalent",
                "mask": int(row["mask"]),
                "order": order,
                "laplacian_eigenvalue": int(eigenvalue),
                "eigenvalue_band": band,
                "coefficient": float(row["coefficient"]),
                "energy": energy,
                "energy_fraction": float(energy / total_variance) if total_variance > 1e-12 and order > 0 else 0.0,
            }
        )
    fractions = {
        f"gft_{band}_energy_fraction": (
            float(energy_by_band[band] / total_variance) if total_variance > 1e-12 else 0.0
        )
        for band in ("low", "mid", "high")
    }
    return rows, fractions


def selected_term_rows(
    rows: list[dict[str, Any]],
    cfg: dict[str, Any],
    *,
    section: str,
) -> list[dict[str, Any]]:
    section_cfg = cfg.get("analysis", {}).get(section, {})
    record_key = "record_all_terms" if section == "walsh" else "record_all_coefficients"
    if bool(section_cfg.get(record_key, True)):
        return rows
    max_rows = int(section_cfg.get("max_recorded_terms_per_graphlet_model", 0))
    if max_rows <= 0:
        return []
    return sorted(
        rows,
        key=lambda row: (
            int(row.get("order", 0)) == 0,
            -abs(float(row.get("coefficient", 0.0))),
            int(row.get("order", 0)),
            int(row.get("mask", 0)),
        ),
    )[:max_rows]


def evaluate_graphlet_model(
    cfg: dict[str, Any],
    *,
    landscape: FitnessLandscapeEnv,
    graphlet: dict[str, Any],
    parameter_set: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    positions = list(graphlet["_positions"])
    amino_acids = list(graphlet["_amino_acids"])
    background_seq = np.asarray(graphlet["_background_seq"], dtype=np.int64)
    dim = int(graphlet["graphlet_dimension"])
    seqs = node_sequences(background_seq, positions, amino_acids)
    base_values = evaluate_batch(landscape, seqs, batch_size=int(cfg.get("batch_size", 4096)))
    scale_context = residual_scale_context(parameter_set, graphlet, base_values=base_values)
    local_scale = float(scale_context["residual_reference_scale"])
    rng = np.random.RandomState(
        stable_seed(
            "residual",
            graphlet["graphlet_id"],
            parameter_set["candidate_model"],
            parameter_set["parameter_set_id"],
            graphlet["landscape_seed"],
        )
    )
    residual = residual_values(
        str(parameter_set["candidate_model"]),
        parameter_set,
        graphlet,
        base_values=base_values,
        local_scale=local_scale,
        rng=rng,
        cfg=cfg,
    )
    if parameter_set["candidate_model"] == "residual_only_diagnostic":
        fitness = residual.copy()
    else:
        fitness = base_values + residual
    total_variance = float(np.var(fitness))
    residual_variance = float(np.var(residual))
    base_variance = float(np.var(base_values))
    walsh_rows, variance_by_order = walsh_decomposition(fitness, dim=dim)
    base_walsh = walsh_variance_summary(base_values, dim=dim)
    residual_walsh = walsh_variance_summary(residual, dim=dim)
    sign_counts = sign_epistasis_counts(fitness, dim=dim)
    extrema = local_extrema_counts(fitness, dim=dim)
    access = accessible_path_metrics(fitness, dim=dim)
    basins = basin_metrics(fitness, dim=dim)
    random_walk = graph_random_walk_metrics(
        fitness,
        dim=dim,
        cfg=cfg,
        rng=np.random.RandomState(stable_seed("rw", graphlet["graphlet_id"], parameter_set["parameter_set_id"])),
    )
    gft_rows, gft_summary = gft_rows_from_walsh(
        walsh_rows,
        graphlet=graphlet,
        parameter_set=parameter_set,
        total_variance=total_variance,
    )
    high_order_fraction = sum(
        float(value) for order, value in variance_by_order.items() if int(order) >= 3
    )
    total_high_order_absolute = float(high_order_fraction * total_variance)
    base_order1_absolute = float(base_walsh["order1_absolute_variance"])
    residual_high_absolute = float(residual_walsh["high_order_absolute_variance"])
    available_mean_scale = float(scale_context["local_available_mean_beneficial_effect_scale"])
    best_available_scale = float(scale_context["local_best_available_effect_scale"])
    metric_row = {
        **{key: value for key, value in graphlet.items() if not key.startswith("_")},
        **parameter_set,
        "fitness_min": float(np.min(fitness)),
        "fitness_mean": float(np.mean(fitness)),
        "fitness_max": float(np.max(fitness)),
        "fitness_sd": float(np.std(fitness)),
        "base_fitness_sd": float(np.std(base_values)),
        "residual_sd": float(np.std(residual)),
        "base_variance": base_variance,
        "residual_variance": residual_variance,
        "total_variance": total_variance,
        "residual_to_base_sd_ratio": float(np.std(residual) / max(np.std(base_values), 1e-12)),
        **scale_context,
        "residual_scale_multiplier": float(
            parameter_set.get("residual_scale_relative_to_local_dfe_sd") or 0.0
        ),
        "residual_to_available_mean_beneficial_effect_ratio": (
            float(np.std(residual) / available_mean_scale)
            if np.isfinite(available_mean_scale) and available_mean_scale > 1e-12
            else math.nan
        ),
        "residual_to_best_available_effect_ratio": (
            float(np.std(residual) / best_available_scale)
            if np.isfinite(best_available_scale) and best_available_scale > 1e-12
            else math.nan
        ),
        "total_high_order_absolute_variance": total_high_order_absolute,
        "base_walsh_high_order_fraction": float(base_walsh["high_order_fraction"]),
        "base_high_order_absolute_variance": float(base_walsh["high_order_absolute_variance"]),
        "base_order1_absolute_variance": base_order1_absolute,
        "residual_walsh_high_order_fraction": float(residual_walsh["high_order_fraction"]),
        "residual_high_order_absolute_variance": residual_high_absolute,
        "residual_order1_absolute_variance": float(residual_walsh["order1_absolute_variance"]),
        "residual_high_order_to_base_order1_variance_ratio": (
            residual_high_absolute / max(base_order1_absolute, 1e-12)
        ),
        "total_high_order_to_base_order1_variance_ratio": (
            total_high_order_absolute / max(base_order1_absolute, 1e-12)
        ),
        "walsh_high_order_fraction": float(high_order_fraction),
        "walsh_high_order_present": bool(high_order_fraction > 1e-6),
        **{f"walsh_variance_order_{order}": float(value) for order, value in variance_by_order.items()},
        **{
            f"base_walsh_variance_order_{order}": float(value)
            for order, value in base_walsh["fraction_by_order"].items()
        },
        **{
            f"base_walsh_absolute_variance_order_{order}": float(value)
            for order, value in base_walsh["absolute_by_order"].items()
        },
        **{
            f"residual_walsh_variance_order_{order}": float(value)
            for order, value in residual_walsh["fraction_by_order"].items()
        },
        **{
            f"residual_walsh_absolute_variance_order_{order}": float(value)
            for order, value in residual_walsh["absolute_by_order"].items()
        },
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
                "split": graphlet["split"],
                "candidate_model": parameter_set["candidate_model"],
                "parameter_set_id": parameter_set["parameter_set_id"],
                "graphlet_dimension": int(dim),
            }
        )
    sign_pair_rows = []
    for row in sign_counts["pair_rows"]:
        out = dict(row)
        out.update(
            {
                "graphlet_id": graphlet["graphlet_id"],
                "split": graphlet["split"],
                "candidate_model": parameter_set["candidate_model"],
                "parameter_set_id": parameter_set["parameter_set_id"],
                "mutation_i": str(graphlet["mutations"]).split(";")[int(row["var_i"])],
                "mutation_j": str(graphlet["mutations"]).split(";")[int(row["var_j"])],
            }
        )
        sign_pair_rows.append(out)
    return (
        metric_row,
        selected_term_rows(walsh_rows, cfg, section="walsh"),
        selected_term_rows(gft_rows, cfg, section="gft"),
        sign_pair_rows,
        access,
    )


def aggregate_metric_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    metrics = [
        "sign_epistasis_fraction",
        "reciprocal_sign_epistasis_fraction",
        "n_local_maxima",
        "n_local_minima",
        "basin_count",
        "basin_entropy",
        "accessible_path_fraction",
        "random_walk_lag1_autocorrelation",
        "random_walk_lag5_autocorrelation",
        "random_walk_autocorrelation_length",
        "walsh_high_order_fraction",
        "gft_low_energy_fraction",
        "gft_mid_energy_fraction",
        "gft_high_energy_fraction",
        "residual_to_base_sd_ratio",
        "residual_to_available_mean_beneficial_effect_ratio",
        "residual_to_best_available_effect_ratio",
        "base_walsh_high_order_fraction",
        "residual_walsh_high_order_fraction",
        "base_order1_absolute_variance",
        "residual_high_order_absolute_variance",
        "residual_high_order_to_base_order1_variance_ratio",
        "total_high_order_to_base_order1_variance_ratio",
    ]
    groups: defaultdict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = (str(row["split"]), str(row["candidate_model"]), str(row["parameter_set_id"]))
        groups[key].append(row)
    out: list[dict[str, Any]] = []
    for (split, model, pset), group_rows in sorted(groups.items()):
        row: dict[str, Any] = {
            "split": split,
            "candidate_model": model,
            "parameter_set_id": pset,
            "n_graphlets": int(len(group_rows)),
        }
        for metric in metrics:
            values = np.asarray(
                [float(item.get(metric, math.nan)) for item in group_rows],
                dtype=np.float64,
            )
            values = values[np.isfinite(values)]
            row[f"{metric}_mean"] = float(np.mean(values)) if values.size else math.nan
            row[f"{metric}_median"] = float(np.median(values)) if values.size else math.nan
            row[f"{metric}_q10"] = float(np.quantile(values, 0.10)) if values.size else math.nan
            row[f"{metric}_q90"] = float(np.quantile(values, 0.90)) if values.size else math.nan
        out.append(row)
    return out


def calibration_score(row: dict[str, Any], baseline: dict[str, float]) -> float:
    rse = float(row.get("reciprocal_sign_epistasis_fraction_mean", 0.0))
    se = float(row.get("sign_epistasis_fraction_mean", 0.0))
    maxima = float(row.get("n_local_maxima_mean", 1.0))
    access = float(row.get("accessible_path_fraction_mean", 1.0))
    lag1 = float(row.get("random_walk_lag1_autocorrelation_mean", 1.0))
    high = float(row.get("walsh_high_order_fraction_mean", 0.0))
    gft_high = float(row.get("gft_high_energy_fraction_mean", 0.0))
    high_penalty = max(0.0, high - 0.35) + max(0.0, gft_high - 0.35)
    return (
        (se - baseline.get("se", 0.0))
        + 2.0 * (rse - baseline.get("rse", 0.0))
        + 0.1 * max(0.0, maxima - baseline.get("maxima", 1.0))
        + max(0.0, baseline.get("access", 1.0) - access)
        + 0.5 * max(0.0, baseline.get("lag1", 1.0) - lag1)
        + 0.5 * high
        - 2.0 * high_penalty
    )


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    aggregate = aggregate_metric_rows(rows)
    smooth = [
        row
        for row in aggregate
        if row["split"] == "calibration" and row["candidate_model"] == "smooth_core_only"
    ]
    baseline = {
        "se": float(smooth[0].get("sign_epistasis_fraction_mean", 0.0)) if smooth else 0.0,
        "rse": float(smooth[0].get("reciprocal_sign_epistasis_fraction_mean", 0.0)) if smooth else 0.0,
        "maxima": float(smooth[0].get("n_local_maxima_mean", 1.0)) if smooth else 1.0,
        "access": float(smooth[0].get("accessible_path_fraction_mean", 1.0)) if smooth else 1.0,
        "lag1": float(smooth[0].get("random_walk_lag1_autocorrelation_mean", 1.0)) if smooth else 1.0,
    }
    for row in aggregate:
        row["calibration_score"] = (
            calibration_score(row, baseline) if row["split"] == "calibration" else math.nan
        )
    selected: dict[str, dict[str, Any]] = {}
    for row in aggregate:
        if row["split"] != "calibration":
            continue
        model = str(row["candidate_model"])
        if model not in selected or float(row["calibration_score"]) > float(selected[model]["calibration_score"]):
            selected[model] = row
    return {
        "aggregate_rows": aggregate,
        "calibration_baseline": baseline,
        "selected_parameter_sets_by_calibration": selected,
    }


def run_experiment(
    cfg: dict[str, Any],
    *,
    config_path: Path,
    project_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    profile = build_alignment_profile(
        cfg,
        project_root=project_root,
        config_path=config_path,
        output_dir=output_dir,
    )
    position_library = build_position_library(profile, cfg)

    start_rows: list[dict[str, Any]] = []
    start_candidate_report: list[dict[str, Any]] = []
    walk_rows: list[dict[str, Any]] = []
    backgrounds: list[dict[str, Any]] = []
    landscape_by_seed: dict[int, FitnessLandscapeEnv] = {}

    for landscape_seed in [int(seed) for seed in cfg["seeds"]["landscape"]]:
        landscape = build_base_landscape(cfg, seed=landscape_seed, profile=profile)
        landscape_by_seed[landscape_seed] = landscape
        starts, candidates = select_starts(
            cfg,
            profile=profile,
            landscape=landscape,
            landscape_seed=landscape_seed,
        )
        start_rows.extend(starts)
        start_candidate_report.extend(candidates)
        for start in starts:
            for raw_mode in cfg["adaptive_walks"]["modes"]:
                mode = adaptive_mode(str(raw_mode))
                for replicate in range(int(cfg["adaptive_walks"]["replicates_per_start"])):
                    walk_id = f"{start['start_id']}__{mode}__walkrep{replicate}"
                    rng = np.random.RandomState(stable_seed("walk", walk_id))
                    walk = run_adaptive_walk(
                        landscape,
                        start["_seq"],
                        mode=mode,
                        rng=rng,
                        max_steps=int(cfg["adaptive_walks"]["max_steps"]),
                        cfg=cfg,
                    )
                    initial_row = walk["trajectory"][0]
                    initial_scan = local_scan_summary(landscape, initial_row["_seq"], cfg)
                    terminal_step = int(walk["trajectory"][-1]["step"])
                    for row in walk["trajectory"]:
                        out = {
                            **{key: value for key, value in row.items() if not key.startswith("_")},
                            "landscape_seed": int(landscape_seed),
                            "walk_id": walk_id,
                            "walk_mode": mode,
                            "start_id": start["start_id"],
                            "start_class": start["start_class"],
                            "start_replicate": int(start["start_replicate"]),
                        }
                        walk_rows.append(out)
                    for step_index, step_class, is_final in graphlet_step_indices(walk["trajectory"], cfg):
                        row = walk["trajectory"][step_index]
                        scan = local_scan_summary(landscape, row["_seq"], cfg)
                        bg = {
                            **{key: value for key, value in row.items() if not key.startswith("_")},
                            **scan,
                            "background_id": f"{walk_id}__step{step_index}",
                            "landscape_seed": int(landscape_seed),
                            "walk_id": walk_id,
                            "walk_mode": mode,
                            "walk_step_class": step_class,
                            "walk_step_is_final": bool(is_final),
                            "walk_terminal_step": terminal_step,
                            "distance_to_walk_terminal_step": int(terminal_step - step_index),
                            "walk_initial_fitness": float(initial_row["fitness"]),
                            "walk_initial_available_beneficial_fraction": float(
                                initial_row["available_beneficial_fraction"]
                            ),
                            "walk_initial_available_mean_beneficial_effect": float(
                                initial_row["available_mean_beneficial_effect"]
                            ),
                            "walk_initial_best_available_effect": float(
                                initial_row["best_available_effect"]
                            ),
                            "walk_initial_local_sd_delta": float(initial_scan["local_sd_delta"]),
                            "start_id": start["start_id"],
                            "start_class": start["start_class"],
                            "start_replicate": int(start["start_replicate"]),
                            "_seq": row["_seq"].copy(),
                        }
                        backgrounds.append(bg)

    rank_quantiles(backgrounds, "fitness", "background_fitness_quantile")
    rank_quantiles(backgrounds, "local_neutral_fraction", "local_neutral_fraction_quantile")
    rank_quantiles(backgrounds, "local_sd_delta", "local_dfe_variance_quantile")
    annotate_global_background_scales(backgrounds)

    ordered_backgrounds = balanced_background_order(backgrounds, cfg)
    graphlets = sample_graphlets(
        cfg,
        profile=profile,
        landscapes_by_seed=landscape_by_seed,
        position_library=position_library,
        backgrounds=ordered_backgrounds,
    )

    parameter_sets = residual_parameter_sets(cfg)
    metric_rows: list[dict[str, Any]] = []
    walsh_rows: list[dict[str, Any]] = []
    gft_rows: list[dict[str, Any]] = []
    sign_pair_rows: list[dict[str, Any]] = []
    accessibility_rows: list[dict[str, Any]] = []

    for graphlet in graphlets:
        landscape = landscape_by_seed[int(graphlet["landscape_seed"])]
        for parameter_set in parameter_sets:
            metric, walsh, gft, sign_pairs, access = evaluate_graphlet_model(
                cfg,
                landscape=landscape,
                graphlet=graphlet,
                parameter_set=parameter_set,
            )
            metric_rows.append(metric)
            walsh_rows.extend(walsh)
            gft_rows.extend(gft)
            if bool(
                cfg.get("analysis", {})
                .get("sign_epistasis", {})
                .get("record_pair_rows", True)
            ):
                sign_pair_rows.extend(sign_pairs)
            accessibility_rows.append(
                {
                    "graphlet_id": graphlet["graphlet_id"],
                    "split": graphlet["split"],
                    "candidate_model": parameter_set["candidate_model"],
                    "parameter_set_id": parameter_set["parameter_set_id"],
                    **access,
                }
            )

    summary = summarize(metric_rows)
    aggregate_rows = summary["aggregate_rows"]
    calibration_rows = [row for row in metric_rows if row["split"] == "calibration"]
    test_rows = [row for row in metric_rows if row["split"] == "test"]

    graphlet_catalog = [
        {key: value for key, value in row.items() if not key.startswith("_")}
        for row in graphlets
    ]
    starts_out = [
        {key: value for key, value in row.items() if not key.startswith("_")}
        for row in start_rows
    ]
    backgrounds_out = [
        {key: value for key, value in row.items() if not key.startswith("_")}
        for row in backgrounds
    ]

    write_csv(output_dir / "adaptive_walk_starts.csv", starts_out)
    write_csv(output_dir / "start_candidate_report.csv", start_candidate_report)
    write_csv(output_dir / "adaptive_walk_trajectories.csv", walk_rows)
    write_csv(output_dir / "adaptive_walk_graphlet_backgrounds.csv", backgrounds_out)
    write_csv(output_dir / "graphlet_catalog.csv", graphlet_catalog)
    write_csv(output_dir / "module_parameter_sweep.csv", aggregate_rows)
    write_csv(output_dir / "calibration_graphlet_metrics.csv", calibration_rows)
    write_csv(output_dir / "test_graphlet_metrics.csv", test_rows)
    write_csv(output_dir / "graphlet_walsh_terms.csv", walsh_rows)
    write_csv(output_dir / "graphlet_gft_coefficients.csv", gft_rows)
    write_csv(output_dir / "graphlet_sign_epistasis_pairs.csv", sign_pair_rows)
    write_csv(output_dir / "graphlet_accessibility.csv", accessibility_rows)
    write_json(
        output_dir / "calibration_target_summary.json",
        {
            "baseline": summary["calibration_baseline"],
            "selected_parameter_sets_by_calibration": summary[
                "selected_parameter_sets_by_calibration"
            ],
            "aggregate_rows": [
                row for row in aggregate_rows if row["split"] == "calibration"
            ],
        },
    )
    write_json(
        output_dir / "heldout_test_summary.json",
        {
            "aggregate_rows": [row for row in aggregate_rows if row["split"] == "test"],
            "selected_parameter_sets_by_calibration": summary[
                "selected_parameter_sets_by_calibration"
            ],
        },
    )
    run_summary = {
        "experiment_id": cfg.get("experiment_id", "EXP-025"),
        "n_alignment_sequences": int(profile.N),
        "n_landscape_seeds": int(len(landscape_by_seed)),
        "n_starts": int(len(start_rows)),
        "n_walk_rows": int(len(walk_rows)),
        "n_graphlet_backgrounds": int(len(backgrounds)),
        "n_graphlets": int(len(graphlets)),
        "n_parameter_sets": int(len(parameter_sets)),
        "n_metric_rows": int(len(metric_rows)),
        "n_walsh_rows": int(len(walsh_rows)),
        "n_gft_rows": int(len(gft_rows)),
        "n_sign_pair_rows": int(len(sign_pair_rows)),
        "calibration_graphlets": int(sum(1 for row in graphlets if row["split"] == "calibration")),
        "test_graphlets": int(sum(1 for row in graphlets if row["split"] == "test")),
        "background_start_class_counts": dict(
            sorted(Counter(str(row["start_class"]) for row in backgrounds).items())
        ),
        "graphlet_start_class_counts": dict(
            sorted(Counter(str(row["start_class"]) for row in graphlets).items())
        ),
        "graphlet_walk_step_class_counts": dict(
            sorted(Counter(str(row["walk_step_class"]) for row in graphlets).items())
        ),
        "record_all_walsh_terms": bool(
            cfg.get("analysis", {}).get("walsh", {}).get("record_all_terms", True)
        ),
        "record_all_gft_coefficients": bool(
            cfg.get("analysis", {}).get("gft", {}).get("record_all_coefficients", True)
        ),
        "adaptive_walk_source_model": "smooth_core_only",
        "gft_basis_note": (
            "Dense binary Hamming graphlets use the Walsh basis as the graph "
            "Fourier basis; Laplacian eigenvalue is 2 * Walsh order."
        ),
        "selected_parameter_sets_by_calibration": summary[
            "selected_parameter_sets_by_calibration"
        ],
    }
    write_json(output_dir / "summary.json", run_summary)
    return run_summary


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    cfg = load_config(config_path, quick=bool(args.quick))
    summary = run_experiment(
        cfg,
        config_path=config_path,
        project_root=project_root,
        output_dir=output_dir,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))


if __name__ == "__main__":
    main()
