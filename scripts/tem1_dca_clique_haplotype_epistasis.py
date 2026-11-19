#!/usr/bin/env python3
"""Run TEM-1 DCA-clique haplotype high-order epistasis diagnostics."""

from __future__ import annotations

import argparse
import csv
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
    json_default,
    merge_dict,
    write_json,
)
from tem1_homologue_walk_parameter_grid import (  # noqa: E402
    build_alignment_profile,
    components_for_seq,
    filled_unique_homologues,
    seq_to_string,
)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = []
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    fieldnames.append(key)
        if not fieldnames:
            fieldnames = ["empty"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


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


def aa_name(index: int) -> str:
    return str(AA_ORDER[int(index)])


def mutation_name(position: int, from_aa: int, to_aa: int) -> str:
    return f"{aa_name(from_aa)}{int(position) + 1}{aa_name(to_aa)}"


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


def evaluate_fitness(landscape: FitnessLandscapeEnv, seq: np.ndarray) -> float:
    return float(landscape.evaluate(np.asarray(seq, dtype=np.int64)))


def select_anchor_sequence(
    cfg: dict[str, Any],
    *,
    profile: AlignmentProfile,
    landscape: FitnessLandscapeEnv,
) -> tuple[np.ndarray, dict[str, Any], list[dict[str, Any]]]:
    anchor_cfg = dict(cfg.get("anchor_selection", {}))
    mode = str(anchor_cfg.get("mode", "closest_realized_margin")).lower()
    reference = np.asarray(landscape.reference, dtype=np.int64)

    candidates: list[tuple[str, int | None, np.ndarray]] = [("reference", None, reference)]
    if mode != "reference":
        homologue_seqs, source_indices = filled_unique_homologues(profile)
        for seq, source_index in zip(homologue_seqs, source_indices):
            candidates.append(
                ("mmseqs_homologue", int(source_index), np.asarray(seq, dtype=np.int64))
            )

    target_margin = float(anchor_cfg.get("target_effective_stability_margin", 1.0))
    rows: list[dict[str, Any]] = []
    best_row: dict[str, Any] | None = None
    best_seq: np.ndarray | None = None
    for source, source_index, seq in candidates:
        comp = components_for_seq(landscape, seq)
        margin = float(comp.get("effective_stability_margin_used", comp["stability_margin_used"]))
        row = {
            "anchor_source": source,
            "source_index": source_index,
            "fitness": evaluate_fitness(landscape, seq),
            "stability": float(comp["stability"]),
            "effective_stability": float(comp["effective_stability"]),
            "stability_margin_used": float(comp["stability_margin_used"]),
            "effective_stability_margin_used": margin,
            "target_effective_stability_margin": target_margin,
            "margin_distance": abs(margin - target_margin),
            "sequence": seq_to_string(seq),
        }
        rows.append(row)
        if mode == "reference":
            best_row = row
            best_seq = seq.copy()
            break
        if best_row is None or row["margin_distance"] < best_row["margin_distance"]:
            best_row = row
            best_seq = seq.copy()

    if best_row is None or best_seq is None:
        raise RuntimeError("No anchor sequence candidates were available")

    max_reported = int(anchor_cfg.get("max_reported_candidates", 5000))
    rows.sort(key=lambda row: float(row["margin_distance"]))
    best_row = dict(best_row)
    best_row["selected_anchor"] = True
    best_row["anchor_selection_mode"] = mode
    return best_seq, best_row, rows[:max_reported]


def weighted_state_counts(
    profile: AlignmentProfile,
    positions: tuple[int, ...],
) -> tuple[Counter[tuple[int, ...]], float]:
    seq_matrix = np.asarray(profile.seq_matrix, dtype=np.int64)
    weights = np.asarray(getattr(profile, "weights", np.ones(seq_matrix.shape[0])), dtype=np.float64)
    pos = np.asarray(positions, dtype=np.int64)
    observed = seq_matrix[:, pos]
    valid = np.all(observed >= 0, axis=1)
    counts: Counter[tuple[int, ...]] = Counter()
    total = 0.0
    for state, weight in zip(observed[valid], weights[valid]):
        counts[tuple(int(x) for x in state.tolist())] += float(weight)
        total += float(weight)
    return counts, total


def marginal_frequency(
    profile: AlignmentProfile,
    position: int,
    aa: int,
) -> float:
    seq_matrix = np.asarray(profile.seq_matrix, dtype=np.int64)
    weights = np.asarray(getattr(profile, "weights", np.ones(seq_matrix.shape[0])), dtype=np.float64)
    observed = seq_matrix[:, int(position)]
    valid = observed >= 0
    if not np.any(valid):
        return 0.0
    total = float(np.sum(weights[valid]))
    if total <= 0:
        return 0.0
    return float(np.sum(weights[valid & (observed == int(aa))]) / total)


def dca_graph(
    profile: AlignmentProfile,
    *,
    scan_top_couplings: int,
    max_positions: int,
) -> tuple[dict[int, set[int]], dict[tuple[int, int], float], list[int]]:
    top = profile.get_top_couplings(n_pairs=int(scan_top_couplings))
    edge_scores: dict[tuple[int, int], float] = {}
    weighted_degree: defaultdict[int, float] = defaultdict(float)
    for pos_i, pos_j, score in top:
        score = float(score)
        if score <= 0:
            continue
        i, j = sorted((int(pos_i), int(pos_j)))
        edge_scores[(i, j)] = max(score, edge_scores.get((i, j), 0.0))
        weighted_degree[i] += score
        weighted_degree[j] += score
    allowed = {
        pos
        for pos, _score in sorted(weighted_degree.items(), key=lambda item: -item[1])[: int(max_positions)]
    }
    adjacency: dict[int, set[int]] = defaultdict(set)
    filtered_scores: dict[tuple[int, int], float] = {}
    for (i, j), score in edge_scores.items():
        if i not in allowed or j not in allowed:
            continue
        adjacency[i].add(j)
        adjacency[j].add(i)
        filtered_scores[(i, j)] = score
    vertices = sorted(adjacency)
    return adjacency, filtered_scores, vertices


def edge_score(edge_scores: dict[tuple[int, int], float], i: int, j: int) -> float:
    a, b = sorted((int(i), int(j)))
    return float(edge_scores.get((a, b), 0.0))


def enumerate_cliques(
    *,
    order: int,
    adjacency: dict[int, set[int]],
    edge_scores: dict[tuple[int, int], float],
    vertices: list[int],
    max_candidates: int,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def extend(prefix: tuple[int, ...], candidates: list[int]) -> None:
        if len(out) >= max_candidates:
            return
        if len(prefix) == order:
            scores = [
                edge_score(edge_scores, i, j)
                for i, j in itertools.combinations(prefix, 2)
            ]
            out.append(
                {
                    "positions": tuple(int(x) for x in prefix),
                    "edge_count": len(scores),
                    "mean_edge_score": float(np.mean(scores)) if scores else 0.0,
                    "min_edge_score": float(np.min(scores)) if scores else 0.0,
                    "max_edge_score": float(np.max(scores)) if scores else 0.0,
                    "clique_score": float(np.mean(scores) * np.min(scores)) if scores else 0.0,
                }
            )
            return
        for idx, vertex in enumerate(candidates):
            if any(vertex not in adjacency[p] for p in prefix):
                continue
            next_candidates = [
                candidate
                for candidate in candidates[idx + 1 :]
                if candidate in adjacency[vertex]
                and all(candidate in adjacency[p] for p in prefix)
            ]
            if len(prefix) + 1 + len(next_candidates) < order:
                continue
            extend((*prefix, int(vertex)), next_candidates)

    extend(tuple(), vertices)
    out.sort(key=lambda row: (-float(row["clique_score"]), row["positions"]))
    return out


def haplotype_for_clique(
    *,
    profile: AlignmentProfile,
    anchor_seq: np.ndarray,
    positions: tuple[int, ...],
    min_frequency: float,
    min_log2_enrichment: float,
) -> dict[str, Any] | None:
    counts, total = weighted_state_counts(profile, positions)
    if total <= 0:
        return None
    candidates: list[dict[str, Any]] = []
    anchor_state = tuple(int(anchor_seq[pos]) for pos in positions)
    for state, count in counts.items():
        if state == anchor_state:
            continue
        if any(int(state[idx]) == int(anchor_state[idx]) for idx in range(len(positions))):
            continue
        frequency = float(count / total)
        if frequency < min_frequency:
            continue
        expected = 1.0
        for pos, aa in zip(positions, state):
            expected *= max(marginal_frequency(profile, int(pos), int(aa)), 1e-12)
        log2_enrichment = float(math.log2(max(frequency, 1e-12) / max(expected, 1e-12)))
        if log2_enrichment < min_log2_enrichment:
            continue
        candidates.append(
            {
                "haplotype_state": tuple(int(x) for x in state),
                "haplotype_frequency": frequency,
                "haplotype_weight": float(count),
                "expected_frequency": float(expected),
                "log2_enrichment": log2_enrichment,
            }
        )
    if not candidates:
        return None
    candidates.sort(
        key=lambda row: (
            -float(row["log2_enrichment"]),
            -float(row["haplotype_frequency"]),
            row["haplotype_state"],
        )
    )
    return candidates[0]


def single_mutant_delta(
    landscape: FitnessLandscapeEnv,
    anchor_seq: np.ndarray,
    parent_fitness: float,
    *,
    position: int,
    alt_aa: int,
    cache: dict[tuple[int, int], float],
) -> float:
    key = (int(position), int(alt_aa))
    if key in cache:
        return cache[key]
    mutant = np.asarray(anchor_seq, dtype=np.int64).copy()
    mutant[int(position)] = int(alt_aa)
    delta = evaluate_fitness(landscape, mutant) - parent_fitness
    cache[key] = float(delta)
    return float(delta)


def context_variable(
    cfg: dict[str, Any],
    *,
    profile: AlignmentProfile,
    landscape: FitnessLandscapeEnv,
    anchor_seq: np.ndarray,
    forbidden_positions: set[int],
    delta_cache: dict[tuple[int, int], float],
) -> dict[str, Any]:
    context_cfg = dict(cfg.get("context_variable", {}))
    min_frequency = float(context_cfg.get("min_frequency", 0.02))
    target_delta = float(context_cfg.get("target_single_delta", -20.0))
    delta_tolerance = float(context_cfg.get("single_delta_tolerance", 80.0))
    max_abs_delta = float(context_cfg.get("max_abs_single_delta", math.inf))
    seq_matrix = np.asarray(profile.seq_matrix, dtype=np.int64)
    weights = np.asarray(getattr(profile, "weights", np.ones(seq_matrix.shape[0])), dtype=np.float64)
    parent_fitness = evaluate_fitness(landscape, anchor_seq)
    candidates: list[dict[str, Any]] = []
    for pos in range(seq_matrix.shape[1]):
        if pos in forbidden_positions:
            continue
        observed = seq_matrix[:, pos]
        valid = observed >= 0
        total = float(np.sum(weights[valid]))
        if total <= 0:
            continue
        for aa in range(len(AA_ORDER)):
            if aa == int(anchor_seq[pos]):
                continue
            frequency = float(np.sum(weights[valid & (observed == aa)]) / total)
            if frequency < min_frequency:
                continue
            delta = single_mutant_delta(
                landscape,
                anchor_seq,
                parent_fitness,
                position=pos,
                alt_aa=aa,
                cache=delta_cache,
            )
            if abs(delta) > max_abs_delta:
                continue
            score = frequency * math.exp(-abs(delta - target_delta) / max(delta_tolerance, 1e-9))
            candidates.append(
                {
                    "context_position": int(pos),
                    "context_anchor_aa": int(anchor_seq[pos]),
                    "context_alt_aa": int(aa),
                    "context_frequency": frequency,
                    "context_single_delta": float(delta),
                    "context_score": float(score),
                }
            )
    if not candidates:
        raise RuntimeError("Could not select a context variable for clique map")
    candidates.sort(key=lambda row: (-float(row["context_score"]), row["context_position"], row["context_alt_aa"]))
    return candidates[0]


def coefficient_sign(term_index: int, seed: int, mode: str) -> int:
    mode = str(mode).lower()
    if mode == "positive":
        return 1
    if mode == "negative":
        return -1
    if mode == "alternating":
        return 1 if (int(term_index) + int(seed)) % 2 == 0 else -1
    raise ValueError(f"Unsupported coefficient_sign_mode: {mode}")


def select_higher_order_terms(
    cfg: dict[str, Any],
    *,
    seed: int,
    profile: AlignmentProfile,
    landscape: FitnessLandscapeEnv,
    anchor_seq: np.ndarray,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    term_cfg = dict(cfg["dca_haplotype_terms"])
    orders = [int(order) for order in term_cfg.get("orders", [2, 3, 4, 5])]
    max_order = max(orders)
    adjacency, edge_scores, vertices = dca_graph(
        profile,
        scan_top_couplings=int(term_cfg.get("scan_top_couplings", 3000)),
        max_positions=int(term_cfg.get("max_clique_positions", 90)),
    )
    max_candidates = int(term_cfg.get("max_clique_candidates_per_order", 5000))
    maps_per_order = int(term_cfg.get("maps_per_order", 4))
    min_frequency_by_order = {
        int(key): float(value)
        for key, value in dict(term_cfg.get("min_haplotype_frequency_by_order", {})).items()
    }
    min_frequency_default = float(term_cfg.get("min_haplotype_frequency", 0.0005))
    min_log2 = float(term_cfg.get("min_haplotype_log2_enrichment", 0.0))
    sign_mode = str(term_cfg.get("coefficient_sign_mode", "alternating"))

    parent_fitness = evaluate_fitness(landscape, anchor_seq)
    delta_cache: dict[tuple[int, int], float] = {}
    terms: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []
    map_specs: list[dict[str, Any]] = []

    for order in orders:
        cliques = enumerate_cliques(
            order=order,
            adjacency=adjacency,
            edge_scores=edge_scores,
            vertices=vertices,
            max_candidates=max_candidates,
        )
        selected_for_order = 0
        for clique_rank, clique in enumerate(cliques):
            positions = tuple(int(pos) for pos in clique["positions"])
            min_freq = min_frequency_by_order.get(int(order), min_frequency_default)
            hap = haplotype_for_clique(
                profile=profile,
                anchor_seq=anchor_seq,
                positions=positions,
                min_frequency=min_freq,
                min_log2_enrichment=min_log2,
            )
            if hap is None:
                if len(rejected_rows) < 500:
                    rejected_rows.append(
                        {
                            "seed": int(seed),
                            "order": int(order),
                            "clique_rank": int(clique_rank),
                            "positions": ";".join(str(pos) for pos in positions),
                            "reason": "no_observed_all-alt_enriched_haplotype",
                            **{k: clique[k] for k in ("mean_edge_score", "min_edge_score", "max_edge_score", "clique_score")},
                        }
                    )
                continue
            states = tuple(int(x) for x in hap["haplotype_state"])
            single_deltas = [
                single_mutant_delta(
                    landscape,
                    anchor_seq,
                    parent_fitness,
                    position=pos,
                    alt_aa=aa,
                    cache=delta_cache,
                )
                for pos, aa in zip(positions, states)
            ]
            context = context_variable(
                cfg,
                profile=profile,
                landscape=landscape,
                anchor_seq=anchor_seq,
                forbidden_positions=set(positions),
                delta_cache=delta_cache,
            )
            term_index = len(terms)
            sign = coefficient_sign(term_index, seed, sign_mode)
            term = {
                "term_id": f"seed{seed}_order{order}_term{selected_for_order}",
                "seed": int(seed),
                "order": int(order),
                "term_index": int(term_index),
                "term_within_order": int(selected_for_order),
                "coefficient_sign": int(sign),
                "positions": positions,
                "amino_acids": states,
                "mutations": tuple(
                    mutation_name(pos, int(anchor_seq[pos]), aa)
                    for pos, aa in zip(positions, states)
                ),
                "single_deltas": tuple(float(x) for x in single_deltas),
                "mean_single_delta": float(np.mean(single_deltas)),
                "min_single_delta": float(np.min(single_deltas)),
                "max_single_delta": float(np.max(single_deltas)),
                **clique,
                **hap,
            }
            terms.append(term)
            map_specs.append(
                {
                    "map_id": f"seed{seed}_order{order}_map{selected_for_order}",
                    "seed": int(seed),
                    "order": int(order),
                    "term_id": term["term_id"],
                    "target_term_index": int(term_index),
                    "positions": tuple((*positions, int(context["context_position"]))),
                    "amino_acids": tuple((*states, int(context["context_alt_aa"]))),
                    "mutations": tuple(
                        (*term["mutations"], mutation_name(int(context["context_position"]), int(context["context_anchor_aa"]), int(context["context_alt_aa"])))
                    ),
                    **context,
                }
            )
            selected_for_order += 1
            if selected_for_order >= maps_per_order:
                break
        if selected_for_order < maps_per_order:
            raise RuntimeError(
                f"Only selected {selected_for_order} order-{order} DCA haplotype terms; "
                f"required {maps_per_order}. Relax clique/haplotype filters."
            )
    return terms, map_specs, rejected_rows


def sequence_from_bits(
    anchor_seq: np.ndarray,
    positions: tuple[int, ...],
    amino_acids: tuple[int, ...],
    bits: np.ndarray,
) -> np.ndarray:
    seq = np.asarray(anchor_seq, dtype=np.int64).copy()
    for pos, aa, bit in zip(positions, amino_acids, bits):
        if int(bit):
            seq[int(pos)] = int(aa)
    return seq


def term_score(seq: np.ndarray, terms: list[dict[str, Any]], *, strength: float) -> tuple[float, list[str]]:
    score = 0.0
    active: list[str] = []
    for term in terms:
        positions = term["positions"]
        aas = term["amino_acids"]
        if all(int(seq[pos]) == int(aa) for pos, aa in zip(positions, aas)):
            contribution = float(strength) * float(term["coefficient_sign"])
            score += contribution
            active.append(str(term["term_id"]))
    return float(score), active


def walsh_decomposition(values: np.ndarray, *, dim: int) -> tuple[list[dict[str, Any]], dict[str, float]]:
    y = np.asarray(values, dtype=np.float64)
    n = 1 << dim
    if y.shape != (n,):
        raise ValueError(f"Expected {n} values for dim={dim}, got {y.shape}")
    rows: list[dict[str, Any]] = []
    variance_by_order: defaultdict[int, float] = defaultdict(float)
    total_variance = float(np.var(y))
    for mask in range(1 << dim):
        signs = np.ones(n, dtype=np.float64)
        for node in range(n):
            parity = (node & mask).bit_count() % 2
            signs[node] = -1.0 if parity else 1.0
        coeff = float(np.mean(y * signs))
        order = int(mask.bit_count())
        variance_by_order[order] += coeff * coeff
        if mask == 0:
            term = "intercept"
        else:
            term = ",".join(str(idx) for idx in range(dim) if (mask >> idx) & 1)
        rows.append(
            {
                "order": order,
                "mask": int(mask),
                "term": term,
                "coefficient": coeff,
                "abs_coefficient": abs(coeff),
            }
        )
    variance_explained = {
        str(order): (float(value / total_variance) if total_variance > 1e-12 and order > 0 else 0.0)
        for order, value in sorted(variance_by_order.items())
        if order > 0
    }
    rows.sort(key=lambda row: (int(row["order"]), -float(row["abs_coefficient"]), int(row["mask"])))
    return rows, variance_explained


def sign_epistasis_counts(fitness: np.ndarray, *, dim: int, eps: float = 1e-9) -> dict[str, Any]:
    sign_count = 0
    reciprocal_count = 0
    total = 0
    pair_rows: list[dict[str, Any]] = []
    for var_i, var_j in itertools.combinations(range(dim), 2):
        other_vars = [idx for idx in range(dim) if idx not in {var_i, var_j}]
        pair_sign = 0
        pair_recip = 0
        pair_total = 0
        for bg_bits in itertools.product([0, 1], repeat=len(other_vars)):
            base = 0
            for idx, bit in zip(other_vars, bg_bits):
                if bit:
                    base |= 1 << idx
            f00 = float(fitness[base])
            f10 = float(fitness[base | (1 << var_i)])
            f01 = float(fitness[base | (1 << var_j)])
            f11 = float(fitness[base | (1 << var_i) | (1 << var_j)])
            effect_i = f10 - f00
            effect_j = f01 - f00
            effect_i_on_j = f11 - f01
            effect_j_on_i = f11 - f10
            sign_i = effect_i * effect_i_on_j < -eps
            sign_j = effect_j * effect_j_on_i < -eps
            recip = bool(sign_i and sign_j)
            sign_epi = bool(sign_i or sign_j)
            pair_sign += int(sign_epi)
            pair_recip += int(recip)
            pair_total += 1
        sign_count += pair_sign
        reciprocal_count += pair_recip
        total += pair_total
        pair_rows.append(
            {
                "var_i": int(var_i),
                "var_j": int(var_j),
                "backgrounds": int(pair_total),
                "sign_epistasis_count": int(pair_sign),
                "reciprocal_sign_epistasis_count": int(pair_recip),
                "sign_epistasis_fraction": float(pair_sign / pair_total) if pair_total else 0.0,
                "reciprocal_sign_epistasis_fraction": float(pair_recip / pair_total) if pair_total else 0.0,
            }
        )
    return {
        "total_pair_backgrounds": int(total),
        "sign_epistasis_count": int(sign_count),
        "reciprocal_sign_epistasis_count": int(reciprocal_count),
        "sign_epistasis_fraction": float(sign_count / total) if total else 0.0,
        "reciprocal_sign_epistasis_fraction": float(reciprocal_count / total) if total else 0.0,
        "pair_rows": pair_rows,
    }


def local_extrema_counts(fitness: np.ndarray, *, dim: int) -> dict[str, Any]:
    local_maxima = []
    local_minima = []
    for node in range(1 << dim):
        neighbors = [node ^ (1 << bit) for bit in range(dim)]
        values = fitness[np.asarray(neighbors, dtype=np.int64)]
        if np.all(fitness[node] >= values):
            local_maxima.append(int(node))
        if np.all(fitness[node] <= values):
            local_minima.append(int(node))
    return {
        "n_local_maxima": int(len(local_maxima)),
        "n_local_minima": int(len(local_minima)),
        "local_maxima": local_maxima,
        "local_minima": local_minima,
    }


def run_hamming_map(
    *,
    cfg: dict[str, Any],
    landscape: FitnessLandscapeEnv,
    anchor_seq: np.ndarray,
    terms: list[dict[str, Any]],
    map_spec: dict[str, Any],
    strength: float,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    positions = tuple(int(x) for x in map_spec["positions"])
    amino_acids = tuple(int(x) for x in map_spec["amino_acids"])
    dim = len(positions)
    target_order = int(map_spec["order"])
    node_count = 1 << dim
    base_fitness = np.zeros(node_count, dtype=np.float64)
    high_order_scores = np.zeros(node_count, dtype=np.float64)
    fitness = np.zeros(node_count, dtype=np.float64)
    node_rows: list[dict[str, Any]] = []
    for mask in range(node_count):
        bits = np.asarray([(mask >> idx) & 1 for idx in range(dim)], dtype=np.int8)
        seq = sequence_from_bits(anchor_seq, positions, amino_acids, bits)
        base = evaluate_fitness(landscape, seq)
        high_order, active_terms = term_score(seq, terms, strength=strength)
        base_fitness[mask] = base
        high_order_scores[mask] = high_order
        fitness[mask] = base + high_order
        mutations = [
            mutation
            for mutation, bit in zip(map_spec["mutations"], bits)
            if int(bit)
        ]
        node_rows.append(
            {
                "map_id": map_spec["map_id"],
                "seed": int(map_spec["seed"]),
                "expected_order": target_order,
                "strength": float(strength),
                "node_index": int(mask),
                "bitstring_lsb_first": "".join(str(int(x)) for x in bits.tolist()),
                "hamming_weight": int(np.sum(bits)),
                "mutations": ";".join(mutations),
                "base_fitness": float(base),
                "high_order_score": float(high_order),
                "fitness": float(base + high_order),
                "active_high_order_terms": ";".join(active_terms),
                "amino_acid_sequence": seq_to_string(seq),
            }
        )

    walsh_rows, variance_explained = walsh_decomposition(fitness, dim=dim)
    base_walsh_rows, base_variance_explained = walsh_decomposition(base_fitness, dim=dim)
    target_mask = (1 << target_order) - 1
    target_coeff = next(row for row in walsh_rows if int(row["mask"]) == target_mask)
    base_target_coeff = next(row for row in base_walsh_rows if int(row["mask"]) == target_mask)
    sign_counts = sign_epistasis_counts(fitness, dim=dim)
    extrema = local_extrema_counts(fitness, dim=dim)
    total_variance = float(np.var(fitness))
    high_order_variance = float(np.var(high_order_scores))
    base_variance = float(np.var(base_fitness))
    summary = {
        "map_id": map_spec["map_id"],
        "seed": int(map_spec["seed"]),
        "expected_order": target_order,
        "grid_dimension": int(dim),
        "node_count": int(node_count),
        "strength": float(strength),
        "term_id": map_spec["term_id"],
        "target_term_index": int(map_spec["target_term_index"]),
        "positions": ";".join(str(pos) for pos in positions),
        "mutations": ";".join(map_spec["mutations"]),
        "context_position": int(map_spec["context_position"]),
        "context_mutation": map_spec["mutations"][-1],
        "fitness_min": float(np.min(fitness)),
        "fitness_mean": float(np.mean(fitness)),
        "fitness_max": float(np.max(fitness)),
        "fitness_sd": float(np.std(fitness)),
        "base_variance": base_variance,
        "high_order_overlay_variance": high_order_variance,
        "total_variance": total_variance,
        "target_order_abs_walsh_coefficient": float(target_coeff["abs_coefficient"]),
        "base_target_order_abs_walsh_coefficient": float(base_target_coeff["abs_coefficient"]),
        "target_order_abs_walsh_lift": float(target_coeff["abs_coefficient"] - base_target_coeff["abs_coefficient"]),
        "target_order_variance_fraction": float(
            (target_coeff["coefficient"] ** 2) / total_variance
        )
        if total_variance > 1e-12
        else 0.0,
        "order_variance_fraction": float(variance_explained.get(str(target_order), 0.0)),
        "base_order_variance_fraction": float(base_variance_explained.get(str(target_order), 0.0)),
        "order_variance_fraction_lift": float(
            variance_explained.get(str(target_order), 0.0)
            - base_variance_explained.get(str(target_order), 0.0)
        ),
        **{f"walsh_variance_order_{order}": float(value) for order, value in variance_explained.items()},
        **{f"base_walsh_variance_order_{order}": float(value) for order, value in base_variance_explained.items()},
        **{key: value for key, value in sign_counts.items() if key != "pair_rows"},
        **extrema,
    }
    for row in walsh_rows:
        row.update(
            {
                "map_id": map_spec["map_id"],
                "seed": int(map_spec["seed"]),
                "expected_order": target_order,
                "strength": float(strength),
            }
        )
    pair_rows = []
    for row in sign_counts["pair_rows"]:
        row.update(
            {
                "map_id": map_spec["map_id"],
                "seed": int(map_spec["seed"]),
                "expected_order": target_order,
                "strength": float(strength),
                "mutation_i": map_spec["mutations"][int(row["var_i"])],
                "mutation_j": map_spec["mutations"][int(row["var_j"])],
            }
        )
        pair_rows.append(row)
    return summary, node_rows, walsh_rows, pair_rows, []


def run_experiment(cfg: dict[str, Any], *, config_path: Path, project_root: Path, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    profile = build_alignment_profile(
        cfg,
        project_root=project_root,
        config_path=config_path,
        output_dir=output_dir,
    )
    seeds = [int(seed) for seed in cfg["seeds"]]
    strengths = [float(x) for x in cfg["dca_haplotype_terms"]["strength_grid"]]

    all_anchor_rows: list[dict[str, Any]] = []
    all_anchor_candidates: list[dict[str, Any]] = []
    all_terms: list[dict[str, Any]] = []
    all_rejected: list[dict[str, Any]] = []
    all_map_specs: list[dict[str, Any]] = []
    all_map_summaries: list[dict[str, Any]] = []
    all_nodes: list[dict[str, Any]] = []
    all_walsh: list[dict[str, Any]] = []
    all_sign_pairs: list[dict[str, Any]] = []

    for seed in seeds:
        landscape = build_base_landscape(cfg, seed=seed, profile=profile)
        anchor_seq, anchor_row, anchor_candidates = select_anchor_sequence(
            cfg,
            profile=profile,
            landscape=landscape,
        )
        anchor_row["seed"] = int(seed)
        all_anchor_rows.append(anchor_row)
        for row in anchor_candidates:
            candidate = dict(row)
            candidate["seed"] = int(seed)
            all_anchor_candidates.append(candidate)

        terms, map_specs, rejected = select_higher_order_terms(
            cfg,
            seed=seed,
            profile=profile,
            landscape=landscape,
            anchor_seq=anchor_seq,
        )
        all_terms.extend(terms)
        all_map_specs.extend(map_specs)
        all_rejected.extend(rejected)
        for strength in strengths:
            for map_spec in map_specs:
                summary, node_rows, walsh_rows, sign_pair_rows, _ = run_hamming_map(
                    cfg=cfg,
                    landscape=landscape,
                    anchor_seq=anchor_seq,
                    terms=terms,
                    map_spec=map_spec,
                    strength=strength,
                )
                all_map_summaries.append(summary)
                all_nodes.extend(node_rows)
                all_walsh.extend(walsh_rows)
                all_sign_pairs.extend(sign_pair_rows)

    term_rows = []
    for term in all_terms:
        term_rows.append(
            {
                "term_id": term["term_id"],
                "seed": term["seed"],
                "order": term["order"],
                "term_index": term["term_index"],
                "term_within_order": term["term_within_order"],
                "coefficient_sign": term["coefficient_sign"],
                "positions": ";".join(str(pos) for pos in term["positions"]),
                "amino_acids": ";".join(aa_name(aa) for aa in term["amino_acids"]),
                "mutations": ";".join(term["mutations"]),
                "haplotype_frequency": term["haplotype_frequency"],
                "haplotype_weight": term["haplotype_weight"],
                "expected_frequency": term["expected_frequency"],
                "log2_enrichment": term["log2_enrichment"],
                "mean_edge_score": term["mean_edge_score"],
                "min_edge_score": term["min_edge_score"],
                "max_edge_score": term["max_edge_score"],
                "clique_score": term["clique_score"],
                "mean_single_delta": term["mean_single_delta"],
                "min_single_delta": term["min_single_delta"],
                "max_single_delta": term["max_single_delta"],
                "single_deltas": ";".join(str(x) for x in term["single_deltas"]),
            }
        )

    map_spec_rows = []
    for spec in all_map_specs:
        map_spec_rows.append(
            {
                "map_id": spec["map_id"],
                "seed": spec["seed"],
                "order": spec["order"],
                "term_id": spec["term_id"],
                "target_term_index": spec["target_term_index"],
                "grid_dimension": len(spec["positions"]),
                "positions": ";".join(str(pos) for pos in spec["positions"]),
                "amino_acids": ";".join(aa_name(aa) for aa in spec["amino_acids"]),
                "mutations": ";".join(spec["mutations"]),
                "context_position": spec["context_position"],
                "context_alt_aa": aa_name(spec["context_alt_aa"]),
                "context_frequency": spec["context_frequency"],
                "context_single_delta": spec["context_single_delta"],
            }
        )

    write_csv(output_dir / "anchor_summary.csv", all_anchor_rows)
    write_csv(output_dir / "anchor_candidates.csv", all_anchor_candidates)
    write_csv(output_dir / "dca_haplotype_terms.csv", term_rows)
    write_csv(output_dir / "dca_clique_rejections.csv", all_rejected)
    write_csv(output_dir / "hamming_map_specs.csv", map_spec_rows)
    write_csv(output_dir / "hamming_map_summary.csv", all_map_summaries)
    write_csv(output_dir / "hamming_nodes.csv", all_nodes)
    write_csv(output_dir / "hamming_walsh_terms.csv", all_walsh)
    write_csv(output_dir / "hamming_sign_epistasis_pairs.csv", all_sign_pairs)

    def grouped(rows: list[dict[str, Any]], *keys: str) -> dict[tuple[Any, ...], list[dict[str, Any]]]:
        out: defaultdict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            out[tuple(row[key] for key in keys)].append(row)
        return dict(out)

    order_strength_summary = []
    for (order, strength), rows in sorted(grouped(all_map_summaries, "expected_order", "strength").items()):
        target_coeffs = [float(row["target_order_abs_walsh_coefficient"]) for row in rows]
        target_lifts = [float(row["target_order_abs_walsh_lift"]) for row in rows]
        order_fracs = [float(row["order_variance_fraction"]) for row in rows]
        sign_fracs = [float(row["sign_epistasis_fraction"]) for row in rows]
        rse_fracs = [float(row["reciprocal_sign_epistasis_fraction"]) for row in rows]
        local_maxima = [float(row["n_local_maxima"]) for row in rows]
        order_strength_summary.append(
            {
                "order": int(order),
                "strength": float(strength),
                "n_maps": int(len(rows)),
                "mean_target_order_abs_walsh_coefficient": float(np.mean(target_coeffs)),
                "mean_target_order_abs_walsh_lift": float(np.mean(target_lifts)),
                "mean_order_variance_fraction": float(np.mean(order_fracs)),
                "mean_sign_epistasis_fraction": float(np.mean(sign_fracs)),
                "mean_reciprocal_sign_epistasis_fraction": float(np.mean(rse_fracs)),
                "mean_local_maxima": float(np.mean(local_maxima)),
            }
        )
    write_csv(output_dir / "order_strength_summary.csv", order_strength_summary)

    max_order = max(int(row["expected_order"]) for row in all_map_summaries)
    max_strength = max(strengths)
    upper = [
        row
        for row in all_map_summaries
        if int(row["expected_order"]) == max_order and float(row["strength"]) == max_strength
    ]
    baseline = [
        row
        for row in all_map_summaries
        if int(row["expected_order"]) == max_order and float(row["strength"]) == min(strengths)
    ]
    summary = {
        "experiment_id": cfg.get("experiment_id"),
        "reference": cfg["reference"],
        "alignment": {
            "n_sequences": int(getattr(profile, "n_sequences", len(profile.seq_matrix))),
            "n_effective_sequences": float(getattr(profile, "n_effective_sequences", np.nan)),
            "length": int(profile.seq_matrix.shape[1]),
        },
        "seeds": seeds,
        "strength_grid": strengths,
        "orders": [int(order) for order in cfg["dca_haplotype_terms"]["orders"]],
        "maps_per_order": int(cfg["dca_haplotype_terms"]["maps_per_order"]),
        "grid_dimension_rule": "dense 2^(n+1) binary Hamming grid for expected DCA order n",
        "n_terms": int(len(all_terms)),
        "n_maps": int(len(all_map_summaries)),
        "n_node_rows": int(len(all_nodes)),
        "n_walsh_rows": int(len(all_walsh)),
        "old_epistasis_disabled": {
            "n_epistatic_terms": 0,
            "epistasis_strength": 0.0,
            "epistasis_covariation_bias": 0.0,
            "n_higher_order_epistatic_terms": 0,
        },
        "upper_order": {
            "order": int(max_order),
            "strength": float(max_strength),
            "mean_target_order_abs_walsh_coefficient": float(np.mean([row["target_order_abs_walsh_coefficient"] for row in upper])),
            "mean_target_order_abs_walsh_lift": float(np.mean([row["target_order_abs_walsh_lift"] for row in upper])),
            "mean_order_variance_fraction": float(np.mean([row["order_variance_fraction"] for row in upper])),
            "mean_sign_epistasis_fraction": float(np.mean([row["sign_epistasis_fraction"] for row in upper])),
            "mean_reciprocal_sign_epistasis_fraction": float(np.mean([row["reciprocal_sign_epistasis_fraction"] for row in upper])),
        },
        "baseline_upper_order": {
            "order": int(max_order),
            "strength": float(min(strengths)),
            "mean_target_order_abs_walsh_coefficient": float(np.mean([row["target_order_abs_walsh_coefficient"] for row in baseline])),
            "mean_order_variance_fraction": float(np.mean([row["order_variance_fraction"] for row in baseline])),
        },
        "outputs": {
            "anchor_summary": str(output_dir / "anchor_summary.csv"),
            "terms": str(output_dir / "dca_haplotype_terms.csv"),
            "map_summary": str(output_dir / "hamming_map_summary.csv"),
            "walsh_terms": str(output_dir / "hamming_walsh_terms.csv"),
            "sign_epistasis_pairs": str(output_dir / "hamming_sign_epistasis_pairs.csv"),
            "order_strength_summary": str(output_dir / "order_strength_summary.csv"),
        },
    }
    write_json(output_dir / "summary.json", summary)
    return summary


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
    print(json.dumps(summary, indent=2, default=json_default))


if __name__ == "__main__":
    main()
