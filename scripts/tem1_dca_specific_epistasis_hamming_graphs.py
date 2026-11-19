#!/usr/bin/env python3
"""Run a TEM-1 DCA-prior specific-epistasis dense Hamming diagnostic."""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
import types
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
LANDSCAPY_SRC = PROJECT_ROOT / "external" / "landscapy" / "src"
SCRIPT_DIR = Path(__file__).resolve().parent
for path in (ADAPT_ENV_ROOT, LANDSCAPY_SRC, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def _install_landscapy_optional_dependency_stubs() -> None:
    """Stub optional landscapy imports that this binary-Hamming path does not use."""

    def module_available(name: str) -> bool:
        try:
            return importlib.util.find_spec(name) is not None
        except (ImportError, ValueError, AttributeError):
            return False

    def missing(*_args: Any, **_kwargs: Any) -> Any:
        raise ImportError(
            "This EXP-023 runner uses only landscapy's dense Hamming APIs; "
            "the requested optional phylo/PLM dependency is not installed."
        )

    if not module_available("torch"):
        torch_mod = types.ModuleType("torch")

        class Tensor:  # noqa: D401 - minimal runtime sentinel
            """Placeholder used only for isinstance checks."""

            def detach(self) -> "Tensor":
                return self

            def cpu(self) -> "Tensor":
                return self

            def numpy(self) -> np.ndarray:
                return np.asarray([])

        torch_mod.Tensor = Tensor
        torch_mod.float32 = "float32"
        torch_mod.long = "long"
        torch_mod.tensor = missing
        torch_mod.zeros = missing
        torch_mod.as_tensor = missing
        torch_mod.stack = missing
        torch_mod.cuda = types.SimpleNamespace(is_available=lambda: False)
        sys.modules["torch"] = torch_mod

    if not module_available("torch_geometric"):
        tg_mod = types.ModuleType("torch_geometric")
        tg_data = types.ModuleType("torch_geometric.data")
        tg_utils = types.ModuleType("torch_geometric.utils")

        class Data:  # noqa: D401 - minimal runtime sentinel
            """Placeholder for unused PyG exports."""

        tg_data.Data = Data
        tg_utils.from_networkx = missing
        tg_mod.data = tg_data
        tg_mod.utils = tg_utils
        sys.modules["torch_geometric"] = tg_mod
        sys.modules["torch_geometric.data"] = tg_data
        sys.modules["torch_geometric.utils"] = tg_utils

    if not module_available("cogent3"):
        cogent_mod = types.ModuleType("cogent3")
        cogent_core = types.ModuleType("cogent3.core")
        cogent_seq = types.ModuleType("cogent3.core.sequence")
        cogent_moltype = types.ModuleType("cogent3.core.moltype")
        cogent_alignment = types.ModuleType("cogent3.core.alignment")
        cogent_tree = types.ModuleType("cogent3.core.tree")

        class Sequence:  # noqa: D401 - minimal runtime sentinel
            """Placeholder for unused cogent3 sequence objects."""

            name = None

            def __str__(self) -> str:
                return getattr(self, "_seq", "")

        class MolType:  # noqa: D401 - minimal runtime sentinel
            """Placeholder for unused cogent3 moltype objects."""

            alphabet = tuple("ACDEFGHIKLMNPQRSTVWYU")

            def make_seq(self, seq: str) -> Sequence:
                out = Sequence()
                out.name = None
                out.moltype = self
                out._seq = str(seq)
                return out

        class Alignment:  # noqa: D401 - minimal runtime sentinel
            """Placeholder for unused cogent3 alignment objects."""

        class PhyloNode:  # noqa: D401 - minimal runtime sentinel
            """Placeholder for unused cogent3 tree objects."""

        cogent_seq.Sequence = Sequence
        cogent_moltype.MolType = MolType
        cogent_alignment.Alignment = Alignment
        cogent_alignment.make_aligned_seqs = missing
        cogent_tree.PhyloNode = PhyloNode
        cogent_mod.Sequence = Sequence
        cogent_mod.MolType = MolType
        cogent_mod.Alignment = Alignment
        cogent_mod.PhyloNode = PhyloNode
        cogent_mod.core = cogent_core
        cogent_core.sequence = cogent_seq
        cogent_core.moltype = cogent_moltype
        cogent_core.alignment = cogent_alignment
        cogent_core.tree = cogent_tree
        cogent_mod.get_moltype = lambda _name=None: MolType()
        cogent_mod.load_unaligned_seqs = missing
        cogent_mod.load_aligned_seqs = missing
        cogent_mod.load_tree = missing
        cogent_mod.get_app = missing
        sys.modules["cogent3"] = cogent_mod
        sys.modules["cogent3.core"] = cogent_core
        sys.modules["cogent3.core.sequence"] = cogent_seq
        sys.modules["cogent3.core.moltype"] = cogent_moltype
        sys.modules["cogent3.core.alignment"] = cogent_alignment
        sys.modules["cogent3.core.tree"] = cogent_tree

    class _UnavailableOptional:
        def __init__(self, *_args: Any, **_kwargs: Any) -> None:
            missing()

    optional_modules = {
        "gudhi": {
            "SimplexTree": type("SimplexTree", (), {}),
            "RipsComplex": _UnavailableOptional,
            "DelaunayComplex": _UnavailableOptional,
        },
        "faiss": {},
        "piqtree": {"model_finder": types.SimpleNamespace()},
        "ray": {
            "is_initialized": lambda: False,
            "init": missing,
            "remote": lambda obj=None, **_kwargs: obj if obj is not None else (lambda inner: inner),
            "get": lambda value: value,
        },
        "softalign": {},
        "softalign.soft_alignment": {"align_soft_sequences": missing},
    }
    for module_name, attrs in optional_modules.items():
        if module_available(module_name):
            continue
        mod = types.ModuleType(module_name)
        if module_name == "softalign":
            mod.__path__ = []  # mark as a package for softalign.soft_alignment
        for attr, value in attrs.items():
            setattr(mod, attr, value)
        sys.modules[module_name] = mod
    if "softalign" in sys.modules and "softalign.soft_alignment" in sys.modules:
        setattr(sys.modules["softalign"], "soft_alignment", sys.modules["softalign.soft_alignment"])

    if "fitness_landscape.embedding" not in sys.modules:
        embedding_pkg = types.ModuleType("fitness_landscape.embedding")
        embedding_pkg.__path__ = []
        particle_mod = types.ModuleType("fitness_landscape.embedding.particle_sampler")
        beam_mod = types.ModuleType("fitness_landscape.embedding.beam_search")
        soft_mod = types.ModuleType("fitness_landscape.embedding.soft_embedding")

        class _UnusedLandscapyEmbedding:
            def __init__(self, *_args: Any, **_kwargs: Any) -> None:
                missing()

        for name in (
            "EvolutionParticleSampler",
            "SequenceGenerator",
            "ParentSelector",
            "TopPSampler",
            "ESMEmbedder",
        ):
            setattr(particle_mod, name, _UnusedLandscapyEmbedding)
        for name in ("PseudoLogLikelihoodScorer", "InterpolationBeamSearch"):
            setattr(beam_mod, name, _UnusedLandscapyEmbedding)
        soft_mod.ESMEmbedder = _UnusedLandscapyEmbedding
        embedding_pkg.particle_sampler = particle_mod
        embedding_pkg.beam_search = beam_mod
        embedding_pkg.soft_embedding = soft_mod
        sys.modules["fitness_landscape.embedding"] = embedding_pkg
        sys.modules["fitness_landscape.embedding.particle_sampler"] = particle_mod
        sys.modules["fitness_landscape.embedding.beam_search"] = beam_mod
        sys.modules["fitness_landscape.embedding.soft_embedding"] = soft_mod


_install_landscapy_optional_dependency_stubs()

from adaptenv import AlignmentProfile, FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from adaptenv.amino_acids import AA_ORDER  # noqa: E402
from fitness_landscape.analysis.dirichlet_energy import (  # noqa: E402
    calculate_ruggedness_dirichlet_energy,
)
from fitness_landscape.analysis.epistasis import calculate_epistasis_walsh  # noqa: E402
from fitness_landscape.analysis.random_walk import (  # noqa: E402
    calculate_ruggedness_autocorrelation_stochastic,
)
from fitness_landscape.core.fitness import NumericFitness  # noqa: E402
from fitness_landscape.core.landscape import FitnessLandscape  # noqa: E402
from fitness_landscape.core.sequence import BinarySequence  # noqa: E402

from adaptive_walk_peak_epistasis_baselines import (  # noqa: E402
    json_default,
    merge_dict,
    write_csv,
    write_json,
)
from tem1_homologue_walk_parameter_grid import (  # noqa: E402
    TEM1_SEQUENCE,
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


def base_landscape_config(
    cfg: dict[str, Any],
    *,
    seed: int,
    empirical_pairwise_strength: float | None,
) -> LandscapeConfig:
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
            "empirical_pairwise_strength": empirical_pairwise_strength,
        }
    )
    return LandscapeConfig(**base)


def build_adaptenv_landscape(
    cfg: dict[str, Any],
    *,
    seed: int,
    profile: AlignmentProfile,
    pairwise_override: dict[str, np.ndarray] | None,
    empirical_pairwise_strength: float | None,
) -> FitnessLandscapeEnv:
    landscape_config = base_landscape_config(
        cfg,
        seed=seed,
        empirical_pairwise_strength=empirical_pairwise_strength,
    )
    return FitnessLandscapeEnv(
        landscape_config,
        alignment_profile=profile,
        reference_sequence=str(cfg["reference"]["sequence"]),
        pairwise_epistasis_override=pairwise_override,
        _empirical_pairwise_target="stability",
    )


def weighted_joint_frequencies(
    profile: AlignmentProfile,
    pos_i: int,
    pos_j: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    observed = np.asarray(profile.seq_matrix[:, [pos_i, pos_j]], dtype=np.int64)
    valid = np.all(observed >= 0, axis=1)
    if not np.any(valid):
        zeros = np.zeros((20, 20), dtype=np.float64)
        return zeros, zeros.sum(axis=1), zeros.sum(axis=0), 0.0

    observed = observed[valid]
    weights = np.asarray(getattr(profile, "weights", np.ones(profile.N)), dtype=np.float64)[valid]
    joint = np.zeros((20, 20), dtype=np.float64)
    np.add.at(joint, (observed[:, 0], observed[:, 1]), weights)
    total = float(np.sum(joint))
    if total <= 0:
        zeros = np.zeros((20, 20), dtype=np.float64)
        return zeros, zeros.sum(axis=1), zeros.sum(axis=0), 0.0
    joint /= total
    return joint, joint.sum(axis=1), joint.sum(axis=0), total


def evaluate_fitness(landscape: FitnessLandscapeEnv, seq: np.ndarray) -> float:
    return float(landscape.evaluate(np.asarray(seq, dtype=np.int64)))


def select_anchor_sequence(
    cfg: dict[str, Any],
    *,
    profile: AlignmentProfile,
    baseline_landscape: FitnessLandscapeEnv,
) -> tuple[np.ndarray, dict[str, Any], list[dict[str, Any]]]:
    anchor_cfg = dict(cfg.get("anchor_selection", {}))
    mode = str(anchor_cfg.get("mode", "closest_realized_margin"))
    target_margin = float(anchor_cfg.get("target_effective_stability_margin", 1.0))
    max_candidates = int(anchor_cfg.get("max_reported_candidates", 5000))

    seqs, source_indices = filled_unique_homologues(profile)
    reference = np.asarray(baseline_landscape.reference, dtype=np.int64)
    consensus = np.asarray(profile.consensus_idx, dtype=np.int64)
    log_odds = profile_log_odds_bits(profile)
    positions = np.arange(profile.L, dtype=np.int64)
    consensus_score = float(np.sum(log_odds[positions, consensus]))

    candidates: list[dict[str, Any]] = []
    candidate_seqs: list[np.ndarray] = []

    def add_candidate(
        seq: np.ndarray,
        *,
        source: str,
        source_alignment_index: int | None,
        unique_homologue_index: int | None,
    ) -> None:
        components = components_for_seq(baseline_landscape, seq)
        row = {
            "anchor_source": source,
            "source_alignment_index": source_alignment_index,
            "unique_homologue_index": unique_homologue_index,
            "sequence": seq_to_string(seq),
            "fitness": evaluate_fitness(baseline_landscape, seq),
            **components,
            **consensus_metrics(
                seq,
                consensus=consensus,
                log_odds_bits=log_odds,
                consensus_score=consensus_score,
            ),
        }
        candidates.append(row)
        candidate_seqs.append(np.asarray(seq, dtype=np.int64).copy())

    if bool(anchor_cfg.get("include_reference_candidate", True)):
        add_candidate(
            reference,
            source="tem1_reference",
            source_alignment_index=None,
            unique_homologue_index=None,
        )

    for unique_idx, seq in enumerate(seqs):
        add_candidate(
            seq,
            source="mmseqs_homologue",
            source_alignment_index=int(source_indices[unique_idx]),
            unique_homologue_index=int(unique_idx),
        )

    if mode == "reference":
        chosen_idx = 0
    elif mode == "min_realized_margin":
        margins = np.asarray(
            [row["effective_stability_margin_used"] for row in candidates],
            dtype=np.float64,
        )
        chosen_idx = int(np.argmin(margins))
    elif mode == "closest_realized_margin":
        margins = np.asarray(
            [row["effective_stability_margin_used"] for row in candidates],
            dtype=np.float64,
        )
        chosen_idx = int(np.argmin(np.abs(margins - target_margin)))
    else:
        raise ValueError(
            "anchor_selection.mode must be one of reference, min_realized_margin, "
            "or closest_realized_margin"
        )

    for idx, row in enumerate(candidates):
        row["selected_anchor"] = bool(idx == chosen_idx)
        row["anchor_selection_mode"] = mode
        row["target_effective_stability_margin"] = target_margin

    candidates.sort(
        key=lambda row: (
            not bool(row["selected_anchor"]),
            abs(float(row["effective_stability_margin_used"]) - target_margin),
        )
    )
    return (
        candidate_seqs[chosen_idx],
        dict(candidates[0]),
        candidates[:max_candidates],
    )


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
    mutant[position] = int(alt_aa)
    delta = evaluate_fitness(landscape, mutant) - parent_fitness
    cache[key] = float(delta)
    return float(delta)


def candidate_state_pairs_for_coupling(
    *,
    profile: AlignmentProfile,
    anchor_seq: np.ndarray,
    reference_seq: np.ndarray,
    baseline_landscape: FitnessLandscapeEnv,
    parent_fitness: float,
    pos_i: int,
    pos_j: int,
    coupling_score: float,
    pair_rank: int,
    dca_cfg: dict[str, Any],
    delta_cache: dict[tuple[int, int], float],
) -> list[dict[str, Any]]:
    joint, fi, fj, effective_pair_weight = weighted_joint_frequencies(profile, pos_i, pos_j)
    min_joint_frequency = float(dca_cfg.get("min_joint_frequency", 0.001))
    min_log2_enrichment = float(dca_cfg.get("min_pair_log2_enrichment", 0.0))
    top_states = int(dca_cfg.get("candidate_state_pairs_per_coupling", 12))
    target_single_delta = float(dca_cfg.get("target_single_delta_for_compensation", -25.0))
    single_delta_tolerance = max(float(dca_cfg.get("single_delta_tolerance", 50.0)), 1e-9)
    anchor_i = int(anchor_seq[pos_i])
    anchor_j = int(anchor_seq[pos_j])
    reference_i = int(reference_seq[pos_i])
    reference_j = int(reference_seq[pos_j])
    eps = 1e-12

    rows: list[dict[str, Any]] = []
    for aa_i in range(20):
        if aa_i == anchor_i or aa_i == reference_i:
            continue
        for aa_j in range(20):
            if aa_j == anchor_j or aa_j == reference_j:
                continue
            observed = float(joint[aa_i, aa_j])
            if observed <= 0:
                continue
            expected = float(fi[aa_i] * fj[aa_j])
            log2_enrichment = float(np.log2((observed + eps) / (expected + eps)))
            delta_i = single_mutant_delta(
                baseline_landscape,
                anchor_seq,
                parent_fitness,
                position=pos_i,
                alt_aa=aa_i,
                cache=delta_cache,
            )
            delta_j = single_mutant_delta(
                baseline_landscape,
                anchor_seq,
                parent_fitness,
                position=pos_j,
                alt_aa=aa_j,
                cache=delta_cache,
            )
            passes_frequency = observed >= min_joint_frequency
            passes_enrichment = log2_enrichment >= min_log2_enrichment
            deleterious_bonus = max(0.0, -delta_i) + max(0.0, -delta_j)
            target_closeness = math.exp(
                -(
                    abs(delta_i - target_single_delta)
                    + abs(delta_j - target_single_delta)
                )
                / single_delta_tolerance
            )
            score = (
                float(coupling_score)
                * max(observed, min_joint_frequency)
                * max(log2_enrichment + 1.0, 0.05)
                * (1.0 + deleterious_bonus / max(abs(parent_fitness), 1.0))
                * (0.25 + 0.75 * target_closeness)
            )
            rows.append(
                {
                    "pair_rank": int(pair_rank),
                    "pos_i": int(pos_i),
                    "pos_j": int(pos_j),
                    "coupling_score": float(coupling_score),
                    "effective_pair_weight": float(effective_pair_weight),
                    "anchor_aa_i": anchor_i,
                    "anchor_aa_j": anchor_j,
                    "reference_aa_i": reference_i,
                    "reference_aa_j": reference_j,
                    "alt_aa_i": int(aa_i),
                    "alt_aa_j": int(aa_j),
                    "joint_frequency": observed,
                    "expected_frequency": expected,
                    "log2_enrichment": log2_enrichment,
                    "single_delta_i": delta_i,
                    "single_delta_j": delta_j,
                    "target_single_delta": target_single_delta,
                    "single_delta_target_closeness": float(target_closeness),
                    "passes_frequency": passes_frequency,
                    "passes_enrichment": passes_enrichment,
                    "selection_score": float(score),
                }
            )

    rows.sort(
        key=lambda row: (
            not (bool(row["passes_frequency"]) and bool(row["passes_enrichment"])),
            -float(row["selection_score"]),
        )
    )
    return rows[:top_states]


def select_dca_specific_pairs(
    cfg: dict[str, Any],
    *,
    profile: AlignmentProfile,
    baseline_landscape: FitnessLandscapeEnv,
    anchor_seq: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, np.ndarray]]:
    dca_cfg = dict(cfg["dca_specific_epistasis"])
    target_mutation_count = int(dca_cfg.get("mutation_count", 10))
    if target_mutation_count < 2:
        raise ValueError("dca_specific_epistasis.mutation_count must be at least 2")
    target_pair_count = int(dca_cfg.get("pair_count", target_mutation_count // 2))
    scan_top_couplings = int(dca_cfg.get("scan_top_couplings", 1000))
    require_disjoint = bool(dca_cfg.get("require_disjoint_pairs", True))
    max_single_delta = float(dca_cfg.get("max_single_delta_for_compensation", math.inf))
    min_selected_pairs = int(dca_cfg.get("min_selected_pairs", target_pair_count))

    parent_fitness = evaluate_fitness(baseline_landscape, anchor_seq)
    reference_seq = np.asarray(baseline_landscape.reference, dtype=np.int64)
    delta_cache: dict[tuple[int, int], float] = {}
    candidate_rows: list[dict[str, Any]] = []
    top_couplings = profile.get_top_couplings(n_pairs=scan_top_couplings)

    for pair_rank, (pos_i, pos_j, coupling_score) in enumerate(top_couplings):
        if float(coupling_score) <= 0:
            continue
        candidate_rows.extend(
            candidate_state_pairs_for_coupling(
                profile=profile,
                anchor_seq=anchor_seq,
                reference_seq=reference_seq,
                baseline_landscape=baseline_landscape,
                parent_fitness=parent_fitness,
                pos_i=int(pos_i),
                pos_j=int(pos_j),
                coupling_score=float(coupling_score),
                pair_rank=int(pair_rank),
                dca_cfg=dca_cfg,
                delta_cache=delta_cache,
            )
        )

    def tier(row: dict[str, Any]) -> int:
        singles_ok = (
            float(row["single_delta_i"]) <= max_single_delta
            and float(row["single_delta_j"]) <= max_single_delta
        )
        dca_ok = bool(row["passes_frequency"]) and bool(row["passes_enrichment"])
        if dca_ok and singles_ok:
            return 0
        if dca_ok:
            return 1
        if singles_ok:
            return 2
        return 3

    candidate_rows.sort(key=lambda row: (tier(row), -float(row["selection_score"])))

    selected: list[dict[str, Any]] = []
    used_positions: set[int] = set()
    seen_pairs: set[tuple[int, int]] = set()
    for row in candidate_rows:
        pair = (int(row["pos_i"]), int(row["pos_j"]))
        if pair in seen_pairs:
            continue
        if require_disjoint and (pair[0] in used_positions or pair[1] in used_positions):
            continue
        selected.append(dict(row))
        seen_pairs.add(pair)
        used_positions.update(pair)
        if len(selected) >= target_pair_count:
            break

    if len(selected) < min_selected_pairs:
        raise RuntimeError(
            f"Only selected {len(selected)} DCA-specific pairs; "
            f"minimum required is {min_selected_pairs}. Relax the DCA filters."
        )

    positions = np.asarray([[row["pos_i"], row["pos_j"]] for row in selected], dtype=np.int64)
    matrices = np.zeros((len(selected), 20, 20), dtype=np.float64)
    raw_sign = float(dca_cfg.get("raw_matrix_entry", 1.0))
    for idx, row in enumerate(selected):
        matrices[idx, int(row["alt_aa_i"]), int(row["alt_aa_j"])] = raw_sign
        row["selected_pair_index"] = int(idx)
        row["selection_tier"] = int(tier(row))
        row["mutation_i"] = mutation_name(
            int(row["pos_i"]),
            int(row["anchor_aa_i"]),
            int(row["alt_aa_i"]),
        )
        row["mutation_j"] = mutation_name(
            int(row["pos_j"]),
            int(row["anchor_aa_j"]),
            int(row["alt_aa_j"]),
        )
        row["term_mode"] = "dca_compensatory_specific"

    return selected, {"positions": positions, "matrices": matrices}


def mutation_variables(selected_pairs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    variables: list[dict[str, Any]] = []
    seen: set[int] = set()
    for pair_index, row in enumerate(selected_pairs):
        for side in ("i", "j"):
            pos = int(row[f"pos_{side}"])
            if pos in seen:
                continue
            seen.add(pos)
            anchor_aa = int(row[f"anchor_aa_{side}"])
            alt_aa = int(row[f"alt_aa_{side}"])
            variables.append(
                {
                    "variable_index": int(len(variables)),
                    "pair_index": int(pair_index),
                    "position": pos,
                    "anchor_aa": anchor_aa,
                    "alt_aa": alt_aa,
                    "mutation": mutation_name(pos, anchor_aa, alt_aa),
                }
            )
    return variables


def sequence_from_bits(
    anchor_seq: np.ndarray,
    variables: list[dict[str, Any]],
    bits: np.ndarray,
) -> np.ndarray:
    seq = np.asarray(anchor_seq, dtype=np.int64).copy()
    for variable, bit in zip(variables, bits):
        if int(bit):
            seq[int(variable["position"])] = int(variable["alt_aa"])
    return seq


def build_dense_hamming_landscape(
    *,
    anchor_seq: np.ndarray,
    variables: list[dict[str, Any]],
    dca_landscape: FitnessLandscapeEnv,
    additive_landscape: FitnessLandscapeEnv,
) -> tuple[FitnessLandscape, list[dict[str, Any]], np.ndarray, np.ndarray]:
    k = len(variables)
    n_nodes = 1 << k
    bits_by_mask = np.zeros((n_nodes, k), dtype=np.int8)
    fitness = np.zeros(n_nodes, dtype=np.float64)
    additive_fitness = np.zeros(n_nodes, dtype=np.float64)
    node_rows: list[dict[str, Any]] = []
    sequences: list[BinarySequence] = []

    for mask in range(n_nodes):
        bits = np.asarray([(mask >> idx) & 1 for idx in range(k)], dtype=np.int8)
        bits_by_mask[mask] = bits
        seq = sequence_from_bits(anchor_seq, variables, bits)
        components = components_for_seq(dca_landscape, seq)
        additive_components = components_for_seq(additive_landscape, seq)
        fitness[mask] = float(components["fitness"] if "fitness" in components else dca_landscape.evaluate(seq))
        additive_fitness[mask] = evaluate_fitness(additive_landscape, seq)
        sequences.append(BinarySequence.from_bits(bits))
        mutations = [
            str(variable["mutation"])
            for variable, bit in zip(variables, bits)
            if int(bit)
        ]
        node_rows.append(
            {
                "node_index": int(mask),
                "bitstring_lsb_first": "".join(str(int(x)) for x in bits.tolist()),
                "hamming_weight": int(np.sum(bits)),
                "amino_acid_sequence": seq_to_string(seq),
                "mutations": ";".join(mutations),
                "fitness": float(fitness[mask]),
                "additive_fitness": float(additive_fitness[mask]),
                "dca_specific_fitness_delta": float(fitness[mask] - additive_fitness[mask]),
                **components,
            }
        )

    layer = NumericFitness(name="default", values=fitness)
    landscape = FitnessLandscape.build(
        sequences,
        fitness_layers={"default": layer},
        graph="hamming",
    )
    return landscape, node_rows, bits_by_mask, fitness


def annotate_local_optima(
    landscape: FitnessLandscape,
    node_rows: list[dict[str, Any]],
    fitness: np.ndarray,
) -> dict[str, Any]:
    local_maxima: list[int] = []
    local_minima: list[int] = []
    for node in landscape.graph.nodes():
        neighbors = list(landscape.graph.neighbors(node))
        neighbor_values = fitness[np.asarray(neighbors, dtype=np.int64)]
        uphill_degree = int(np.sum(neighbor_values > fitness[node]))
        downhill_degree = int(np.sum(neighbor_values < fitness[node]))
        node_rows[int(node)]["uphill_degree"] = uphill_degree
        node_rows[int(node)]["downhill_degree"] = downhill_degree
        node_rows[int(node)]["is_local_maximum"] = uphill_degree == 0
        node_rows[int(node)]["is_local_minimum"] = downhill_degree == 0
        if uphill_degree == 0:
            local_maxima.append(int(node))
        if downhill_degree == 0:
            local_minima.append(int(node))
    return {
        "n_local_maxima": int(len(local_maxima)),
        "n_local_minima": int(len(local_minima)),
        "local_maxima": local_maxima,
        "local_minima": local_minima,
        "global_maximum_node": int(np.argmax(fitness)),
        "global_minimum_node": int(np.argmin(fitness)),
    }


def accessible_shortest_paths(fitness: np.ndarray, *, k: int) -> dict[str, Any]:
    target = (1 << k) - 1
    paths = np.zeros(1 << k, dtype=object)
    paths[0] = 1
    for mask in range(1, 1 << k):
        total = 0
        for bit in range(k):
            if not ((mask >> bit) & 1):
                continue
            prev = mask ^ (1 << bit)
            if fitness[prev] < fitness[mask]:
                total += paths[prev]
        paths[mask] = total
    total_shortest = math.factorial(k)
    return {
        "target_node": int(target),
        "accessible_shortest_paths_to_all_ones": int(paths[target]),
        "total_shortest_paths_to_all_ones": int(total_shortest),
        "accessible_shortest_path_fraction_to_all_ones": float(paths[target] / total_shortest),
    }


def greedy_walk_summary(
    landscape: FitnessLandscape,
    fitness: np.ndarray,
    *,
    start_node: int = 0,
) -> dict[str, Any]:
    current = int(start_node)
    path = [current]
    gains: list[float] = []
    while True:
        neighbors = list(landscape.graph.neighbors(current))
        better = [node for node in neighbors if fitness[node] > fitness[current]]
        if not better:
            break
        next_node = int(max(better, key=lambda node: fitness[node]))
        gains.append(float(fitness[next_node] - fitness[current]))
        current = next_node
        path.append(current)
        if len(path) > landscape.graph.number_of_nodes():
            raise RuntimeError("Greedy walk exceeded graph size; graph may contain a cycle bug")
    return {
        "greedy_start_node": int(start_node),
        "greedy_final_node": int(current),
        "greedy_steps": int(len(path) - 1),
        "greedy_total_gain": float(fitness[current] - fitness[start_node]),
        "greedy_path": path,
        "greedy_step_gains": gains,
    }


def sign_flags(effect_i: float, effect_j: float, effect_i_on_j: float, effect_j_on_i: float) -> dict[str, bool]:
    sign_i = np.sign(effect_i)
    sign_j = np.sign(effect_j)
    sign_i_on_j = np.sign(effect_i_on_j)
    sign_j_on_i = np.sign(effect_j_on_i)
    sign_i_changed = bool(sign_i != 0 and sign_i_on_j != 0 and sign_i != sign_i_on_j)
    sign_j_changed = bool(sign_j != 0 and sign_j_on_i != 0 and sign_j != sign_j_on_i)
    return {
        "sign_epistasis": bool(sign_i_changed or sign_j_changed),
        "reciprocal_sign_epistasis": bool(sign_i_changed and sign_j_changed),
    }


def pair_epistasis_rows(
    *,
    selected_pairs: list[dict[str, Any]],
    variables: list[dict[str, Any]],
    fitness: np.ndarray,
    k: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    var_by_position = {int(variable["position"]): int(variable["variable_index"]) for variable in variables}
    rows: list[dict[str, Any]] = []
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for pair_index, pair in enumerate(selected_pairs):
        var_i = var_by_position[int(pair["pos_i"])]
        var_j = var_by_position[int(pair["pos_j"])]
        for background in range(1 << k):
            if ((background >> var_i) & 1) or ((background >> var_j) & 1):
                continue
            m00 = background
            m10 = background | (1 << var_i)
            m01 = background | (1 << var_j)
            m11 = background | (1 << var_i) | (1 << var_j)
            f00 = float(fitness[m00])
            f10 = float(fitness[m10])
            f01 = float(fitness[m01])
            f11 = float(fitness[m11])
            effect_i = f10 - f00
            effect_j = f01 - f00
            effect_i_on_j = f11 - f01
            effect_j_on_i = f11 - f10
            epistasis = f11 - f10 - f01 + f00
            flags = sign_flags(effect_i, effect_j, effect_i_on_j, effect_j_on_i)
            row = {
                "selected_pair_index": int(pair_index),
                "background_node": int(background),
                "var_i": int(var_i),
                "var_j": int(var_j),
                "pos_i": int(pair["pos_i"]),
                "pos_j": int(pair["pos_j"]),
                "mutation_i": str(pair["mutation_i"]),
                "mutation_j": str(pair["mutation_j"]),
                "coupling_score": float(pair["coupling_score"]),
                "joint_frequency": float(pair["joint_frequency"]),
                "log2_enrichment": float(pair["log2_enrichment"]),
                "f00": f00,
                "f10": f10,
                "f01": f01,
                "f11": f11,
                "effect_i": effect_i,
                "effect_j": effect_j,
                "effect_i_on_j": effect_i_on_j,
                "effect_j_on_i": effect_j_on_i,
                "epistasis": float(epistasis),
                "abs_epistasis": float(abs(epistasis)),
                **flags,
            }
            rows.append(row)
            grouped[pair_index].append(row)

    summary_rows = []
    for pair_index, pair_rows in grouped.items():
        epis = np.asarray([row["epistasis"] for row in pair_rows], dtype=np.float64)
        signs = np.asarray([row["sign_epistasis"] for row in pair_rows], dtype=bool)
        recip = np.asarray([row["reciprocal_sign_epistasis"] for row in pair_rows], dtype=bool)
        summary_rows.append(
            {
                "selected_pair_index": int(pair_index),
                "mean_epistasis": float(np.mean(epis)),
                "mean_abs_epistasis": float(np.mean(np.abs(epis))),
                "min_epistasis": float(np.min(epis)),
                "max_epistasis": float(np.max(epis)),
                "sign_epistasis_fraction": float(np.mean(signs)),
                "reciprocal_sign_epistasis_fraction": float(np.mean(recip)),
                "n_backgrounds": int(len(pair_rows)),
            }
        )
    return rows, {"pair_summaries": summary_rows}


def walsh_rows(landscape: FitnessLandscape, *, max_order: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    result = calculate_epistasis_walsh(landscape, order=max_order)
    rows: list[dict[str, Any]] = []
    variance_by_order = {
        int(order): float(value)
        for order, value in result.get("variance_explained", {}).items()
    }
    for order, terms in sorted(result.get("by_order", {}).items()):
        for term, value in terms.items():
            rows.append(
                {
                    "order": int(order),
                    "term": str(term),
                    "coefficient": float(value),
                    "abs_coefficient": float(abs(value)),
                }
            )
    rows.sort(key=lambda row: (int(row["order"]), -float(row["abs_coefficient"])))
    return rows, {
        "walsh_variance_explained": variance_by_order,
        "walsh_statistics": result.get("statistics", {}),
    }


def run_experiment(cfg: dict[str, Any], *, config_path: Path, project_root: Path, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    seed = int(cfg.get("seed", 9401))
    profile = build_alignment_profile(
        cfg,
        project_root=project_root,
        config_path=config_path,
        output_dir=output_dir,
    )

    additive_landscape = build_adaptenv_landscape(
        cfg,
        seed=seed,
        profile=profile,
        pairwise_override=None,
        empirical_pairwise_strength=0.0,
    )
    anchor_seq, anchor_row, anchor_candidates = select_anchor_sequence(
        cfg,
        profile=profile,
        baseline_landscape=additive_landscape,
    )
    selected_pairs, pairwise_override = select_dca_specific_pairs(
        cfg,
        profile=profile,
        baseline_landscape=additive_landscape,
        anchor_seq=anchor_seq,
    )
    variables = mutation_variables(selected_pairs)
    pairwise_strength = float(cfg["dca_specific_epistasis"]["pairwise_strength"])
    dca_landscape = build_adaptenv_landscape(
        cfg,
        seed=seed,
        profile=profile,
        pairwise_override=pairwise_override,
        empirical_pairwise_strength=pairwise_strength,
    )

    hamming_landscape, node_rows, bits_by_mask, fitness = build_dense_hamming_landscape(
        anchor_seq=anchor_seq,
        variables=variables,
        dca_landscape=dca_landscape,
        additive_landscape=additive_landscape,
    )
    optima_summary = annotate_local_optima(hamming_landscape, node_rows, fitness)
    accessibility = accessible_shortest_paths(fitness, k=len(variables))
    greedy = greedy_walk_summary(hamming_landscape, fitness, start_node=0)
    pair_rows, pair_summary = pair_epistasis_rows(
        selected_pairs=selected_pairs,
        variables=variables,
        fitness=fitness,
        k=len(variables),
    )
    walsh_order = min(int(cfg.get("analysis", {}).get("walsh_max_order", 4)), len(variables))
    walsh_table, walsh_summary = walsh_rows(hamming_landscape, max_order=walsh_order)
    random_cfg = dict(cfg.get("analysis", {}).get("random_walk", {}))
    random_walk = calculate_ruggedness_autocorrelation_stochastic(
        hamming_landscape,
        n_walks=int(random_cfg.get("n_walks", 200)),
        steps=int(random_cfg.get("steps", 200)),
        lag_max=int(random_cfg.get("lag_max", 8)),
        seed=seed,
    )
    dirichlet = calculate_ruggedness_dirichlet_energy(hamming_landscape)

    bundle_dir = output_dir / "landscapy_hamming_landscape"
    hamming_landscape.save_bundle_dir(
        bundle_dir,
        metadata={
            "dataset_name": "tem1-dca-specific-epistasis-hamming-graph",
            "source_name": "project-adapt-env EXP-023",
            "protein_gene": "TEM-1 beta-lactamase",
            "assay_type": "synthetic stability-margin local Hamming graph",
            "version": "RUN pending",
            "tags": ["TEM-1", "DCA", "specific-epistasis", "hamming-cube"],
        },
        include_embeddings=False,
    )

    variable_rows = [
        {
            **variable,
            "anchor_aa_name": aa_name(int(variable["anchor_aa"])),
            "alt_aa_name": aa_name(int(variable["alt_aa"])),
        }
        for variable in variables
    ]
    for row in selected_pairs:
        row["anchor_aa_i_name"] = aa_name(int(row["anchor_aa_i"]))
        row["anchor_aa_j_name"] = aa_name(int(row["anchor_aa_j"]))
        row["alt_aa_i_name"] = aa_name(int(row["alt_aa_i"]))
        row["alt_aa_j_name"] = aa_name(int(row["alt_aa_j"]))

    np.savez_compressed(
        output_dir / "dca_pairwise_override.npz",
        positions=pairwise_override["positions"],
        matrices=pairwise_override["matrices"],
        variables=np.asarray(
            [
                [
                    int(variable["variable_index"]),
                    int(variable["position"]),
                    int(variable["anchor_aa"]),
                    int(variable["alt_aa"]),
                ]
                for variable in variables
            ],
            dtype=np.int64,
        ),
        bits_by_mask=bits_by_mask,
    )

    write_csv(
        output_dir / "anchor_candidates.csv",
        anchor_candidates,
        [
            "selected_anchor",
            "anchor_selection_mode",
            "target_effective_stability_margin",
            "anchor_source",
            "source_alignment_index",
            "unique_homologue_index",
            "sequence",
            "fitness",
            "stability",
            "effective_stability",
            "stability_margin_used",
            "effective_stability_margin_used",
            "stability_gate",
            "epistatic_score",
            "generic_pairwise_score",
            "empirical_pairwise_score",
            "latent_readout",
            "consensus_bitscore",
            "consensus_distance_bits",
            "consensus_hamming_distance",
            "consensus_identity",
        ],
    )
    write_csv(
        output_dir / "selected_dca_pairs.csv",
        selected_pairs,
        [
            "selected_pair_index",
            "selection_tier",
            "term_mode",
            "pair_rank",
            "pos_i",
            "pos_j",
            "mutation_i",
            "mutation_j",
            "coupling_score",
            "effective_pair_weight",
            "joint_frequency",
            "expected_frequency",
            "log2_enrichment",
            "single_delta_i",
            "single_delta_j",
            "target_single_delta",
            "single_delta_target_closeness",
            "anchor_aa_i_name",
            "anchor_aa_j_name",
            "alt_aa_i_name",
            "alt_aa_j_name",
            "reference_aa_i",
            "reference_aa_j",
            "selection_score",
        ],
    )
    write_csv(
        output_dir / "hamming_variables.csv",
        variable_rows,
        [
            "variable_index",
            "pair_index",
            "position",
            "anchor_aa",
            "anchor_aa_name",
            "alt_aa",
            "alt_aa_name",
            "mutation",
        ],
    )
    write_csv(
        output_dir / "hamming_nodes.csv",
        node_rows,
        [
            "node_index",
            "bitstring_lsb_first",
            "hamming_weight",
            "mutations",
            "fitness",
            "additive_fitness",
            "dca_specific_fitness_delta",
            "stability",
            "effective_stability",
            "stability_margin_used",
            "effective_stability_margin_used",
            "stability_gate",
            "epistatic_score",
            "generic_pairwise_score",
            "empirical_pairwise_score",
            "latent_readout",
            "uphill_degree",
            "downhill_degree",
            "is_local_maximum",
            "is_local_minimum",
            "amino_acid_sequence",
        ],
    )
    write_csv(
        output_dir / "hamming_pair_epistasis.csv",
        pair_rows,
        [
            "selected_pair_index",
            "background_node",
            "var_i",
            "var_j",
            "pos_i",
            "pos_j",
            "mutation_i",
            "mutation_j",
            "coupling_score",
            "joint_frequency",
            "log2_enrichment",
            "f00",
            "f10",
            "f01",
            "f11",
            "effect_i",
            "effect_j",
            "effect_i_on_j",
            "effect_j_on_i",
            "epistasis",
            "abs_epistasis",
            "sign_epistasis",
            "reciprocal_sign_epistasis",
        ],
    )
    write_csv(
        output_dir / "hamming_walsh_terms.csv",
        walsh_table,
        ["order", "term", "coefficient", "abs_coefficient"],
    )

    random_autocorr = np.asarray(random_walk["autocorrelation"], dtype=np.float64)
    graph_metrics = {
        "n_variables": int(len(variables)),
        "n_nodes": int(hamming_landscape.graph.number_of_nodes()),
        "n_edges": int(hamming_landscape.graph.number_of_edges()),
        "fitness_min": float(np.min(fitness)),
        "fitness_max": float(np.max(fitness)),
        "fitness_mean": float(np.mean(fitness)),
        "fitness_sd": float(np.std(fitness)),
        "random_walk_lag1_autocorrelation": (
            float(random_autocorr[1]) if random_autocorr.size > 1 else None
        ),
        "random_walk_correlation_length": float(random_walk["correlation_length"]),
        "dirichlet_energy_per_node": float(dirichlet["total_dirichlet_energy"]),
        **optima_summary,
        **accessibility,
        **greedy,
        **pair_summary,
        **walsh_summary,
    }
    write_json(output_dir / "hamming_graph_metrics.json", graph_metrics)

    summary = {
        "experiment_id": str(cfg.get("experiment_id", "EXP-023")),
        "reference": cfg["reference"],
        "seed": seed,
        "alignment": {
            "n_sequences": int(profile.N),
            "n_effective_sequences": float(profile.N_eff),
            "length": int(profile.L),
            "consensus_sequence": str(profile.consensus_seq),
        },
        "old_epistasis_disabled": {
            "n_epistatic_terms": 0,
            "epistasis_strength": 0.0,
            "epistasis_covariation_bias": 0.0,
            "n_higher_order_epistatic_terms": 0,
        },
        "dca_specific_epistasis": {
            **cfg["dca_specific_epistasis"],
            "selected_pair_count": int(len(selected_pairs)),
            "selected_mutation_count": int(len(variables)),
            "pairwise_strength": pairwise_strength,
        },
        "anchor": anchor_row,
        "variables": variable_rows,
        "graph_metrics": graph_metrics,
        "outputs": {
            "landscapy_bundle": str(bundle_dir),
            "node_table": str(output_dir / "hamming_nodes.csv"),
            "pair_epistasis_table": str(output_dir / "hamming_pair_epistasis.csv"),
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
    print(json.dumps(summary["graph_metrics"], indent=2, sort_keys=True, default=json_default))


if __name__ == "__main__":
    main()
