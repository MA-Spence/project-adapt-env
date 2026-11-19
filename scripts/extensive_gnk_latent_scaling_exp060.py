#!/usr/bin/env python3
"""EXP-060 paired latent scale-law test for an extensive GNK candidate.

The experiment ends at additive stability plus the centered GNK field.  No
GNK score is routed through AdaptEnv's sigmoid and no observable-fitness,
accessibility, or optimizer metric is calculated.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPTENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
for _path in (SRC_ROOT, ADAPTENV_ROOT):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from adaptenv import FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from adaptenv.amino_acids import AA_TO_IDX  # noqa: E402
from adaptenv.landscape import _stable_standard_normal as adaptenv_standard_normal  # noqa: E402
from project_adapt_env.extensive_gnk import (  # noqa: E402
    ExtensiveGNKField,
    additive_mutation_energy,
    stable_uint64,
    uniform_alternative_substitution_energy,
)


ARM_LEGACY = "legacy_public_global_fixed"
ARM_LOCAL = "legacy_support_local_normalized"
ARM_EXTENSIVE = "extensive_balanced_bounded"
ARMS = (ARM_LEGACY, ARM_LOCAL, ARM_EXTENSIVE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def deep_merge(base: dict[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = deep_merge(dict(out[key]), value)
        else:
            out[key] = value
    return out


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("config root must be a mapping")
    if quick:
        override = payload.get("quick", {})
        if not isinstance(override, Mapping):
            raise ValueError("quick config must be a mapping")
        payload = deep_merge(payload, override)
    payload["_runtime"] = {"quick": bool(quick)}
    return payload


def json_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def array_hash(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(json.dumps(array.shape).encode("ascii"))
    digest.update(array.view(np.uint8))
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv_atomic(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    *,
    columns: Sequence[str] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame = pd.DataFrame(list(rows), columns=list(columns) if columns else None)
    if frame.empty and not len(frame.columns):
        frame = pd.DataFrame(columns=["status"])
    elif columns is None:
        frame = frame.reindex(sorted(frame.columns), axis=1)
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def safe_name(value: str) -> str:
    return "".join(character if character.isalnum() or character in "._-=" else "_" for character in str(value))


def resolve_processed_root(project_root: Path, config: Mapping[str, Any]) -> Path:
    configured = Path(str(config.get("data", {}).get("processed_root", "data/processed/ensemble-landscape-validation-datasets")))
    if not configured.is_absolute():
        configured = project_root / configured
    override = os.environ.get("STAGE4_PROCESSED_ROOT")
    if override:
        configured = Path(override)
    return configured


def _valid_protein_sequence(value: Any) -> bool:
    sequence = str(value or "").strip().upper()
    return len(sequence) > 0 and all(character in AA_TO_IDX for character in sequence)


def load_reference_catalog(project_root: Path, config: Mapping[str, Any]) -> list[dict[str, Any]]:
    root = resolve_processed_root(project_root, config)
    rows: list[dict[str, Any]] = []
    reference_config = config.get("references", {})
    synthetic = reference_config.get("synthetic", [])
    for item in synthetic:
        reference_id = str(item["id"])
        length = int(item["length"])
        seed = int(item.get("seed", stable_uint64("synthetic-reference", reference_id, length)))
        rng = np.random.default_rng(seed)
        indices = rng.integers(0, 20, size=length, dtype=np.int64)
        sequence = "".join("ACDEFGHIKLMNPQRSTVWY"[int(index)] for index in indices)
        rows.append(
            {
                "reference_id": reference_id,
                "dataset_id": "",
                "source": "synthetic_bridge",
                "source_field": "generated_uniform_q20",
                "length": length,
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            }
        )

    for item in reference_config.get("inline", []):
        reference_id = str(item["id"])
        dataset_id = str(item.get("dataset_id", reference_id))
        sequence = str(item["sequence"]).strip().upper()
        if not _valid_protein_sequence(sequence):
            raise ValueError(f"invalid inline protein reference: {reference_id}")
        expected_length = item.get("expected_length")
        if expected_length is not None and len(sequence) != int(expected_length):
            raise ValueError(
                f"{reference_id} inline reference length {len(sequence)} does not match "
                f"{expected_length}"
            )
        sequence_sha256 = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        declared_sha256 = item.get("sequence_sha256")
        if declared_sha256 is not None and sequence_sha256 != str(declared_sha256):
            raise ValueError(f"inline reference hash mismatch: {reference_id}")
        rows.append(
            {
                "reference_id": reference_id,
                "dataset_id": dataset_id,
                "source": str(item.get("source", "config.yaml:inline")),
                "source_field": str(item.get("source_field", "inline_sequence")),
                "length": len(sequence),
                "sequence": sequence,
                "sequence_sha256": sequence_sha256,
            }
        )

    preferred_fields = tuple(reference_config.get("field_priority", ["target_sequence", "wildtype_sequence", "sequence"]))
    for item in reference_config.get("datasets", []):
        dataset_id = str(item["dataset_id"])
        reference_id = str(item.get("id", dataset_id))
        path = root / "datasets" / dataset_id / "variants.csv"
        if not path.exists():
            raise FileNotFoundError(f"missing declared reference source: {path}")
        with path.open("r", encoding="utf-8", newline="") as handle:
            first = next(csv.DictReader(handle), None)
        if first is None:
            raise ValueError(f"empty reference source: {path}")
        sequence = ""
        source_field = ""
        for field in preferred_fields:
            if field in first and _valid_protein_sequence(first[field]):
                sequence = str(first[field]).strip().upper()
                source_field = str(field)
                break
        if not sequence:
            raise ValueError(f"no valid protein reference field found in {path}")
        expected_length = item.get("expected_length")
        if expected_length is not None and len(sequence) != int(expected_length):
            raise ValueError(
                f"{dataset_id} reference length {len(sequence)} does not match {expected_length}"
            )
        rows.append(
            {
                "reference_id": reference_id,
                "dataset_id": dataset_id,
                "source": str(path.relative_to(project_root)) if path.is_relative_to(project_root) else str(path),
                "source_field": source_field,
                "length": len(sequence),
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            }
        )
    scientific_decision_enabled = not bool(
        config.get("experiment", {}).get("disable_scientific_decision", False)
    )
    minimum_lengths = 3 if scientific_decision_enabled else 2
    if len(rows) < minimum_lengths:
        raise ValueError(
            f"EXP-060 requires at least {minimum_lengths} distinct sequence lengths "
            f"when scientific_decision_enabled={scientific_decision_enabled}"
        )
    lengths = [int(row["length"]) for row in rows]
    if len(lengths) != len(set(lengths)):
        raise ValueError(f"reference lengths must be unique for scale-law fitting: {lengths}")
    return sorted(rows, key=lambda row: int(row["length"]))


def sequence_to_indices(sequence: str) -> np.ndarray:
    return np.asarray([AA_TO_IDX[character] for character in sequence], dtype=np.int64)


def distance_specs(length: int, values: Sequence[Any]) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    for value in values:
        if isinstance(value, str) and value.endswith("L"):
            fraction = float(value[:-1])
            distance = max(1, min(int(round(fraction * length)), length))
            label = f"fraction_{fraction:g}"
        else:
            distance = max(1, min(int(value), length))
            label = f"fixed_{int(value)}"
        out.append((label, distance))
    labels = [label for label, _ in out]
    if len(labels) != len(set(labels)):
        raise ValueError(f"distance labels must be unique: {labels}")
    return out


def validate_config(config: Mapping[str, Any]) -> None:
    """Fail before computation when the registered design is incomplete."""

    required_mappings = (
        "experiment",
        "data",
        "references",
        "replication",
        "legacy_control",
        "candidate",
        "intervention",
        "probe",
        "acceptance",
        "outputs",
    )
    for key in required_mappings:
        if not isinstance(config.get(key), Mapping):
            raise ValueError(f"config.{key} must be a mapping")

    references = config["references"]
    reference_items = (
        list(references.get("synthetic", []))
        + list(references.get("inline", []))
        + list(references.get("datasets", []))
    )
    reference_ids = [str(item.get("id", item.get("dataset_id", ""))) for item in reference_items]
    if any(not value for value in reference_ids):
        raise ValueError("every reference must declare id or dataset_id")
    if len(reference_ids) != len(set(reference_ids)):
        raise ValueError("reference ids must be unique")
    safe_ids = [safe_name(value) for value in reference_ids]
    if len(safe_ids) != len(set(safe_ids)):
        raise ValueError("reference ids collide after filesystem-safe normalization")

    seeds = [int(value) for value in config["replication"].get("seeds", [])]
    if not seeds or len(seeds) != len(set(seeds)):
        raise ValueError("replication.seeds must be non-empty and unique")

    candidate = config["candidate"]
    orders = [int(value) for value in candidate.get("orders", [])]
    order_weights = [float(value) for value in candidate.get("order_weights", [])]
    target_degrees = [float(value) for value in candidate.get("target_degrees", [])]
    if not orders or len({*orders}) != len(orders):
        raise ValueError("candidate.orders must be non-empty and unique")
    if not (len(orders) == len(order_weights) == len(target_degrees)):
        raise ValueError("candidate orders, order_weights, and target_degrees must align")
    if any(value < 0.0 or not np.isfinite(value) for value in order_weights) or sum(order_weights) <= 0.0:
        raise ValueError("candidate.order_weights must be finite, non-negative, and nonzero")
    if any(value <= 0.0 or not np.isfinite(value) for value in target_degrees):
        raise ValueError("candidate.target_degrees must be finite and positive")

    legacy = config["legacy_control"]
    if int(legacy.get("terms", 0)) <= 0:
        raise ValueError("legacy_control.terms must be positive")
    if int(legacy.get("min_order", 0)) <= 0 or int(legacy.get("max_order", 0)) < int(legacy.get("min_order", 0)):
        raise ValueError("legacy_control order range is invalid")

    weights = [float(value) for value in config["intervention"].get("weights", [])]
    if not weights or len(weights) != len(set(weights)) or any(value < 0.0 for value in weights):
        raise ValueError("intervention.weights must be unique, non-negative, and non-empty")
    if not any(np.isclose(value, 0.0) for value in weights) or not any(np.isclose(value, 1.0) for value in weights):
        raise ValueError("intervention.weights must include the causal null 0 and working weight 1")

    probe = config["probe"]
    regimes = {str(value) for value in probe.get("regimes", [])}
    required_regimes = {"reference_shell", "uniform_background"}
    if not required_regimes.issubset(regimes):
        raise ValueError(f"probe.regimes must include {sorted(required_regimes)}")
    distance_labels = {
        label for label, _ in distance_specs(1000, probe.get("hamming_distances", []))
    }
    required_distances = {
        "fixed_1",
        "fixed_2",
        "fixed_4",
        "fixed_8",
        "fraction_0.1",
        "fraction_0.5",
    }
    if not required_distances.issubset(distance_labels):
        raise ValueError(
            f"probe.hamming_distances must generate {sorted(required_distances)}"
        )
    for key in ("pairs_per_distance", "eval_chunk_size", "supported_pair_curvature_samples"):
        if int(probe.get(key, 0)) <= 0:
            raise ValueError(f"probe.{key} must be positive")

    acceptance_keys = {
        "bootstrap_replicates",
        "bootstrap_seed",
        "reference_center_atol",
        "support_max_degree_range",
        "support_max_cv",
        "support_max_duplicate_fraction",
        "contrast_moment_atol",
        "local_ratio_median_band",
        "local_ratio_seed_band",
        "local_ratio_min_seed_fraction",
        "local_slope_abs_max",
        "local_slope_ci_abs_max",
        "min_slope_improvement",
        "field_variance_slope_band",
        "field_variance_per_length_cv_max",
        "fractional_distance_slope_band",
        "distance_ratio_slope_abs_max",
        "fixed_distance_slope_abs_max",
        "weight_linearity_relative_error_max",
        "weight_linearity_r2_min",
        "order_variance_fraction_atol",
        "curvature_slope_abs_max",
        "q20_median_pearson_kurtosis_max",
        "legacy_local_slope_expected_band",
    }
    missing_acceptance = sorted(acceptance_keys - set(config["acceptance"]))
    if missing_acceptance:
        raise ValueError(f"acceptance is missing required keys: {missing_acceptance}")
    if int(config["acceptance"]["bootstrap_replicates"]) <= 0:
        raise ValueError("acceptance.bootstrap_replicates must be positive")


def generate_pairs(
    *,
    reference: np.ndarray,
    seed: int,
    regime: str,
    distance_label: str,
    hamming_distance: int,
    n_pairs: int,
    alphabet_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    length = int(reference.size)
    rng = np.random.default_rng(
        stable_uint64("EXP-060-pairs", int(seed), regime, distance_label, length, n_pairs)
    )
    if regime == "reference_shell":
        sources = np.repeat(reference[None, :], n_pairs, axis=0)
    elif regime == "uniform_background":
        sources = rng.integers(0, alphabet_size, size=(n_pairs, length), dtype=np.int64)
    else:
        raise ValueError(f"unsupported pair regime: {regime}")
    targets = sources.copy()
    for row in range(n_pairs):
        sites = rng.choice(length, size=hamming_distance, replace=False)
        current = targets[row, sites]
        alternatives = rng.integers(0, alphabet_size - 1, size=hamming_distance, dtype=np.int64)
        alternatives += alternatives >= current
        targets[row, sites] = alternatives
    observed = np.sum(sources != targets, axis=1)
    if not np.all(observed == hamming_distance):
        raise AssertionError("pair generator did not preserve exact Hamming distance")
    return sources, targets


def build_additive_core(reference: np.ndarray, seed: int) -> FitnessLandscapeEnv:
    """Build AdaptEnv's additive core with public GNK disabled."""

    config = LandscapeConfig(
        L=int(reference.size),
        alphabet_size=20,
        seed=int(seed),
        gnk_epistasis_weight=0.0,
        gnk_terms=0,
        n_epistatic_terms=0,
        n_higher_order_epistatic_terms=0,
        noise_amplitude=0.0,
        peak_search_steps=0,
        peak_search_top_single_starts=0,
        peak_search_random_starts=0,
        strict_peak_lock=False,
    )
    return FitnessLandscapeEnv(config, reference_sequence=reference)


def retrofit_current_latent_gnk(
    landscape: FitnessLandscapeEnv,
    config: Mapping[str, Any],
) -> None:
    """Build the current public GNK field after core construction.

    The initial environment has ``gnk_terms=0``.  Retrofitting the field here
    means no GNK value is ever passed through the landscape readout.
    """

    legacy = config["legacy_control"]
    landscape.cfg.gnk_terms = int(legacy["terms"])
    landscape.cfg.gnk_min_order = int(legacy["min_order"])
    landscape.cfg.gnk_max_order = int(legacy["max_order"])
    landscape.cfg.gnk_order_weights = tuple(float(value) for value in legacy["order_weights"])
    landscape.cfg.gnk_support_mode = str(legacy.get("support_mode", "structured"))
    landscape.cfg.gnk_epistasis_weight = 1.0
    landscape.cfg.gnk_seed = int(legacy.get("seed_offset", 10000)) + int(landscape.cfg.seed)
    landscape._build_hashed_gnk_epistatic_field()


def term_operator_energy(
    *,
    length: int,
    support: Sequence[int],
    contrasts: Sequence[np.ndarray],
    frequencies: np.ndarray,
) -> float:
    second_moments = [
        float(np.sum(frequencies[site] * contrast * contrast))
        for site, contrast in zip(support, contrasts)
    ]
    delta_moments = [
        uniform_alternative_substitution_energy(contrast, frequencies[site])
        for site, contrast in zip(support, contrasts)
    ]
    total = 0.0
    for mutated_index in range(len(support)):
        contribution = float(delta_moments[mutated_index])
        for other_index, moment in enumerate(second_moments):
            if other_index != mutated_index:
                contribution *= float(moment)
        total += contribution
    return float(total / length)


def compile_legacy_catalog(landscape: FitnessLandscapeEnv) -> tuple[list[dict[str, Any]], float]:
    length = int(landscape.L)
    alphabet_size = int(landscape.A)
    frequencies = np.full((length, alphabet_size), 1.0 / alphabet_size, dtype=np.float64)
    seed = int(landscape._gnk_seed())
    rows: list[dict[str, Any]] = []
    term_id = 0
    for order, count in sorted(landscape._gnk_order_counts.items()):
        order_weight = float(landscape._gnk_order_variance[order])
        coefficient_sd = math.sqrt(order_weight / int(count))
        for _ in range(int(count)):
            support = tuple(int(site) for site in landscape._gnk_support_for_term(term_id, order))
            contrasts = tuple(
                np.asarray(landscape._gnk_contrast_vector(term_id, site), dtype=np.float64)
                for site in support
            )
            coefficient = coefficient_sd * adaptenv_standard_normal(seed, "coeff", term_id, order)
            operator_energy = term_operator_energy(
                length=length,
                support=support,
                contrasts=contrasts,
                frequencies=frequencies,
            )
            rows.append(
                {
                    "term_id": int(term_id),
                    "order": int(order),
                    "support": support,
                    "contrasts": contrasts,
                    "coefficient": float(coefficient),
                    "operator_energy": float(operator_energy),
                }
            )
            term_id += 1
    total_energy = float(sum(row["coefficient"] ** 2 * row["operator_energy"] for row in rows))
    if total_energy <= 0.0:
        raise ValueError("legacy field has zero mutation-operator energy")
    return rows, total_energy


def legacy_catalog_hash(rows: Sequence[Mapping[str, Any]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(np.asarray([row["term_id"], row["order"]], dtype=np.int64).tobytes())
        digest.update(np.asarray(row["support"], dtype=np.int64).tobytes())
        digest.update(np.asarray([row["coefficient"]], dtype=np.float64).tobytes())
        for contrast in row["contrasts"]:
            digest.update(np.asarray(contrast, dtype=np.float64).tobytes())
    return digest.hexdigest()


def score_in_chunks(
    scorer: Callable[[np.ndarray], np.ndarray],
    sequences: np.ndarray,
    *,
    chunk_size: int,
) -> np.ndarray:
    chunks = []
    for start in range(0, sequences.shape[0], chunk_size):
        chunks.append(np.asarray(scorer(sequences[start : start + chunk_size]), dtype=np.float64))
    return np.concatenate(chunks) if chunks else np.empty(0, dtype=np.float64)


def additive_scores(matrix: np.ndarray, sequences: np.ndarray) -> np.ndarray:
    positions = np.arange(matrix.shape[0], dtype=np.int64)
    return np.sum(matrix[positions[None, :], sequences], axis=1)


def raw_moments(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    rms = float(np.sqrt(np.mean(array * array))) if array.size else float("nan")
    mean = float(np.mean(array)) if array.size else float("nan")
    sd = float(np.std(array)) if array.size else float("nan")
    if array.size and sd > 0.0:
        kurtosis = float(np.mean(((array - mean) / sd) ** 4))
    else:
        kurtosis = float("nan")
    p99 = float(np.quantile(np.abs(array), 0.99)) if array.size else float("nan")
    return {
        "rms": rms,
        "mean": mean,
        "sd": sd,
        "pearson_kurtosis": kurtosis,
        "p99_abs_over_rms": float(p99 / rms) if rms > 0.0 else float("nan"),
    }


def support_audit_rows(
    *,
    reference_id: str,
    length: int,
    seed: int,
    arm: str,
    supports_by_order: Mapping[int, Sequence[Sequence[int]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for order, supports in sorted(supports_by_order.items()):
        degrees = np.zeros(length, dtype=np.int64)
        normalized = [tuple(int(site) for site in support) for support in supports]
        for support in normalized:
            degrees[np.asarray(support, dtype=np.int64)] += 1
        mean = float(np.mean(degrees))
        rows.append(
            {
                "reference_id": reference_id,
                "length": int(length),
                "seed": int(seed),
                "arm": arm,
                "order": int(order),
                "term_count": int(len(normalized)),
                "degree_mean": mean,
                "degree_sd": float(np.std(degrees)),
                "degree_cv": float(np.std(degrees) / mean) if mean > 0.0 else float("inf"),
                "degree_min": int(np.min(degrees)),
                "degree_max": int(np.max(degrees)),
                "uncovered_sites": int(np.sum(degrees == 0)),
                "duplicate_support_fraction": float(1.0 - len(set(normalized)) / len(normalized)) if normalized else 0.0,
            }
        )
    return rows


def contrast_audit_rows(
    *,
    reference_id: str,
    length: int,
    seed: int,
    arm: str,
    terms: Sequence[Mapping[str, Any]],
    alphabet_size: int,
) -> list[dict[str, Any]]:
    grouped: dict[int, list[np.ndarray]] = {}
    for term in terms:
        grouped.setdefault(int(term["order"]), []).extend(
            np.asarray(contrast, dtype=np.float64) for contrast in term["contrasts"]
        )
    uniform = np.full(alphabet_size, 1.0 / alphabet_size, dtype=np.float64)
    rows: list[dict[str, Any]] = []
    for order, contrasts in sorted(grouped.items()):
        means = np.asarray([np.mean(contrast) for contrast in contrasts], dtype=np.float64)
        seconds = np.asarray([np.mean(contrast * contrast) for contrast in contrasts], dtype=np.float64)
        maxima = np.asarray([np.max(np.abs(contrast)) for contrast in contrasts], dtype=np.float64)
        deltas = np.asarray(
            [uniform_alternative_substitution_energy(contrast, uniform) for contrast in contrasts],
            dtype=np.float64,
        )
        rows.append(
            {
                "reference_id": reference_id,
                "length": int(length),
                "seed": int(seed),
                "arm": arm,
                "order": int(order),
                "n_contrasts": int(len(contrasts)),
                "max_abs_mean": float(np.max(np.abs(means))),
                "mean_second_moment": float(np.mean(seconds)),
                "max_abs_value": float(np.max(maxima)),
                "mean_substitution_energy": float(np.mean(deltas)),
                "expected_q20_substitution_energy": float(40.0 / 19.0),
            }
        )
    return rows


def legacy_terms_for_audit(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "order": int(row["order"]),
            "support": tuple(row["support"]),
            "contrasts": tuple(row["contrasts"]),
        }
        for row in rows
    ]


def extensive_terms_for_audit(field: ExtensiveGNKField) -> list[dict[str, Any]]:
    return [
        {
            "order": int(term.order),
            "support": tuple(term.support),
            "contrasts": tuple(term.contrasts),
        }
        for term in field.terms
    ]


def sample_supported_curvature(
    field: ExtensiveGNKField,
    *,
    order: int,
    seed: int,
    n_samples: int,
    chunk_size: int,
) -> dict[str, float]:
    pairs = field.support_pairs(order)
    if not pairs:
        return {"curvature_rms": float("nan"), "curvature_kurtosis": float("nan"), "n": 0}
    rng = np.random.default_rng(
        stable_uint64("EXP-060-curvature", int(seed), field.length, int(order))
    )
    base = rng.integers(0, field.alphabet_size, size=(n_samples, field.length), dtype=np.int64)
    first = base.copy()
    second = base.copy()
    both = base.copy()
    for row in range(n_samples):
        left, right = pairs[int(rng.integers(0, len(pairs)))]
        left_alt = int(rng.integers(0, field.alphabet_size - 1))
        right_alt = int(rng.integers(0, field.alphabet_size - 1))
        if left_alt >= int(base[row, left]):
            left_alt += 1
        if right_alt >= int(base[row, right]):
            right_alt += 1
        first[row, left] = left_alt
        second[row, right] = right_alt
        both[row, left] = left_alt
        both[row, right] = right_alt
    combined = np.concatenate([base, first, second, both], axis=0)
    values = score_in_chunks(
        lambda sequences: field.score_batch_order(sequences, order),
        combined,
        chunk_size=chunk_size,
    )
    n = n_samples
    curvature = values[3 * n :] - values[n : 2 * n] - values[2 * n : 3 * n] + values[:n]
    moments = raw_moments(curvature)
    return {
        "curvature_rms": float(moments["rms"]),
        "curvature_kurtosis": float(moments["pearson_kurtosis"]),
        "n": int(n_samples),
    }


def run_task(payload: Mapping[str, Any]) -> dict[str, Any]:
    config = dict(payload["config"])
    reference_row = dict(payload["reference"])
    seed = int(payload["seed"])
    task_signature = str(payload["task_signature"])
    output_path = Path(str(payload["checkpoint_path"]))
    if output_path.exists():
        cached = json.loads(output_path.read_text(encoding="utf-8"))
        if str(cached.get("task", {}).get("task_signature", "")) == task_signature:
            return cached

    reference_id = str(reference_row["reference_id"])
    reference = sequence_to_indices(str(reference_row["sequence"]))
    length = int(reference.size)
    probe = config["probe"]
    chunk_size = int(probe["eval_chunk_size"])
    core = build_additive_core(reference, seed)
    retrofit_current_latent_gnk(core, config)
    additive_matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    additive_energy = additive_mutation_energy(additive_matrix)
    legacy_catalog, legacy_energy = compile_legacy_catalog(core)
    local_scale = math.sqrt(additive_energy / legacy_energy)
    legacy_reference_raw = float(core._gnk_epistatic_components_batch(reference[None, :])[0])
    legacy_hash = legacy_catalog_hash(legacy_catalog)

    candidate = config["candidate"]
    field = ExtensiveGNKField(
        additive_matrix=additive_matrix,
        reference_sequence=reference,
        seed=int(seed) + int(candidate.get("seed_offset", 30000)),
        orders=tuple(int(value) for value in candidate["orders"]),
        order_weights=tuple(float(value) for value in candidate["order_weights"]),
        target_degrees=tuple(float(value) for value in candidate["target_degrees"]),
    )
    replay = ExtensiveGNKField(
        additive_matrix=additive_matrix,
        reference_sequence=reference,
        seed=int(seed) + int(candidate.get("seed_offset", 30000)),
        orders=tuple(int(value) for value in candidate["orders"]),
        order_weights=tuple(float(value) for value in candidate["order_weights"]),
        target_degrees=tuple(float(value) for value in candidate["target_degrees"]),
    )
    replay_match = bool(field.fingerprint == replay.fingerprint)

    def legacy_score(sequences: np.ndarray) -> np.ndarray:
        return np.asarray(core._gnk_epistatic_components_batch(sequences), dtype=np.float64)

    def local_score(sequences: np.ndarray) -> np.ndarray:
        return local_scale * (legacy_score(sequences) - legacy_reference_raw)

    scorers: dict[str, Callable[[np.ndarray], np.ndarray]] = {
        ARM_LEGACY: legacy_score,
        ARM_LOCAL: local_score,
        ARM_EXTENSIVE: field.score_batch,
    }
    reference_scores = {
        arm: float(score_in_chunks(scorer, reference[None, :], chunk_size=chunk_size)[0])
        for arm, scorer in scorers.items()
    }
    field_hashes = {
        ARM_LEGACY: legacy_hash,
        ARM_LOCAL: json_hash({"legacy_hash": legacy_hash, "local_scale": local_scale}),
        ARM_EXTENSIVE: field.fingerprint,
    }

    legacy_supports: dict[int, list[tuple[int, ...]]] = {}
    for row in legacy_catalog:
        legacy_supports.setdefault(int(row["order"]), []).append(tuple(row["support"]))
    candidate_supports = {
        order: [term.support for term in field.terms if term.order == order]
        for order in field.orders
    }
    support_rows = []
    for arm, supports in (
        (ARM_LEGACY, legacy_supports),
        (ARM_LOCAL, legacy_supports),
        (ARM_EXTENSIVE, candidate_supports),
    ):
        support_rows.extend(
            support_audit_rows(
                reference_id=reference_id,
                length=length,
                seed=seed,
                arm=arm,
                supports_by_order=supports,
            )
        )

    contrast_rows = []
    for arm, terms in (
        (ARM_LEGACY, legacy_terms_for_audit(legacy_catalog)),
        (ARM_LOCAL, legacy_terms_for_audit(legacy_catalog)),
        (ARM_EXTENSIVE, extensive_terms_for_audit(field)),
    ):
        contrast_rows.extend(
            contrast_audit_rows(
                reference_id=reference_id,
                length=length,
                seed=seed,
                arm=arm,
                terms=terms,
                alphabet_size=20,
            )
        )

    weights = [float(value) for value in config["intervention"]["weights"]]
    regimes = [str(value) for value in probe["regimes"]]
    n_pairs = int(probe["pairs_per_distance"])
    metric_rows: list[dict[str, Any]] = []
    pairing_rows: list[dict[str, Any]] = []
    order_rows: list[dict[str, Any]] = []
    for regime in regimes:
        for distance_label, hamming in distance_specs(length, probe["hamming_distances"]):
            sources, targets = generate_pairs(
                reference=reference,
                seed=seed,
                regime=regime,
                distance_label=distance_label,
                hamming_distance=hamming,
                n_pairs=n_pairs,
                alphabet_size=20,
            )
            pair_hash = json_hash(
                {
                    "sources": array_hash(sources),
                    "targets": array_hash(targets),
                    "regime": regime,
                    "distance_label": distance_label,
                    "hamming": hamming,
                }
            )
            pairing_rows.append(
                {
                    "reference_id": reference_id,
                    "length": length,
                    "seed": seed,
                    "regime": regime,
                    "distance_label": distance_label,
                    "hamming_distance": hamming,
                    "n_pairs": n_pairs,
                    "pair_hash": pair_hash,
                    "n_arms": len(ARMS),
                    "arm_or_weight_in_sampling_key": False,
                }
            )
            additive_delta = additive_scores(additive_matrix, targets) - additive_scores(additive_matrix, sources)
            additive_moments = raw_moments(additive_delta)
            for arm, scorer in scorers.items():
                combined = np.concatenate([sources, targets], axis=0)
                values = score_in_chunks(scorer, combined, chunk_size=chunk_size)
                source_values = values[:n_pairs]
                base_delta = values[n_pairs:] - values[:n_pairs]
                for weight in weights:
                    delta = weight * base_delta
                    moments = raw_moments(delta)
                    metric_rows.append(
                        {
                            "reference_id": reference_id,
                            "length": length,
                            "seed": seed,
                            "arm": arm,
                            "weight": weight,
                            "weight_semantics": (
                                "raw_public_multiplier" if arm == ARM_LEGACY else "mutation_rms_ratio_to_additive"
                            ),
                            "regime": regime,
                            "distance_label": distance_label,
                            "hamming_distance": hamming,
                            "n_pairs": n_pairs,
                            "pair_hash": pair_hash,
                            "field_hash": field_hashes[arm],
                            "additive_delta_rms": additive_moments["rms"],
                            "gnk_delta_rms": moments["rms"],
                            "gnk_to_additive_rms_ratio": (
                                moments["rms"] / additive_moments["rms"]
                                if additive_moments["rms"] > 0.0
                                else float("nan")
                            ),
                            "gnk_delta_mean": moments["mean"],
                            "gnk_delta_sd": moments["sd"],
                            "gnk_delta_pearson_kurtosis": moments["pearson_kurtosis"],
                            "gnk_delta_p99_abs_over_rms": moments["p99_abs_over_rms"],
                            "source_field_variance": float(np.var(weight * source_values)),
                            "expected_field_variance": (
                                float(weight * weight * field.expected_field_variance())
                                if arm == ARM_EXTENSIVE
                                else float("nan")
                            ),
                            "reference_gnk_score": weight * reference_scores[arm],
                            "all_finite": bool(
                                np.all(np.isfinite(delta))
                                and np.all(np.isfinite(additive_delta))
                            ),
                        }
                    )

                if arm == ARM_EXTENSIVE and regime == "uniform_background" and distance_label == "fixed_1":
                    source_by_order = field.score_batch_by_order(sources)
                    target_by_order = field.score_batch_by_order(targets)
                    energies: dict[int, float] = {}
                    for order in field.orders:
                        order_delta = target_by_order[order] - source_by_order[order]
                        energies[order] = float(np.mean(order_delta * order_delta))
                    total_order_energy = float(sum(energies.values()))
                    for order in field.orders:
                        order_rows.append(
                            {
                                "reference_id": reference_id,
                                "length": length,
                                "seed": seed,
                                "arm": arm,
                                "order": int(order),
                                "target_variance_fraction": float(field.order_weights[order]),
                                "expected_mutation_energy": float(field.expected_order_mutation_energy[order]),
                                "expected_variance_fraction": float(
                                    field.expected_order_mutation_energy[order]
                                    / sum(field.expected_order_mutation_energy.values())
                                ),
                                "sample_mutation_energy": energies[order],
                                "sample_variance_fraction": (
                                    energies[order] / total_order_energy if total_order_energy > 0.0 else float("nan")
                                ),
                                "additive_mutation_energy": additive_energy,
                                "expected_field_variance": field.expected_field_variance(),
                            }
                        )

    curvature_by_order = {
        int(order): sample_supported_curvature(
            field,
            order=int(order),
            seed=seed,
            n_samples=int(probe["supported_pair_curvature_samples"]),
            chunk_size=chunk_size,
        )
        for order in field.orders
    }
    for row in order_rows:
        row.update(curvature_by_order[int(row["order"])])

    result = {
        "task": {
            "reference_id": reference_id,
            "length": length,
            "seed": seed,
            "task_signature": task_signature,
            "replay_match": replay_match,
            "additive_mutation_energy": additive_energy,
            "legacy_expected_mutation_energy": legacy_energy,
            "legacy_local_scale": local_scale,
            "candidate_expected_mutation_energy": float(sum(field.expected_order_mutation_energy.values())),
            "candidate_expected_field_variance": field.expected_field_variance(),
            "candidate_fingerprint": field.fingerprint,
        },
        "pairing_rows": pairing_rows,
        "support_rows": support_rows,
        "contrast_rows": contrast_rows,
        "metric_rows": metric_rows,
        "order_rows": order_rows,
    }
    write_json_atomic(output_path, result)
    return result


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2 or np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        return float("nan"), float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return float(slope), float(intercept), float(r2)


def scale_fit(
    frame: pd.DataFrame,
    *,
    fit_id: str,
    arm: str,
    regime: str,
    distance_label: str,
    metric: str,
    bootstrap_replicates: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    selected = frame[
        (frame["arm"] == arm)
        & (frame["regime"] == regime)
        & (frame["distance_label"] == distance_label)
        & np.isclose(frame["weight"].astype(float), 1.0)
    ].copy()
    selected = selected[np.isfinite(selected[metric].astype(float)) & (selected[metric].astype(float) > 0.0)]
    aggregate = selected.groupby(["reference_id", "length"], as_index=False)[metric].median()
    slope, intercept, r2 = linear_fit(
        np.log(aggregate["length"].to_numpy(dtype=float)),
        np.log(aggregate[metric].to_numpy(dtype=float)),
    )
    seed_values = sorted(int(value) for value in selected["seed"].unique())
    rng = np.random.default_rng(int(bootstrap_seed))
    bootstraps: list[float] = []
    if seed_values and len(aggregate) >= 3:
        by_length = {
            int(length): selected[selected["length"] == length].set_index("seed")[metric].to_dict()
            for length in sorted(selected["length"].unique())
        }
        for _ in range(int(bootstrap_replicates)):
            sampled = rng.choice(seed_values, size=len(seed_values), replace=True)
            lengths: list[int] = []
            medians: list[float] = []
            for length, mapping in by_length.items():
                values = [float(mapping[int(seed)]) for seed in sampled if int(seed) in mapping]
                if values and float(np.median(values)) > 0.0:
                    lengths.append(length)
                    medians.append(float(np.median(values)))
            if len(lengths) >= 3:
                boot_slope, _, _ = linear_fit(np.log(lengths), np.log(medians))
                if np.isfinite(boot_slope):
                    bootstraps.append(float(boot_slope))
    ci_low = float(np.quantile(bootstraps, 0.025)) if bootstraps else float("nan")
    ci_high = float(np.quantile(bootstraps, 0.975)) if bootstraps else float("nan")
    return {
        "fit_id": fit_id,
        "arm": arm,
        "regime": regime,
        "distance_label": distance_label,
        "metric": metric,
        "n_lengths": int(len(aggregate)),
        "n_seeds": int(len(seed_values)),
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "bootstrap_ci_low": ci_low,
        "bootstrap_ci_high": ci_high,
        "bootstrap_replicates": int(len(bootstraps)),
    }


def make_scale_fits(metrics: pd.DataFrame, config: Mapping[str, Any]) -> list[dict[str, Any]]:
    acceptance = config["acceptance"]
    replicates = int(acceptance["bootstrap_replicates"])
    seed = int(acceptance["bootstrap_seed"])
    requests: list[tuple[str, str, str, str, str]] = []
    for arm in ARMS:
        requests.append((f"{arm}_local", arm, "uniform_background", "fixed_1", "gnk_delta_rms"))
    requests.append(("candidate_field_variance", ARM_EXTENSIVE, "uniform_background", "fixed_1", "expected_field_variance"))
    requests.append(
        (
            "candidate_empirical_field_variance",
            ARM_EXTENSIVE,
            "uniform_background",
            "fixed_1",
            "source_field_variance",
        )
    )
    for label in ("fixed_2", "fixed_4", "fixed_8", "fraction_0.1", "fraction_0.5"):
        requests.append((f"candidate_delta_{label}", ARM_EXTENSIVE, "uniform_background", label, "gnk_delta_rms"))
        requests.append(
            (
                f"candidate_ratio_{label}",
                ARM_EXTENSIVE,
                "uniform_background",
                label,
                "gnk_to_additive_rms_ratio",
            )
        )
        requests.append(
            (
                f"candidate_reference_shell_delta_{label}",
                ARM_EXTENSIVE,
                "reference_shell",
                label,
                "gnk_delta_rms",
            )
        )
        requests.append(
            (
                f"candidate_reference_shell_ratio_{label}",
                ARM_EXTENSIVE,
                "reference_shell",
                label,
                "gnk_to_additive_rms_ratio",
            )
        )
    return [
        scale_fit(
            metrics,
            fit_id=fit_id,
            arm=arm,
            regime=regime,
            distance_label=distance,
            metric=metric,
            bootstrap_replicates=replicates,
            bootstrap_seed=seed + index,
        )
        for index, (fit_id, arm, regime, distance, metric) in enumerate(requests)
    ]


def acceptance_checks(
    *,
    tasks: pd.DataFrame,
    pairing: pd.DataFrame,
    supports: pd.DataFrame,
    contrasts: pd.DataFrame,
    metrics: pd.DataFrame,
    orders: pd.DataFrame,
    fits: pd.DataFrame,
    config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    gates = config["acceptance"]
    checks: list[dict[str, Any]] = []

    def add(check_id: str, passed: bool, value: Any, criterion: str, *, required: bool = True) -> None:
        checks.append(
            {
                "check_id": check_id,
                "passed": bool(passed),
                "required": bool(required),
                "value": value,
                "criterion": criterion,
            }
        )

    add(
        "deterministic_replay",
        bool(tasks["replay_match"].all()),
        int(tasks["replay_match"].sum()),
        f"all {len(tasks)} candidate fields reproduce their fingerprint",
    )
    add(
        "pairing_keys_exclude_treatment",
        bool((~pairing["arm_or_weight_in_sampling_key"].astype(bool)).all()),
        int(pairing["pair_hash"].nunique()),
        "all sequence-pair designs are generated once without arm or weight in the key",
    )
    add(
        "finite_latent_outputs",
        bool(metrics["all_finite"].astype(bool).all()),
        int(metrics["all_finite"].astype(bool).sum()),
        "all additive and GNK differences are finite",
    )
    zero = metrics[np.isclose(metrics["weight"].astype(float), 0.0)]
    max_zero = float(zero["gnk_delta_rms"].abs().max()) if len(zero) else float("inf")
    add("weight_zero_exact", max_zero == 0.0, max_zero, "maximum weight-zero GNK RMS equals 0 exactly")
    centered = metrics[
        metrics["arm"].isin([ARM_LOCAL, ARM_EXTENSIVE])
        & np.isclose(metrics["weight"].astype(float), 1.0)
    ]
    max_reference = float(centered["reference_gnk_score"].abs().max())
    center_tolerance = float(gates["reference_center_atol"])
    add(
        "reference_centering",
        max_reference <= center_tolerance,
        max_reference,
        f"max |G(reference)| <= {center_tolerance}",
    )

    primary_support = supports[supports["arm"] == ARM_EXTENSIVE]
    max_degree_range = int((primary_support["degree_max"] - primary_support["degree_min"]).max())
    max_degree_cv = float(primary_support["degree_cv"].max())
    uncovered = int(primary_support["uncovered_sites"].sum())
    max_duplicate = float(primary_support["duplicate_support_fraction"].max())
    add("support_complete", uncovered == 0, uncovered, "zero uncovered sites at every active order")
    add(
        "support_nearly_regular",
        max_degree_range <= int(gates["support_max_degree_range"])
        and max_degree_cv <= float(gates["support_max_cv"]),
        json.dumps({"max_degree_range": max_degree_range, "max_cv": max_degree_cv}),
        f"degree range <= {gates['support_max_degree_range']} and CV <= {gates['support_max_cv']}",
    )
    add(
        "support_unique",
        max_duplicate <= float(gates["support_max_duplicate_fraction"]),
        max_duplicate,
        f"duplicate support fraction <= {gates['support_max_duplicate_fraction']}",
    )

    primary_contrast = contrasts[contrasts["arm"] == ARM_EXTENSIVE]
    mean_error = float(primary_contrast["max_abs_mean"].max())
    second_error = float(np.max(np.abs(primary_contrast["mean_second_moment"].astype(float) - 1.0)))
    max_abs = float(primary_contrast["max_abs_value"].max())
    delta_error = float(
        np.max(
            np.abs(
                primary_contrast["mean_substitution_energy"].astype(float)
                - primary_contrast["expected_q20_substitution_energy"].astype(float)
            )
        )
    )
    contrast_tolerance = float(gates["contrast_moment_atol"])
    add(
        "bounded_q20_contrasts",
        mean_error <= contrast_tolerance
        and second_error <= contrast_tolerance
        and max_abs <= 1.0 + contrast_tolerance
        and delta_error <= contrast_tolerance,
        json.dumps(
            {
                "mean_error": mean_error,
                "second_error": second_error,
                "max_abs": max_abs,
                "substitution_energy_error": delta_error,
            }
        ),
        f"uniform q20 contrast moments exact within {contrast_tolerance} and max |contrast| <= 1",
    )

    local = metrics[
        (metrics["arm"] == ARM_EXTENSIVE)
        & (metrics["regime"] == "uniform_background")
        & (metrics["distance_label"] == "fixed_1")
        & np.isclose(metrics["weight"].astype(float), 1.0)
    ]
    length_medians = local.groupby("length")["gnk_to_additive_rms_ratio"].median()
    median_low, median_high = map(float, gates["local_ratio_median_band"])
    medians_pass = bool(((length_medians >= median_low) & (length_medians <= median_high)).all())
    seed_low, seed_high = map(float, gates["local_ratio_seed_band"])
    seed_fraction = float(
        np.mean(
            (local["gnk_to_additive_rms_ratio"].astype(float) >= seed_low)
            & (local["gnk_to_additive_rms_ratio"].astype(float) <= seed_high)
        )
    )
    add(
        "local_ratio_length_medians",
        medians_pass,
        json.dumps({str(index): float(value) for index, value in length_medians.items()}),
        f"every length median RMS(DeltaG)/RMS(DeltaA) lies in [{median_low}, {median_high}]",
    )
    add(
        "local_ratio_seed_replication",
        seed_fraction >= float(gates["local_ratio_min_seed_fraction"]),
        seed_fraction,
        f"at least {gates['local_ratio_min_seed_fraction']} of seed rows lie in [{seed_low}, {seed_high}]",
    )

    fit_map = fits.set_index("fit_id").to_dict(orient="index")
    candidate_local = fit_map["extensive_balanced_bounded_local"]
    legacy_local = fit_map["legacy_public_global_fixed_local"]
    normalized_control = fit_map["legacy_support_local_normalized_local"]
    normalized_rows = metrics[
        (metrics["arm"] == ARM_LOCAL)
        & (metrics["regime"] == "uniform_background")
        & (metrics["distance_label"] == "fixed_1")
        & np.isclose(metrics["weight"].astype(float), 1.0)
    ]
    normalized_medians = normalized_rows.groupby("length")["gnk_to_additive_rms_ratio"].median()
    normalized_control_pass = bool(
        ((normalized_medians >= median_low) & (normalized_medians <= median_high)).all()
        and abs(float(normalized_control["slope"])) <= float(gates["local_slope_abs_max"])
    )
    add(
        "local_normalized_control_scale",
        normalized_control_pass,
        json.dumps(
            {
                "slope": float(normalized_control["slope"]),
                "length_medians": {
                    str(index): float(value) for index, value in normalized_medians.items()
                },
            }
        ),
        "fixed-support local-normalization control meets the candidate median band and local slope limit",
        required=False,
    )
    slope_limit = float(gates["local_slope_abs_max"])
    ci_limit = float(gates["local_slope_ci_abs_max"])
    add(
        "local_scale_slope",
        abs(float(candidate_local["slope"])) <= slope_limit,
        float(candidate_local["slope"]),
        f"|log-log single-substitution slope| <= {slope_limit}",
    )
    add(
        "local_scale_slope_ci",
        float(candidate_local["bootstrap_ci_low"]) >= -ci_limit
        and float(candidate_local["bootstrap_ci_high"]) <= ci_limit,
        json.dumps(
            {
                "low": float(candidate_local["bootstrap_ci_low"]),
                "high": float(candidate_local["bootstrap_ci_high"]),
            }
        ),
        f"bootstrap slope CI lies within [-{ci_limit}, {ci_limit}]",
    )
    improvement = abs(float(legacy_local["slope"])) - abs(float(candidate_local["slope"]))
    add(
        "slope_repair_over_legacy",
        improvement >= float(gates["min_slope_improvement"]),
        improvement,
        f"candidate slope is at least {gates['min_slope_improvement']} closer to zero than legacy",
    )

    field_fit = fit_map["candidate_field_variance"]
    field_selected = metrics[
        (metrics["arm"] == ARM_EXTENSIVE)
        & (metrics["regime"] == "uniform_background")
        & (metrics["distance_label"] == "fixed_1")
        & np.isclose(metrics["weight"].astype(float), 1.0)
    ].copy()
    field_selected["variance_per_length"] = field_selected["expected_field_variance"] / field_selected["length"]
    var_by_length = field_selected.groupby("length")["variance_per_length"].median()
    var_cv = float(np.std(var_by_length) / np.mean(var_by_length))
    variance_low, variance_high = map(float, gates["field_variance_slope_band"])
    add(
        "extensive_field_variance",
        variance_low <= float(field_fit["slope"]) <= variance_high
        and var_cv <= float(gates["field_variance_per_length_cv_max"]),
        json.dumps({"slope": float(field_fit["slope"]), "variance_per_length_cv": var_cv}),
        f"Var(G) slope in [{variance_low}, {variance_high}] and CV[Var(G)/L] <= {gates['field_variance_per_length_cv_max']}",
    )

    fractional_low, fractional_high = map(float, gates["fractional_distance_slope_band"])
    ratio_slope_limit = float(gates["distance_ratio_slope_abs_max"])
    for label in ("fraction_0.1", "fraction_0.5"):
        delta_fit = fit_map[f"candidate_reference_shell_delta_{label}"]
        ratio_fit = fit_map[f"candidate_ratio_{label}"]
        add(
            f"full_sequence_delta_scaling_{label}",
            fractional_low <= float(delta_fit["slope"]) <= fractional_high,
            float(delta_fit["slope"]),
            f"reference-shell fractional-distance GNK RMS slope lies in [{fractional_low}, {fractional_high}]",
        )
        add(
            f"background_operator_ratio_scaling_{label}",
            abs(float(ratio_fit["slope"])) <= ratio_slope_limit,
            float(ratio_fit["slope"]),
            f"uniform-background fractional-distance GNK/additive ratio slope has absolute value <= {ratio_slope_limit}",
        )
    fixed_slope_limit = float(gates["fixed_distance_slope_abs_max"])
    for label in ("fixed_2", "fixed_4", "fixed_8"):
        fit = fit_map[f"candidate_reference_shell_delta_{label}"]
        add(
            f"fixed_distance_scaling_{label}",
            abs(float(fit["slope"])) <= fixed_slope_limit,
            float(fit["slope"]),
            f"reference-shell fixed-distance GNK RMS slope has absolute value <= {fixed_slope_limit}",
        )

    positive_weights = sorted(float(value) for value in metrics["weight"].unique() if float(value) > 0.0)
    linearity_errors: list[float] = []
    linearity_r2: list[float] = []
    grouped_columns = ["reference_id", "seed", "arm", "regime", "distance_label"]
    for _, group in metrics.groupby(grouped_columns):
        ordered = group.sort_values("weight")
        x = ordered["weight"].to_numpy(dtype=float)
        y = ordered["gnk_delta_rms"].to_numpy(dtype=float)
        slope, _, r2 = linear_fit(x, y)
        predicted = slope * x
        denominator = np.maximum(np.abs(predicted), 1e-12)
        mask = x > 0.0
        if np.any(mask):
            linearity_errors.append(float(np.max(np.abs(y[mask] - predicted[mask]) / denominator[mask])))
        linearity_r2.append(float(r2))
    max_linearity_error = max(linearity_errors) if linearity_errors else float("inf")
    min_linearity_r2 = min(linearity_r2) if linearity_r2 else float("-inf")
    add(
        "weight_linearity",
        max_linearity_error <= float(gates["weight_linearity_relative_error_max"])
        and min_linearity_r2 >= float(gates["weight_linearity_r2_min"]),
        json.dumps({"max_relative_error": max_linearity_error, "min_r2": min_linearity_r2, "weights": positive_weights}),
        f"relative error <= {gates['weight_linearity_relative_error_max']} and R2 >= {gates['weight_linearity_r2_min']}",
    )

    order_error = float(
        np.max(
            np.abs(
                orders.groupby(["length", "order"])["sample_variance_fraction"].median().to_numpy(dtype=float)
                - orders.groupby(["length", "order"])["target_variance_fraction"].first().to_numpy(dtype=float)
            )
        )
    )
    add(
        "order_variance_allocation",
        order_error <= float(gates["order_variance_fraction_atol"]),
        order_error,
        f"maximum median order-share error <= {gates['order_variance_fraction_atol']}",
    )
    curvature_by_length = (
        orders.groupby(["length", "seed", "order"])["curvature_rms"]
        .first()
        .reset_index()
    )
    curvature_slopes: dict[int, float] = {}
    for order, group in curvature_by_length.groupby("order"):
        aggregate = group.groupby("length")["curvature_rms"].median()
        slope, _, _ = linear_fit(
            np.log(aggregate.index.to_numpy(dtype=float)),
            np.log(aggregate.to_numpy(dtype=float)),
        )
        curvature_slopes[int(order)] = float(slope)
    maximum_curvature_slope = max(
        (abs(value) for value in curvature_slopes.values()),
        default=float("inf"),
    )
    add(
        "supported_pair_curvature",
        bool((curvature_by_length["curvature_rms"] > 0.0).all())
        and maximum_curvature_slope <= float(gates["curvature_slope_abs_max"]),
        json.dumps(
            {
                "slope_by_order": curvature_slopes,
                "minimum_rms": float(curvature_by_length["curvature_rms"].min()),
            }
        ),
        f"order-stratified co-supported mixed differences are positive and every length slope <= {gates['curvature_slope_abs_max']} in absolute value",
    )
    tail_median = float(local["gnk_delta_pearson_kurtosis"].median())
    add(
        "q20_local_tail",
        tail_median <= float(gates["q20_median_pearson_kurtosis_max"]),
        tail_median,
        f"median q20 local-effect Pearson kurtosis <= {gates['q20_median_pearson_kurtosis_max']}",
        required=False,
    )

    legacy_band = tuple(float(value) for value in gates["legacy_local_slope_expected_band"])
    add(
        "legacy_scaling_defect_reproduced",
        legacy_band[0] <= float(legacy_local["slope"]) <= legacy_band[1],
        float(legacy_local["slope"]),
        f"legacy local-effect slope lies in [{legacy_band[0]}, {legacy_band[1]}]",
        required=bool(gates.get("require_legacy_defect_reproduction", True)),
    )
    return checks


def classify_scientific_decision(
    checks: Sequence[Mapping[str, Any]],
    config: Mapping[str, Any],
) -> tuple[bool, str]:
    """Map preregistered gate failures to an interpretable latent-stage decision."""

    required = [row for row in checks if bool(row["required"])]
    passed = bool(required and all(bool(row["passed"]) for row in required))
    labels = config.get("decision_rule", {})
    if passed:
        return True, str(labels.get("pass", "LATENT_GNK_EXTENSIVE_SCALE_LAW_VALIDATED"))

    check_map = {str(row["check_id"]): bool(row["passed"]) for row in checks}
    failed = {check_id for check_id, value in check_map.items() if not value}
    technical = {
        "deterministic_replay",
        "pairing_keys_exclude_treatment",
        "finite_latent_outputs",
        "weight_zero_exact",
        "reference_centering",
        "weight_linearity",
    }
    normalization = {
        "local_ratio_length_medians",
        "local_ratio_seed_replication",
        "local_scale_slope",
        "local_scale_slope_ci",
        "slope_repair_over_legacy",
    }
    support = {"support_complete", "support_nearly_regular", "support_unique"}
    contrast = {"bounded_q20_contrasts"}
    full_sequence = {
        value
        for value in failed
        if value.startswith("full_sequence_")
        or value.startswith("background_operator_ratio_")
        or value.startswith("fixed_distance_")
        or value in {"extensive_field_variance", "order_variance_allocation", "supported_pair_curvature"}
    }
    if failed & technical:
        return False, str(labels.get("technical_fail", "EXECUTION_OR_PAIRING_INVALID"))
    if failed & normalization:
        return False, str(labels.get("local_scale_fail", "OPERATOR_NORMALIZATION_FAILED"))
    if bool(check_map.get("local_normalized_control_scale")) and failed & (support | contrast):
        return False, str(
            labels.get(
                "local_normalization_pass_primary_fail",
                "SUPPORT_OR_CONTRAST_REPAIR_INCOMPLETE",
            )
        )
    if failed & support:
        return False, str(
            labels.get("support_fail", "EXTENSIVE_SUPPORT_CONSTRUCTION_FAILED")
        )
    if failed & contrast:
        return False, str(
            labels.get("contrast_fail", "BOUNDED_Q20_CONTRAST_IMPLEMENTATION_FAILED")
        )
    if full_sequence:
        return False, str(
            labels.get(
                "full_sequence_fail",
                "LOCAL_SCALE_DOES_NOT_GENERALIZE_TO_FULL_SEQUENCE_DISTANCE",
            )
        )
    return False, str(labels.get("other_fail", "LATENT_GNK_SCALE_LAW_FAILED_OTHER_GATE"))


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, quick=bool(args.quick))
    validate_config(config)
    references = load_reference_catalog(project_root, config)
    seeds = [int(value) for value in config["replication"]["seeds"]]
    if len(set(seeds)) != len(seeds):
        raise ValueError("replication seeds must be unique")

    scientific_decision_enabled = not bool(
        config.get("experiment", {}).get("disable_scientific_decision", False)
    )
    implementation_files = {
        "runner": Path(__file__).resolve(),
        "field": SRC_ROOT / "project_adapt_env" / "extensive_gnk.py",
    }
    implementation_sha256 = json_hash(
        {
            name: hashlib.sha256(path.read_bytes()).hexdigest()
            for name, path in implementation_files.items()
        }
    )
    preflight = {
        "experiment_id": str(config.get("experiment_id", "EXP-060")),
        "quick": bool(args.quick),
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "config_runtime_hash": json_hash(config),
        "implementation_sha256": implementation_sha256,
        "reference_count": len(references),
        "lengths": [int(row["length"]) for row in references],
        "seed_count": len(seeds),
        "scientific_decision_enabled": scientific_decision_enabled,
        "arms": list(ARMS),
        "gnk_transmitted_through_sigmoid": False,
        "observable_fitness_metrics": False,
        "accessibility_or_optimizer_metrics": False,
        "additive_core_gnk_terms_during_construction": 0,
        "current_gnk_retrofitted_after_core_construction": True,
    }
    write_json_atomic(output_dir / "preflight.json", preflight)
    reference_output = [
        {key: value for key, value in row.items() if key != "sequence"}
        for row in references
    ]
    write_csv_atomic(output_dir / "reference_catalog.csv", reference_output)

    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    task_payloads: list[dict[str, Any]] = []
    for reference in references:
        for seed in seeds:
            task_signature = json_hash(
                {
                    "config": config,
                    "reference": reference,
                    "seed": int(seed),
                    "implementation_sha256": implementation_sha256,
                }
            )
            task_payloads.append(
                {
                    "config": config,
                    "reference": reference,
                    "seed": seed,
                    "task_signature": task_signature,
                    "checkpoint_path": str(
                        checkpoint_dir
                        / (
                            f"{safe_name(reference['reference_id'])}"
                            f"__L-{int(reference['length'])}"
                            f"__seed-{seed}"
                            f"__{task_signature[:16]}.json"
                        )
                    ),
                }
            )
    worker_default = os.environ.get("SLURM_CPUS_PER_TASK", "1")
    workers = max(1, int(os.environ.get("LABPROJ_WORKERS", worker_default)))
    results: list[dict[str, Any]] = []
    if workers == 1:
        for index, payload in enumerate(task_payloads, start=1):
            results.append(run_task(payload))
            print(f"[EXP-060] completed {index}/{len(task_payloads)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_map = {executor.submit(run_task, payload): payload for payload in task_payloads}
            for index, future in enumerate(as_completed(future_map), start=1):
                results.append(future.result())
                print(f"[EXP-060] completed {index}/{len(task_payloads)}", flush=True)

    results.sort(
        key=lambda result: (
            int(result["task"]["length"]),
            str(result["task"]["reference_id"]),
            int(result["task"]["seed"]),
        )
    )

    task_rows = [result["task"] for result in results]
    pairing_rows = [row for result in results for row in result["pairing_rows"]]
    support_rows = [row for result in results for row in result["support_rows"]]
    contrast_rows = [row for result in results for row in result["contrast_rows"]]
    metric_rows = [row for result in results for row in result["metric_rows"]]
    order_rows = [row for result in results for row in result["order_rows"]]

    tasks = pd.DataFrame(task_rows)
    pairing = pd.DataFrame(pairing_rows)
    supports = pd.DataFrame(support_rows)
    contrasts = pd.DataFrame(contrast_rows)
    metrics = pd.DataFrame(metric_rows)
    orders = pd.DataFrame(order_rows)
    fit_rows = make_scale_fits(metrics, config)
    fits = pd.DataFrame(fit_rows)
    checks = acceptance_checks(
        tasks=tasks,
        pairing=pairing,
        supports=supports,
        contrasts=contrasts,
        metrics=metrics,
        orders=orders,
        fits=fits,
        config=config,
    )
    required_checks = [row for row in checks if bool(row["required"])]
    if scientific_decision_enabled:
        passed, decision = classify_scientific_decision(checks, config)
    else:
        passed = False
        decision = "QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION"

    write_csv_atomic(output_dir / "pairing_audit.csv", pairing_rows)
    write_csv_atomic(output_dir / "support_audit.csv", support_rows)
    write_csv_atomic(output_dir / "contrast_audit.csv", contrast_rows)
    write_csv_atomic(output_dir / "latent_scale_metrics.csv", metric_rows)
    write_csv_atomic(output_dir / "order_component_metrics.csv", order_rows)
    write_csv_atomic(output_dir / "scale_law_fits.csv", fit_rows)
    write_csv_atomic(output_dir / "acceptance_checks.csv", checks)
    summary = {
        **preflight,
        "task_count": len(task_rows),
        "metric_row_count": len(metric_rows),
        "required_checks_passed": int(sum(bool(row["passed"]) for row in required_checks)),
        "required_check_count": len(required_checks),
        "scientific_pass": passed if scientific_decision_enabled else None,
        "decision": decision,
        "public_model_changed": False,
        "sigmoid_transmission_tested": False,
        "next_stage_blocked_until_latent_pass": not (
            scientific_decision_enabled and passed
        ),
    }
    write_json_atomic(output_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
