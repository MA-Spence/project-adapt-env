#!/usr/bin/env python3
"""EXP-062 global additive-core calibration and held-out GNK validation.

The experiment has a strict two-phase boundary.  Phase 1 estimates one
mutation-unit scale and selects a sigmoid operating point from additive-core
outcomes only.  Phase 2 freezes that choice before scoring GNK outcomes on
disjoint proteins and crossed additive-core/GNK seeds.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _path in (
    PROJECT_ROOT / "src",
    PROJECT_ROOT / "external" / "Adapt-Env",
    PROJECT_ROOT / "scripts",
):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from project_adapt_env.extensive_gnk import (  # noqa: E402
    ExtensiveGNKField,
    additive_mutation_energy,
    stable_uint64,
)
import extensive_gnk_latent_scaling_exp060 as latent  # noqa: E402
import extensive_gnk_sigmoid_transmission_exp061 as previous  # noqa: E402


PRIMARY_GLOBAL = "combined_global"
EXACT_BETA = "exact_per_landscape_beta"
STEEPNESS_ONLY = "steepness_only_margin8"
THRESHOLD_ONLY = "threshold_only_k20"
RAW_MARGIN8 = "candidate_raw_margin8_w0.6"
RAW_MARGIN10 = "candidate_raw_margin10_w0.6"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def sigmoid(values: np.ndarray | float) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    return 1.0 / (1.0 + np.exp(-np.clip(values, -500.0, 500.0)))


def dimensionless_fitness(
    additive_scores: np.ndarray | float,
    gnk_scores: np.ndarray | float,
    *,
    additive_reference_score: float,
    sigma_a: float,
    sigma_g: float,
    beta: float,
    mu: float,
    amplitude: float,
    output_scale: float,
) -> np.ndarray:
    """Evaluate C sigmoid(mu + beta Z_A + a Z_G)."""

    if float(sigma_a) <= 0.0 or float(sigma_g) <= 0.0:
        raise ValueError("sigma_a and sigma_g must be positive")
    additive = np.asarray(additive_scores, dtype=np.float64)
    gnk = np.asarray(gnk_scores, dtype=np.float64)
    logit = (
        float(mu)
        + float(beta) * (additive - float(additive_reference_score)) / float(sigma_a)
        + float(amplitude) * gnk / float(sigma_g)
    )
    return float(output_scale) * sigmoid(logit)


def paired_residual(
    additive_scores: np.ndarray | float,
    gnk_scores: np.ndarray | float,
    *,
    additive_reference_score: float,
    sigma_a: float,
    sigma_g: float,
    beta: float,
    mu: float,
    amplitude: float,
    output_scale: float,
) -> np.ndarray:
    """Return F(A+aG)-F(A), exactly zero at a=0."""

    additive = np.asarray(additive_scores, dtype=np.float64)
    if float(amplitude) == 0.0:
        return np.zeros_like(additive, dtype=np.float64)
    full = dimensionless_fitness(
        additive,
        gnk_scores,
        additive_reference_score=additive_reference_score,
        sigma_a=sigma_a,
        sigma_g=sigma_g,
        beta=beta,
        mu=mu,
        amplitude=amplitude,
        output_scale=output_scale,
    )
    core = dimensionless_fitness(
        additive,
        np.zeros_like(np.asarray(gnk_scores, dtype=np.float64)),
        additive_reference_score=additive_reference_score,
        sigma_a=sigma_a,
        sigma_g=sigma_g,
        beta=beta,
        mu=mu,
        amplitude=0.0,
        output_scale=output_scale,
    )
    return np.asarray(full - core, dtype=np.float64)


def raw_parameterization(
    *,
    additive_reference_score: float,
    sigma_a: float,
    sigma_g: float,
    beta: float,
    mu: float,
    amplitude: float,
) -> dict[str, float]:
    """Map the dimensionless form to k(A+wG-theta)."""

    steepness = float(beta) / float(sigma_a)
    if steepness <= 0.0 or float(sigma_g) <= 0.0:
        raise ValueError("positive beta, sigma_a, and sigma_g are required")
    return {
        "steepness": steepness,
        "threshold": float(additive_reference_score) - float(mu) / steepness,
        "gnk_weight": float(amplitude) / (steepness * float(sigma_g)),
    }


def realized_additive_logit_rms(
    *, beta: float, task_sigma_a: float, scale_sigma_a: float
) -> float:
    return float(beta) * float(task_sigma_a) / float(scale_sigma_a)


def select_core_calibration(
    rows: Sequence[Mapping[str, Any]],
    *,
    beta_order: Sequence[float],
    mu_order: Sequence[float],
) -> dict[str, Any] | None:
    """Choose the largest eligible beta, then mu, without reading GNK fields."""

    eligible = [row for row in rows if bool(row.get("core_pass", False))]
    if not eligible:
        return None
    beta_rank = {float(value): index for index, value in enumerate(beta_order)}
    mu_rank = {float(value): index for index, value in enumerate(mu_order)}
    return dict(
        max(
            eligible,
            key=lambda row: (
                beta_rank[float(row["beta"])],
                mu_rank[float(row["mu"])],
            ),
        )
    )


def classify_decision(
    *,
    quick: bool,
    technical_pass: bool,
    selection_exists: bool,
    confirmation_core_pass: bool,
    global_transmission_pass: bool,
    global_directionality_pass: bool,
    global_length_stability_pass: bool,
    exact_full_pass: bool,
    frozen_texture_pass: bool,
) -> str:
    if quick:
        return "QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION"
    if not technical_pass:
        return "EXECUTION_OR_PAIRING_INVALID"
    if not selection_exists:
        return "NO_CORE_OPERATING_POINT_QUALIFIES"
    if not confirmation_core_pass:
        return "CORE_CALIBRATION_DOES_NOT_GENERALIZE"
    if global_transmission_pass and frozen_texture_pass:
        if not global_directionality_pass:
            return "CALIBRATED_GNK_TRANSMISSION_NOT_DIRECTIONAL"
        if not global_length_stability_pass:
            return "CALIBRATED_GNK_TRANSMISSION_NOT_LENGTH_STABLE"
        return "GLOBAL_OPERATOR_CALIBRATED_GNK_TRANSMISSION_VALIDATED_NOISELESS"
    if exact_full_pass and frozen_texture_pass:
        return "PER_LANDSCAPE_CORE_NORMALIZATION_REQUIRED"
    if frozen_texture_pass:
        return "GNK_TEXTURE_REMAINS_FROZEN_CORE_ONLY"
    return "CALIBRATED_SIGMOID_TRANSMISSION_FAILED"


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("config root must be a mapping")
    if quick:
        payload = latent.deep_merge(payload, payload.get("quick", {}))
    payload["_runtime"] = {"quick": bool(quick)}
    return payload


def validate_config(config: Mapping[str, Any]) -> None:
    required = {
        "experiment", "references", "replication", "additive_core", "candidate",
        "readout", "calibration", "intervention", "probe", "acceptance", "outputs",
    }
    missing = sorted(required - set(config))
    if missing:
        raise ValueError(f"config missing required sections: {missing}")
    if str(config.get("experiment_id")) != "EXP-062":
        raise ValueError("EXP-062 runner requires experiment_id EXP-062")
    beta = [float(value) for value in config["calibration"]["beta_grid"]]
    occupancies = [float(value) for value in config["calibration"]["reference_occupancies"]]
    if not beta or beta != sorted(set(beta)) or any(value <= 0.0 for value in beta):
        raise ValueError("calibration.beta_grid must be positive, unique, and increasing")
    if not occupancies or occupancies != sorted(set(occupancies)):
        raise ValueError("reference occupancies must be unique and increasing")
    if any(not 0.0 < value < 1.0 for value in occupancies):
        raise ValueError("reference occupancies must lie strictly between zero and one")
    amplitudes = [float(value) for value in config["intervention"]["dimensionless_amplitudes"]]
    for required_amplitude in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
        if required_amplitude not in amplitudes and not bool(config.get("_runtime", {}).get("quick")):
            raise ValueError(f"full amplitude grid is missing {required_amplitude}")
    groups = config["replication"]
    seed_sets = [
        set(int(x) for x in groups["calibration_core_seeds"]),
        set(int(x) for x in groups["confirmation_core_seeds"]),
        set(int(x) for x in groups["confirmation_gnk_seeds"]),
    ]
    if any(not values for values in seed_sets) or any(
        seed_sets[left] & seed_sets[right]
        for left in range(3) for right in range(left + 1, 3)
    ):
        raise ValueError("calibration/core/GNK seed sets must be nonempty and disjoint")


def source_hashes(project_root: Path, config: Mapping[str, Any]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for relative in config["outputs"]["implementation_hash_components"]:
        path = project_root / str(relative)
        if not path.is_file():
            raise FileNotFoundError(f"missing implementation component: {relative}")
        hashes[str(relative)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def load_references(
    project_root: Path, config: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], str]:
    path = project_root / str(config["references"]["catalog"])
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    rows = payload.get("references", [])
    if not isinstance(rows, list) or not rows:
        raise ValueError("reference catalog must contain a nonempty references list")
    selected = set(str(value) for value in config["references"].get("include_ids", []))
    out: list[dict[str, Any]] = []
    for raw in rows:
        row = dict(raw)
        reference_id = str(row["id"])
        if selected and reference_id not in selected:
            continue
        sequence = str(row["sequence"]).strip().upper()
        if len(sequence) != int(row["length"]):
            raise ValueError(f"reference length mismatch: {reference_id}")
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if digest != str(row["sequence_sha256"]):
            raise ValueError(f"reference hash mismatch: {reference_id}")
        if not set(sequence).issubset(set("ACDEFGHIKLMNPQRSTVWY")):
            raise ValueError(f"noncanonical residue in reference: {reference_id}")
        out.append({
            "reference_id": reference_id,
            "split": str(row["split"]),
            "length_bin": str(row["length_bin"]),
            "length": len(sequence),
            "sequence": sequence,
            "sequence_sha256": digest,
            "source_dataset": str(row["source_dataset"]),
            "source_record": str(row["source_record"]),
            "selection_rule": str(row["selection_rule"]),
        })
    if not out:
        raise ValueError("reference filtering produced an empty panel")
    hashes = [row["sequence_sha256"] for row in out]
    if len(hashes) != len(set(hashes)):
        raise ValueError("reference sequences must be unique")
    splits = {row["split"] for row in out}
    if splits != {"calibration", "confirmation"}:
        raise ValueError("both calibration and confirmation references are required")
    if set(row["sequence_sha256"] for row in out if row["split"] == "calibration") & set(
        row["sequence_sha256"] for row in out if row["split"] == "confirmation"
    ):
        raise ValueError("calibration and confirmation sequences overlap")
    if not bool(config.get("_runtime", {}).get("quick")):
        expected_bins = int(config["acceptance"]["expected_length_bins"])
        expected_calibration = int(config["references"]["calibration_per_length_bin"])
        expected_confirmation = int(config["references"]["confirmation_per_length_bin"])
        counts = (
            pd.DataFrame(out)
            .groupby(["length_bin", "split"])
            .size()
            .unstack(fill_value=0)
        )
        if len(counts) != expected_bins:
            raise ValueError(
                f"full panel requires {expected_bins} length bins, found {len(counts)}"
            )
        if not {"calibration", "confirmation"}.issubset(counts.columns):
            raise ValueError("every full-panel length bin requires both splits")
        if not bool((counts["calibration"] == expected_calibration).all()):
            raise ValueError("unexpected calibration reference count in a length bin")
        if not bool((counts["confirmation"] == expected_confirmation).all()):
            raise ValueError("unexpected confirmation reference count in a length bin")
    return sorted(out, key=lambda row: (row["split"], row["length_bin"], row["reference_id"])), hashlib.sha256(path.read_bytes()).hexdigest()


def build_core(reference: np.ndarray, seed: int, config: Mapping[str, Any]):
    return previous.build_additive_core(reference, int(seed), config)


def complete_h1(reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return previous.complete_single_mutants(reference, alphabet_size=20)


def reference_shell_pairs(
    reference: np.ndarray,
    *,
    hamming: int,
    count: int,
    seed: int,
    label: str,
) -> tuple[np.ndarray, np.ndarray, str]:
    rng = np.random.default_rng(
        stable_uint64("EXP-062-reference-shell", int(seed), int(reference.size), label, hamming, count)
    )
    sources = np.repeat(reference[None, :], int(count), axis=0)
    targets = sources.copy()
    for row in range(int(count)):
        sites = rng.choice(reference.size, size=int(hamming), replace=False)
        current = targets[row, sites]
        alternative = rng.integers(0, 19, size=int(hamming), dtype=np.int64)
        alternative += alternative >= current
        targets[row, sites] = alternative
    if not np.all(np.sum(sources != targets, axis=1) == int(hamming)):
        raise AssertionError("reference-shell sampler changed the requested distance")
    design_hash = latent.json_hash({
        "sources": latent.array_hash(sources),
        "targets": latent.array_hash(targets),
        "label": label,
        "hamming": int(hamming),
    })
    return sources, targets, design_hash


def core_scale_task(payload: Mapping[str, Any]) -> dict[str, Any]:
    row = dict(payload["reference"])
    seed = int(payload["core_seed"])
    config = payload["config"]
    reference = latent.sequence_to_indices(str(row["sequence"]))
    core = build_core(reference, seed, config)
    matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    sources, targets = complete_h1(reference)
    a_source = latent.additive_scores(matrix, sources)
    a_target = latent.additive_scores(matrix, targets)
    sigma_a = float(np.sqrt(np.mean((a_target - a_source) ** 2)))
    analytic = math.sqrt(additive_mutation_energy(matrix))
    return {
        "reference_id": row["reference_id"],
        "split": row["split"],
        "length_bin": row["length_bin"],
        "length": int(reference.size),
        "core_seed": seed,
        "sigma_a_complete_h1": sigma_a,
        "sigma_a_analytic_q20": analytic,
        "sigma_ratio_complete_to_analytic": sigma_a / analytic,
        "additive_reference_score": float(latent.additive_scores(matrix, reference[None, :])[0]),
        "public_threshold_margin8": float(core.stability_threshold),
        "additive_matrix_hash": latent.array_hash(matrix),
        "finite": bool(np.isfinite(matrix).all() and np.isfinite(sigma_a) and sigma_a > 0.0),
    }


def core_metrics_for_scale(
    *,
    matrix: np.ndarray,
    reference: np.ndarray,
    reference_score: float,
    scale_sigma_a: float,
    beta: float,
    mu: float,
    output_scale: float,
    probe_seed: int,
    config: Mapping[str, Any],
) -> dict[str, float | bool | str]:
    _, h1_targets = complete_h1(reference)
    h1_a = latent.additive_scores(matrix, h1_targets)
    h1_z = float(mu) + float(beta) * (h1_a - reference_score) / float(scale_sigma_a)
    h1_p = sigmoid(h1_z)
    ref_fitness = float(output_scale) * float(sigmoid(float(mu)))
    h1_fitness = float(output_scale) * h1_p
    fractions = config["probe"]["core_plateau_distances"]
    plateau: dict[str, float] = {}
    hashes: dict[str, str] = {}
    for label, fraction in fractions.items():
        hamming = max(1, min(int(reference.size), int(round(float(fraction) * int(reference.size)))))
        _, targets, design_hash = reference_shell_pairs(
            reference,
            hamming=hamming,
            count=int(config["probe"]["core_plateau_pairs"]),
            seed=probe_seed,
            label=f"core_{label}",
        )
        values = latent.additive_scores(matrix, targets)
        occupancy = sigmoid(float(mu) + float(beta) * (values - reference_score) / float(scale_sigma_a))
        plateau[label] = float(np.mean(occupancy <= 0.05))
        hashes[label] = design_hash
    gates = config["acceptance"]["core"]
    passed = bool(
        ref_fitness >= float(gates["reference_fitness_min"]) - 1e-12
        and float(np.mean((h1_p >= 0.05) & (h1_p <= 0.95))) >= float(gates["h1_dynamic_fraction_min"])
        and float(np.sqrt(np.mean((h1_fitness - ref_fitness) ** 2))) >= float(gates["h1_effect_rms_min"])
        and plateau["fraction_0.1"] >= float(gates["fraction_0.1_low_plateau_min"])
        and plateau["fraction_0.5"] >= float(gates["fraction_0.5_low_plateau_min"])
    )
    return {
        "reference_fitness": ref_fitness,
        "h1_target_dynamic_fraction": float(np.mean((h1_p >= 0.05) & (h1_p <= 0.95))),
        "h1_core_effect_rms": float(np.sqrt(np.mean((h1_fitness - ref_fitness) ** 2))),
        "fraction_0.1_low_plateau_fraction": plateau["fraction_0.1"],
        "fraction_0.5_low_plateau_fraction": plateau["fraction_0.5"],
        "fraction_0.1_pairing_hash": hashes["fraction_0.1"],
        "fraction_0.5_pairing_hash": hashes["fraction_0.5"],
        "core_pass": passed,
    }


def calibration_task(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    row = dict(payload["reference"])
    seed = int(payload["core_seed"])
    config = payload["config"]
    global_sigma_a = float(payload["global_sigma_a"])
    reference = latent.sequence_to_indices(str(row["sequence"]))
    core = build_core(reference, seed, config)
    matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    reference_score = float(latent.additive_scores(matrix, reference[None, :])[0])
    sigma_task = complete_h1_rms(matrix, reference)
    output_scale = float(config["readout"]["output_scale"])
    out: list[dict[str, Any]] = []
    for beta in config["calibration"]["beta_grid"]:
        for occupancy in config["calibration"]["reference_occupancies"]:
            mu = float(math.log(float(occupancy) / (1.0 - float(occupancy))))
            metrics = core_metrics_for_scale(
                matrix=matrix,
                reference=reference,
                reference_score=reference_score,
                scale_sigma_a=global_sigma_a,
                beta=float(beta),
                mu=mu,
                output_scale=output_scale,
                probe_seed=stable_uint64("EXP-062-calibration-core", row["reference_id"], seed),
                config=config,
            )
            out.append({
                "panel": "calibration",
                "scale_mode": PRIMARY_GLOBAL,
                "reference_id": row["reference_id"],
                "length_bin": row["length_bin"],
                "length": int(reference.size),
                "core_seed": seed,
                "beta": float(beta),
                "reference_occupancy": float(occupancy),
                "mu": mu,
                "scale_sigma_a": global_sigma_a,
                "task_sigma_a": sigma_task,
                "realized_beta": realized_additive_logit_rms(
                    beta=float(beta), task_sigma_a=sigma_task, scale_sigma_a=global_sigma_a
                ),
                **metrics,
            })
    return out


def confirmation_core_task(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    row = dict(payload["reference"])
    seed = int(payload["core_seed"])
    config = payload["config"]
    beta = float(payload["beta"])
    mu = float(payload["mu"])
    global_sigma = float(payload["global_sigma_a"])
    reference = latent.sequence_to_indices(str(row["sequence"]))
    core = build_core(reference, seed, config)
    matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    reference_score = float(latent.additive_scores(matrix, reference[None, :])[0])
    sigma_task = complete_h1_rms(matrix, reference)
    rows: list[dict[str, Any]] = []
    for mode, scale in ((PRIMARY_GLOBAL, global_sigma), (EXACT_BETA, sigma_task)):
        metrics = core_metrics_for_scale(
            matrix=matrix,
            reference=reference,
            reference_score=reference_score,
            scale_sigma_a=scale,
            beta=beta,
            mu=mu,
            output_scale=float(config["readout"]["output_scale"]),
            probe_seed=stable_uint64("EXP-062-confirmation-core", row["reference_id"], seed),
            config=config,
        )
        rows.append({
            "panel": "confirmation",
            "scale_mode": mode,
            "reference_id": row["reference_id"],
            "length_bin": row["length_bin"],
            "length": int(reference.size),
            "core_seed": seed,
            "beta": beta,
            "mu": mu,
            "scale_sigma_a": scale,
            "task_sigma_a": sigma_task,
            "realized_beta": realized_additive_logit_rms(
                beta=beta, task_sigma_a=sigma_task, scale_sigma_a=scale
            ),
            **metrics,
        })
    return rows


def raw_paired_residual(
    additive_scores: np.ndarray,
    gnk_scores: np.ndarray,
    *,
    steepness: float,
    threshold: float,
    gnk_weight: float,
    output_scale: float,
) -> np.ndarray:
    additive = np.asarray(additive_scores, dtype=np.float64)
    gnk = np.asarray(gnk_scores, dtype=np.float64)
    if float(gnk_weight) == 0.0:
        return np.zeros_like(additive)
    base_logit = float(steepness) * (additive - float(threshold))
    shift = float(steepness) * float(gnk_weight) * gnk
    return float(output_scale) * (sigmoid(base_logit + shift) - sigmoid(base_logit))


def prefixed_stats(x: np.ndarray, y: np.ndarray, *, prefix: str, floor: float) -> dict[str, Any]:
    return {f"{prefix}{key}": value for key, value in previous.sufficient_stats(
        np.asarray(x, dtype=np.float64),
        np.asarray(y, dtype=np.float64),
        resolution_floor=float(floor),
    ).items()}


def additive_surrogate_matrix(additive_matrix: np.ndarray, *, seed: int) -> np.ndarray:
    matrix = np.asarray(additive_matrix, dtype=np.float64)
    rng = np.random.default_rng(stable_uint64("EXP-062-additive-surrogate", int(seed), matrix.shape))
    surrogate = rng.standard_normal(matrix.shape)
    surrogate -= np.mean(surrogate, axis=1, keepdims=True)
    return surrogate * math.sqrt(additive_mutation_energy(matrix) / additive_mutation_energy(surrogate))


def score_additive_field(
    matrix: np.ndarray, reference: np.ndarray, sequences: np.ndarray
) -> np.ndarray:
    return previous.score_additive_field(matrix, reference, sequences)


def generate_squares(
    field: ExtensiveGNKField,
    reference: np.ndarray,
    *,
    core_seed: int,
    gnk_seed: int,
    stratum: str,
    background: str,
    count: int,
) -> dict[str, Any]:
    """Generate a treatment-independent deterministic plaquette panel."""

    rng = np.random.default_rng(stable_uint64(
        "EXP-062-squares", int(core_seed), int(gnk_seed), int(reference.size),
        str(stratum), str(background), int(count),
    ))
    if stratum.startswith("candidate_supported_order_"):
        order = int(stratum.rsplit("_", 1)[1])
        pair_pool = field.support_pairs(order)
        if not pair_pool:
            raise ValueError(f"no candidate support pairs at order {order}")
    elif stratum == "uniform_site_pair":
        pair_pool = ()
    else:
        raise ValueError(f"unsupported square stratum: {stratum}")
    if background == "reference":
        base = np.repeat(reference[None, :], int(count), axis=0)
    elif background == "uniform_q20":
        base = rng.integers(0, field.alphabet_size, size=(int(count), field.length), dtype=np.int64)
    else:
        raise ValueError(f"unsupported square background: {background}")
    first = base.copy()
    second = base.copy()
    both = base.copy()
    sites = np.empty((int(count), 2), dtype=np.int64)
    for row in range(int(count)):
        if pair_pool:
            left, right = pair_pool[int(rng.integers(0, len(pair_pool)))]
        else:
            left, right = sorted(rng.choice(field.length, size=2, replace=False).tolist())
        sites[row] = (int(left), int(right))
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
    square_hash = latent.json_hash({
        "base": latent.array_hash(base),
        "first": latent.array_hash(first),
        "second": latent.array_hash(second),
        "both": latent.array_hash(both),
        "sites": latent.array_hash(sites),
        "stratum": stratum,
        "background": background,
    })
    return {
        "base": base,
        "first": first,
        "second": second,
        "both": both,
        "sites": sites,
        "square_hash": square_hash,
    }


def score_nodes(
    scorer: Callable[[np.ndarray], np.ndarray],
    arrays: Sequence[np.ndarray],
    *,
    chunk_size: int,
) -> list[np.ndarray]:
    return previous.score_nodes(scorer, arrays, chunk_size=int(chunk_size))


def dimensionless_residual_and_stable(
    additive: np.ndarray,
    gnk: np.ndarray,
    *,
    reference_score: float,
    scale_sigma_a: float,
    sigma_g: float,
    beta: float,
    mu: float,
    amplitude: float,
    output_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    additive = np.asarray(additive, dtype=np.float64)
    gnk = np.asarray(gnk, dtype=np.float64)
    direct = paired_residual(
        additive,
        gnk,
        additive_reference_score=reference_score,
        sigma_a=scale_sigma_a,
        sigma_g=sigma_g,
        beta=beta,
        mu=mu,
        amplitude=amplitude,
        output_scale=output_scale,
    )
    base_logit = float(mu) + float(beta) * (additive - reference_score) / float(scale_sigma_a)
    shift = float(amplitude) * gnk / float(sigma_g)
    if float(amplitude) == 0.0:
        stable = np.zeros_like(direct)
    else:
        # The inputs are already logits.  A steepness of one avoids the EXP-061
        # unit pitfall of multiplying a dimensionless shift by the public k=20.
        stable = previous.stable_sigmoid_residual(
            base_logit, shift, output_scale=output_scale, steepness=1.0
        )
    return direct, stable, base_logit, shift


def raw_residual_and_stable(
    additive: np.ndarray,
    gnk: np.ndarray,
    *,
    steepness: float,
    threshold: float,
    gnk_weight: float,
    output_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    direct = raw_paired_residual(
        additive,
        gnk,
        steepness=steepness,
        threshold=threshold,
        gnk_weight=gnk_weight,
        output_scale=output_scale,
    )
    base_logit = float(steepness) * (np.asarray(additive, dtype=np.float64) - float(threshold))
    shift = float(steepness) * float(gnk_weight) * np.asarray(gnk, dtype=np.float64)
    stable = previous.stable_sigmoid_residual(
        base_logit, shift, output_scale=output_scale, steepness=1.0
    )
    return direct, stable, base_logit, shift


def edge_metric_row(
    *,
    ids: Mapping[str, Any],
    sample_id: str,
    pairing_hash: str,
    hamming_distance: int,
    a_source: np.ndarray,
    a_target: np.ndarray,
    g_source: np.ndarray,
    g_target: np.ndarray,
    mode: str,
    amplitude: float | None,
    scale_sigma_a: float,
    sigma_g: float,
    beta: float,
    mu: float,
    reference_score: float,
    public_threshold: float,
    output_scale: float,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    if mode in {PRIMARY_GLOBAL, EXACT_BETA}:
        if amplitude is None:
            raise ValueError("dimensionless modes require an amplitude")
        r_source, stable_source, _, shift_source = dimensionless_residual_and_stable(
            a_source, g_source,
            reference_score=reference_score,
            scale_sigma_a=scale_sigma_a,
            sigma_g=sigma_g,
            beta=beta,
            mu=mu,
            amplitude=float(amplitude),
            output_scale=output_scale,
        )
        r_target, stable_target, _, shift_target = dimensionless_residual_and_stable(
            a_target, g_target,
            reference_score=reference_score,
            scale_sigma_a=scale_sigma_a,
            sigma_g=sigma_g,
            beta=beta,
            mu=mu,
            amplitude=float(amplitude),
            output_scale=output_scale,
        )
        raw = raw_parameterization(
            additive_reference_score=reference_score,
            sigma_a=scale_sigma_a,
            sigma_g=sigma_g,
            beta=beta,
            mu=mu,
            amplitude=float(amplitude),
        )
    else:
        global_k = beta / scale_sigma_a
        if mode == STEEPNESS_ONLY:
            k = global_k
            theta = public_threshold
            weight = 1.0 / (k * sigma_g)
            reported_amplitude = 1.0
        elif mode == THRESHOLD_ONLY:
            k = float(config["readout"]["public_steepness"])
            theta = reference_score - mu / k
            weight = 1.0 / (k * sigma_g)
            reported_amplitude = 1.0
        elif mode == RAW_MARGIN8:
            k = float(config["readout"]["public_steepness"])
            theta = public_threshold
            weight = float(config["readout"]["public_gnk_weight"])
            reported_amplitude = k * weight * sigma_g
        elif mode == RAW_MARGIN10:
            k = float(config["readout"]["public_steepness"])
            theta = public_threshold - float(config["readout"]["constructor_margin_increment"])
            weight = float(config["readout"]["public_gnk_weight"])
            reported_amplitude = k * weight * sigma_g
        else:
            raise ValueError(f"unsupported edge mode: {mode}")
        r_source, stable_source, _, shift_source = raw_residual_and_stable(
            a_source, g_source, steepness=k, threshold=theta,
            gnk_weight=weight, output_scale=output_scale,
        )
        r_target, stable_target, _, shift_target = raw_residual_and_stable(
            a_target, g_target, steepness=k, threshold=theta,
            gnk_weight=weight, output_scale=output_scale,
        )
        amplitude = reported_amplitude
        raw = {"steepness": k, "threshold": theta, "gnk_weight": weight}
    x = np.asarray(shift_target - shift_source, dtype=np.float64)
    direct = np.asarray(r_target - r_source, dtype=np.float64)
    stable = np.asarray(stable_target - stable_source, dtype=np.float64)
    floor = float(config["acceptance"]["transmission"]["operational_detection_threshold"])
    return {
        **ids,
        "field_arm": "extensive_balanced_bounded",
        "scale_mode": mode,
        "additive_core_mode": "actual_additive_core",
        "sample_id": sample_id,
        "hamming_distance": int(hamming_distance),
        "pairing_hash": pairing_hash,
        "amplitude": float(amplitude),
        "scale_sigma_a": float(scale_sigma_a),
        "sigma_g": float(sigma_g),
        "realized_beta": realized_additive_logit_rms(
            beta=beta,
            task_sigma_a=float(ids["task_sigma_a"]),
            scale_sigma_a=scale_sigma_a,
        ),
        "raw_steepness": float(raw["steepness"]),
        "raw_threshold": float(raw["threshold"]),
        "raw_gnk_weight": float(raw["gnk_weight"]),
        "delta_a_rms": float(np.sqrt(np.mean((a_target - a_source) ** 2))),
        **prefixed_stats(x, direct, prefix="direct_", floor=floor),
        **prefixed_stats(x, stable, prefix="stable_", floor=floor),
        "max_direct_stable_abs_difference": float(np.max(np.abs(direct - stable))),
        "all_finite": bool(
            np.isfinite(x).all() and np.isfinite(direct).all() and np.isfinite(stable).all()
        ),
    }


def texture_arrays(
    *,
    additive_nodes: Sequence[np.ndarray],
    gnk_nodes: Sequence[np.ndarray],
    core_mode: str,
    reference_score: float,
    scale_sigma_a: float,
    sigma_g: float,
    beta: float,
    mu: float,
    amplitude: float,
    output_scale: float,
) -> dict[str, np.ndarray]:
    a0, aa, ab, aab = [np.asarray(value, dtype=np.float64) for value in additive_nodes]
    g0, ga, gb, gab = [np.asarray(value, dtype=np.float64) for value in gnk_nodes]
    if core_mode == "frozen_additive_core":
        core_nodes = [a0, a0, a0, a0]
    elif core_mode == "actual_additive_core":
        core_nodes = [a0, aa, ab, aab]
    else:
        raise ValueError(f"unsupported core mode: {core_mode}")
    residual_nodes: list[np.ndarray] = []
    for additive, gnk in zip(core_nodes, (g0, ga, gb, gab)):
        direct, _, _, _ = dimensionless_residual_and_stable(
            additive,
            gnk,
            reference_score=reference_score,
            scale_sigma_a=scale_sigma_a,
            sigma_g=sigma_g,
            beta=beta,
            mu=mu,
            amplitude=amplitude,
            output_scale=output_scale,
        )
        residual_nodes.append(direct)
    epsilon_r = residual_nodes[3] - residual_nodes[1] - residual_nodes[2] + residual_nodes[0]
    latent_epsilon = gab - ga - gb + g0
    additivized_gab = ga + gb - g0
    z_ab = float(mu) + float(beta) * (core_nodes[3] - reference_score) / float(scale_sigma_a)
    base = z_ab + float(amplitude) * additivized_gab / float(sigma_g)
    shift = float(amplitude) * latent_epsilon / float(sigma_g)
    if float(amplitude) == 0.0:
        tau = np.zeros_like(shift)
        tau_stable = np.zeros_like(shift)
    else:
        tau = float(output_scale) * (sigmoid(base + shift) - sigmoid(base))
        tau_stable = previous.stable_sigmoid_residual(
            base, shift, output_scale=output_scale, steepness=1.0
        )
    return {
        "x": np.asarray(shift, dtype=np.float64),
        "tau": np.asarray(tau, dtype=np.float64),
        "tau_stable": np.asarray(tau_stable, dtype=np.float64),
        "epsilon_r": np.asarray(epsilon_r, dtype=np.float64),
        "latent_epsilon": np.asarray(latent_epsilon, dtype=np.float64),
    }


def stratified_permutation_projection(
    blocks: Sequence[tuple[np.ndarray, np.ndarray]],
    *,
    replicates: int,
    seed: int,
) -> dict[str, float | int]:
    """Project after permuting the latent predictor within support-order blocks."""

    parsed = [
        (np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
        for x, y in blocks
    ]
    denominator = float(sum(np.dot(x, x) for x, _ in parsed))
    if not parsed or denominator <= 0.0 or any(x.size < 2 for x, _ in parsed):
        return {
            "observed_projection": float("nan"),
            "null_mean": float("nan"),
            "null_q95": float("nan"),
            "one_sided_pvalue": float("nan"),
            "permutation_seed": int(seed),
        }
    observed = float(sum(np.dot(x, y) for x, y in parsed) / denominator)
    rng = np.random.default_rng(int(seed))
    null = np.empty(int(replicates), dtype=np.float64)
    for index in range(int(replicates)):
        numerator = float(sum(np.dot(rng.permutation(x), y) for x, y in parsed))
        null[index] = numerator / denominator
    return {
        "observed_projection": observed,
        "null_mean": float(np.mean(null)),
        "null_q95": float(np.quantile(null, 0.95)),
        "one_sided_pvalue": float((1 + np.sum(null >= observed)) / (int(replicates) + 1)),
        "permutation_seed": int(seed),
    }


def texture_metric_row(
    *,
    ids: Mapping[str, Any],
    scale_mode: str,
    core_mode: str,
    field_arm: str,
    amplitude: float,
    stratum: str,
    background: str,
    square_hash: str,
    arrays: Mapping[str, np.ndarray],
    scale_sigma_a: float,
    sigma_g: float,
    beta: float,
    config: Mapping[str, Any],
    pooled: bool,
    permutation_blocks: Sequence[tuple[np.ndarray, np.ndarray]] | None = None,
) -> dict[str, Any]:
    x = np.asarray(arrays["x"], dtype=np.float64)
    tau = np.asarray(arrays["tau"], dtype=np.float64)
    stable = np.asarray(arrays["tau_stable"], dtype=np.float64)
    floor = float(config["acceptance"]["transmission"]["operational_detection_threshold"])
    permutation_required = bool(
        pooled
        and field_arm == "extensive_balanced_bounded"
        and float(amplitude) == 1.0
        and scale_mode in {PRIMARY_GLOBAL, EXACT_BETA}
    )
    if permutation_required and float(np.var(x)) > 0.0:
        permutation_seed = stable_uint64(
            "EXP-062-permutation", ids["reference_id"], ids["core_seed"],
            ids["gnk_seed"], scale_mode, core_mode,
        )
        if permutation_blocks is None:
            raise ValueError("pooled primary permutation requires support-order blocks")
        permutation = stratified_permutation_projection(
            permutation_blocks,
            replicates=int(config["acceptance"]["permutation_replicates"]),
            seed=int(permutation_seed),
        )
    else:
        permutation = {
            "observed_projection": float("nan"),
            "null_mean": float("nan"),
            "null_q95": float("nan"),
            "one_sided_pvalue": float("nan"),
            "permutation_seed": -1,
        }
    return {
        **ids,
        "scale_mode": scale_mode,
        "additive_core_mode": core_mode,
        "field_arm": field_arm,
        "amplitude": float(amplitude),
        "square_stratum": stratum,
        "square_background": background,
        "square_hash": square_hash,
        "pooled_supported_reference": bool(pooled),
        "scale_sigma_a": float(scale_sigma_a),
        "sigma_g": float(sigma_g),
        "realized_beta": realized_additive_logit_rms(
            beta=beta,
            task_sigma_a=float(ids["task_sigma_a"]),
            scale_sigma_a=scale_sigma_a,
        ),
        "latent_epsilon_rms": float(np.sqrt(np.mean(np.asarray(arrays["latent_epsilon"]) ** 2))),
        "latent_epsilon_max_abs": float(np.max(np.abs(np.asarray(arrays["latent_epsilon"])))),
        "epsilon_r_rms": float(np.sqrt(np.mean(np.asarray(arrays["epsilon_r"]) ** 2))),
        **prefixed_stats(x, tau, prefix="direct_tau_", floor=floor),
        **prefixed_stats(x, stable, prefix="stable_tau_", floor=floor),
        **permutation,
        "permutation_required": permutation_required,
        "permutation_pass": bool(
            permutation_required
            and float(permutation["one_sided_pvalue"])
            <= float(config["acceptance"]["texture_permutation_alpha"])
        ) if permutation_required else None,
        "max_direct_stable_abs_difference": float(np.max(np.abs(tau - stable))),
        "all_finite": bool(
            np.isfinite(x).all()
            and np.isfinite(tau).all()
            and np.isfinite(stable).all()
            and np.isfinite(np.asarray(arrays["epsilon_r"])).all()
        ),
    }


def complete_h1_rms(matrix: np.ndarray, reference: np.ndarray) -> float:
    sources, targets = complete_h1(reference)
    delta = latent.additive_scores(matrix, targets) - latent.additive_scores(matrix, sources)
    return float(np.sqrt(np.mean(delta * delta)))


def validation_task(payload: Mapping[str, Any]) -> dict[str, Any]:
    config = payload["config"]
    row = dict(payload["reference"])
    core_seed = int(payload["core_seed"])
    gnk_seed = int(payload["gnk_seed"])
    beta = float(payload["beta"])
    mu = float(payload["mu"])
    global_sigma_a = float(payload["global_sigma_a"])
    signature = str(payload["task_signature"])
    checkpoint = Path(str(payload["checkpoint_path"]))
    if checkpoint.exists():
        cached = json.loads(checkpoint.read_text(encoding="utf-8"))
        if str(cached.get("technical", {}).get("task_signature")) == signature:
            return cached

    reference = latent.sequence_to_indices(str(row["sequence"]))
    core = build_core(reference, core_seed, config)
    matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    reference_score = float(latent.additive_scores(matrix, reference[None, :])[0])
    task_sigma_a = complete_h1_rms(matrix, reference)
    sigma_g = math.sqrt(additive_mutation_energy(matrix))
    candidate_cfg = config["candidate"]
    field_seed = int(gnk_seed) + int(candidate_cfg.get("seed_offset", 30000))
    candidate = ExtensiveGNKField(
        additive_matrix=matrix,
        reference_sequence=reference,
        seed=field_seed,
        orders=tuple(int(value) for value in candidate_cfg["orders"]),
        order_weights=tuple(float(value) for value in candidate_cfg["order_weights"]),
        target_degrees=tuple(float(value) for value in candidate_cfg["target_degrees"]),
    )
    replay = ExtensiveGNKField(
        additive_matrix=matrix,
        reference_sequence=reference,
        seed=field_seed,
        orders=tuple(int(value) for value in candidate_cfg["orders"]),
        order_weights=tuple(float(value) for value in candidate_cfg["order_weights"]),
        target_degrees=tuple(float(value) for value in candidate_cfg["target_degrees"]),
    )
    surrogate = additive_surrogate_matrix(
        matrix,
        seed=gnk_seed + int(config["candidate"].get("surrogate_seed_offset", 50000)),
    )
    surrogate_sigma_g = math.sqrt(additive_mutation_energy(surrogate))

    def score_candidate(sequences: np.ndarray) -> np.ndarray:
        return candidate.score_batch(sequences)

    def score_surrogate(sequences: np.ndarray) -> np.ndarray:
        return score_additive_field(surrogate, reference, sequences)

    ids = {
        "reference_id": row["reference_id"],
        "length_bin": row["length_bin"],
        "length": int(reference.size),
        "core_seed": core_seed,
        "gnk_seed": gnk_seed,
        "task_sigma_a": task_sigma_a,
    }
    output_scale = float(config["readout"]["output_scale"])
    public_threshold = float(core.stability_threshold)
    amplitudes = [float(value) for value in config["intervention"]["dimensionless_amplitudes"]]
    edge_rows: list[dict[str, Any]] = []
    pair_replay = True
    pair_panels: list[tuple[str, int, np.ndarray, np.ndarray, str]] = []
    h1_sources, h1_targets = complete_h1(reference)
    h1_hash = latent.json_hash({
        "sources": latent.array_hash(h1_sources),
        "targets": latent.array_hash(h1_targets),
        "label": "complete_h1",
    })
    pair_panels.append(("reference_shell_h1_complete", 1, h1_sources, h1_targets, h1_hash))
    for label, fraction in config["probe"]["validation_distances"].items():
        hamming = max(1, min(int(reference.size), int(round(float(fraction) * int(reference.size)))))
        sources, targets, design_hash = reference_shell_pairs(
            reference,
            hamming=hamming,
            count=int(config["probe"]["validation_pairs_per_distance"]),
            seed=stable_uint64("EXP-062-validation-pairs", row["reference_id"], core_seed, gnk_seed),
            label=str(label),
        )
        replay_sources, replay_targets, replay_hash = reference_shell_pairs(
            reference,
            hamming=hamming,
            count=int(config["probe"]["validation_pairs_per_distance"]),
            seed=stable_uint64("EXP-062-validation-pairs", row["reference_id"], core_seed, gnk_seed),
            label=str(label),
        )
        pair_replay = bool(
            pair_replay
            and design_hash == replay_hash
            and np.array_equal(sources, replay_sources)
            and np.array_equal(targets, replay_targets)
        )
        pair_panels.append((str(label), hamming, sources, targets, design_hash))

    max_formula_parity = 0.0
    raw_spec = raw_parameterization(
        additive_reference_score=reference_score,
        sigma_a=global_sigma_a,
        sigma_g=sigma_g,
        beta=beta,
        mu=mu,
        amplitude=1.0,
    )
    for sample_id, hamming, sources, targets, pairing_hash in pair_panels:
        a_source = latent.additive_scores(matrix, sources)
        a_target = latent.additive_scores(matrix, targets)
        g_source, g_target = score_nodes(
            score_candidate, [sources, targets], chunk_size=int(config["probe"]["eval_chunk_size"])
        )
        dimensionless_probe = dimensionless_fitness(
            a_target[: min(64, len(a_target))],
            g_target[: min(64, len(g_target))],
            additive_reference_score=reference_score,
            sigma_a=global_sigma_a,
            sigma_g=sigma_g,
            beta=beta,
            mu=mu,
            amplitude=1.0,
            output_scale=output_scale,
        )
        raw_probe = output_scale * sigmoid(
            raw_spec["steepness"]
            * (
                a_target[: len(dimensionless_probe)]
                + raw_spec["gnk_weight"] * g_target[: len(dimensionless_probe)]
                - raw_spec["threshold"]
            )
        )
        max_formula_parity = max(
            max_formula_parity,
            float(np.max(np.abs(dimensionless_probe - raw_probe))),
        )
        for mode, scale in ((PRIMARY_GLOBAL, global_sigma_a), (EXACT_BETA, task_sigma_a)):
            for amplitude in amplitudes:
                edge_rows.append(edge_metric_row(
                    ids=ids,
                    sample_id=sample_id,
                    pairing_hash=pairing_hash,
                    hamming_distance=hamming,
                    a_source=a_source,
                    a_target=a_target,
                    g_source=g_source,
                    g_target=g_target,
                    mode=mode,
                    amplitude=amplitude,
                    scale_sigma_a=scale,
                    sigma_g=sigma_g,
                    beta=beta,
                    mu=mu,
                    reference_score=reference_score,
                    public_threshold=public_threshold,
                    output_scale=output_scale,
                    config=config,
                ))
        for mode in (STEEPNESS_ONLY, THRESHOLD_ONLY, RAW_MARGIN8, RAW_MARGIN10):
            edge_rows.append(edge_metric_row(
                ids=ids,
                sample_id=sample_id,
                pairing_hash=pairing_hash,
                hamming_distance=hamming,
                a_source=a_source,
                a_target=a_target,
                g_source=g_source,
                g_target=g_target,
                mode=mode,
                amplitude=None,
                scale_sigma_a=global_sigma_a,
                sigma_g=sigma_g,
                beta=beta,
                mu=mu,
                reference_score=reference_score,
                public_threshold=public_threshold,
                output_scale=output_scale,
                config=config,
            ))

    texture_rows: list[dict[str, Any]] = []
    pooled: dict[tuple[str, str, float, str], list[dict[str, np.ndarray]]] = {}
    pooled_hashes: dict[tuple[str, str, float, str], list[str]] = {}
    square_replay = True
    square_hashes: list[str] = []
    surrogate_latent_epsilon_max = 0.0
    for stratum in config["probe"]["square_strata"]:
        for background in config["probe"]["square_backgrounds"]:
            square = generate_squares(
                candidate,
                reference,
                core_seed=core_seed,
                gnk_seed=gnk_seed,
                stratum=str(stratum),
                background=str(background),
                count=int(config["probe"]["squares_per_stratum"]),
            )
            square_again = generate_squares(
                candidate,
                reference,
                core_seed=core_seed,
                gnk_seed=gnk_seed,
                stratum=str(stratum),
                background=str(background),
                count=int(config["probe"]["squares_per_stratum"]),
            )
            square_replay = bool(square_replay and square["square_hash"] == square_again["square_hash"])
            square_hashes.append(str(square["square_hash"]))
            nodes = [square["base"], square["first"], square["second"], square["both"]]
            additive_nodes = [latent.additive_scores(matrix, values) for values in nodes]
            gnk_nodes = score_nodes(
                score_candidate, nodes, chunk_size=int(config["probe"]["eval_chunk_size"])
            )
            for mode, scale in ((PRIMARY_GLOBAL, global_sigma_a), (EXACT_BETA, task_sigma_a)):
                for core_mode in ("actual_additive_core", "frozen_additive_core"):
                    for amplitude in amplitudes:
                        arrays = texture_arrays(
                            additive_nodes=additive_nodes,
                            gnk_nodes=gnk_nodes,
                            core_mode=core_mode,
                            reference_score=reference_score,
                            scale_sigma_a=scale,
                            sigma_g=sigma_g,
                            beta=beta,
                            mu=mu,
                            amplitude=amplitude,
                            output_scale=output_scale,
                        )
                        texture_rows.append(texture_metric_row(
                            ids=ids,
                            scale_mode=mode,
                            core_mode=core_mode,
                            field_arm="extensive_balanced_bounded",
                            amplitude=amplitude,
                            stratum=str(stratum),
                            background=str(background),
                            square_hash=str(square["square_hash"]),
                            arrays=arrays,
                            scale_sigma_a=scale,
                            sigma_g=sigma_g,
                            beta=beta,
                            config=config,
                            pooled=False,
                            permutation_blocks=None,
                        ))
                        if str(background) == "reference" and str(stratum).startswith("candidate_supported_order_"):
                            key = (mode, core_mode, amplitude, "extensive_balanced_bounded")
                            pooled.setdefault(key, []).append(arrays)
                            pooled_hashes.setdefault(key, []).append(str(square["square_hash"]))
            if str(background) == "reference" and str(stratum).startswith("candidate_supported_order_"):
                surrogate_nodes = score_nodes(
                    score_surrogate, nodes, chunk_size=int(config["probe"]["eval_chunk_size"])
                )
                surrogate_arrays = texture_arrays(
                    additive_nodes=additive_nodes,
                    gnk_nodes=surrogate_nodes,
                    core_mode="actual_additive_core",
                    reference_score=reference_score,
                    scale_sigma_a=global_sigma_a,
                    sigma_g=surrogate_sigma_g,
                    beta=beta,
                    mu=mu,
                    amplitude=1.0,
                    output_scale=output_scale,
                )
                surrogate_latent_epsilon_max = max(
                    surrogate_latent_epsilon_max,
                    float(np.max(np.abs(surrogate_arrays["latent_epsilon"]))),
                )
                key = (PRIMARY_GLOBAL, "actual_additive_core", 1.0, "additive_surrogate")
                pooled.setdefault(key, []).append(surrogate_arrays)
                pooled_hashes.setdefault(key, []).append(str(square["square_hash"]))

    for (mode, core_mode, amplitude, field_arm), pieces in pooled.items():
        combined = {
            key: np.concatenate([piece[key] for piece in pieces])
            for key in ("x", "tau", "tau_stable", "epsilon_r", "latent_epsilon")
        }
        scale = global_sigma_a if mode == PRIMARY_GLOBAL else task_sigma_a
        reported_sigma_g = (
            sigma_g if field_arm == "extensive_balanced_bounded" else surrogate_sigma_g
        )
        texture_rows.append(texture_metric_row(
            ids=ids,
            scale_mode=mode,
            core_mode=core_mode,
            field_arm=field_arm,
            amplitude=amplitude,
            stratum="pooled_candidate_supported_orders",
            background="reference",
            square_hash=latent.json_hash(sorted(pooled_hashes[(mode, core_mode, amplitude, field_arm)])),
            arrays=combined,
            scale_sigma_a=scale,
            sigma_g=reported_sigma_g,
            beta=beta,
            config=config,
            pooled=True,
            permutation_blocks=[
                (np.asarray(piece["x"]), np.asarray(piece["tau"])) for piece in pieces
            ],
        ))

    centered_reference_abs = abs(float(score_candidate(reference[None, :])[0]))
    finite = bool(
        np.isfinite(matrix).all()
        and all(bool(row_value["all_finite"]) for row_value in edge_rows)
        and all(bool(row_value["all_finite"]) for row_value in texture_rows)
    )
    technical = {
        **ids,
        "task_signature": signature,
        "candidate_fingerprint": candidate.fingerprint,
        "candidate_replay_match": bool(candidate.fingerprint == replay.fingerprint),
        "candidate_reference_center_abs": centered_reference_abs,
        "formula_parity_max_abs": max_formula_parity,
        "pair_design_replay_match": pair_replay,
        "square_design_replay_match": square_replay,
        "pair_design_hash": latent.json_hash([value[-1] for value in pair_panels]),
        "square_design_hash": latent.json_hash(sorted(square_hashes)),
        "surrogate_latent_epsilon_max_abs": surrogate_latent_epsilon_max,
        "all_outputs_finite": finite,
    }
    result = {"technical": technical, "edge_rows": edge_rows, "texture_rows": texture_rows}
    latent.write_json_atomic(checkpoint, result)
    return result


def parallel_execute(
    function: Callable[[Mapping[str, Any]], Any],
    payloads: Sequence[Mapping[str, Any]],
    *,
    workers: int,
    label: str,
) -> list[Any]:
    results: list[Any] = []
    total = len(payloads)
    if int(workers) <= 1:
        for index, payload in enumerate(payloads, start=1):
            results.append(function(payload))
            print(f"[EXP-062] {label} {index}/{total}", flush=True)
        return results
    with ProcessPoolExecutor(max_workers=int(workers)) as executor:
        futures = {executor.submit(function, payload): payload for payload in payloads}
        for index, future in enumerate(as_completed(futures), start=1):
            results.append(future.result())
            print(f"[EXP-062] {label} {index}/{total}", flush=True)
    return results


def calibration_selection_rows(
    calibration_rows: Sequence[Mapping[str, Any]], config: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    frame = pd.DataFrame(calibration_rows)
    minimum = float(config["acceptance"]["core"]["minimum_cell_pass_fraction_per_length_bin"])
    summaries: list[dict[str, Any]] = []
    for (beta, mu, occupancy), cell in frame.groupby(
        ["beta", "mu", "reference_occupancy"], sort=True
    ):
        fractions = cell.groupby("length_bin")["core_pass"].mean().sort_index()
        summaries.append({
            "beta": float(beta),
            "mu": float(mu),
            "reference_occupancy": float(occupancy),
            "n_cells": int(len(cell)),
            "n_length_bins": int(len(fractions)),
            "minimum_length_bin_pass_fraction": float(fractions.min()),
            "per_length_bin_pass_fraction": json.dumps(
                {str(key): float(value) for key, value in fractions.items()}, sort_keys=True
            ),
            "core_pass": bool(len(fractions) > 0 and float(fractions.min()) >= minimum),
            "selection_rule": "largest eligible beta then largest eligible mu; GNK outcomes unavailable",
            "gnk_outcomes_used": False,
        })
    beta_order = [float(value) for value in config["calibration"]["beta_grid"]]
    mu_order = [
        float(math.log(float(value) / (1.0 - float(value))))
        for value in config["calibration"]["reference_occupancies"]
    ]
    selected = select_core_calibration(summaries, beta_order=beta_order, mu_order=mu_order)
    for row in summaries:
        row["selected"] = bool(
            selected is not None
            and float(row["beta"]) == float(selected["beta"])
            and float(row["mu"]) == float(selected["mu"])
        )
    return summaries, selected


def minimum_bin_pass_fraction(frame: pd.DataFrame, pass_column: str) -> tuple[float, dict[str, float]]:
    if frame.empty:
        return 0.0, {}
    fractions = frame.groupby("length_bin")[pass_column].mean().sort_index()
    return float(fractions.min()), {str(key): float(value) for key, value in fractions.items()}


def hierarchical_bootstrap_lower(
    frame: pd.DataFrame,
    *,
    metric: str,
    replicates: int,
    seed: int,
) -> float:
    selected = frame[np.isfinite(frame[metric])].copy().sort_values(
        ["reference_id", "core_seed", "gnk_seed"], kind="mergesort"
    )
    references = sorted(str(value) for value in selected["reference_id"].unique())
    if not references:
        return float("nan")
    matrices: dict[str, np.ndarray] = {}
    for reference_id in references:
        ref = selected[selected["reference_id"] == reference_id]
        cores = sorted(int(value) for value in ref["core_seed"].unique())
        gnks = sorted(int(value) for value in ref["gnk_seed"].unique())
        if len(ref) != len(cores) * len(gnks):
            raise ValueError("crossed bootstrap requires a complete core/GNK design")
        matrix = ref.pivot(
            index="core_seed", columns="gnk_seed", values=metric
        ).reindex(index=cores, columns=gnks).to_numpy(dtype=np.float64)
        if matrix.shape != (len(cores), len(gnks)) or not np.isfinite(matrix).all():
            raise ValueError("crossed bootstrap requires one finite row per core/GNK cell")
        matrices[reference_id] = matrix
    rng = np.random.default_rng(int(seed))
    boot = np.empty(int(replicates), dtype=np.float64)
    for replicate in range(int(replicates)):
        sampled_references = rng.choice(references, size=len(references), replace=True)
        reference_values: list[float] = []
        for reference_id in sampled_references:
            matrix = matrices[str(reference_id)]
            sampled_cores = rng.integers(0, matrix.shape[0], size=matrix.shape[0])
            sampled_gnks = rng.integers(0, matrix.shape[1], size=matrix.shape[1])
            reference_values.append(
                float(np.median(matrix[np.ix_(sampled_cores, sampled_gnks)]))
            )
        boot[replicate] = float(np.median(reference_values))
    return float(np.quantile(boot, 0.025))


def core_confirmation_gate(
    rows: pd.DataFrame, *, mode: str, config: Mapping[str, Any]
) -> tuple[bool, dict[str, Any]]:
    selected = rows[(rows["panel"] == "confirmation") & (rows["scale_mode"] == mode)].copy()
    minimum, fractions = minimum_bin_pass_fraction(selected, "core_pass")
    threshold = float(config["acceptance"]["core"]["minimum_cell_pass_fraction_per_length_bin"])
    return bool(fractions and minimum >= threshold), {
        "mode": mode,
        "minimum_length_bin_pass_fraction": minimum,
        "per_length_bin_pass_fraction": fractions,
        "threshold": threshold,
        "n_cells": int(len(selected)),
    }


def transmission_gate(
    edge: pd.DataFrame,
    texture: pd.DataFrame,
    *,
    mode: str,
    core_mode: str,
    config: Mapping[str, Any],
) -> tuple[bool, dict[str, Any]]:
    gates = config["acceptance"]["transmission"]
    edges = edge[
        (edge["scale_mode"] == mode)
        & (edge["sample_id"] == "reference_shell_h1_complete")
        & np.isclose(edge["amplitude"], 1.0)
    ].copy()
    edges["gate_pass"] = (
        (edges["direct_practical_resolution_fraction"] >= float(gates["edge_practical_fraction_min"]))
        & (edges["direct_correlation"] >= float(gates["correlation_min"]))
    )
    edge_minimum, edge_fractions = minimum_bin_pass_fraction(edges, "gate_pass")
    if core_mode == "frozen_additive_core":
        edge_minimum = 1.0
        edge_fractions = {"not_required_for_frozen_texture_control": 1.0}

    selected = texture[
        (texture["scale_mode"] == mode)
        & (texture["additive_core_mode"] == core_mode)
        & (texture["field_arm"] == "extensive_balanced_bounded")
        & texture["pooled_supported_reference"].eq(True)
        & np.isclose(texture["amplitude"], 1.0)
    ].copy()
    if selected.empty:
        return False, {"mode": mode, "core_mode": core_mode, "rows": 0}
    ratio_minimum = float(gates["actual_to_frozen_rms_ratio_min"])
    if core_mode == "actual_additive_core":
        frozen = texture[
            (texture["scale_mode"] == mode)
            & (texture["additive_core_mode"] == "frozen_additive_core")
            & (texture["field_arm"] == "extensive_balanced_bounded")
            & texture["pooled_supported_reference"].eq(True)
            & np.isclose(texture["amplitude"], 1.0)
        ][["reference_id", "core_seed", "gnk_seed", "direct_tau_y_rms"]].rename(
            columns={"direct_tau_y_rms": "frozen_y_rms"}
        )
        selected = selected.merge(
            frozen,
            on=["reference_id", "core_seed", "gnk_seed"],
            how="left",
            validate="one_to_one",
        )
        selected["actual_to_frozen_rms_ratio"] = (
            selected["direct_tau_y_rms"] / selected["frozen_y_rms"].replace(0.0, np.nan)
        )
        ratio_pass = selected["actual_to_frozen_rms_ratio"] >= ratio_minimum
    else:
        selected["actual_to_frozen_rms_ratio"] = np.nan
        ratio_pass = pd.Series(True, index=selected.index)
    selected["texture_gate_pass"] = (
        (selected["direct_tau_practical_resolution_fraction"] >= float(gates["texture_practical_fraction_min"]))
        & (selected["direct_tau_projection"] > 0.0)
        & (selected["direct_tau_correlation"] >= float(gates["correlation_min"]))
        & selected["permutation_pass"].eq(True)
        & ratio_pass
    )
    texture_minimum, texture_fractions = minimum_bin_pass_fraction(selected, "texture_gate_pass")
    if core_mode == "actual_additive_core":
        edge_cells = edges[[
            "reference_id", "length_bin", "core_seed", "gnk_seed", "gate_pass"
        ]].rename(columns={"gate_pass": "edge_gate_pass"})
        joint = selected.merge(
            edge_cells,
            on=["reference_id", "length_bin", "core_seed", "gnk_seed"],
            how="left",
            validate="one_to_one",
        )
        joint["joint_gate_pass"] = (
            joint["texture_gate_pass"].eq(True) & joint["edge_gate_pass"].eq(True)
        )
    else:
        joint = selected.copy()
        joint["joint_gate_pass"] = joint["texture_gate_pass"].eq(True)
    joint_minimum, joint_fractions = minimum_bin_pass_fraction(joint, "joint_gate_pass")
    bootstrap_lower = hierarchical_bootstrap_lower(
        selected,
        metric="direct_tau_projection",
        replicates=int(config["acceptance"]["bootstrap_replicates"]),
        seed=int(config["acceptance"]["bootstrap_seed"])
        + (0 if mode == PRIMARY_GLOBAL else 100)
        + (0 if core_mode == "actual_additive_core" else 10),
    )
    cell_fraction_minimum = float(gates["minimum_cell_pass_fraction_per_length_bin"])
    passed = bool(
        edge_fractions
        and texture_fractions
        and joint_fractions
        and edge_minimum >= cell_fraction_minimum
        and texture_minimum >= cell_fraction_minimum
        and joint_minimum >= cell_fraction_minimum
        and np.isfinite(bootstrap_lower)
        and bootstrap_lower > float(gates["projection_bootstrap_lower_min"])
    )
    return passed, {
        "mode": mode,
        "core_mode": core_mode,
        "edge_rows": int(len(edges)),
        "texture_rows": int(len(selected)),
        "edge_minimum_length_bin_pass_fraction": edge_minimum,
        "edge_per_length_bin_pass_fraction": edge_fractions,
        "texture_minimum_length_bin_pass_fraction": texture_minimum,
        "texture_per_length_bin_pass_fraction": texture_fractions,
        "joint_minimum_length_bin_pass_fraction": joint_minimum,
        "joint_per_length_bin_pass_fraction": joint_fractions,
        "projection_hierarchical_bootstrap_lower": bootstrap_lower,
        "minimum_actual_to_frozen_rms_ratio": (
            float(selected["actual_to_frozen_rms_ratio"].min())
            if core_mode == "actual_additive_core" else None
        ),
    }


def rank_correlation(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=np.float64)
    if array.size < 2 or not np.isfinite(array).all():
        return float("nan")
    ranks = pd.Series(array).rank(method="average").to_numpy(dtype=float)
    if float(np.std(ranks)) == 0.0:
        return 0.0
    return float(np.corrcoef(np.arange(array.size, dtype=float), ranks)[0, 1])


def dose_response_gate(
    texture: pd.DataFrame,
    *,
    mode: str,
    config: Mapping[str, Any],
) -> tuple[bool, list[dict[str, Any]], dict[str, Any]]:
    gates = config["acceptance"]["dose_response"]
    doses = [float(value) for value in config["intervention"]["directionality_amplitudes"]]
    selected = texture[
        (texture["scale_mode"] == mode)
        & (texture["additive_core_mode"] == "actual_additive_core")
        & (texture["field_arm"] == "extensive_balanced_bounded")
        & texture["pooled_supported_reference"].eq(True)
        & texture["amplitude"].isin(doses)
    ].copy()
    rows: list[dict[str, Any]] = []

    def summarize_curve(group: pd.DataFrame) -> dict[str, Any]:
        medians = group.groupby("amplitude")[[
            "direct_tau_y_rms", "direct_tau_practical_resolution_fraction"
        ]].median().reindex(doses)
        if medians.notna().all().all():
            rms_values = medians["direct_tau_y_rms"].to_numpy(dtype=float)
            practical_values = medians["direct_tau_practical_resolution_fraction"].to_numpy(dtype=float)
            rms_rho = rank_correlation(rms_values)
            practical_rho = rank_correlation(practical_values)
            rms_delta = float(rms_values[-1] - rms_values[0])
            practical_delta = float(practical_values[-1] - practical_values[0])
        else:
            rms_values = np.full(len(doses), np.nan)
            practical_values = np.full(len(doses), np.nan)
            rms_rho = practical_rho = rms_delta = practical_delta = float("nan")
        passed = bool(
            np.isfinite(rms_rho)
            and np.isfinite(practical_rho)
            and rms_rho >= float(gates["spearman_min"])
            and practical_rho >= float(gates["spearman_min"])
            and rms_delta >= float(gates["rms_growth_floor"])
            and practical_delta >= float(gates["practical_fraction_growth_floor"])
        )
        return {
            "amplitudes": json.dumps(doses),
            "tau_rms_medians": json.dumps([float(value) for value in rms_values]),
            "practical_fraction_medians": json.dumps([float(value) for value in practical_values]),
            "tau_rms_spearman": rms_rho,
            "practical_fraction_spearman": practical_rho,
            "tau_rms_a2_minus_a0.25": rms_delta,
            "practical_fraction_a2_minus_a0.25": practical_delta,
            "directionality_pass": passed,
        }

    for (reference_id, length_bin, length, core_seed, gnk_seed), group in selected.groupby(
        ["reference_id", "length_bin", "length", "core_seed", "gnk_seed"], sort=True
    ):
        rows.append({
            "aggregation_unit": "protein_core_gnk_cell",
            "scale_mode": mode,
            "reference_id": reference_id,
            "length_bin": length_bin,
            "length": int(length),
            "core_seed": int(core_seed),
            "gnk_seed": int(gnk_seed),
            **summarize_curve(group),
        })
    for (reference_id, length_bin, length), group in selected.groupby(
        ["reference_id", "length_bin", "length"], sort=True
    ):
        rows.append({
            "aggregation_unit": "protein_median_over_cells",
            "scale_mode": mode,
            "reference_id": reference_id,
            "length_bin": length_bin,
            "length": int(length),
            "core_seed": None,
            "gnk_seed": None,
            **summarize_curve(group),
        })
    frame = pd.DataFrame(rows)
    protein = frame[frame["aggregation_unit"] == "protein_median_over_cells"].copy()
    cells = frame[frame["aggregation_unit"] == "protein_core_gnk_cell"].copy()
    if protein.empty or cells.empty:
        return False, rows, {"mode": mode, "protein_count": 0, "cell_count": 0}
    protein_pass_fraction = float(protein["directionality_pass"].mean())
    protein_bin_min, protein_bin_fractions = minimum_bin_pass_fraction(
        protein, "directionality_pass"
    )
    cell_bin_min, cell_bin_fractions = minimum_bin_pass_fraction(
        cells, "directionality_pass"
    )
    rms_lower = simple_bootstrap_lower(
        protein["tau_rms_a2_minus_a0.25"].to_numpy(dtype=float),
        replicates=int(config["acceptance"]["bootstrap_replicates"]),
        seed=int(config["acceptance"]["bootstrap_seed"]) + (200 if mode == PRIMARY_GLOBAL else 300),
    )
    practical_lower = simple_bootstrap_lower(
        protein["practical_fraction_a2_minus_a0.25"].to_numpy(dtype=float),
        replicates=int(config["acceptance"]["bootstrap_replicates"]),
        seed=int(config["acceptance"]["bootstrap_seed"]) + (201 if mode == PRIMARY_GLOBAL else 301),
    )
    passed = bool(
        protein_pass_fraction >= float(gates["minimum_protein_pass_fraction"])
        and protein_bin_min >= float(gates["minimum_protein_pass_fraction_per_length_bin"])
        and cell_bin_min >= float(gates["minimum_cell_pass_fraction_per_length_bin"])
        and rms_lower > 0.0
        and practical_lower > 0.0
    )
    return passed, rows, {
        "mode": mode,
        "protein_count": int(len(protein)),
        "cell_count": int(len(cells)),
        "protein_pass_fraction": protein_pass_fraction,
        "minimum_protein_pass_fraction_per_length_bin": protein_bin_min,
        "protein_pass_fraction_per_length_bin": protein_bin_fractions,
        "minimum_cell_pass_fraction_per_length_bin": cell_bin_min,
        "cell_pass_fraction_per_length_bin": cell_bin_fractions,
        "tau_rms_growth_bootstrap_lower": rms_lower,
        "practical_fraction_growth_bootstrap_lower": practical_lower,
    }


def simple_bootstrap_lower(values: np.ndarray, *, replicates: int, seed: int) -> float:
    array = np.asarray(values, dtype=np.float64)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return float("nan")
    rng = np.random.default_rng(int(seed))
    boot = np.empty(int(replicates), dtype=np.float64)
    for index in range(int(replicates)):
        boot[index] = float(np.median(rng.choice(array, size=array.size, replace=True)))
    return float(np.quantile(boot, 0.025))


def length_stability_gate(
    texture: pd.DataFrame,
    *,
    mode: str,
    config: Mapping[str, Any],
) -> tuple[bool, dict[str, Any]]:
    gates = config["acceptance"]["length_stability"]
    selected = texture[
        (texture["scale_mode"] == mode)
        & (texture["additive_core_mode"] == "actual_additive_core")
        & (texture["field_arm"] == "extensive_balanced_bounded")
        & texture["pooled_supported_reference"].eq(True)
        & np.isclose(texture["amplitude"], 1.0)
    ]
    protein = selected.groupby(
        ["reference_id", "length_bin", "length"], as_index=False
    )["direct_tau_y_rms"].median()
    protein = protein[np.isfinite(protein["direct_tau_y_rms"]) & (protein["direct_tau_y_rms"] > 0.0)]
    if len(protein) < 3:
        return False, {
            "fit_id": f"{mode}_actual_tau_a1",
            "n_proteins": int(len(protein)),
            "slope": float("nan"),
            "slope_ci_low": float("nan"),
            "slope_ci_high": float("nan"),
        }
    x = np.log(protein["length"].to_numpy(dtype=float))
    y = np.log(protein["direct_tau_y_rms"].to_numpy(dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    rng = np.random.default_rng(
        int(config["acceptance"]["bootstrap_seed"]) + (400 if mode == PRIMARY_GLOBAL else 500)
    )
    boot: list[float] = []
    bins = sorted(str(value) for value in protein["length_bin"].unique())
    for _ in range(int(config["acceptance"]["bootstrap_replicates"])):
        pieces = []
        for length_bin in bins:
            block = protein[protein["length_bin"] == length_bin]
            indices = rng.integers(0, len(block), size=len(block))
            pieces.append(block.iloc[indices])
        replicate = pd.concat(pieces, ignore_index=True)
        bx = np.log(replicate["length"].to_numpy(dtype=float))
        by = np.log(replicate["direct_tau_y_rms"].to_numpy(dtype=float))
        if float(np.std(bx)) > 0.0:
            boot.append(float(np.polyfit(bx, by, 1)[0]))
    low = float(np.quantile(boot, 0.025)) if boot else float("nan")
    high = float(np.quantile(boot, 0.975)) if boot else float("nan")
    passed = bool(
        len(set(protein["length_bin"])) == int(config["acceptance"]["expected_length_bins"])
        and abs(float(slope)) <= float(gates["slope_abs_max"])
        and low >= -float(gates["slope_ci_abs_max"])
        and high <= float(gates["slope_ci_abs_max"])
    )
    return passed, {
        "fit_id": f"{mode}_actual_tau_a1",
        "metric": "protein_median_direct_tau_y_rms",
        "n_proteins": int(len(protein)),
        "n_length_bins": int(len(set(protein["length_bin"]))),
        "slope": float(slope),
        "intercept": float(intercept),
        "slope_ci_low": low,
        "slope_ci_high": high,
        "bootstrap_replicates_realized": int(len(boot)),
    }


def acceptance_checks(
    *,
    references: Sequence[Mapping[str, Any]],
    scale_rows: Sequence[Mapping[str, Any]],
    selection: Mapping[str, Any] | None,
    core_rows: pd.DataFrame,
    technical_rows: pd.DataFrame,
    edge: pd.DataFrame,
    texture: pd.DataFrame,
    config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    checks: list[dict[str, Any]] = []

    def add(
        check_id: str,
        passed: bool,
        value: Any,
        criterion: str,
        *,
        family: str,
        required: bool = True,
    ) -> None:
        checks.append({
            "check_id": check_id,
            "family": family,
            "required": bool(required),
            "passed": bool(passed),
            "value": json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value,
            "criterion": criterion,
        })

    calibration_hashes = {
        str(row["sequence_sha256"]) for row in references if row["split"] == "calibration"
    }
    confirmation_hashes = {
        str(row["sequence_sha256"]) for row in references if row["split"] == "confirmation"
    }
    split_pass = bool(calibration_hashes and confirmation_hashes and calibration_hashes.isdisjoint(confirmation_hashes))
    add(
        "reference_split_disjoint",
        split_pass,
        {"calibration": len(calibration_hashes), "confirmation": len(confirmation_hashes)},
        "calibration and confirmation contain distinct sequence hashes",
        family="technical",
    )
    scale_finite = bool(scale_rows and all(bool(row["finite"]) for row in scale_rows))
    add(
        "calibration_scale_finite",
        scale_finite,
        {"rows": len(scale_rows)},
        "every complete-h1 calibration scale is finite and positive",
        family="technical",
    )

    heldout_core_failed = False
    if selection is not None and not core_rows.empty and "panel" in core_rows.columns:
        confirmation_preview = core_rows[core_rows["panel"] == "confirmation"]
        if not confirmation_preview.empty:
            heldout_core_failed = not core_confirmation_gate(
                core_rows, mode=PRIMARY_GLOBAL, config=config
            )[0]
    if technical_rows.empty:
        task_technical = bool(selection is None or heldout_core_failed)
        technical_values: dict[str, Any] = {
            "validation_tasks": 0,
            "stage_skipped_no_selection": selection is None,
            "stage_skipped_core_non_generalization": heldout_core_failed,
        }
    else:
        tolerance = config["acceptance"]["technical"]
        zero_edge = edge[
            edge["scale_mode"].isin([PRIMARY_GLOBAL, EXACT_BETA])
            & np.isclose(edge["amplitude"], 0.0)
        ]
        zero_texture = texture[
            texture["field_arm"].eq("extensive_balanced_bounded")
            & np.isclose(texture["amplitude"], 0.0)
        ]
        zero_max = max(
            float(zero_edge["direct_y_rms"].max()) if not zero_edge.empty else float("inf"),
            float(zero_texture["direct_tau_y_rms"].max()) if not zero_texture.empty else float("inf"),
        )
        technical_values = {
            "validation_tasks": int(len(technical_rows)),
            "candidate_replay_all": bool(technical_rows["candidate_replay_match"].all()),
            "pair_replay_all": bool(technical_rows["pair_design_replay_match"].all()),
            "square_replay_all": bool(technical_rows["square_design_replay_match"].all()),
            "finite_all": bool(technical_rows["all_outputs_finite"].all()),
            "max_reference_center_abs": float(technical_rows["candidate_reference_center_abs"].max()),
            "max_formula_parity_abs": float(technical_rows["formula_parity_max_abs"].max()),
            "max_surrogate_latent_epsilon_abs": float(
                technical_rows["surrogate_latent_epsilon_max_abs"].max()
            ),
            "max_zero_amplitude_residual_rms": zero_max,
        }
        task_technical = bool(
            technical_values["candidate_replay_all"]
            and technical_values["pair_replay_all"]
            and technical_values["square_replay_all"]
            and technical_values["finite_all"]
            and technical_values["max_reference_center_abs"]
            <= float(tolerance["reference_center_atol"])
            and technical_values["max_formula_parity_abs"]
            <= float(tolerance["formula_parity_atol"])
            and technical_values["max_surrogate_latent_epsilon_abs"]
            <= float(tolerance["additive_surrogate_epistasis_atol"])
            and technical_values["max_zero_amplitude_residual_rms"]
            <= float(tolerance["zero_amplitude_atol"])
        )
    add(
        "validation_execution_and_algebra",
        task_technical,
        technical_values,
        "deterministic designs, finite arrays, centering, formula parity, additive null, and a=0 exactness all pass",
        family="technical",
    )
    technical_pass = bool(split_pass and scale_finite and task_technical)
    add(
        "core_operating_point_selected",
        selection is not None,
        None if selection is None else {
            "beta": float(selection["beta"]),
            "mu": float(selection["mu"]),
            "reference_occupancy": float(selection["reference_occupancy"]),
        },
        "at least one core-only beta/mu cell passes every per-length calibration gate",
        family="calibration",
    )

    if selection is None or core_rows.empty:
        global_core_pass = exact_core_pass = False
        global_core_values = exact_core_values = {"stage_skipped": True}
        global_transmission = frozen_texture = exact_transmission = False
        global_transmission_values = frozen_values = exact_transmission_values = {"stage_skipped": True}
        global_direction = exact_direction = False
        global_dose_rows: list[dict[str, Any]] = []
        exact_dose_rows: list[dict[str, Any]] = []
        global_direction_values = exact_direction_values = {"stage_skipped": True}
        global_length = exact_length = False
        global_fit = {"fit_id": f"{PRIMARY_GLOBAL}_actual_tau_a1", "stage_skipped": True}
        exact_fit = {"fit_id": f"{EXACT_BETA}_actual_tau_a1", "stage_skipped": True}
    else:
        global_core_pass, global_core_values = core_confirmation_gate(
            core_rows, mode=PRIMARY_GLOBAL, config=config
        )
        exact_core_pass, exact_core_values = core_confirmation_gate(
            core_rows, mode=EXACT_BETA, config=config
        )
        if not global_core_pass:
            skipped = {"stage_skipped_core_non_generalization": True}
            global_transmission = frozen_texture = exact_transmission = False
            global_transmission_values = dict(skipped)
            frozen_values = dict(skipped)
            exact_transmission_values = dict(skipped)
            global_direction = exact_direction = False
            global_dose_rows = []
            exact_dose_rows = []
            global_direction_values = dict(skipped)
            exact_direction_values = dict(skipped)
            global_length = exact_length = False
            global_fit = {
                "fit_id": f"{PRIMARY_GLOBAL}_actual_tau_a1", **skipped
            }
            exact_fit = {"fit_id": f"{EXACT_BETA}_actual_tau_a1", **skipped}
        else:
            global_transmission, global_transmission_values = transmission_gate(
                edge, texture, mode=PRIMARY_GLOBAL,
                core_mode="actual_additive_core", config=config,
            )
            frozen_texture, frozen_values = transmission_gate(
                edge, texture, mode=PRIMARY_GLOBAL,
                core_mode="frozen_additive_core", config=config,
            )
            exact_transmission, exact_transmission_values = transmission_gate(
                edge, texture, mode=EXACT_BETA,
                core_mode="actual_additive_core", config=config,
            )
            global_direction, global_dose_rows, global_direction_values = dose_response_gate(
                texture, mode=PRIMARY_GLOBAL, config=config
            )
            exact_direction, exact_dose_rows, exact_direction_values = dose_response_gate(
                texture, mode=EXACT_BETA, config=config
            )
            global_length, global_fit = length_stability_gate(
                texture, mode=PRIMARY_GLOBAL, config=config
            )
            exact_length, exact_fit = length_stability_gate(
                texture, mode=EXACT_BETA, config=config
            )

    add(
        "global_core_calibration_heldout",
        global_core_pass,
        global_core_values,
        "the frozen global core calibration independently passes in each confirmation length bin",
        family="confirmation_core",
    )
    add(
        "global_actual_core_gnk_transmission",
        global_transmission,
        global_transmission_values,
        "h1 edges and candidate-additivized actual-core texture pass practical, correlation, projection, permutation, and ratio gates",
        family="transmission",
    )
    add(
        "global_actual_core_gnk_directionality",
        global_direction,
        global_direction_values,
        "protein-level tau RMS and practical fraction increase monotonically over a=0.25..2 with positive paired growth",
        family="directionality",
    )
    add(
        "global_actual_core_gnk_length_stability",
        global_length,
        global_fit,
        "protein-level a=1 tau RMS log-length slope and bootstrap interval lie inside preregistered bounds",
        family="length_stability",
    )
    add(
        "frozen_core_texture_positive_control",
        frozen_texture,
        frozen_values,
        "the frozen-additive-core candidate texture positive control passes",
        family="control",
    )
    exact_full = bool(
        exact_core_pass
        and exact_transmission
        and exact_direction
        and exact_length
        and frozen_texture
    )
    add(
        "exact_per_landscape_beta_diagnostic",
        exact_full,
        {
            "core": exact_core_values,
            "transmission": exact_transmission_values,
            "directionality": exact_direction_values,
            "length_stability": exact_fit,
        },
        "the diagnostic exact-beta arm passes all core and GNK gates",
        family="diagnostic",
        required=False,
    )

    ablation_values: dict[str, Any] = {}
    ablation_complete = not edge.empty
    if not edge.empty:
        def finite_median(values: pd.Series) -> float | None:
            array = values.to_numpy(dtype=float)
            array = array[np.isfinite(array)]
            return float(np.median(array)) if array.size else None

        for mode in (STEEPNESS_ONLY, THRESHOLD_ONLY, RAW_MARGIN8, RAW_MARGIN10):
            selected_edges = edge[
                (edge["scale_mode"] == mode)
                & (edge["sample_id"] == "reference_shell_h1_complete")
            ]
            if selected_edges.empty:
                ablation_complete = False
            ablation_values[mode] = {
                "n_rows": int(len(selected_edges)),
                "median_practical_fraction": (
                    finite_median(selected_edges["direct_practical_resolution_fraction"])
                    if not selected_edges.empty else None
                ),
                "median_correlation": (
                    finite_median(selected_edges["direct_correlation"])
                    if not selected_edges.empty else None
                ),
            }
    add(
        "steepness_threshold_margin_ablation_reported",
        ablation_complete and len(ablation_values) == 4,
        ablation_values,
        "steepness-only, threshold-only, raw margin-8, and constructor margin-10 controls are reported",
        family="diagnostic",
        required=False,
    )
    flags = {
        "technical_pass": technical_pass,
        "selection_exists": selection is not None,
        "confirmation_core_pass": global_core_pass,
        "global_transmission_pass": bool(global_transmission and frozen_texture),
        "global_directionality_pass": global_direction,
        "global_length_stability_pass": global_length,
        "exact_full_pass": exact_full,
        "frozen_texture_pass": frozen_texture,
        "exact_core_pass": exact_core_pass,
        "exact_transmission_pass": exact_transmission,
        "exact_directionality_pass": exact_direction,
        "exact_length_stability_pass": exact_length,
    }
    return checks, flags, global_dose_rows + exact_dose_rows, [global_fit, exact_fit]


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, quick=bool(args.quick))
    validate_config(config)
    references, catalog_hash = load_references(project_root, config)
    implementation_hashes = source_hashes(project_root, config)
    implementation_hash = latent.json_hash(implementation_hashes)
    workers = max(1, int(os.environ.get(
        "LABPROJ_WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "1")
    )))
    calibration_references = [row for row in references if row["split"] == "calibration"]
    confirmation_references = [row for row in references if row["split"] == "confirmation"]
    calibration_seeds = [int(value) for value in config["replication"]["calibration_core_seeds"]]
    confirmation_core_seeds = [
        int(value) for value in config["replication"]["confirmation_core_seeds"]
    ]
    confirmation_gnk_seeds = [
        int(value) for value in config["replication"]["confirmation_gnk_seeds"]
    ]
    preflight = {
        "experiment_id": "EXP-062",
        "quick": bool(args.quick),
        "scientific_decision_enabled": not bool(args.quick),
        "config_file_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "effective_config_hash": latent.json_hash(config),
        "reference_catalog_sha256": catalog_hash,
        "implementation_hash": implementation_hash,
        "implementation_hash_components": implementation_hashes,
        "reference_count": len(references),
        "calibration_reference_count": len(calibration_references),
        "confirmation_reference_count": len(confirmation_references),
        "length_bins": sorted(set(row["length_bin"] for row in references)),
        "calibration_core_seed_count": len(calibration_seeds),
        "confirmation_core_seed_count": len(confirmation_core_seeds),
        "confirmation_gnk_seed_count": len(confirmation_gnk_seeds),
        "confirmation_design": "crossed additive-core and GNK seeds nested within held-out protein",
        "gnk_outcomes_available_during_selection": False,
        "observation_noise_enabled": False,
        "candidate_is_experiment_local": True,
        "public_model_changed": False,
        "empirical_fitness_realism_tested": False,
        "optimizer_or_accessibility_tested": False,
        "stage_4c_4d_4e_hard_stop_retained": True,
    }
    latent.write_json_atomic(output_dir / "preflight.json", preflight)
    latent.write_csv_atomic(
        output_dir / "reference_catalog.csv",
        [{key: value for key, value in row.items() if key != "sequence"} for row in references],
    )

    scale_payloads = [
        {"reference": reference, "core_seed": seed, "config": config}
        for reference in calibration_references for seed in calibration_seeds
    ]
    scale_rows = parallel_execute(
        core_scale_task, scale_payloads, workers=workers, label="calibration-scale"
    )
    scale_rows.sort(key=lambda row: (
        str(row["length_bin"]), str(row["reference_id"]), int(row["core_seed"])
    ))
    global_sigma_a = float(np.median([
        float(row["sigma_a_complete_h1"]) for row in scale_rows
    ]))
    calibration_payloads = [
        {
            "reference": reference,
            "core_seed": seed,
            "config": config,
            "global_sigma_a": global_sigma_a,
        }
        for reference in calibration_references for seed in calibration_seeds
    ]
    calibration_nested = parallel_execute(
        calibration_task, calibration_payloads, workers=workers, label="core-calibration"
    )
    calibration_rows = [row for rows in calibration_nested for row in rows]
    calibration_rows.sort(key=lambda row: (
        float(row["beta"]), float(row["mu"]), str(row["length_bin"]),
        str(row["reference_id"]), int(row["core_seed"]),
    ))
    selection_rows, selection = calibration_selection_rows(calibration_rows, config)
    latent.write_csv_atomic(output_dir / "calibration_selection.csv", selection_rows)

    confirmation_core_rows: list[dict[str, Any]] = []
    validation_results: list[dict[str, Any]] = []
    if selection is not None:
        beta = float(selection["beta"])
        mu = float(selection["mu"])
        confirmation_core_payloads = [
            {
                "reference": reference,
                "core_seed": seed,
                "config": config,
                "global_sigma_a": global_sigma_a,
                "beta": beta,
                "mu": mu,
            }
            for reference in confirmation_references for seed in confirmation_core_seeds
        ]
        confirmation_nested = parallel_execute(
            confirmation_core_task,
            confirmation_core_payloads,
            workers=workers,
            label="heldout-core",
        )
        confirmation_core_rows = [row for rows in confirmation_nested for row in rows]
        confirmation_core_rows.sort(key=lambda row: (
            str(row["scale_mode"]), str(row["length_bin"]),
            str(row["reference_id"]), int(row["core_seed"]),
        ))
        heldout_core_ready, heldout_core_values = core_confirmation_gate(
            pd.DataFrame(confirmation_core_rows), mode=PRIMARY_GLOBAL, config=config
        )
        if heldout_core_ready:
            checkpoint_dir = output_dir / "checkpoints"
            checkpoint_dir.mkdir(parents=True, exist_ok=True)
            validation_payloads: list[dict[str, Any]] = []
            for reference in confirmation_references:
                for core_seed in confirmation_core_seeds:
                    for gnk_seed in confirmation_gnk_seeds:
                        signature = latent.json_hash({
                            "reference": reference,
                            "core_seed": core_seed,
                            "gnk_seed": gnk_seed,
                            "global_sigma_a": global_sigma_a,
                            "beta": beta,
                            "mu": mu,
                            "effective_config_hash": preflight["effective_config_hash"],
                            "implementation_hash": implementation_hash,
                        })
                        validation_payloads.append({
                            "reference": reference,
                            "core_seed": core_seed,
                            "gnk_seed": gnk_seed,
                            "config": config,
                            "global_sigma_a": global_sigma_a,
                            "beta": beta,
                            "mu": mu,
                            "task_signature": signature,
                            "checkpoint_path": str(
                                checkpoint_dir
                                / f"{latent.safe_name(reference['reference_id'])}__core-{core_seed}__gnk-{gnk_seed}__{signature[:16]}.json"
                            ),
                        })
            validation_results = parallel_execute(
                validation_task,
                validation_payloads,
                workers=workers,
                label="heldout-gnk",
            )
            validation_results.sort(key=lambda result: (
                str(result["technical"]["length_bin"]),
                str(result["technical"]["reference_id"]),
                int(result["technical"]["core_seed"]),
                int(result["technical"]["gnk_seed"]),
            ))
        else:
            print(
                "[EXP-062] heldout GNK stage skipped because global core calibration "
                f"did not generalize: {json.dumps(heldout_core_values, sort_keys=True)}",
                flush=True,
            )

    core_output_rows: list[dict[str, Any]] = []
    for row in scale_rows:
        core_output_rows.append({"row_type": "scale_estimate", **row})
    for row in calibration_rows:
        core_output_rows.append({"row_type": "core_grid_cell", **row})
    for row in confirmation_core_rows:
        core_output_rows.append({"row_type": "heldout_core_confirmation", **row})
    latent.write_csv_atomic(output_dir / "core_calibration_metrics.csv", core_output_rows)

    technical_rows = [result["technical"] for result in validation_results]
    edge_rows = [row for result in validation_results for row in result["edge_rows"]]
    texture_rows = [row for result in validation_results for row in result["texture_rows"]]
    edge_rows.sort(key=lambda row: (
        str(row["length_bin"]), str(row["reference_id"]), int(row["core_seed"]),
        int(row["gnk_seed"]), str(row["scale_mode"]), str(row["sample_id"]),
        float(row["amplitude"]),
    ))
    texture_rows.sort(key=lambda row: (
        str(row["length_bin"]), str(row["reference_id"]), int(row["core_seed"]),
        int(row["gnk_seed"]), str(row["scale_mode"]), str(row["additive_core_mode"]),
        str(row["field_arm"]), str(row["square_background"]),
        str(row["square_stratum"]), float(row["amplitude"]),
    ))
    latent.write_csv_atomic(output_dir / "validation_edge_metrics.csv", edge_rows)
    latent.write_csv_atomic(output_dir / "validation_texture_metrics.csv", texture_rows)
    core_frame = pd.DataFrame(calibration_rows + confirmation_core_rows)
    technical_frame = pd.DataFrame(technical_rows)
    edge_frame = pd.DataFrame(edge_rows)
    texture_frame = pd.DataFrame(texture_rows)
    checks, flags, dose_rows, fit_rows = acceptance_checks(
        references=references,
        scale_rows=scale_rows,
        selection=selection,
        core_rows=core_frame,
        technical_rows=technical_frame,
        edge=edge_frame,
        texture=texture_frame,
        config=config,
    )
    latent.write_csv_atomic(output_dir / "dose_response_metrics.csv", dose_rows)
    latent.write_csv_atomic(output_dir / "scale_law_fits.csv", fit_rows)
    latent.write_csv_atomic(output_dir / "acceptance_checks.csv", checks)
    decision = classify_decision(quick=bool(args.quick), **{
        key: bool(flags[key]) for key in (
            "technical_pass", "selection_exists", "confirmation_core_pass",
            "global_transmission_pass", "global_directionality_pass",
            "global_length_stability_pass", "exact_full_pass", "frozen_texture_pass",
        )
    })
    required = [row for row in checks if bool(row["required"])]
    summary = {
        **preflight,
        "global_sigma_a_complete_h1_median": global_sigma_a,
        "selected_calibration": None if selection is None else {
            "beta": float(selection["beta"]),
            "mu": float(selection["mu"]),
            "reference_occupancy": float(selection["reference_occupancy"]),
            "global_raw_steepness": float(selection["beta"]) / global_sigma_a,
        },
        "calibration_scale_task_count": len(scale_rows),
        "calibration_grid_cell_count": len(calibration_rows),
        "confirmation_core_cell_count": len(confirmation_core_rows),
        "confirmation_gnk_task_count": len(validation_results),
        "validation_edge_row_count": len(edge_rows),
        "validation_texture_row_count": len(texture_rows),
        "required_checks_passed": int(sum(bool(row["passed"]) for row in required)),
        "required_check_count": len(required),
        "decision_flags": flags,
        "scientific_pass": (
            None if bool(args.quick)
            else decision == "GLOBAL_OPERATOR_CALIBRATED_GNK_TRANSMISSION_VALIDATED_NOISELESS"
        ),
        "decision": decision,
        "public_model_changed": False,
        "candidate_is_experiment_local": True,
        "empirical_realism_tested": False,
        "optimizer_or_accessibility_tested": False,
    }
    latent.write_json_atomic(output_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
