#!/usr/bin/env python3
"""EXP-061 paired transmission test for GNK texture through the public sigmoid.

The runner scores additive and GNK latent values once, then streams readout
conditions into seed-level sufficient statistics.  It deliberately separates
genotype-wise transfer-function calibration, common-threshold mechanistic
landscapes, and the direct float64 deployment calculation.
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
for _path in (PROJECT_ROOT / "src", PROJECT_ROOT / "external" / "Adapt-Env", PROJECT_ROOT / "scripts"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from adaptenv import FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from project_adapt_env.extensive_gnk import (  # noqa: E402
    ExtensiveGNKField,
    additive_mutation_energy,
    stable_uint64,
)
import extensive_gnk_latent_scaling_exp060 as latent  # noqa: E402


ARM_CANDIDATE = "extensive_balanced_bounded"
ARM_LOCAL = "legacy_support_local_normalized"
ARM_PUBLIC = "legacy_public_global_fixed"
ARM_ADDITIVE = "additive_surrogate"
ARMS = (ARM_CANDIDATE, ARM_LOCAL, ARM_PUBLIC, ARM_ADDITIVE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


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
        "references", "replication", "candidate", "legacy_control", "readout",
        "intervention", "probe", "acceptance", "decision_rule", "outputs",
    }
    missing = sorted(required - set(config))
    if missing:
        raise ValueError(f"config is missing required sections: {missing}")
    if str(config.get("experiment_id")) != "EXP-061":
        raise ValueError("EXP-061 runner requires experiment_id EXP-061")
    readout = config["readout"]
    if float(readout["threshold_steepness"]) <= 0 or float(readout["output_scale"]) <= 0:
        raise ValueError("readout scale and steepness must be positive")
    if 0.5 not in [float(x) for x in readout["mechanistic_reference_occupancies"]]:
        raise ValueError("reference occupancy 0.5 is required")
    if 0.5 not in [float(x) for x in readout["pointwise_occupancies"]]:
        raise ValueError("pointwise occupancy 0.5 is required")
    scopes = readout.get("operating_point_scopes", {})
    if not {"pointwise_occupancy", "reference_anchored", "deployment"}.issubset(scopes):
        raise ValueError("all three operating-point scopes are required")
    probe = config["probe"]
    if not bool(probe.get("complete_reference_single_mutant_panel")):
        raise ValueError("complete reference single-mutant panel is required")
    if int(probe.get("uniform_background_single_mutant_pairs", 0)) <= 0:
        raise ValueError("uniform-background h=1 pairs are required")
    if not bool(probe.get("candidate_additivized_counterfactual", {}).get("enabled")):
        raise ValueError("candidate-additivized plaquette control is required")
    if not {"reference", "uniform_q20"}.issubset(set(probe.get("square_backgrounds", []))):
        raise ValueError("reference and uniform-q20 square backgrounds are required")
    acceptance = config["acceptance"]
    for key in (
        "identity_atol", "zero_weight_atol", "reference_center_atol",
        "operational_detection_threshold", "texture_permutation_alpha",
        "topology_non_tie_agreement_min", "texture_projection_bootstrap_lower_min",
        "texture_permutation_pass_fraction_min",
    ):
        if key not in acceptance:
            raise ValueError(f"acceptance is missing {key}")
    if float(acceptance["topology_non_tie_agreement_min"]) != 1.0:
        raise ValueError("non-tie topology agreement must be exact")
    amplitudes = [float(x) for x in config["intervention"]["dimensionless_amplitudes"]]
    if 0.0 not in amplitudes or 0.5 not in amplitudes:
        raise ValueError("dimensionless amplitudes must include 0 and 0.5")


def source_hashes(project_root: Path, config: Mapping[str, Any]) -> dict[str, str]:
    out: dict[str, str] = {}
    for value in config["outputs"].get("implementation_hash_components", []):
        path = project_root / str(value)
        if not path.is_file():
            raise FileNotFoundError(f"declared implementation component missing: {value}")
        out[str(value)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return out


def load_references(project_root: Path, config: Mapping[str, Any]) -> tuple[list[dict[str, Any]], str]:
    ref_spec = config["references"]
    source = project_root / str(ref_spec["source_config"])
    source_payload = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    reference_config = {"references": source_payload[str(ref_spec.get("source_section", "references"))]}
    rows = latent.load_reference_catalog(project_root, reference_config)
    if bool(config.get("_runtime", {}).get("quick")):
        permitted = set(str(x) for x in ref_spec.get("quick_reference_ids", []))
        rows = [row for row in rows if str(row["reference_id"]) in permitted]
    if not rows:
        raise ValueError("reference panel is empty")
    return rows, hashlib.sha256(source.read_bytes()).hexdigest()


def sigmoid(logit: np.ndarray | float) -> np.ndarray:
    """AdaptEnv's direct float64 sigmoid, including its public +/-500 clip."""

    x = np.asarray(logit, dtype=np.float64)
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500.0, 500.0)))


def analytic_sigmoid(logit: np.ndarray | float) -> np.ndarray:
    """Unclipped logistic used only for the stable attenuation diagnostic."""

    x = np.asarray(logit, dtype=np.float64)
    out = np.empty_like(x)
    positive = x >= 0.0
    out[positive] = 1.0 / (1.0 + np.exp(-np.minimum(x[positive], 745.0)))
    exp_x = np.exp(np.maximum(x[~positive], -745.0))
    out[~positive] = exp_x / (1.0 + exp_x)
    return out


def direct_sigmoid_residual(
    base_logit: np.ndarray | float,
    latent_shift: np.ndarray | float,
    *,
    output_scale: float,
    steepness: float,
) -> np.ndarray:
    base = np.asarray(base_logit, dtype=np.float64)
    shift = float(steepness) * np.asarray(latent_shift, dtype=np.float64)
    return float(output_scale) * (sigmoid(base + shift) - sigmoid(base))


def stable_sigmoid_residual(
    base_logit: np.ndarray | float,
    latent_shift: np.ndarray | float,
    *,
    output_scale: float,
    steepness: float,
) -> np.ndarray:
    """Evaluate C[sigma(z+k*u)-sigma(z)] without plateau subtraction."""

    z, u = np.broadcast_arrays(
        np.asarray(base_logit, dtype=np.float64),
        float(steepness) * np.asarray(latent_shift, dtype=np.float64),
    )
    target = z + u
    out = np.empty_like(z)
    high = (z >= 0.0) & (target >= 0.0)
    low = (z <= 0.0) & (target <= 0.0)
    cross = ~(high | low)
    if np.any(high):
        q = np.exp(np.maximum(-z[high], -745.0))
        qt = np.exp(np.maximum(-target[high], -745.0))
        delta = u[high]
        numerator = q - qt
        close = np.abs(delta) < 50.0
        numerator[close] = q[close] * (-np.expm1(-delta[close]))
        out[high] = numerator / ((1.0 + q) * (1.0 + qt))
    if np.any(low):
        e = np.exp(np.maximum(z[low], -745.0))
        et = np.exp(np.maximum(target[low], -745.0))
        delta = u[low]
        numerator = et - e
        close = np.abs(delta) < 50.0
        numerator[close] = e[close] * np.expm1(delta[close])
        out[low] = numerator / ((1.0 + e) * (1.0 + et))
    if np.any(cross):
        out[cross] = analytic_sigmoid(target[cross]) - analytic_sigmoid(z[cross])
    return float(output_scale) * out


def normalized_secant_gain(
    base_logit: np.ndarray,
    latent_shift: np.ndarray,
    stable_residual: np.ndarray,
    *,
    output_scale: float,
    steepness: float,
    floor: float,
) -> np.ndarray:
    u = np.asarray(latent_shift, dtype=np.float64)
    residual = np.asarray(stable_residual, dtype=np.float64)
    denominator = (float(output_scale) * float(steepness) / 4.0) * u
    out = np.empty_like(u)
    regular = np.abs(u) > float(floor)
    out[regular] = residual[regular] / denominator[regular]
    p = sigmoid(np.asarray(base_logit, dtype=np.float64))
    out[~regular] = 4.0 * p[~regular] * (1.0 - p[~regular])
    return out


def finite(values: np.ndarray) -> bool:
    return bool(np.isfinite(np.asarray(values, dtype=np.float64)).all())


def require_finite(label: str, *values: np.ndarray | float) -> None:
    if not all(finite(np.asarray(value)) for value in values):
        raise FloatingPointError(f"non-finite values in {label}")


def identity_readout_residual(latent_shift: np.ndarray | float) -> np.ndarray:
    """Explicit identity-readout control used to audit residual algebra."""

    return np.asarray(latent_shift, dtype=np.float64)


def safe_correlation(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2 or float(np.std(x)) <= 0.0 or float(np.std(y)) <= 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def sufficient_stats(
    x: np.ndarray,
    y: np.ndarray,
    *,
    resolution_floor: float,
) -> dict[str, Any]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    n = int(x.size)
    if n == 0:
        return {
            "n": 0, "sum_x": 0.0, "sum_y": 0.0, "sum_x2": 0.0,
            "sum_y2": 0.0, "sum_xy": 0.0, "x_rms": float("nan"),
            "y_rms": float("nan"), "correlation": float("nan"),
            "projection": float("nan"), "practical_resolution_fraction": float("nan"),
            "exact_zero_fraction": float("nan"), "sign_agreement_fraction": float("nan"),
        }
    sx = float(np.sum(x)); sy = float(np.sum(y))
    sx2 = float(np.dot(x, x)); sy2 = float(np.dot(y, y)); sxy = float(np.dot(x, y))
    nonzero = (x != 0.0) | (y != 0.0)
    sign_agree = float(np.mean(np.signbit(x[nonzero]) == np.signbit(y[nonzero]))) if np.any(nonzero) else 1.0
    return {
        "n": n,
        "sum_x": sx,
        "sum_y": sy,
        "sum_x2": sx2,
        "sum_y2": sy2,
        "sum_xy": sxy,
        "x_mean": float(np.mean(x)),
        "y_mean": float(np.mean(y)),
        "x_rms": float(math.sqrt(sx2 / n)),
        "y_rms": float(math.sqrt(sy2 / n)),
        "y_abs_median": float(np.median(np.abs(y))),
        "y_abs_q95": float(np.quantile(np.abs(y), 0.95)),
        "correlation": safe_correlation(x, y),
        "projection": float(sxy / sx2) if sx2 > 0.0 else float("nan"),
        "practical_resolution_fraction": float(np.mean(np.abs(y) >= float(resolution_floor))),
        "exact_zero_fraction": float(np.mean(y == 0.0)),
        "sign_agreement_fraction": sign_agree,
    }


def linearization_stats(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    observed = np.asarray(observed, dtype=np.float64)
    predicted = np.asarray(predicted, dtype=np.float64)
    denominator = float(np.sqrt(np.mean(observed * observed)))
    rmse = float(np.sqrt(np.mean((observed - predicted) ** 2)))
    return {
        "linear_correlation": safe_correlation(observed, predicted),
        "linear_nrmse": float(rmse / denominator) if denominator > 0.0 else float("nan"),
    }


def permutation_projection(
    x: np.ndarray,
    y: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, float | int]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    denom = float(np.dot(x, x))
    if x.size < 2 or denom <= 0.0 or float(np.var(x)) <= 0.0:
        return {
            "observed_projection": float("nan"), "null_mean": float("nan"),
            "null_q95": float("nan"), "one_sided_pvalue": float("nan"),
            "permutation_seed": int(seed),
        }
    observed = float(np.dot(x, y) / denom)
    rng = np.random.default_rng(int(seed))
    null = np.empty(int(replicates), dtype=np.float64)
    for index in range(int(replicates)):
        null[index] = float(np.dot(rng.permutation(x), y) / denom)
    pvalue = float((1 + np.sum(null >= observed)) / (int(replicates) + 1))
    return {
        "observed_projection": observed,
        "null_mean": float(np.mean(null)),
        "null_q95": float(np.quantile(null, 0.95)),
        "one_sided_pvalue": pvalue,
        "permutation_seed": int(seed),
    }


def build_additive_core(reference: np.ndarray, seed: int, config: Mapping[str, Any]) -> FitnessLandscapeEnv:
    core_cfg = config["additive_core"]
    readout = config["readout"]
    landscape_config = LandscapeConfig(
        model_system=str(core_cfg["model_system"]),
        L=int(reference.size),
        alphabet_size=int(core_cfg.get("alphabet_size", 20)),
        seed=int(seed),
        stability_margin=float(readout["public_base"]["stability_margin"]),
        stability_scale=float(core_cfg.get("stability_scale", 1.0)),
        use_blosum_kernel=bool(core_cfg.get("use_blosum_kernel", True)),
        blosum_blend=float(core_cfg.get("blosum_blend", 0.5)),
        threshold_steepness=float(readout["threshold_steepness"]),
        fitness_output_scale=float(readout["output_scale"]),
        observed_fitness_combine_mode=str(core_cfg["observed_fitness_combine_mode"]),
        generic_epistasis_target=str(core_cfg["generic_epistasis_target"]),
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
    return FitnessLandscapeEnv(landscape_config, reference_sequence=reference)


def additive_surrogate_matrix(additive_matrix: np.ndarray, *, seed: int) -> np.ndarray:
    length, alphabet_size = additive_matrix.shape
    rng = np.random.default_rng(stable_uint64("EXP-061-additive-surrogate", int(seed), length))
    matrix = rng.standard_normal((length, alphabet_size))
    matrix -= np.mean(matrix, axis=1, keepdims=True)
    target = additive_mutation_energy(additive_matrix)
    realized = additive_mutation_energy(matrix)
    return matrix * math.sqrt(target / realized)


def score_additive_field(matrix: np.ndarray, reference: np.ndarray, sequences: np.ndarray) -> np.ndarray:
    positions = np.arange(matrix.shape[0], dtype=np.int64)
    raw = np.sum(matrix[positions[None, :], sequences], axis=1)
    center = float(np.sum(matrix[positions, reference]))
    return np.asarray(raw - center, dtype=np.float64)


def complete_single_mutants(reference: np.ndarray, alphabet_size: int = 20) -> tuple[np.ndarray, np.ndarray]:
    length = int(reference.size)
    sources = np.repeat(reference[None, :], length * (alphabet_size - 1), axis=0)
    targets = sources.copy()
    row = 0
    for site in range(length):
        for allele in range(alphabet_size):
            if allele == int(reference[site]):
                continue
            targets[row, site] = int(allele)
            row += 1
    return sources, targets


def sample_pairs(
    reference: np.ndarray,
    *,
    seed: int,
    regime: str,
    distance_label: str,
    hamming: int,
    count: int,
) -> tuple[np.ndarray, np.ndarray]:
    return latent.generate_pairs(
        reference=reference,
        seed=int(seed),
        regime=str(regime),
        distance_label=str(distance_label),
        hamming_distance=int(hamming),
        n_pairs=int(count),
        alphabet_size=20,
    )


def pair_hash(sources: np.ndarray, targets: np.ndarray, **labels: Any) -> str:
    return latent.json_hash(
        {"sources": latent.array_hash(sources), "targets": latent.array_hash(targets), **labels}
    )


def occupancy_strata(source: np.ndarray, target: np.ndarray, band: Sequence[float]) -> np.ndarray:
    low, high = map(float, band)
    source_dynamic = (source >= low) & (source <= high)
    target_dynamic = (target >= low) & (target <= high)
    out = np.full(source.shape, "one_crossing", dtype=object)
    out[source_dynamic & target_dynamic] = "both_dynamic"
    out[(source < low) & (target < low)] = "same_low"
    out[(source > high) & (target > high)] = "same_high"
    return out


def weight_grid(config: Mapping[str, Any], candidate_mutation_rms: float) -> list[dict[str, Any]]:
    k = float(config["readout"]["threshold_steepness"])
    rows: list[dict[str, Any]] = []
    for amplitude in config["intervention"]["dimensionless_amplitudes"]:
        amplitude = float(amplitude)
        rows.append({
            "weight_label": f"a_{amplitude:g}",
            "dimensionless_amplitude": amplitude,
            "weight": float(amplitude / (k * candidate_mutation_rms)),
            "weight_kind": "operator_scaled",
        })
    for weight in config["intervention"]["fixed_weights"]:
        weight = float(weight)
        rows.append({
            "weight_label": f"fixed_{weight:g}",
            "dimensionless_amplitude": float(k * weight * candidate_mutation_rms),
            "weight": weight,
            "weight_kind": "fixed",
        })
    deduped: dict[tuple[str, float], dict[str, Any]] = {}
    for row in rows:
        deduped[(str(row["weight_label"]), float(row["weight"]))] = row
    return list(deduped.values())


def prefixed(values: Mapping[str, Any], prefix: str) -> dict[str, Any]:
    return {f"{prefix}{key}": value for key, value in values.items()}


def build_pair_panels(
    reference: np.ndarray,
    *,
    seed: int,
    config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    probe = config["probe"]
    panels: list[dict[str, Any]] = []
    sources, targets = complete_single_mutants(reference)
    panels.append({
        "sample_id": "reference_shell_h1_complete",
        "regime": "reference_shell",
        "distance_label": "fixed_1",
        "hamming_distance": 1,
        "sources": sources,
        "targets": targets,
    })
    sources, targets = sample_pairs(
        reference,
        seed=seed,
        regime="uniform_background",
        distance_label="fixed_1",
        hamming=1,
        count=int(probe["uniform_background_single_mutant_pairs"]),
    )
    panels.append({
        "sample_id": "uniform_background_h1",
        "regime": "uniform_background",
        "distance_label": "fixed_1",
        "hamming_distance": 1,
        "sources": sources,
        "targets": targets,
    })
    for regime in ("reference_shell", "uniform_background"):
        for label in probe["fractional_hamming_distances"]:
            fraction = float(str(label).removesuffix("L"))
            hamming = max(1, min(int(reference.size), int(round(fraction * int(reference.size)))))
            sources, targets = sample_pairs(
                reference,
                seed=seed,
                regime=regime,
                distance_label=f"frac_{fraction:g}",
                hamming=hamming,
                count=int(probe["pairs_per_fractional_distance"]),
            )
            panels.append({
                "sample_id": f"{regime}_frac_{fraction:g}",
                "regime": regime,
                "distance_label": f"frac_{fraction:g}",
                "hamming_distance": hamming,
                "sources": sources,
                "targets": targets,
            })
    for panel in panels:
        panel["pairing_hash"] = pair_hash(
            panel["sources"], panel["targets"],
            regime=panel["regime"], distance_label=panel["distance_label"],
            hamming_distance=panel["hamming_distance"],
        )
    return panels


def generate_squares(
    field: ExtensiveGNKField,
    reference: np.ndarray,
    *,
    generator_seed: int,
    stratum: str,
    background: str,
    count: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(
        stable_uint64("EXP-061-squares", int(generator_seed), int(reference.size), stratum, background, count)
    )
    if stratum.startswith("candidate_supported_order_"):
        order = int(stratum.rsplit("_", 1)[1])
        pair_pool = field.support_pairs(order)
        if not pair_pool:
            raise ValueError(f"candidate has no supported pairs at order {order}")
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
    first = base.copy(); second = base.copy(); both = base.copy()
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
    design_hash = latent.json_hash({
        "base": latent.array_hash(base), "first": latent.array_hash(first),
        "second": latent.array_hash(second), "both": latent.array_hash(both),
        "sites": latent.array_hash(sites), "stratum": stratum, "background": background,
    })
    return {
        "base": base, "first": first, "second": second, "both": both,
        "sites": sites, "square_hash": design_hash,
    }


def score_nodes(
    scorer: Callable[[np.ndarray], np.ndarray],
    arrays: Sequence[np.ndarray],
    *,
    chunk_size: int,
) -> list[np.ndarray]:
    sizes = [int(array.shape[0]) for array in arrays]
    values = latent.score_in_chunks(scorer, np.concatenate(arrays, axis=0), chunk_size=chunk_size)
    out: list[np.ndarray] = []
    start = 0
    for size in sizes:
        out.append(values[start : start + size])
        start += size
    return out


def logit_for_occupancy(occupancy: float) -> float:
    p = float(occupancy)
    if not 0.0 < p < 1.0:
        raise ValueError("occupancy must lie strictly between zero and one")
    return float(math.log(p / (1.0 - p)))


def scope_specs(
    *,
    additive_reference_score: float,
    public_theta_margin8: float,
    config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    k = float(config["readout"]["threshold_steepness"])
    rows: list[dict[str, Any]] = []
    for p in config["readout"]["mechanistic_reference_occupancies"]:
        p = float(p)
        rows.append({
            "readout_scope": "reference_anchored",
            "operating_point": f"p_ref_{p:g}",
            "occupancy": p,
            "theta": float(additive_reference_score - logit_for_occupancy(p) / k),
            "deployment_margin": float("nan"),
        })
    base_margin = float(config["readout"]["deployment_margins"]["base_template_margin8"])
    constructor_margin = float(config["readout"]["deployment_margins"]["python_constructor_default"])
    rows.extend([
        {
            "readout_scope": "deployment",
            "operating_point": "public_base_margin8",
            "occupancy": float("nan"),
            "theta": float(public_theta_margin8),
            "deployment_margin": base_margin,
        },
        {
            "readout_scope": "deployment",
            "operating_point": "constructor_default_margin10",
            "occupancy": float("nan"),
            "theta": float(public_theta_margin8 - (constructor_margin - base_margin)),
            "deployment_margin": constructor_margin,
        },
    ])
    return rows


def evaluate_pair_panel(
    *,
    reference_id: str,
    length: int,
    generator_seed: int,
    panel: Mapping[str, Any],
    additive_matrix: np.ndarray,
    additive_reference_score: float,
    public_theta_margin8: float,
    field_scores: Mapping[str, tuple[np.ndarray, np.ndarray]],
    weights: Sequence[Mapping[str, Any]],
    config: Mapping[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    sources = np.asarray(panel["sources"], dtype=np.int64)
    targets = np.asarray(panel["targets"], dtype=np.int64)
    a_source = latent.additive_scores(additive_matrix, sources)
    a_target = latent.additive_scores(additive_matrix, targets)
    require_finite("pair additive scores", a_source, a_target)
    delta_a = a_target - a_source
    k = float(config["readout"]["threshold_steepness"])
    c = float(config["readout"]["output_scale"])
    resolution = float(config["acceptance"]["operational_detection_threshold"])
    gamma_floor = float(config["readout"]["pointwise_latent_shift_abs_floor"])
    dynamic_band = config["readout"]["dynamic_occupancy_band"]
    base_ids = {
        "reference_id": reference_id, "length": int(length),
        "generator_seed": int(generator_seed), "sample_id": str(panel["sample_id"]),
        "regime": str(panel["regime"]), "distance_label": str(panel["distance_label"]),
        "hamming_distance": int(panel["hamming_distance"]),
        "pairing_hash": str(panel["pairing_hash"]),
    }
    paired_rows: list[dict[str, Any]] = []
    mutation_rows: list[dict[str, Any]] = []
    occupancy_rows: list[dict[str, Any]] = []
    topology_rows: list[dict[str, Any]] = []

    # Genotype-wise transfer calibration: theta changes per genotype and is not a
    # landscape.  Local panels suffice; repeating the transfer-function grid at
    # extensive Hamming distances would add rows but no new transfer estimand.
    if int(panel["hamming_distance"]) == 1:
        for occupancy in config["readout"]["pointwise_occupancies"]:
            occupancy = float(occupancy)
            # Use each target genotype once.  Repeating the reference source 19L
            # times would turn G(reference)=0 into an artificial 50% zero mass.
            z = np.full(a_target.size, logit_for_occupancy(occupancy), dtype=np.float64)
            for arm, (g_source, g_target) in field_scores.items():
                g_all = g_target
                for weight_row in weights:
                    w = float(weight_row["weight"])
                    u = w * g_all
                    direct = direct_sigmoid_residual(z, u, output_scale=c, steepness=k)
                    stable = stable_sigmoid_residual(z, u, output_scale=c, steepness=k)
                    gamma = normalized_secant_gain(
                        z, u, stable, output_scale=c, steepness=k, floor=gamma_floor,
                    )
                    linear = c * k * occupancy * (1.0 - occupancy) * u
                    require_finite("pointwise residuals", u, direct, stable, gamma, linear)
                    row = {
                        **base_ids, "arm": arm, "readout_scope": "pointwise_occupancy",
                        "operating_point": f"p_point_{occupancy:g}", "occupancy": occupancy,
                        "theta": "per_genotype", **weight_row,
                        **prefixed(sufficient_stats(u, stable, resolution_floor=resolution), "stable_"),
                        **prefixed(sufficient_stats(u, direct, resolution_floor=resolution), "direct_"),
                        **linearization_stats(stable, linear),
                        "gamma_median": float(np.median(gamma)),
                        "gamma_q05": float(np.quantile(gamma, 0.05)),
                        "max_direct_stable_abs_difference": float(np.max(np.abs(direct - stable))),
                        "scientific_landscape_scope": False,
                    }
                    paired_rows.append(row)

    for scope in scope_specs(
        additive_reference_score=additive_reference_score,
        public_theta_margin8=public_theta_margin8,
        config=config,
    ):
        fractional_panel = int(panel["hamming_distance"]) != 1
        if (
            fractional_panel
            and scope["readout_scope"] == "reference_anchored"
            and str(scope["operating_point"]) != "p_ref_0.5"
        ):
            continue
        theta = float(scope["theta"])
        z_source = k * (a_source - theta)
        z_target = k * (a_target - theta)
        p_source = sigmoid(z_source); p_target = sigmoid(z_target)
        require_finite("edge occupancy calibration", theta, z_source, z_target, p_source, p_target)
        strata = occupancy_strata(p_source, p_target, dynamic_band)
        occupancy_rows.append({
            **base_ids, **scope, "n": int(a_source.size),
            "source_occupancy_median": float(np.median(p_source)),
            "target_occupancy_median": float(np.median(p_target)),
            "both_dynamic_fraction": float(np.mean(strata == "both_dynamic")),
            "one_crossing_fraction": float(np.mean(strata == "one_crossing")),
            "same_low_fraction": float(np.mean(strata == "same_low")),
            "same_high_fraction": float(np.mean(strata == "same_high")),
        })
        for arm, (g_source, g_target) in field_scores.items():
            for weight_row in weights:
                if fractional_panel and str(weight_row["weight_label"]) not in {
                    "a_0", "a_0.5", "a_1", "fixed_0.6", "fixed_1",
                }:
                    continue
                w = float(weight_row["weight"])
                us = w * g_source; ut = w * g_target
                direct_s = direct_sigmoid_residual(z_source, us, output_scale=c, steepness=k)
                direct_t = direct_sigmoid_residual(z_target, ut, output_scale=c, steepness=k)
                stable_s = stable_sigmoid_residual(z_source, us, output_scale=c, steepness=k)
                stable_t = stable_sigmoid_residual(z_target, ut, output_scale=c, steepness=k)
                delta_direct = direct_t - direct_s
                delta_stable = stable_t - stable_s
                x = ut - us
                linear_s = c * k * p_source * (1.0 - p_source) * us
                linear_t = c * k * p_target * (1.0 - p_target) * ut
                delta_linear = linear_t - linear_s
                require_finite(
                    "edge residuals", us, ut, direct_s, direct_t, stable_s,
                    stable_t, delta_direct, delta_stable, delta_linear,
                )
                labels = (
                    ["all", "both_dynamic", "one_crossing", "same_low", "same_high"]
                    if str(panel["sample_id"]) == "reference_shell_h1_complete"
                    else ["all"]
                )
                for label in labels:
                    mask = np.ones(x.shape, dtype=bool) if label == "all" else strata == label
                    if not np.any(mask):
                        continue
                    mutation_rows.append({
                        **base_ids, **scope, "arm": arm, **weight_row,
                        "occupancy_stratum": label,
                        "n": int(np.sum(mask)),
                        "delta_A_mean": float(np.mean(delta_a[mask])),
                        "delta_A_rms": float(np.sqrt(np.mean(delta_a[mask] ** 2))),
                        "core_occupancy_source_mean": float(np.mean(p_source[mask])),
                        "core_occupancy_target_mean": float(np.mean(p_target[mask])),
                        **prefixed(sufficient_stats(x[mask], delta_stable[mask], resolution_floor=resolution), "stable_"),
                        **prefixed(sufficient_stats(x[mask], delta_direct[mask], resolution_floor=resolution), "direct_"),
                        **linearization_stats(delta_stable[mask], delta_linear[mask]),
                        "max_direct_stable_abs_difference": float(np.max(np.abs(delta_direct[mask] - delta_stable[mask]))),
                    })
                latent_delta = delta_a + x
                full_source = c * sigmoid(z_source + k * us)
                full_target = c * sigmoid(z_target + k * ut)
                observed_delta = full_target - full_source
                non_tie = observed_delta != 0.0
                agreement = (
                    float(np.mean(np.signbit(observed_delta[non_tie]) == np.signbit(latent_delta[non_tie])))
                    if np.any(non_tie) else 1.0
                )
                topology_rows.append({
                    **base_ids, **scope, "arm": arm, **weight_row,
                    "n": int(observed_delta.size), "non_tie_count": int(np.sum(non_tie)),
                    "direct_tie_fraction": float(np.mean(~non_tie)),
                    "non_tie_rank_agreement": agreement,
                    "non_tie_mismatch_count": int(np.sum(
                        np.signbit(observed_delta[non_tie]) != np.signbit(latent_delta[non_tie])
                    )),
                })
    return {
        "paired_rows": paired_rows,
        "mutation_rows": mutation_rows,
        "occupancy_rows": occupancy_rows,
        "topology_rows": topology_rows,
    }


def evaluate_square_panel(
    *,
    reference_id: str,
    length: int,
    generator_seed: int,
    stratum: str,
    background: str,
    square: Mapping[str, Any],
    additive_matrix: np.ndarray,
    additive_reference_score: float,
    public_theta_margin8: float,
    field_scores: Mapping[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    weights: Sequence[Mapping[str, Any]],
    config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    arrays = [square["base"], square["first"], square["second"], square["both"]]
    additive = [latent.additive_scores(additive_matrix, np.asarray(array)) for array in arrays]
    require_finite("square additive scores", *additive)
    k = float(config["readout"]["threshold_steepness"])
    c = float(config["readout"]["output_scale"])
    resolution = float(config["acceptance"]["operational_detection_threshold"])
    permutation_reps = int(config["acceptance"]["permutation_replicates"])
    alpha = float(config["acceptance"]["texture_permutation_alpha"])
    rows: list[dict[str, Any]] = []
    support_order = (
        int(stratum.rsplit("_", 1)[1])
        if stratum.startswith("candidate_supported_order_") else 0
    )
    ids = {
        "reference_id": reference_id, "length": int(length),
        "generator_seed": int(generator_seed), "square_stratum": stratum,
        "support_order": support_order, "square_background": background,
        "square_hash": str(square["square_hash"]),
    }
    for scope in scope_specs(
        additive_reference_score=additive_reference_score,
        public_theta_margin8=public_theta_margin8,
        config=config,
    ):
        theta = float(scope["theta"])
        for core_mode in config["probe"]["square_core_modes"]:
            if core_mode == "actual_additive_core":
                a_nodes = additive
            elif core_mode == "frozen_additive_core":
                a_nodes = [additive[0], additive[0], additive[0], additive[0]]
            else:
                raise ValueError(f"unsupported square core mode: {core_mode}")
            z_nodes = [k * (values - theta) for values in a_nodes]
            for arm, g_nodes in field_scores.items():
                g0, ga, gb, gab = g_nodes
                latent_epsilon = gab - ga - gb + g0
                gab_additivized = ga + gb - g0
                for weight_row in weights:
                    if scope["readout_scope"] == "reference_anchored":
                        if str(scope["operating_point"]) == "p_ref_0.5":
                            keep_condition = str(weight_row["weight_kind"]) == "operator_scaled"
                        else:
                            keep_condition = (
                                arm == ARM_CANDIDATE
                                and str(weight_row["weight_label"]) == "a_0.5"
                            )
                    else:
                        keep_condition = str(weight_row["weight_label"]) in {
                            "a_0", "fixed_0.6", "fixed_1",
                        }
                    if not keep_condition:
                        continue
                    w = float(weight_row["weight"])
                    direct_nodes = [
                        direct_sigmoid_residual(z, w * g, output_scale=c, steepness=k)
                        for z, g in zip(z_nodes, g_nodes)
                    ]
                    stable_nodes = [
                        stable_sigmoid_residual(z, w * g, output_scale=c, steepness=k)
                        for z, g in zip(z_nodes, g_nodes)
                    ]
                    epsilon_direct = direct_nodes[3] - direct_nodes[1] - direct_nodes[2] + direct_nodes[0]
                    epsilon_stable = stable_nodes[3] - stable_nodes[1] - stable_nodes[2] + stable_nodes[0]
                    x = w * latent_epsilon
                    # tau_texture is exactly a single secant at the double-mutant
                    # endpoint.  Evaluate that secant directly; subtracting two R
                    # values would reintroduce the cancellation this experiment is
                    # designed to diagnose.
                    ab_null_logit = z_nodes[3] + k * w * gab_additivized
                    tau_direct = direct_sigmoid_residual(
                        ab_null_logit, x, output_scale=c, steepness=k,
                    )
                    tau_stable = stable_sigmoid_residual(
                        ab_null_logit, x, output_scale=c, steepness=k,
                    )
                    p_ab_null = sigmoid(ab_null_logit)
                    tau_linear = c * k * p_ab_null * (1.0 - p_ab_null) * x
                    require_finite(
                        "square residuals", *direct_nodes, *stable_nodes,
                        epsilon_direct, epsilon_stable, tau_direct, tau_stable,
                        tau_linear, latent_epsilon,
                    )
                    stable_stats = sufficient_stats(x, tau_stable, resolution_floor=resolution)
                    direct_stats = sufficient_stats(x, tau_direct, resolution_floor=resolution)
                    epsilon_stats = sufficient_stats(x, epsilon_stable, resolution_floor=resolution)
                    primary_permutation = bool(
                        arm == ARM_CANDIDATE
                        and support_order > 0
                        and (
                            (
                                scope["readout_scope"] == "reference_anchored"
                                and abs(float(scope["occupancy"]) - 0.5) <= 1e-12
                                and str(weight_row["weight_label"]) == "a_0.5"
                            )
                            or (
                                scope["operating_point"] == "public_base_margin8"
                                and str(weight_row["weight_label"]) == "fixed_0.6"
                            )
                        )
                    )
                    if primary_permutation and float(np.var(x)) > 0.0:
                        perm_seed = stable_uint64(
                            "EXP-061-permutation", generator_seed, reference_id, stratum,
                            background, core_mode, scope["operating_point"], weight_row["weight_label"],
                        )
                        permutation = permutation_projection(
                            x, tau_stable if scope["readout_scope"] != "deployment" else tau_direct,
                            replicates=permutation_reps, seed=perm_seed,
                        )
                    else:
                        permutation = {
                            "observed_projection": float("nan"), "null_mean": float("nan"),
                            "null_q95": float("nan"), "one_sided_pvalue": float("nan"),
                            "permutation_seed": -1,
                        }
                    rows.append({
                        **ids, **scope, "additive_core_mode": core_mode, "arm": arm,
                        **weight_row, "n": int(x.size),
                        "latent_epsilon_rms": float(np.sqrt(np.mean(latent_epsilon ** 2))),
                        "latent_epsilon_max_abs": float(np.max(np.abs(latent_epsilon))),
                        "epsilon_R_stable_rms": float(np.sqrt(np.mean(epsilon_stable ** 2))),
                        "epsilon_R_direct_rms": float(np.sqrt(np.mean(epsilon_direct ** 2))),
                        "tau_texture_definition": "epsilon_R_actual_minus_candidate_additivized",
                        **prefixed(stable_stats, "stable_tau_"),
                        **prefixed(direct_stats, "direct_tau_"),
                        **prefixed(epsilon_stats, "stable_epsilon_"),
                        **linearization_stats(tau_stable, tau_linear),
                        **permutation,
                        "permutation_required_primary": primary_permutation,
                        "permutation_pass": bool(
                            primary_permutation
                            and np.isfinite(float(permutation["one_sided_pvalue"]))
                            and float(permutation["one_sided_pvalue"]) <= alpha
                        ) if primary_permutation else None,
                        "max_direct_stable_tau_abs_difference": float(np.max(np.abs(tau_direct - tau_stable))),
                    })
    return rows


def make_field_bundle(
    reference: np.ndarray,
    *,
    generator_seed: int,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    core = build_additive_core(reference, generator_seed, config)
    additive_matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    additive_energy = additive_mutation_energy(additive_matrix)
    additive_reference_score = float(latent.additive_scores(additive_matrix, reference[None, :])[0])
    public_theta_margin8 = float(core.stability_threshold)

    # Retrofitting occurs only after the additive matrix/threshold are frozen.
    latent.retrofit_current_latent_gnk(core, config)
    legacy_catalog, legacy_energy = latent.compile_legacy_catalog(core)
    legacy_reference_raw = float(core._gnk_epistatic_components_batch(reference[None, :])[0])
    local_scale = math.sqrt(additive_energy / legacy_energy)

    candidate_cfg = config["candidate"]
    candidate = ExtensiveGNKField(
        additive_matrix=additive_matrix,
        reference_sequence=reference,
        seed=int(generator_seed) + int(candidate_cfg.get("seed_offset", 30000)),
        orders=tuple(int(x) for x in candidate_cfg["orders"]),
        order_weights=tuple(float(x) for x in candidate_cfg["order_weights"]),
        target_degrees=tuple(float(x) for x in candidate_cfg["target_degrees"]),
    )
    replay = ExtensiveGNKField(
        additive_matrix=additive_matrix,
        reference_sequence=reference,
        seed=int(generator_seed) + int(candidate_cfg.get("seed_offset", 30000)),
        orders=tuple(int(x) for x in candidate_cfg["orders"]),
        order_weights=tuple(float(x) for x in candidate_cfg["order_weights"]),
        target_degrees=tuple(float(x) for x in candidate_cfg["target_degrees"]),
    )
    surrogate_matrix = additive_surrogate_matrix(
        additive_matrix,
        seed=int(generator_seed) + int(config["additive_surrogate"].get("seed_offset", 50000)),
    )

    def score_public(sequences: np.ndarray) -> np.ndarray:
        return np.asarray(core._gnk_epistatic_components_batch(sequences), dtype=np.float64)

    def score_local(sequences: np.ndarray) -> np.ndarray:
        return local_scale * (score_public(sequences) - legacy_reference_raw)

    def score_candidate(sequences: np.ndarray) -> np.ndarray:
        return candidate.score_batch(sequences)

    def score_surrogate(sequences: np.ndarray) -> np.ndarray:
        return score_additive_field(surrogate_matrix, reference, sequences)

    scorers: dict[str, Callable[[np.ndarray], np.ndarray]] = {
        ARM_CANDIDATE: score_candidate,
        ARM_LOCAL: score_local,
        ARM_PUBLIC: score_public,
        ARM_ADDITIVE: score_surrogate,
    }
    reference_scores = {
        arm: float(scorer(reference[None, :])[0]) for arm, scorer in scorers.items()
    }
    field_hashes = {
        ARM_CANDIDATE: candidate.fingerprint,
        ARM_LOCAL: latent.json_hash({"legacy": latent.legacy_catalog_hash(legacy_catalog), "scale": local_scale}),
        ARM_PUBLIC: latent.legacy_catalog_hash(legacy_catalog),
        ARM_ADDITIVE: latent.array_hash(surrogate_matrix),
    }
    probe = reference[None, :]
    public_formula = float(config["readout"]["output_scale"]) * sigmoid(
        float(config["readout"]["threshold_steepness"])
        * (latent.additive_scores(additive_matrix, probe) - public_theta_margin8)
    )
    public_implementation = float(config["readout"]["output_scale"]) * core._stability_gates_batch(
        latent.additive_scores(additive_matrix, probe)
    )
    return {
        "core": core,
        "additive_matrix": additive_matrix,
        "additive_energy": float(additive_energy),
        "additive_reference_score": additive_reference_score,
        "public_theta_margin8": public_theta_margin8,
        "candidate": candidate,
        "scorers": scorers,
        "field_hashes": field_hashes,
        "reference_scores": reference_scores,
        "candidate_replay_match": bool(candidate.fingerprint == replay.fingerprint),
        "local_scale": float(local_scale),
        "legacy_energy": float(legacy_energy),
        "public_formula_parity_max_abs": float(np.max(np.abs(public_formula - public_implementation))),
    }


def run_task(payload: Mapping[str, Any]) -> dict[str, Any]:
    config = dict(payload["config"])
    reference_row = dict(payload["reference"])
    generator_seed = int(payload["generator_seed"])
    task_signature = str(payload["task_signature"])
    checkpoint = Path(str(payload["checkpoint_path"]))
    if checkpoint.exists():
        cached = json.loads(checkpoint.read_text(encoding="utf-8"))
        if str(cached.get("task", {}).get("task_signature")) == task_signature:
            return cached

    reference_id = str(reference_row["reference_id"])
    reference = latent.sequence_to_indices(str(reference_row["sequence"]))
    length = int(reference.size)
    chunk_size = int(config["probe"]["eval_chunk_size"])
    bundle = make_field_bundle(reference, generator_seed=generator_seed, config=config)
    candidate: ExtensiveGNKField = bundle["candidate"]
    weights = weight_grid(config, math.sqrt(float(bundle["additive_energy"])))

    design_rows: list[dict[str, Any]] = []
    for arm in ARMS:
        design_rows.append({
            "row_type": "field", "reference_id": reference_id, "length": length,
            "generator_seed": generator_seed, "arm": arm,
            "field_hash": bundle["field_hashes"][arm],
            "reference_score": bundle["reference_scores"][arm],
            "candidate_replay_match": bundle["candidate_replay_match"] if arm == ARM_CANDIDATE else None,
            "additive_mutation_energy": bundle["additive_energy"],
            "candidate_mutation_rms": math.sqrt(float(bundle["additive_energy"])),
            "legacy_operator_energy": bundle["legacy_energy"] if arm in {ARM_LOCAL, ARM_PUBLIC} else None,
            "legacy_local_scale": bundle["local_scale"] if arm == ARM_LOCAL else None,
            "public_formula_parity_max_abs": bundle["public_formula_parity_max_abs"],
            "task_signature": task_signature,
        })

    paired_rows: list[dict[str, Any]] = []
    mutation_rows: list[dict[str, Any]] = []
    occupancy_rows: list[dict[str, Any]] = []
    topology_rows: list[dict[str, Any]] = []
    pair_panels = build_pair_panels(reference, seed=generator_seed, config=config)
    replay_panels = build_pair_panels(reference, seed=generator_seed, config=config)
    pair_replay = {
        str(panel["sample_id"]): str(panel["pairing_hash"]) for panel in replay_panels
    }
    public_end_to_end_parity = 0.0
    identity_point_max_error = 0.0
    identity_edge_max_error = 0.0
    identity_square_max_error = 0.0
    raw_arrays_finite = bool(
        finite(bundle["additive_matrix"])
        and np.isfinite(float(bundle["additive_energy"]))
        and np.isfinite(float(bundle["additive_reference_score"]))
        and np.isfinite(float(bundle["public_theta_margin8"]))
    )
    for panel_index, panel in enumerate(pair_panels):
        sources = np.asarray(panel["sources"], dtype=np.int64)
        targets = np.asarray(panel["targets"], dtype=np.int64)
        arrays = [sources, targets]
        field_scores: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for arm, scorer in bundle["scorers"].items():
            scored = score_nodes(scorer, arrays, chunk_size=chunk_size)
            field_scores[arm] = (scored[0], scored[1])
            raw_arrays_finite = raw_arrays_finite and finite(scored[0]) and finite(scored[1])
            for weight_row in weights:
                w = float(weight_row["weight"])
                identity_source = identity_readout_residual(w * scored[0])
                identity_target = identity_readout_residual(w * scored[1])
                identity_point_max_error = max(
                    identity_point_max_error,
                    float(np.max(np.abs(identity_source - w * scored[0]))),
                    float(np.max(np.abs(identity_target - w * scored[1]))),
                )
                identity_edge_max_error = max(
                    identity_edge_max_error,
                    float(np.max(np.abs(
                        (identity_target - identity_source) - w * (scored[1] - scored[0])
                    ))),
                )
        design_rows.append({
            "row_type": "pair_panel", "reference_id": reference_id, "length": length,
            "generator_seed": generator_seed, "sample_id": panel["sample_id"],
            "regime": panel["regime"], "distance_label": panel["distance_label"],
            "hamming_distance": panel["hamming_distance"], "n": int(sources.shape[0]),
            "pairing_hash": panel["pairing_hash"], "sampling_key_excludes_treatment": True,
            "design_replay_match": str(panel["pairing_hash"]) == pair_replay[str(panel["sample_id"])],
            "task_signature": task_signature,
        })
        if panel_index == 0:
            sample = np.concatenate([sources[:16], targets[:16]], axis=0)
            a_sample = latent.additive_scores(bundle["additive_matrix"], sample)
            raw_sample = bundle["scorers"][ARM_PUBLIC](sample)
            deployment_w = float(config["intervention"]["deployment_weight"])
            expected = float(config["readout"]["output_scale"]) * sigmoid(
                float(config["readout"]["threshold_steepness"])
                * (a_sample + deployment_w * raw_sample - bundle["public_theta_margin8"])
            )
            original_weight = float(bundle["core"].cfg.gnk_epistasis_weight)
            bundle["core"].cfg.gnk_epistasis_weight = deployment_w
            observed = np.asarray(bundle["core"]._evaluate_components_batch(sample)["fitness"], dtype=np.float64)
            bundle["core"].cfg.gnk_epistasis_weight = original_weight
            public_end_to_end_parity = float(np.max(np.abs(expected - observed)))
        evaluated = evaluate_pair_panel(
            reference_id=reference_id,
            length=length,
            generator_seed=generator_seed,
            panel=panel,
            additive_matrix=bundle["additive_matrix"],
            additive_reference_score=bundle["additive_reference_score"],
            public_theta_margin8=bundle["public_theta_margin8"],
            field_scores=field_scores,
            weights=weights,
            config=config,
        )
        paired_rows.extend(evaluated["paired_rows"])
        mutation_rows.extend(evaluated["mutation_rows"])
        occupancy_rows.extend(evaluated["occupancy_rows"])
        topology_rows.extend(evaluated["topology_rows"])

    epistasis_rows: list[dict[str, Any]] = []
    for stratum in config["probe"]["square_strata"]:
        for background in config["probe"]["square_backgrounds"]:
            square = generate_squares(
                candidate,
                reference,
                generator_seed=generator_seed,
                stratum=str(stratum),
                background=str(background),
                count=int(config["probe"]["squares_per_stratum"]),
            )
            square_replay = generate_squares(
                candidate,
                reference,
                generator_seed=generator_seed,
                stratum=str(stratum),
                background=str(background),
                count=int(config["probe"]["squares_per_stratum"]),
            )
            arrays = [square["base"], square["first"], square["second"], square["both"]]
            field_scores4: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
            for arm, scorer in bundle["scorers"].items():
                values = score_nodes(scorer, arrays, chunk_size=chunk_size)
                field_scores4[arm] = (values[0], values[1], values[2], values[3])
                raw_arrays_finite = raw_arrays_finite and all(finite(value) for value in values)
                for weight_row in weights:
                    w = float(weight_row["weight"])
                    identity_nodes = [identity_readout_residual(w * value) for value in values]
                    identity_epsilon = (
                        identity_nodes[3] - identity_nodes[1] - identity_nodes[2] + identity_nodes[0]
                    )
                    expected_epsilon = w * (values[3] - values[1] - values[2] + values[0])
                    identity_square_max_error = max(
                        identity_square_max_error,
                        float(np.max(np.abs(identity_epsilon - expected_epsilon))),
                    )
            design_rows.append({
                "row_type": "square_panel", "reference_id": reference_id, "length": length,
                "generator_seed": generator_seed, "square_stratum": stratum,
                "square_background": background, "n": int(square["base"].shape[0]),
                "square_hash": square["square_hash"], "sampling_key_excludes_treatment": True,
                "design_replay_match": str(square["square_hash"]) == str(square_replay["square_hash"]),
                "task_signature": task_signature,
            })
            epistasis_rows.extend(evaluate_square_panel(
                reference_id=reference_id,
                length=length,
                generator_seed=generator_seed,
                stratum=str(stratum),
                background=str(background),
                square=square,
                additive_matrix=bundle["additive_matrix"],
                additive_reference_score=bundle["additive_reference_score"],
                public_theta_margin8=bundle["public_theta_margin8"],
                field_scores=field_scores4,
                weights=weights,
                config=config,
            ))

    # Task-wide execution audits are duplicated on the task row for compact checking.
    centered_arms = (ARM_CANDIDATE, ARM_LOCAL, ARM_ADDITIVE)
    reference_occupancy_errors = [
        abs(float(row["source_occupancy_median"]) - float(row["occupancy"]))
        for row in occupancy_rows
        if row["readout_scope"] == "reference_anchored"
        and row["sample_id"] == "reference_shell_h1_complete"
    ]
    task = {
        "reference_id": reference_id,
        "length": length,
        "generator_seed": generator_seed,
        "task_signature": task_signature,
        "candidate_replay_match": bool(bundle["candidate_replay_match"]),
        "max_centered_reference_abs": float(max(abs(bundle["reference_scores"][arm]) for arm in centered_arms)),
        "public_formula_parity_max_abs": float(bundle["public_formula_parity_max_abs"]),
        "public_end_to_end_parity_max_abs": public_end_to_end_parity,
        "identity_point_max_error": identity_point_max_error,
        "identity_edge_max_error": identity_edge_max_error,
        "identity_square_max_error": identity_square_max_error,
        "reference_occupancy_calibration_max_abs": float(max(reference_occupancy_errors, default=float("inf"))),
        "constructor_margin10_comparator_present": bool(any(
            row["operating_point"] == "constructor_default_margin10" for row in occupancy_rows
        )),
        "all_outputs_finite": bool(
            raw_arrays_finite
            and
            all(np.isfinite(float(row.get("stable_y_rms", 0.0))) for row in paired_rows)
            and all(np.isfinite(float(row.get("stable_y_rms", 0.0))) for row in mutation_rows)
            and all(np.isfinite(float(row.get("stable_tau_y_rms", 0.0))) for row in epistasis_rows)
        ),
    }
    result = {
        "task": task,
        "design_rows": design_rows,
        "occupancy_rows": occupancy_rows,
        "paired_rows": paired_rows,
        "mutation_rows": mutation_rows,
        "epistasis_rows": epistasis_rows,
        "topology_rows": topology_rows,
    }
    latent.write_json_atomic(checkpoint, result)
    return result


def log_length_fit(
    frame: pd.DataFrame,
    *,
    metric: str,
    fit_id: str,
    bootstrap_replicates: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    selected = frame[["length", "generator_seed", metric]].dropna()
    selected = selected[np.isfinite(selected[metric]) & (selected[metric] > 0.0)]
    medians = selected.groupby("length")[metric].median().sort_index()
    if len(medians) < 3:
        return {
            "fit_id": fit_id, "metric": metric, "n_lengths": int(len(medians)),
            "slope": float("nan"), "intercept": float("nan"),
            "slope_ci_low": float("nan"), "slope_ci_high": float("nan"),
        }
    x = np.log(medians.index.to_numpy(dtype=float))
    y = np.log(medians.to_numpy(dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    seeds = sorted(int(value) for value in selected["generator_seed"].unique())
    rng = np.random.default_rng(int(bootstrap_seed))
    boot = []
    for _ in range(int(bootstrap_replicates)):
        sampled = rng.choice(seeds, size=len(seeds), replace=True)
        pieces = [selected[selected["generator_seed"] == seed] for seed in sampled]
        replicate = pd.concat(pieces, ignore_index=True)
        replicate_medians = replicate.groupby("length")[metric].median().reindex(medians.index)
        if replicate_medians.notna().all() and (replicate_medians > 0.0).all():
            boot.append(float(np.polyfit(x, np.log(replicate_medians.to_numpy(dtype=float)), 1)[0]))
    return {
        "fit_id": fit_id, "metric": metric, "n_lengths": int(len(medians)),
        "lengths": json.dumps([int(value) for value in medians.index]),
        "per_length_medians": json.dumps({str(int(key)): float(value) for key, value in medians.items()}),
        "slope": float(slope), "intercept": float(intercept),
        "slope_ci_low": float(np.quantile(boot, 0.025)) if boot else float("nan"),
        "slope_ci_high": float(np.quantile(boot, 0.975)) if boot else float("nan"),
        "bootstrap_replicates_realized": int(len(boot)),
    }


def make_scale_fits(
    paired: pd.DataFrame,
    mutation: pd.DataFrame,
    epistasis: pd.DataFrame,
    config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    b = int(config["acceptance"]["bootstrap_replicates"])
    seed = int(config["acceptance"]["bootstrap_seed"])
    minimum_n = int(config["acceptance"]["minimum_stratum_observations"])
    point = paired[
        (paired["arm"] == ARM_CANDIDATE)
        & (paired["readout_scope"] == "pointwise_occupancy")
        & (paired["operating_point"] == "p_point_0.5")
        & (paired["weight_label"] == "a_0.5")
        & (paired["sample_id"] == "reference_shell_h1_complete")
    ]
    anchored = mutation[
        (mutation["arm"] == ARM_CANDIDATE)
        & (mutation["operating_point"] == "p_ref_0.5")
        & (mutation["weight_label"] == "a_0.5")
        & (mutation["sample_id"] == "reference_shell_h1_complete")
        & (mutation["occupancy_stratum"] == "both_dynamic")
        & (mutation["n"] >= minimum_n)
    ]
    deployment = mutation[
        (mutation["arm"] == ARM_CANDIDATE)
        & (mutation["operating_point"] == "public_base_margin8")
        & (mutation["weight_label"] == "fixed_0.6")
        & (mutation["sample_id"] == "reference_shell_h1_complete")
        & (mutation["occupancy_stratum"] == "all")
    ]
    texture_frozen = epistasis[
        (epistasis["arm"] == ARM_CANDIDATE)
        & (epistasis["support_order"] > 0)
        & (epistasis["square_background"] == "reference")
        & (epistasis["additive_core_mode"] == "frozen_additive_core")
        & (epistasis["operating_point"] == "p_ref_0.5")
        & (epistasis["weight_label"] == "a_0.5")
    ]
    texture_actual = epistasis[
        (epistasis["arm"] == ARM_CANDIDATE)
        & (epistasis["support_order"] > 0)
        & (epistasis["square_background"] == "reference")
        & (epistasis["additive_core_mode"] == "actual_additive_core")
        & (epistasis["operating_point"] == "p_ref_0.5")
        & (epistasis["weight_label"] == "a_0.5")
    ]
    return [
        log_length_fit(point, metric="stable_y_rms", fit_id="pointwise_p0.5_a0.5_h1", bootstrap_replicates=b, bootstrap_seed=seed),
        log_length_fit(anchored, metric="stable_y_rms", fit_id="reference_anchored_p0.5_a0.5_h1_both_dynamic", bootstrap_replicates=b, bootstrap_seed=seed + 1),
        log_length_fit(deployment, metric="direct_y_rms", fit_id="deployment_margin8_w0.6_h1", bootstrap_replicates=b, bootstrap_seed=seed + 2),
        log_length_fit(texture_frozen, metric="stable_tau_y_rms", fit_id="texture_frozen_p0.5_a0.5", bootstrap_replicates=b, bootstrap_seed=seed + 3),
        log_length_fit(texture_actual, metric="stable_tau_y_rms", fit_id="texture_actual_p0.5_a0.5", bootstrap_replicates=b, bootstrap_seed=seed + 4),
    ]


def bootstrap_median_lower(values: Sequence[float], *, replicates: int, seed: int) -> float:
    array = np.asarray([value for value in values if np.isfinite(value)], dtype=np.float64)
    if array.size == 0:
        return float("nan")
    rng = np.random.default_rng(int(seed))
    boot = np.empty(int(replicates), dtype=np.float64)
    for index in range(int(replicates)):
        boot[index] = float(np.median(rng.choice(array, size=array.size, replace=True)))
    return float(np.quantile(boot, 0.025))


def acceptance_checks(
    *,
    tasks: pd.DataFrame,
    design: pd.DataFrame,
    paired: pd.DataFrame,
    mutation: pd.DataFrame,
    epistasis: pd.DataFrame,
    topology: pd.DataFrame,
    fits: pd.DataFrame,
    config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, bool]]:
    gates = config["acceptance"]
    checks: list[dict[str, Any]] = []

    def add(check_id: str, passed: bool, value: Any, criterion: str, *, family: str, required: bool = True) -> None:
        checks.append({
            "check_id": check_id, "passed": bool(passed), "required": bool(required),
            "family": family, "value": value, "criterion": criterion,
        })

    parity_tol = float(gates["public_formula_parity_atol"])
    add(
        "deterministic_field_replay",
        bool(tasks["candidate_replay_match"].all()),
        int(tasks["candidate_replay_match"].sum()),
        "candidate fingerprint replay matches in every task",
        family="technical",
    )
    max_center = float(tasks["max_centered_reference_abs"].max())
    add(
        "reference_centering",
        max_center <= float(gates["reference_center_atol"]),
        max_center,
        f"maximum centered-arm reference score <= {gates['reference_center_atol']}",
        family="technical",
    )
    max_formula = float(tasks["public_formula_parity_max_abs"].max())
    max_end_to_end = float(tasks["public_end_to_end_parity_max_abs"].max())
    add(
        "public_float64_formula_parity",
        max_formula <= parity_tol and max_end_to_end <= parity_tol,
        json.dumps({"core": max_formula, "end_to_end_256_term": max_end_to_end}),
        f"direct formula and AdaptEnv implementation agree within {parity_tol}",
        family="technical",
    )
    add(
        "finite_outputs",
        bool(tasks["all_outputs_finite"].all()),
        int(tasks["all_outputs_finite"].sum()),
        "all task summaries are finite",
        family="technical",
    )
    identity_max = float(tasks[[
        "identity_point_max_error", "identity_edge_max_error", "identity_square_max_error",
    ]].to_numpy(dtype=float).max())
    add(
        "identity_readout_algebra",
        identity_max <= float(gates["identity_atol"]),
        identity_max,
        f"R=wG, DeltaR=wDeltaG, and epsilon_R=w epsilon_G within {gates['identity_atol']}",
        family="technical",
    )
    occupancy_error = float(tasks["reference_occupancy_calibration_max_abs"].max())
    add(
        "reference_occupancy_calibration",
        occupancy_error <= float(gates["occupancy_atol"]),
        occupancy_error,
        f"reference-anchored core occupancy agrees with its target within {gates['occupancy_atol']}",
        family="technical",
    )
    comparator_required = bool(gates.get("constructor_margin10_comparator_required_for_execution", True))
    comparator_present = bool(tasks["constructor_margin10_comparator_present"].all())
    add(
        "constructor_margin10_comparator_present",
        comparator_present or not comparator_required,
        int(tasks["constructor_margin10_comparator_present"].sum()),
        "margin-10 constructor comparator is present in every task",
        family="technical",
        required=comparator_required,
    )
    sampling = design[design["row_type"].isin(["pair_panel", "square_panel"])]
    sampling_flags = sampling["sampling_key_excludes_treatment"].eq(True)
    replay_flags = sampling["design_replay_match"].eq(True)
    add(
        "treatment_independent_sampling",
        bool(sampling_flags.all() and replay_flags.all()),
        json.dumps({"treatment_independent": int(sampling_flags.sum()), "replay_matches": int(replay_flags.sum())}),
        "every pair/square hash excludes treatment and independently replays",
        family="technical",
    )
    zero_paired = paired[paired["weight_label"] == "a_0"]
    zero_mutation = mutation[mutation["weight_label"] == "a_0"]
    zero_epi = epistasis[epistasis["weight_label"] == "a_0"]
    zero_max = max(
        float(zero_paired["stable_y_rms"].max()),
        float(zero_paired["direct_y_rms"].max()),
        float(zero_mutation["stable_y_rms"].max()),
        float(zero_mutation["direct_y_rms"].max()),
        float(zero_epi["stable_tau_y_rms"].max()),
        float(zero_epi["direct_tau_y_rms"].max()),
    )
    add(
        "weight_zero_exact",
        zero_max <= float(gates["zero_weight_atol"]),
        zero_max,
        f"all paired residual estimands at a=0 are <= {gates['zero_weight_atol']}",
        family="technical",
    )
    surrogate = epistasis[epistasis["arm"] == ARM_ADDITIVE]
    surrogate_max = float(surrogate["latent_epsilon_max_abs"].max())
    surrogate_observable_rows = int(np.sum(surrogate["epsilon_R_stable_rms"] > 0.0))
    add(
        "additive_surrogate_latent_epistasis_zero",
        surrogate_max <= float(gates["additive_surrogate_latent_epistasis_atol"]),
        json.dumps({"max_latent_abs": surrogate_max, "observable_nonzero_rows": surrogate_observable_rows}),
        "order-1 surrogate latent epsilon is zero and observable epsilon is explicitly reported",
        family="technical",
    )
    minimum_topology = float(topology["non_tie_rank_agreement"].min())
    mismatch_count = int(topology["non_tie_mismatch_count"].sum())
    add(
        "strict_monotone_non_tie_topology",
        minimum_topology >= float(gates["topology_non_tie_agreement_min"]) and mismatch_count == 0,
        json.dumps({"minimum_agreement": minimum_topology, "mismatches": mismatch_count}),
        "zero rank mismatches among direct float64 non-ties",
        family="technical",
    )

    point = paired[
        (paired["arm"] == ARM_CANDIDATE)
        & (paired["sample_id"] == "reference_shell_h1_complete")
        & (paired["operating_point"] == "p_point_0.5")
        & (paired["weight_label"] == "a_0.5")
    ]
    gamma_by_length = point.groupby("length")["gamma_median"].median()
    resolution_by_length = point.groupby("length")["stable_practical_resolution_fraction"].median()
    point_linear_corr = float(point["linear_correlation"].median())
    point_linear_nrmse = float(point["linear_nrmse"].median())
    shoulder = paired[
        (paired["arm"] == ARM_CANDIDATE)
        & (paired["sample_id"] == "reference_shell_h1_complete")
        & (paired["operating_point"].isin(["p_point_0.1", "p_point_0.9"]))
        & (paired["weight_label"] == "a_0.5")
    ]
    shoulder_min = float(shoulder.groupby("operating_point")["gamma_median"].median().min())
    point_fit = fits[fits["fit_id"] == "pointwise_p0.5_a0.5_h1"].iloc[0]
    point_pass = bool(
        len(gamma_by_length) >= 3
        and (gamma_by_length >= float(gates["pointwise_gamma_median_min"])).all()
        and (resolution_by_length >= float(gates["pointwise_detectable_fraction_min"])).all()
        and shoulder_min >= float(gates["shoulder_gamma_median_min"])
        and point_linear_corr >= float(gates["linearization_correlation_min"])
        and point_linear_nrmse <= float(gates["linearization_nrmse_max"])
        and abs(float(point_fit["slope"])) <= float(gates["local_slope_abs_max"])
        and float(point_fit["slope_ci_low"]) >= -float(gates["local_slope_ci_abs_max"])
        and float(point_fit["slope_ci_high"]) <= float(gates["local_slope_ci_abs_max"])
    )
    add(
        "pointwise_transfer_calibration",
        point_pass,
        json.dumps({
            "minimum_gamma": float(gamma_by_length.min()) if len(gamma_by_length) else None,
            "minimum_practical_fraction": float(resolution_by_length.min()) if len(resolution_by_length) else None,
            "shoulder_gamma": shoulder_min, "linear_correlation": point_linear_corr,
            "linear_nrmse": point_linear_nrmse, "slope": float(point_fit["slope"]),
            "slope_ci": [float(point_fit["slope_ci_low"]), float(point_fit["slope_ci_high"])],
        }),
        "technical scalar-transfer calibration passes all pointwise bounds",
        family="technical",
    )

    controlled_edges_all = mutation[
        (mutation["arm"] == ARM_CANDIDATE)
        & (mutation["readout_scope"] == "reference_anchored")
        & (mutation["operating_point"] == "p_ref_0.5")
        & (mutation["weight_label"] == "a_0.5")
        & (mutation["sample_id"] == "reference_shell_h1_complete")
        & (mutation["occupancy_stratum"] == "both_dynamic")
    ]
    minimum_stratum_n = int(gates["minimum_stratum_observations"])
    controlled_edges = controlled_edges_all[
        controlled_edges_all["n"] >= minimum_stratum_n
    ]
    expected_index = pd.MultiIndex.from_frame(
        tasks[["length", "generator_seed"]].drop_duplicates().sort_values(["length", "generator_seed"])
    )
    observed_counts = controlled_edges.groupby(["length", "generator_seed"])["n"].sum()
    complete_counts = observed_counts.reindex(expected_index, fill_value=0)
    eligible_counts = complete_counts.groupby(level="length").sum()
    seed_coverage_by_length = (complete_counts > 0).groupby(level="length").mean()
    eligible_seed_fraction = float(seed_coverage_by_length.min()) if len(seed_coverage_by_length) else 0.0
    controlled_corr_by_length = controlled_edges.groupby("length")["stable_correlation"].median()
    controlled_resolution_by_length = controlled_edges.groupby("length")["stable_practical_resolution_fraction"].median()
    controlled_corr = float(controlled_corr_by_length.min()) if len(controlled_corr_by_length) else float("nan")
    controlled_resolution = float(controlled_resolution_by_length.min()) if len(controlled_resolution_by_length) else float("nan")
    anchored_fit = fits[fits["fit_id"] == "reference_anchored_p0.5_a0.5_h1_both_dynamic"].iloc[0]
    anchored_scale_pass = bool(
        np.isfinite(float(anchored_fit["slope"]))
        and abs(float(anchored_fit["slope"])) <= float(gates["local_slope_abs_max"])
        and float(anchored_fit["slope_ci_low"]) >= -float(gates["local_slope_ci_abs_max"])
        and float(anchored_fit["slope_ci_high"]) <= float(gates["local_slope_ci_abs_max"])
    )
    controlled_edge_pass = bool(
        len(eligible_counts) > 0
        and not controlled_edges.empty
        and int(controlled_edges["n"].min()) >= minimum_stratum_n
        and eligible_seed_fraction >= float(gates["minimum_stratum_seed_coverage_fraction_per_length"])
        and controlled_corr >= float(gates["transmission_correlation_min"])
        and controlled_resolution >= float(gates["texture_detectable_fraction_min"])
        and anchored_scale_pass
    )
    add(
        "reference_anchored_edge_transmission",
        controlled_edge_pass,
        json.dumps({
            "minimum_retained_row_n": int(controlled_edges["n"].min()) if not controlled_edges.empty else 0,
            "minimum_pooled_retained_n_per_length": int(eligible_counts.min()) if len(eligible_counts) else 0,
            "raw_rows_below_minimum_n": int(np.sum(controlled_edges_all["n"] < minimum_stratum_n)),
            "eligible_seed_fraction": eligible_seed_fraction,
            "minimum_per_length_median_correlation": controlled_corr,
            "minimum_per_length_median_practical_fraction": controlled_resolution,
            "length_slope": float(anchored_fit["slope"]),
            "length_slope_ci": [float(anchored_fit["slope_ci_low"]), float(anchored_fit["slope_ci_high"])],
        }),
        "common-threshold h=1 edges retain >=32-observation seed strata with adequate coverage, correlation, resolution, and scale stability",
        family="controlled_transfer",
    )

    def texture_summary(core_mode: str, deployment: bool) -> tuple[bool, dict[str, Any]]:
        frame = epistasis[
            (epistasis["arm"] == ARM_CANDIDATE)
            & (epistasis["support_order"] > 0)
            & (epistasis["additive_core_mode"] == core_mode)
            & (epistasis["square_background"] == "reference")
        ]
        if deployment:
            frame = frame[
                (frame["operating_point"] == "public_base_margin8")
                & (frame["weight_label"] == "fixed_0.6")
            ]
            prefix = "direct_tau_"
        else:
            frame = frame[
                (frame["operating_point"] == "p_ref_0.5")
                & (frame["weight_label"] == "a_0.5")
            ]
            prefix = "stable_tau_"
        if frame.empty:
            return False, {"rows": 0}
        expected_lengths = sorted(int(value) for value in tasks["length"].unique())
        seed_projection = (
            frame.groupby(["length", "generator_seed"])[f"{prefix}projection"]
            .median()
            .reindex(expected_index)
        )
        seed_coverage_by_length = (
            seed_projection.notna().groupby(level="length").mean()
            .reindex(expected_lengths, fill_value=0.0)
        )
        positive_fraction_by_length = (
            seed_projection.fillna(-np.inf).gt(0.0).groupby(level="length").mean()
            .reindex(expected_lengths, fill_value=0.0)
        )
        ci_lower_by_length: dict[int, float] = {}
        bootstrap_offset = (11 if deployment else 7) + (1 if core_mode.startswith("frozen") else 0)
        for length_index, length in enumerate(expected_lengths):
            values = seed_projection.xs(length, level="length").dropna().to_numpy(dtype=float)
            ci_lower_by_length[length] = bootstrap_median_lower(
                values,
                replicates=int(gates["bootstrap_replicates"]),
                seed=int(gates["bootstrap_seed"]) + bootstrap_offset + 101 * length_index,
            )
        if deployment:
            practical_frame = frame
        else:
            practical_frame = epistasis[
                (epistasis["arm"] == ARM_CANDIDATE)
                & (epistasis["support_order"] > 0)
                & (epistasis["additive_core_mode"] == core_mode)
                & (epistasis["square_background"] == "reference")
                & (epistasis["operating_point"] == "p_ref_0.5")
                & (epistasis["weight_label"] == "a_1")
            ]
        practical_by_length = (
            practical_frame.groupby("length")[f"{prefix}practical_resolution_fraction"].median()
            .reindex(expected_lengths)
            if not practical_frame.empty else pd.Series(index=expected_lengths, dtype=float)
        )
        predictor_energy_by_length = (
            frame.assign(
                _predictor_energy=frame[f"{prefix}sum_x2"] / frame[f"{prefix}n"].clip(lower=1)
            )
            .groupby("length")["_predictor_energy"].min()
            .reindex(expected_lengths)
        )
        minimum_n_by_length = frame.groupby("length")["n"].min().reindex(expected_lengths)
        permutation_required = frame[frame["permutation_required_primary"].eq(True)]
        permutation_fraction_by_length = (
            permutation_required.assign(_passed=permutation_required["permutation_pass"].eq(True))
            .groupby("length")["_passed"].mean()
            .reindex(expected_lengths, fill_value=0.0)
            if not permutation_required.empty
            else pd.Series(0.0, index=expected_lengths, dtype=float)
        )
        linear_corr_by_length = (
            frame.groupby("length")["linear_correlation"].median().reindex(expected_lengths)
        )
        projection_threshold = float(
            gates[
                "deployment_projection_positive_seed_fraction_min"
                if deployment else "transmission_projection_positive_seed_fraction_min"
            ]
        )
        all_per_length_finite = bool(
            practical_by_length.notna().all()
            and predictor_energy_by_length.notna().all()
            and minimum_n_by_length.notna().all()
            and (deployment or linear_corr_by_length.notna().all())
            and all(np.isfinite(value) for value in ci_lower_by_length.values())
        )
        if deployment:
            scale_pass = True
            scale_values: dict[str, Any] = {"length_slope_required": False}
        else:
            fit_id = (
                "texture_frozen_p0.5_a0.5"
                if core_mode == "frozen_additive_core"
                else "texture_actual_p0.5_a0.5"
            )
            matching_fit = fits[fits["fit_id"] == fit_id]
            if len(matching_fit) == 1:
                texture_fit = matching_fit.iloc[0]
                scale_pass = bool(
                    int(texture_fit["n_lengths"]) == len(expected_lengths)
                    and np.isfinite(float(texture_fit["slope"]))
                    and abs(float(texture_fit["slope"])) <= float(gates["local_slope_abs_max"])
                    and float(texture_fit["slope_ci_low"]) >= -float(gates["local_slope_ci_abs_max"])
                    and float(texture_fit["slope_ci_high"]) <= float(gates["local_slope_ci_abs_max"])
                )
                scale_values = {
                    "length_slope_required": True,
                    "length_slope": float(texture_fit["slope"]),
                    "length_slope_ci": [
                        float(texture_fit["slope_ci_low"]),
                        float(texture_fit["slope_ci_high"]),
                    ],
                    "n_lengths_in_fit": int(texture_fit["n_lengths"]),
                }
            else:
                scale_pass = False
                scale_values = {"length_slope_required": True, "fit_missing": fit_id}
        passed = bool(
            all_per_length_finite
            and len(expected_lengths) > 0
            and float(seed_coverage_by_length.min()) >= float(gates["minimum_stratum_seed_coverage_fraction_per_length"])
            and int(minimum_n_by_length.min()) >= int(gates["minimum_nonzero_texture_count"])
            and float(predictor_energy_by_length.min()) >= float(gates["minimum_texture_predictor_variance"])
            and float(positive_fraction_by_length.min()) >= projection_threshold
            and min(ci_lower_by_length.values()) > float(gates["texture_projection_bootstrap_lower_min"])
            and float(practical_by_length.min()) >= float(gates["texture_detectable_fraction_min"])
            and float(permutation_fraction_by_length.min()) >= float(gates["texture_permutation_pass_fraction_min"])
            and (deployment or float(linear_corr_by_length.min()) >= float(gates["transmission_correlation_min"]))
            and scale_pass
        )
        return passed, {
            "rows": int(len(frame)),
            "minimum_seed_coverage": float(seed_coverage_by_length.min()),
            "minimum_per_length_positive_seed_fraction": float(positive_fraction_by_length.min()),
            "minimum_per_length_projection_bootstrap_lower": float(min(ci_lower_by_length.values())),
            "minimum_per_length_practical_fraction": float(practical_by_length.min()),
            "minimum_per_length_permutation_pass_fraction": float(permutation_fraction_by_length.min()),
            "minimum_per_length_median_linear_correlation": (
                None if deployment else float(linear_corr_by_length.min())
            ),
            "minimum_per_length_predictor_energy": float(predictor_energy_by_length.min()),
            "minimum_per_length_n": int(minimum_n_by_length.min()),
            "per_length_positive_seed_fraction": {
                str(int(key)): float(value) for key, value in positive_fraction_by_length.items()
            },
            "per_length_projection_bootstrap_lower": {
                str(int(key)): float(value) for key, value in ci_lower_by_length.items()
            },
            "per_length_practical_fraction": {
                str(int(key)): float(value) for key, value in practical_by_length.items()
            },
            "per_length_permutation_pass_fraction": {
                str(int(key)): float(value) for key, value in permutation_fraction_by_length.items()
            },
            **scale_values,
        }

    frozen_pass, frozen_values = texture_summary("frozen_additive_core", False)
    actual_pass, actual_values = texture_summary("actual_additive_core", False)
    add(
        "frozen_core_texture_transmission",
        frozen_pass,
        json.dumps(frozen_values),
        "candidate-additivized tau_texture passes in the frozen-core mechanistic control",
        family="controlled_texture",
    )
    add(
        "actual_core_texture_transmission",
        actual_pass,
        json.dumps(actual_values),
        "candidate-additivized tau_texture passes with the real additive core",
        family="controlled_texture",
    )

    public_edges = mutation[
        (mutation["arm"] == ARM_CANDIDATE)
        & (mutation["operating_point"] == "public_base_margin8")
        & (mutation["weight_label"] == "fixed_0.6")
        & (mutation["sample_id"] == "reference_shell_h1_complete")
        & (mutation["occupancy_stratum"] == "all")
    ]
    public_edge_fraction_by_length = (
        public_edges.groupby("length")["direct_practical_resolution_fraction"]
        .median()
        .reindex(sorted(int(value) for value in tasks["length"].unique()))
    )
    public_edge_fraction = float(public_edge_fraction_by_length.min())
    public_edge_pass = bool(
        public_edge_fraction_by_length.notna().all()
        and public_edge_fraction >= float(gates["deployment_detectable_fraction_min"])
    )
    public_texture_pass, public_texture_values = texture_summary("actual_additive_core", True)
    public_pass = bool(public_edge_pass and public_texture_pass)
    add(
        "margin8_direct_float64_deployment",
        public_pass,
        json.dumps({
            "minimum_per_length_edge_practical_fraction": public_edge_fraction,
            "per_length_edge_practical_fraction": {
                str(int(key)): float(value) for key, value in public_edge_fraction_by_length.items()
            },
            **public_texture_values,
        }),
        "direct float64 margin-8 edges and candidate-additivized texture pass every per-length deployment bound",
        family="deployment",
    )
    return checks, {
        "technical": bool(all(row["passed"] for row in checks if row["family"] == "technical")),
        "controlled_edges": controlled_edge_pass,
        "frozen_texture": frozen_pass,
        "actual_texture": actual_pass,
        "public": public_pass,
    }


def classify_decision(flags: Mapping[str, bool], config: Mapping[str, Any]) -> tuple[bool, str]:
    labels = config["decision_rule"]
    if not bool(flags["technical"]):
        return False, str(labels["technical_fail"])
    if not bool(flags["controlled_edges"]):
        return False, str(labels["transfer_fail"])
    if bool(flags["frozen_texture"]) and not bool(flags["actual_texture"]):
        return False, str(labels["frozen_core_only"])
    if not bool(flags["actual_texture"]):
        return False, str(labels["texture_fail"])
    if bool(flags["public"]):
        return True, str(labels["public_pass"])
    return False, str(labels["transition_only"])


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, quick=bool(args.quick))
    validate_config(config)
    references, reference_source_sha256 = load_references(project_root, config)
    seeds = [int(value) for value in config["replication"]["seeds"]]
    if len(seeds) != len(set(seeds)):
        raise ValueError("generator seeds must be unique")
    implementation_hashes = source_hashes(project_root, config)
    implementation_hash = latent.json_hash(implementation_hashes)
    scientific_decision_enabled = not bool(config["experiment"].get("disable_scientific_decision", False))
    preflight = {
        "experiment_id": "EXP-061",
        "quick": bool(args.quick),
        "scientific_decision_enabled": scientific_decision_enabled,
        "config_file_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "effective_config_hash": latent.json_hash(config),
        "reference_source_sha256": reference_source_sha256,
        "implementation_hash": implementation_hash,
        "implementation_hash_components": implementation_hashes,
        "reference_count": len(references),
        "lengths": [int(row["length"]) for row in references],
        "generator_seed_count": len(seeds),
        "task_count": len(references) * len(seeds),
        "arms": list(ARMS),
        "readout_scopes": ["pointwise_occupancy", "reference_anchored", "deployment"],
        "pointwise_scope_is_one_landscape": False,
        "deployment_primary": "public_base_margin8_direct_float64",
        "deployment_comparator": "constructor_default_margin10_direct_float64",
        "stable_deployment_residual_can_satisfy_gate": False,
        "observation_noise_enabled": False,
        "candidate_replaces_public_generator_in_primary_arm": True,
        "public_model_changed": False,
        "candidate_is_experiment_local": True,
        "stage_4c_4d_4e_hard_stop_retained": True,
    }
    latent.write_json_atomic(output_dir / "preflight.json", preflight)

    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    payloads: list[dict[str, Any]] = []
    for reference in references:
        for generator_seed in seeds:
            signature = latent.json_hash({
                "config": config,
                "reference": reference,
                "generator_seed": generator_seed,
                "reference_source_sha256": reference_source_sha256,
                "implementation_hash": implementation_hash,
            })
            payloads.append({
                "config": config,
                "reference": reference,
                "generator_seed": generator_seed,
                "task_signature": signature,
                "checkpoint_path": str(
                    checkpoint_dir
                    / f"{latent.safe_name(reference['reference_id'])}__L-{reference['length']}__seed-{generator_seed}__{signature[:16]}.json"
                ),
            })

    workers = max(1, int(os.environ.get("LABPROJ_WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "1"))))
    results: list[dict[str, Any]] = []
    if workers == 1:
        for index, payload in enumerate(payloads, start=1):
            results.append(run_task(payload))
            print(f"[EXP-061] completed {index}/{len(payloads)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(run_task, payload): payload for payload in payloads}
            for index, future in enumerate(as_completed(futures), start=1):
                results.append(future.result())
                print(f"[EXP-061] completed {index}/{len(payloads)}", flush=True)
    results.sort(key=lambda result: (
        int(result["task"]["length"]),
        str(result["task"]["reference_id"]),
        int(result["task"]["generator_seed"]),
    ))

    task_rows = [result["task"] for result in results]
    design_rows = [row for result in results for row in result["design_rows"]]
    occupancy_rows = [row for result in results for row in result["occupancy_rows"]]
    paired_rows = [row for result in results for row in result["paired_rows"]]
    mutation_rows = [row for result in results for row in result["mutation_rows"]]
    epistasis_rows = [row for result in results for row in result["epistasis_rows"]]
    topology_rows = [row for result in results for row in result["topology_rows"]]
    tasks = pd.DataFrame(task_rows)
    design = pd.DataFrame(design_rows)
    paired = pd.DataFrame(paired_rows)
    mutation = pd.DataFrame(mutation_rows)
    epistasis = pd.DataFrame(epistasis_rows)
    topology = pd.DataFrame(topology_rows)
    fit_rows = make_scale_fits(paired, mutation, epistasis, config)
    fits = pd.DataFrame(fit_rows)
    checks, flags = acceptance_checks(
        tasks=tasks,
        design=design,
        paired=paired,
        mutation=mutation,
        epistasis=epistasis,
        topology=topology,
        fits=fits,
        config=config,
    )
    if scientific_decision_enabled:
        scientific_pass, decision = classify_decision(flags, config)
    else:
        scientific_pass = None
        decision = str(config["decision_rule"]["quick"])

    latent.write_csv_atomic(output_dir / "design_catalog.csv", design_rows)
    latent.write_csv_atomic(output_dir / "occupancy_calibration.csv", occupancy_rows)
    latent.write_csv_atomic(output_dir / "paired_residual_metrics.csv", paired_rows)
    latent.write_csv_atomic(output_dir / "mutation_transmission_metrics.csv", mutation_rows)
    latent.write_csv_atomic(output_dir / "epistasis_transmission_metrics.csv", epistasis_rows)
    latent.write_csv_atomic(output_dir / "topology_audit.csv", topology_rows)
    latent.write_csv_atomic(output_dir / "scale_law_fits.csv", fit_rows)
    latent.write_csv_atomic(output_dir / "acceptance_checks.csv", checks)
    required = [row for row in checks if bool(row["required"])]
    summary = {
        **preflight,
        "completed_task_count": len(task_rows),
        "design_row_count": len(design_rows),
        "paired_residual_row_count": len(paired_rows),
        "mutation_transmission_row_count": len(mutation_rows),
        "epistasis_transmission_row_count": len(epistasis_rows),
        "topology_audit_row_count": len(topology_rows),
        "required_checks_passed": int(sum(bool(row["passed"]) for row in required)),
        "required_check_count": len(required),
        "decision_flags": flags,
        "scientific_pass": scientific_pass,
        "decision": decision,
        "public_model_changed": False,
        "empirical_realism_tested": False,
        "optimizer_or_accessibility_tested": False,
        "stable_residual_is_deployment_gate": False,
    }
    latent.write_json_atomic(output_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
