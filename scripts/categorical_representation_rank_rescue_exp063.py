#!/usr/bin/env python3
"""EXP-063 paired categorical-code and representation-rank intervention.

The runner freezes the EXP-062 additive coordinate and sigmoid, reuses its
fixed failure-challenge protein panel, and changes only the categorical
representation carried by the extensive GNK supports.  All four arms are
evaluated on identical supports and probe genotypes:

* binary_r1
* dense_q20_r1
* binary_r8
* dense_q20_r8

The public AdaptEnv implementation is not modified.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
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

import dimensionless_sigmoid_calibration_exp062 as exp062  # noqa: E402
import extensive_gnk_latent_scaling_exp060 as latent  # noqa: E402
from project_adapt_env.extensive_gnk import (  # noqa: E402
    ExtensiveGNKField,
    additive_mutation_energy,
    stable_uint64,
)
from project_adapt_env.ranked_categorical_gnk import (  # noqa: E402
    RankedCategoricalGNKField,
)


EXP062_CONFIG = PROJECT_ROOT / (
    "experiments/2026-07-23_EXP-062_dimensionless-additive-core-sigmoid-"
    "calibration-and-held-out-gnk-dose-response-validation/config.yaml"
)
ARM_COMPLEXITY = (
    "binary_r1",
    "dense_q20_r1",
    "binary_r8",
    "dense_q20_r8",
)
RESCUE_ARMS = ARM_COMPLEXITY[1:]


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


def helper_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Create the EXP-062-compatible subset used by its audited helpers."""

    canonical = exp062.load_config(EXP062_CONFIG, quick=False)
    out = copy.deepcopy(canonical)
    out["experiment_id"] = "EXP-063"
    out["_runtime"] = dict(config.get("_runtime", {}))
    out["additive_core"] = copy.deepcopy(config["additive_core"])
    out["candidate"].update({
        key: copy.deepcopy(config["candidate"][key])
        for key in (
            "seed_offset",
            "surrogate_seed_offset",
            "orders",
            "order_weights",
            "target_degrees",
        )
    })
    out["intervention"]["dimensionless_amplitudes"] = copy.deepcopy(
        config["intervention"]["dimensionless_amplitudes"]
    )
    out["intervention"]["directionality_amplitudes"] = copy.deepcopy(
        config["intervention"]["directionality_amplitudes"]
    )
    for key in (
        "square_strata",
        "square_backgrounds",
        "squares_per_stratum",
        "eval_chunk_size",
    ):
        out["probe"][key] = copy.deepcopy(config["probe"][key])
    out["acceptance"]["expected_length_bins"] = int(
        config["acceptance"]["expected_length_bins"]
    )
    for key in (
        "bootstrap_replicates",
        "bootstrap_seed",
        "permutation_replicates",
        "texture_permutation_alpha",
    ):
        out["acceptance"][key] = copy.deepcopy(config["acceptance"][key])
    for section in ("core", "dose_response", "length_stability"):
        out["acceptance"][section].update(
            copy.deepcopy(config["acceptance"][section])
        )
    # The inherited EXP-062 helper performs a direct >= comparison.  Shift its
    # private threshold so EXP-063 implements the registered tolerance as
    # rho + atol >= spearman_min without modifying the prior experiment.
    out["acceptance"]["dose_response"]["spearman_min"] = (
        float(config["acceptance"]["dose_response"]["spearman_min"])
        - float(config["acceptance"]["dose_response"]["rho_comparison_atol"])
    )
    out["acceptance"]["transmission"].update(
        copy.deepcopy(config["acceptance"]["transmission"])
    )
    out["readout"]["output_scale"] = float(
        config["frozen_operator"]["output_scale"]
    )
    return out


def validate_config(config: Mapping[str, Any]) -> None:
    required = {
        "experiment",
        "references",
        "replication",
        "additive_core",
        "frozen_operator",
        "candidate",
        "readout",
        "intervention",
        "probe",
        "acceptance",
        "decision",
        "outputs",
    }
    missing = sorted(required - set(config))
    if missing:
        raise ValueError(f"config missing sections: {missing}")
    if str(config.get("experiment_id")) != "EXP-063":
        raise ValueError("EXP-063 runner requires experiment_id EXP-063")
    arms = [str(row["id"]) for row in config["candidate"]["arms"]]
    if tuple(arms) != ARM_COMPLEXITY:
        raise ValueError(f"candidate arms must be {ARM_COMPLEXITY}")
    for row in config["candidate"]["arms"]:
        family = str(row["contrast_family"])
        rank = int(row["representation_rank"])
        expected_family, expected_rank = {
            "binary_r1": ("binary", 1),
            "dense_q20_r1": ("dense_q20", 1),
            "binary_r8": ("binary", 8),
            "dense_q20_r8": ("dense_q20", 8),
        }[str(row["id"])]
        if (family, rank) != (expected_family, expected_rank):
            raise ValueError(f"incorrect representation specification for {row['id']}")
    frozen = config["frozen_operator"]
    if not bool(frozen.get("recalibration_forbidden")):
        raise ValueError("the EXP-062 operating point must remain frozen")
    if not math.isclose(
        float(frozen["mu"]),
        math.log(float(frozen["reference_occupancy"]) / (1.0 - float(frozen["reference_occupancy"]))),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise ValueError("frozen mu and reference occupancy disagree")
    amplitudes = [float(value) for value in config["intervention"]["dimensionless_amplitudes"]]
    if 0.0 not in amplitudes or 1.0 not in amplitudes:
        raise ValueError("amplitude grid must contain 0 and 1")
    if not bool(config.get("_runtime", {}).get("quick")):
        required_amplitudes = {0.0, 0.25, 0.5, 1.0, 2.0, 4.0}
        if set(amplitudes) != required_amplitudes:
            raise ValueError("full run requires the preregistered amplitude grid")
    core_seeds = {int(value) for value in config["replication"]["core_seeds"]}
    gnk_seeds = {int(value) for value in config["replication"]["gnk_seeds"]}
    if not core_seeds or not gnk_seeds or core_seeds & gnk_seeds:
        raise ValueError("core and GNK seeds must be nonempty and disjoint")


def source_hashes(
    project_root: Path, config: Mapping[str, Any]
) -> dict[str, str]:
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
    catalog_path = project_root / str(config["references"]["catalog"])
    payload = yaml.safe_load(catalog_path.read_text(encoding="utf-8")) or {}
    selected = {str(value) for value in config["references"].get("include_ids", [])}
    required_split = str(config["references"]["split"])
    rows: list[dict[str, Any]] = []
    for raw in payload.get("references", []):
        if str(raw["split"]) != required_split:
            continue
        if selected and str(raw["id"]) not in selected:
            continue
        sequence = str(raw["sequence"]).strip().upper()
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if digest != str(raw["sequence_sha256"]):
            raise ValueError(f"reference hash mismatch: {raw['id']}")
        if len(sequence) != int(raw["length"]):
            raise ValueError(f"reference length mismatch: {raw['id']}")
        if not set(sequence).issubset(set("ACDEFGHIKLMNPQRSTVWY")):
            raise ValueError(f"noncanonical reference: {raw['id']}")
        rows.append({
            "reference_id": str(raw["id"]),
            "split": str(raw["split"]),
            "length_bin": str(raw["length_bin"]),
            "length": int(raw["length"]),
            "sequence": sequence,
            "sequence_sha256": digest,
            "source_dataset": str(raw["source_dataset"]),
            "source_record": str(raw["source_record"]),
            "selection_rule": str(raw["selection_rule"]),
        })
    if not rows:
        raise ValueError("reference filtering produced an empty panel")
    counts = pd.DataFrame(rows).groupby("length_bin").size()
    expected_bins = int(config["references"]["expected_length_bins"])
    expected_per_bin = int(config["references"]["expected_references_per_length_bin"])
    if len(counts) != expected_bins or not bool((counts == expected_per_bin).all()):
        raise ValueError(
            "reference panel does not match expected length-bin replication"
        )
    if selected != {row["reference_id"] for row in rows}:
        raise ValueError("one or more configured reference IDs were not loaded")
    return (
        sorted(rows, key=lambda row: (row["length_bin"], row["reference_id"])),
        hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
    )


def _clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean_json(item) for item in value]
    if isinstance(value, np.generic):
        return _clean_json(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def strict_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_clean_json(dict(payload)), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def parallel_execute(
    function: Callable[[Mapping[str, Any]], Any],
    payloads: Sequence[Mapping[str, Any]],
    *,
    workers: int,
    label: str,
) -> list[Any]:
    results: list[Any] = []
    total = len(payloads)
    if workers <= 1:
        for index, payload in enumerate(payloads, start=1):
            results.append(function(payload))
            print(f"[EXP-063] {label} {index}/{total}", flush=True)
        return results
    with ProcessPoolExecutor(max_workers=int(workers)) as executor:
        futures = {
            executor.submit(function, payload): payload for payload in payloads
        }
        for index, future in enumerate(as_completed(futures), start=1):
            results.append(future.result())
            print(f"[EXP-063] {label} {index}/{total}", flush=True)
    return results


def core_guard_task(payload: Mapping[str, Any]) -> dict[str, Any]:
    reference_row = dict(payload["reference"])
    core_seed = int(payload["core_seed"])
    config = payload["helper_config"]
    frozen = payload["frozen_operator"]
    reference = latent.sequence_to_indices(str(reference_row["sequence"]))
    core = exp062.build_core(reference, core_seed, config)
    matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    reference_score = float(latent.additive_scores(matrix, reference[None, :])[0])
    metrics = exp062.core_metrics_for_scale(
        matrix=matrix,
        reference=reference,
        reference_score=reference_score,
        scale_sigma_a=float(frozen["global_sigma_a"]),
        beta=float(frozen["beta"]),
        mu=float(frozen["mu"]),
        output_scale=float(frozen["output_scale"]),
        probe_seed=stable_uint64(
            "EXP-063-core-guard",
            reference_row["reference_id"],
            core_seed,
        ),
        config=config,
    )
    return {
        "reference_id": reference_row["reference_id"],
        "length_bin": reference_row["length_bin"],
        "length": int(reference.size),
        "core_seed": core_seed,
        "global_sigma_a": float(frozen["global_sigma_a"]),
        "task_sigma_a": exp062.complete_h1_rms(matrix, reference),
        "additive_reference_score": reference_score,
        "additive_matrix_hash": latent.array_hash(matrix),
        **metrics,
    }


def core_guard_pass(
    rows: Sequence[Mapping[str, Any]], config: Mapping[str, Any]
) -> tuple[bool, dict[str, Any]]:
    frame = pd.DataFrame(rows)
    if frame.empty:
        return False, {"n_cells": 0, "per_length_bin": {}}
    fractions = frame.groupby("length_bin")["core_pass"].mean().sort_index()
    threshold = float(
        config["acceptance"]["core"]["minimum_cell_pass_fraction_per_length_bin"]
    )
    values = {str(key): float(value) for key, value in fractions.items()}
    passed = bool(
        len(fractions) == int(config["acceptance"]["expected_length_bins"])
        and float(fractions.min()) >= threshold
    )
    return passed, {
        "n_cells": int(len(frame)),
        "minimum_length_bin_pass_fraction": float(fractions.min()),
        "per_length_bin_pass_fraction": values,
        "threshold": threshold,
    }


def arm_spec(config: Mapping[str, Any], arm: str) -> tuple[str, int]:
    for row in config["candidate"]["arms"]:
        if str(row["id"]) == str(arm):
            return str(row["contrast_family"]), int(row["representation_rank"])
    raise KeyError(arm)


def build_field(
    *,
    matrix: np.ndarray,
    reference: np.ndarray,
    field_seed: int,
    arm: str,
    config: Mapping[str, Any],
) -> RankedCategoricalGNKField:
    family, rank = arm_spec(config, arm)
    candidate = config["candidate"]
    return RankedCategoricalGNKField(
        additive_matrix=matrix,
        reference_sequence=reference,
        seed=int(field_seed),
        orders=tuple(int(value) for value in candidate["orders"]),
        order_weights=tuple(float(value) for value in candidate["order_weights"]),
        target_degrees=tuple(float(value) for value in candidate["target_degrees"]),
        contrast_family=family,
        representation_rank=rank,
    )


def _component_zero_coefficient_hash(
    field: RankedCategoricalGNKField,
) -> str:
    digest = hashlib.sha256()
    for term in field.terms:
        if int(term.component_id) != 0:
            continue
        digest.update(str(term.term_id).encode("ascii"))
        digest.update(np.asarray([term.base_coefficient], dtype=np.float64).tobytes())
    return digest.hexdigest()


def _pair_pool_hash(field: RankedCategoricalGNKField) -> str:
    return latent.json_hash({
        str(order): [list(pair) for pair in field.support_pairs(order)]
        for order in field.orders
    })


def representation_audit(
    field: RankedCategoricalGNKField,
    *,
    arm: str,
    ids: Mapping[str, Any],
    rank_relative_threshold: float,
) -> dict[str, Any]:
    energy_errors: dict[str, float] = {}
    effective_ranks: list[int] = []
    minimum_relative_eigenvalues: list[float] = []
    participation_ratios: list[float] = []
    for order in field.orders:
        target = float(
            field.order_weights[order] * field.additive_mutation_energy
        )
        realized = float(field.expected_order_mutation_energy[order])
        energy_errors[str(order)] = abs(realized - target) / max(abs(target), 1.0e-300)
        for support_id in field.support_ids(order)[: min(24, len(field.support_ids(order)))]:
            gram = field.support_operator_gram(support_id)
            eigenvalues = np.maximum(
                np.linalg.eigvalsh(0.5 * (gram + gram.T)),
                0.0,
            )
            maximum = float(np.max(eigenvalues))
            tolerance = max(
                np.finfo(np.float64).eps * max(1, len(eigenvalues)) * maximum,
                float(rank_relative_threshold) * maximum,
            )
            effective_ranks.append(int(np.sum(eigenvalues > tolerance)))
            if maximum > 0.0:
                minimum_relative_eigenvalues.append(
                    float(np.min(eigenvalues) / maximum)
                )
            denominator = float(np.dot(eigenvalues, eigenvalues))
            participation_ratios.append(
                float(np.sum(eigenvalues) ** 2 / denominator)
                if denominator > 0.0
                else 0.0
            )

    sampled_terms = []
    for order in field.orders:
        sampled_terms.extend(
            [term for term in field.terms if term.order == order][:32]
        )
    maximum_mean = 0.0
    maximum_second_error = 0.0
    tie_numerator = 0
    tie_denominator = 0
    binary_balanced = True
    for term in sampled_terms:
        for contrast in term.contrasts:
            values = np.asarray(contrast, dtype=np.float64)
            maximum_mean = max(maximum_mean, abs(float(np.mean(values))))
            maximum_second_error = max(
                maximum_second_error,
                abs(float(np.mean(values * values)) - 1.0),
            )
            _, counts = np.unique(values, return_counts=True)
            tie_numerator += int(np.sum(counts * (counts - 1)))
            tie_denominator += int(values.size * (values.size - 1))
            if field.contrast_family == "binary":
                unique, binary_counts = np.unique(values, return_counts=True)
                binary_balanced = bool(
                    binary_balanced
                    and len(unique) == 2
                    and sorted(binary_counts.tolist()) == [10, 10]
                )
    tie_fraction = (
        float(tie_numerator / tie_denominator) if tie_denominator else float("nan")
    )
    return {
        **ids,
        "representation_arm": arm,
        "contrast_family": field.contrast_family,
        "representation_rank": int(field.representation_rank),
        "support_count": int(sum(field.term_counts.values())),
        "component_count": int(sum(field.component_counts.values())),
        "support_fingerprint": field.support_fingerprint,
        "pair_pool_hash": _pair_pool_hash(field),
        "component_zero_base_coefficient_hash": _component_zero_coefficient_hash(field),
        "field_fingerprint": field.fingerprint,
        "order_energy_relative_errors": json.dumps(energy_errors, sort_keys=True),
        "maximum_operator_energy_relative_error": float(max(energy_errors.values())),
        "minimum_sampled_operator_gram_rank": int(min(effective_ranks)),
        "minimum_sampled_operator_gram_relative_eigenvalue": float(
            min(minimum_relative_eigenvalues)
        ),
        "median_sampled_operator_gram_participation_ratio": float(
            np.median(participation_ratios)
        ),
        "sampled_operator_gram_relative_eigenvalue_quantiles": json.dumps(
            {
                "minimum": float(np.min(minimum_relative_eigenvalues)),
                "q05": float(np.quantile(minimum_relative_eigenvalues, 0.05)),
                "median": float(np.median(minimum_relative_eigenvalues)),
                "maximum": float(np.max(minimum_relative_eigenvalues)),
            },
            sort_keys=True,
        ),
        "expected_field_variance": float(field.expected_field_variance()),
        "maximum_sampled_contrast_mean_abs": maximum_mean,
        "maximum_sampled_contrast_second_moment_error": maximum_second_error,
        "sampled_nonself_tie_fraction": tie_fraction,
        "sampled_binary_10_10_balance_pass": bool(binary_balanced),
    }


def _augment_texture_row(
    row: dict[str, Any],
    arrays: Mapping[str, np.ndarray],
    *,
    arm: str,
    core_label: str,
) -> dict[str, Any]:
    latent_epsilon = np.asarray(arrays["latent_epsilon"], dtype=np.float64)
    row["representation_arm"] = arm
    row["additive_core_mode"] = core_label
    row["latent_epsilon_exact_zero_fraction"] = float(
        np.mean(latent_epsilon == 0.0)
    )
    row["latent_epsilon_near_zero_fraction"] = float(
        np.mean(np.abs(latent_epsilon) <= 1.0e-12)
    )
    row["latent_epsilon_abs_median"] = float(
        np.median(np.abs(latent_epsilon))
    )
    row["latent_epsilon_abs_q95"] = float(
        np.quantile(np.abs(latent_epsilon), 0.95)
    )
    centered = latent_epsilon - float(np.mean(latent_epsilon))
    variance = float(np.mean(centered * centered))
    row["latent_epsilon_excess_kurtosis"] = (
        float(np.mean(centered ** 4) / (variance * variance) - 3.0)
        if variance > 0.0
        else float("nan")
    )
    return row


def generate_squares(
    field: RankedCategoricalGNKField,
    reference: np.ndarray,
    *,
    core_seed: int,
    gnk_seed: int,
    stratum: str,
    background: str,
    count: int,
) -> dict[str, Any]:
    """Generate paired squares, with a genuinely unsupported uniform control."""

    if stratum.startswith("candidate_supported_order_"):
        return exp062.generate_squares(
            field,
            reference,
            core_seed=core_seed,
            gnk_seed=gnk_seed,
            stratum=stratum,
            background=background,
            count=count,
        )
    if stratum != "uniform_unsupported_site_pair":
        raise ValueError(f"unsupported square stratum: {stratum}")

    supported = {
        tuple(pair)
        for order in field.orders
        for pair in field.support_pairs(order)
    }
    pair_pool = tuple(
        (left, right)
        for left in range(field.length)
        for right in range(left + 1, field.length)
        if (left, right) not in supported
    )
    if not pair_pool:
        raise ValueError("no unsupported site pairs remain for the control")
    rng = np.random.default_rng(stable_uint64(
        "EXP-063-squares",
        int(core_seed),
        int(gnk_seed),
        int(reference.size),
        str(stratum),
        str(background),
        int(count),
    ))
    if background == "reference":
        base = np.repeat(reference[None, :], int(count), axis=0)
    elif background == "uniform_q20":
        base = rng.integers(
            0,
            field.alphabet_size,
            size=(int(count), field.length),
            dtype=np.int64,
        )
    else:
        raise ValueError(f"unsupported square background: {background}")
    first = base.copy()
    second = base.copy()
    both = base.copy()
    sites = np.empty((int(count), 2), dtype=np.int64)
    for row in range(int(count)):
        left, right = pair_pool[int(rng.integers(0, len(pair_pool)))]
        sites[row] = (left, right)
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
        "candidate_support_union_hash": latent.json_hash(
            [list(pair) for pair in sorted(supported)]
        ),
    })
    return {
        "base": base,
        "first": first,
        "second": second,
        "both": both,
        "sites": sites,
        "square_hash": square_hash,
    }


def _pooled_texture_row(
    *,
    pieces: Sequence[Mapping[str, np.ndarray]],
    hashes: Sequence[str],
    ids: Mapping[str, Any],
    arm: str,
    core_label: str,
    amplitude: float,
    sigma_g: float,
    frozen: Mapping[str, Any],
    helper: Mapping[str, Any],
) -> dict[str, Any]:
    arrays = {
        key: np.concatenate([np.asarray(piece[key]) for piece in pieces])
        for key in ("x", "tau", "tau_stable", "epsilon_r", "latent_epsilon")
    }
    metric_core = (
        "actual_additive_core"
        if core_label == "actual_additive_core"
        else "frozen_additive_core"
    )
    row = exp062.texture_metric_row(
        ids=ids,
        scale_mode=exp062.PRIMARY_GLOBAL,
        core_mode=metric_core,
        field_arm="extensive_balanced_bounded",
        amplitude=float(amplitude),
        stratum="pooled_candidate_supported_orders",
        background="reference",
        square_hash=latent.json_hash(sorted(str(value) for value in hashes)),
        arrays=arrays,
        scale_sigma_a=float(frozen["global_sigma_a"]),
        sigma_g=float(sigma_g),
        beta=float(frozen["beta"]),
        config=helper,
        pooled=True,
        permutation_blocks=[
            (np.asarray(piece["x"]), np.asarray(piece["tau"])) for piece in pieces
        ],
    )
    return _augment_texture_row(
        row, arrays, arm=arm, core_label=core_label
    )


def context_task(payload: Mapping[str, Any]) -> dict[str, Any]:
    config = payload["config"]
    helper = payload["helper_config"]
    frozen = config["frozen_operator"]
    reference_row = dict(payload["reference"])
    core_seed = int(payload["core_seed"])
    gnk_seed = int(payload["gnk_seed"])
    signature = str(payload["task_signature"])
    checkpoint = Path(str(payload["checkpoint_path"]))
    if checkpoint.exists():
        cached = json.loads(checkpoint.read_text(encoding="utf-8"))
        if str(cached.get("task_signature")) == signature:
            return cached

    reference = latent.sequence_to_indices(str(reference_row["sequence"]))
    core = exp062.build_core(reference, core_seed, helper)
    matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    reference_score = float(latent.additive_scores(matrix, reference[None, :])[0])
    task_sigma_a = exp062.complete_h1_rms(matrix, reference)
    sigma_g = math.sqrt(additive_mutation_energy(matrix))
    field_seed = int(gnk_seed) + int(config["candidate"]["seed_offset"])
    output_scale = float(frozen["output_scale"])
    ids = {
        "reference_id": reference_row["reference_id"],
        "length_bin": reference_row["length_bin"],
        "length": int(reference.size),
        "core_seed": core_seed,
        "gnk_seed": gnk_seed,
        "task_sigma_a": task_sigma_a,
    }

    h1_sources, h1_targets = exp062.complete_h1(reference)
    h1_hash = latent.json_hash({
        "sources": latent.array_hash(h1_sources),
        "targets": latent.array_hash(h1_targets),
        "label": "complete_h1",
    })
    h1_a_source = latent.additive_scores(matrix, h1_sources)
    h1_a_target = latent.additive_scores(matrix, h1_targets)

    baseline = build_field(
        matrix=matrix,
        reference=reference,
        field_seed=field_seed,
        arm="binary_r1",
        config=config,
    )
    square_panels: list[dict[str, Any]] = []
    square_hashes: list[str] = []
    for stratum in config["probe"]["square_strata"]:
        for background in config["probe"]["square_backgrounds"]:
            square = generate_squares(
                baseline,
                reference,
                core_seed=core_seed,
                gnk_seed=gnk_seed,
                stratum=str(stratum),
                background=str(background),
                count=int(config["probe"]["squares_per_stratum"]),
            )
            nodes = [
                square["base"],
                square["first"],
                square["second"],
                square["both"],
            ]
            square_panels.append({
                "stratum": str(stratum),
                "background": str(background),
                "hash": str(square["square_hash"]),
                "nodes": nodes,
                "additive_nodes": [
                    latent.additive_scores(matrix, values) for values in nodes
                ],
            })
            square_hashes.append(str(square["square_hash"]))
    common_square_hash = latent.json_hash(sorted(square_hashes))

    edge_rows: list[dict[str, Any]] = []
    texture_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    technical_rows: list[dict[str, Any]] = []
    amplitudes = [
        float(value) for value in config["intervention"]["dimensionless_amplitudes"]
    ]

    for arm in ARM_COMPLEXITY:
        field = (
            baseline
            if arm == "binary_r1"
            else build_field(
                matrix=matrix,
                reference=reference,
                field_seed=field_seed,
                arm=arm,
                config=config,
            )
        )
        audit = representation_audit(
            field,
            arm=arm,
            ids=ids,
            rank_relative_threshold=float(
                config["acceptance"]["technical"][
                    "rank_relative_singular_value_min"
                ]
            ),
        )
        audit_rows.append(audit)

        g_source, g_target = exp062.score_nodes(
            field.score_batch,
            [h1_sources, h1_targets],
            chunk_size=int(config["probe"]["eval_chunk_size"]),
        )
        maximum_formula_parity = 0.0
        raw = exp062.raw_parameterization(
            additive_reference_score=reference_score,
            sigma_a=float(frozen["global_sigma_a"]),
            sigma_g=sigma_g,
            beta=float(frozen["beta"]),
            mu=float(frozen["mu"]),
            amplitude=1.0,
        )
        probe_count = min(64, len(h1_a_target))
        dimensionless = exp062.dimensionless_fitness(
            h1_a_target[:probe_count],
            g_target[:probe_count],
            additive_reference_score=reference_score,
            sigma_a=float(frozen["global_sigma_a"]),
            sigma_g=sigma_g,
            beta=float(frozen["beta"]),
            mu=float(frozen["mu"]),
            amplitude=1.0,
            output_scale=output_scale,
        )
        raw_values = output_scale * exp062.sigmoid(
            raw["steepness"]
            * (
                h1_a_target[:probe_count]
                + raw["gnk_weight"] * g_target[:probe_count]
                - raw["threshold"]
            )
        )
        maximum_formula_parity = float(
            np.max(np.abs(dimensionless - raw_values))
        )
        for amplitude in amplitudes:
            row = exp062.edge_metric_row(
                ids=ids,
                sample_id="reference_shell_h1_complete",
                pairing_hash=h1_hash,
                hamming_distance=1,
                a_source=h1_a_source,
                a_target=h1_a_target,
                g_source=g_source,
                g_target=g_target,
                mode=exp062.PRIMARY_GLOBAL,
                amplitude=amplitude,
                scale_sigma_a=float(frozen["global_sigma_a"]),
                sigma_g=sigma_g,
                beta=float(frozen["beta"]),
                mu=float(frozen["mu"]),
                reference_score=reference_score,
                public_threshold=float(core.stability_threshold),
                output_scale=output_scale,
                config=helper,
            )
            row["representation_arm"] = arm
            edge_rows.append(row)

        pooled_actual: dict[float, list[dict[str, np.ndarray]]] = {
            amplitude: [] for amplitude in amplitudes
        }
        pooled_actual_hashes: dict[float, list[str]] = {
            amplitude: [] for amplitude in amplitudes
        }
        pooled_midpoint: list[dict[str, np.ndarray]] = []
        pooled_midpoint_hashes: list[str] = []
        pooled_saturation: list[dict[str, np.ndarray]] = []
        pooled_saturation_hashes: list[str] = []
        for panel in square_panels:
            gnk_nodes = exp062.score_nodes(
                field.score_batch,
                panel["nodes"],
                chunk_size=int(config["probe"]["eval_chunk_size"]),
            )
            for amplitude in amplitudes:
                arrays = exp062.texture_arrays(
                    additive_nodes=panel["additive_nodes"],
                    gnk_nodes=gnk_nodes,
                    core_mode="actual_additive_core",
                    reference_score=reference_score,
                    scale_sigma_a=float(frozen["global_sigma_a"]),
                    sigma_g=sigma_g,
                    beta=float(frozen["beta"]),
                    mu=float(frozen["mu"]),
                    amplitude=amplitude,
                    output_scale=output_scale,
                )
                metric = exp062.texture_metric_row(
                    ids=ids,
                    scale_mode=exp062.PRIMARY_GLOBAL,
                    core_mode="actual_additive_core",
                    field_arm="extensive_balanced_bounded",
                    amplitude=amplitude,
                    stratum=panel["stratum"],
                    background=panel["background"],
                    square_hash=panel["hash"],
                    arrays=arrays,
                    scale_sigma_a=float(frozen["global_sigma_a"]),
                    sigma_g=sigma_g,
                    beta=float(frozen["beta"]),
                    config=helper,
                    pooled=False,
                    permutation_blocks=None,
                )
                texture_rows.append(
                    _augment_texture_row(
                        metric,
                        arrays,
                        arm=arm,
                        core_label="actual_additive_core",
                    )
                )
                if (
                    panel["background"] == "reference"
                    and panel["stratum"].startswith("candidate_supported_order_")
                ):
                    pooled_actual[amplitude].append(arrays)
                    pooled_actual_hashes[amplitude].append(panel["hash"])

            midpoint_score = (
                reference_score
                - float(frozen["mu"])
                * float(frozen["global_sigma_a"])
                / float(frozen["beta"])
            )
            midpoint_nodes = [
                np.full_like(values, midpoint_score, dtype=np.float64)
                for values in panel["additive_nodes"]
            ]
            midpoint_arrays = exp062.texture_arrays(
                additive_nodes=midpoint_nodes,
                gnk_nodes=gnk_nodes,
                core_mode="frozen_additive_core",
                reference_score=reference_score,
                scale_sigma_a=float(frozen["global_sigma_a"]),
                sigma_g=sigma_g,
                beta=float(frozen["beta"]),
                mu=float(frozen["mu"]),
                amplitude=1.0,
                output_scale=output_scale,
            )
            midpoint_metric = exp062.texture_metric_row(
                ids=ids,
                scale_mode=exp062.PRIMARY_GLOBAL,
                core_mode="frozen_additive_core",
                field_arm="extensive_balanced_bounded",
                amplitude=1.0,
                stratum=panel["stratum"],
                background=panel["background"],
                square_hash=panel["hash"],
                arrays=midpoint_arrays,
                scale_sigma_a=float(frozen["global_sigma_a"]),
                sigma_g=sigma_g,
                beta=float(frozen["beta"]),
                config=helper,
                pooled=False,
                permutation_blocks=None,
            )
            texture_rows.append(
                _augment_texture_row(
                    midpoint_metric,
                    midpoint_arrays,
                    arm=arm,
                    core_label="frozen_midpoint_occupancy_0.5",
                )
            )

            saturation_nodes = [
                np.full_like(values, reference_score, dtype=np.float64)
                for values in panel["additive_nodes"]
            ]
            saturation_arrays = exp062.texture_arrays(
                additive_nodes=saturation_nodes,
                gnk_nodes=gnk_nodes,
                core_mode="frozen_additive_core",
                reference_score=reference_score,
                scale_sigma_a=float(frozen["global_sigma_a"]),
                sigma_g=sigma_g,
                beta=float(frozen["beta"]),
                mu=float(frozen["mu"]),
                amplitude=1.0,
                output_scale=output_scale,
            )
            saturation_metric = exp062.texture_metric_row(
                ids=ids,
                scale_mode=exp062.PRIMARY_GLOBAL,
                core_mode="frozen_additive_core",
                field_arm="extensive_balanced_bounded",
                amplitude=1.0,
                stratum=panel["stratum"],
                background=panel["background"],
                square_hash=panel["hash"],
                arrays=saturation_arrays,
                scale_sigma_a=float(frozen["global_sigma_a"]),
                sigma_g=sigma_g,
                beta=float(frozen["beta"]),
                config=helper,
                pooled=False,
                permutation_blocks=None,
            )
            texture_rows.append(
                _augment_texture_row(
                    saturation_metric,
                    saturation_arrays,
                    arm=arm,
                    core_label="frozen_saturation_occupancy_0.99",
                )
            )
            if (
                panel["background"] == "reference"
                and panel["stratum"].startswith("candidate_supported_order_")
            ):
                pooled_midpoint.append(midpoint_arrays)
                pooled_midpoint_hashes.append(panel["hash"])
                pooled_saturation.append(saturation_arrays)
                pooled_saturation_hashes.append(panel["hash"])

        for amplitude in amplitudes:
            texture_rows.append(_pooled_texture_row(
                pieces=pooled_actual[amplitude],
                hashes=pooled_actual_hashes[amplitude],
                ids=ids,
                arm=arm,
                core_label="actual_additive_core",
                amplitude=amplitude,
                sigma_g=sigma_g,
                frozen=frozen,
                helper=helper,
            ))
        texture_rows.append(_pooled_texture_row(
            pieces=pooled_midpoint,
            hashes=pooled_midpoint_hashes,
            ids=ids,
            arm=arm,
            core_label="frozen_midpoint_occupancy_0.5",
            amplitude=1.0,
            sigma_g=sigma_g,
            frozen=frozen,
            helper=helper,
        ))
        texture_rows.append(_pooled_texture_row(
            pieces=pooled_saturation,
            hashes=pooled_saturation_hashes,
            ids=ids,
            arm=arm,
            core_label="frozen_saturation_occupancy_0.99",
            amplitude=1.0,
            sigma_g=sigma_g,
            frozen=frozen,
            helper=helper,
        ))

        legacy_parity = 0.0
        if arm == "binary_r1":
            legacy = ExtensiveGNKField(
                additive_matrix=matrix,
                reference_sequence=reference,
                seed=field_seed,
                orders=tuple(int(value) for value in config["candidate"]["orders"]),
                order_weights=tuple(
                    float(value) for value in config["candidate"]["order_weights"]
                ),
                target_degrees=tuple(
                    float(value) for value in config["candidate"]["target_degrees"]
                ),
            )
            legacy_values = legacy.score_batch(h1_targets)
            ranked_values = field.score_batch(h1_targets)
            legacy_parity = float(np.max(np.abs(legacy_values - ranked_values)))

        arm_edges = [
            row
            for row in edge_rows
            if row["representation_arm"] == arm
        ]
        arm_textures = [
            row
            for row in texture_rows
            if row["representation_arm"] == arm
        ]
        actual_primary = next(
            row
            for row in arm_textures
            if row["additive_core_mode"] == "actual_additive_core"
            and bool(row["pooled_supported_reference"])
            and math.isclose(float(row["amplitude"]), 1.0)
        )
        frozen_saturation = next(
            row
            for row in arm_textures
            if row["additive_core_mode"] == "frozen_saturation_occupancy_0.99"
            and bool(row["pooled_supported_reference"])
            and math.isclose(float(row["amplitude"]), 1.0)
        )
        saturation_rms = float(frozen_saturation["direct_tau_y_rms"])
        actual_to_saturation_ratio = (
            float(actual_primary["direct_tau_y_rms"]) / saturation_rms
            if saturation_rms > 0.0
            else float("inf")
        )
        zero_rows = [
            row
            for row in arm_edges
            if float(row["amplitude"]) == 0.0
        ]
        technical_rows.append({
            **ids,
            "representation_arm": arm,
            "support_fingerprint": audit["support_fingerprint"],
            "pair_pool_hash": audit["pair_pool_hash"],
            "component_zero_base_coefficient_hash": audit[
                "component_zero_base_coefficient_hash"
            ],
            "pair_design_hash": h1_hash,
            "square_design_hash": common_square_hash,
            "binary_rank1_legacy_max_abs_difference": legacy_parity,
            "maximum_formula_parity_abs_difference": maximum_formula_parity,
            "complete_h1_gnk_mutation_rms": float(
                np.sqrt(np.mean((g_target - g_source) ** 2))
            ),
            "complete_h1_gnk_mutation_rms_ratio": float(
                np.sqrt(np.mean((g_target - g_source) ** 2)) / sigma_g
            ),
            "maximum_zero_amplitude_residual_rms": (
                max(float(row["direct_y_rms"]) for row in zero_rows)
                if zero_rows
                else float("nan")
            ),
            "candidate_reference_center_abs": abs(
                float(field.score_batch(reference[None, :])[0])
            ),
            "maximum_operator_energy_relative_error": audit[
                "maximum_operator_energy_relative_error"
            ],
            "minimum_sampled_operator_gram_rank": audit[
                "minimum_sampled_operator_gram_rank"
            ],
            "minimum_sampled_operator_gram_relative_eigenvalue": audit[
                "minimum_sampled_operator_gram_relative_eigenvalue"
            ],
            "median_sampled_operator_gram_participation_ratio": audit[
                "median_sampled_operator_gram_participation_ratio"
            ],
            "actual_to_frozen_saturation_tau_rms_ratio": (
                actual_to_saturation_ratio
            ),
            "sampled_nonself_tie_fraction": audit[
                "sampled_nonself_tie_fraction"
            ],
            "sampled_binary_10_10_balance_pass": audit[
                "sampled_binary_10_10_balance_pass"
            ],
            "maximum_sampled_contrast_mean_abs": audit[
                "maximum_sampled_contrast_mean_abs"
            ],
            "maximum_sampled_contrast_second_moment_error": audit[
                "maximum_sampled_contrast_second_moment_error"
            ],
            "all_outputs_finite": bool(
                all(bool(row["all_finite"]) for row in arm_edges)
                and all(bool(row["all_finite"]) for row in arm_textures)
            ),
        })
        if arm != "binary_r1":
            del field

    surrogate = exp062.additive_surrogate_matrix(
        matrix,
        seed=gnk_seed + int(config["candidate"]["surrogate_seed_offset"]),
    )
    surrogate_sigma = math.sqrt(additive_mutation_energy(surrogate))
    surrogate_pieces: list[dict[str, np.ndarray]] = []
    surrogate_hashes: list[str] = []
    for panel in square_panels:
        if not (
            panel["background"] == "reference"
            and panel["stratum"].startswith("candidate_supported_order_")
        ):
            continue
        surrogate_nodes = exp062.score_nodes(
            lambda values: exp062.score_additive_field(
                surrogate, reference, values
            ),
            panel["nodes"],
            chunk_size=int(config["probe"]["eval_chunk_size"]),
        )
        arrays = exp062.texture_arrays(
            additive_nodes=panel["additive_nodes"],
            gnk_nodes=surrogate_nodes,
            core_mode="actual_additive_core",
            reference_score=reference_score,
            scale_sigma_a=float(frozen["global_sigma_a"]),
            sigma_g=surrogate_sigma,
            beta=float(frozen["beta"]),
            mu=float(frozen["mu"]),
            amplitude=1.0,
            output_scale=output_scale,
        )
        surrogate_pieces.append(arrays)
        surrogate_hashes.append(panel["hash"])
    surrogate_row = _pooled_texture_row(
        pieces=surrogate_pieces,
        hashes=surrogate_hashes,
        ids=ids,
        arm="shared_additive_surrogate",
        core_label="actual_additive_core",
        amplitude=1.0,
        sigma_g=surrogate_sigma,
        frozen=frozen,
        helper=helper,
    )
    surrogate_max = float(surrogate_row["latent_epsilon_max_abs"])
    texture_rows.append(surrogate_row)
    for row in technical_rows:
        row["surrogate_latent_epsilon_max_abs"] = surrogate_max

    result = {
        "task_signature": signature,
        "technical_rows": technical_rows,
        "audit_rows": audit_rows,
        "edge_rows": edge_rows,
        "texture_rows": texture_rows,
    }
    strict_json_atomic(checkpoint, result)
    return result


def pairing_audit(
    technical: pd.DataFrame,
    config: Mapping[str, Any],
) -> tuple[bool, list[dict[str, Any]]]:
    keys = ["reference_id", "length_bin", "length", "core_seed", "gnk_seed"]
    rows: list[dict[str, Any]] = []
    tolerance = config["acceptance"]["technical"]
    for values, group in technical.groupby(keys, sort=True):
        arms = sorted(str(value) for value in group["representation_arm"])
        rank_checks = []
        for _, row in group.iterrows():
            expected_rank = int(
                arm_spec(config, str(row["representation_arm"]))[1]
            )
            rank_checks.append(
                int(row["minimum_sampled_operator_gram_rank"])
                >= (1 if expected_rank == 1 else int(config["acceptance"]["manipulation"]["rank8_effective_rank_min"]))
            )
        dense = group[group["representation_arm"].str.startswith("dense_q20")]
        passed = bool(
            arms == sorted(ARM_COMPLEXITY)
            and int(group["support_fingerprint"].nunique()) == 1
            and int(group["pair_pool_hash"].nunique()) == 1
            and int(group["component_zero_base_coefficient_hash"].nunique()) == 1
            and int(group["pair_design_hash"].nunique()) == 1
            and int(group["square_design_hash"].nunique()) == 1
            and bool(group["all_outputs_finite"].all())
            and float(group["maximum_operator_energy_relative_error"].max())
            <= float(tolerance["operator_energy_relative_error_max"])
            and float(group["maximum_formula_parity_abs_difference"].max())
            <= float(tolerance["formula_parity_atol"])
            and float(group["candidate_reference_center_abs"].max())
            <= float(tolerance["reference_center_atol"])
            and float(group["maximum_sampled_contrast_mean_abs"].max())
            <= float(tolerance["contrast_mean_atol"])
            and float(
                group["maximum_sampled_contrast_second_moment_error"].max()
            )
            <= float(tolerance["contrast_second_moment_atol"])
            and float(group["surrogate_latent_epsilon_max_abs"].max())
            <= float(tolerance["additive_surrogate_epistasis_atol"])
            and float(group["maximum_zero_amplitude_residual_rms"].max())
            <= float(tolerance["zero_amplitude_atol"])
            and float(
                group.loc[
                    group["representation_arm"] == "binary_r1",
                    "binary_rank1_legacy_max_abs_difference",
                ].max()
            )
            <= float(tolerance["binary_rank1_score_parity_atol"])
            and (
                dense.empty
                or float(dense["sampled_nonself_tie_fraction"].max())
                <= float(
                    config["acceptance"]["manipulation"][
                        "dense_q20_nonself_tie_fraction_max"
                    ]
                )
            )
            and all(rank_checks)
            and bool(
                group.loc[
                    group["representation_arm"].str.startswith("binary"),
                    "sampled_binary_10_10_balance_pass",
                ].all()
            )
        )
        rows.append({
            **dict(zip(keys, values)),
            "row_type": "cross_arm_pairing",
            "representation_arm": "all_arms",
            "arms": json.dumps(arms),
            "support_hash_count": int(group["support_fingerprint"].nunique()),
            "pair_pool_hash_count": int(group["pair_pool_hash"].nunique()),
            "component_zero_coefficient_hash_count": int(
                group["component_zero_base_coefficient_hash"].nunique()
            ),
            "pair_design_hash_count": int(group["pair_design_hash"].nunique()),
            "square_design_hash_count": int(group["square_design_hash"].nunique()),
            "pairing_and_manipulation_pass": passed,
        })
    return bool(rows and all(row["pairing_and_manipulation_pass"] for row in rows)), rows


def _minimum_bin_fraction(frame: pd.DataFrame, column: str) -> tuple[float, dict[str, float]]:
    if frame.empty:
        return 0.0, {}
    fractions = frame.groupby("length_bin")[column].mean().sort_index()
    return float(fractions.min()), {
        str(key): float(value) for key, value in fractions.items()
    }


def arm_gate(
    *,
    arm: str,
    edge: pd.DataFrame,
    texture: pd.DataFrame,
    helper: Mapping[str, Any],
    config: Mapping[str, Any],
    technical_pass: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    transmission = config["acceptance"]["transmission"]
    edge_selected = edge[
        (edge["representation_arm"] == arm)
        & np.isclose(edge["amplitude"], float(transmission["primary_amplitude"]))
        & (edge["sample_id"] == "reference_shell_h1_complete")
    ].copy()
    edge_selected["cell_pass"] = (
        (edge_selected["direct_practical_resolution_fraction"] >= float(transmission["edge_practical_fraction_min"]))
        & (edge_selected["direct_correlation"] >= float(transmission["correlation_min"]))
    )
    edge_min, edge_bins = _minimum_bin_fraction(edge_selected, "cell_pass")
    edge_pass = bool(
        edge_bins
        and edge_min
        >= float(transmission["minimum_cell_pass_fraction_per_length_bin"])
    )

    texture_flags: dict[str, Any] = {}
    checks: list[dict[str, Any]] = []
    for core_label, flag_name in (
        ("actual_additive_core", "actual_texture"),
        ("frozen_midpoint_occupancy_0.5", "midpoint_texture"),
    ):
        selected = texture[
            (texture["representation_arm"] == arm)
            & (texture["additive_core_mode"] == core_label)
            & texture["pooled_supported_reference"].eq(True)
            & np.isclose(texture["amplitude"], float(transmission["primary_amplitude"]))
        ].copy()
        selected["cell_pass"] = (
            (selected["direct_tau_practical_resolution_fraction"] >= float(transmission["texture_practical_fraction_min"]))
            & (selected["direct_tau_correlation"] >= float(transmission["correlation_min"]))
            & (selected["direct_tau_projection"] > 0.0)
            & (selected["one_sided_pvalue"] <= float(config["acceptance"]["texture_permutation_alpha"]))
        )
        minimum, bins = _minimum_bin_fraction(selected, "cell_pass")
        projection_lower = exp062.hierarchical_bootstrap_lower(
            selected,
            metric="observed_projection",
            replicates=int(config["acceptance"]["bootstrap_replicates"]),
            seed=int(config["acceptance"]["bootstrap_seed"])
            + (10 if flag_name == "actual_texture" else 20)
            + 100 * ARM_COMPLEXITY.index(arm),
        )
        passed = bool(
            bins
            and minimum
            >= float(transmission["minimum_cell_pass_fraction_per_length_bin"])
            and projection_lower
            > float(transmission["projection_bootstrap_lower_min"])
        )
        texture_flags[flag_name] = passed
        texture_flags[f"{flag_name}_values"] = {
            "n_cells": int(len(selected)),
            "minimum_length_bin_pass_fraction": minimum,
            "per_length_bin_pass_fraction": bins,
            "projection_bootstrap_lower": projection_lower,
        }

    arm_texture = texture[
        (texture["representation_arm"] == arm)
        & (texture["additive_core_mode"] == "actual_additive_core")
    ].copy()
    direction_pass, dose_rows, direction_values = exp062.dose_response_gate(
        arm_texture,
        mode=exp062.PRIMARY_GLOBAL,
        config=helper,
    )
    length_pass, length_values = exp062.length_stability_gate(
        arm_texture,
        mode=exp062.PRIMARY_GLOBAL,
        config=helper,
    )
    for row in dose_rows:
        row["representation_arm"] = arm
    length_values = {**length_values, "representation_arm": arm}

    checks.extend([
        {
            "check_id": f"{arm}_technical_and_pairing",
            "representation_arm": arm,
            "family": "technical",
            "required": False,
            "decision_relevant": True,
            "passed": bool(technical_pass),
            "criterion": "paired supports/probes, exact energy, replay, and null controls pass",
            "value": json.dumps({"technical_pass": bool(technical_pass)}),
        },
        {
            "check_id": f"{arm}_complete_h1_edge_transfer",
            "representation_arm": arm,
            "family": "edge",
            "required": False,
            "decision_relevant": True,
            "passed": edge_pass,
            "criterion": "at least 80% of cells per length bin pass complete-h1 edge gates",
            "value": json.dumps({
                "minimum_length_bin_pass_fraction": edge_min,
                "per_length_bin_pass_fraction": edge_bins,
            }, sort_keys=True),
        },
        {
            "check_id": f"{arm}_actual_core_texture",
            "representation_arm": arm,
            "family": "texture",
            "required": False,
            "decision_relevant": True,
            "passed": bool(texture_flags["actual_texture"]),
            "criterion": "actual-core supported texture passes prevalence, alignment, and permutation gates",
            "value": json.dumps(texture_flags["actual_texture_values"], sort_keys=True),
        },
        {
            "check_id": f"{arm}_midpoint_frozen_positive_control",
            "representation_arm": arm,
            "family": "control",
            "required": False,
            "decision_relevant": True,
            "passed": bool(texture_flags["midpoint_texture"]),
            "criterion": "occupancy-0.5 frozen texture positive control passes",
            "value": json.dumps(texture_flags["midpoint_texture_values"], sort_keys=True),
        },
        {
            "check_id": f"{arm}_dose_directionality",
            "representation_arm": arm,
            "family": "directionality",
            "required": False,
            "decision_relevant": True,
            "passed": bool(direction_pass),
            "criterion": "actual-core texture is directional over a=0.25..2",
            "value": json.dumps(_clean_json(direction_values), sort_keys=True),
        },
        {
            "check_id": f"{arm}_length_stability",
            "representation_arm": arm,
            "family": "length_stability",
            "required": False,
            "decision_relevant": True,
            "passed": bool(length_pass),
            "criterion": "a=1 protein-median texture RMS passes length-slope bounds",
            "value": json.dumps(_clean_json(length_values), sort_keys=True),
        },
    ])
    flags = {
        "technical_pass": bool(technical_pass),
        "edge_pass": edge_pass,
        "actual_texture_pass": bool(texture_flags["actual_texture"]),
        "midpoint_control_pass": bool(texture_flags["midpoint_texture"]),
        "directionality_pass": bool(direction_pass),
        "length_stability_pass": bool(length_pass),
    }
    flags["package_pass_before_paired_superiority"] = bool(all(flags.values()))
    return checks, flags, dose_rows, length_values


def paired_effects(
    texture: pd.DataFrame,
    config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, bool]]:
    selected = texture[
        texture["representation_arm"].isin(ARM_COMPLEXITY)
        & (texture["additive_core_mode"] == "actual_additive_core")
        & texture["pooled_supported_reference"].eq(True)
        & np.isclose(texture["amplitude"], 1.0)
    ][
        [
            "reference_id",
            "length_bin",
            "length",
            "core_seed",
            "gnk_seed",
            "representation_arm",
            "direct_tau_practical_resolution_fraction",
        ]
    ].copy()
    pivot = selected.pivot(
        index=["reference_id", "length_bin", "length", "core_seed", "gnk_seed"],
        columns="representation_arm",
        values="direct_tau_practical_resolution_fraction",
    ).reset_index()
    if any(arm not in pivot.columns for arm in ARM_COMPLEXITY):
        return [], [], {arm: False for arm in RESCUE_ARMS}

    effect_rows: list[dict[str, Any]] = []
    superiority: dict[str, bool] = {}
    gain_min = float(
        config["acceptance"]["manipulation"]["median_practical_fraction_gain_min"]
    )
    lower_min = float(
        config["acceptance"]["manipulation"][
            "candidate_minus_incumbent_practical_fraction_bootstrap_lower_min"
        ]
    )
    for arm in RESCUE_ARMS:
        frame = pivot[
            ["reference_id", "length_bin", "length", "core_seed", "gnk_seed"]
        ].copy()
        frame["difference"] = pivot[arm] - pivot["binary_r1"]
        lower = exp062.hierarchical_bootstrap_lower(
            frame,
            metric="difference",
            replicates=int(config["acceptance"]["bootstrap_replicates"]),
            seed=int(config["acceptance"]["bootstrap_seed"])
            + 600
            + ARM_COMPLEXITY.index(arm),
        )
        protein = frame.groupby(
            ["reference_id", "length_bin", "length"], as_index=False
        )["difference"].median()
        median_gain = float(np.median(protein["difference"]))
        passed = bool(median_gain >= gain_min and lower > lower_min)
        superiority[arm] = passed
        for _, row in protein.iterrows():
            effect_rows.append({
                "row_type": "protein_median_paired_difference",
                "comparison_arm": arm,
                "reference_id": row["reference_id"],
                "length_bin": row["length_bin"],
                "length": int(row["length"]),
                "practical_fraction_difference_vs_binary_r1": float(row["difference"]),
                "median_gain": None,
                "bootstrap_lower": None,
                "paired_superiority_pass": None,
            })
        effect_rows.append({
            "row_type": "comparison_summary",
            "comparison_arm": arm,
            "reference_id": None,
            "length_bin": None,
            "length": None,
            "practical_fraction_difference_vs_binary_r1": None,
            "median_gain": median_gain,
            "bootstrap_lower": lower,
            "paired_superiority_pass": passed,
        })

    factorial = pd.DataFrame({
        "reference_id": pivot["reference_id"],
        "length_bin": pivot["length_bin"],
        "length": pivot["length"],
        "core_seed": pivot["core_seed"],
        "gnk_seed": pivot["gnk_seed"],
        "dense_main_effect": 0.5
        * (
            (pivot["dense_q20_r1"] - pivot["binary_r1"])
            + (pivot["dense_q20_r8"] - pivot["binary_r8"])
        ),
        "rank_main_effect": 0.5
        * (
            (pivot["binary_r8"] - pivot["binary_r1"])
            + (pivot["dense_q20_r8"] - pivot["dense_q20_r1"])
        ),
        "dense_by_rank_interaction": (
            (pivot["dense_q20_r8"] - pivot["binary_r8"])
            - (pivot["dense_q20_r1"] - pivot["binary_r1"])
        ),
    })
    factorial_rows: list[dict[str, Any]] = []
    for metric in (
        "dense_main_effect",
        "rank_main_effect",
        "dense_by_rank_interaction",
    ):
        protein = factorial.groupby(
            ["reference_id", "length_bin", "length"], as_index=False
        )[metric].median()
        for _, row in protein.iterrows():
            factorial_rows.append({
                "row_type": "protein_median_factorial_effect",
                "effect": metric,
                "reference_id": row["reference_id"],
                "length_bin": row["length_bin"],
                "length": int(row["length"]),
                "value": float(row[metric]),
                "median": None,
                "bootstrap_lower": None,
                "bootstrap_upper": None,
            })
        values = protein[metric].to_numpy(dtype=float)
        rng = np.random.default_rng(
            int(config["acceptance"]["bootstrap_seed"])
            + 700
            + list(
                (
                    "dense_main_effect",
                    "rank_main_effect",
                    "dense_by_rank_interaction",
                )
            ).index(metric)
        )
        boot = np.asarray([
            np.median(rng.choice(values, size=len(values), replace=True))
            for _ in range(int(config["acceptance"]["bootstrap_replicates"]))
        ])
        factorial_rows.append({
            "row_type": "factorial_effect_summary",
            "effect": metric,
            "reference_id": None,
            "length_bin": None,
            "length": None,
            "value": None,
            "median": float(np.median(values)),
            "bootstrap_lower": float(np.quantile(boot, 0.025)),
            "bootstrap_upper": float(np.quantile(boot, 0.975)),
        })
    return effect_rows, factorial_rows, superiority


def classify_decision(
    *,
    quick: bool,
    pairing_pass: bool,
    core_pass: bool,
    arm_flags: Mapping[str, Mapping[str, Any]],
    superiority: Mapping[str, bool],
) -> str:
    if quick:
        return "QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION"
    if not pairing_pass:
        return "EXECUTION_PAIRING_OR_MANIPULATION_INVALID"
    if not core_pass:
        return "FROZEN_EXP062_CORE_OPERATING_POINT_NOT_REPRODUCED"
    if bool(arm_flags["binary_r1"]["package_pass_before_paired_superiority"]):
        return "INCUMBENT_PREVALENCE_FAILURE_NOT_REPLICATED"
    eligible = {
        arm
        for arm in RESCUE_ARMS
        if bool(arm_flags[arm]["package_pass_before_paired_superiority"])
        and bool(superiority.get(arm, False))
    }
    if "dense_q20_r1" in eligible and "binary_r8" in eligible:
        return "ALTERNATIVE_SINGLE_FACTOR_REPAIRS_BOTH_SUFFICIENT"
    if "dense_q20_r1" in eligible:
        return "DENSE_Q20_RANK1_SUFFICIENT"
    if "binary_r8" in eligible:
        return "BINARY_RANK8_SUFFICIENT"
    if "dense_q20_r8" in eligible:
        return "COMBINED_DENSE_Q20_RANK8_REQUIRED"
    if any(
        bool(arm_flags[arm]["package_pass_before_paired_superiority"])
        for arm in RESCUE_ARMS
    ):
        return "APPARENT_RESCUE_NOT_PAIRED_SUPERIOR_TO_BINARY_R1"
    if not any(bool(arm_flags[arm]["midpoint_control_pass"]) for arm in RESCUE_ARMS):
        return "NO_INTERVENTION_PASSES_MIDPOINT_CONTROL"
    if not any(bool(arm_flags[arm]["edge_pass"]) for arm in RESCUE_ARMS):
        return "NO_INTERVENTION_PASSES_EDGE_TRANSFER"
    if not any(bool(arm_flags[arm]["actual_texture_pass"]) for arm in RESCUE_ARMS):
        return "NO_INTERVENTION_RESCUES_TEXTURE_PREVALENCE"
    if not any(bool(arm_flags[arm]["directionality_pass"]) for arm in RESCUE_ARMS):
        return "REPRESENTATION_RESCUE_NOT_DIRECTIONAL"
    if not any(bool(arm_flags[arm]["length_stability_pass"]) for arm in RESCUE_ARMS):
        return "REPRESENTATION_RESCUE_NOT_LENGTH_STABLE"
    return "REPRESENTATION_REPAIR_INSUFFICIENT_AT_FIXED_SIGMOID"


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, quick=bool(args.quick))
    validate_config(config)
    helper = helper_config(config)
    references, catalog_hash = load_references(project_root, config)
    implementation_hashes = source_hashes(project_root, config)
    implementation_hash = latent.json_hash(implementation_hashes)
    workers = max(
        1,
        int(os.environ.get(
            "LABPROJ_WORKERS",
            os.environ.get("SLURM_CPUS_PER_TASK", "1"),
        )),
    )
    core_seeds = [int(value) for value in config["replication"]["core_seeds"]]
    gnk_seeds = [int(value) for value in config["replication"]["gnk_seeds"]]
    preflight = {
        "experiment_id": "EXP-063",
        "quick": bool(args.quick),
        "scientific_decision_enabled": not bool(args.quick),
        "config_file_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "effective_config_hash": latent.json_hash(config),
        "reference_catalog_sha256": catalog_hash,
        "implementation_hash": implementation_hash,
        "implementation_hash_components": implementation_hashes,
        "reference_count": len(references),
        "length_bins": sorted({row["length_bin"] for row in references}),
        "core_seed_count": len(core_seeds),
        "gnk_seed_count": len(gnk_seeds),
        "paired_context_count": len(references) * len(core_seeds) * len(gnk_seeds),
        "representation_arms": list(ARM_COMPLEXITY),
        "frozen_operator": dict(config["frozen_operator"]),
        "panel_role": config["references"]["reuse_role"],
        "candidate_is_experiment_local": True,
        "public_model_changed": False,
        "empirical_outcomes_loaded": False,
        "optimizer_or_accessibility_tested": False,
        "stage_4c_4d_4e_hard_stop_retained": True,
    }
    strict_json_atomic(output_dir / "preflight.json", preflight)
    latent.write_csv_atomic(
        output_dir / "reference_catalog.csv",
        [
            {key: value for key, value in row.items() if key != "sequence"}
            for row in references
        ],
    )
    latent.write_csv_atomic(
        output_dir / "representation_catalog.csv",
        [
            {
                "representation_arm": str(row["id"]),
                "contrast_family": str(row["contrast_family"]),
                "representation_rank": int(row["representation_rank"]),
                "role": str(row["role"]),
                "orders": json.dumps(config["candidate"]["orders"]),
                "order_weights": json.dumps(config["candidate"]["order_weights"]),
                "target_degrees": json.dumps(config["candidate"]["target_degrees"]),
            }
            for row in config["candidate"]["arms"]
        ],
    )

    core_payloads = [
        {
            "reference": reference,
            "core_seed": core_seed,
            "helper_config": helper,
            "frozen_operator": config["frozen_operator"],
        }
        for reference in references
        for core_seed in core_seeds
    ]
    core_rows = parallel_execute(
        core_guard_task,
        core_payloads,
        workers=workers,
        label="core-guard",
    )
    core_rows.sort(
        key=lambda row: (
            str(row["length_bin"]),
            str(row["reference_id"]),
            int(row["core_seed"]),
        )
    )
    latent.write_csv_atomic(output_dir / "core_guardrail_metrics.csv", core_rows)
    shared_core_pass, core_values = core_guard_pass(core_rows, config)

    results: list[dict[str, Any]] = []
    if shared_core_pass or bool(args.quick):
        checkpoint_dir = output_dir / "checkpoints"
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        payloads: list[dict[str, Any]] = []
        for reference in references:
            for core_seed in core_seeds:
                for gnk_seed in gnk_seeds:
                    signature = latent.json_hash({
                        "reference": reference,
                        "core_seed": core_seed,
                        "gnk_seed": gnk_seed,
                        "effective_config_hash": preflight["effective_config_hash"],
                        "implementation_hash": implementation_hash,
                    })
                    payloads.append({
                        "reference": reference,
                        "core_seed": core_seed,
                        "gnk_seed": gnk_seed,
                        "config": config,
                        "helper_config": helper,
                        "task_signature": signature,
                        "checkpoint_path": str(
                            checkpoint_dir
                            / (
                                f"{latent.safe_name(reference['reference_id'])}"
                                f"__core-{core_seed}__gnk-{gnk_seed}"
                                f"__{signature[:16]}.json"
                            )
                        ),
                    })
        results = parallel_execute(
            context_task,
            payloads,
            workers=workers,
            label="paired-representation",
        )

    technical_rows = [
        row for result in results for row in result["technical_rows"]
    ]
    audit_rows = [row for result in results for row in result["audit_rows"]]
    edge_rows = [row for result in results for row in result["edge_rows"]]
    texture_rows = [row for result in results for row in result["texture_rows"]]
    technical = pd.DataFrame(technical_rows)
    edge = pd.DataFrame(edge_rows)
    texture = pd.DataFrame(texture_rows)
    if technical.empty:
        paired_pass, pairing_rows = False, []
    else:
        paired_pass, pairing_rows = pairing_audit(technical, config)
        paired_pass = bool(
            paired_pass
            and len(pairing_rows) == int(preflight["paired_context_count"])
        )
    for row in audit_rows:
        row["row_type"] = "arm_context_audit"
    for row in technical_rows:
        row["row_type"] = "technical_context_audit"
    def context_sort_key(row: Mapping[str, Any]) -> tuple[str, str, int, int, str]:
        return (
            str(row.get("length_bin", "")),
            str(row.get("reference_id", "")),
            int(row.get("core_seed", -1)),
            int(row.get("gnk_seed", -1)),
            str(row.get("representation_arm", "")),
        )
    audit_rows.sort(key=context_sort_key)
    technical_rows.sort(key=context_sort_key)
    pairing_rows.sort(key=context_sort_key)
    latent.write_csv_atomic(
        output_dir / "representation_audit.csv",
        audit_rows + technical_rows + pairing_rows,
    )
    edge_rows.sort(
        key=lambda row: (
            str(row["length_bin"]),
            str(row["reference_id"]),
            int(row["core_seed"]),
            int(row["gnk_seed"]),
            str(row["representation_arm"]),
            float(row["amplitude"]),
        )
    )
    texture_rows.sort(
        key=lambda row: (
            str(row["length_bin"]),
            str(row["reference_id"]),
            int(row["core_seed"]),
            int(row["gnk_seed"]),
            str(row["representation_arm"]),
            str(row["additive_core_mode"]),
            str(row["square_background"]),
            str(row["square_stratum"]),
            float(row["amplitude"]),
        )
    )
    latent.write_csv_atomic(output_dir / "validation_edge_metrics.csv", edge_rows)
    latent.write_csv_atomic(output_dir / "validation_texture_metrics.csv", texture_rows)
    latent_geometry = [
        {
            key: row.get(key)
            for key in (
                "reference_id",
                "length_bin",
                "length",
                "core_seed",
                "gnk_seed",
                "representation_arm",
                "additive_core_mode",
                "amplitude",
                "square_stratum",
                "square_background",
                "latent_epsilon_rms",
                "latent_epsilon_max_abs",
                "latent_epsilon_exact_zero_fraction",
                "latent_epsilon_near_zero_fraction",
                "latent_epsilon_abs_median",
                "latent_epsilon_abs_q95",
                "latent_epsilon_excess_kurtosis",
                "direct_tau_y_rms",
                "direct_tau_practical_resolution_fraction",
            )
        }
        for row in texture_rows
        if bool(row.get("pooled_supported_reference", False))
    ]
    latent.write_csv_atomic(
        output_dir / "latent_geometry_metrics.csv", latent_geometry
    )

    checks: list[dict[str, Any]] = [
        {
            "check_id": "shared_frozen_core_guardrail",
            "representation_arm": "shared",
            "family": "core",
            "required": True,
            "decision_relevant": True,
            "passed": bool(shared_core_pass),
            "criterion": "frozen EXP-062 core gates pass in every length stratum",
            "value": json.dumps(core_values, sort_keys=True),
        },
        {
            "check_id": "cross_arm_pairing_and_manipulation",
            "representation_arm": "all_arms",
            "family": "technical",
            "required": True,
            "decision_relevant": True,
            "passed": bool(paired_pass),
            "criterion": "all four arms share supports, base coefficient draws, and probes and pass manipulation checks",
            "value": json.dumps({
                "paired_context_count": len(pairing_rows),
                "all_contexts_pass": bool(paired_pass),
            }, sort_keys=True),
        },
    ]
    arm_flags: dict[str, dict[str, Any]] = {}
    dose_rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    if not edge.empty and not texture.empty:
        for arm in ARM_COMPLEXITY:
            arm_checks, flags, arm_doses, fit = arm_gate(
                arm=arm,
                edge=edge,
                texture=texture,
                helper=helper,
                config=config,
                technical_pass=paired_pass,
            )
            checks.extend(arm_checks)
            arm_flags[arm] = flags
            dose_rows.extend(arm_doses)
            fit_rows.append(fit)
        paired_rows, factorial_rows, superiority = paired_effects(texture, config)
    else:
        paired_rows, factorial_rows = [], []
        superiority = {arm: False for arm in RESCUE_ARMS}
        arm_flags = {
            arm: {
                "technical_pass": False,
                "edge_pass": False,
                "actual_texture_pass": False,
                "midpoint_control_pass": False,
                "directionality_pass": False,
                "length_stability_pass": False,
                "package_pass_before_paired_superiority": False,
            }
            for arm in ARM_COMPLEXITY
        }
    for arm in RESCUE_ARMS:
        summary_rows = [
            row
            for row in paired_rows
            if row["row_type"] == "comparison_summary"
            and row["comparison_arm"] == arm
        ]
        checks.append({
            "check_id": f"{arm}_paired_superiority",
            "representation_arm": arm,
            "family": "causal_rescue",
            "required": False,
            "decision_relevant": True,
            "passed": bool(superiority.get(arm, False)),
            "criterion": "median practical-fraction gain >=0.05 and paired bootstrap lower >0",
            "value": json.dumps(
                _clean_json(summary_rows[0] if summary_rows else {}),
                sort_keys=True,
            ),
        })

    latent.write_csv_atomic(output_dir / "paired_arm_effects.csv", paired_rows)
    latent.write_csv_atomic(output_dir / "factorial_effects.csv", factorial_rows)
    latent.write_csv_atomic(output_dir / "dose_response_metrics.csv", dose_rows)
    latent.write_csv_atomic(output_dir / "scale_law_fits.csv", fit_rows)
    decision = classify_decision(
        quick=bool(args.quick),
        pairing_pass=paired_pass,
        core_pass=shared_core_pass,
        arm_flags=arm_flags,
        superiority=superiority,
    )
    decision_registered = decision in {
        str(value) for value in config["decision"]["precedence"]
    }
    checks.append({
        "check_id": "registered_terminal_decision_resolved",
        "representation_arm": "all_arms",
        "family": "decision_contract",
        "required": True,
        "decision_relevant": True,
        "passed": bool(decision_registered),
        "criterion": "runner terminal decision is preregistered in config.yaml",
        "value": decision,
    })
    latent.write_csv_atomic(output_dir / "acceptance_checks.csv", checks)
    decision_rows = [
        {
            "representation_arm": arm,
            **flags,
            "paired_superiority_pass": (
                None if arm == "binary_r1" else bool(superiority.get(arm, False))
            ),
            "selected_decision": decision,
        }
        for arm, flags in arm_flags.items()
    ]
    latent.write_csv_atomic(output_dir / "decision_table.csv", decision_rows)
    required_checks = [row for row in checks if bool(row["required"])]
    rescue_decisions = {
        "DENSE_Q20_RANK1_SUFFICIENT",
        "BINARY_RANK8_SUFFICIENT",
        "ALTERNATIVE_SINGLE_FACTOR_REPAIRS_BOTH_SUFFICIENT",
        "COMBINED_DENSE_Q20_RANK8_REQUIRED",
    }
    summary = {
        **preflight,
        "core_guardrail": core_values,
        "shared_core_pass": bool(shared_core_pass),
        "pairing_and_manipulation_pass": bool(paired_pass),
        "context_task_count": len(results),
        "validation_edge_row_count": len(edge_rows),
        "validation_texture_row_count": len(texture_rows),
        "arm_decision_flags": arm_flags,
        "paired_superiority_flags": superiority,
        "required_checks_passed": int(
            sum(bool(row["passed"]) for row in required_checks)
        ),
        "required_check_count": len(required_checks),
        "decision": decision,
        "scientific_pass": (
            None if bool(args.quick) else decision in rescue_decisions
        ),
        "candidate_is_experiment_local": True,
        "public_model_changed": False,
        "independent_protein_generalization_tested": False,
        "empirical_realism_tested": False,
        "optimizer_or_accessibility_tested": False,
    }
    strict_json_atomic(output_dir / "summary.json", summary)
    print(json.dumps(_clean_json(summary), indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
