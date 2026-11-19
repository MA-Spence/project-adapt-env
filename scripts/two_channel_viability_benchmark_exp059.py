#!/usr/bin/env python3
"""EXP-059 saturated-gate mechanism calibration and confirmation.

The runner deliberately treats sequence/graphlet sampling as part of the
experimental design.  A sampling key never contains architecture, GNK routing,
or GNK weight; all treatments are therefore evaluated on the same sequences
and with the same sampled edges, plaquettes, and walk starts.

Calibration enumerates complete four-site, 20-state graphlets embedded in an
L=48 landscape.  Confirmation uses the empirical design cardinalities but not
their fitness values.  Both arms use the same sampled greedy-basin metric;
calibration additionally audits exact Hamming-shortest paths for fixed starts.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
ADAPTENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
for _path in (SRC_ROOT, SCRIPT_ROOT, ADAPTENV_ROOT):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from design_matched_support_aware_probe_exp057 import (  # noqa: E402
    allele_choices,
    choose_support_sites,
    factorial_sequences,
    load_templates,
    stable_uint64,
    support_catalog,
)
from hashed_gnk_benchmark_knob_validation import spearman  # noqa: E402

from adaptenv import FitnessLandscapeEnv, LandscapeConfig, PUBLIC_MODEL_SYSTEM  # noqa: E402


ARCH_PRODUCT = "two_channel_product"
ARCH_MIXED = "mixed_score_stability"
ARCH_GATE = "single_stability_gate"
ARCH_CLOSED = "two_channel_closed_gate"
VALID_ARCHITECTURES = {ARCH_PRODUCT, ARCH_MIXED, ARCH_GATE, ARCH_CLOSED}
VALID_ROUTINGS = {"structured", "shuffled", "random"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--arm", choices=("calibration", "confirmation"))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(dict(out[key]), value)
        else:
            out[key] = value
    return out


def load_config(path: Path, *, quick: bool = False) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("config root must be a mapping")
    if quick and isinstance(payload.get("quick"), dict):
        payload = deep_merge(payload, dict(payload["quick"]))
    payload["_runtime"] = {"quick": bool(quick)}
    return payload


def json_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def array_hash(values: np.ndarray) -> str:
    arr = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(arr.dtype).encode("ascii"))
    digest.update(json.dumps(arr.shape).encode("ascii"))
    digest.update(arr.view(np.uint8))
    return digest.hexdigest()


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._=-" else "_" for ch in str(value))


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_csv_atomic(path: Path, rows: Sequence[dict[str, Any]], columns: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    effective_columns = list(columns) if columns is not None else None
    if not rows and effective_columns is None:
        # Scientific hard stops are normal outcomes. Keep every declared CSV
        # parseable even when its stage was correctly skipped.
        effective_columns = ["status"]
    frame = pd.DataFrame(list(rows), columns=effective_columns)
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def sigmoid(values: np.ndarray | float) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    return 1.0 / (1.0 + np.exp(-np.clip(arr, -500.0, 500.0)))


def clipped_logit(values: np.ndarray) -> np.ndarray:
    arr = np.clip(np.asarray(values, dtype=np.float64), 1e-6, 1.0 - 1e-6)
    return np.log(arr) - np.log1p(-arr)


def architecture_name(value: Any) -> str:
    aliases = {
        "C": ARCH_PRODUCT,
        "C'": ARCH_PRODUCT,
        "product": ARCH_PRODUCT,
        "B": ARCH_MIXED,
        "mixed": ARCH_MIXED,
        "A": ARCH_GATE,
        "gate": ARCH_GATE,
        "closed": ARCH_CLOSED,
        "closed_gate": ARCH_CLOSED,
    }
    normalized = aliases.get(str(value), str(value))
    if normalized not in VALID_ARCHITECTURES:
        raise ValueError(f"unsupported architecture: {value}")
    return normalized


def arm_name(config: dict[str, Any], override: str | None = None) -> str:
    arm = override or config.get("arm") or config.get("experiment", {}).get("arm")
    if arm not in {"calibration", "confirmation"}:
        raise ValueError("set experiment.arm (or --arm) to calibration or confirmation")
    return str(arm)


def configured_stage_order(config: dict[str, Any], override: str | None = None) -> list[str]:
    """Return the explicit stage or the preregistered calibration->confirmation order."""

    if override is not None:
        return [arm_name(config, override)]
    configured = config.get("experiment", {}).get("stage_order", ["calibration", "confirmation"])
    stages = [str(value) for value in configured]
    if not stages or any(stage not in {"calibration", "confirmation"} for stage in stages):
        raise ValueError("experiment.stage_order must contain calibration/confirmation stages")
    if "confirmation" in stages and "calibration" in stages:
        if stages.index("confirmation") < stages.index("calibration"):
            raise ValueError("confirmation cannot precede calibration")
    return stages


def arm_config(config: dict[str, Any], arm: str) -> dict[str, Any]:
    section = config.get(arm, {})
    return dict(section) if isinstance(section, dict) else {}


def configured_weights(config: dict[str, Any]) -> list[float]:
    gen = config.get("generator", {})
    values = gen.get("weights", gen.get("gnk_epistasis_weight_grid", [0.0, 0.3, 0.6, 1.0]))
    weights = sorted({float(value) for value in values})
    if not weights or weights[0] != 0.0:
        weights = [0.0, *weights]
    return weights


def configured_architectures(config: dict[str, Any], arm: str) -> list[str]:
    section = arm_config(config, arm)
    default = [ARCH_PRODUCT, ARCH_MIXED, ARCH_GATE] if arm == "calibration" else [ARCH_PRODUCT, ARCH_MIXED]
    values = section.get("architectures", config.get("generator", {}).get("architectures", default))
    architectures = [architecture_name(value) for value in values]
    negative_control = config.get("generator", {}).get("direct_negative_control")
    if negative_control:
        control = architecture_name(negative_control)
        if control not in architectures:
            architectures.append(control)
    return architectures


def configured_routings(config: dict[str, Any], arm: str) -> list[str]:
    section = arm_config(config, arm)
    generator = config.get("generator", {})
    values = section.get(
        "support_modes",
        section.get(
            "routing_modes",
            generator.get("support_modes", generator.get("routing_modes", ["structured", "shuffled", "random"])),
        ),
    )
    routes = [str(value) for value in values]
    invalid = sorted(set(routes) - VALID_ROUTINGS)
    if invalid:
        raise ValueError(f"unsupported routing modes: {invalid}")
    return routes


def configured_seeds(config: dict[str, Any], arm: str) -> list[int]:
    defaults = [701, 702, 703, 704] if arm == "calibration" else list(range(1701, 1709))
    section = arm_config(config, arm)
    values = section.get("seeds", config.get("generator", {}).get("seeds", defaults))
    seeds = [int(value) for value in values]
    if arm == "confirmation":
        calibration_seeds = {int(value) for value in arm_config(config, "calibration").get("seeds", [701, 702, 703, 704])}
        overlap = sorted(calibration_seeds.intersection(seeds))
        if overlap:
            raise ValueError(f"confirmation seeds must be disjoint from calibration seeds: {overlap}")
    return seeds


def validate_global_config(config: dict[str, Any], arm: str, *, quick: bool = False) -> None:
    length = int(config.get("generator", {}).get("sequence_length", 48))
    if length != 48:
        raise ValueError(f"EXP-059 requires L=48 in both arms, got {length}")
    if arm == "calibration":
        section = arm_config(config, arm)
        dimension = int(section.get("dimension", 4))
        alphabet = int(section.get("alphabet_size", 20))
        if alphabet != 20:
            raise ValueError("calibration must use the native q=20 alphabet")
        if dimension != 4 and not (quick and 1 <= dimension < 4):
            raise ValueError("full calibration must enumerate complete q=20, four-site graphlets")


def _normalise_order_weights(values: Any) -> tuple[float, ...]:
    if isinstance(values, dict):
        return tuple(float(values[key]) for key in sorted(values, key=lambda item: int(item)))
    return tuple(float(value) for value in values)


def build_environment(
    config: dict[str, Any],
    *,
    seed: int,
    operating_point: dict[str, Any],
    routing_mode: str,
    gnk_weight: float,
    architecture: str = ARCH_PRODUCT,
) -> FitnessLandscapeEnv:
    """Build an AdaptEnv instance used for component evaluation or validation."""

    gen = config.get("generator", {})
    architecture = architecture_name(architecture)
    mode = "product" if architecture in {ARCH_PRODUCT, ARCH_CLOSED} else "stability_gate"
    target = "function" if architecture in {ARCH_PRODUCT, ARCH_CLOSED} else "stability"
    kwargs: dict[str, Any] = {
        "model_system": PUBLIC_MODEL_SYSTEM,
        "observed_fitness_combine_mode": mode,
        "L": 48,
        "alphabet_size": 20,
        "seed": int(seed),
        "stability_margin": float(operating_point["stability_margin"]),
        "stability_scale": float(gen.get("stability_scale", 1.0)),
        "threshold_steepness": float(operating_point["threshold_steepness"]),
        "fitness_output_scale": float(gen.get("fitness_output_scale", 10.0)),
        "gnk_epistasis_weight": float(gnk_weight),
        "gnk_terms": int(gen.get("gnk_terms", 384)),
        "gnk_min_order": int(gen.get("gnk_min_order", 2)),
        "gnk_max_order": int(gen.get("gnk_max_order", 4)),
        "gnk_order_weights": _normalise_order_weights(gen.get("gnk_order_weights", [0.35, 0.45, 0.20])),
        "gnk_support_mode": str(routing_mode),
        "gnk_seed": int(seed) + int(gen.get("gnk_seed_offset", 10000)),
        "n_epistatic_terms": 0,
        "noise_amplitude": 0.0,
        "strict_peak_lock": False,
        "generic_epistasis_target": target,
    }
    if mode == "product":
        kwargs.update(
            {
                "observed_fitness_terms": (
                    {"source": "stability_gate"},
                    {"source": "function_capacity"},
                ),
                "n_functional_dims": int(gen.get("n_functional_dims", gen.get("functional_dims", 1))),
                "functional_sigma_base": float(gen.get("functional_sigma_base", 2.0)),
                "functional_sigma_anisotropy": float(gen.get("functional_sigma_anisotropy", 0.0)),
                "functional_weight": 1.0,
            }
        )
    return FitnessLandscapeEnv(LandscapeConfig(**kwargs))


def evaluate_components(
    env: FitnessLandscapeEnv,
    sequences: np.ndarray,
    *,
    keys: Sequence[str],
    chunk_size: int,
) -> dict[str, np.ndarray]:
    chunks: dict[str, list[np.ndarray]] = {str(key): [] for key in keys}
    for start in range(0, int(sequences.shape[0]), int(chunk_size)):
        components = env.evaluate_batch_components(sequences[start : start + int(chunk_size)])
        for key in chunks:
            if key not in components:
                raise KeyError(f"AdaptEnv component path did not return {key}")
            chunks[key].append(np.asarray(components[key], dtype=np.float64))
    return {key: np.concatenate(parts) if parts else np.empty(0, dtype=np.float64) for key, parts in chunks.items()}


def derive_observable(
    architecture: str,
    *,
    weight: float,
    base_margin: np.ndarray,
    base_gate: np.ndarray,
    base_function_capacity: np.ndarray,
    raw_gnk_score: np.ndarray,
    threshold_steepness: float,
    fitness_output_scale: float,
    lambda_score: float = 1.0,
) -> dict[str, np.ndarray]:
    """Derive all treatment observables from paired base components.

    ``raw_gnk_score`` must come from an otherwise matched weight-1 environment.
    The formulas mirror ``FitnessLandscapeEnv._evaluate_components_batch`` and
    the EXP-058 mixed score/stability construction.
    """

    architecture = architecture_name(architecture)
    margin = np.asarray(base_margin, dtype=np.float64)
    gate0 = np.asarray(base_gate, dtype=np.float64)
    capacity0 = np.asarray(base_function_capacity, dtype=np.float64)
    gnk = float(weight) * np.asarray(raw_gnk_score, dtype=np.float64)
    scale = float(fitness_output_scale)
    steepness = float(threshold_steepness)
    if architecture in {ARCH_PRODUCT, ARCH_CLOSED}:
        capacity = sigmoid(clipped_logit(capacity0) + gnk)
        gate = gate0
        fitness = scale * gate * capacity
        effective_margin = margin
    elif architecture == ARCH_GATE:
        effective_margin = margin + gnk
        gate = sigmoid(steepness * effective_margin)
        capacity = capacity0
        fitness = scale * gate
    else:
        lambda_score = float(lambda_score)
        stability_weight = 1.0 - lambda_score
        effective_margin = margin + stability_weight * gnk
        gate = sigmoid(steepness * effective_margin)
        capacity = capacity0
        fitness = scale * gate + scale * lambda_score * gnk
    return {
        "fitness": np.asarray(fitness, dtype=np.float64),
        "stability_gate": np.asarray(gate, dtype=np.float64),
        "effective_stability_margin": np.asarray(effective_margin, dtype=np.float64),
        "effective_functional_capacity": np.asarray(capacity, dtype=np.float64),
        "weighted_gnk_score": np.asarray(gnk, dtype=np.float64),
    }


@dataclass(frozen=True)
class OperatingPoint:
    operating_point_id: str
    stability_margin: float
    threshold_steepness: float

    def as_mapping(self) -> dict[str, Any]:
        return {
            "operating_point_id": self.operating_point_id,
            "stability_margin": self.stability_margin,
            "threshold_steepness": self.threshold_steepness,
        }


def _point_from_mapping(row: dict[str, Any], index: int = 0) -> OperatingPoint:
    margin = float(row.get("stability_margin", row.get("margin", 10.0)))
    steepness = float(row.get("threshold_steepness", row.get("steepness", 4.0)))
    point_id = str(row.get("operating_point_id", row.get("id", f"op{index:03d}_m{margin:g}_k{steepness:g}")))
    return OperatingPoint(point_id, margin, steepness)


def operating_points(config: dict[str, Any], arm: str, project_root: Path) -> list[OperatingPoint]:
    section = arm_config(config, arm)
    if arm == "calibration":
        rows = section.get("operating_points", config.get("generator", {}).get("operating_points"))
        if rows is None:
            margins = section.get("stability_margin_grid", config.get("generator", {}).get("stability_margin_grid", [14.0, 18.0, 22.0]))
            steepnesses = section.get("threshold_steepness_grid", config.get("generator", {}).get("threshold_steepness_grid", [2.0, 4.0]))
            rows = [
                {"stability_margin": margin, "threshold_steepness": steepness}
                for margin in margins
                for steepness in steepnesses
            ]
        points = [_point_from_mapping(dict(row), index) for index, row in enumerate(rows)]
        return points

    frozen = section.get("frozen_operating_point")
    frozen_path = section.get("frozen_operating_point_path")
    if frozen is None and frozen_path:
        path = Path(str(frozen_path))
        path = path if path.is_absolute() else project_root / path
        frozen = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(frozen, dict) and isinstance(frozen.get("selected_operating_point"), dict):
        frozen = frozen["selected_operating_point"]
    if frozen == "from_calibration":
        raise ValueError("confirmation frozen point must be supplied by the preceding calibration stage")
    if not isinstance(frozen, dict):
        raise ValueError("confirmation requires confirmation.frozen_operating_point or _path")
    return [_point_from_mapping(dict(frozen), 0)]


@dataclass(frozen=True)
class GraphletDesign:
    arm: str
    operating_point_id: str
    seed: int
    dataset_id: str
    graphlet_index: int
    anchor_regime: str
    sampling_group_id: str
    graphlet_id: str
    anchor_sequence: tuple[int, ...]
    sites: tuple[int, ...]
    choices: tuple[tuple[int, ...], ...]
    shape: tuple[int, ...]
    target_margin: float
    actual_margin: float
    gate_occupancy: float
    effective_capacity_sd: float
    effective_capacity_interior_fraction: float
    dead_gate_fraction: float
    product_reconstruction_error: float
    raw_gnk_score_sd: float
    closed_control_low_gate_fraction: float
    preflight_pass: bool
    rejection_reason: str

    @property
    def design_hash(self) -> str:
        return json_hash(asdict(self))


def _template_specs(config: dict[str, Any], arm: str, project_root: Path) -> list[dict[str, Any]]:
    if arm == "calibration":
        section = arm_config(config, arm)
        dimension = int(section.get("dimension", 4))
        alphabet = int(section.get("alphabet_size", 20))
        return [
            {
                "dataset_id": f"embedded_q{alphabet}_d{dimension}",
                "dimension": dimension,
                "alphabet_cardinalities": [alphabet] * dimension,
                "max_alphabet_cardinality": alphabet,
                "n_expected_nodes": alphabet**dimension,
            }
        ]
    return load_templates(config, project_root)


def _graphlet_count(config: dict[str, Any], arm: str, template: dict[str, Any]) -> int:
    section = arm_config(config, arm)
    if arm == "calibration":
        return int(section.get("graphlets_per_seed", section.get("graphlets", 1)))
    configured = section.get("graphlets_per_template", 1)
    if isinstance(configured, dict):
        return int(configured.get(str(template["dataset_id"]), configured.get("default", 1)))
    return int(configured)


def prepare_designs_for_operating_point(
    config: dict[str, Any],
    *,
    arm: str,
    point: OperatingPoint,
    seeds: Sequence[int],
    project_root: Path,
) -> tuple[list[GraphletDesign], list[dict[str, Any]]]:
    """Construct one frozen viable graphlet per preregistered base-design unit.

    The basis environment is deliberately independent of the candidate operating
    point.  In particular, the sampling key omits the operating-point ID, and the
    designed-peak background, sites, and allele choices cannot change during
    calibration selection.
    """

    templates = _template_specs(config, arm, project_root)
    designs: list[GraphletDesign] = []
    support_rows: list[dict[str, Any]] = []
    for seed in seeds:
        basis_point = {"stability_margin": 0.0, "threshold_steepness": 1.0}
        base_env = build_environment(
            config,
            seed=int(seed),
            operating_point=basis_point,
            routing_mode="structured",
            gnk_weight=0.0,
        )
        structured_env = build_environment(
            config,
            seed=int(seed),
            operating_point=basis_point,
            routing_mode="structured",
            gnk_weight=1.0,
        )
        catalog_id = f"{arm}__seed{int(seed)}__structured"
        terms = support_catalog(structured_env, catalog_id)
        support_rows.extend(
            row
            | {
                "arm": arm,
                "operating_point_id": point.operating_point_id,
                "seed": int(seed),
                "routing_mode": "structured",
            }
            for row in terms
        )
        for template in templates:
            dataset_id = str(template["dataset_id"])
            dimension = int(template["dimension"])
            cardinalities = [int(value) for value in template["alphabet_cardinalities"]]
            for graphlet_index in range(_graphlet_count(config, arm, template)):
                sampling_group_id = (
                    f"{arm}__seed{int(seed)}"
                    f"__{dataset_id}__g{int(graphlet_index):03d}"
                )
                site_seed = (
                    int(seed)
                    + int(config.get("probe", {}).get("site_seed_offset", 51000))
                    + int(graphlet_index)
                ) % (2**32)
                sites = choose_support_sites(
                    terms,
                    dimension=dimension,
                    graphlet_index=int(graphlet_index),
                    sequence_length=48,
                    seed=site_seed,
                )
                anchor = np.asarray(getattr(base_env, "_designed_peak"), dtype=np.int64)
                allele_seed = (
                    int(seed)
                    + int(config.get("probe", {}).get("allele_seed_offset", 61000))
                    + int(graphlet_index)
                ) % (2**32)
                if arm == "calibration":
                    choices_list = [list(range(20)) for _ in sites]
                else:
                    selection = str(config.get("probe", {}).get("allele_selection", "hashed_including_anchor"))
                    if selection == "stability_ranked_including_anchor":
                        choices_list = []
                        for offset, site in enumerate(sites):
                            count = min(int(cardinalities[offset]), 20)
                            anchor_allele = int(anchor[int(site)])
                            if count >= 20:
                                choices_list.append(list(range(20)))
                                continue
                            alternatives = sorted(
                                [allele for allele in range(20) if allele != anchor_allele],
                                key=lambda allele: (
                                    -float(base_env.stability_matrix[int(site), int(allele)]),
                                    stable_uint64("stability-allele-tie", allele_seed, int(site), int(allele)),
                                ),
                            )
                            choices_list.append([anchor_allele, *alternatives[: max(count - 1, 0)]])
                    elif selection == "hashed_including_anchor":
                        choices_list = allele_choices(
                            anchor,
                            sites,
                            cardinalities,
                            alphabet_size=20,
                            seed=allele_seed,
                        )
                    else:
                        raise ValueError(f"unsupported probe.allele_selection: {selection}")
                choices = tuple(tuple(int(value) for value in axis) for axis in choices_list)
                shape = tuple(len(axis) for axis in choices)
                graphlet_id = f"{arm}__{point.operating_point_id}__seed{int(seed)}__{dataset_id}__g{int(graphlet_index):03d}"
                designs.append(
                    GraphletDesign(
                        arm=arm,
                        operating_point_id=point.operating_point_id,
                        seed=int(seed),
                        dataset_id=dataset_id,
                        graphlet_index=int(graphlet_index),
                        anchor_regime="viable_designed_peak",
                        sampling_group_id=sampling_group_id,
                        graphlet_id=graphlet_id,
                        anchor_sequence=tuple(int(value) for value in anchor),
                        sites=tuple(int(value) for value in sites),
                        choices=choices,
                        shape=shape,
                        target_margin=float(point.stability_margin),
                        actual_margin=math.nan,
                        gate_occupancy=math.nan,
                        effective_capacity_sd=math.nan,
                        effective_capacity_interior_fraction=math.nan,
                        dead_gate_fraction=math.nan,
                        product_reconstruction_error=math.nan,
                        raw_gnk_score_sd=math.nan,
                        closed_control_low_gate_fraction=math.nan,
                        preflight_pass=True,
                        rejection_reason="",
                    )
                )
    return designs, support_rows


def materialize_design(design: GraphletDesign) -> np.ndarray:
    sequences, shape = factorial_sequences(
        np.asarray(design.anchor_sequence, dtype=np.int64),
        design.sites,
        [list(axis) for axis in design.choices],
    )
    if tuple(shape) != tuple(design.shape):
        raise RuntimeError(f"shape drift for {design.graphlet_id}: {shape} != {design.shape}")
    return sequences


def dead_control_sequences(config: dict[str, Any], design: GraphletDesign) -> np.ndarray:
    """Draw the frozen, treatment-independent dead-state control panel."""

    count = int(config.get("probe", {}).get("dead_sequence_samples", 4096))
    alphabet = int(config.get("generator", {}).get("alphabet_size", 20))
    rng = np.random.default_rng(
        stable_uint64("dead-controls", design.sampling_group_id) % (2**63 - 1)
    )
    return rng.integers(0, alphabet, size=(count, 48), dtype=np.int64)


def closed_control_point(config: dict[str, Any], point: OperatingPoint) -> OperatingPoint:
    """The direct closed-gate control has designed-peak achieved margin -2."""

    preflight = config.get("gates", {}).get("preflight", {})
    target = float(
        preflight.get(
            "closed_control_target_effective_margin",
            config.get("generator", {}).get("closed_gate_target_effective_margin", -2.0),
        )
    )
    return OperatingPoint(
        operating_point_id=f"{point.operating_point_id}__closed_m{target:g}",
        stability_margin=target,
        threshold_steepness=point.threshold_steepness,
    )


def component_preflight(
    config: dict[str, Any],
    designs: Sequence[GraphletDesign],
    *,
    point: OperatingPoint,
) -> list[GraphletDesign]:
    """Hard-check primary product gate and effective function-capacity paths.

    Only component quantities are used here; no ruggedness/topology outcome can
    influence operating-point selection.
    """

    probe = config.get("probe", {})
    gates = config.get("gates", {}).get("preflight", {})
    chunk_size = int(probe.get("eval_chunk_size", 32768))
    weights = configured_weights(config)
    viable_threshold = float(gates.get("viable_gate_threshold", 0.90))
    occupancy_min = float(gates.get("viable_graphlet_min_occupancy", 0.90))
    dead_threshold = float(gates.get("dead_gate_threshold", 0.10))
    dead_fraction_min = float(gates.get("dead_sample_min_fraction", 0.90))
    capacity_sd_min = float(
        gates.get("effective_capacity_min_sd", probe.get("effective_capacity_sd_min", 0.05))
    )
    interior_min = float(
        gates.get(
            "effective_capacity_min_interior_fraction",
            probe.get("effective_capacity_interior_fraction_min", 0.70),
        )
    )
    if bool(config.get("_runtime", {}).get("quick", False)):
        # Explicit smoke-only relaxation; full runs retain the preregistered floor.
        interior_min = min(
            interior_min,
            float(config.get("quick", {}).get("preflight_capacity_min_interior_fraction", 0.60)),
        )
    interior_low = float(
        gates.get("effective_capacity_interior_low", probe.get("effective_capacity_interior_low", 0.05))
    )
    interior_high = float(
        gates.get("effective_capacity_interior_high", probe.get("effective_capacity_interior_high", 0.95))
    )
    reconstruction_atol = float(gates.get("product_reconstruction_atol", 1e-8))
    raw_gnk_sd_min = float(gates.get("raw_gnk_score_min_sd", 1e-3))
    out: list[GraphletDesign] = []
    for design in designs:
        sequences = materialize_design(design)
        dead_sequences = dead_control_sequences(config, design)
        reasons: list[str] = []
        capacity_sds: list[float] = []
        interior_fractions: list[float] = []
        reconstruction_errors: list[float] = []
        gate_occupancies: list[float] = []
        dead_fractions: list[float] = []
        raw_gnk_sds: list[float] = []
        actual_margin = math.nan
        closed_control_low_fraction = math.nan
        for weight in weights:
            env = build_environment(
                config,
                seed=design.seed,
                operating_point=point.as_mapping(),
                routing_mode="structured",
                gnk_weight=float(weight),
            )
            components = evaluate_components(
                env,
                sequences,
                keys=(
                    "fitness",
                    "effective_stability_margin_used",
                    "stability_gate",
                    "effective_functional_capacity",
                    "gnk_epistatic_score",
                ),
                chunk_size=chunk_size,
            )
            gate = components["stability_gate"]
            viable = gate >= viable_threshold
            occupancy = float(np.mean(viable))
            capacity = components["effective_functional_capacity"]
            capacity_sd = float(np.std(capacity[viable])) if np.any(viable) else math.nan
            interior_fraction = (
                float(np.mean((capacity[viable] > interior_low) & (capacity[viable] < interior_high)))
                if np.any(viable)
                else math.nan
            )
            reconstruction = (
                float(config.get("generator", {}).get("fitness_output_scale", 10.0))
                * gate
                * capacity
            )
            reconstruction_error = float(np.max(np.abs(components["fitness"] - reconstruction)))
            dead_gate = evaluate_components(
                env,
                dead_sequences,
                keys=("stability_gate",),
                chunk_size=chunk_size,
            )["stability_gate"]
            dead_fraction = float(np.mean(dead_gate <= dead_threshold))
            if weight == 0.0:
                anchor_components = evaluate_components(
                    env,
                    np.asarray(design.anchor_sequence, dtype=np.int64)[None, :],
                    keys=("effective_stability_margin_used",),
                    chunk_size=1,
                )
                actual_margin = float(anchor_components["effective_stability_margin_used"][0])
            gate_occupancies.append(occupancy)
            capacity_sds.append(capacity_sd)
            interior_fractions.append(interior_fraction)
            reconstruction_errors.append(reconstruction_error)
            dead_fractions.append(dead_fraction)
            if np.isclose(float(weight), 1.0):
                raw_gnk_sd = float(np.std(components["gnk_epistatic_score"]))
                raw_gnk_sds.append(raw_gnk_sd)
                if raw_gnk_sd < raw_gnk_sd_min:
                    reasons.append("raw_gnk_score_sd_below_floor_w1")
            suffix = f"_w{float(weight):g}"
            if occupancy < occupancy_min:
                reasons.append(f"viable_gate_occupancy_below_floor{suffix}")
            if not np.isfinite(capacity_sd) or capacity_sd < capacity_sd_min:
                reasons.append(f"effective_capacity_sd_below_floor{suffix}")
            if not np.isfinite(interior_fraction) or interior_fraction < interior_min:
                reasons.append(f"effective_capacity_interior_fraction_below_floor{suffix}")
            if reconstruction_error > reconstruction_atol:
                reasons.append(f"product_reconstruction_error_above_tolerance{suffix}")
            if dead_fraction < dead_fraction_min:
                reasons.append(f"dead_gate_fraction_below_floor{suffix}")
        passed = bool(not reasons)
        closed_point = closed_control_point(config, point)
        closed_env = build_environment(
            config,
            seed=design.seed,
            operating_point=closed_point.as_mapping(),
            routing_mode="structured",
            gnk_weight=0.0,
            architecture=ARCH_PRODUCT,
        )
        closed_components = evaluate_components(
            closed_env,
            sequences,
            keys=("stability_gate", "effective_stability_margin_used"),
            chunk_size=chunk_size,
        )
        closed_cfg = config.get("gates", {}).get("preflight", {})
        closed_threshold = float(
            closed_cfg.get("closed_control_gate_threshold", closed_cfg.get("dead_gate_threshold", 0.10))
        )
        closed_min_fraction = float(closed_cfg.get("closed_control_min_low_gate_fraction", 0.90))
        closed_control_low_fraction = float(
            np.mean(closed_components["stability_gate"] <= closed_threshold)
        )
        closed_anchor_margin = float(
            evaluate_components(
                closed_env,
                np.asarray(design.anchor_sequence, dtype=np.int64)[None, :],
                keys=("effective_stability_margin_used",),
                chunk_size=1,
            )["effective_stability_margin_used"][0]
        )
        if abs(closed_anchor_margin - float(closed_point.stability_margin)) > 1e-8:
            reasons.append("closed_control_target_margin_not_achieved")
        if closed_control_low_fraction < closed_min_fraction:
            reasons.append("closed_control_low_gate_fraction_below_floor")
        passed = bool(not reasons)
        out.append(
            replace(
                design,
                actual_margin=actual_margin,
                gate_occupancy=float(min(gate_occupancies)) if gate_occupancies else math.nan,
                effective_capacity_sd=float(min(capacity_sds)) if capacity_sds else math.nan,
                effective_capacity_interior_fraction=(
                    float(min(interior_fractions)) if interior_fractions else math.nan
                ),
                dead_gate_fraction=float(min(dead_fractions)) if dead_fractions else math.nan,
                product_reconstruction_error=(
                    float(max(reconstruction_errors)) if reconstruction_errors else math.nan
                ),
                raw_gnk_score_sd=float(min(raw_gnk_sds)) if raw_gnk_sds else math.nan,
                closed_control_low_gate_fraction=closed_control_low_fraction,
                preflight_pass=passed,
                rejection_reason=";".join(sorted(set(reasons))),
            )
        )
    return out


def select_operating_point(
    points: Sequence[OperatingPoint],
    designs_by_point: dict[str, list[GraphletDesign]],
) -> OperatingPoint | None:
    """Choose deterministically using preflight/component checks only."""

    for point in points:
        designs = designs_by_point.get(point.operating_point_id, [])
        if designs and all(design.preflight_pass for design in designs):
            return point
    return None


@dataclass(frozen=True)
class MetricPlan:
    edge_source: tuple[int, ...]
    edge_target: tuple[int, ...]
    pair_nodes: tuple[tuple[int, int, int, int], ...]
    basin_starts: tuple[int, ...]
    path_starts: tuple[int, ...]

    @property
    def plan_hash(self) -> str:
        return json_hash(asdict(self))


def _coordinates(indices: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    return np.column_stack(np.unravel_index(np.asarray(indices, dtype=np.int64), shape))


def _indices(coordinates: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    return np.ravel_multi_index(np.asarray(coordinates, dtype=np.int64).T, shape)


def make_metric_plan(config: dict[str, Any], *, shape: tuple[int, ...], sampling_group_id: str) -> MetricPlan:
    descriptors = config.get("descriptors", {})
    probe = config.get("probe", {})
    n_nodes = int(np.prod(shape, dtype=np.int64))
    rng = np.random.default_rng(stable_uint64("metric-plan", sampling_group_id) % (2**63 - 1))

    n_edges = min(
        int(descriptors.get("single_step_effect_samples", probe.get("single_step_effect_samples", 20000))),
        n_nodes * sum(k - 1 for k in shape),
    )
    edge_source = rng.integers(0, n_nodes, size=n_edges, dtype=np.int64)
    edge_coords = _coordinates(edge_source, shape)
    edge_target_coords = edge_coords.copy()
    edge_axes = rng.integers(0, len(shape), size=n_edges)
    for row, axis in enumerate(edge_axes):
        old = int(edge_target_coords[row, axis])
        alternative = int(rng.integers(0, shape[int(axis)] - 1))
        edge_target_coords[row, axis] = alternative + int(alternative >= old)
    edge_target = _indices(edge_target_coords, shape)

    n_pairs = int(descriptors.get("pairwise_epistasis_samples", probe.get("pairwise_epistasis_samples", 8000)))
    pair_rows: list[tuple[int, int, int, int]] = []
    if len(shape) >= 2:
        for _ in range(n_pairs):
            left, right = [int(value) for value in rng.choice(len(shape), size=2, replace=False)]
            base_coord = np.asarray([int(rng.integers(0, size)) for size in shape], dtype=np.int64)
            left_alt = int(rng.integers(0, shape[left] - 1))
            left_alt += int(left_alt >= int(base_coord[left]))
            right_alt = int(rng.integers(0, shape[right] - 1))
            right_alt += int(right_alt >= int(base_coord[right]))
            c00 = base_coord.copy()
            c10 = base_coord.copy()
            c01 = base_coord.copy()
            c11 = base_coord.copy()
            c10[left] = left_alt
            c01[right] = right_alt
            c11[left] = left_alt
            c11[right] = right_alt
            nodes = np.ravel_multi_index(np.vstack((c00, c10, c01, c11)).T, shape)
            pair_rows.append(tuple(int(value) for value in nodes))

    basin_count = min(
        int(
            descriptors.get(
                "basin_starts",
                descriptors.get("reachability_starts", probe.get("reachability_starts", 256)),
            )
        ),
        n_nodes,
    )
    path_count = min(
        int(
            descriptors.get(
                "exact_path_starts",
                probe.get("exact_shortest_path_starts", probe.get("exact_shortest_path_starts_per_distance", 64)),
            )
        ),
        n_nodes,
    )
    basin_starts = tuple(int(value) for value in rng.choice(n_nodes, size=basin_count, replace=False))
    path_starts = tuple(int(value) for value in rng.choice(n_nodes, size=path_count, replace=False))
    return MetricPlan(
        edge_source=tuple(int(value) for value in edge_source),
        edge_target=tuple(int(value) for value in edge_target),
        pair_nodes=tuple(pair_rows),
        basin_starts=basin_starts,
        path_starts=path_starts,
    )


def exact_interaction_fraction(values: np.ndarray, shape: tuple[int, ...]) -> float:
    """Uniform factorial variance not explained by intercept + main effects.

    This is the q-ary interaction fraction of order >=2.  It is intentionally
    not labelled a binary Walsh order>=3 fraction.
    """

    tensor = np.asarray(values, dtype=np.float64).reshape(shape)
    centered = tensor - float(np.mean(tensor))
    total = float(np.sum(centered * centered))
    if total <= 1e-18:
        return 0.0
    grand = float(np.mean(tensor))
    main_ss = 0.0
    for axis in range(len(shape)):
        reduction_axes = tuple(index for index in range(len(shape)) if index != axis)
        effect = np.mean(tensor, axis=reduction_axes) - grand
        repeats = int(np.prod([shape[index] for index in reduction_axes], dtype=np.int64))
        main_ss += repeats * float(np.sum(effect * effect))
    return float(np.clip(1.0 - main_ss / total, 0.0, 1.0))


def masked_strict_local_maxima(
    values: np.ndarray,
    shape: tuple[int, ...],
    viable: np.ndarray,
    *,
    tolerance: float,
) -> tuple[int, float]:
    arr = np.asarray(values, dtype=np.float64).reshape(shape)
    mask = np.asarray(viable, dtype=bool).reshape(shape)
    is_max = mask.copy()
    # Each axis is a collection of complete one-mutant fibres.  Comparing each
    # node with the largest *other* viable allele vectorises the complete
    # neighbourhood scan and preserves the fixed absolute tolerance.
    for axis, size in enumerate(shape):
        fibre_values = np.moveaxis(arr, axis, 0).reshape(int(size), -1)
        fibre_mask = np.moveaxis(mask, axis, 0).reshape(int(size), -1)
        masked_values = np.where(fibre_mask, fibre_values, -np.inf)
        if int(size) == 1:
            other_max = np.full_like(masked_values, -np.inf)
        else:
            order = np.argsort(masked_values, axis=0)
            top_index = order[-1]
            top_value = np.take_along_axis(masked_values, top_index[None, :], axis=0)[0]
            second_value = np.take_along_axis(masked_values, order[-2][None, :], axis=0)[0]
            allele_index = np.arange(int(size), dtype=np.int64)[:, None]
            other_max = np.where(allele_index == top_index[None, :], second_value[None, :], top_value[None, :])
        axis_is_max = fibre_mask & (fibre_values > other_max + float(tolerance))
        axis_is_max = np.moveaxis(axis_is_max.reshape((int(size),) + np.moveaxis(arr, axis, 0).shape[1:]), 0, axis)
        is_max &= axis_is_max
    count = int(np.sum(is_max))
    denominator = int(np.sum(mask))
    return count, float(count / denominator) if denominator else math.nan


def sampled_effect_metrics(
    values: np.ndarray,
    viable: np.ndarray,
    plan: MetricPlan,
    *,
    tolerance: float,
) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    mask = np.asarray(viable, dtype=bool)
    source = np.asarray(plan.edge_source, dtype=np.int64)
    target = np.asarray(plan.edge_target, dtype=np.int64)
    keep = mask[source] & mask[target]
    effects = arr[target[keep]] - arr[source[keep]]
    pair_array = np.asarray(plan.pair_nodes, dtype=np.int64)
    sign_count = reciprocal_count = 0
    pair_count = 0
    if pair_array.size:
        keep_pairs = np.all(mask[pair_array], axis=1)
        for n00, n10, n01, n11 in pair_array[keep_pairs]:
            left_low = float(arr[n10] - arr[n00])
            left_high = float(arr[n11] - arr[n01])
            right_low = float(arr[n01] - arr[n00])
            right_high = float(arr[n11] - arr[n10])
            left_flip = (
                abs(left_low) > tolerance
                and abs(left_high) > tolerance
                and left_low * left_high < 0.0
            )
            right_flip = (
                abs(right_low) > tolerance
                and abs(right_high) > tolerance
                and right_low * right_high < 0.0
            )
            sign_count += int(left_flip or right_flip)
            reciprocal_count += int(left_flip and right_flip)
            pair_count += 1
    return {
        "single_step_effect_sd": float(np.std(effects)) if effects.size else math.nan,
        "beneficial_fraction": float(np.mean(effects > tolerance)) if effects.size else math.nan,
        "neutral_fraction": float(np.mean(np.abs(effects) <= tolerance)) if effects.size else math.nan,
        "deleterious_fraction": float(np.mean(effects < -tolerance)) if effects.size else math.nan,
        "sign_epistasis_fraction": float(sign_count / pair_count) if pair_count else math.nan,
        "reciprocal_sign_epistasis_fraction": float(reciprocal_count / pair_count) if pair_count else math.nan,
        "n_sampled_viable_edges": float(effects.size),
        "n_sampled_viable_plaquettes": float(pair_count),
    }


def greedy_basin_metrics(
    values: np.ndarray,
    shape: tuple[int, ...],
    viable: np.ndarray,
    starts: Sequence[int],
    *,
    tolerance: float,
) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    mask = np.asarray(viable, dtype=bool)
    viable_indices = np.flatnonzero(mask)
    if viable_indices.size == 0:
        return {
            "basin_convergence_fraction": math.nan,
            "mean_greedy_walk_length": math.nan,
            "basin_target_node": math.nan,
            "n_viable_basin_starts": 0.0,
        }
    target = int(viable_indices[int(np.argmax(arr[viable_indices]))])
    successes = 0
    lengths: list[int] = []
    used = 0
    max_steps = max(4 * len(shape) * max(shape), 1)
    for start in starts:
        if not mask[int(start)]:
            continue
        used += 1
        current = int(start)
        steps = 0
        for _ in range(max_steps):
            coordinate = list(np.unravel_index(current, shape))
            best = current
            best_value = float(arr[current])
            for axis, size in enumerate(shape):
                old = coordinate[axis]
                for allele in range(size):
                    if allele == old:
                        continue
                    coordinate[axis] = allele
                    neighbor = int(np.ravel_multi_index(tuple(coordinate), shape))
                    if mask[neighbor] and float(arr[neighbor]) > best_value + tolerance:
                        best = neighbor
                        best_value = float(arr[neighbor])
                coordinate[axis] = old
            if best == current:
                break
            current = best
            steps += 1
        successes += int(current == target)
        lengths.append(steps)
    return {
        "basin_convergence_fraction": float(successes / used) if used else math.nan,
        "mean_greedy_walk_length": float(np.mean(lengths)) if lengths else math.nan,
        "basin_target_node": float(target),
        "n_viable_basin_starts": float(used),
    }


def exact_shortest_path_audit(
    values: np.ndarray,
    shape: tuple[int, ...],
    viable: np.ndarray,
    starts: Sequence[int],
    *,
    tolerance: float,
) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    mask = np.asarray(viable, dtype=bool)
    viable_indices = np.flatnonzero(mask)
    if viable_indices.size == 0:
        return {
            "accessible_shortest_path_fraction": math.nan,
            "accessible_shortest_path_count": 0.0,
            "shortest_path_count": 0.0,
            "n_path_starts": 0.0,
        }
    target = int(viable_indices[int(np.argmax(arr[viable_indices]))])
    target_coordinate = np.asarray(np.unravel_index(target, shape), dtype=np.int64)
    fractions: list[float] = []
    accessible_total = 0
    path_total = 0
    for start in starts:
        start = int(start)
        if start == target or not mask[start]:
            continue
        coordinate = np.asarray(np.unravel_index(start, shape), dtype=np.int64)
        axes = [axis for axis in range(len(shape)) if int(coordinate[axis]) != int(target_coordinate[axis])]
        total = math.factorial(len(axes))
        accessible = 0
        for order in itertools.permutations(axes):
            node_coordinate = coordinate.copy()
            node = start
            previous = float(arr[node])
            ok = True
            for axis in order:
                node_coordinate[int(axis)] = target_coordinate[int(axis)]
                node = int(np.ravel_multi_index(tuple(node_coordinate), shape))
                current = float(arr[node])
                if not mask[node] or current <= previous + tolerance:
                    ok = False
                    break
                previous = current
            accessible += int(ok)
        accessible_total += accessible
        path_total += total
        fractions.append(float(accessible / total) if total else math.nan)
    return {
        "accessible_shortest_path_fraction": float(np.mean(fractions)) if fractions else math.nan,
        "accessible_shortest_path_count": float(accessible_total),
        "shortest_path_count": float(path_total),
        "n_path_starts": float(len(fractions)),
    }


def metric_bundle(
    values: np.ndarray,
    shape: tuple[int, ...],
    gate: np.ndarray,
    plan: MetricPlan,
    *,
    tolerance: float,
    exact_paths: bool,
) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    viable = np.asarray(gate, dtype=np.float64) >= 0.9
    maxima_count, maxima_density = masked_strict_local_maxima(
        arr,
        shape,
        viable,
        tolerance=tolerance,
    )
    interaction_fraction = exact_interaction_fraction(arr, shape)
    out = {
        "node_count": float(arr.size),
        "viable_node_count": float(np.sum(viable)),
        "gate_occupancy": float(np.mean(viable)) if viable.size else math.nan,
        "fitness_mean": float(np.mean(arr)) if arr.size else math.nan,
        "fitness_sd": float(np.std(arr)) if arr.size else math.nan,
        "fitness_min": float(np.min(arr)) if arr.size else math.nan,
        "fitness_max": float(np.max(arr)) if arr.size else math.nan,
        "fitness_median": float(np.median(arr)) if arr.size else math.nan,
        "fitness_q95": float(np.quantile(arr, 0.95)) if arr.size else math.nan,
        "viable_fitness_median": float(np.median(arr[viable])) if np.any(viable) else math.nan,
        "viable_fitness_sd": float(np.std(arr[viable])) if np.any(viable) else math.nan,
        "local_maxima_count": float(maxima_count),
        "local_maxima_density": maxima_density,
        "high_order_fraction": interaction_fraction,
        "nonadditive_variance_fraction": interaction_fraction,
        "interaction_fraction_ge2": interaction_fraction,
        "fixed_tolerance": float(tolerance),
        **sampled_effect_metrics(arr, viable, plan, tolerance=tolerance),
        **greedy_basin_metrics(arr, shape, viable, plan.basin_starts, tolerance=tolerance),
    }
    if exact_paths:
        out.update(exact_shortest_path_audit(arr, shape, viable, plan.path_starts, tolerance=tolerance))
    return out


def fixed_fitness_tolerance(config: dict[str, Any]) -> float:
    descriptors = config.get("descriptors", {})
    return float(descriptors.get("fixed_fitness_tolerance", descriptors.get("fixed_tolerance", 1e-6)))


def working_weight(config: dict[str, Any]) -> float:
    generator = config.get("generator", {})
    return float(generator.get("working_weight", config.get("gates", {}).get("working_weight", 1.0)))


def _titration_summary(
    config: dict[str, Any],
    *,
    architecture: str,
    weight: float,
    point: OperatingPoint,
    capacity_at_background: float,
    raw_gnk_at_background: float,
) -> dict[str, float | bool]:
    """Fixed-background stability titration with the function value held fixed."""

    gate_cfg = config.get("gates", {}).get("g1_viability", {})
    titration_cfg = gate_cfg.get("stability_titration", {})
    margins = np.asarray(
        titration_cfg.get(
            "achieved_margins",
            titration_cfg.get("offsets", [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0]),
        ),
        dtype=np.float64,
    )
    transition_window = titration_cfg.get("transition_window", [-1.0, 1.0])
    saturated_window = titration_cfg.get("saturated_window", [2.0, 4.0])
    transition_mask = (margins >= float(transition_window[0])) & (margins <= float(transition_window[1]))
    saturated_mask = (margins >= float(saturated_window[0])) & (margins <= float(saturated_window[1]))
    gate = sigmoid(float(point.threshold_steepness) * margins)
    scale = float(config.get("generator", {}).get("fitness_output_scale", 10.0))
    raw = float(raw_gnk_at_background)
    if architecture in {ARCH_PRODUCT, ARCH_CLOSED}:
        capacity = float(sigmoid(clipped_logit(np.asarray([capacity_at_background])) + float(weight) * raw)[0])
        values = scale * gate * capacity
    elif architecture == ARCH_MIXED:
        lambda_score = float(config.get("generator", {}).get("lambda_score", 1.0))
        stability_weight = 1.0 - lambda_score
        shifted_gate = sigmoid(float(point.threshold_steepness) * (margins + stability_weight * float(weight) * raw))
        values = scale * shifted_gate + scale * lambda_score * float(weight) * raw
    else:
        values = scale * sigmoid(float(point.threshold_steepness) * (margins + float(weight) * raw))

    def slope(mask: np.ndarray) -> float:
        if int(np.sum(mask)) < 2:
            return math.nan
        return float(np.polyfit(margins[mask], np.asarray(values)[mask], 1)[0])

    transition_slope = slope(transition_mask)
    saturated_slope = slope(saturated_mask)
    ratio = (
        float(saturated_slope / transition_slope)
        if np.isfinite(transition_slope) and transition_slope > 0.0 and np.isfinite(saturated_slope)
        else math.nan
    )
    min_transition = float(titration_cfg.get("min_transition_slope", 1e-4))
    min_saturated = float(titration_cfg.get("min_saturated_slope", 0.0))
    max_ratio = float(titration_cfg.get("max_saturated_to_transition_slope_ratio", 0.20))
    passed = bool(
        np.isfinite(transition_slope)
        and np.isfinite(saturated_slope)
        and np.isfinite(ratio)
        and transition_slope >= min_transition
        and saturated_slope >= min_saturated - 1e-12
        and ratio <= max_ratio
    )
    return {
        "titration_transition_slope": transition_slope,
        "titration_saturated_slope": saturated_slope,
        "titration_saturated_to_transition_ratio": ratio,
        "titration_pass": passed,
    }


def _treatment_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        row["stage"],
        row["operating_point_id"],
        int(row["seed"]),
        row["dataset_id"],
        int(row["graphlet_index"]),
        row["architecture"],
        row["support_mode"],
        float(row["gnk_weight"]),
    )


def evaluate_design(
    config: dict[str, Any],
    design: GraphletDesign,
    *,
    point: OperatingPoint,
    output_dir: Path,
    run_config_hash: str,
) -> dict[str, Any]:
    """Evaluate all treatments for one base design and checkpoint atomically."""

    checkpoint_path = output_dir / "checkpoints" / "designs" / f"{safe_name(design.graphlet_id)}.json"
    checkpoint_hash = json_hash(
        {
            "schema": 4,
            "config_hash": run_config_hash,
            "runner_source_hash": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "design": asdict(design),
            "point": point.as_mapping(),
        }
    )
    if checkpoint_path.exists():
        try:
            cached = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            if cached.get("checkpoint_hash") == checkpoint_hash:
                cached["resumed"] = True
                return cached
        except (OSError, json.JSONDecodeError):
            pass

    if not design.preflight_pass:
        payload = {
            "checkpoint_hash": checkpoint_hash,
            "graphlet_id": design.graphlet_id,
            "resumed": False,
            "component_rows": [],
            "descriptor_rows": [],
        }
        write_json_atomic(checkpoint_path, payload)
        return payload

    chunk_size = int(config.get("probe", {}).get("eval_chunk_size", 32768))
    tolerance = fixed_fitness_tolerance(config)
    scale = float(config.get("generator", {}).get("fitness_output_scale", 10.0))
    lambda_score = float(config.get("generator", {}).get("lambda_score", 1.0))
    sequences = materialize_design(design)
    dead_sequences = dead_control_sequences(config, design)
    plan = make_metric_plan(config, shape=design.shape, sampling_group_id=design.sampling_group_id)
    sequence_digest = array_hash(sequences)
    dead_digest = array_hash(dead_sequences)

    base_by_architecture: dict[str, tuple[OperatingPoint, dict[str, np.ndarray], dict[str, np.ndarray]]] = {}
    base_keys = (
        "effective_stability_margin_used",
        "stability_gate",
        "functional_capacity",
        "effective_functional_capacity",
    )
    for architecture in configured_architectures(config, design.arm):
        base_point = closed_control_point(config, point) if architecture == ARCH_CLOSED else point
        cache_key = ARCH_CLOSED if architecture == ARCH_CLOSED else "open"
        if cache_key in base_by_architecture:
            continue
        base_env = build_environment(
            config,
            seed=design.seed,
            operating_point=base_point.as_mapping(),
            routing_mode="structured",
            gnk_weight=0.0,
            architecture=ARCH_PRODUCT,
        )
        base_by_architecture[cache_key] = (
            base_point,
            evaluate_components(base_env, sequences, keys=base_keys, chunk_size=chunk_size),
            evaluate_components(base_env, dead_sequences, keys=base_keys, chunk_size=chunk_size),
        )

    texture_by_route: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for route in configured_routings(config, design.arm):
        texture_env = build_environment(
            config,
            seed=design.seed,
            operating_point=point.as_mapping(),
            routing_mode=route,
            gnk_weight=1.0,
            architecture=ARCH_PRODUCT,
        )
        texture_by_route[route] = (
            evaluate_components(
                texture_env,
                sequences,
                keys=("gnk_epistatic_score",),
                chunk_size=chunk_size,
            )["gnk_epistatic_score"],
            evaluate_components(
                texture_env,
                dead_sequences,
                keys=("gnk_epistatic_score",),
                chunk_size=chunk_size,
            )["gnk_epistatic_score"],
        )

    component_rows: list[dict[str, Any]] = []
    descriptor_rows: list[dict[str, Any]] = []
    for architecture in configured_architectures(config, design.arm):
        base_key = ARCH_CLOSED if architecture == ARCH_CLOSED else "open"
        architecture_point, base, dead_base = base_by_architecture[base_key]
        anchor = np.asarray(design.anchor_sequence, dtype=np.int64)
        anchor_index_candidates = np.flatnonzero(np.all(sequences == anchor[None, :], axis=1))
        anchor_index = int(anchor_index_candidates[0]) if anchor_index_candidates.size else 0
        for route in configured_routings(config, design.arm):
            texture, dead_texture = texture_by_route[route]
            for weight in configured_weights(config):
                observed = derive_observable(
                    architecture,
                    weight=float(weight),
                    base_margin=base["effective_stability_margin_used"],
                    base_gate=base["stability_gate"],
                    base_function_capacity=base["functional_capacity"],
                    raw_gnk_score=texture,
                    threshold_steepness=architecture_point.threshold_steepness,
                    fitness_output_scale=scale,
                    lambda_score=lambda_score,
                )
                dead_observed = derive_observable(
                    architecture,
                    weight=float(weight),
                    base_margin=dead_base["effective_stability_margin_used"],
                    base_gate=dead_base["stability_gate"],
                    base_function_capacity=dead_base["functional_capacity"],
                    raw_gnk_score=dead_texture,
                    threshold_steepness=architecture_point.threshold_steepness,
                    fitness_output_scale=scale,
                    lambda_score=lambda_score,
                )
                base_row = {
                    "stage": design.arm,
                    "operating_point_id": design.operating_point_id,
                    "seed": int(design.seed),
                    "dataset_id": design.dataset_id,
                    "graphlet_index": int(design.graphlet_index),
                    "graphlet_id": design.graphlet_id,
                    "sampling_group_id": design.sampling_group_id,
                    "architecture": architecture,
                    "support_mode": route,
                    "gnk_weight": float(weight),
                    "sequence_hash": sequence_digest,
                    "dead_control_hash": dead_digest,
                    "metric_plan_hash": plan.plan_hash,
                    "fixed_fitness_tolerance": tolerance,
                    "preflight_pass": bool(design.preflight_pass),
                }
                gate = observed["stability_gate"]
                viable = gate >= float(config.get("gates", {}).get("preflight", {}).get("viable_gate_threshold", 0.90))
                capacity = observed["effective_functional_capacity"]
                reconstruction_error = (
                    float(np.max(np.abs(observed["fitness"] - scale * gate * capacity)))
                    if architecture in {ARCH_PRODUCT, ARCH_CLOSED}
                    else math.nan
                )
                titration = _titration_summary(
                    config,
                    architecture=architecture,
                    weight=float(weight),
                    point=architecture_point,
                    capacity_at_background=float(base["functional_capacity"][anchor_index]),
                    raw_gnk_at_background=float(texture[anchor_index]),
                )
                component_row = base_row | {
                    "actual_anchor_margin": float(base["effective_stability_margin_used"][anchor_index]),
                    "gate_occupancy": float(np.mean(viable)),
                    "dead_gate_low_fraction": float(
                        np.mean(
                            dead_observed["stability_gate"]
                            <= float(config.get("gates", {}).get("preflight", {}).get("dead_gate_threshold", 0.10))
                        )
                    ),
                    "effective_functional_capacity_sd": (
                        float(np.std(capacity[viable])) if np.any(viable) else math.nan
                    ),
                    "effective_functional_capacity_interior_fraction": (
                        float(np.mean((capacity[viable] > 0.05) & (capacity[viable] < 0.95)))
                        if np.any(viable)
                        else math.nan
                    ),
                    "product_reconstruction_max_abs_error": reconstruction_error,
                    "raw_gnk_score_sd": float(np.std(texture)),
                    "dead_fitness_median": float(np.median(dead_observed["fitness"])),
                    "dead_fitness_q95": float(np.quantile(dead_observed["fitness"], 0.95)),
                    **titration,
                }
                if architecture == ARCH_CLOSED:
                    closed_target = float(architecture_point.stability_margin)
                    preflight_cfg = config.get("gates", {}).get("preflight", {})
                    closed_gate_threshold = float(
                        preflight_cfg.get(
                            "closed_control_gate_threshold",
                            preflight_cfg.get("dead_gate_threshold", 0.10),
                        )
                    )
                    closed_min_low_fraction = float(
                        preflight_cfg.get(
                            "closed_control_min_low_gate_fraction",
                            1.0 - float(preflight_cfg.get("closed_control_max_gate_occupancy", 0.10)),
                        )
                    )
                    closed_low_fraction = float(np.mean(gate <= closed_gate_threshold))
                    component_row.update(
                        {
                            "closed_control_target_effective_margin": closed_target,
                            "closed_control_gate_threshold": closed_gate_threshold,
                            "closed_control_min_low_gate_fraction": closed_min_low_fraction,
                            "closed_control_margin_pass": bool(
                                abs(float(component_row["actual_anchor_margin"]) - closed_target) <= 1e-8
                            ),
                            "closed_control_gate_low_fraction": closed_low_fraction,
                            "closed_control_occupancy_pass": bool(
                                closed_low_fraction >= closed_min_low_fraction
                            ),
                            "closed_control_preflight_pass": bool(
                                abs(float(component_row["actual_anchor_margin"]) - closed_target) <= 1e-8
                                and closed_low_fraction >= closed_min_low_fraction
                            ),
                        }
                    )
                metrics = metric_bundle(
                    observed["fitness"],
                    design.shape,
                    gate,
                    plan,
                    tolerance=tolerance,
                    exact_paths=design.arm == "calibration",
                )
                descriptor_rows.append(
                    base_row
                    | metrics
                    | {
                        "dead_fitness_median": component_row["dead_fitness_median"],
                        "dead_fitness_q95": component_row["dead_fitness_q95"],
                        **titration,
                    }
                )
                component_rows.append(component_row)

    payload = {
        "checkpoint_hash": checkpoint_hash,
        "graphlet_id": design.graphlet_id,
        "resumed": False,
        "component_rows": component_rows,
        "descriptor_rows": descriptor_rows,
    }
    write_json_atomic(checkpoint_path, payload)
    return payload


def _finite_median(values: Iterable[Any]) -> float:
    array = np.asarray(list(values), dtype=np.float64)
    array = array[np.isfinite(array)]
    return float(np.median(array)) if array.size else math.nan


def _all_finite_true(values: Iterable[Any]) -> bool:
    items = list(values)
    return bool(items) and all(bool(value) for value in items)


def build_viability_checks(config: dict[str, Any], descriptor_rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    gate_cfg = config.get("gates", {}).get("g1_viability", {})
    scale = float(config.get("generator", {}).get("fitness_output_scale", 10.0))
    range_tolerance = float(gate_cfg.get("range_tolerance", 1e-8))
    observable_min = float(gate_cfg.get("observable_min", 0.0))
    observable_max = scale * float(gate_cfg.get("observable_max_fraction_of_scale", 1.0))
    dead_q95_max = scale * float(gate_cfg.get("dead_fitness_95th_max_fraction_of_scale", 0.10))
    viable_median_min = scale * float(gate_cfg.get("viable_median_min_fraction_of_scale", 0.30))
    viable_sd_min = scale * float(gate_cfg.get("viable_sd_min_fraction_of_scale", 0.05))
    rows: list[dict[str, Any]] = []
    for metric in descriptor_rows:
        range_pass = bool(
            np.isfinite(float(metric.get("fitness_min", math.nan)))
            and np.isfinite(float(metric.get("fitness_max", math.nan)))
            and float(metric["fitness_min"]) >= observable_min - range_tolerance
            and float(metric["fitness_max"]) <= observable_max + range_tolerance
        )
        dead_pass = bool(
            np.isfinite(float(metric.get("dead_fitness_q95", math.nan)))
            and float(metric["dead_fitness_q95"]) <= dead_q95_max
        )
        viable_pass = bool(
            np.isfinite(float(metric.get("viable_fitness_median", math.nan)))
            and np.isfinite(float(metric.get("viable_fitness_sd", math.nan)))
            and float(metric["viable_fitness_median"]) >= viable_median_min
            and float(metric["viable_fitness_sd"]) >= viable_sd_min
        )
        titration_pass = bool(metric.get("titration_pass", False))
        rows.append(
            {key: metric[key] for key in (
                "stage",
                "operating_point_id",
                "seed",
                "dataset_id",
                "graphlet_index",
                "graphlet_id",
                "architecture",
                "support_mode",
                "gnk_weight",
            )}
            | {
                "observable_min": float(metric.get("fitness_min", math.nan)),
                "observable_max": float(metric.get("fitness_max", math.nan)),
                "dead_fitness_q95": float(metric.get("dead_fitness_q95", math.nan)),
                "viable_fitness_median": float(metric.get("viable_fitness_median", math.nan)),
                "viable_fitness_sd": float(metric.get("viable_fitness_sd", math.nan)),
                "titration_transition_slope": float(metric.get("titration_transition_slope", math.nan)),
                "titration_saturated_slope": float(metric.get("titration_saturated_slope", math.nan)),
                "titration_saturated_to_transition_ratio": float(
                    metric.get("titration_saturated_to_transition_ratio", math.nan)
                ),
                "observable_range_pass": range_pass,
                "dead_stays_dead_pass": dead_pass,
                "viable_expressed_pass": viable_pass,
                "titration_implementation_pass": titration_pass,
                "g1_pass": bool(range_pass and dead_pass and viable_pass and titration_pass),
            }
        )
    return rows


TRANSMISSION_METRICS = ("nonadditive_variance_fraction", "sign_epistasis_fraction")


def build_transmission_checks(
    config: dict[str, Any],
    descriptor_rows: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    working = working_weight(config)
    frame = pd.DataFrame(list(descriptor_rows))
    if frame.empty:
        return []
    identity = [
        "stage",
        "operating_point_id",
        "seed",
        "dataset_id",
        "graphlet_index",
        "graphlet_id",
        "architecture",
        "support_mode",
    ]
    graphlet_rows: list[dict[str, Any]] = []
    for keys, group in frame.groupby(identity, dropna=False, sort=False):
        base = dict(zip(identity, keys))
        zero = group[np.isclose(group["gnk_weight"].astype(float), 0.0)]
        endpoint = group[np.isclose(group["gnk_weight"].astype(float), working)]
        for metric in TRANSMISSION_METRICS:
            value0 = _finite_median(zero[metric]) if metric in zero else math.nan
            value1 = _finite_median(endpoint[metric]) if metric in endpoint else math.nan
            graphlet_rows.append(
                base
                | {
                    "aggregation_level": "graphlet",
                    "metric": metric,
                    "weight_zero_value": value0,
                    "working_weight_value": value1,
                    "signed_delta": float(value1 - value0) if np.isfinite(value0) and np.isfinite(value1) else math.nan,
                }
            )
    graphlet_frame = pd.DataFrame(graphlet_rows)
    template_identity = [
        "stage",
        "operating_point_id",
        "seed",
        "dataset_id",
        "architecture",
        "support_mode",
        "metric",
    ]
    template_rows: list[dict[str, Any]] = []
    if not graphlet_frame.empty:
        for keys, group in graphlet_frame.groupby(template_identity, dropna=False, sort=False):
            template_rows.append(
                dict(zip(template_identity, keys))
                | {
                    "aggregation_level": "template_seed",
                    "weight_zero_value": _finite_median(group["weight_zero_value"]),
                    "working_weight_value": _finite_median(group["working_weight_value"]),
                    "signed_delta": _finite_median(group["signed_delta"]),
                    "graphlet_count": int(len(group)),
                }
            )
    return graphlet_rows + template_rows


def exact_one_sided_sign_p(values: Sequence[float]) -> float:
    """Exact P(Binomial(n,.5)>=positive); zero/tied effects count nonpositive."""

    array = np.asarray(values, dtype=np.float64)
    if array.size == 0 or np.any(~np.isfinite(array)):
        return math.nan
    positives = int(np.sum(array > 0.0))
    n = int(array.size)
    return float(sum(math.comb(n, count) for count in range(positives, n + 1)) / (2**n))


def fixed_seed_bootstrap_lower(
    values: Sequence[float],
    *,
    replicates: int,
    confidence: float,
    seed_key: str,
    base_seed: int = 59059,
) -> float:
    array = np.asarray(values, dtype=np.float64)
    if array.size == 0 or np.any(~np.isfinite(array)):
        return math.nan
    rng = np.random.default_rng(stable_uint64("seed-bootstrap", int(base_seed), seed_key) % (2**63 - 1))
    draws = rng.integers(0, array.size, size=(int(replicates), array.size), dtype=np.int64)
    statistics = np.median(array[draws], axis=1)
    # One-sided lower confidence bound: 5th percentile for a 95% interval.
    return float(np.quantile(statistics, 1.0 - float(confidence)))


def transmission_gate(
    config: dict[str, Any],
    checks: Sequence[dict[str, Any]],
    *,
    stage: str,
    architecture: str,
    support_mode: str,
    expected_seeds: Sequence[int],
    expected_templates: Sequence[str],
) -> tuple[bool, list[dict[str, Any]]]:
    cfg = config.get("gates", {}).get("g2_transmission", {})
    frame = pd.DataFrame(
        [
            row
            for row in checks
            if row.get("aggregation_level") == "template_seed"
            and row.get("stage") == stage
            and row.get("architecture") == architecture
            and row.get("support_mode") == support_mode
        ]
    )
    details: list[dict[str, Any]] = []
    overall_pass = True
    floors = {
        "nonadditive_variance_fraction": float(
            cfg.get("metrics", {}).get("nonadditive_variance_fraction", {}).get("min_median_delta", 0.02)
        ),
        "sign_epistasis_fraction": float(
            cfg.get("metrics", {}).get("sign_epistasis_fraction", {}).get("min_median_delta", 0.02)
        ),
    }
    for metric in TRANSMISSION_METRICS:
        subset = frame[frame["metric"] == metric] if not frame.empty else pd.DataFrame()
        if subset.empty:
            result = {"metric": metric, "g2_metric_pass": False, "reason": "missing_metric_rows"}
            details.append(result)
            overall_pass = False
            continue
        seed_values = (
            subset.groupby("seed", sort=False)["signed_delta"].median().reindex(list(expected_seeds))
        )
        template_values = (
            subset.groupby("dataset_id", sort=False)["signed_delta"].median().reindex(list(expected_templates))
        )
        seed_array = seed_values.to_numpy(dtype=np.float64)
        template_array = template_values.to_numpy(dtype=np.float64)
        median_delta = _finite_median(seed_array)
        positive_seeds = int(np.sum(seed_array > 0.0))
        finite_complete = bool(np.all(np.isfinite(seed_array)) and np.all(np.isfinite(template_array)))
        if stage == "calibration":
            required_positive = len(expected_seeds)
            sign_p = math.nan
            bootstrap_lower = math.nan
            metric_pass = bool(
                finite_complete
                and median_delta >= floors[metric]
                and positive_seeds == required_positive
            )
        else:
            required_positive = min(
                int(cfg.get("confirmation_min_positive_seeds", 7)), len(expected_seeds)
            )
            sign_p = exact_one_sided_sign_p(seed_array)
            bootstrap_lower = fixed_seed_bootstrap_lower(
                seed_array,
                replicates=int(cfg.get("bootstrap_replicates", 10000)),
                confidence=float(cfg.get("bootstrap_ci", 0.95)),
                seed_key=f"g2|{stage}|{architecture}|{support_mode}|{metric}",
                base_seed=int(cfg.get("bootstrap_seed", 59059)),
            )
            metric_pass = bool(
                finite_complete
                and median_delta >= floors[metric]
                and np.all(template_array > 0.0)
                and positive_seeds >= required_positive
                and np.isfinite(sign_p)
                and sign_p <= float(cfg.get("sign_test_max_p", 0.0352))
                and np.isfinite(bootstrap_lower)
                and bootstrap_lower > 0.0
            )
        details.append(
            {
                "metric": metric,
                "overall_median_delta": median_delta,
                "absolute_floor": floors[metric],
                "positive_seed_count": positive_seeds,
                "required_positive_seed_count": required_positive,
                "every_template_positive": bool(
                    template_array.size and np.all(np.isfinite(template_array)) and np.all(template_array > 0.0)
                ),
                "one_sided_sign_p": sign_p,
                "bootstrap_one_sided_lower": bootstrap_lower,
                "g2_metric_pass": metric_pass,
            }
        )
        overall_pass = bool(overall_pass and metric_pass)
    return overall_pass, details


TOPOLOGY_RHO_METRICS = (
    "nonadditive_variance_fraction",
    "sign_epistasis_fraction",
    "local_maxima_density",
)


def build_topography_checks(
    config: dict[str, Any],
    descriptor_rows: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    frame = pd.DataFrame(list(descriptor_rows))
    if frame.empty:
        return []
    working = working_weight(config)
    identity = [
        "stage",
        "operating_point_id",
        "seed",
        "dataset_id",
        "architecture",
        "support_mode",
    ]
    rows: list[dict[str, Any]] = []
    for keys, group in frame.groupby(identity, dropna=False, sort=False):
        base = dict(zip(identity, keys))
        curves = group.groupby("gnk_weight", sort=True)[
            [*TOPOLOGY_RHO_METRICS, "local_maxima_count", "basin_convergence_fraction"]
        ].median(numeric_only=True)
        weights = [float(value) for value in curves.index]

        def curve_value(metric: str, weight: float) -> float:
            matches = [index for index in curves.index if np.isclose(float(index), float(weight))]
            if not matches:
                return math.nan
            return float(curves.loc[matches[0], metric])

        density0 = curve_value("local_maxima_density", 0.0)
        density1 = curve_value("local_maxima_density", working)
        count0 = curve_value("local_maxima_count", 0.0)
        count1 = curve_value("local_maxima_count", working)
        basin0 = curve_value("basin_convergence_fraction", 0.0)
        basin1 = curve_value("basin_convergence_fraction", working)
        row = base | {
            "aggregation_level": "template_seed",
            "weight_count": int(len(curves)),
            "local_maxima_density_w0": density0,
            "local_maxima_density_working": density1,
            "local_maxima_density_delta": (
                float(density1 - density0) if np.isfinite(density0) and np.isfinite(density1) else math.nan
            ),
            "local_maxima_count_w0": count0,
            "local_maxima_count_working": count1,
            "local_maxima_count_delta": (
                float(count1 - count0) if np.isfinite(count0) and np.isfinite(count1) else math.nan
            ),
            "basin_convergence_w0": basin0,
            "basin_convergence_working": basin1,
            "basin_convergence_delta": (
                float(basin1 - basin0) if np.isfinite(basin0) and np.isfinite(basin1) else math.nan
            ),
        }
        for metric in TOPOLOGY_RHO_METRICS:
            row[f"rho_{metric}"] = spearman(weights, [float(value) for value in curves[metric]])
        rows.append(row)
    return rows


def topography_gate(
    config: dict[str, Any],
    checks: Sequence[dict[str, Any]],
    *,
    stage: str,
    architecture: str,
    support_mode: str,
    expected_seeds: Sequence[int],
    expected_templates: Sequence[str],
) -> tuple[bool, dict[str, Any]]:
    cfg = config.get("gates", {}).get("g3_topography", {})
    rows = [
        row
        for row in checks
        if row.get("stage") == stage
        and row.get("architecture") == architecture
        and row.get("support_mode") == support_mode
    ]
    frame = pd.DataFrame(rows)
    if frame.empty:
        return False, {"g3_pass": False, "reason": "missing_topography_rows"}
    seed_frame = frame.groupby("seed", sort=False).median(numeric_only=True).reindex(list(expected_seeds))
    template_frame = frame.groupby("dataset_id", sort=False).median(numeric_only=True).reindex(list(expected_templates))
    complete = bool(
        len(frame) >= len(expected_seeds) * len(expected_templates)
        and set(int(value) for value in frame["seed"]) >= set(int(value) for value in expected_seeds)
        and set(str(value) for value in frame["dataset_id"]) >= set(str(value) for value in expected_templates)
    )
    median_rho_floor = float(cfg.get("min_overall_median_rho", 0.80))
    template_rho_floor = float(cfg.get("min_template_median_rho", 0.50))
    required_correct = min(int(cfg.get("confirmation_min_correct_direction_seeds", 7)), len(expected_seeds))
    rho_details: dict[str, Any] = {}
    rho_pass = True
    for metric in TOPOLOGY_RHO_METRICS:
        column = f"rho_{metric}"
        seed_values = seed_frame[column].to_numpy(dtype=np.float64)
        template_values = template_frame[column].to_numpy(dtype=np.float64)
        median_rho = _finite_median(seed_values)
        correct = int(np.sum(seed_values > 0.0))
        if stage == "calibration":
            metric_pass = bool(
                np.all(np.isfinite(seed_values))
                and correct == len(expected_seeds)
                and median_rho >= median_rho_floor
            )
        else:
            metric_pass = bool(
                np.all(np.isfinite(seed_values))
                and np.all(np.isfinite(template_values))
                and correct >= required_correct
                and median_rho >= median_rho_floor
                and np.all(template_values >= template_rho_floor)
            )
        rho_details[f"{metric}_median_rho"] = median_rho
        rho_details[f"{metric}_correct_direction_seed_count"] = correct
        rho_details[f"{metric}_rho_pass"] = metric_pass
        rho_pass = bool(rho_pass and metric_pass)

    density_floor = float(cfg.get("local_maxima_density_min_delta", 0.005))
    count_floor = float(cfg.get("local_maxima_count_min_increase", 1.0))
    density_seed = seed_frame["local_maxima_density_delta"].to_numpy(dtype=np.float64)
    count_seed = seed_frame["local_maxima_count_delta"].to_numpy(dtype=np.float64)
    density_template = template_frame["local_maxima_density_delta"].to_numpy(dtype=np.float64)
    count_template = template_frame["local_maxima_count_delta"].to_numpy(dtype=np.float64)
    endpoint_seed_pass = np.isfinite(density_seed) & np.isfinite(count_seed) & (density_seed >= density_floor) & (count_seed >= count_floor)
    if stage == "calibration":
        endpoint_pass = bool(
            np.all(np.isfinite(density_seed))
            and np.all(np.isfinite(count_seed))
            and _finite_median(density_seed) >= density_floor
            and _finite_median(count_seed) >= count_floor
            and np.all(density_seed > 0.0)
            and np.all(count_seed > 0.0)
        )
        endpoint_required = len(expected_seeds)
    else:
        endpoint_required = required_correct
        endpoint_pass = bool(
            np.all(np.isfinite(density_seed))
            and np.all(np.isfinite(count_seed))
            and _finite_median(density_seed) >= density_floor
            and _finite_median(count_seed) >= count_floor
            and int(np.sum(endpoint_seed_pass)) >= endpoint_required
            and np.all(np.isfinite(density_template))
            and np.all(np.isfinite(count_template))
            and np.all(density_template > 0.0)
            and np.all(count_template > 0.0)
        )

    basin_min = float(cfg.get("working_weight_basin_min", 0.05))
    basin_max = float(cfg.get("working_weight_basin_max", 0.80))
    basin_delta_max = float(cfg.get("basin_max_delta_from_zero", -0.05))
    basin_seed_value = seed_frame["basin_convergence_working"].to_numpy(dtype=np.float64)
    basin_seed_delta = seed_frame["basin_convergence_delta"].to_numpy(dtype=np.float64)
    basin_template_value = template_frame["basin_convergence_working"].to_numpy(dtype=np.float64)
    basin_template_delta = template_frame["basin_convergence_delta"].to_numpy(dtype=np.float64)
    basin_seed_pass = (
        np.isfinite(basin_seed_value)
        & np.isfinite(basin_seed_delta)
        & (basin_seed_value >= basin_min)
        & (basin_seed_value <= basin_max)
        & (basin_seed_delta <= basin_delta_max)
    )
    if stage == "calibration":
        basin_pass = bool(np.all(basin_seed_pass))
        basin_required = len(expected_seeds)
    else:
        basin_required = required_correct
        basin_pass = bool(
            int(np.sum(basin_seed_pass)) >= basin_required
            and np.all(np.isfinite(basin_template_value))
            and np.all(np.isfinite(basin_template_delta))
            and np.all((basin_template_value >= basin_min) & (basin_template_value <= basin_max))
            and np.all(basin_template_delta <= basin_delta_max)
        )
    details = {
        **rho_details,
        "endpoint_density_median_delta": _finite_median(density_seed),
        "endpoint_count_median_delta": _finite_median(count_seed),
        "endpoint_seed_pass_count": int(np.sum(endpoint_seed_pass)),
        "endpoint_required_seed_count": endpoint_required,
        "endpoint_pass": endpoint_pass,
        "basin_working_median": _finite_median(basin_seed_value),
        "basin_delta_median": _finite_median(basin_seed_delta),
        "basin_seed_pass_count": int(np.sum(basin_seed_pass)),
        "basin_required_seed_count": basin_required,
        "basin_pass": basin_pass,
        "complete_expected_units": complete,
    }
    passed = bool(complete and rho_pass and endpoint_pass and basin_pass)
    details["g3_pass"] = passed
    return passed, details


def paired_seed_bootstrap_interval(
    values: Sequence[float],
    *,
    interval: float,
    replicates: int,
    seed_key: str,
    base_seed: int = 59060,
) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    if array.size == 0 or np.any(~np.isfinite(array)):
        return math.nan, math.nan
    rng = np.random.default_rng(stable_uint64("g4-bootstrap", int(base_seed), seed_key) % (2**63 - 1))
    draws = rng.integers(0, array.size, size=(int(replicates), array.size), dtype=np.int64)
    statistics = np.median(array[draws], axis=1)
    alpha = 0.5 * (1.0 - float(interval))
    return float(np.quantile(statistics, alpha)), float(np.quantile(statistics, 1.0 - alpha))


def build_routing_equivalence(
    config: dict[str, Any],
    descriptor_rows: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    cfg = config.get("gates", {}).get("g4_routing_equivalence", {})
    working = working_weight(config)
    interval = float(cfg.get("interval", 0.90))
    replicates = int(cfg.get("bootstrap_replicates", 10000))
    margins = cfg.get("margins", {})
    metrics = (
        "nonadditive_variance_fraction",
        "sign_epistasis_fraction",
        "local_maxima_density",
        "basin_convergence_fraction",
    )
    frame = pd.DataFrame(list(descriptor_rows))
    if frame.empty:
        return []
    frame = frame[np.isclose(frame["gnk_weight"].astype(float), working)]
    rows: list[dict[str, Any]] = []
    for (stage, point_id, architecture), arch_frame in frame.groupby(
        ["stage", "operating_point_id", "architecture"], sort=False
    ):
        structured = arch_frame[arch_frame["support_mode"] == "structured"]
        if structured.empty:
            continue
        for control in ("shuffled", "random"):
            control_frame = arch_frame[arch_frame["support_mode"] == control]
            if control_frame.empty:
                continue
            join_keys = ["seed", "dataset_id", "graphlet_index"]
            for metric in metrics:
                paired = structured[join_keys + [metric]].merge(
                    control_frame[join_keys + [metric]],
                    on=join_keys,
                    how="inner",
                    suffixes=("_structured", "_control"),
                )
                if paired.empty:
                    continue
                paired["difference"] = paired[f"{metric}_structured"] - paired[f"{metric}_control"]
                seed_values = paired.groupby("seed", sort=False)["difference"].median().to_numpy(dtype=np.float64)
                lower, upper = paired_seed_bootstrap_interval(
                    seed_values,
                    interval=interval,
                    replicates=replicates,
                    seed_key=f"{stage}|{point_id}|{architecture}|{control}|{metric}",
                    base_seed=int(cfg.get("bootstrap_seed", 59060)),
                )
                margin = float(margins.get(metric, 0.0))
                equivalent = bool(
                    np.isfinite(lower)
                    and np.isfinite(upper)
                    and lower >= -margin
                    and upper <= margin
                )
                rows.append(
                    {
                        "stage": stage,
                        "operating_point_id": point_id,
                        "architecture": architecture,
                        "structured_vs": control,
                        "metric": metric,
                        "paired_seed_count": int(seed_values.size),
                        "paired_seed_median_difference": _finite_median(seed_values),
                        "bootstrap_interval": interval,
                        "bootstrap_lower": lower,
                        "bootstrap_upper": upper,
                        "equivalence_margin": margin,
                        "equivalent": equivalent,
                        "diagnostic_only": True,
                    }
                )
    return rows


def viability_gate(
    checks: Sequence[dict[str, Any]],
    *,
    stage: str,
    architecture: str,
    support_mode: str,
    expected_designs: Sequence[GraphletDesign],
    weights: Sequence[float],
) -> tuple[bool, dict[str, Any]]:
    rows = [
        row
        for row in checks
        if row.get("stage") == stage
        and row.get("architecture") == architecture
        and row.get("support_mode") == support_mode
    ]
    observed = {
        (str(row["graphlet_id"]), float(row["gnk_weight"]))
        for row in rows
        if bool(row.get("g1_pass", False))
    }
    expected = {
        (str(design.graphlet_id), float(weight))
        for design in expected_designs
        for weight in weights
    }
    passed = bool(expected and observed >= expected)
    return passed, {
        "expected_treatment_count": int(len(expected)),
        "passing_treatment_count": int(len(observed.intersection(expected))),
        "g1_pass": passed,
    }


def confirmation_seed_summary(
    transmission_checks: Sequence[dict[str, Any]],
    topography_checks: Sequence[dict[str, Any]],
    *,
    architecture: str,
    support_mode: str,
) -> list[dict[str, Any]]:
    transmission = pd.DataFrame(
        [
            row
            for row in transmission_checks
            if row.get("stage") == "confirmation"
            and row.get("aggregation_level") == "template_seed"
            and row.get("architecture") == architecture
            and row.get("support_mode") == support_mode
        ]
    )
    topography = pd.DataFrame(
        [
            row
            for row in topography_checks
            if row.get("stage") == "confirmation"
            and row.get("architecture") == architecture
            and row.get("support_mode") == support_mode
        ]
    )
    seeds = sorted(
        set(int(value) for value in transmission.get("seed", pd.Series(dtype=int)))
        | set(int(value) for value in topography.get("seed", pd.Series(dtype=int)))
    )
    rows: list[dict[str, Any]] = []
    for seed in seeds:
        row: dict[str, Any] = {
            "stage": "confirmation",
            "seed": seed,
            "architecture": architecture,
            "support_mode": support_mode,
        }
        if not transmission.empty:
            subset = transmission[transmission["seed"] == seed]
            for metric in TRANSMISSION_METRICS:
                row[f"delta_{metric}"] = _finite_median(
                    subset.loc[subset["metric"] == metric, "signed_delta"]
                )
        if not topography.empty:
            subset = topography[topography["seed"] == seed]
            for metric in TOPOLOGY_RHO_METRICS:
                row[f"rho_{metric}"] = _finite_median(subset[f"rho_{metric}"])
            for column in (
                "local_maxima_density_delta",
                "local_maxima_count_delta",
                "basin_convergence_working",
                "basin_convergence_delta",
            ):
                row[column] = _finite_median(subset[column])
        rows.append(row)
    return rows


def design_catalog_rows(config: dict[str, Any], designs: Sequence[GraphletDesign]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    stability_hashes: dict[int, str] = {}
    for design in designs:
        if design.seed not in stability_hashes:
            basis_env = build_environment(
                config,
                seed=design.seed,
                operating_point={"stability_margin": 0.0, "threshold_steepness": 1.0},
                routing_mode="structured",
                gnk_weight=0.0,
                architecture=ARCH_PRODUCT,
            )
            stability_hashes[design.seed] = array_hash(np.asarray(basis_env.stability_matrix))
        rows.append(
            asdict(design)
            | {
                "anchor_sequence": ";".join(str(value) for value in design.anchor_sequence),
                "sites": ";".join(str(value) for value in design.sites),
                "choices": "|".join(",".join(str(value) for value in axis) for axis in design.choices),
                "shape": ";".join(str(value) for value in design.shape),
                "sequence_hash": array_hash(materialize_design(design)),
                "dead_control_hash": array_hash(dead_control_sequences(config, design)),
                "dead_control_count": int(config.get("probe", {}).get("dead_sequence_samples", 4096)),
                "allele_selection": str(
                    config.get("probe", {}).get("allele_selection", "hashed_including_anchor")
                ),
                "stability_matrix_hash": stability_hashes[design.seed],
            }
        )
    return rows


def expected_template_ids(config: dict[str, Any], stage: str, designs: Sequence[GraphletDesign]) -> list[str]:
    if stage == "calibration":
        return sorted({design.dataset_id for design in designs})
    configured = [str(value) for value in config.get("data", {}).get("dataset_ids", [])]
    return configured or sorted({design.dataset_id for design in designs})


def stage_decisions(
    config: dict[str, Any],
    *,
    stage: str,
    designs: Sequence[GraphletDesign],
    viability_rows: Sequence[dict[str, Any]],
    transmission_rows: Sequence[dict[str, Any]],
    topography_rows: Sequence[dict[str, Any]],
    component_rows: Sequence[dict[str, Any]] = (),
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    seeds = configured_seeds(config, stage)
    templates = expected_template_ids(config, stage, designs)
    architectures = [
        architecture
        for architecture in configured_architectures(config, stage)
        if architecture != ARCH_CLOSED
    ]
    rows: list[dict[str, Any]] = []
    details: dict[str, dict[str, Any]] = {}
    closed_required = ARCH_CLOSED in configured_architectures(config, stage)
    closed_rows = [
        row
        for row in component_rows
        if row.get("stage") == stage and row.get("architecture") == ARCH_CLOSED
    ]
    closed_control_pass = bool(
        not closed_required
        or (
            closed_rows
            and all(bool(row.get("closed_control_preflight_pass", False)) for row in closed_rows)
        )
    )
    for architecture in architectures:
        g1_pass, g1_detail = viability_gate(
            viability_rows,
            stage=stage,
            architecture=architecture,
            support_mode="structured",
            expected_designs=designs,
            weights=configured_weights(config),
        )
        g2_pass, g2_detail = transmission_gate(
            config,
            transmission_rows,
            stage=stage,
            architecture=architecture,
            support_mode="structured",
            expected_seeds=seeds,
            expected_templates=templates,
        )
        g3_pass, g3_detail = topography_gate(
            config,
            topography_rows,
            stage=stage,
            architecture=architecture,
            support_mode="structured",
            expected_seeds=seeds,
            expected_templates=templates,
        )
        preflight_pass = bool(
            designs and all(design.preflight_pass for design in designs) and closed_control_pass
        )
        mechanism_pass = bool(preflight_pass and g1_pass and g2_pass and g3_pass)
        status = (
            "LEGACY_TWO_CHANNEL_MECHANISM_DEMONSTRATED"
            if architecture == ARCH_PRODUCT and mechanism_pass
            else "MECHANISM_CHECKS_PASS"
            if mechanism_pass
            else "MECHANISM_FAIL"
        )
        rows.append(
            {
                "stage": stage,
                "architecture": architecture,
                "primary_architecture": architecture == ARCH_PRODUCT,
                "support_mode": "structured",
                "preflight_pass": preflight_pass,
                "closed_control_preflight_pass": closed_control_pass,
                "g1_pass": g1_pass,
                "g2_pass": g2_pass,
                "g3_pass": g3_pass,
                "mechanism_pass": mechanism_pass,
                "decision": status,
                "scope": "internally_controllable_candidate_mechanism_only",
            }
        )
        details[architecture] = {
            "closed_control_preflight_pass": closed_control_pass,
            "g1": g1_detail,
            "g2": g2_detail,
            "g3": g3_detail,
        }
    return rows, details


def _deduplicate_rows(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    for row in rows:
        digest = json_hash(row)
        if digest not in seen:
            seen.add(digest)
            out.append(row)
    return out


def _primary_stage_pass(rows: Sequence[dict[str, Any]], stage: str) -> bool:
    matches = [
        row
        for row in rows
        if row.get("stage") == stage and row.get("architecture") == ARCH_PRODUCT
    ]
    return bool(matches and bool(matches[0].get("mechanism_pass", False)))


def resolve_confirmation_point(
    config: dict[str, Any],
    *,
    project_root: Path,
    output_dir: Path,
    selected_calibration_point: OperatingPoint | None,
) -> OperatingPoint:
    if selected_calibration_point is not None:
        return selected_calibration_point
    section = arm_config(config, "confirmation")
    frozen = section.get("frozen_operating_point")
    if isinstance(frozen, dict):
        if isinstance(frozen.get("selected_operating_point"), dict):
            frozen = frozen["selected_operating_point"]
        return _point_from_mapping(dict(frozen), 0)
    frozen_path = section.get("frozen_operating_point_path")
    candidates: list[Path] = []
    if frozen_path:
        path = Path(str(frozen_path))
        candidates.append(path if path.is_absolute() else project_root / path)
    candidates.append(output_dir / "frozen_operating_point.json")
    for path in candidates:
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        row = payload.get("selected_operating_point", payload)
        if isinstance(row, dict):
            return _point_from_mapping(dict(row), 0)
    raise ValueError("confirmation requires a frozen point from calibration")


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config).resolve(), quick=bool(args.quick))
    stages = configured_stage_order(config, args.arm)
    for stage in stages:
        validate_global_config(config, stage, quick=bool(args.quick))
    run_config_hash = json_hash(config)

    all_designs: list[GraphletDesign] = []
    all_support_rows: list[dict[str, Any]] = []
    component_rows: list[dict[str, Any]] = []
    descriptor_rows: list[dict[str, Any]] = []
    calibration_rows: list[dict[str, Any]] = []
    mechanism_rows: list[dict[str, Any]] = []
    decision_details: dict[str, Any] = {}
    preflight_points: list[dict[str, Any]] = []
    resumed_checkpoints = 0
    computed_checkpoints = 0
    selected_point: OperatingPoint | None = None
    calibration_mechanism_pass: bool | None = None
    confirmation_skipped_reason = ""

    if "calibration" in stages:
        points = operating_points(config, "calibration", project_root)
        seeds = configured_seeds(config, "calibration")
        designs_by_point: dict[str, list[GraphletDesign]] = {}
        for point in points:
            print(f"[EXP-059] calibration preflight {point.operating_point_id}", flush=True)
            designs, support_rows = prepare_designs_for_operating_point(
                config,
                arm="calibration",
                point=point,
                seeds=seeds,
                project_root=project_root,
            )
            checked = component_preflight(config, designs, point=point)
            designs_by_point[point.operating_point_id] = checked
            all_designs.extend(checked)
            all_support_rows.extend(support_rows)
            eligible = bool(checked and all(design.preflight_pass for design in checked))
            calibration_rows.append(
                {
                    "operating_point_id": point.operating_point_id,
                    "stability_margin": point.stability_margin,
                    "threshold_steepness": point.threshold_steepness,
                    "order_index": len(calibration_rows),
                    "design_count": len(checked),
                    "preflight_pass_count": sum(int(design.preflight_pass) for design in checked),
                    "component_eligible": eligible,
                    "selected": False,
                    "selection_uses_topology": False,
                }
            )
            preflight_points.append(
                {
                    "stage": "calibration",
                    "operating_point": point.as_mapping(),
                    "component_eligible": eligible,
                    "designs": [asdict(design) for design in checked],
                }
            )
        selected_point = select_operating_point(points, designs_by_point)
        if selected_point is not None:
            for row in calibration_rows:
                row["selected"] = row["operating_point_id"] == selected_point.operating_point_id
            selected_designs = designs_by_point[selected_point.operating_point_id]
            print(f"[EXP-059] selected {selected_point.operating_point_id} using component preflight only", flush=True)
            for index, design in enumerate(selected_designs, start=1):
                print(f"[EXP-059] calibration graphlet {index}/{len(selected_designs)}: {design.graphlet_id}", flush=True)
                payload = evaluate_design(
                    config,
                    design,
                    point=selected_point,
                    output_dir=output_dir,
                    run_config_hash=run_config_hash,
                )
                resumed_checkpoints += int(bool(payload.get("resumed", False)))
                computed_checkpoints += int(not bool(payload.get("resumed", False)))
                component_rows.extend(payload.get("component_rows", []))
                descriptor_rows.extend(payload.get("descriptor_rows", []))
            calibration_viability = build_viability_checks(
                config, [row for row in descriptor_rows if row.get("stage") == "calibration"]
            )
            calibration_transmission = build_transmission_checks(
                config, [row for row in descriptor_rows if row.get("stage") == "calibration"]
            )
            calibration_topography = build_topography_checks(
                config, [row for row in descriptor_rows if row.get("stage") == "calibration"]
            )
            rows, details = stage_decisions(
                config,
                stage="calibration",
                designs=selected_designs,
                viability_rows=calibration_viability,
                transmission_rows=calibration_transmission,
                topography_rows=calibration_topography,
                component_rows=[row for row in component_rows if row.get("stage") == "calibration"],
            )
            mechanism_rows.extend(rows)
            decision_details["calibration"] = details
            calibration_mechanism_pass = _primary_stage_pass(rows, "calibration")
            for row in calibration_rows:
                if row["selected"]:
                    row["g1_g2_g3_pass"] = calibration_mechanism_pass
        else:
            calibration_mechanism_pass = False
            confirmation_skipped_reason = "NO_COMPONENT_QUALIFIED_POINT"

    should_run_confirmation = "confirmation" in stages
    if "calibration" in stages and not bool(calibration_mechanism_pass):
        should_run_confirmation = False
        if not confirmation_skipped_reason:
            confirmation_skipped_reason = "SELECTED_COMPONENT_QUALIFIED_POINT_FAIL"

    confirmation_designs: list[GraphletDesign] = []
    if should_run_confirmation:
        confirmation_point = resolve_confirmation_point(
            config,
            project_root=project_root,
            output_dir=output_dir,
            selected_calibration_point=selected_point,
        )
        if selected_point is None:
            selected_point = confirmation_point
        seeds = configured_seeds(config, "confirmation")
        print(f"[EXP-059] confirmation frozen at {confirmation_point.operating_point_id}", flush=True)
        confirmation_designs, support_rows = prepare_designs_for_operating_point(
            config,
            arm="confirmation",
            point=confirmation_point,
            seeds=seeds,
            project_root=project_root,
        )
        confirmation_designs = component_preflight(config, confirmation_designs, point=confirmation_point)
        all_designs.extend(confirmation_designs)
        all_support_rows.extend(support_rows)
        preflight_points.append(
            {
                "stage": "confirmation",
                "operating_point": confirmation_point.as_mapping(),
                "component_eligible": bool(
                    confirmation_designs and all(design.preflight_pass for design in confirmation_designs)
                ),
                "designs": [asdict(design) for design in confirmation_designs],
            }
        )
        for index, design in enumerate(confirmation_designs, start=1):
            print(f"[EXP-059] confirmation graphlet {index}/{len(confirmation_designs)}: {design.graphlet_id}", flush=True)
            payload = evaluate_design(
                config,
                design,
                point=confirmation_point,
                output_dir=output_dir,
                run_config_hash=run_config_hash,
            )
            resumed_checkpoints += int(bool(payload.get("resumed", False)))
            computed_checkpoints += int(not bool(payload.get("resumed", False)))
            component_rows.extend(payload.get("component_rows", []))
            descriptor_rows.extend(payload.get("descriptor_rows", []))

    viability_rows = build_viability_checks(config, descriptor_rows)
    transmission_rows = build_transmission_checks(config, descriptor_rows)
    topography_rows = build_topography_checks(config, descriptor_rows)
    routing_rows = build_routing_equivalence(config, descriptor_rows)

    if should_run_confirmation:
        rows, details = stage_decisions(
            config,
            stage="confirmation",
            designs=confirmation_designs,
            viability_rows=viability_rows,
            transmission_rows=transmission_rows,
            topography_rows=topography_rows,
            component_rows=component_rows,
        )
        mechanism_rows.extend(rows)
        decision_details["confirmation"] = details

    seed_rows: list[dict[str, Any]] = []
    for architecture in (ARCH_PRODUCT, ARCH_MIXED):
        seed_rows.extend(
            confirmation_seed_summary(
                transmission_rows,
                topography_rows,
                architecture=architecture,
                support_mode="structured",
            )
        )

    write_json_atomic(
        output_dir / "frozen_operating_point.json",
        {
            "selection_basis": "component_preflight_only",
            "selected_operating_point": selected_point.as_mapping() if selected_point else None,
            "reason": None if selected_point else "NO_COMPONENT_QUALIFIED_POINT",
            "accessibility_gate": {
                "metric": "basin_convergence_fraction",
                "working_weight_min": float(
                    config.get("gates", {}).get("g3_topography", {}).get("working_weight_basin_min", 0.05)
                ),
                "working_weight_max": float(
                    config.get("gates", {}).get("g3_topography", {}).get("working_weight_basin_max", 0.80)
                ),
            },
        },
    )

    preflight_payload = {
        "config_hash": run_config_hash,
        "quick": bool(args.quick),
        "sequence_length": 48,
        "sampling_key_fields": ["stage", "dataset_id", "seed", "graphlet_index"],
        "sampling_key_excludes": ["operating_point_id", "architecture", "support_mode", "gnk_weight"],
        "weights_checked_for_capacity": configured_weights(config),
        "selected_operating_point": selected_point.as_mapping() if selected_point else None,
        "points": preflight_points,
    }
    write_json_atomic(output_dir / "preflight.json", preflight_payload)
    write_csv_atomic(output_dir / "design_catalog.csv", design_catalog_rows(config, all_designs))
    write_csv_atomic(output_dir / "component_metrics.csv", component_rows)
    write_csv_atomic(output_dir / "descriptor_metrics.csv", descriptor_rows)
    write_csv_atomic(output_dir / "calibration_decision.csv", calibration_rows)
    write_csv_atomic(output_dir / "confirmation_seed_metrics.csv", seed_rows)
    write_csv_atomic(output_dir / "viability_checks.csv", viability_rows)
    write_csv_atomic(output_dir / "gnk_transmission_checks.csv", transmission_rows)
    write_csv_atomic(output_dir / "topography_checks.csv", topography_rows)
    write_csv_atomic(output_dir / "routing_equivalence.csv", routing_rows)
    write_csv_atomic(output_dir / "mechanism_decision.csv", mechanism_rows)
    write_csv_atomic(
        output_dir / "enumerable_landscape_metrics.csv",
        [row for row in descriptor_rows if row.get("stage") == "calibration"],
    )
    write_csv_atomic(output_dir / "support_catalog.csv", _deduplicate_rows(all_support_rows))

    confirmation_pass = _primary_stage_pass(mechanism_rows, "confirmation")
    if selected_point is None and "calibration" in stages:
        decision = "NO_COMPONENT_QUALIFIED_POINT"
    elif "calibration" in stages and not bool(calibration_mechanism_pass):
        decision = "SELECTED_COMPONENT_QUALIFIED_POINT_FAIL"
    elif should_run_confirmation and confirmation_pass:
        decision = "LEGACY_TWO_CHANNEL_MECHANISM_DEMONSTRATED"
    elif should_run_confirmation:
        decision = "CONFIRMATION_MECHANISM_FAIL"
    elif "confirmation" in stages:
        decision = confirmation_skipped_reason or "CONFIRMATION_NOT_RUN"
    else:
        decision = "CALIBRATION_MECHANISM_PASS"
    closed_rows = [row for row in component_rows if row.get("architecture") == ARCH_CLOSED]
    closed_control_pass = (
        bool(all(bool(row.get("closed_control_preflight_pass", False)) for row in closed_rows))
        if closed_rows
        else None
    )
    summary = {
        "experiment_id": config.get("experiment_id", "EXP-059"),
        "decision": decision,
        "scope": "internally_controllable_candidate_mechanism_only",
        "not_evidence_for": [
            "empirical_realism",
            "biological_routing_privilege",
            "optimizer_performance",
            "public_generator_readiness",
            "functional_dimensionality",
        ],
        "selected_operating_point": selected_point.as_mapping() if selected_point else None,
        "calibration_mechanism_pass": calibration_mechanism_pass,
        "confirmation_mechanism_pass": confirmation_pass if should_run_confirmation else None,
        "confirmation_skipped_reason": confirmation_skipped_reason,
        "closed_control_preflight_pass": closed_control_pass,
        "preflight_design_count": len(all_designs),
        "preflight_failure_count": sum(int(not design.preflight_pass) for design in all_designs),
        "descriptor_row_count": len(descriptor_rows),
        "checkpoint_resumed_count": resumed_checkpoints,
        "checkpoint_computed_count": computed_checkpoints,
        "fixed_fitness_tolerance": fixed_fitness_tolerance(config),
        "routing_equivalence_is_diagnostic_only": True,
        "decision_details": decision_details,
    }
    write_json_atomic(output_dir / "summary.json", summary)
    print(f"[EXP-059] decision: {decision}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
