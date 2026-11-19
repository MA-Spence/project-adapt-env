#!/usr/bin/env python3
"""EXP-064 observable-ruggedness coverage and specificity probe.

This is a forward-generation experiment.  It deliberately does not fit any
empirical fitness values.  The only empirical inputs are frozen descriptor
tables from RUN-135.  Generated graphlets are evaluated with the exact
Stage-4A binary-cube descriptor path.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import sys
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
EXTERNAL_ADAPTENV = PROJECT_ROOT / "external" / "Adapt-Env"
for _path in (SRC_ROOT, SCRIPT_DIR, EXTERNAL_ADAPTENV):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from design_matched_support_aware_probe_exp057 import (  # noqa: E402
    choose_support_sites,
    dataset_template,
    load_dataset,
    resolve_processed_root,
)
import dimensionless_sigmoid_calibration_exp062 as exp062  # noqa: E402
from ensemble_landscape_validation_stage4 import metrics_from_cube  # noqa: E402
import extensive_gnk_latent_scaling_exp060 as latent  # noqa: E402
from phot_coupled_landscape_grid_search import load_config  # noqa: E402
from project_adapt_env.extensive_gnk import (  # noqa: E402
    ExtensiveGNKField,
    additive_mutation_energy,
)
from project_adapt_env.ranked_categorical_gnk import (  # noqa: E402
    RankedCategoricalGNKField,
)


CANONICAL_DESCRIPTOR_SOURCES: dict[str, str] = {
    "local_maxima_density": "n_local_maxima",
    "accessible_path_fraction": "accessible_path_fraction_dp",
    "walsh_high_order_fraction": "walsh_high_order_fraction",
    "sign_epistasis_fraction": "sign_epistasis_fraction_local",
    "reciprocal_sign_epistasis_fraction": "reciprocal_sign_epistasis_fraction_local",
    "random_walk_lag1_autocorrelation": "random_walk_lag1_autocorrelation",
    "roughness_to_slope_ratio": "magnitude_epistasis_to_single_mutant_sd_ratio",
    "dirichlet_energy_localisation": "edge_energy_top10_fraction",
}

PUBLISHED_MEDIAN_COLUMNS: dict[str, str] = {
    "local_maxima_density": "n_local_maxima",
    "accessible_path_fraction": "accessible_path_fraction_dp",
    "walsh_high_order_fraction": "walsh_high_order_fraction",
    "sign_epistasis_fraction": "sign_epistasis_fraction_local",
    "reciprocal_sign_epistasis_fraction": "reciprocal_sign_epistasis_fraction_local",
    "random_walk_lag1_autocorrelation": "random_walk_lag1_autocorrelation",
    "roughness_to_slope_ratio": "magnitude_epistasis_to_single_mutant_sd_ratio",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def stable_uint64(*parts: Any, seed: int = 0) -> int:
    payload = json.dumps(
        [int(seed), *parts],
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little")


def json_hash(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()
    ).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: Any, default: float = math.nan) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if np.isfinite(result) else default


def strict_json_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): strict_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [strict_json_value(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if np.isfinite(number) else None
    return value


def write_json_strict(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(
            strict_json_value(payload),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv_atomic(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    pd.DataFrame(list(rows)).to_csv(temporary, index=False)
    temporary.replace(path)


def design_key(
    *,
    template: str,
    support_mode: str,
    representation: str,
    beta: float,
    p0: float,
    seed: int,
) -> str:
    """Return the paired design identity; amplitude is intentionally absent."""

    return json_hash(
        {
            "template": str(template),
            "support_mode": str(support_mode),
            "representation": str(representation),
            "beta": float(beta),
            "p0": float(p0),
            "seed": int(seed),
        }
    )


def logit(probability: float) -> float:
    probability = float(probability)
    if not 0.0 < probability < 1.0:
        raise ValueError("reference occupancy p0 must lie strictly between zero and one")
    return float(math.log(probability / (1.0 - probability)))


def observable_fitness(
    additive: np.ndarray,
    gnk: np.ndarray,
    *,
    additive_reference: float,
    sigma_a: float,
    sigma_g: float,
    beta: float,
    p0: float,
    amplitude: float,
    output_scale: float,
) -> np.ndarray:
    """Evaluate the experiment-local stability-gate operating-point override."""

    mu = logit(p0)
    logits = (
        mu
        + float(beta)
        * (np.asarray(additive, dtype=np.float64) - float(additive_reference))
        / float(sigma_a)
        + float(amplitude) * np.asarray(gnk, dtype=np.float64) / float(sigma_g)
    )
    return float(output_scale) * exp062.sigmoid(logits)


def canonical_descriptors(
    metrics: Mapping[str, Any],
    *,
    dimension: int,
) -> dict[str, float]:
    """Translate the frozen Stage-4A names to the eight EXP-064 descriptors."""

    out: dict[str, float] = {}
    for canonical, source in CANONICAL_DESCRIPTOR_SOURCES.items():
        value = finite(metrics.get(source))
        if canonical == "local_maxima_density" and np.isfinite(value):
            value /= float(1 << int(dimension))
        out[canonical] = value
    return out


def descriptors_for_values(values: np.ndarray, *, dimension: int) -> dict[str, float]:
    return canonical_descriptors(
        metrics_from_cube(np.asarray(values, dtype=np.float64), dim=int(dimension)),
        dimension=int(dimension),
    )


def _expected_input_paths(
    config: Mapping[str, Any],
    project_root: Path,
) -> dict[str, Path]:
    empirical = config["empirical"]
    return {
        key: project_root / str(empirical[key]["path"])
        for key in (
            "published_key_bands",
            "published_dataset_medians",
            "frozen_graphlet_metrics",
            "frozen_scope_audit",
        )
    }


def _verify_input_hashes(
    config: Mapping[str, Any],
    paths: Mapping[str, Path],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, path in paths.items():
        expected = str(config["empirical"][key]["sha256"])
        observed = sha256_file(path) if path.exists() else ""
        rows.append(
            {
                "check": f"sha256_{key}",
                "path": str(path),
                "expected": expected,
                "observed": observed,
                "pass": bool(path.exists() and observed == expected),
            }
        )
    return rows


def _descriptor_series(frame: pd.DataFrame, canonical: str) -> pd.Series:
    source = CANONICAL_DESCRIPTOR_SOURCES[canonical]
    values = pd.to_numeric(frame[source], errors="coerce")
    if canonical == "local_maxima_density":
        dimensions = pd.to_numeric(frame["dimension"], errors="coerce")
        values = values / np.power(2.0, dimensions)
    return values


def resolve_empirical_bands(
    config: Mapping[str, Any],
    project_root: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    """Resolve per-template bands and reproduce the frozen RES-050 summaries."""

    paths = _expected_input_paths(config, project_root)
    checks = _verify_input_hashes(config, paths)
    if not all(path.exists() for path in paths.values()):
        return [], checks, []

    raw = pd.read_csv(paths["frozen_graphlet_metrics"])
    scope = pd.read_csv(paths["frozen_scope_audit"])
    key_bands = pd.read_csv(paths["published_key_bands"])
    medians = pd.read_csv(paths["published_dataset_medians"])
    requested = [str(value) for value in config["data"]["tier1_dataset_ids"]]
    eligible = scope[
        scope["higher_order_allowed"].astype(bool)
        & scope["topography_allowed"].astype(bool)
    ]["dataset_id"].astype(str).tolist()
    quick = bool(config.get("_runtime", {}).get("quick", False))
    scope_pass = (
        set(requested).issubset(set(eligible))
        if quick
        else sorted(requested) == sorted(eligible)
    )
    checks.append(
        {
            "check": "tier1_scope_catalog",
            "expected": sorted(requested),
            "observed": sorted(eligible),
            "mode": "quick_subset" if quick else "full_exact_catalog",
            "pass": scope_pass,
        }
    )

    tolerance = float(config["empirical"]["reproduction_atol"])
    raw_all = raw[raw["dataset_id"].astype(str).isin(eligible)].copy()
    key_graph = key_bands[key_bands["source_table"] == "graphlet_metrics"].copy()
    for row in key_graph.to_dict(orient="records"):
        source = str(row["descriptor_name"])
        if source not in raw_all:
            continue
        values = pd.to_numeric(raw_all[source], errors="coerce").dropna()
        observed = {
            "q05": float(values.quantile(0.05)),
            "median": float(values.quantile(0.50)),
            "q95": float(values.quantile(0.95)),
        }
        expected = {
            "q05": float(row["empirical_q05"]),
            "median": float(row["empirical_median"]),
            "q95": float(row["empirical_q95"]),
        }
        maximum_error = max(abs(observed[key] - expected[key]) for key in expected)
        checks.append(
            {
                "check": f"published_pooled_band_{source}",
                "maximum_absolute_error": maximum_error,
                "tolerance": tolerance,
                "pass": bool(maximum_error <= tolerance),
            }
        )

    raw = raw_all[raw_all["dataset_id"].astype(str).isin(requested)].copy()
    published_by_id = medians.set_index("dataset_id")
    for dataset_id in requested:
        group = raw[raw["dataset_id"].astype(str) == dataset_id]
        if dataset_id not in published_by_id.index:
            checks.append(
                {
                    "check": f"published_dataset_median_{dataset_id}",
                    "pass": False,
                    "reason": "missing published row",
                }
            )
            continue
        errors: list[float] = []
        for canonical, column in PUBLISHED_MEDIAN_COLUMNS.items():
            if column not in published_by_id or column not in group:
                continue
            observed_values = pd.to_numeric(group[column], errors="coerce").dropna()
            if observed_values.empty:
                continue
            observed = float(observed_values.median())
            expected = finite(published_by_id.loc[dataset_id, column])
            if np.isfinite(observed) and np.isfinite(expected):
                errors.append(abs(observed - expected))
        maximum_error = max(errors) if errors else math.inf
        checks.append(
            {
                "check": f"published_dataset_median_{dataset_id}",
                "maximum_absolute_error": maximum_error,
                "tolerance": tolerance,
                "pass": bool(maximum_error <= tolerance),
            }
        )

    directions = {
        str(key): str(value)
        for key, value in config["descriptors"]["rugged_direction"].items()
    }
    bands: list[dict[str, Any]] = []
    for dataset_id in requested:
        group = raw[raw["dataset_id"].astype(str) == dataset_id]
        for canonical in CANONICAL_DESCRIPTOR_SOURCES:
            values = _descriptor_series(group, canonical).dropna()
            if values.empty:
                bands.append(
                    {
                        "template": dataset_id,
                        "descriptor": canonical,
                        "source_column": CANONICAL_DESCRIPTOR_SOURCES[canonical],
                        "rugged_direction": directions[canonical],
                        "n_graphlets": 0,
                        "empirical_q05": math.nan,
                        "empirical_median": math.nan,
                        "empirical_q90": math.nan,
                        "empirical_q95": math.nan,
                        "empirical_rugged_q90": math.nan,
                        "descriptor_available": False,
                    }
                )
                continue
            direction = directions[canonical]
            rugged_probability = 0.90 if direction == "increase" else 0.10
            bands.append(
                {
                    "template": dataset_id,
                    "descriptor": canonical,
                    "source_column": CANONICAL_DESCRIPTOR_SOURCES[canonical],
                    "rugged_direction": direction,
                    "n_graphlets": int(len(values)),
                    "empirical_q05": float(values.quantile(0.05)),
                    "empirical_median": float(values.quantile(0.50)),
                    "empirical_q90": float(values.quantile(0.90)),
                    "empirical_q95": float(values.quantile(0.95)),
                    "empirical_rugged_q90": float(values.quantile(rugged_probability)),
                    "descriptor_available": True,
                }
            )

    local = pd.DataFrame(bands)
    local = local[local["descriptor"] == "local_maxima_density"].sort_values(
        ["empirical_median", "template"], ascending=[False, True]
    )
    rugged_count = max(
        1,
        int(math.ceil(float(config["gates"]["rugged_end_top_fraction"]) * len(requested))),
    )
    rugged_templates = local.head(rugged_count)["template"].astype(str).tolist()
    rugged_set = set(rugged_templates)
    for row in bands:
        row["rugged_end_template"] = bool(row["template"] in rugged_set)

    synthetic = np.asarray(
        [float(node.bit_count()) for node in range(1 << 4)],
        dtype=np.float64,
    )
    synthetic_metrics = metrics_from_cube(synthetic, dim=4)
    synthetic_pass = bool(
        abs(finite(synthetic_metrics.get("n_local_maxima")) - 1.0) <= 1e-12
        and abs(finite(synthetic_metrics.get("accessible_path_fraction_dp")) - 1.0)
        <= 1e-12
        and abs(finite(synthetic_metrics.get("walsh_high_order_fraction"), 0.0))
        <= 1e-12
    )
    checks.append(
        {
            "check": "known_additive_cube_descriptor_control",
            "n_local_maxima": finite(synthetic_metrics.get("n_local_maxima")),
            "accessible_path_fraction": finite(
                synthetic_metrics.get("accessible_path_fraction_dp")
            ),
            "walsh_high_order_fraction": finite(
                synthetic_metrics.get("walsh_high_order_fraction")
            ),
            "pass": synthetic_pass,
        }
    )
    return bands, checks, rugged_templates


def load_template_catalog(
    config: Mapping[str, Any],
    project_root: Path,
) -> list[dict[str, Any]]:
    root = resolve_processed_root(config["data"], project_root)
    raw_graphlets = pd.read_csv(
        project_root / str(config["empirical"]["frozen_graphlet_metrics"]["path"])
    )
    rows: list[dict[str, Any]] = []
    for dataset_id in config["data"]["tier1_dataset_ids"]:
        frame = load_dataset(root, str(dataset_id))
        template = dataset_template(frame, str(dataset_id))
        empirical = raw_graphlets[
            raw_graphlets["dataset_id"].astype(str) == str(dataset_id)
        ]
        dimensions = sorted(
            int(value)
            for value in pd.to_numeric(empirical["dimension"], errors="coerce")
            .dropna()
            .unique()
        )
        if len(dimensions) != 1:
            raise ValueError(f"{dataset_id} does not have one frozen graphlet dimension")
        rows.append(
            {
                **template,
                "frozen_graphlet_dimension": dimensions[0],
                "frozen_graphlet_count": int(len(empirical)),
            }
        )
    return rows


def arm_parameters(
    config: Mapping[str, Any],
    arm: str,
) -> tuple[str, int]:
    for row in config["generator"]["representations"]:
        if str(row["id"]) == str(arm):
            return str(row["contrast_family"]), int(row["rank"])
    raise KeyError(arm)


def build_ranked_field(
    *,
    matrix: np.ndarray,
    reference: np.ndarray,
    seed: int,
    arm: str,
    config: Mapping[str, Any],
) -> RankedCategoricalGNKField:
    family, rank = arm_parameters(config, arm)
    generator = config["generator"]
    return RankedCategoricalGNKField(
        additive_matrix=matrix,
        reference_sequence=reference,
        seed=int(seed),
        orders=tuple(int(value) for value in generator["orders"]),
        order_weights=tuple(float(value) for value in generator["order_weights"]),
        target_degrees=tuple(float(value) for value in generator["target_degrees"]),
        contrast_family=family,
        representation_rank=rank,
    )


def _support_rows(
    field: RankedCategoricalGNKField,
    *,
    template: str,
    seed: int,
    arm: str,
    support_mode: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for order in field.orders:
        degrees = field.support_degrees(order)
        rows.append(
            {
                "template": template,
                "seed": int(seed),
                "representation": arm,
                "support_mode": support_mode,
                "order": int(order),
                "n_supports": int(field.term_counts[order]),
                "degree_min": int(np.min(degrees)),
                "degree_max": int(np.max(degrees)),
                "degree_mean": float(np.mean(degrees)),
                "expected_mutation_energy": float(
                    field.expected_order_mutation_energy[order]
                ),
                "support_fingerprint": field.support_fingerprint,
                "field_fingerprint": field.fingerprint,
            }
        )
    return rows


def _reference_for_task(template: str, seed: int, length: int) -> np.ndarray:
    rng = np.random.default_rng(
        stable_uint64("EXP-064-reference", template, int(seed)) % (2**32)
    )
    return rng.integers(0, 20, size=int(length), dtype=np.int64)


def _route_sequences(
    sequences: np.ndarray,
    reference: np.ndarray,
    permutation: np.ndarray,
) -> np.ndarray:
    """Route categorical mutation offsets through a bijective site map."""

    array = np.asarray(sequences, dtype=np.int64)
    routed = np.repeat(reference[None, :], array.shape[0], axis=0)
    for source, target in enumerate(np.asarray(permutation, dtype=np.int64)):
        offsets = (array[:, source] - int(reference[source])) % 20
        routed[:, int(target)] = (int(reference[int(target)]) + offsets) % 20
    return routed


def _graphlet_panel(
    field: RankedCategoricalGNKField,
    reference: np.ndarray,
    *,
    template: str,
    seed: int,
    dimension: int,
    count: int,
) -> tuple[list[np.ndarray], list[dict[str, Any]]]:
    support_terms = [
        {
            "term_id": int(support_id),
            "support_positions": ";".join(str(value) for value in support),
        }
        for support_id, support in zip(field.support_ids(), field.supports())
    ]
    panels: list[np.ndarray] = []
    audits: list[dict[str, Any]] = []
    for index in range(int(count)):
        sites = choose_support_sites(
            support_terms,
            dimension=int(dimension),
            graphlet_index=int(index),
            sequence_length=int(reference.size),
            seed=stable_uint64("EXP-064-sites", template, int(seed)) % (2**32),
        )
        alternatives: list[int] = []
        for axis, site in enumerate(sites):
            raw = stable_uint64(
                "EXP-064-allele",
                template,
                int(seed),
                int(index),
                int(axis),
                int(site),
            ) % 19
            alternative = int(raw)
            if alternative >= int(reference[int(site)]):
                alternative += 1
            alternatives.append(alternative)
        sequences = np.repeat(
            reference[None, :],
            1 << int(dimension),
            axis=0,
        )
        for node in range(1 << int(dimension)):
            for axis, site in enumerate(sites):
                if (node >> axis) & 1:
                    sequences[node, int(site)] = int(alternatives[axis])
        panel_hash = latent.array_hash(sequences)
        panels.append(sequences)
        audits.append(
            {
                "template": template,
                "seed": int(seed),
                "graphlet_index": int(index),
                "dimension": int(dimension),
                "n_nodes": int(len(sequences)),
                "varied_sites": ";".join(str(value) for value in sites),
                "alternative_alleles": ";".join(str(value) for value in alternatives),
                "graphlet_hash": panel_hash,
                "support_aware": True,
            }
        )
    return panels, audits


def _field_audit(
    field: RankedCategoricalGNKField,
    *,
    arm: str,
    template: str,
    seed: int,
    tolerance: float,
) -> dict[str, Any]:
    target_by_order = {
        order: field.order_weights[order] * field.additive_mutation_energy
        for order in field.orders
    }
    errors = [
        abs(field.expected_order_mutation_energy[order] - target_by_order[order])
        / target_by_order[order]
        for order in field.orders
    ]
    effective_ranks: list[int] = []
    centered: list[float] = []
    finite_terms = True
    dense_ties: list[int] = []
    for support_id in field.support_ids():
        eigenvalues = np.linalg.eigvalsh(field.support_operator_gram(support_id))
        threshold = max(float(np.max(eigenvalues)), 1.0) * 1e-10
        effective_ranks.append(int(np.sum(eigenvalues > threshold)))
    for term in field.terms:
        for contrast in term.contrasts:
            centered.append(abs(float(np.mean(contrast))))
            finite_terms = finite_terms and bool(np.isfinite(contrast).all())
            if field.contrast_family == "dense_q20":
                dense_ties.append(int(20 - len(np.unique(contrast))))
    return {
        "template": template,
        "seed": int(seed),
        "representation": arm,
        "support_fingerprint": field.support_fingerprint,
        "field_fingerprint": field.fingerprint,
        "maximum_operator_energy_relative_error": max(errors),
        "operator_energy_pass": bool(max(errors) <= float(tolerance)),
        "maximum_contrast_mean_abs": max(centered) if centered else math.nan,
        "centered_and_finite": bool(
            finite_terms and centered and max(centered) <= 1e-12
        ),
        "minimum_effective_rank": min(effective_ranks) if effective_ranks else 0,
        "effective_rank_pass": bool(
            field.representation_rank < 8
            or (effective_ranks and min(effective_ranks) >= 4)
        ),
        "maximum_dense_q20_ties": max(dense_ties) if dense_ties else 0,
        "dense_q20_tie_pass": bool(not dense_ties or max(dense_ties) == 0),
    }


def _checkpoint_name(template: str, seed: int) -> str:
    safe = "".join(char if char.isalnum() or char in "._-" else "_" for char in template)
    return f"{safe}__seed-{int(seed)}"


def _task_checkpoint_dir(output_dir: Path, template: str, seed: int) -> Path:
    return output_dir / "checkpoints" / _checkpoint_name(template, seed)


def _load_task_checkpoint(path: Path) -> dict[str, list[dict[str, Any]]] | None:
    done = path / "done.json"
    if not done.exists():
        return None
    return {
        "metrics": pd.read_csv(path / "generated_descriptor_metrics.csv").to_dict(
            orient="records"
        ),
        "supports": pd.read_csv(path / "support_catalog.csv").to_dict(orient="records"),
        "graphlets": pd.read_csv(path / "graphlet_audit.csv").to_dict(orient="records"),
        "audits": pd.read_csv(path / "technical_audit.csv").to_dict(orient="records"),
    }


def _write_task_checkpoint(
    path: Path,
    payload: Mapping[str, Sequence[Mapping[str, Any]]],
) -> None:
    path.mkdir(parents=True, exist_ok=True)
    done = path / "done.json"
    if done.exists():
        done.unlink()
    write_csv_atomic(path / "generated_descriptor_metrics.csv", payload["metrics"])
    write_csv_atomic(path / "support_catalog.csv", payload["supports"])
    write_csv_atomic(path / "graphlet_audit.csv", payload["graphlets"])
    write_csv_atomic(path / "technical_audit.csv", payload["audits"])
    write_json_strict(
        done,
        {
            "n_metric_rows": len(payload["metrics"]),
            "n_support_rows": len(payload["supports"]),
            "n_graphlet_rows": len(payload["graphlets"]),
            "n_audit_rows": len(payload["audits"]),
        },
    )


def _copy_resume_checkpoints(output_dir: Path) -> list[str]:
    source = os.environ.get("LABPROJ_RESUME_OUTPUT_DIR")
    if not source:
        return []
    source_root = Path(source) / "checkpoints"
    destination_root = output_dir / "checkpoints"
    if not source_root.exists() or source_root.resolve() == destination_root.resolve():
        return []
    destination_root.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for child in sorted(source_root.iterdir()):
        if not child.is_dir() or not (child / "done.json").exists():
            continue
        destination = destination_root / child.name
        if destination.exists():
            continue
        shutil.copytree(child, destination)
        copied.append(child.name)
    return copied


def generation_task(payload: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    config = payload["config"]
    template = dict(payload["template"])
    template_id = str(template["dataset_id"])
    seed = int(payload["seed"])
    output_dir = Path(str(payload["output_dir"]))
    checkpoint = _task_checkpoint_dir(output_dir, template_id, seed)
    cached = _load_task_checkpoint(checkpoint)
    if cached is not None:
        return cached

    generator = config["generator"]
    reference = _reference_for_task(
        template_id,
        seed,
        int(generator["sequence_length"]),
    )
    core_seed = int(generator["core_seed_offset"]) + seed
    core = exp062.build_core(reference, core_seed, config["core_helper"])
    matrix = np.asarray(core.stability_matrix, dtype=np.float64)
    additive_reference = float(latent.additive_scores(matrix, reference[None, :])[0])
    sigma_g = math.sqrt(additive_mutation_energy(matrix))
    sigma_a = float(config["readout"]["sigma_a"])
    if not np.isfinite(sigma_g) or sigma_g <= 0.0:
        raise ValueError("non-positive GNK mutation scale")

    arms = [str(row["id"]) for row in generator["representations"]]
    field_seed = int(generator["gnk_seed_offset"]) + seed
    random_seed = field_seed + int(generator["random_support_seed_offset"])
    fields = {
        (arm, "structured"): build_ranked_field(
            matrix=matrix,
            reference=reference,
            seed=field_seed,
            arm=arm,
            config=config,
        )
        for arm in arms
    }
    fields.update(
        {
            (arm, "random"): build_ranked_field(
                matrix=matrix,
                reference=reference,
                seed=random_seed,
                arm=arm,
                config=config,
            )
            for arm in arms
        }
    )
    baseline = fields[("binary_r1", "structured")]
    panels, graphlet_rows = _graphlet_panel(
        baseline,
        reference,
        template=template_id,
        seed=seed,
        dimension=int(template["frozen_graphlet_dimension"]),
        count=int(generator["graphlets_per_template_seed"]),
    )
    flattened = np.concatenate(panels, axis=0)
    panel_hash = json_hash([row["graphlet_hash"] for row in graphlet_rows])
    additive_flat = latent.additive_scores(matrix, flattened)
    additive_panels = np.split(
        additive_flat,
        np.cumsum([len(panel) for panel in panels])[:-1],
    )

    permutation_rng = np.random.default_rng(
        stable_uint64("EXP-064-shuffled-route", template_id, seed) % (2**32)
    )
    shuffled_permutation = permutation_rng.permutation(reference.size)
    shuffled_flattened = _route_sequences(
        flattened,
        reference,
        shuffled_permutation,
    )

    support_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    energy_tolerance = float(config["technical"]["operator_energy_relative_error_max"])
    for arm in arms:
        structured = fields[(arm, "structured")]
        random = fields[(arm, "random")]
        support_rows.extend(
            _support_rows(
                structured,
                template=template_id,
                seed=seed,
                arm=arm,
                support_mode="structured",
            )
        )
        support_rows.extend(
            _support_rows(
                random,
                template=template_id,
                seed=seed,
                arm=arm,
                support_mode="random",
            )
        )
        audit_rows.append(
            _field_audit(
                structured,
                arm=arm,
                template=template_id,
                seed=seed,
                tolerance=energy_tolerance,
            )
        )

    support_match = len(
        {fields[(arm, "structured")].support_fingerprint for arm in arms}
    ) == 1
    audit_rows.append(
        {
            "template": template_id,
            "seed": seed,
            "representation": "all",
            "check": "support_hash_match_across_representation_arms",
            "pass": bool(support_match),
            "graphlet_panel_hash": panel_hash,
        }
    )

    legacy = ExtensiveGNKField(
        additive_matrix=matrix,
        reference_sequence=reference,
        seed=field_seed,
        orders=tuple(int(value) for value in generator["orders"]),
        order_weights=tuple(float(value) for value in generator["order_weights"]),
        target_degrees=tuple(float(value) for value in generator["target_degrees"]),
    )
    replay_difference = float(
        np.max(
            np.abs(
                baseline.score_batch(flattened)
                - legacy.score_batch(flattened)
            )
        )
    )
    audit_rows.append(
        {
            "template": template_id,
            "seed": seed,
            "representation": "binary_r1",
            "check": "binary_r1_legacy_replay",
            "maximum_absolute_difference": replay_difference,
            "pass": bool(replay_difference == 0.0),
        }
    )

    dimensions = int(template["frozen_graphlet_dimension"])
    operating_betas = [float(value) for value in generator["beta_grid"]]
    operating_p0 = [float(value) for value in generator["p0_grid"]]
    amplitudes = [float(value) for value in generator["amplitudes"]]
    output_scale = float(config["readout"]["output_scale"])
    support_modes = [str(value) for value in generator["support_modes"]]
    metric_rows: list[dict[str, Any]] = []
    for arm in arms:
        for support_mode in support_modes:
            if support_mode == "structured":
                score_input = flattened
                field = fields[(arm, "structured")]
            elif support_mode == "shuffled":
                score_input = shuffled_flattened
                field = fields[(arm, "structured")]
            elif support_mode == "random":
                score_input = flattened
                field = fields[(arm, "random")]
            else:
                raise ValueError(f"unsupported support mode: {support_mode}")
            gnk_flat = field.score_batch(score_input)
            gnk_panels = np.split(
                gnk_flat,
                np.cumsum([len(panel) for panel in panels])[:-1],
            )
            for beta, p0 in itertools.product(operating_betas, operating_p0):
                identity = design_key(
                    template=template_id,
                    support_mode=support_mode,
                    representation=arm,
                    beta=beta,
                    p0=p0,
                    seed=seed,
                )
                for amplitude in amplitudes:
                    descriptor_pieces: dict[str, list[float]] = {
                        descriptor: [] for descriptor in CANONICAL_DESCRIPTOR_SOURCES
                    }
                    for additive_values, gnk_values in zip(
                        additive_panels,
                        gnk_panels,
                    ):
                        fitness = observable_fitness(
                            additive_values,
                            gnk_values,
                            additive_reference=additive_reference,
                            sigma_a=sigma_a,
                            sigma_g=sigma_g,
                            beta=beta,
                            p0=p0,
                            amplitude=amplitude,
                            output_scale=output_scale,
                        )
                        descriptors = descriptors_for_values(
                            fitness,
                            dimension=dimensions,
                        )
                        for descriptor, value in descriptors.items():
                            if np.isfinite(value):
                                descriptor_pieces[descriptor].append(float(value))
                    collapsed = {
                        descriptor: (
                            float(np.median(values)) if values else math.nan
                        )
                        for descriptor, values in descriptor_pieces.items()
                    }
                    metric_rows.append(
                        {
                            "template": template_id,
                            "seed": seed,
                            "core_seed": core_seed,
                            "gnk_seed": field_seed,
                            "design_key": identity,
                            "graphlet_panel_hash": panel_hash,
                            "representation": arm,
                            "support_mode": support_mode,
                            "beta": beta,
                            "p0": p0,
                            "mu": logit(p0),
                            "amplitude": amplitude,
                            "n_graphlets": len(panels),
                            "graphlet_dimension": dimensions,
                            "sigma_a": sigma_a,
                            "sigma_g": sigma_g,
                            "additive_reference": additive_reference,
                            **collapsed,
                        }
                    )

    payload_out = {
        "metrics": metric_rows,
        "supports": support_rows,
        "graphlets": graphlet_rows,
        "audits": audit_rows,
    }
    _write_task_checkpoint(checkpoint, payload_out)
    return payload_out


def run_generation(
    config: Mapping[str, Any],
    templates: Sequence[Mapping[str, Any]],
    *,
    output_dir: Path,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    copied = _copy_resume_checkpoints(output_dir)
    payloads = [
        {
            "config": config,
            "template": dict(template),
            "seed": int(seed),
            "output_dir": str(output_dir),
        }
        for template in templates
        for seed in config["generator"]["seeds"]
    ]
    workers = max(
        1,
        min(
            int(os.environ.get("LABPROJ_WORKERS", "1")),
            len(payloads),
        ),
    )
    results: list[dict[str, list[dict[str, Any]]]] = []
    if workers == 1:
        for index, payload in enumerate(payloads, start=1):
            results.append(generation_task(payload))
            print(f"[EXP-064] generation {index}/{len(payloads)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(generation_task, payload) for payload in payloads]
            for index, future in enumerate(as_completed(futures), start=1):
                results.append(future.result())
                print(f"[EXP-064] generation {index}/{len(payloads)}", flush=True)
    metrics = [row for result in results for row in result["metrics"]]
    supports = [row for result in results for row in result["supports"]]
    graphlets = [row for result in results for row in result["graphlets"]]
    audits = [row for result in results for row in result["audits"]]
    resume = {
        "copied_checkpoint_directories": copied,
        "n_tasks": len(payloads),
        "n_checkpoint_directories": len(
            list((output_dir / "checkpoints").glob("*/done.json"))
        ),
        "workers": workers,
    }
    return metrics, supports, graphlets, audits, resume


def _row_within_band(
    row: Mapping[str, Any],
    band_by_descriptor: Mapping[str, Mapping[str, Any]],
) -> tuple[int, list[str]]:
    in_band: list[str] = []
    for descriptor, band in band_by_descriptor.items():
        value = finite(row.get(descriptor))
        if (
            np.isfinite(value)
            and value >= finite(band["empirical_q05"]) - 1e-12
            and value <= finite(band["empirical_q95"]) + 1e-12
        ):
            in_band.append(descriptor)
    return len(in_band), in_band


def coverage_and_decision(
    config: Mapping[str, Any],
    generated_rows: Sequence[Mapping[str, Any]],
    empirical_bands: Sequence[Mapping[str, Any]],
    technical_checks: Sequence[Mapping[str, Any]],
    *,
    quick: bool,
) -> dict[str, Any]:
    generated = pd.DataFrame(generated_rows)
    bands = pd.DataFrame(empirical_bands)
    primary = generated[
        generated["representation"].eq("binary_r8")
        & generated["support_mode"].eq("structured")
        & pd.to_numeric(generated["amplitude"], errors="coerce").gt(0.0)
    ].copy()
    rugged_templates = sorted(
        set(
            bands.loc[
                bands["rugged_end_template"].astype(bool),
                "template",
            ].astype(str)
        )
    )
    coverage_rows: list[dict[str, Any]] = []
    for band in bands.to_dict(orient="records"):
        template = str(band["template"])
        descriptor = str(band["descriptor"])
        values = pd.to_numeric(
            primary.loc[primary["template"].eq(template), descriptor],
            errors="coerce",
        ).dropna()
        minimum = float(values.min()) if not values.empty else math.nan
        maximum = float(values.max()) if not values.empty else math.nan
        median = finite(band["empirical_median"])
        rugged_target = finite(band["empirical_rugged_q90"])
        direction = str(band["rugged_direction"])
        contains_median = bool(
            np.isfinite(minimum)
            and minimum <= median + 1e-12
            and maximum >= median - 1e-12
        )
        reaches_rugged = bool(
            np.isfinite(minimum)
            and (
                maximum >= rugged_target - 1e-12
                if direction == "increase"
                else minimum <= rugged_target + 1e-12
            )
        )
        if not np.isfinite(rugged_target):
            shortfall = math.nan
        elif direction == "increase" and np.isfinite(maximum):
            shortfall = max(0.0, rugged_target - maximum)
        elif direction == "decrease" and np.isfinite(minimum):
            shortfall = max(0.0, minimum - rugged_target)
        else:
            shortfall = math.inf
        coverage_rows.append(
            {
                **band,
                "reachable_min": minimum,
                "reachable_max": maximum,
                "contains_empirical_median": contains_median,
                "reaches_empirical_rugged_q90": reaches_rugged,
                "cover": bool(contains_median and reaches_rugged),
                "rugged_shortfall": shortfall,
            }
        )

    minimum_joint = int(config["gates"]["joint_descriptors_min"])
    accessibility_floor = float(config["gates"]["accessibility_floor"])
    reachability_rows: list[dict[str, Any]] = []
    for template in rugged_templates:
        template_bands = {
            str(row["descriptor"]): row
            for row in bands[bands["template"].eq(template)].to_dict(orient="records")
        }
        candidates: list[dict[str, Any]] = []
        for row in primary[primary["template"].eq(template)].to_dict(orient="records"):
            count, names = _row_within_band(row, template_bands)
            candidates.append(
                {
                    **row,
                    "n_descriptors_in_band": count,
                    "descriptors_in_band": ";".join(sorted(names)),
                    "joint_reach_pass": bool(count >= minimum_joint),
                    "navigability_pass": bool(
                        finite(row.get("accessible_path_fraction"))
                        >= accessibility_floor
                    ),
                }
            )
        candidates.sort(
            key=lambda row: (
                int(row["n_descriptors_in_band"]),
                int(row["navigability_pass"]),
                finite(row.get("accessible_path_fraction"), -math.inf),
            ),
            reverse=True,
        )
        best = candidates[0] if candidates else {
            "template": template,
            "n_descriptors_in_band": 0,
            "joint_reach_pass": False,
            "navigability_pass": False,
        }
        navigable = [
            row
            for row in candidates
            if row["joint_reach_pass"] and row["navigability_pass"]
        ]
        reachability_rows.append(
            {
                **best,
                "any_joint_reach_pass": bool(
                    any(row["joint_reach_pass"] for row in candidates)
                ),
                "any_navigable_joint_reach_pass": bool(navigable),
            }
        )

    rugged_fraction = (
        float(
            np.mean(
                [bool(row["any_joint_reach_pass"]) for row in reachability_rows]
            )
        )
        if reachability_rows
        else 0.0
    )
    navigable_fraction = (
        float(
            np.mean(
                [
                    bool(row["any_navigable_joint_reach_pass"])
                    for row in reachability_rows
                ]
            )
        )
        if reachability_rows
        else 0.0
    )
    landscape_fraction_floor = float(
        config["gates"]["rugged_end_landscape_fraction_min"]
    )
    g2_pass = bool(rugged_fraction >= landscape_fraction_floor)
    g4_pass = bool(navigable_fraction >= landscape_fraction_floor)

    cover_frame = pd.DataFrame(coverage_rows)
    cover_fractions = (
        cover_frame[cover_frame["template"].isin(rugged_templates)]
        .groupby("template")["cover"]
        .mean()
    )
    g1_pass = bool(
        len(cover_fractions) == len(rugged_templates)
        and float(
            np.mean(
                cover_fractions
                >= float(config["gates"]["range_descriptor_fraction_min"])
            )
        )
        >= landscape_fraction_floor
    )

    deltas = {
        str(key): float(value)
        for key, value in config["gates"]["specificity_deltas"].items()
    }
    increasing = [
        descriptor
        for descriptor, direction in config["descriptors"]["rugged_direction"].items()
        if str(direction) == "increase"
    ]
    specificity_rows: list[dict[str, Any]] = []
    for template in rugged_templates:
        template_generated = generated[generated["template"].eq(template)]
        composite = template_generated[
            template_generated["representation"].eq("binary_r8")
            & template_generated["support_mode"].eq("structured")
            & pd.to_numeric(
                template_generated["amplitude"], errors="coerce"
            ).gt(0.0)
        ]
        additive = template_generated[
            template_generated["representation"].eq("binary_r8")
            & template_generated["support_mode"].eq("structured")
            & pd.to_numeric(
                template_generated["amplitude"], errors="coerce"
            ).eq(0.0)
        ]
        shuffled = template_generated[
            template_generated["representation"].eq("binary_r8")
            & template_generated["support_mode"].eq("shuffled")
            & pd.to_numeric(
                template_generated["amplitude"], errors="coerce"
            ).gt(0.0)
        ]
        random = template_generated[
            template_generated["representation"].eq("binary_r8")
            & template_generated["support_mode"].eq("random")
            & pd.to_numeric(
                template_generated["amplitude"], errors="coerce"
            ).gt(0.0)
        ]
        for descriptor in increasing:
            comp_max = finite(pd.to_numeric(composite[descriptor], errors="coerce").max())
            additive_max = finite(pd.to_numeric(additive[descriptor], errors="coerce").max())
            shuffled_max = finite(pd.to_numeric(shuffled[descriptor], errors="coerce").max())
            random_max = finite(pd.to_numeric(random[descriptor], errors="coerce").max())
            delta = deltas[descriptor]
            rows = {
                "template": template,
                "descriptor": descriptor,
                "composite_max": comp_max,
                "additive_core_max": additive_max,
                "shuffled_support_max": shuffled_max,
                "random_support_max": random_max,
                "required_delta": delta,
                "composite_minus_additive": comp_max - additive_max,
                "composite_minus_shuffled": comp_max - shuffled_max,
                "composite_minus_random": comp_max - random_max,
            }
            rows["specificity_pass"] = bool(
                all(
                    finite(rows[key], -math.inf) >= delta - 1e-12
                    for key in (
                        "composite_minus_additive",
                        "composite_minus_shuffled",
                        "composite_minus_random",
                    )
                )
            )
            specificity_rows.append(rows)

    specificity = pd.DataFrame(specificity_rows)
    if specificity.empty:
        specificity_fraction = 0.0
    else:
        per_descriptor = specificity.groupby("descriptor")["specificity_pass"].mean()
        descriptor_pass = (
            per_descriptor
            >= float(config["gates"]["specificity_landscape_fraction_min"])
        )
        specificity_fraction = float(np.mean(descriptor_pass))
    g3_pass = bool(
        specificity_fraction
        > float(config["gates"]["specificity_descriptor_majority_fraction"])
    )

    selected_rows: list[dict[str, Any]] = []
    if not primary.empty:
        selections: list[dict[str, Any]] = []
        band_maps = {
            template: {
                str(row["descriptor"]): row
                for row in bands[bands["template"].eq(template)].to_dict(
                    orient="records"
                )
            }
            for template in rugged_templates
        }
        for row in primary[primary["template"].isin(rugged_templates)].to_dict(
            orient="records"
        ):
            count, _names = _row_within_band(
                row,
                band_maps[str(row["template"])],
            )
            selections.append(
                {
                    "template": row["template"],
                    "seed": int(row["seed"]),
                    "beta": float(row["beta"]),
                    "p0": float(row["p0"]),
                    "amplitude": float(row["amplitude"]),
                    "n_descriptors_in_band": count,
                    "navigable": bool(
                        finite(row["accessible_path_fraction"])
                        >= accessibility_floor
                    ),
                }
            )
        selected = pd.DataFrame(selections)
        if not selected.empty:
            grouped = (
                selected.groupby(["beta", "p0", "amplitude"], as_index=False)
                .agg(
                    mean_descriptors_in_band=("n_descriptors_in_band", "mean"),
                    maximum_descriptors_in_band=("n_descriptors_in_band", "max"),
                    rugged_template_fraction_joint=(
                        "n_descriptors_in_band",
                        lambda values: float(np.mean(np.asarray(values) >= minimum_joint)),
                    ),
                    navigable_fraction=("navigable", "mean"),
                    n_rows=("template", "size"),
                )
                .sort_values(
                    [
                        "rugged_template_fraction_joint",
                        "mean_descriptors_in_band",
                        "navigable_fraction",
                    ],
                    ascending=False,
                )
            )
            selected_rows = grouped.to_dict(orient="records")

    empirical_prefixes = (
        "sha256_",
        "tier1_scope_catalog",
        "published_pooled_band_",
        "published_dataset_median_",
    )
    empirical_check_rows = [
        row
        for row in technical_checks
        if str(row.get("check", "")).startswith(empirical_prefixes)
    ]
    execution_check_rows = [
        row
        for row in technical_checks
        if row not in empirical_check_rows
    ]
    empirical_pass = bool(
        empirical_check_rows
        and all(bool(row.get("pass", False)) for row in empirical_check_rows)
    )
    execution_pass = bool(
        execution_check_rows
        and all(bool(row.get("pass", False)) for row in execution_check_rows)
    )
    technical_pass = bool(empirical_pass and execution_pass)
    coverage_package = bool(g1_pass and g2_pass)
    if quick:
        decision = "QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION"
    elif not execution_pass:
        decision = "EXECUTION_OR_PAIRING_INVALID"
    elif not empirical_pass:
        decision = "EMPIRICAL_BANDS_NOT_REPRODUCED"
    elif not coverage_package:
        decision = "OBSERVABLE_CONSTRUCTION_LIMIT"
    elif not g3_pass:
        decision = "COVERAGE_NONSPECIFIC"
    elif not g4_pass:
        decision = "NAVIGABILITY_DEGENERATE"
    else:
        decision = "PROCEED_TO_FITTING"

    decision_rows = [
        {"gate": "execution_and_pairing", "pass": execution_pass},
        {"gate": "empirical_band_reproduction", "pass": empirical_pass},
        {"gate": "G1_range_coverage", "pass": g1_pass},
        {
            "gate": "G2_joint_rugged_end_reachability",
            "pass": g2_pass,
            "observed_fraction": rugged_fraction,
            "required_fraction": landscape_fraction_floor,
        },
        {
            "gate": "G3_specificity",
            "pass": g3_pass,
            "descriptor_pass_fraction": specificity_fraction,
        },
        {
            "gate": "G4_navigable_rugged_sanity",
            "pass": g4_pass,
            "observed_fraction": navigable_fraction,
            "required_fraction": landscape_fraction_floor,
        },
        {"gate": "terminal_decision", "pass": True, "decision": decision},
    ]
    return {
        "coverage_rows": coverage_rows,
        "reachability_rows": reachability_rows,
        "specificity_rows": specificity_rows,
        "selected_rows": selected_rows,
        "decision_rows": decision_rows,
        "decision": decision,
        "gate_summary": {
            "technical_pass": technical_pass,
            "execution_and_pairing_pass": execution_pass,
            "empirical_band_reproduction_pass": empirical_pass,
            "g1_range_coverage_pass": g1_pass,
            "g2_joint_rugged_end_pass": g2_pass,
            "g3_specificity_pass": g3_pass,
            "g4_navigability_pass": g4_pass,
            "rugged_end_joint_fraction": rugged_fraction,
            "rugged_end_navigable_fraction": navigable_fraction,
            "specificity_descriptor_pass_fraction": specificity_fraction,
        },
    }


def _technical_rows(
    empirical_checks: Sequence[Mapping[str, Any]],
    generation_audits: Sequence[Mapping[str, Any]],
    graphlet_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows = [dict(row) for row in empirical_checks]
    for row in generation_audits:
        if "pass" in row:
            rows.append(dict(row))
        else:
            passed = bool(
                row.get("operator_energy_pass")
                and row.get("centered_and_finite")
                and row.get("effective_rank_pass")
                and row.get("dense_q20_tie_pass")
            )
            rows.append(
                {
                    **row,
                    "check": "representation_manipulation_package",
                    "pass": passed,
                }
            )
    graphlets = pd.DataFrame(graphlet_rows)
    rows.append(
        {
            "check": "graphlet_hashes_unique_by_template_seed_index",
            "n_rows": int(len(graphlets)),
            "n_unique_keys": int(
                graphlets[
                    ["template", "seed", "graphlet_index"]
                ].drop_duplicates().shape[0]
            )
            if not graphlets.empty
            else 0,
            "pass": bool(
                not graphlets.empty
                and len(graphlets)
                == graphlets[
                    ["template", "seed", "graphlet_index"]
                ].drop_duplicates().shape[0]
            ),
        }
    )
    return rows


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=bool(args.quick))
    config["_runtime"] = {"quick": bool(args.quick)}

    bands, empirical_checks, rugged_templates = resolve_empirical_bands(
        config,
        project_root,
    )
    templates = load_template_catalog(config, project_root)
    write_json_strict(
        output_dir / "preflight.json",
        {
            "experiment_id": config["experiment_id"],
            "quick": bool(args.quick),
            "uses_exp050_transfer_function": False,
            "empirical_fitness_values_loaded": False,
            "descriptor_source": "RUN-135_frozen_operator_free_graphlet_metrics",
            "descriptor_implementation": (
                "scripts/ensemble_landscape_validation_stage4.py::metrics_from_cube"
            ),
            "observable_operating_point_override": (
                "experiment_local_beta_and_p0_no_public_preset_change"
            ),
            "tier1_template_count": len(templates),
            "rugged_end_templates": rugged_templates,
        },
    )
    write_csv_atomic(output_dir / "empirical_bands_resolved.csv", bands)
    write_csv_atomic(output_dir / "template_catalog.csv", templates)

    generated, supports, graphlets, audits, resume = run_generation(
        config,
        templates,
        output_dir=output_dir,
    )
    technical = _technical_rows(empirical_checks, audits, graphlets)
    result = coverage_and_decision(
        config,
        generated,
        bands,
        technical,
        quick=bool(args.quick),
    )

    write_csv_atomic(output_dir / "support_catalog.csv", supports)
    write_csv_atomic(output_dir / "support_aware_graphlet_audit.csv", graphlets)
    write_csv_atomic(output_dir / "generated_descriptor_metrics.csv", generated)
    write_csv_atomic(output_dir / "coverage_matrix.csv", result["coverage_rows"])
    write_csv_atomic(
        output_dir / "rugged_end_reachability.csv",
        result["reachability_rows"],
    )
    write_csv_atomic(output_dir / "specificity_table.csv", result["specificity_rows"])
    write_csv_atomic(
        output_dir / "selected_operating_region.csv",
        result["selected_rows"],
    )
    write_csv_atomic(output_dir / "decision_table.csv", result["decision_rows"])
    write_csv_atomic(output_dir / "technical_checks.csv", technical)
    write_json_strict(output_dir / "checkpoint_resume_summary.json", resume)
    write_json_strict(
        output_dir / "summary.json",
        {
            "experiment_id": config["experiment_id"],
            "quick": bool(args.quick),
            "decision": result["decision"],
            "gate_summary": result["gate_summary"],
            "n_templates": len(templates),
            "n_rugged_end_templates": len(rugged_templates),
            "n_generated_rows": len(generated),
            "n_coverage_rows": len(result["coverage_rows"]),
            "n_specificity_rows": len(result["specificity_rows"]),
            "uses_exp050_transfer_function": False,
            "stage_4c_4d_4e_hard_stop": "retained",
            "public_model_changed": False,
        },
    )


if __name__ == "__main__":
    main()
