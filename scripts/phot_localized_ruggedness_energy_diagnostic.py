#!/usr/bin/env python3
"""Localized ruggedness diagnostics for PHOT fitted/control landscapes.

This runner does not refit any parameters. It rebuilds the dense PHOT graphlet
context from the frozen RUN-077 source, loads the fitted ABC parameter sets from
RUN-091, and measures whether Dirichlet energy is locally clustered or diffuse.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_ROOT, ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import tem1_aggregate_graphlet_residual_ruggedness as graphlet_runner  # noqa: E402
from phot_coupled_family_discriminating_abc import (  # noqa: E402
    ASSAY_ID,
    build_graphlets,
)
from phot_coupled_landscape_grid_search import (  # noqa: E402
    apply_readout_model,
    feature_label,
    load_config,
    write_csv,
    write_json,
)
from project_adapt_env.proteingym_panel import merge_dict  # noqa: E402
from project_adapt_env.coupled_landscape import (  # noqa: E402
    coupled_values,
    plaquette_energy_by_pair,
    safe_corr,
)
from proteingym_system_specific_current_graphlet_refit import gnk_values  # noqa: E402


ABC_FAMILIES = {
    "no_module",
    "empirical_module",
    "random_module",
    "shuffled_module",
    "uncoupled_module",
}

STATIC_CANDIDATES = [
    {
        "candidate_family": "empirical_target",
        "candidate_model": "empirical_observed__latent_phi",
        "parameter_set_id": "latent_phi",
        "base_anchor": "",
        "module_source": "",
        "value_key": "empirical",
        "include_total": True,
        "include_innovation": False,
        "interpretation": "PHOT empirical latent_phi landscape",
    },
    {
        "candidate_family": "empirical_readout_residual",
        "candidate_model": "empirical_global_readout_residual",
        "parameter_set_id": "latent_phi__minus_poly_fitted_current_best",
        "base_anchor": "fitted_current_best",
        "module_source": "empirical_calibration_plaquette",
        "value_key": "empirical_readout_residual",
        "include_total": True,
        "include_innovation": False,
        "interpretation": "PHOT empirical residual after polynomial readout of fitted a(x)",
    },
    {
        "candidate_family": "fitted_current",
        "candidate_model": "fitted_current_best",
        "parameter_set_id": f"{ASSAY_ID}__fitted_current_best",
        "base_anchor": "fitted_current_best",
        "module_source": "baseline_control",
        "value_key": "fitted_current_best",
        "include_total": True,
        "include_innovation": False,
        "interpretation": "ABC-refit current core a(x)",
    },
    {
        "candidate_family": "fitted_current_posterior_mean",
        "candidate_model": "fitted_current_posterior_mean",
        "parameter_set_id": f"{ASSAY_ID}__fitted_current_posterior_mean",
        "base_anchor": "fitted_current_posterior_mean",
        "module_source": "baseline_control",
        "value_key": "fitted_current_posterior_mean",
        "include_total": True,
        "include_innovation": False,
        "interpretation": "ABC-refit current core posterior mean",
    },
    {
        "candidate_family": "smooth_core",
        "candidate_model": "fitted_smooth_core_only",
        "parameter_set_id": f"{ASSAY_ID}__fitted_smooth_core_only",
        "base_anchor": "fitted_smooth_core_only",
        "module_source": "baseline_control",
        "value_key": "fitted_smooth_core_only",
        "include_total": True,
        "include_innovation": False,
        "interpretation": "smooth core-only control",
    },
    {
        "candidate_family": "residual_only_diagnostic",
        "candidate_model": "fitted_residual_only_diagnostic",
        "parameter_set_id": f"{ASSAY_ID}__fitted_residual_only_diagnostic",
        "base_anchor": "fitted_smooth_core_only",
        "module_source": "baseline_control",
        "value_key": "fitted_residual_only_diagnostic",
        "include_total": True,
        "include_innovation": False,
        "interpretation": "current fitted a(x) minus smooth core diagnostic",
    },
    {
        "candidate_family": "gnk_local_null",
        "candidate_model": "gnk_local_null",
        "parameter_set_id": f"{ASSAY_ID}__gnk_k2",
        "base_anchor": "",
        "module_source": "gnk_control",
        "value_key": "gnk_local_null",
        "include_total": True,
        "include_innovation": False,
        "interpretation": "matched-mean/sd GNK local ruggedness null",
    },
]

KEY_LOCALIZATION_METRICS = [
    "edge_energy_adjacent_moran_i",
    "edge_energy_adjacent_pearson",
    "edge_energy_top10_fraction",
    "edge_energy_effective_fraction",
    "node_incident_top10_fraction",
    "node_incident_effective_fraction",
    "site_energy_top1_fraction",
    "site_energy_effective_fraction",
    "plaquette_energy_top1_fraction",
    "plaquette_energy_effective_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def resolve_path(path_value: str, project_root: Path) -> Path:
    path = Path(str(path_value))
    if not path.is_absolute():
        path = project_root / path
    return path


def load_diagnostic_config(path: Path, *, project_root: Path, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    base_config: dict[str, Any] = {}
    base_path_value = dict(raw.get("source_abc_run", {})).get("config_path")
    if base_path_value:
        base_path = resolve_path(str(base_path_value), project_root)
        base_config = load_config(base_path, quick=quick)
    quick_overrides = raw.pop("quick", {}) if quick else {}
    merged = merge_dict(base_config, raw)
    if quick_overrides:
        merged = merge_dict(merged, quick_overrides)
    return merged


def coerce_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def coerce_float(value: Any, default: float = math.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def parse_orders(value: Any) -> list[int]:
    if isinstance(value, list):
        return [int(item) for item in value]
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return [2, 3]
    text = text.strip("[]")
    text = text.replace(";", ",")
    return [int(float(part.strip())) for part in text.split(",") if part.strip()]


def finite_or_nan(value: float) -> float:
    return float(value) if np.isfinite(float(value)) else math.nan


def concentration_metrics(values: np.ndarray, *, prefix: str) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    arr = np.where(np.isfinite(arr), np.maximum(arr, 0.0), 0.0)
    total = float(np.sum(arr))
    n = int(arr.size)
    if n == 0 or total <= 1e-18:
        return {
            f"{prefix}_count": float(n),
            f"{prefix}_total": 0.0,
            f"{prefix}_mean": math.nan,
            f"{prefix}_sd": math.nan,
            f"{prefix}_hhi": math.nan,
            f"{prefix}_effective_count": math.nan,
            f"{prefix}_effective_fraction": math.nan,
            f"{prefix}_normalized_entropy": math.nan,
            f"{prefix}_gini": math.nan,
            f"{prefix}_top1_fraction": math.nan,
            f"{prefix}_top5_fraction": math.nan,
            f"{prefix}_top10_fraction": math.nan,
            f"{prefix}_top25_fraction": math.nan,
        }
    probs = arr / total
    sorted_probs = np.sort(probs)[::-1]
    hhi = float(np.sum(probs * probs))
    entropy = -float(np.sum(probs * np.log(np.maximum(probs, 1e-18))))
    sorted_arr = np.sort(arr)
    denom = n * float(np.sum(sorted_arr))
    gini = (
        float((2.0 * np.sum((np.arange(1, n + 1) * sorted_arr)) / denom) - ((n + 1.0) / n))
        if denom > 1e-18
        else math.nan
    )

    def top_fraction(frac: float) -> float:
        count = max(1, int(math.ceil(frac * n)))
        return float(np.sum(sorted_probs[:count]))

    return {
        f"{prefix}_count": float(n),
        f"{prefix}_total": total,
        f"{prefix}_mean": float(np.mean(arr)),
        f"{prefix}_sd": float(np.std(arr)),
        f"{prefix}_hhi": hhi,
        f"{prefix}_effective_count": float(1.0 / hhi) if hhi > 1e-18 else math.nan,
        f"{prefix}_effective_fraction": float((1.0 / hhi) / n) if hhi > 1e-18 else math.nan,
        f"{prefix}_normalized_entropy": float(entropy / math.log(n)) if n > 1 else 0.0,
        f"{prefix}_gini": gini,
        f"{prefix}_top1_fraction": top_fraction(1.0 / n),
        f"{prefix}_top5_fraction": top_fraction(0.05),
        f"{prefix}_top10_fraction": top_fraction(0.10),
        f"{prefix}_top25_fraction": top_fraction(0.25),
    }


def hypercube_edges(dim: int) -> list[dict[str, int]]:
    edges: list[dict[str, int]] = []
    for node in range(1 << int(dim)):
        for bit in range(int(dim)):
            neighbor = node ^ (1 << bit)
            if node < neighbor:
                edges.append(
                    {
                        "edge_index": len(edges),
                        "source_node": int(node),
                        "target_node": int(neighbor),
                        "bit_index": int(bit),
                    }
                )
    return edges


def adjacent_edge_pairs(edges: list[dict[str, int]]) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    endpoints = [
        {int(edge["source_node"]), int(edge["target_node"])}
        for edge in edges
    ]
    for i, j in itertools.combinations(range(len(edges)), 2):
        if endpoints[i] & endpoints[j]:
            pairs.append((i, j))
    return pairs


def edge_energy_table(values: np.ndarray, *, dim: int, features: tuple[Any, ...]) -> list[dict[str, Any]]:
    arr = np.asarray(values, dtype=np.float64)
    rows: list[dict[str, Any]] = []
    edges = hypercube_edges(dim)
    raw_energy: list[float] = []
    for edge in edges:
        delta = float(arr[edge["target_node"]] - arr[edge["source_node"]])
        raw_energy.append(delta * delta)
    total = float(np.sum(raw_energy))
    for edge, energy in zip(edges, raw_energy):
        bit = int(edge["bit_index"])
        rows.append(
            {
                **edge,
                "edge_id": f"{edge['source_node']:0{dim}b}-{edge['target_node']:0{dim}b}",
                "mutation_feature": feature_label(features[bit]),
                "delta": float(arr[edge["target_node"]] - arr[edge["source_node"]]),
                "edge_energy": float(energy),
                "edge_energy_fraction": float(energy / total) if total > 1e-18 else math.nan,
            }
        )
    return rows


def line_graph_autocorrelation(edge_energy: np.ndarray, edges: list[dict[str, int]]) -> dict[str, float]:
    pairs = adjacent_edge_pairs(edges)
    if not pairs:
        return {
            "edge_adjacency_pair_count": 0.0,
            "edge_energy_adjacent_pearson": math.nan,
            "edge_energy_adjacent_moran_i": math.nan,
        }
    arr = np.asarray(edge_energy, dtype=np.float64)
    left = np.asarray([arr[i] for i, _j in pairs], dtype=np.float64)
    right = np.asarray([arr[j] for _i, j in pairs], dtype=np.float64)
    mean = float(np.mean(arr))
    denom = float(np.sum((arr - mean) ** 2))
    moran = (
        float(len(arr) / len(pairs) * np.sum((left - mean) * (right - mean)) / denom)
        if denom > 1e-18
        else math.nan
    )
    return {
        "edge_adjacency_pair_count": float(len(pairs)),
        "edge_energy_adjacent_pearson": safe_corr(left, right),
        "edge_energy_adjacent_moran_i": moran,
    }


def ruggedness_metrics_for_values(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    edges = hypercube_edges(dim)
    edge_energy = np.asarray(
        [
            (float(arr[edge["target_node"]]) - float(arr[edge["source_node"]])) ** 2
            for edge in edges
        ],
        dtype=np.float64,
    )
    node_incident = np.zeros(1 << int(dim), dtype=np.float64)
    site_energy = np.zeros(int(dim), dtype=np.float64)
    for edge, energy in zip(edges, edge_energy):
        node_incident[edge["source_node"]] += float(energy)
        node_incident[edge["target_node"]] += float(energy)
        site_energy[edge["bit_index"]] += float(energy)
    plaquettes = plaquette_energy_by_pair(arr, dim=int(dim))
    pair_energy = np.asarray([float(value) for value in plaquettes.values()], dtype=np.float64)
    return {
        "fitness_mean": float(np.mean(arr)),
        "fitness_sd": float(np.std(arr)),
        "fitness_range": float(np.max(arr) - np.min(arr)),
        **line_graph_autocorrelation(edge_energy, edges),
        **concentration_metrics(edge_energy, prefix="edge_energy"),
        **concentration_metrics(node_incident, prefix="node_incident"),
        **concentration_metrics(site_energy, prefix="site_energy"),
        **concentration_metrics(pair_energy, prefix="plaquette_energy"),
    }


def module_enrichment_rows(
    *,
    values: np.ndarray,
    dim: int,
    features: tuple[Any, ...],
    modules_by_source: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    arr = np.asarray(values, dtype=np.float64)
    edges = hypercube_edges(dim)
    edge_energy = np.asarray(
        [
            (float(arr[edge["target_node"]]) - float(arr[edge["source_node"]])) ** 2
            for edge in edges
        ],
        dtype=np.float64,
    )
    edge_total = float(np.sum(edge_energy))
    plaquettes = plaquette_energy_by_pair(arr, dim=int(dim))
    plaquette_total = float(sum(float(value) for value in plaquettes.values()))
    index_by_feature = {feature: idx for idx, feature in enumerate(features)}
    rows: list[dict[str, Any]] = []
    for source in ["empirical_calibration_plaquette", "random_module", "shuffled_module", "uncoupled_empirical_module"]:
        covered_bits: set[int] = set()
        covered_pairs: set[tuple[int, int]] = set()
        modules = modules_by_source.get(source, [])
        active_modules = 0
        for module in modules:
            local = sorted(
                index_by_feature[feature]
                for feature in module.get("features", [])
                if feature in index_by_feature
            )
            if local:
                active_modules += 1
                covered_bits.update(local)
            if len(local) >= 2:
                covered_pairs.update(tuple(pair) for pair in itertools.combinations(local, 2))
        edge_covered = float(
            sum(float(energy) for edge, energy in zip(edges, edge_energy) if edge["bit_index"] in covered_bits)
        )
        pair_covered = float(
            sum(
                float(energy)
                for pair, energy in plaquettes.items()
                if tuple(sorted(pair)) in covered_pairs
            )
        )
        bit_coverage = float(len(covered_bits) / max(int(dim), 1))
        pair_count = max(len(plaquettes), 1)
        pair_coverage = float(len(covered_pairs) / pair_count)
        edge_fraction = float(edge_covered / edge_total) if edge_total > 1e-18 else math.nan
        pair_fraction = float(pair_covered / plaquette_total) if plaquette_total > 1e-18 else math.nan
        rows.append(
            {
                "module_source_tested": source,
                "n_modules_total": int(len(modules)),
                "n_modules_overlapping_graphlet": int(active_modules),
                "covered_site_count": int(len(covered_bits)),
                "covered_site_fraction": bit_coverage,
                "covered_pair_count": int(len(covered_pairs)),
                "covered_pair_fraction": pair_coverage,
                "edge_energy_in_module_sites": edge_covered,
                "edge_energy_module_site_fraction": edge_fraction,
                "edge_energy_module_site_enrichment": (
                    float(edge_fraction / bit_coverage)
                    if np.isfinite(edge_fraction) and bit_coverage > 1e-18
                    else math.nan
                ),
                "plaquette_energy_in_module_pairs": pair_covered,
                "plaquette_energy_module_pair_fraction": pair_fraction,
                "plaquette_energy_module_pair_enrichment": (
                    float(pair_fraction / pair_coverage)
                    if np.isfinite(pair_fraction) and pair_coverage > 1e-18
                    else math.nan
                ),
            }
        )
    return rows


def selected_parameter_ids(abc_output_dir: Path) -> set[str]:
    selected_path = abc_output_dir / "selected_family_distances.csv"
    if not selected_path.exists():
        return set()
    df = pd.read_csv(selected_path)
    if "graphlet_split" in df:
        df = df[df["graphlet_split"] == "calibration"]
    return {str(value) for value in df.get("parameter_set_id", []) if str(value) and str(value) != "nan"}


def selected_distance_rows(abc_output_dir: Path) -> list[dict[str, Any]]:
    path = abc_output_dir / "selected_family_distances.csv"
    if not path.exists():
        return []
    return pd.read_csv(path).to_dict(orient="records")


def load_evaluated_parameters(abc_output_dir: Path, *, config: dict[str, Any], quick: bool) -> list[dict[str, Any]]:
    path = abc_output_dir / "evaluated_parameter_sets.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing RUN-091 fitted parameter table: {path}")
    df = pd.read_csv(path)
    keep_families = [str(item) for item in config["diagnostic"]["abc_families"]]
    df = df[df["candidate_family"].astype(str).isin(keep_families)].copy()
    selected_ids = selected_parameter_ids(abc_output_dir)
    if bool(config["diagnostic"].get("selected_only", False)):
        df = df[df["parameter_set_id"].astype(str).isin(selected_ids)].copy()
    if quick:
        max_per_family = int(config["diagnostic"].get("quick_max_parameters_per_family", 2))
        df = df.sort_values(["candidate_family", "posterior_rank"]).groupby("candidate_family").head(max_per_family)

    rows: list[dict[str, Any]] = []
    for row in df.to_dict(orient="records"):
        params = {
            "candidate_model": str(row["candidate_model"]),
            "candidate_family": str(row["candidate_family"]),
            "parameter_set_id": str(row["parameter_set_id"]),
            "base_anchor": str(row["base_anchor"]),
            "module_source": str(row["module_source"]),
            "use_modules": coerce_bool(row.get("use_modules", True)),
            "uncoupled_auxiliary": coerce_bool(row.get("uncoupled_auxiliary", False)),
            "auxiliary_correlation": coerce_float(row.get("auxiliary_correlation"), 0.0),
            "operator_lambda": coerce_float(row.get("operator_lambda"), 0.0),
            "operator_gate_midpoint": coerce_float(row.get("operator_gate_midpoint"), 0.0),
            "operator_gate_steepness": coerce_float(row.get("operator_gate_steepness"), 1.0),
            "operator_gate_direction": str(row.get("operator_gate_direction", "high")),
            "module_residual_scale": coerce_float(row.get("module_residual_scale"), 0.0),
            "module_residual_density": coerce_float(row.get("module_residual_density"), 0.0),
            "module_residual_orders": parse_orders(row.get("module_residual_orders", "[2, 3]")),
            "posterior_rank": int(coerce_float(row.get("posterior_rank"), 0.0)),
            "posterior_objective_distance": coerce_float(row.get("posterior_objective_distance")),
            "is_calibration_selected": str(row["parameter_set_id"]) in selected_ids,
            "selection_role": (
                "calibration_selected"
                if str(row["parameter_set_id"]) in selected_ids
                else "posterior_predictive"
            ),
            "interpretation": "RUN-091 ABC fitted coupled-family candidate",
        }
        rows.append(params)
    return rows


def candidate_metadata(candidate: dict[str, Any]) -> dict[str, Any]:
    return {
        "candidate_family": candidate.get("candidate_family", ""),
        "candidate_model": candidate.get("candidate_model", ""),
        "parameter_set_id": candidate.get("parameter_set_id", ""),
        "base_anchor": candidate.get("base_anchor", ""),
        "module_source": candidate.get("module_source", ""),
        "posterior_rank": candidate.get("posterior_rank", ""),
        "posterior_objective_distance": candidate.get("posterior_objective_distance", ""),
        "is_calibration_selected": bool(candidate.get("is_calibration_selected", False)),
        "selection_role": candidate.get("selection_role", "static_control"),
        "interpretation": candidate.get("interpretation", ""),
    }


def add_diagnostic_rows(
    *,
    graphlet: dict[str, Any],
    candidate: dict[str, Any],
    target_kind: str,
    values: np.ndarray,
    rows: list[dict[str, Any]],
    edge_rows: list[dict[str, Any]],
    module_rows: list[dict[str, Any]],
    modules_by_source: dict[str, list[dict[str, Any]]],
) -> None:
    dim = int(graphlet["graphlet_dimension"])
    features = tuple(graphlet["features"])
    base = {
        "graphlet_id": str(graphlet["graphlet_id"]),
        "assay_id": ASSAY_ID,
        "graphlet_split": str(graphlet["graphlet_split"]),
        "graphlet_dimension": dim,
        "mutation_features": str(graphlet["mutation_features"]),
        "target_kind": target_kind,
        **candidate_metadata(candidate),
    }
    metric = ruggedness_metrics_for_values(np.asarray(values, dtype=np.float64), dim=dim)
    rows.append({**base, **metric})
    for row in edge_energy_table(np.asarray(values, dtype=np.float64), dim=dim, features=features):
        edge_rows.append({**base, **row})
    for row in module_enrichment_rows(
        values=np.asarray(values, dtype=np.float64),
        dim=dim,
        features=features,
        modules_by_source=modules_by_source,
    ):
        module_rows.append({**base, **row})


def graphlet_value_for_static(
    *,
    graphlet: dict[str, Any],
    candidate: dict[str, Any],
    primary_score: str,
    empirical_innovation: dict[str, np.ndarray],
) -> np.ndarray:
    key = str(candidate["value_key"])
    if key == "empirical":
        return np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
    if key == "empirical_readout_residual":
        return np.asarray(empirical_innovation[str(graphlet["graphlet_id"])], dtype=np.float64)
    if key == "gnk_local_null":
        empirical_values = np.asarray(graphlet["empirical_values"][primary_score], dtype=np.float64)
        return gnk_values(
            dim=int(graphlet["graphlet_dimension"]),
            seed=graphlet_runner.stable_seed("phot-localized-ruggedness-gnk", graphlet["graphlet_id"]),
            target_mean=float(np.mean(empirical_values)),
            target_sd=float(np.std(empirical_values)),
            k=2,
        )
    return np.asarray(graphlet[key], dtype=np.float64)


def aggregate_metric_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    group_cols = [
        "target_kind",
        "graphlet_split",
        "candidate_family",
        "candidate_model",
        "parameter_set_id",
        "selection_role",
        "is_calibration_selected",
    ]
    numeric_cols = [
        col
        for col in df.columns
        if col not in group_cols
        and col not in {"graphlet_id", "assay_id", "mutation_features", "interpretation", "base_anchor", "module_source"}
        and pd.api.types.is_numeric_dtype(df[col])
    ]
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {col: value for col, value in zip(group_cols, keys)}
        row["n_graphlets"] = int(len(sub))
        first = sub.iloc[0].to_dict()
        row["base_anchor"] = first.get("base_anchor", "")
        row["module_source"] = first.get("module_source", "")
        row["posterior_rank"] = first.get("posterior_rank", "")
        row["posterior_objective_distance"] = first.get("posterior_objective_distance", "")
        row["interpretation"] = first.get("interpretation", "")
        for col in numeric_cols:
            values = pd.to_numeric(sub[col], errors="coerce").dropna()
            if values.empty:
                continue
            row[f"{col}_mean"] = float(values.mean())
            row[f"{col}_median"] = float(values.median())
        out.append(row)
    return out


def aggregate_module_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    group_cols = [
        "target_kind",
        "graphlet_split",
        "candidate_family",
        "candidate_model",
        "parameter_set_id",
        "selection_role",
        "is_calibration_selected",
        "module_source_tested",
    ]
    numeric_cols = [
        col
        for col in [
            "edge_energy_module_site_fraction",
            "edge_energy_module_site_enrichment",
            "plaquette_energy_module_pair_fraction",
            "plaquette_energy_module_pair_enrichment",
            "covered_site_fraction",
            "covered_pair_fraction",
            "n_modules_overlapping_graphlet",
        ]
        if col in df.columns
    ]
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {col: value for col, value in zip(group_cols, keys)}
        row["n_graphlets"] = int(len(sub))
        first = sub.iloc[0].to_dict()
        row["base_anchor"] = first.get("base_anchor", "")
        row["module_source"] = first.get("module_source", "")
        row["posterior_rank"] = first.get("posterior_rank", "")
        for col in numeric_cols:
            values = pd.to_numeric(sub[col], errors="coerce").dropna()
            if values.empty:
                continue
            row[f"{col}_mean"] = float(values.mean())
            row[f"{col}_median"] = float(values.median())
        out.append(row)
    return out


def selected_vs_empirical(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary = pd.DataFrame(aggregate_metric_rows(rows))
    if summary.empty:
        return []
    selected = summary[
        (summary["selection_role"].isin(["calibration_selected", "static_control"]))
        & (summary["target_kind"] == "total_landscape")
    ].copy()
    out: list[dict[str, Any]] = []
    for split in sorted(selected["graphlet_split"].dropna().unique()):
        ref = selected[
            (selected["graphlet_split"] == split)
            & (selected["candidate_family"] == "empirical_target")
        ]
        if ref.empty:
            continue
        ref_row = ref.iloc[0]
        for _, row in selected[selected["graphlet_split"] == split].iterrows():
            if row["candidate_family"] == "empirical_target":
                continue
            base = {
                "graphlet_split": split,
                "target_kind": "total_landscape",
                "candidate_family": row["candidate_family"],
                "candidate_model": row["candidate_model"],
                "parameter_set_id": row["parameter_set_id"],
                "selection_role": row["selection_role"],
            }
            for metric in KEY_LOCALIZATION_METRICS:
                col = f"{metric}_median"
                if col not in row or col not in ref_row:
                    continue
                base[f"{metric}_candidate_median"] = finite_or_nan(row[col])
                base[f"{metric}_empirical_median"] = finite_or_nan(ref_row[col])
                base[f"{metric}_delta"] = finite_or_nan(row[col] - ref_row[col])
            out.append(base)
    return out


def metric_rank_rows(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    df = pd.DataFrame(summary_rows)
    if df.empty:
        return []
    df = df[df["selection_role"].isin(["calibration_selected", "static_control"])].copy()
    directions = {
        "edge_energy_adjacent_moran_i_median": "descending",
        "edge_energy_adjacent_pearson_median": "descending",
        "edge_energy_top10_fraction_median": "descending",
        "edge_energy_effective_fraction_median": "ascending",
        "node_incident_top10_fraction_median": "descending",
        "node_incident_effective_fraction_median": "ascending",
        "site_energy_top1_fraction_median": "descending",
        "site_energy_effective_fraction_median": "ascending",
        "plaquette_energy_top1_fraction_median": "descending",
        "plaquette_energy_effective_fraction_median": "ascending",
    }
    rows: list[dict[str, Any]] = []
    for (target_kind, split), sub in df.groupby(["target_kind", "graphlet_split"], dropna=False):
        for metric, direction in directions.items():
            if metric not in sub.columns:
                continue
            ranked = sub[np.isfinite(pd.to_numeric(sub[metric], errors="coerce"))].copy()
            if ranked.empty:
                continue
            ranked[metric] = pd.to_numeric(ranked[metric], errors="coerce")
            ranked = ranked.sort_values(metric, ascending=(direction == "ascending"))
            for rank, row in enumerate(ranked.to_dict(orient="records"), start=1):
                rows.append(
                    {
                        "target_kind": target_kind,
                        "graphlet_split": split,
                        "metric": metric.replace("_median", ""),
                        "localized_direction": direction,
                        "rank": int(rank),
                        "candidate_family": row["candidate_family"],
                        "candidate_model": row["candidate_model"],
                        "parameter_set_id": row["parameter_set_id"],
                        "selection_role": row["selection_role"],
                        "value": float(row[metric]),
                    }
                )
    return rows


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_diagnostic_config(Path(args.config), project_root=project_root, quick=args.quick)

    source_abc_dir = resolve_path(config["source_abc_run"]["output_dir"], project_root)
    if not source_abc_dir.exists():
        raise FileNotFoundError(
            f"RUN-091 ABC output directory is missing: {source_abc_dir}. Restore it with DVC before running."
        )

    graphlets, context, modules_by_source, _readout = build_graphlets(
        config=config,
        project_root=project_root,
        output_dir=output_dir,
    )
    if args.quick:
        max_graphlets = int(config["diagnostic"].get("quick_max_graphlets", 4))
        graphlets = graphlets[:max_graphlets]

    primary_score = str(config["graphlet_evaluation"]["primary_score_column"])
    empirical_innovation = context["empirical_innovation"]
    abc_candidates = load_evaluated_parameters(source_abc_dir, config=config, quick=args.quick)
    selected_ids = selected_parameter_ids(source_abc_dir)

    static_candidates: list[dict[str, Any]] = []
    for row in STATIC_CANDIDATES:
        candidate = dict(row)
        candidate["is_calibration_selected"] = str(candidate["parameter_set_id"]) in selected_ids
        candidate["selection_role"] = "static_control"
        static_candidates.append(candidate)

    graphlet_rows: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    module_rows: list[dict[str, Any]] = []

    for graphlet in graphlets:
        for candidate in static_candidates:
            values = graphlet_value_for_static(
                graphlet=graphlet,
                candidate=candidate,
                primary_score=primary_score,
                empirical_innovation=empirical_innovation,
            )
            add_diagnostic_rows(
                graphlet=graphlet,
                candidate=candidate,
                target_kind="total_landscape",
                values=values,
                rows=graphlet_rows,
                edge_rows=edge_rows,
                module_rows=module_rows,
                modules_by_source=modules_by_source,
            )

        for candidate in abc_candidates:
            base_values = np.asarray(graphlet[str(candidate["base_anchor"])], dtype=np.float64)
            values, innovation, active, coeff_energy = coupled_values(
                base_values=base_values,
                graphlet_features=tuple(graphlet["features"]),
                modules=modules_by_source[str(candidate["module_source"])],
                parameter_set=candidate,
                seed=graphlet_runner.stable_seed(candidate["parameter_set_id"], graphlet["graphlet_id"]),
            )
            enriched_candidate = {
                **candidate,
                **active,
                "module_coefficient_energy_total": float(sum(coeff_energy.values())),
            }
            add_diagnostic_rows(
                graphlet=graphlet,
                candidate=enriched_candidate,
                target_kind="total_landscape",
                values=values,
                rows=graphlet_rows,
                edge_rows=edge_rows,
                module_rows=module_rows,
                modules_by_source=modules_by_source,
            )
            add_diagnostic_rows(
                graphlet=graphlet,
                candidate=enriched_candidate,
                target_kind="innovation_from_base",
                values=innovation,
                rows=graphlet_rows,
                edge_rows=edge_rows,
                module_rows=module_rows,
                modules_by_source=modules_by_source,
            )

    model_summary_rows = aggregate_metric_rows(graphlet_rows)
    module_summary_rows = aggregate_module_rows(module_rows)
    delta_rows = selected_vs_empirical(graphlet_rows)
    rank_rows = metric_rank_rows(model_summary_rows)

    write_csv(output_dir / "localized_ruggedness_graphlet_metrics.csv", graphlet_rows)
    write_csv(output_dir / "edge_dirichlet_energy_rows.csv", edge_rows)
    write_csv(output_dir / "module_enrichment_by_graphlet.csv", module_rows)
    write_csv(output_dir / "model_localization_summary.csv", model_summary_rows)
    write_csv(output_dir / "module_enrichment_summary.csv", module_summary_rows)
    write_csv(output_dir / "selected_vs_empirical_localization_delta.csv", delta_rows)
    write_csv(output_dir / "localization_metric_ranks.csv", rank_rows)
    write_csv(output_dir / "run091_selected_family_distances.csv", selected_distance_rows(source_abc_dir))
    write_csv(output_dir / "run091_reused_parameter_sets.csv", abc_candidates)

    preflight = {
        "assay_id": ASSAY_ID,
        "experiment_id": config["experiment"]["id"],
        "mode": "localized_dirichlet_energy_diagnostic_reusing_run091_abc_particles",
        "source_run": config["source_run"],
        "source_abc_run": config["source_abc_run"],
        "primary_score_column": primary_score,
        "n_graphlets": int(len(graphlets)),
        "n_calibration_graphlets": int(sum(g["graphlet_split"] == "calibration" for g in graphlets)),
        "n_test_graphlets": int(sum(g["graphlet_split"] == "test" for g in graphlets)),
        "n_static_candidates": int(len(static_candidates)),
        "n_abc_parameter_sets_reused": int(len(abc_candidates)),
        "selected_parameter_ids": sorted(selected_ids),
        "module_counts": {source: int(len(modules)) for source, modules in modules_by_source.items()},
        "diagnostic_metrics": {
            "edge_autocorrelation": [
                "edge_energy_adjacent_pearson",
                "edge_energy_adjacent_moran_i",
            ],
            "energy_concentration": [
                "edge_energy_effective_fraction",
                "edge_energy_top10_fraction",
                "node_incident_effective_fraction",
                "site_energy_effective_fraction",
                "plaquette_energy_effective_fraction",
            ],
            "module_enrichment": [
                "edge_energy_module_site_enrichment",
                "plaquette_energy_module_pair_enrichment",
            ],
        },
    }
    write_json(output_dir / "preflight.json", preflight)
    write_json(
        output_dir / "summary.json",
        {
            "experiment_id": config["experiment"]["id"],
            "preflight": preflight,
            "output_files": {
                "preflight": str(output_dir / "preflight.json"),
                "localized_ruggedness_graphlet_metrics": str(output_dir / "localized_ruggedness_graphlet_metrics.csv"),
                "edge_dirichlet_energy_rows": str(output_dir / "edge_dirichlet_energy_rows.csv"),
                "module_enrichment_by_graphlet": str(output_dir / "module_enrichment_by_graphlet.csv"),
                "model_localization_summary": str(output_dir / "model_localization_summary.csv"),
                "module_enrichment_summary": str(output_dir / "module_enrichment_summary.csv"),
                "selected_vs_empirical_localization_delta": str(output_dir / "selected_vs_empirical_localization_delta.csv"),
                "localization_metric_ranks": str(output_dir / "localization_metric_ranks.csv"),
                "run091_selected_family_distances": str(output_dir / "run091_selected_family_distances.csv"),
                "run091_reused_parameter_sets": str(output_dir / "run091_reused_parameter_sets.csv"),
            },
            "selection_policy": "Primary comparisons use calibration-selected RUN-091 fitted/control candidates; posterior rows are retained to show fit-particle spread.",
        },
    )


if __name__ == "__main__":
    main()
