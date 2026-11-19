#!/usr/bin/env python3
"""PHOT P0 target and measurement-operator isolation diagnostic."""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import random
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_ROOT, ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import tem1_aggregate_graphlet_residual_ruggedness as graphlet_runner  # noqa: E402
from phot_coupled_landscape_grid_search import load_config, write_csv, write_json  # noqa: E402
from phot_nested_fit_protocol_benchmark import epistasis_property_metrics  # noqa: E402
from phot_trait_space_bounded_gate_fast_fit import (  # noqa: E402
    sparse_walsh_config,
    sparse_walsh_order_weights,
    sparse_walsh_orders,
)
from project_adapt_env.dense_descriptive import (  # noqa: E402
    CompleteGraphlet,
    build_score_map,
    enumerate_complete_graphlets,
    exact_summary,
    values_for_graphlet,
)
from project_adapt_env.proteingym_panel import fit_mavenn_model  # noqa: E402


ASSAY_ID = "PHOT_CHLRE_Chen_2023"

P0_1_TARGET_METRICS = [
    "walsh_high_order_fraction",
    "n_local_maxima",
    "accessible_path_fraction",
    "edge_energy_top10_fraction",
    "edge_energy_effective_fraction",
    "site_energy_top1_fraction",
    "site_energy_effective_fraction",
]

P0_2_METRICS = [
    *P0_1_TARGET_METRICS,
    "magnitude_epistasis_to_single_mutant_sd_ratio",
    "single_mutant_effect_sd",
    "sign_epistasis_fraction_local",
    "reciprocal_sign_epistasis_fraction_local",
    "diminishing_returns_epistasis_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def resolve_project_path(project_root: Path, raw: str | Path) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        path = project_root / path
    return path


def source_output_dir(config: dict[str, Any], project_root: Path) -> Path:
    path = resolve_project_path(project_root, str(config["source_run"]["output_dir"]))
    if not path.exists():
        raise FileNotFoundError(f"source RUN output directory is missing: {path}")
    return path


def load_source_payload(config: dict[str, Any], project_root: Path) -> dict[str, Any]:
    source_dir = source_output_dir(config, project_root)
    system_dir = source_dir / "systems" / ASSAY_ID
    fit_dir = system_dir / "fit"
    assay_dir = fit_dir / "checkpoints" / "panel_preparation" / ASSAY_ID
    assay_path = assay_dir / "prepared_assay.csv"
    metrics_path = assay_dir / "mavenn_metrics.json"
    panel_path = fit_dir / "selected_panel.csv"
    if not assay_path.is_file() or not metrics_path.is_file() or not panel_path.is_file():
        raise FileNotFoundError(f"missing PHOT prepared assay payload under {fit_dir}")
    assay = pd.read_csv(assay_path)
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    panel = pd.read_csv(panel_path)
    wildtype = str(panel.loc[panel["DMS_id"] == ASSAY_ID, "target_seq"].iloc[0])
    return {
        "source_dir": source_dir,
        "system_dir": system_dir,
        "fit_dir": fit_dir,
        "assay": assay,
        "mavenn_metrics": metrics,
        "panel": panel,
        "wildtype": wildtype,
        "wildtype_phi": float(metrics["wildtype_phi"]),
        "wildtype_yhat": float(metrics["wildtype_yhat"]),
    }


def graphlets_from_assay(
    *,
    assay: pd.DataFrame,
    wildtype: str,
    wt_scores: dict[str, float],
    score_columns: list[str],
    config: dict[str, Any],
) -> list[CompleteGraphlet]:
    cfg = dict(config["phot_graphlets"])
    score_map = build_score_map(
        assay,
        wildtype=wildtype,
        score_columns=score_columns,
        wt_scores=wt_scores,
    )
    return enumerate_complete_graphlets(
        score_map,
        assay_id=ASSAY_ID,
        score_columns=score_columns,
        dimensions=[int(value) for value in cfg["dimensions"]],
        max_graphlets_per_dimension=int(cfg["max_graphlets_per_dimension"]),
        top_feature_pool=int(cfg["top_feature_pool"]),
        max_combinations_per_dimension=cfg.get("max_combinations_per_dimension"),
    )


def graphlet_catalog_rows(graphlets: list[CompleteGraphlet]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, graphlet in enumerate(graphlets):
        split = "calibration" if index % 2 == 0 else "test"
        rows.append(
            {
                "assay_id": ASSAY_ID,
                "graphlet_id": graphlet.graphlet_id,
                "graphlet_index": int(index),
                "graphlet_split": split,
                "graphlet_dimension": int(graphlet.dimension),
                "mutation_features": graphlet.mutation_features,
            }
        )
    return rows


def split_by_graphlet_id(graphlets: list[CompleteGraphlet]) -> dict[str, str]:
    return {graphlet.graphlet_id: ("calibration" if index % 2 == 0 else "test") for index, graphlet in enumerate(graphlets)}


def metric_bundle(values: np.ndarray, *, dim: int) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    metrics = exact_summary(arr, dim=int(dim))
    metrics.update(epistasis_property_metrics(arr, dim=int(dim)))
    return {
        key: float(value)
        for key, value in metrics.items()
        if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(float(value))
    }


def metric_rows_for_score(
    *,
    assay: pd.DataFrame,
    wildtype: str,
    wt_score: float,
    score_column: str,
    graphlets: list[CompleteGraphlet],
    candidate_model: str,
    parameter_set_id: str,
    phase: str,
    extra: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    score_map = build_score_map(
        assay,
        wildtype=wildtype,
        score_columns=[score_column],
        wt_scores={score_column: float(wt_score)},
    )
    split_map = split_by_graphlet_id(graphlets)
    rows: list[dict[str, Any]] = []
    for graphlet in graphlets:
        values = values_for_graphlet(score_map, graphlet, score_column=score_column)
        rows.append(
            {
                "phase": phase,
                "assay_id": ASSAY_ID,
                "graphlet_id": graphlet.graphlet_id,
                "graphlet_split": split_map[graphlet.graphlet_id],
                "graphlet_dimension": int(graphlet.dimension),
                "mutation_features": graphlet.mutation_features,
                "candidate_model": candidate_model,
                "parameter_set_id": parameter_set_id,
                "score_column": score_column,
                **(extra or {}),
                **metric_bundle(values, dim=int(graphlet.dimension)),
            }
        )
    return rows


def bootstrap_metric_bands(
    rows: list[dict[str, Any]],
    *,
    metrics: list[str],
    iterations: int,
    seed: int,
) -> list[dict[str, Any]]:
    df = pd.DataFrame(rows)
    if df.empty:
        return []
    rng = np.random.RandomState(int(seed))
    group_cols = ["phase", "candidate_model", "parameter_set_id", "score_column"]
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        group = {col: value for col, value in zip(group_cols, keys)}
        n = int(len(sub))
        if n == 0:
            continue
        for metric in metrics:
            if metric not in sub:
                continue
            values = pd.to_numeric(sub[metric], errors="coerce").dropna().to_numpy(dtype=np.float64)
            if values.size == 0:
                continue
            means = []
            for _ in range(max(int(iterations), 1)):
                sample = rng.choice(values, size=values.size, replace=True)
                means.append(float(np.mean(sample)))
            boot = np.asarray(means, dtype=np.float64)
            out.append(
                {
                    **group,
                    "metric": metric,
                    "n_graphlets": n,
                    "graphlet_mean": float(np.mean(values)),
                    "graphlet_sd": float(np.std(values, ddof=0)),
                    "bootstrap_mean_q05": float(np.quantile(boot, 0.05)),
                    "bootstrap_mean_q50": float(np.quantile(boot, 0.50)),
                    "bootstrap_mean_q95": float(np.quantile(boot, 0.95)),
                }
            )
    return out


def set_training_seed(seed: int) -> None:
    random.seed(int(seed))
    np.random.seed(int(seed) % (2**32))
    try:
        import tensorflow as tf  # type: ignore

        tf.random.set_seed(int(seed))
    except ModuleNotFoundError:
        pass


def mavenn_available() -> bool:
    return importlib.util.find_spec("mavenn") is not None


def quick_surrogate_phi(values: np.ndarray, *, seed: int, profile_index: int) -> np.ndarray:
    rng = np.random.RandomState(int(seed) + 1009 * int(profile_index))
    arr = np.asarray(values, dtype=np.float64)
    centered = arr - float(np.mean(arr))
    sd = max(float(np.std(centered)), 1e-8)
    nonlinear = np.tanh(centered / sd)
    return centered / sd + 0.15 * nonlinear + rng.normal(0.0, 0.02, size=arr.shape[0])


def fit_or_surrogate_mavenn(
    *,
    assay: pd.DataFrame,
    wildtype: str,
    cfg: dict[str, Any],
    seed: int,
    allow_quick_surrogate: bool,
    profile_index: int,
) -> tuple[np.ndarray, float, dict[str, Any]]:
    if allow_quick_surrogate and bool(cfg.get("force_quick_surrogate", False)):
        phi = quick_surrogate_phi(assay["DMS_score"].to_numpy(dtype=np.float64), seed=seed, profile_index=profile_index)
        return phi, 0.0, {
            "fit_variant_count": int(len(assay)),
            "train_variant_count": 0,
            "test_variant_count": 0,
            "wildtype_phi": 0.0,
            "wildtype_yhat": 0.0,
            "quick_surrogate": True,
            "quick_surrogate_forced": True,
        }
    if mavenn_available():
        set_training_seed(seed)
        model, metrics = fit_mavenn_model(
            assay_frame=assay,
            wildtype_sequence=wildtype,
            cfg=cfg,
        )
        phi = np.asarray(model.x_to_phi(assay["mutated_sequence"].to_numpy(dtype=object)), dtype=np.float64).reshape(-1)
        wt_phi = float(np.asarray(model.x_to_phi(np.asarray([wildtype], dtype=object))).reshape(-1)[0])
        return phi, wt_phi, {**metrics, "quick_surrogate": False}
    if not allow_quick_surrogate:
        raise RuntimeError(
            "mavenn is required for full P0.1/P0.2 measurement-operator tests; "
            "install the project overlay or run through labproj."
        )
    phi = quick_surrogate_phi(assay["DMS_score"].to_numpy(dtype=np.float64), seed=seed, profile_index=profile_index)
    return phi, 0.0, {
        "fit_variant_count": int(len(assay)),
        "train_variant_count": 0,
        "test_variant_count": 0,
        "wildtype_phi": 0.0,
        "wildtype_yhat": 0.0,
        "quick_surrogate": True,
    }


def run_p0_1(
    *,
    config: dict[str, Any],
    payload: dict[str, Any],
    graphlets: list[CompleteGraphlet],
    quick: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    assay = payload["assay"].copy()
    wildtype = str(payload["wildtype"])
    metric_rows = metric_rows_for_score(
        assay=assay,
        wildtype=wildtype,
        wt_score=float(payload["wildtype_phi"]),
        score_column="latent_phi",
        graphlets=graphlets,
        candidate_model="empirical_cached_mavenn_latent_phi",
        parameter_set_id="RUN-077_cached_latent_phi",
        phase="P0.1",
        extra={"inference_profile": "cached", "refit_seed": math.nan, "quick_surrogate": False},
    )
    refit_metric_rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    refit_cfg = dict(config["p0_1"]["mavenn_refits"])
    seeds = [int(seed) for seed in refit_cfg["seeds"]]
    profiles = [dict(profile) for profile in refit_cfg["profiles"]]
    allow_surrogate = bool(quick and refit_cfg.get("quick_surrogate_without_mavenn", True))
    for profile_index, profile in enumerate(profiles, start=1):
        profile_name = str(profile.pop("name"))
        for seed in seeds:
            cfg = copy.deepcopy(config["mavenn"])
            cfg.update(profile)
            cfg["sample_seed"] = int(seed)
            phi, wt_phi, metrics = fit_or_surrogate_mavenn(
                assay=assay,
                wildtype=wildtype,
                cfg=cfg,
                seed=seed,
                allow_quick_surrogate=allow_surrogate,
                profile_index=profile_index,
            )
            score_column = f"latent_phi__{profile_name}__seed{seed}"
            assay[score_column] = phi
            parameter_set_id = f"{profile_name}__seed{seed}"
            rows = metric_rows_for_score(
                assay=assay,
                wildtype=wildtype,
                wt_score=wt_phi,
                score_column=score_column,
                graphlets=graphlets,
                candidate_model="empirical_refit_mavenn_latent_phi",
                parameter_set_id=parameter_set_id,
                phase="P0.1",
                extra={
                    "inference_profile": profile_name,
                    "refit_seed": int(seed),
                    "quick_surrogate": bool(metrics.get("quick_surrogate", False)),
                },
            )
            refit_metric_rows.extend(rows)
            fit_rows.append(
                {
                    "phase": "P0.1",
                    "inference_profile": profile_name,
                    "refit_seed": int(seed),
                    "parameter_set_id": parameter_set_id,
                    "gpmap_type": cfg.get("gpmap_type"),
                    "theta_regularization": cfg.get("theta_regularization"),
                    "eta_regularization": cfg.get("eta_regularization"),
                    "epochs": cfg.get("epochs"),
                    "batch_size": cfg.get("batch_size"),
                    **{f"mavenn_{key}": value for key, value in metrics.items() if not isinstance(value, dict)},
                    **{
                        f"mavenn_test_{key}": value
                        for key, value in dict(metrics.get("test_metrics", {})).items()
                    },
                    **{
                        f"mavenn_train_{key}": value
                        for key, value in dict(metrics.get("train_metrics", {})).items()
                    },
                }
            )
    all_rows = metric_rows + refit_metric_rows
    band_rows = bootstrap_metric_bands(
        all_rows,
        metrics=P0_1_TARGET_METRICS,
        iterations=int(config["p0_1"].get("bootstrap_iterations", 500)),
        seed=int(config["p0_1"].get("bootstrap_seed", 20260618)),
    )
    replicate_rows = replicate_audit(assay=payload["assay"], config=config)
    return all_rows, band_rows, fit_rows, replicate_rows


def replicate_audit(*, assay: pd.DataFrame, config: dict[str, Any]) -> list[dict[str, Any]]:
    configured = [str(col) for col in config["p0_1"].get("replicate_columns", []) if str(col) in assay.columns]
    detected = [
        str(col)
        for col in assay.columns
        if "replicate" in str(col).lower() or str(col).lower().startswith("dms_score_rep")
    ]
    columns = sorted(set(configured + detected))
    return [
        {
            "phase": "P0.1",
            "assay_id": ASSAY_ID,
            "has_assay_replicate_columns": bool(columns),
            "replicate_columns": ";".join(columns),
            "n_replicate_columns": int(len(columns)),
            "action": "propagate_if_present",
            "note": (
                "No replicate-valued columns were present in the prepared PHOT assay."
                if not columns
                else "Replicate columns detected; downstream ANA should inspect replicate-specific target bands."
            ),
        }
    ]


def mutation_mask_for_sequence(wildtype: str, sequence: str) -> frozenset[tuple[int, str]]:
    return frozenset(
        (idx, aa)
        for idx, (wt, aa) in enumerate(zip(str(wildtype), str(sequence)))
        if aa != wt
    )


def stable_uint64(*parts: Any) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little")


def synthetic_terms(
    *,
    feature_pool: list[tuple[int, str]],
    profile: dict[str, Any],
) -> list[dict[str, Any]]:
    seed = int(profile.get("seed", 20260618))
    rng = np.random.RandomState(seed)
    order_weights = {int(k): float(v) for k, v in dict(profile["order_weights"]).items()}
    n_terms = int(profile.get("n_terms", 96))
    total_weight = max(float(sum(max(v, 0.0) for v in order_weights.values())), 1e-12)
    expected = {order: n_terms * max(weight, 0.0) / total_weight for order, weight in order_weights.items()}
    counts = {order: int(math.floor(value)) for order, value in expected.items()}
    remaining = n_terms - sum(counts.values())
    ranked = sorted(expected, key=lambda order: (-(expected[order] - counts[order]), order))
    for order in ranked[: max(remaining, 0)]:
        counts[order] += 1
    terms: list[dict[str, Any]] = []
    for order, count in sorted(counts.items()):
        if order <= 0 or count <= 0:
            continue
        coeff_sd = math.sqrt(max(order_weights[order], 0.0) / max(count, 1))
        for term_index in range(count):
            if len(feature_pool) < order:
                support = tuple(feature_pool)
            else:
                indices = rng.choice(np.arange(len(feature_pool)), size=order, replace=False)
                support = tuple(sorted((feature_pool[int(idx)] for idx in indices), key=lambda x: (x[0], x[1])))
            terms.append(
                {
                    "profile": str(profile["name"]),
                    "term_id": f"{profile['name']}__order{order}__{term_index:04d}",
                    "order": int(order),
                    "coefficient": float(rng.normal(0.0, coeff_sd)),
                    "support": support,
                    "support_features": ";".join(f"{pos}:{aa}" for pos, aa in support),
                }
            )
    return terms


def raw_term_values_for_masks(
    masks: Iterable[frozenset[tuple[int, str]]],
    *,
    terms: list[dict[str, Any]],
    additive_weight: float,
) -> np.ndarray:
    values: list[float] = []
    for mask in masks:
        mask_set = set(mask)
        value = additive_weight * float(len(mask_set))
        for term in terms:
            product = 1.0
            for feature in term["support"]:
                product *= 1.0 if feature in mask_set else -1.0
            value += float(term["coefficient"]) * product
        values.append(value)
    return np.asarray(values, dtype=np.float64)


def normalize_assay_scale(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    sd = float(np.std(arr))
    return (arr - float(np.mean(arr))) / max(sd, 1e-12)


def run_p0_2(
    *,
    config: dict[str, Any],
    payload: dict[str, Any],
    graphlets: list[CompleteGraphlet],
    quick: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    assay_base = payload["assay"].copy()
    wildtype = str(payload["wildtype"])
    feature_pool = sorted(
        {feature for graphlet in graphlets for feature in graphlet.features},
        key=lambda feature: (int(feature[0]), str(feature[1])),
    )
    all_terms: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    transfer_rows: list[dict[str, Any]] = []
    mavenn_rows: list[dict[str, Any]] = []
    masks = [mutation_mask_for_sequence(wildtype, seq) for seq in assay_base["mutated_sequence"].astype(str)]
    wt_mask = frozenset()
    allow_surrogate = bool(quick and config["p0_2"].get("quick_surrogate_without_mavenn", True))
    for profile_index, profile in enumerate(config["p0_2"]["known_spectrum_profiles"], start=1):
        profile = dict(profile)
        profile_name = str(profile["name"])
        terms = synthetic_terms(feature_pool=feature_pool, profile=profile)
        all_terms.extend(
            {
                key: (";".join(f"{pos}:{aa}" for pos, aa in value) if key == "support" else value)
                for key, value in term.items()
            }
            for term in terms
        )
        raw_score = raw_term_values_for_masks(
            masks,
            terms=terms,
            additive_weight=float(profile.get("additive_weight", 0.0)),
        )
        synthetic_score = normalize_assay_scale(raw_score)
        raw_mean = float(np.mean(raw_score))
        raw_sd = max(float(np.std(raw_score)), 1e-12)
        wt_raw = float(
            raw_term_values_for_masks(
                [wt_mask],
                terms=terms,
                additive_weight=float(profile.get("additive_weight", 0.0)),
            )[0]
        )
        wt_score = float((wt_raw - raw_mean) / raw_sd)
        assay = assay_base.copy()
        assay["synthetic_dms_score"] = synthetic_score
        assay["DMS_score"] = synthetic_score
        injected_rows = metric_rows_for_score(
            assay=assay,
            wildtype=wildtype,
            wt_score=wt_score,
            score_column="synthetic_dms_score",
            graphlets=graphlets,
            candidate_model="known_injected_synthetic_landscape",
            parameter_set_id=profile_name,
            phase="P0.2",
            extra={
                "spectrum_profile": profile_name,
                "measurement_stage": "injected_direct",
                "quick_surrogate": False,
            },
        )
        metric_rows.extend(injected_rows)

        cfg = copy.deepcopy(config["mavenn"])
        cfg.update(dict(config["p0_2"].get("mavenn_overrides", {})))
        seed = int(profile.get("mavenn_seed", config["p0_2"].get("mavenn_seed", 20260618)))
        phi, wt_phi, metrics = fit_or_surrogate_mavenn(
            assay=assay,
            wildtype=wildtype,
            cfg=cfg,
            seed=seed,
            allow_quick_surrogate=allow_surrogate,
            profile_index=profile_index,
        )
        recovered_column = f"roundtrip_latent_phi__{profile_name}"
        assay[recovered_column] = phi
        recovered_rows = metric_rows_for_score(
            assay=assay,
            wildtype=wildtype,
            wt_score=wt_phi,
            score_column=recovered_column,
            graphlets=graphlets,
            candidate_model="recovered_mavenn_latent_phi_after_roundtrip",
            parameter_set_id=profile_name,
            phase="P0.2",
            extra={
                "spectrum_profile": profile_name,
                "measurement_stage": "recovered_after_mavenn",
                "quick_surrogate": bool(metrics.get("quick_surrogate", False)),
            },
        )
        metric_rows.extend(recovered_rows)
        mavenn_rows.append(
            {
                "phase": "P0.2",
                "spectrum_profile": profile_name,
                "mavenn_seed": seed,
                **{f"mavenn_{key}": value for key, value in metrics.items() if not isinstance(value, dict)},
                **{
                    f"mavenn_test_{key}": value
                    for key, value in dict(metrics.get("test_metrics", {})).items()
                },
                **{
                    f"mavenn_train_{key}": value
                    for key, value in dict(metrics.get("train_metrics", {})).items()
                },
            }
        )

        injected_df = pd.DataFrame(injected_rows)
        recovered_df = pd.DataFrame(recovered_rows)
        for metric in P0_2_METRICS:
            if metric not in injected_df or metric not in recovered_df:
                continue
            injected_mean = float(pd.to_numeric(injected_df[metric], errors="coerce").mean())
            recovered_mean = float(pd.to_numeric(recovered_df[metric], errors="coerce").mean())
            transfer_rows.append(
                {
                    "phase": "P0.2",
                    "spectrum_profile": profile_name,
                    "metric": metric,
                    "injected_mean": injected_mean,
                    "recovered_mean": recovered_mean,
                    "recovered_minus_injected": float(recovered_mean - injected_mean),
                    "recovered_to_injected_ratio": (
                        float(recovered_mean / injected_mean) if abs(injected_mean) > 1e-12 else math.nan
                    ),
                }
            )
    return all_terms, metric_rows, transfer_rows, mavenn_rows


def parse_feature_list(raw: str) -> list[tuple[int, str]]:
    features: list[tuple[int, str]] = []
    for piece in str(raw).split(";"):
        if not piece:
            continue
        pos, aa = piece.split(":", 1)
        features.append((int(pos), aa))
    return features


def module_clusters_from_catalog(path: Path) -> dict[str, list[list[int]]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path)
    out: dict[str, list[list[int]]] = {}
    for row in df.to_dict(orient="records"):
        source = str(row.get("module_source", "empirical_calibration_plaquette"))
        features = parse_feature_list(str(row.get("module_features", "")))
        positions = sorted({int(pos) for pos, _aa in features})
        if len(positions) >= 2:
            out.setdefault(source, []).append(positions)
    return out


def sparse_cluster_positions(
    *,
    config: dict[str, Any],
    graphlets: list[CompleteGraphlet],
    project_root: Path,
) -> list[list[int]]:
    cfg = sparse_walsh_config(config)
    if cfg.get("cluster_positions"):
        return [
            sorted({int(pos) for pos in cluster})
            for cluster in cfg.get("cluster_positions", [])
            if len({int(pos) for pos in cluster}) >= 2
        ]
    module_path = resolve_project_path(project_root, str(config["p0_3"].get("module_catalog_path", "")))
    modules_by_source = module_clusters_from_catalog(module_path)
    source = str(cfg.get("support_source", "empirical_calibration_plaquette"))
    clusters = list(modules_by_source.get(source, []))
    if bool(cfg.get("include_graphlet_feature_clusters", True)):
        splits = {str(value) for value in cfg.get("support_graphlet_splits", ["calibration", "test"])}
        split_map = split_by_graphlet_id(graphlets)
        seen = {tuple(cluster) for cluster in clusters}
        for graphlet in graphlets:
            if split_map[graphlet.graphlet_id] not in splits:
                continue
            positions = tuple(sorted({int(pos) for pos, _aa in graphlet.features}))
            if len(positions) >= 2 and positions not in seen:
                clusters.append(list(positions))
                seen.add(positions)
    return clusters


def allocate_sparse_counts(weights: dict[int, float], *, n_terms: int) -> dict[int, int]:
    expected = {order: float(weight) * int(n_terms) for order, weight in weights.items()}
    counts = {order: int(math.floor(value)) for order, value in expected.items()}
    remaining = int(n_terms) - int(sum(counts.values()))
    ranked = sorted(expected, key=lambda order: (-(expected[order] - counts[order]), -weights[order], order))
    for order in ranked[: max(remaining, 0)]:
        counts[order] += 1
    return {order: count for order, count in counts.items() if count > 0}


def support_for_term(
    *,
    term_id: int,
    order: int,
    seed: int,
    clusters: list[list[int]],
    sequence_length: int,
) -> tuple[int, ...]:
    order = int(np.clip(int(order), 1, int(sequence_length)))
    if clusters:
        cluster = tuple(clusters[stable_uint64(seed, "support-cluster", int(term_id)) % len(clusters)])
        candidates = list(cluster)
    else:
        candidates = []
    neighbor_map: dict[int, set[int]] = {}
    for cluster in clusters:
        for pos in cluster:
            neighbor_map.setdefault(int(pos), set()).update(int(other) for other in cluster if int(other) != int(pos))
    if len(candidates) < order and candidates:
        expanded = set(candidates)
        for pos in list(candidates):
            expanded.update(neighbor_map.get(int(pos), set()))
        candidates = sorted(expanded)
    if len(candidates) < order:
        expanded = set(candidates)
        ranked = sorted(
            range(int(sequence_length)),
            key=lambda pos: stable_uint64(seed, "support-fill", int(term_id), int(pos)),
        )
        for pos in ranked:
            expanded.add(int(pos))
            if len(expanded) >= order:
                break
        candidates = sorted(expanded)
    chosen = sorted(
        candidates,
        key=lambda pos: stable_uint64(seed, "support-pick", int(term_id), int(pos)),
    )[:order]
    return tuple(sorted(int(pos) for pos in chosen))


def run_p0_3(
    *,
    config: dict[str, Any],
    payload: dict[str, Any],
    graphlets: list[CompleteGraphlet],
    project_root: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    cfg = sparse_walsh_config(config)
    clusters = sparse_cluster_positions(config=config, graphlets=graphlets, project_root=project_root)
    graphlet_position_sets = {
        graphlet.graphlet_id: set(int(pos) for pos, _aa in graphlet.features)
        for graphlet in graphlets
    }
    split_map = split_by_graphlet_id(graphlets)
    all_positions = set().union(*graphlet_position_sets.values()) if graphlet_position_sets else set()
    cluster_positions = set(pos for cluster in clusters for pos in cluster)
    position_rows = [
        {
            "phase": "P0.3",
            "assay_id": ASSAY_ID,
            "n_dense_graphlets": int(len(graphlets)),
            "n_graphlet_positions": int(len(all_positions)),
            "n_generator_clusters": int(len(clusters)),
            "n_generator_cluster_positions": int(len(cluster_positions)),
            "n_shared_positions": int(len(all_positions & cluster_positions)),
            "generator_cluster_position_fraction_in_graphlet_pool": (
                float(len(all_positions & cluster_positions) / len(cluster_positions)) if cluster_positions else math.nan
            ),
            "graphlet_position_fraction_in_generator_cluster_pool": (
                float(len(all_positions & cluster_positions) / len(all_positions)) if all_positions else math.nan
            ),
        }
    ]
    orders = sparse_walsh_orders(config)
    weights = {
        order: weight
        for order, weight in zip(
            orders,
            sparse_walsh_order_weights(raw={}, config=config),
        )
    }
    counts = allocate_sparse_counts(weights, n_terms=int(cfg.get("terms", 1024)))
    term_rows: list[dict[str, Any]] = []
    term_id = 0
    for order, count in sorted(counts.items()):
        for _ in range(int(count)):
            support = set(
                support_for_term(
                    term_id=term_id,
                    order=int(order),
                    seed=int(cfg.get("seed", 20260617)),
                    clusters=clusters,
                    sequence_length=len(str(payload["wildtype"])),
                )
            )
            contained_graphlets = [
                graphlet_id
                for graphlet_id, positions in graphlet_position_sets.items()
                if support.issubset(positions)
            ]
            calibration_hits = [
                gid for gid in contained_graphlets if split_map.get(gid) == "calibration"
            ]
            test_hits = [gid for gid in contained_graphlets if split_map.get(gid) == "test"]
            term_rows.append(
                {
                    "phase": "P0.3",
                    "term_id": int(term_id),
                    "order": int(order),
                    "support_positions": ";".join(str(pos) for pos in sorted(support)),
                    "support_size": int(len(support)),
                    "support_fully_in_graphlet_position_pool": bool(support.issubset(all_positions)),
                    "contained_in_any_sampled_graphlet": bool(contained_graphlets),
                    "contained_in_calibration_graphlet": bool(calibration_hits),
                    "contained_in_test_graphlet": bool(test_hits),
                    "n_containing_graphlets": int(len(contained_graphlets)),
                    "n_containing_calibration_graphlets": int(len(calibration_hits)),
                    "n_containing_test_graphlets": int(len(test_hits)),
                    "containing_graphlets": ";".join(contained_graphlets[:12]),
                }
            )
            term_id += 1
    term_df = pd.DataFrame(term_rows)
    summary_rows: list[dict[str, Any]] = []
    if not term_df.empty:
        for order, sub in term_df.groupby("order", dropna=False):
            summary_rows.append(
                {
                    "phase": "P0.3",
                    "order": int(order),
                    "n_terms": int(len(sub)),
                    "fraction_support_fully_in_graphlet_position_pool": float(
                        pd.to_numeric(sub["support_fully_in_graphlet_position_pool"]).mean()
                    ),
                    "fraction_contained_in_any_sampled_graphlet": float(
                        pd.to_numeric(sub["contained_in_any_sampled_graphlet"]).mean()
                    ),
                    "fraction_contained_in_calibration_graphlet": float(
                        pd.to_numeric(sub["contained_in_calibration_graphlet"]).mean()
                    ),
                    "fraction_contained_in_test_graphlet": float(
                        pd.to_numeric(sub["contained_in_test_graphlet"]).mean()
                    ),
                    "mean_n_containing_graphlets": float(
                        pd.to_numeric(sub["n_containing_graphlets"]).mean()
                    ),
                }
            )
    return position_rows, term_rows, summary_rows


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=args.quick)
    payload = load_source_payload(config, project_root)

    graphlets = graphlets_from_assay(
        assay=payload["assay"],
        wildtype=str(payload["wildtype"]),
        wt_scores={"latent_phi": float(payload["wildtype_phi"])},
        score_columns=["latent_phi"],
        config=config,
    )
    write_csv(output_dir / "graphlet_catalog.csv", graphlet_catalog_rows(graphlets))
    preflight = {
        "experiment_id": config["experiment"]["id"],
        "run_mode": config["experiment"]["mode"],
        "scratch_source": config["experiment"]["scientific_boundary"]["advisory_source"],
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_output_dir": str(source_output_dir(config, project_root)),
        "n_assay_rows": int(len(payload["assay"])),
        "n_graphlets": int(len(graphlets)),
        "graphlet_dimensions": sorted({int(graphlet.dimension) for graphlet in graphlets}),
        "p0_1_refit_count": int(
            len(config["p0_1"]["mavenn_refits"]["seeds"])
            * len(config["p0_1"]["mavenn_refits"]["profiles"])
        ),
        "p0_2_known_spectrum_profiles": [
            str(profile["name"]) for profile in config["p0_2"]["known_spectrum_profiles"]
        ],
        "p0_3_sparse_walsh_terms": int(sparse_walsh_config(config).get("terms", 0)),
        "quick": bool(args.quick),
    }
    write_json(output_dir / "preflight.json", preflight)

    p01_metric_rows, p01_bands, p01_fit_rows, p01_replicates = run_p0_1(
        config=config,
        payload=payload,
        graphlets=graphlets,
        quick=bool(args.quick),
    )
    p02_terms, p02_metric_rows, p02_transfer, p02_mavenn = run_p0_2(
        config=config,
        payload=payload,
        graphlets=graphlets,
        quick=bool(args.quick),
    )
    p03_positions, p03_terms, p03_summary = run_p0_3(
        config=config,
        payload=payload,
        graphlets=graphlets,
        project_root=project_root,
    )

    write_csv(output_dir / "p0_1_target_graphlet_metrics.csv", p01_metric_rows)
    write_csv(output_dir / "p0_1_target_bootstrap.csv", p01_bands)
    write_csv(output_dir / "p0_1_mavenn_refit_metrics.csv", p01_fit_rows)
    write_csv(output_dir / "p0_1_replicate_audit.csv", p01_replicates)
    write_csv(output_dir / "p0_2_injected_terms.csv", p02_terms)
    write_csv(output_dir / "p0_2_roundtrip_graphlet_metrics.csv", p02_metric_rows)
    write_csv(output_dir / "p0_2_roundtrip_transfer.csv", p02_transfer)
    write_csv(output_dir / "p0_2_mavenn_metrics.csv", p02_mavenn)
    write_csv(output_dir / "p0_3_position_pool_audit.csv", p03_positions)
    write_csv(output_dir / "p0_3_support_terms.csv", p03_terms)
    write_csv(output_dir / "p0_3_support_containment_summary.csv", p03_summary)
    write_json(
        output_dir / "summary.json",
        {
            **preflight,
            "n_p0_1_graphlet_metric_rows": int(len(p01_metric_rows)),
            "n_p0_1_bootstrap_rows": int(len(p01_bands)),
            "n_p0_1_mavenn_refits": int(len(p01_fit_rows)),
            "n_p0_2_injected_terms": int(len(p02_terms)),
            "n_p0_2_graphlet_metric_rows": int(len(p02_metric_rows)),
            "n_p0_2_transfer_rows": int(len(p02_transfer)),
            "n_p0_3_support_terms": int(len(p03_terms)),
            "output_files": {
                "preflight": str(output_dir / "preflight.json"),
                "graphlet_catalog": str(output_dir / "graphlet_catalog.csv"),
                "p0_1_target_bootstrap": str(output_dir / "p0_1_target_bootstrap.csv"),
                "p0_1_mavenn_refit_metrics": str(output_dir / "p0_1_mavenn_refit_metrics.csv"),
                "p0_2_roundtrip_transfer": str(output_dir / "p0_2_roundtrip_transfer.csv"),
                "p0_3_support_containment_summary": str(output_dir / "p0_3_support_containment_summary.csv"),
            },
        },
    )


if __name__ == "__main__":
    main()
