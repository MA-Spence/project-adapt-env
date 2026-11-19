#!/usr/bin/env python3
"""Certify non-additive PHOT MAVE-NN measurement operators by round-trip tests."""

from __future__ import annotations

import argparse
import copy
import math
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_ROOT, ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from phot_coupled_landscape_grid_search import load_config, write_csv, write_json  # noqa: E402
from phot_miscalibration_p0_isolation import (  # noqa: E402
    ASSAY_ID,
    P0_2_METRICS,
    fit_or_surrogate_mavenn,
    graphlet_catalog_rows,
    graphlets_from_assay,
    load_source_payload,
    metric_rows_for_score,
    mutation_mask_for_sequence,
    normalize_assay_scale,
    raw_term_values_for_masks,
    source_output_dir,
    synthetic_terms,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def profile_catalog(config: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for profile in config["operator_profiles"]:
        cfg = copy.deepcopy(config["mavenn"])
        cfg.update(dict(profile.get("mavenn_overrides", {})))
        rows.append(
            {
                "operator_profile": str(profile["name"]),
                "description": str(profile.get("description", "")),
                "gpmap_type": cfg.get("gpmap_type"),
                "ge_nonlinearity_type": cfg.get("ge_nonlinearity_type"),
                "theta_regularization": cfg.get("theta_regularization"),
                "eta_regularization": cfg.get("eta_regularization"),
                "epochs": cfg.get("epochs"),
                "batch_size": cfg.get("batch_size"),
                "max_variants_per_assay": cfg.get("max_variants_per_assay"),
                "sample_seed": cfg.get("sample_seed"),
            }
        )
    return rows


def alphabet_size(name: str) -> int:
    if str(name).lower() == "protein":
        return 20
    if str(name).lower() in {"dna", "rna"}:
        return 4
    return 20


def pairwise_parameter_bytes(wildtype: str, cfg: dict[str, Any]) -> int:
    length = len(str(wildtype))
    states = max(alphabet_size(str(cfg.get("alphabet", "protein"))) - 1, 1)
    n_additive = length * states
    n_pairwise = (length * max(length - 1, 0) // 2) * (states**2)
    return int((n_additive + n_pairwise) * 8)


def preflight_skip_reason(config: dict[str, Any], wildtype: str, cfg: dict[str, Any]) -> str:
    preflight = dict(config.get("operator_preflight", {}))
    if not bool(preflight.get("enabled", True)):
        return ""
    gpmap = str(cfg.get("gpmap_type", "")).lower()
    if gpmap != "pairwise":
        return ""
    max_length = int(preflight.get("max_pairwise_sequence_length", 0))
    if max_length > 0 and len(str(wildtype)) > max_length:
        return (
            "pairwise_gpmap_preflight_skip: sequence_length="
            f"{len(str(wildtype))} exceeds max_pairwise_sequence_length={max_length}"
        )
    limit = int(preflight.get("max_pairwise_parameter_bytes", 1_000_000_000))
    estimate = pairwise_parameter_bytes(wildtype, cfg)
    if estimate > limit:
        return (
            "pairwise_gpmap_preflight_skip: estimated_parameter_bytes="
            f"{estimate} exceeds max_pairwise_parameter_bytes={limit}"
        )
    return ""


def ratio_pass(
    *,
    injected: float,
    recovered: float,
    ratio_min: float,
    ratio_max: float,
    near_zero_tolerance: float,
) -> tuple[float, bool]:
    if not np.isfinite(injected) or not np.isfinite(recovered):
        return math.nan, False
    if abs(float(injected)) <= float(near_zero_tolerance):
        return math.nan, abs(float(recovered) - float(injected)) <= float(near_zero_tolerance)
    ratio = float(recovered) / float(injected)
    return ratio, bool(float(ratio_min) <= ratio <= float(ratio_max))


def run_roundtrip(config: dict[str, Any], project_root: Path, output_dir: Path, *, quick: bool) -> dict[str, Any]:
    payload = load_source_payload(config, project_root)
    assay_base = payload["assay"].copy()
    wildtype = str(payload["wildtype"])
    graphlets = graphlets_from_assay(
        assay=assay_base,
        wildtype=wildtype,
        wt_scores={"latent_phi": float(payload["wildtype_phi"])},
        score_columns=["latent_phi"],
        config=config,
    )
    write_csv(output_dir / "graphlet_catalog.csv", graphlet_catalog_rows(graphlets))
    write_csv(output_dir / "operator_profile_catalog.csv", profile_catalog(config))

    feature_pool = sorted(
        {feature for graphlet in graphlets for feature in graphlet.features},
        key=lambda feature: (int(feature[0]), str(feature[1])),
    )
    masks = [mutation_mask_for_sequence(wildtype, seq) for seq in assay_base["mutated_sequence"].astype(str)]
    wt_mask = frozenset()
    all_terms: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    transfer_rows: list[dict[str, Any]] = []
    mavenn_rows: list[dict[str, Any]] = []
    selected_metrics = [str(metric) for metric in config["certification"]["metrics"]]
    ratio_min = float(config["certification"].get("ratio_min", 0.7))
    ratio_max = float(config["certification"].get("ratio_max", 1.4))
    near_zero_tolerance = float(config["certification"].get("near_zero_absolute_tolerance", 1e-6))
    allow_surrogate = bool(quick and config.get("quick_surrogate_without_mavenn", True))

    profile_rows_cache: dict[str, list[dict[str, Any]]] = {}
    profile_wt_scores: dict[str, float] = {}
    for spectrum_profile in config["known_spectrum_profiles"]:
        spectrum_profile = dict(spectrum_profile)
        spectrum_name = str(spectrum_profile["name"])
        terms = synthetic_terms(feature_pool=feature_pool, profile=spectrum_profile)
        # Build term rows explicitly to avoid relying on dictionary-comprehension side effects.
        for term in terms:
            row = {
                key: (";".join(f"{pos}:{aa}" for pos, aa in value) if key == "support" else value)
                for key, value in term.items()
            }
            row["spectrum_profile"] = spectrum_name
            all_terms.append(row)

        raw_score = raw_term_values_for_masks(
            masks,
            terms=terms,
            additive_weight=float(spectrum_profile.get("additive_weight", 0.0)),
        )
        synthetic_score = normalize_assay_scale(raw_score)
        raw_mean = float(np.mean(raw_score))
        raw_sd = max(float(np.std(raw_score)), 1e-12)
        wt_raw = float(
            raw_term_values_for_masks(
                [wt_mask],
                terms=terms,
                additive_weight=float(spectrum_profile.get("additive_weight", 0.0)),
            )[0]
        )
        profile_wt_scores[spectrum_name] = float((wt_raw - raw_mean) / raw_sd)
        assay = assay_base.copy()
        assay["synthetic_dms_score"] = synthetic_score
        assay["DMS_score"] = synthetic_score
        injected_rows = metric_rows_for_score(
            assay=assay,
            wildtype=wildtype,
            wt_score=profile_wt_scores[spectrum_name],
            score_column="synthetic_dms_score",
            graphlets=graphlets,
            candidate_model="known_injected_synthetic_landscape",
            parameter_set_id=spectrum_name,
            phase="WS-A",
            extra={
                "operator_profile": "direct",
                "spectrum_profile": spectrum_name,
                "measurement_stage": "injected_direct",
                "quick_surrogate": False,
            },
        )
        profile_rows_cache[spectrum_name] = injected_rows

    for operator_index, operator_profile in enumerate(config["operator_profiles"], start=1):
        operator_name = str(operator_profile["name"])
        operator_cfg = copy.deepcopy(config["mavenn"])
        operator_cfg.update(dict(operator_profile.get("mavenn_overrides", {})))
        for spectrum_index, spectrum_profile in enumerate(config["known_spectrum_profiles"], start=1):
            spectrum_profile = dict(spectrum_profile)
            spectrum_name = str(spectrum_profile["name"])
            injected_rows = [
                {**row, "operator_profile": operator_name}
                for row in profile_rows_cache[spectrum_name]
            ]
            metric_rows.extend(injected_rows)

            raw_score = raw_term_values_for_masks(
                masks,
                terms=synthetic_terms(feature_pool=feature_pool, profile=spectrum_profile),
                additive_weight=float(spectrum_profile.get("additive_weight", 0.0)),
            )
            assay = assay_base.copy()
            assay["DMS_score"] = normalize_assay_scale(raw_score)
            seed = int(operator_profile.get("mavenn_seed", spectrum_profile.get("mavenn_seed", operator_cfg.get("sample_seed", 0))))
            fit_error = ""
            metrics: dict[str, Any] = {}
            recovered_rows: list[dict[str, Any]] = []
            try:
                skip_reason = preflight_skip_reason(config, wildtype, operator_cfg)
                if skip_reason:
                    fit_error = skip_reason
                    metrics = {
                        "quick_surrogate": False,
                        "fit_failed": True,
                        "fit_skipped_preflight": True,
                        "fit_error": fit_error,
                        "estimated_pairwise_parameter_bytes": pairwise_parameter_bytes(wildtype, operator_cfg),
                    }
                else:
                    phi, wt_phi, metrics = fit_or_surrogate_mavenn(
                        assay=assay,
                        wildtype=wildtype,
                        cfg=operator_cfg,
                        seed=seed,
                        allow_quick_surrogate=allow_surrogate,
                        profile_index=operator_index * 100 + spectrum_index,
                    )
                    recovered_column = f"roundtrip_latent_phi__{operator_name}__{spectrum_name}"
                    assay[recovered_column] = phi
                    recovered_rows = metric_rows_for_score(
                        assay=assay,
                        wildtype=wildtype,
                        wt_score=wt_phi,
                        score_column=recovered_column,
                        graphlets=graphlets,
                        candidate_model="recovered_mavenn_latent_phi_after_roundtrip",
                        parameter_set_id=spectrum_name,
                        phase="WS-A",
                        extra={
                            "operator_profile": operator_name,
                            "spectrum_profile": spectrum_name,
                            "measurement_stage": "recovered_after_mavenn",
                            "quick_surrogate": bool(metrics.get("quick_surrogate", False)),
                        },
                    )
                    metric_rows.extend(recovered_rows)
            except Exception as exc:  # pragma: no cover - runtime certification failure path
                fit_error = "".join(traceback.format_exception_only(type(exc), exc)).strip()
                metrics = {"quick_surrogate": False, "fit_failed": True, "fit_error": fit_error}

            mavenn_rows.append(
                {
                    "phase": "WS-A",
                    "operator_profile": operator_name,
                    "spectrum_profile": spectrum_name,
                    "mavenn_seed": seed,
                    "fit_failed": bool(fit_error),
                    "fit_error": fit_error,
                    "gpmap_type": operator_cfg.get("gpmap_type"),
                    "theta_regularization": operator_cfg.get("theta_regularization"),
                    "eta_regularization": operator_cfg.get("eta_regularization"),
                    "epochs": operator_cfg.get("epochs"),
                    "batch_size": operator_cfg.get("batch_size"),
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
                injected_mean = (
                    float(pd.to_numeric(injected_df.get(metric), errors="coerce").mean())
                    if metric in injected_df
                    else math.nan
                )
                recovered_mean = (
                    float(pd.to_numeric(recovered_df.get(metric), errors="coerce").mean())
                    if metric in recovered_df
                    else math.nan
                )
                ratio, passed = ratio_pass(
                    injected=injected_mean,
                    recovered=recovered_mean,
                    ratio_min=ratio_min,
                    ratio_max=ratio_max,
                    near_zero_tolerance=near_zero_tolerance,
                )
                transfer_rows.append(
                    {
                        "phase": "WS-A",
                        "operator_profile": operator_name,
                        "spectrum_profile": spectrum_name,
                        "metric": metric,
                        "injected_mean": injected_mean,
                        "recovered_mean": recovered_mean,
                        "recovered_minus_injected": (
                            float(recovered_mean - injected_mean)
                            if np.isfinite(injected_mean) and np.isfinite(recovered_mean)
                            else math.nan
                        ),
                        "recovered_to_injected_ratio": ratio,
                        "certification_metric": metric in selected_metrics,
                        "ratio_within_tolerance": bool(passed) if metric in selected_metrics else math.nan,
                        "fit_failed": bool(fit_error),
                    }
                )

    return {
        "payload": payload,
        "graphlets": graphlets,
        "terms": all_terms,
        "metrics": metric_rows,
        "transfer": transfer_rows,
        "mavenn": mavenn_rows,
    }


def summarize_certification(config: dict[str, Any], transfer_rows: list[dict[str, Any]], mavenn_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    transfer = pd.DataFrame(transfer_rows)
    mavenn = pd.DataFrame(mavenn_rows)
    selected_metrics = [str(metric) for metric in config["certification"]["metrics"]]
    rows: list[dict[str, Any]] = []
    for profile in config["operator_profiles"]:
        operator_name = str(profile["name"])
        sub = transfer[(transfer["operator_profile"] == operator_name) & (transfer["metric"].isin(selected_metrics))]
        fit_sub = mavenn[mavenn["operator_profile"] == operator_name] if not mavenn.empty else pd.DataFrame()
        n_fit_failures = int(pd.to_numeric(fit_sub.get("fit_failed", pd.Series(dtype=bool)), errors="coerce").fillna(False).sum())
        pass_values = pd.to_numeric(sub.get("ratio_within_tolerance", pd.Series(dtype=float)), errors="coerce").dropna()
        n_required = int(len(selected_metrics) * len(config["known_spectrum_profiles"]))
        rows.append(
            {
                "phase": "WS-A",
                "operator_profile": operator_name,
                "gpmap_type": dict(profile.get("mavenn_overrides", {})).get("gpmap_type", config["mavenn"].get("gpmap_type")),
                "n_spectrum_profiles": int(len(config["known_spectrum_profiles"])),
                "n_certification_metrics": int(len(selected_metrics)),
                "n_required_metric_profile_checks": n_required,
                "n_ratio_checks_reported": int(len(pass_values)),
                "n_ratio_checks_passed": int(pass_values.sum()) if len(pass_values) else 0,
                "n_fit_failures": n_fit_failures,
                "operator_certified": bool(n_fit_failures == 0 and len(pass_values) == n_required and int(pass_values.sum()) == n_required),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=bool(args.quick))
    result = run_roundtrip(config, project_root, output_dir, quick=bool(args.quick))

    write_csv(output_dir / "injected_terms.csv", result["terms"])
    write_csv(output_dir / "roundtrip_graphlet_metrics.csv", result["metrics"])
    write_csv(output_dir / "roundtrip_transfer.csv", result["transfer"])
    write_csv(output_dir / "mavenn_metrics.csv", result["mavenn"])
    certification_summary = summarize_certification(config, result["transfer"], result["mavenn"])
    write_csv(output_dir / "operator_certification_summary.csv", certification_summary)

    payload = result["payload"]
    graphlets = result["graphlets"]
    preflight = {
        "experiment_id": config["experiment"]["id"],
        "mode": config["experiment"]["mode"],
        "scratch_source": config["experiment"]["scientific_boundary"]["advisory_source"],
        "assay_id": ASSAY_ID,
        "source_run_id": config["source_run"]["run_id"],
        "source_output_dir": str(source_output_dir(config, project_root)),
        "n_assay_rows": int(len(payload["assay"])),
        "n_graphlets": int(len(graphlets)),
        "operator_profiles": [str(profile["name"]) for profile in config["operator_profiles"]],
        "known_spectrum_profiles": [str(profile["name"]) for profile in config["known_spectrum_profiles"]],
        "certification_metrics": [str(metric) for metric in config["certification"]["metrics"]],
        "quick": bool(args.quick),
    }
    write_json(output_dir / "preflight.json", preflight)
    write_json(
        output_dir / "summary.json",
        {
            **preflight,
            "n_injected_terms": int(len(result["terms"])),
            "n_roundtrip_graphlet_metric_rows": int(len(result["metrics"])),
            "n_roundtrip_transfer_rows": int(len(result["transfer"])),
            "n_mavenn_metric_rows": int(len(result["mavenn"])),
            "operator_certification_summary": certification_summary,
        },
    )


if __name__ == "__main__":
    main()
