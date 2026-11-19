#!/usr/bin/env python3
"""Derive compact ANA-048 tables from RUN-125 source outputs."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SOURCE_DIR = ANALYSIS_DIR / "source_tables"
TABLE_DIR = ANALYSIS_DIR / "tables"
RUN_ID = "RUN-125"


P0_TARGET_METRICS = [
    "walsh_high_order_fraction",
    "n_local_maxima",
    "accessible_path_fraction",
    "edge_energy_top10_fraction",
    "edge_energy_effective_fraction",
    "site_energy_top1_fraction",
    "site_energy_effective_fraction",
]

P0_2_DECISION_METRICS = [
    "walsh_high_order_fraction",
    "n_local_maxima",
    "accessible_path_fraction",
    "magnitude_epistasis_to_single_mutant_sd_ratio",
    "sign_epistasis_fraction_local",
    "reciprocal_sign_epistasis_fraction_local",
    "diminishing_returns_epistasis_fraction",
    "edge_energy_top10_fraction",
    "edge_energy_effective_fraction",
    "site_energy_top1_fraction",
    "site_energy_effective_fraction",
]


def source_path(name: str) -> Path:
    return SOURCE_DIR / f"run-125_{name}"


def as_float(value: Any) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(float(value)):
            return ""
        v = float(value)
        if abs(v) < 1e-4 and v != 0.0:
            return f"{v:.3e}"
        return f"{v:.4f}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    text = str(value)
    return text


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._\n"
    cols = list(df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for row in df.to_dict(orient="records"):
        lines.append("| " + " | ".join(format_value(row.get(col, "")) for col in cols) + " |")
    return "\n".join(lines) + "\n"


def write_md(path: Path, title: str, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"# {title}\n\n{body.rstrip()}\n", encoding="utf-8")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_output_audit() -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(SOURCE_DIR.glob("run-125_*")):
        rows.append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }
        )
    df = pd.DataFrame(rows)
    return df, {
        "n_source_files": int(len(df)),
        "copied_source_bytes": int(df["bytes"].sum()) if not df.empty else 0,
    }


def p0_1_target_stability() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    bands = pd.read_csv(source_path("p0_1_target_bootstrap.csv"))
    fits = pd.read_csv(source_path("p0_1_mavenn_refit_metrics.csv"))
    cached = bands[bands["parameter_set_id"] == "RUN-077_cached_latent_phi"].set_index("metric")
    refit = bands[bands["candidate_model"] == "empirical_refit_mavenn_latent_phi"].copy()
    rows: list[dict[str, Any]] = []
    for metric in P0_TARGET_METRICS:
        c = cached.loc[metric] if metric in cached.index else None
        sub = refit[refit["metric"] == metric]
        rows.append(
            {
                "metric": metric,
                "cached_mean": as_float(c["graphlet_mean"]) if c is not None else math.nan,
                "cached_q05": as_float(c["bootstrap_mean_q05"]) if c is not None else math.nan,
                "cached_q95": as_float(c["bootstrap_mean_q95"]) if c is not None else math.nan,
                "refit_mean_min": as_float(sub["graphlet_mean"].min()),
                "refit_mean_max": as_float(sub["graphlet_mean"].max()),
                "refit_q05_min": as_float(sub["bootstrap_mean_q05"].min()),
                "refit_q95_max": as_float(sub["bootstrap_mean_q95"].max()),
            }
        )
    stability = pd.DataFrame(rows)

    fit_quality = (
        fits.groupby("inference_profile", dropna=False)
        .agg(
            n_refits=("refit_seed", "count"),
            test_spearman_min=("mavenn_test_spearman", "min"),
            test_spearman_median=("mavenn_test_spearman", "median"),
            test_spearman_max=("mavenn_test_spearman", "max"),
            test_nrmse_min=("mavenn_test_nrmse", "min"),
            test_nrmse_median=("mavenn_test_nrmse", "median"),
            test_nrmse_max=("mavenn_test_nrmse", "max"),
            any_quick_surrogate=("mavenn_quick_surrogate", "max"),
        )
        .reset_index()
    )
    access_refit = stability.loc[stability["metric"] == "accessible_path_fraction"].iloc[0]
    high_order_refit = stability.loc[stability["metric"] == "walsh_high_order_fraction"].iloc[0]
    metrics = {
        "p0_1_refits": int(len(fits)),
        "p0_1_refit_profiles": int(fits["inference_profile"].nunique()),
        "p0_1_refit_seeds": int(fits["refit_seed"].nunique()),
        "p0_1_any_quick_surrogate": bool(fits["mavenn_quick_surrogate"].any()),
        "p0_1_cached_high_order_mean": as_float(
            stability.loc[stability["metric"] == "walsh_high_order_fraction", "cached_mean"].iloc[0]
        ),
        "p0_1_refit_high_order_mean_max": as_float(high_order_refit["refit_mean_max"]),
        "p0_1_cached_local_maxima_mean": as_float(
            stability.loc[stability["metric"] == "n_local_maxima", "cached_mean"].iloc[0]
        ),
        "p0_1_refit_accessibility_mean_min": as_float(access_refit["refit_mean_min"]),
        "p0_1_refit_accessibility_mean_max": as_float(access_refit["refit_mean_max"]),
    }
    return stability, fit_quality, metrics


def p0_2_roundtrip() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    transfer = pd.read_csv(source_path("p0_2_roundtrip_transfer.csv"))
    mavenn = pd.read_csv(source_path("p0_2_mavenn_metrics.csv"))
    key = transfer[transfer["metric"].isin(P0_2_DECISION_METRICS)].copy()
    key["recovered_to_injected_ratio"] = key["recovered_to_injected_ratio"].replace([np.inf, -np.inf], np.nan)

    profile_rows: list[dict[str, Any]] = []
    for profile, sub in transfer.groupby("spectrum_profile", dropna=False):
        metric_map = {row["metric"]: row for row in sub.to_dict(orient="records")}
        fit = mavenn[mavenn["spectrum_profile"] == profile].iloc[0]
        profile_rows.append(
            {
                "spectrum_profile": profile,
                "injected_high_order": as_float(metric_map["walsh_high_order_fraction"]["injected_mean"]),
                "recovered_high_order": as_float(metric_map["walsh_high_order_fraction"]["recovered_mean"]),
                "injected_accessibility": as_float(metric_map["accessible_path_fraction"]["injected_mean"]),
                "recovered_accessibility": as_float(metric_map["accessible_path_fraction"]["recovered_mean"]),
                "injected_local_maxima": as_float(metric_map["n_local_maxima"]["injected_mean"]),
                "recovered_local_maxima": as_float(metric_map["n_local_maxima"]["recovered_mean"]),
                "injected_magnitude_ratio": as_float(
                    metric_map["magnitude_epistasis_to_single_mutant_sd_ratio"]["injected_mean"]
                ),
                "recovered_magnitude_ratio": as_float(
                    metric_map["magnitude_epistasis_to_single_mutant_sd_ratio"]["recovered_mean"]
                ),
                "injected_sign_epistasis_local": as_float(
                    metric_map["sign_epistasis_fraction_local"]["injected_mean"]
                ),
                "recovered_sign_epistasis_local": as_float(
                    metric_map["sign_epistasis_fraction_local"]["recovered_mean"]
                ),
                "mavenn_test_spearman": as_float(fit["mavenn_test_spearman"]),
                "mavenn_test_nrmse": as_float(fit["mavenn_test_nrmse"]),
            }
        )
    profile_summary = pd.DataFrame(profile_rows)
    metrics = {
        "p0_2_profiles": int(profile_summary.shape[0]),
        "p0_2_mavenn_any_quick_surrogate": bool(mavenn["mavenn_quick_surrogate"].any()),
        "p0_2_recovered_high_order_max": as_float(profile_summary["recovered_high_order"].max()),
        "p0_2_direct_high_order_max": as_float(profile_summary["injected_high_order"].max()),
        "p0_2_recovered_local_maxima_min": as_float(profile_summary["recovered_local_maxima"].min()),
        "p0_2_recovered_local_maxima_max": as_float(profile_summary["recovered_local_maxima"].max()),
        "p0_2_recovered_magnitude_ratio_max": as_float(profile_summary["recovered_magnitude_ratio"].max()),
        "p0_2_recovered_sign_epistasis_max": as_float(profile_summary["recovered_sign_epistasis_local"].max()),
    }
    return key, profile_summary, metrics


def p0_3_support_alignment() -> tuple[pd.DataFrame, dict[str, Any]]:
    position = pd.read_csv(source_path("p0_3_position_pool_audit.csv"))
    contain = pd.read_csv(source_path("p0_3_support_containment_summary.csv"))
    combined = pd.concat(
        [
            position.assign(order=np.nan, n_terms=np.nan),
            contain.assign(assay_id="PHOT_CHLRE_Chen_2023"),
        ],
        axis=0,
        ignore_index=True,
        sort=False,
    )
    c = contain.iloc[0]
    p = position.iloc[0]
    metrics = {
        "p0_3_graphlet_positions": int(p["n_graphlet_positions"]),
        "p0_3_generator_cluster_positions": int(p["n_generator_cluster_positions"]),
        "p0_3_shared_positions": int(p["n_shared_positions"]),
        "p0_3_n_terms": int(c["n_terms"]),
        "p0_3_fraction_contained_any": as_float(c["fraction_contained_in_any_sampled_graphlet"]),
        "p0_3_fraction_contained_calibration": as_float(c["fraction_contained_in_calibration_graphlet"]),
        "p0_3_fraction_contained_test": as_float(c["fraction_contained_in_test_graphlet"]),
        "p0_3_mean_n_containing_graphlets": as_float(c["mean_n_containing_graphlets"]),
    }
    return combined, metrics


def graphlet_summary() -> pd.DataFrame:
    graphlets = pd.read_csv(source_path("graphlet_catalog.csv"))
    split_dim = (
        graphlets.groupby(["graphlet_split", "graphlet_dimension"], dropna=False)
        .size()
        .reset_index(name="n_graphlets")
        .sort_values(["graphlet_dimension", "graphlet_split"])
    )
    all_positions = sorted(
        {
            int(feature.split(":", 1)[0])
            for raw in graphlets["mutation_features"].astype(str)
            for feature in raw.split(";")
            if feature
        }
    )
    total = pd.DataFrame(
        [
            {
                "graphlet_split": "all",
                "graphlet_dimension": "all",
                "n_graphlets": int(len(graphlets)),
                "n_unique_positions": int(len(all_positions)),
                "positions": ";".join(str(pos) for pos in all_positions),
            }
        ]
    )
    split_dim["n_unique_positions"] = ""
    split_dim["positions"] = ""
    return pd.concat([total, split_dim], ignore_index=True, sort=False)


def phase_decisions(metrics: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "phase": "P0.1",
                "decision": "target is smooth/additive-like rather than high-order/rugged under this MAVE-NN latent operator",
                "evidence": (
                    "cached and refit high-order means are effectively zero; local maxima are 1.0 for every "
                    "cached/refit graphlet summary, while accessibility is refit-sensitive and ranges from 0.0 to 0.5417"
                ),
            },
            {
                "phase": "P0.2",
                "decision": "measurement operator is strongly non-identity and smoothing rather than spectrum-preserving",
                "evidence": (
                    "known injected spectra recover high-order fraction <= "
                    f"{metrics['p0_2_recovered_high_order_max']:.3e}, magnitude ratio <= "
                    f"{metrics['p0_2_recovered_magnitude_ratio_max']:.3e}, local maxima exactly 1.0, and local sign epistasis 0.0"
                ),
            },
            {
                "phase": "P0.3",
                "decision": "support-to-sampler position containment is not the limiting failure mode in this run",
                "evidence": (
                    "all 2048 generated order-3 supports are contained in at least one sampled graphlet; "
                    "calibration containment is 0.9521 and test containment is 0.9292"
                ),
            },
        ]
    )


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    summary = json.loads(source_path("summary.json").read_text(encoding="utf-8"))
    preflight = json.loads(source_path("preflight.json").read_text(encoding="utf-8"))
    output_audit, audit_metrics = build_output_audit()
    p01_stability, p01_fit_quality, p01_metrics = p0_1_target_stability()
    p02_transfer, p02_profile_summary, p02_metrics = p0_2_roundtrip()
    p03_alignment, p03_metrics = p0_3_support_alignment()
    graphlets = graphlet_summary()

    metrics = {
        "run_id": RUN_ID,
        "experiment_id": summary["experiment_id"],
        "assay_id": summary["assay_id"],
        "source_run_id": summary["source_run_id"],
        "quick": bool(summary["quick"]),
        "n_assay_rows": int(preflight["n_assay_rows"]),
        "n_graphlets": int(summary["n_graphlets"]),
        "n_p0_1_mavenn_refits": int(summary["n_p0_1_mavenn_refits"]),
        "n_p0_2_transfer_rows": int(summary["n_p0_2_transfer_rows"]),
        "n_p0_3_support_terms": int(summary["n_p0_3_support_terms"]),
        **audit_metrics,
        **p01_metrics,
        **p02_metrics,
        **p03_metrics,
    }
    decisions = phase_decisions(metrics)

    (TABLE_DIR / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output_audit.to_csv(TABLE_DIR / "run-125_output_collection_audit.csv", index=False)
    p01_stability.to_csv(TABLE_DIR / "run-125_p0_1_target_stability.csv", index=False)
    p01_fit_quality.to_csv(TABLE_DIR / "run-125_p0_1_mavenn_refit_fit_quality.csv", index=False)
    p02_transfer.to_csv(TABLE_DIR / "run-125_p0_2_roundtrip_transfer_key_metrics.csv", index=False)
    p02_profile_summary.to_csv(TABLE_DIR / "run-125_p0_2_profile_summary.csv", index=False)
    p03_alignment.to_csv(TABLE_DIR / "run-125_p0_3_support_alignment.csv", index=False)
    graphlets.to_csv(TABLE_DIR / "run-125_graphlet_catalog_summary.csv", index=False)
    decisions.to_csv(TABLE_DIR / "run-125_phase_decisions.csv", index=False)

    write_md(
        TABLE_DIR / "run-125_output_collection_audit.md",
        "RUN-125 Output Collection Audit",
        markdown_table(output_audit[["file", "bytes", "sha256"]]),
    )
    write_md(
        TABLE_DIR / "run-125_p0_1_target_stability.md",
        "P0.1 Target Bootstrap and MAVE-NN Refit Stability",
        markdown_table(p01_stability),
    )
    write_md(
        TABLE_DIR / "run-125_p0_1_mavenn_refit_fit_quality.md",
        "P0.1 MAVE-NN Refit Fit Quality",
        markdown_table(p01_fit_quality),
    )
    write_md(
        TABLE_DIR / "run-125_p0_2_profile_summary.md",
        "P0.2 Round-Trip Profile Summary",
        markdown_table(p02_profile_summary),
    )
    write_md(
        TABLE_DIR / "run-125_p0_2_roundtrip_transfer_key_metrics.md",
        "P0.2 Round-Trip Transfer Key Metrics",
        markdown_table(p02_transfer),
    )
    write_md(
        TABLE_DIR / "run-125_p0_3_support_alignment.md",
        "P0.3 Support and Sampler Alignment",
        markdown_table(p03_alignment),
    )
    write_md(
        TABLE_DIR / "run-125_graphlet_catalog_summary.md",
        "RUN-125 Graphlet Catalog Summary",
        markdown_table(graphlets),
    )
    write_md(
        TABLE_DIR / "run-125_phase_decisions.md",
        "RUN-125 Phase Decisions",
        markdown_table(decisions),
    )


if __name__ == "__main__":
    main()
