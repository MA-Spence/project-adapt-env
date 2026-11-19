from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SOURCE_DIR = ANALYSIS_DIR / "source_tables"
TABLE_DIR = ANALYSIS_DIR / "tables"


def read_csv(name: str) -> list[dict[str, str]]:
    with (SOURCE_DIR / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in {"", None} else float("nan")


def as_bool(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip().lower() == "true"


def write_csv(name: str, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with (TABLE_DIR / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def format_value(value: object) -> str:
    if isinstance(value, float):
        if value != value:
            return ""
        if abs(value) >= 100:
            return f"{value:.3f}"
        if abs(value) >= 1:
            return f"{value:.4f}"
        return f"{value:.6g}"
    return str(value)


def write_md_table(name: str, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    lines = []
    lines.append("| " + " | ".join(fieldnames) + " |")
    lines.append("| " + " | ".join(["---"] * len(fieldnames)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(format_value(row.get(field, "")) for field in fieldnames) + " |")
    (TABLE_DIR / name).write_text("\n".join(lines) + "\n")


def audit_source_tables() -> tuple[list[dict[str, object]], int]:
    rows = []
    total = 0
    for path in sorted(SOURCE_DIR.iterdir()):
        if not path.is_file():
            continue
        payload = path.read_bytes()
        total += len(payload)
        rows.append(
            {
                "file": path.name,
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        )
    return rows, total


def grouped(rows: list[dict[str, str]], *keys: str) -> dict[tuple[str, ...], list[dict[str, str]]]:
    out: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        out[tuple(row[key] for key in keys)].append(row)
    return out


def main() -> None:
    TABLE_DIR.mkdir(exist_ok=True)

    summary = json.loads((SOURCE_DIR / "run-127_summary.json").read_text())
    monotonicity = read_csv("run-127_monotonicity_checks.csv")
    knob = read_csv("run-127_knob_summary.csv")
    core = read_csv("run-127_core_readout_guardrail_checks.csv")
    threshold = read_csv("run-127_threshold_epistasis_metrics.csv")
    hamming = read_csv("run-127_hamming_layer_metrics.csv")
    support = read_csv("run-127_support_catalog.csv")

    audit_rows, copied_bytes = audit_source_tables()
    write_csv("run-127_output_collection_audit.csv", audit_rows)
    write_md_table("run-127_output_collection_audit.md", audit_rows)

    mono_fields = [
        "support_mode",
        "metric",
        "expected_direction",
        "spearman_rho",
        "direction_pass",
        "n_parameter_sets",
    ]
    mono_rows = [
        {
            "support_mode": row["support_mode"],
            "metric": row["metric"],
            "expected_direction": row["expected_direction"],
            "spearman_rho": round(float(row["spearman_rho"]), 6),
            "direction_pass": row["direction_pass"],
            "n_parameter_sets": int(row["n_parameter_sets"]),
        }
        for row in monotonicity
    ]
    write_csv("run-127_monotonicity_summary.csv", mono_rows, mono_fields)
    write_md_table("run-127_monotonicity_summary.md", mono_rows, mono_fields)

    metric_fields = [
        "n_local_maxima",
        "accessible_path_fraction",
        "walsh_high_order_fraction",
        "sign_epistasis_fraction",
        "reciprocal_sign_epistasis_fraction",
        "magnitude_epistasis_to_single_mutant_sd_ratio",
        "random_walk_lag1_autocorrelation",
    ]
    epi_rows = []
    for (epi_weight,), group in sorted(grouped(knob, "epi_weight").items(), key=lambda item: float(item[0][0])):
        row = {"epi_weight": float(epi_weight), "n_parameter_sets": len(group)}
        for field in metric_fields:
            values = [as_float(item, field) for item in group]
            row[f"{field}_mean"] = mean(values)
            row[f"{field}_min"] = min(values)
            row[f"{field}_max"] = max(values)
        epi_rows.append(row)
    epi_fields = ["epi_weight", "n_parameter_sets"]
    for field in metric_fields:
        epi_fields.extend([f"{field}_mean", f"{field}_min", f"{field}_max"])
    write_csv("run-127_epi_weight_metric_summary.csv", epi_rows, epi_fields)
    write_md_table("run-127_epi_weight_metric_summary.md", epi_rows, epi_fields)

    mode_rows = []
    for (mode, epi_weight), group in sorted(
        grouped(knob, "support_mode", "epi_weight").items(),
        key=lambda item: (item[0][0], float(item[0][1])),
    ):
        row = {"support_mode": mode, "epi_weight": float(epi_weight), "n_parameter_sets": len(group)}
        for field in metric_fields:
            row[f"{field}_mean"] = mean(as_float(item, field) for item in group)
        row["joint_multipeak_navigable_rows"] = sum(
            as_float(item, "n_local_maxima") > 1.0 and as_float(item, "accessible_path_fraction") > 0.0
            for item in group
        )
        mode_rows.append(row)
    mode_fields = ["support_mode", "epi_weight", "n_parameter_sets"]
    mode_fields.extend(f"{field}_mean" for field in metric_fields)
    mode_fields.append("joint_multipeak_navigable_rows")
    write_csv("run-127_support_mode_knob_summary.csv", mode_rows, mode_fields)
    write_md_table("run-127_support_mode_knob_summary.md", mode_rows, mode_fields)

    joint_rows = []
    for epi_weight in sorted({float(row["epi_weight"]) for row in knob}):
        group = [row for row in knob if float(row["epi_weight"]) == epi_weight]
        joint = [
            row
            for row in group
            if as_float(row, "n_local_maxima") > 1.0 and as_float(row, "accessible_path_fraction") > 0.0
        ]
        joint_rows.append(
            {
                "epi_weight": epi_weight,
                "joint_rows": len(joint),
                "n_parameter_sets": len(group),
                "mean_local_maxima": mean(as_float(row, "n_local_maxima") for row in group),
                "mean_accessible_path_fraction": mean(as_float(row, "accessible_path_fraction") for row in group),
                "mean_walsh_high_order_fraction": mean(as_float(row, "walsh_high_order_fraction") for row in group),
                "mean_sign_epistasis_fraction": mean(as_float(row, "sign_epistasis_fraction") for row in group),
            }
        )
    joint_fields = [
        "epi_weight",
        "joint_rows",
        "n_parameter_sets",
        "mean_local_maxima",
        "mean_accessible_path_fraction",
        "mean_walsh_high_order_fraction",
        "mean_sign_epistasis_fraction",
    ]
    write_csv("run-127_joint_regime_summary.csv", joint_rows, joint_fields)
    write_md_table("run-127_joint_regime_summary.md", joint_rows, joint_fields)

    core_rows = [
        {
            "n_core_rows": len(core),
            "beneficial_fraction_declines_rows": sum(as_bool(row, "beneficial_fraction_declines_toward_peak") for row in core),
            "plateau_rows": sum(as_bool(row, "has_low_and_high_readout_plateaus") for row in core),
            "transition_epistasis_enrichment_rows": sum(
                as_bool(row, "transition_epistasis_enriched_over_deep_regimes") for row in core
            ),
            "mean_beneficial_fraction_spearman": mean(
                as_float(row, "beneficial_fraction_vs_fraction_of_peak_spearman") for row in core
            ),
            "mean_high_plateau_fraction": mean(as_float(row, "high_plateau_fraction") for row in core),
            "mean_low_plateau_fraction": mean(as_float(row, "low_plateau_fraction") for row in core),
        }
    ]
    core_fields = list(core_rows[0])
    write_csv("run-127_core_readout_guardrail_summary.csv", core_rows, core_fields)
    write_md_table("run-127_core_readout_guardrail_summary.md", core_rows, core_fields)

    hamming_epi0 = [
        row
        for row in hamming
        if row["score_layer"] == "observable_fitness" and abs(as_float(row, "epi_weight")) < 1e-12
    ]
    hamming_rows = []
    for (distance,), group in sorted(grouped(hamming_epi0, "hamming_distance_from_peak").items(), key=lambda item: int(item[0][0])):
        hamming_rows.append(
            {
                "hamming_distance_from_peak": int(distance),
                "fraction_of_peak": float(group[0]["fraction_of_peak"]),
                "mean_beneficial_fraction": mean(as_float(row, "beneficial_fraction") for row in group),
                "mean_beneficial_effect": mean(as_float(row, "mean_beneficial_effect") for row in group),
                "n_rows": len(group),
            }
        )
    hamming_fields = [
        "hamming_distance_from_peak",
        "fraction_of_peak",
        "mean_beneficial_fraction",
        "mean_beneficial_effect",
        "n_rows",
    ]
    write_csv("run-127_epi0_hamming_beneficial_summary.csv", hamming_rows, hamming_fields)
    write_md_table("run-127_epi0_hamming_beneficial_summary.md", hamming_rows, hamming_fields)

    threshold_obs = [row for row in threshold if row["score_layer"] == "observable_fitness"]
    threshold_rows = []
    for (epi_weight, band), group in sorted(
        grouped(threshold_obs, "epi_weight", "threshold_band").items(),
        key=lambda item: (float(item[0][0]), item[0][1]),
    ):
        threshold_rows.append(
            {
                "epi_weight": float(epi_weight),
                "threshold_band": band,
                "n_rows": len(group),
                "mean_abs_pairwise_epistasis": mean(as_float(row, "mean_abs_pairwise_epistasis") for row in group),
                "mean_sign_epistasis_fraction": mean(as_float(row, "sign_epistasis_fraction") for row in group),
                "mean_reciprocal_sign_epistasis_fraction": mean(
                    as_float(row, "reciprocal_sign_epistasis_fraction") for row in group
                ),
            }
        )
    threshold_fields = [
        "epi_weight",
        "threshold_band",
        "n_rows",
        "mean_abs_pairwise_epistasis",
        "mean_sign_epistasis_fraction",
        "mean_reciprocal_sign_epistasis_fraction",
    ]
    write_csv("run-127_threshold_epistasis_summary.csv", threshold_rows, threshold_fields)
    write_md_table("run-127_threshold_epistasis_summary.md", threshold_rows, threshold_fields)

    support_rows = []
    for (mode, order), group in sorted(grouped(support, "support_mode", "order").items(), key=lambda item: (item[0][0], int(item[0][1]))):
        support_rows.append(
            {
                "support_mode": mode,
                "order": int(order),
                "n_terms": len(group),
                "n_unique_supports": len({row["support_positions"] for row in group}),
            }
        )
    support_fields = ["support_mode", "order", "n_terms", "n_unique_supports"]
    write_csv("run-127_support_order_summary.csv", support_rows, support_fields)
    write_md_table("run-127_support_order_summary.md", support_rows, support_fields)

    all_direction_passed = all(as_bool(row, "direction_pass") for row in monotonicity)
    core_decline_passed = all(as_bool(row, "beneficial_fraction_declines_toward_peak") for row in core)
    core_plateau_passed = all(as_bool(row, "has_low_and_high_readout_plateaus") for row in core)
    core_transition_passed = all(as_bool(row, "transition_epistasis_enriched_over_deep_regimes") for row in core)
    epi06_rows = [row for row in knob if abs(as_float(row, "epi_weight") - 0.6) < 1e-12]
    epi1_rows = [row for row in knob if abs(as_float(row, "epi_weight") - 1.0) < 1e-12]
    metrics = {
        **summary,
        "run_id": "RUN-127",
        "scheduler_job_id": "1935",
        "scheduler_exit_code": "0",
        "copied_source_payload_bytes": copied_bytes,
        "n_copied_source_files": len(audit_rows),
        "all_direction_checks_passed": all_direction_passed,
        "core_beneficial_decline_passed": core_decline_passed,
        "core_plateau_guardrail_passed": core_plateau_passed,
        "core_transition_epistasis_enrichment_passed": core_transition_passed,
        "n_core_transition_epistasis_enrichment_rows": sum(
            as_bool(row, "transition_epistasis_enriched_over_deep_regimes") for row in core
        ),
        "joint_multipeak_navigable_rows": sum(
            as_float(row, "n_local_maxima") > 1.0 and as_float(row, "accessible_path_fraction") > 0.0
            for row in knob
        ),
        "joint_multipeak_navigable_rows_epi_weight_0_6": sum(
            as_float(row, "n_local_maxima") > 1.0 and as_float(row, "accessible_path_fraction") > 0.0
            for row in epi06_rows
        ),
        "mean_local_maxima_epi_weight_0": mean(as_float(row, "n_local_maxima") for row in knob if as_float(row, "epi_weight") == 0.0),
        "mean_accessible_path_fraction_epi_weight_0": mean(
            as_float(row, "accessible_path_fraction") for row in knob if as_float(row, "epi_weight") == 0.0
        ),
        "mean_local_maxima_epi_weight_0_6": mean(as_float(row, "n_local_maxima") for row in epi06_rows),
        "mean_accessible_path_fraction_epi_weight_0_6": mean(
            as_float(row, "accessible_path_fraction") for row in epi06_rows
        ),
        "mean_walsh_high_order_fraction_epi_weight_0_6": mean(
            as_float(row, "walsh_high_order_fraction") for row in epi06_rows
        ),
        "mean_sign_epistasis_fraction_epi_weight_0_6": mean(
            as_float(row, "sign_epistasis_fraction") for row in epi06_rows
        ),
        "mean_local_maxima_epi_weight_1": mean(as_float(row, "n_local_maxima") for row in epi1_rows),
        "mean_accessible_path_fraction_epi_weight_1": mean(
            as_float(row, "accessible_path_fraction") for row in epi1_rows
        ),
        "mean_walsh_high_order_fraction_epi_weight_1": mean(
            as_float(row, "walsh_high_order_fraction") for row in epi1_rows
        ),
        "mean_sign_epistasis_fraction_epi_weight_1": mean(
            as_float(row, "sign_epistasis_fraction") for row in epi1_rows
        ),
        "min_direction_spearman_abs": min(abs(float(row["spearman_rho"])) for row in monotonicity),
    }
    (TABLE_DIR / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
