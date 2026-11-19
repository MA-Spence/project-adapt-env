from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


ANA_DIR = Path(__file__).resolve().parents[1]
SOURCE_DIR = ANA_DIR / "source_tables"
TABLE_DIR = ANA_DIR / "tables"


def read_csv(name: str) -> list[dict[str, str]]:
    path = SOURCE_DIR / name
    if path.stat().st_size == 0 or path.read_text().strip() == "":
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f(value: str | float | int | None) -> float:
    if value in (None, ""):
        return math.nan
    return float(value)


def bool_str(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() == "true"


def fmt(value: float, digits: int = 4) -> str:
    if value is None or math.isnan(value):
        return "NA"
    if value == 0:
        return "0"
    if abs(value) < 0.001 or abs(value) >= 1000:
        return f"{value:.3e}"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |"]
    out.append("| " + " | ".join(["---"] * len(headers)) + " |")
    out.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(out) + "\n"


def write_text(name: str, text: str) -> None:
    (TABLE_DIR / name).write_text(text)


def main() -> None:
    TABLE_DIR.mkdir(exist_ok=True)
    summary = json.loads((SOURCE_DIR / "run-123_summary.json").read_text())
    spectrum = read_csv("run-123_sparse_walsh_spectrum_calibration.csv")
    probe = read_csv("run-123_sparse_walsh_forward_probe.csv")
    comparator = read_csv("run-123_abc_comparator_selected_family_distances.csv")
    graphlets = read_csv("run-123_graphlet_catalog.csv")
    modules = read_csv("run-123_module_catalog.csv")
    manifest_rows = []
    manifest_path = SOURCE_DIR / "run-123_remote_output_manifest.tsv"
    if manifest_path.exists():
        for line in manifest_path.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) == 3:
                manifest_rows.append(
                    {"file": parts[0], "bytes": int(parts[1]), "sha256": parts[2]}
                )

    passing_orders = [
        int(row["order"]) for row in spectrum if bool_str(row["passes_noise_floor"])
    ]
    order_weights = {
        int(row["order"]): f(row["global_order_weight_normalized"]) for row in spectrum
    }

    probe_metrics: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "n": 0,
            "passes": 0,
            "empirical": math.nan,
            "model_means": [],
            "qlows": [],
            "qhighs": [],
        }
    )
    probe_candidates: dict[int, dict[str, object]] = defaultdict(
        lambda: {
            "n": 0,
            "passes": 0,
            "epi_weight": math.nan,
            "high_order_weight_multiplier": math.nan,
        }
    )
    for row in probe:
        metric = row["metric"]
        metric_info = probe_metrics[metric]
        metric_info["n"] = int(metric_info["n"]) + 1
        metric_info["passes"] = int(metric_info["passes"]) + int(
            bool_str(row["metric_pass"])
        )
        metric_info["empirical"] = f(row["empirical_mean"])
        metric_info["model_means"].append(f(row["model_mean"]))
        metric_info["qlows"].append(f(row["model_central_q_low"]))
        metric_info["qhighs"].append(f(row["model_central_q_high"]))

        candidate_idx = int(row["probe_candidate_index"])
        candidate_info = probe_candidates[candidate_idx]
        candidate_info["n"] = int(candidate_info["n"]) + 1
        candidate_info["passes"] = int(candidate_info["passes"]) + int(
            bool_str(row["metric_pass"])
        )
        candidate_info["epi_weight"] = f(row["epi_weight"])
        candidate_info["high_order_weight_multiplier"] = f(
            row["high_order_weight_multiplier"]
        )

    failed_metrics_all_candidates = sorted(
        metric
        for metric, info in probe_metrics.items()
        if int(info["passes"]) == 0 and int(info["n"]) > 0
    )
    joint_pass_rows = sum(1 for row in probe if bool_str(row["joint_pass"]))
    n_probe_candidates = len(probe_candidates)
    n_probe_metrics = len(probe_metrics)

    test_comparator = [
        row for row in comparator if row.get("graphlet_split") == "test"
    ]
    test_comparator.sort(key=lambda row: f(row["distance_to_empirical_primary"]))
    best_available = test_comparator[0] if test_comparator else {}
    best_nested = next(
        (
            row
            for row in test_comparator
            if row.get("candidate_family") != "empirical_readout_residual"
        ),
        {},
    )

    stage_rows = []
    for row in summary["fix_stage_summaries"]:
        stage_rows.append(
            [
                str(row["fix_stage"]),
                str(row["fixes_applied"]),
                str(row["fast_fit_enabled"]),
                str(row["sparse_walsh_forward_probe_passed"]),
                str(row["n_candidate_evaluations"]),
                str(row["n_selected_candidates"]),
                fmt(float(row["best_objective_distance"]))
                if not math.isnan(float(row["best_objective_distance"]))
                else "NA",
            ]
        )
    write_text(
        "run-123_stage_gate_and_fit_decision.md",
        "# RUN-123 Stage Gate And Fit Decision\n\n"
        + markdown_table(
            [
                "stage",
                "fixes",
                "fast_fit_enabled",
                "H1 passed",
                "candidate evaluations",
                "selected candidates",
                "best objective",
            ],
            stage_rows,
        )
        + "\nThe sequence-conditioned sparse-Walsh FIT stage produced no candidate "
        "rows because the preregistered H1 forward gate did not pass.\n",
    )

    spectrum_rows = []
    for row in spectrum:
        spectrum_rows.append(
            [
                row["order"],
                row["n_graphlets"],
                fmt(f(row["local_variance_mean"])),
                fmt(f(row["local_variance_bootstrap_q10"])),
                fmt(f(row["local_variance_bootstrap_q90"])),
                fmt(f(row["local_fraction_mean"])),
                row["passes_noise_floor"],
                fmt(f(row["support_hit_probability"])),
                fmt(f(row["global_sigma2_unscaled"])),
                fmt(f(row["global_order_weight_normalized"])),
            ]
        )
    write_text(
        "run-123_sparse_walsh_spectrum_calibration.md",
        "# RUN-123 Sparse Walsh Spectrum Calibration\n\n"
        + markdown_table(
            [
                "order",
                "calibration graphlets",
                "local variance mean",
                "bootstrap q10",
                "bootstrap q90",
                "local fraction mean",
                "passes floor",
                "support hit probability",
                "global sigma2 unscaled",
                "normalized global weight",
            ],
            spectrum_rows,
        )
        + "\nOnly the preregistered pass flag is used for promotion to H1 weights; "
        "here it leaves all normalized global weight on order 3.\n",
    )

    metric_rows = []
    for metric, info in sorted(probe_metrics.items()):
        model_means = info["model_means"]
        qlows = info["qlows"]
        qhighs = info["qhighs"]
        empirical = float(info["empirical"])
        envelope_low = min(qlows)
        envelope_high = max(qhighs)
        if empirical < envelope_low:
            interpretation = "empirical below model interval"
        elif empirical > envelope_high:
            interpretation = "empirical above model interval"
        else:
            interpretation = "empirical inside at least one interval"
        metric_rows.append(
            [
                metric,
                f"{info['passes']}/{info['n']}",
                fmt(empirical),
                f"{fmt(min(model_means))}-{fmt(max(model_means))}",
                f"{fmt(envelope_low)}-{fmt(envelope_high)}",
                interpretation,
            ]
        )
    write_text(
        "run-123_forward_probe_metric_gate.md",
        "# RUN-123 Forward Sparse-Walsh H1 Metric Gate\n\n"
        + markdown_table(
            [
                "metric",
                "pass count",
                "empirical mean",
                "model mean range",
                "central interval envelope",
                "interpretation",
            ],
            metric_rows,
        )
        + "\nThe joint gate requires all target metrics for a probe candidate to pass. "
        "No candidate passed jointly.\n",
    )

    candidate_rows = []
    for idx, info in sorted(probe_candidates.items()):
        candidate_rows.append(
            [
                str(idx),
                fmt(float(info["epi_weight"])),
                fmt(float(info["high_order_weight_multiplier"])),
                f"{info['passes']}/{info['n']}",
                "False",
            ]
        )
    write_text(
        "run-123_forward_probe_candidate_gate.md",
        "# RUN-123 Forward Sparse-Walsh H1 Candidate Gate\n\n"
        + markdown_table(
            [
                "probe candidate",
                "epi weight",
                "high-order multiplier",
                "metric passes",
                "joint pass",
            ],
            candidate_rows,
        )
        + "\nAll 12 H1 candidates passed exactly four of seven metrics and failed "
        "the joint gate.\n",
    )

    comparator_rows = []
    for row in test_comparator[:10]:
        comparator_rows.append(
            [
                row["candidate_family"],
                row["candidate_model"],
                row["parameter_set_id"],
                fmt(f(row["distance_to_empirical_primary"])),
            ]
        )
    write_text(
        "run-123_comparator_context.md",
        "# RUN-123 Comparator Context\n\n"
        + markdown_table(
            ["candidate family", "candidate model", "parameter set", "test distance"],
            comparator_rows,
        )
        + "\nRUN-123 generated no sparse-Walsh FIT candidate rows, so these values "
        "are context only; no operational equivalence comparison was possible.\n",
    )

    split_counts = Counter(row["graphlet_split"] for row in graphlets)
    dim_counts = Counter(row["graphlet_dimension"] for row in graphlets)
    module_source_counts = Counter(row["module_source"] for row in modules)
    catalog_rows = [
        ["graphlets total", str(len(graphlets))],
        ["graphlets by split", ", ".join(f"{k}: {v}" for k, v in sorted(split_counts.items()))],
        ["graphlets by dimension", ", ".join(f"{k}: {v}" for k, v in sorted(dim_counts.items()))],
        ["modules total", str(len(modules))],
        ["modules by source", ", ".join(f"{k}: {v}" for k, v in sorted(module_source_counts.items()))],
    ]
    write_text(
        "run-123_graphlet_and_module_catalog_summary.md",
        "# RUN-123 Graphlet And Module Catalog Summary\n\n"
        + markdown_table(["quantity", "value"], catalog_rows),
    )

    copied_files = sorted(SOURCE_DIR.glob("run-123_*"))
    copied_bytes = sum(path.stat().st_size for path in copied_files)
    manifest_table_rows = []
    for row in sorted(manifest_rows, key=lambda row: row["file"]):
        copied = SOURCE_DIR / f"run-123_{row['file']}"
        manifest_table_rows.append(
            [
                row["file"],
                str(row["bytes"]),
                row["sha256"][:16],
                "yes" if copied.exists() else "no",
            ]
        )
    write_text(
        "run-123_output_collection_audit.md",
        "# RUN-123 Output Collection Audit\n\n"
        f"Copied source payload: `{copied_bytes}` bytes across `{len(copied_files)}` files.\n\n"
        + markdown_table(
            ["remote file", "remote bytes", "sha256 prefix", "copied to ANA"],
            manifest_table_rows,
        )
        + "\nThe PHOT alignment FASTA was not copied because it is not needed for "
        "the gate decision; it remains on the runtime volume.\n",
    )

    metrics = {
        "run_id": "RUN-123",
        "experiment_id": "EXP-048",
        "scheduler_job_id": "1930",
        "scheduler_exit_code": "0",
        "n_graphlets": int(summary["n_graphlets"]),
        "n_sparse_walsh_spectrum_rows": int(summary["n_sparse_walsh_spectrum_rows"]),
        "n_sparse_walsh_forward_probe_rows": int(
            summary["n_sparse_walsh_forward_probe_rows"]
        ),
        "sparse_walsh_forward_probe_passed": bool(
            summary["sparse_walsh_forward_probe_passed"]
        ),
        "n_probe_candidates": n_probe_candidates,
        "n_probe_metrics": n_probe_metrics,
        "n_joint_pass_rows": joint_pass_rows,
        "failed_metrics_all_candidates": failed_metrics_all_candidates,
        "passing_spectrum_orders": passing_orders,
        "global_order_weights": order_weights,
        "n_candidate_evaluations": int(summary["n_candidate_evaluations"]),
        "n_selected_candidates": int(summary["n_selected_candidates"]),
        "n_equivalence_rows": int(summary["n_equivalence_rows"]),
        "n_operationally_equivalent": int(summary["n_operationally_equivalent"]),
        "best_available_run098_test_distance": f(
            best_available.get("distance_to_empirical_primary")
        )
        if best_available
        else math.nan,
        "best_available_run098_family": best_available.get("candidate_family"),
        "best_nested_run098_test_distance": f(
            best_nested.get("distance_to_empirical_primary")
        )
        if best_nested
        else math.nan,
        "best_nested_run098_family": best_nested.get("candidate_family"),
        "copied_source_payload_bytes": copied_bytes,
    }
    (TABLE_DIR / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
