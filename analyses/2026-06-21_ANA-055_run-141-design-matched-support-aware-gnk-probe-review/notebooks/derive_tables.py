from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ANALYSIS_DIR.parents[1]
RUN_DIR = PROJECT_ROOT / "data/processed/stage-4b-design-matched-support-aware-gnk-probe/RUN-141"
OUTPUT_DIR = RUN_DIR / "outputs"
LOG_DIR = RUN_DIR / "logs"
TABLE_DIR = ANALYSIS_DIR / "tables"
SOURCE_DIR = ANALYSIS_DIR / "source_tables"


def write_markdown_table(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        path.write_text("(empty)\n")
        return
    rendered = df.copy()
    for col in rendered.columns:
        if pd.api.types.is_float_dtype(rendered[col]):
            rendered[col] = rendered[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
    headers = list(rendered.columns)
    rows = rendered.astype(str).values.tolist()
    widths = [
        max(len(str(header)), *(len(row[i]) for row in rows))
        for i, header in enumerate(headers)
    ]
    lines = [
        "| " + " | ".join(str(header).ljust(widths[i]) for i, header in enumerate(headers)) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(headers))) + " |")
    path.write_text("\n".join(lines) + "\n")


def save_table(df: pd.DataFrame, stem: str) -> None:
    csv_path = TABLE_DIR / f"{stem}.csv"
    md_path = TABLE_DIR / f"{stem}.md"
    df.to_csv(csv_path, index=False)
    write_markdown_table(df, md_path)


def grouped_pass_fraction(df: pd.DataFrame, group_cols: list[str], pass_col: str) -> pd.DataFrame:
    return (
        df.groupby(group_cols, dropna=False)
        .agg(
            n=(pass_col, "size"),
            pass_n=(pass_col, "sum"),
            pass_fraction=(pass_col, "mean"),
            median_spearman_rho=("spearman_rho", "median")
            if "spearman_rho" in df.columns
            else (pass_col, "mean"),
        )
        .reset_index()
    )


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)

    summary = json.loads((OUTPUT_DIR / "summary.json").read_text())
    preflight = json.loads((OUTPUT_DIR / "preflight.json").read_text())
    checkpoint = json.loads((OUTPUT_DIR / "checkpoint_resume_summary.json").read_text())

    for name in ["summary.json", "preflight.json", "checkpoint_resume_summary.json"]:
        (SOURCE_DIR / f"run-141_{name}").write_text((OUTPUT_DIR / name).read_text())

    output_rows = []
    for path in sorted(OUTPUT_DIR.glob("*")):
        if path.is_file():
            output_rows.append(
                {
                    "relative_path": str(path.relative_to(PROJECT_ROOT)),
                    "bytes": path.stat().st_size,
                    "record_role": "run_output",
                }
            )
    for path in sorted(LOG_DIR.glob("*")):
        if path.is_file():
            output_rows.append(
                {
                    "relative_path": str(path.relative_to(PROJECT_ROOT)),
                    "bytes": path.stat().st_size,
                    "record_role": "run_log",
                }
            )
    output_audit = pd.DataFrame(output_rows)
    save_table(output_audit, "run-141_output_collection_audit")

    mono = pd.read_csv(OUTPUT_DIR / "generator_knob_monotonicity.csv")
    mono_by_metric = grouped_pass_fraction(mono, ["metric", "expected_direction"], "direction_pass")
    mono_by_metric = mono_by_metric.sort_values(["pass_fraction", "metric"], ascending=[False, True])
    save_table(mono_by_metric, "run-141_monotonicity_by_metric")

    mono_by_support_anchor = grouped_pass_fraction(
        mono, ["support_mode", "anchor_regime"], "direction_pass"
    ).sort_values(["support_mode", "anchor_regime"])
    save_table(mono_by_support_anchor, "run-141_monotonicity_by_support_anchor")

    mono_by_dataset = grouped_pass_fraction(mono, ["dataset_id"], "direction_pass")
    mono_by_dataset = mono_by_dataset.sort_values(["pass_fraction", "dataset_id"], ascending=[False, True])
    save_table(mono_by_dataset, "run-141_monotonicity_by_dataset")

    sep = pd.read_csv(OUTPUT_DIR / "core_gnk_separation_gates.csv")
    positive_sep = sep[sep["epi_weight"] > 0].copy()
    gate_cols = [
        "high_order_fraction_separation_pass",
        "sign_epistasis_fraction_separation_pass",
        "local_maxima_density_separation_pass",
        "core_gnk_separation_pass",
    ]
    core_gnk_by_metric = pd.DataFrame(
        [
            {
                "gate": col.replace("_separation_pass", "").replace("core_gnk", "any_preregistered_metric"),
                "n_positive_epi_weight_rows": len(positive_sep),
                "pass_n": int(positive_sep[col].sum()),
                "pass_fraction": float(positive_sep[col].mean()),
            }
            for col in gate_cols
        ]
    )
    save_table(core_gnk_by_metric, "run-141_core_gnk_gate_by_metric")

    core_gnk_by_epi = (
        positive_sep.groupby("epi_weight", dropna=False)
        .agg(
            n=("core_gnk_separation_pass", "size"),
            pass_n=("core_gnk_separation_pass", "sum"),
            pass_fraction=("core_gnk_separation_pass", "mean"),
            mean_delta_high_order=("high_order_fraction_full_minus_core", "mean"),
            mean_abs_delta_high_order=(
                "high_order_fraction_full_minus_core",
                lambda x: x.abs().mean(),
            ),
            mean_delta_sign=("sign_epistasis_fraction_full_minus_core", "mean"),
            mean_delta_local_maxima=("local_maxima_density_full_minus_core", "mean"),
        )
        .reset_index()
    )
    save_table(core_gnk_by_epi, "run-141_core_gnk_gate_by_epi_weight")

    core_gnk_by_dataset = (
        positive_sep.groupby("dataset_id", dropna=False)
        .agg(
            n=("core_gnk_separation_pass", "size"),
            pass_n=("core_gnk_separation_pass", "sum"),
            pass_fraction=("core_gnk_separation_pass", "mean"),
            mean_abs_delta_high_order=(
                "high_order_fraction_full_minus_core",
                lambda x: x.abs().mean(),
            ),
            max_abs_delta_high_order=(
                "high_order_fraction_full_minus_core",
                lambda x: x.abs().max(),
            ),
        )
        .reset_index()
        .sort_values(["pass_fraction", "dataset_id"], ascending=[False, True])
    )
    save_table(core_gnk_by_dataset, "run-141_core_gnk_gate_by_dataset")

    support = pd.read_csv(OUTPUT_DIR / "support_aware_graphlet_audit.csv")
    templates = pd.read_csv(OUTPUT_DIR / "empirical_design_templates.csv")
    empirical_descriptor_status = pd.read_csv(OUTPUT_DIR / "empirical_template_descriptor_metrics.csv")
    support_summary = (
        support.groupby(["dataset_id", "dimension", "alphabet_cardinalities"], dropna=False)
        .agg(
            n_graphlets=("graphlet_id", "nunique"),
            support_aware_fraction=("support_aware", "mean"),
            mean_node_count=("n_nodes", "mean"),
        )
        .reset_index()
    )
    template_audit = support_summary.merge(
        templates[
            [
                "dataset_id",
                "n_expected_nodes",
                "n_observed_unique_states",
                "observed_state_fraction",
                "complete_factorial",
                "alphabet",
            ]
        ],
        on="dataset_id",
        how="left",
    ).merge(
        empirical_descriptor_status[["dataset_id", "descriptor_status", "skip_reason"]],
        on="dataset_id",
        how="left",
    )
    save_table(template_audit, "run-141_support_template_audit")

    anchors = pd.read_csv(OUTPUT_DIR / "anchor_regime_audit.csv")
    anchor_summary = (
        anchors.groupby("anchor_regime", dropna=False)
        .agg(
            n=("actual_effective_stability_margin", "size"),
            target_margin=("target_effective_stability_margin", "median"),
            actual_margin_mean=("actual_effective_stability_margin", "mean"),
            actual_margin_sd=("actual_effective_stability_margin", "std"),
            mean_abs_margin_error=(
                "actual_effective_stability_margin",
                lambda x: (
                    x - anchors.loc[x.index, "target_effective_stability_margin"]
                ).abs().mean(),
            ),
            mean_distance_from_designed_peak=("distance_from_designed_peak", "mean"),
            max_distance_from_designed_peak=("distance_from_designed_peak", "max"),
        )
        .reset_index()
    )
    save_table(anchor_summary, "run-141_anchor_regime_audit")

    gen_cols = [
        "model_family",
        "dataset_id",
        "support_mode",
        "anchor_regime",
        "stability_margin",
        "epi_weight",
        "local_maxima_density",
        "peak_targeted_reachability_fraction",
        "random_walk_lag1_autocorrelation",
        "high_order_fraction",
        "sign_epistasis_fraction",
        "reciprocal_sign_epistasis_fraction",
    ]
    gen = pd.read_csv(OUTPUT_DIR / "generator_descriptor_metrics.csv", usecols=gen_cols)
    wide = gen[gen["model_family"].isin(["full_public_observable", "additive_core_only_observable"])].pivot_table(
        index=["dataset_id", "support_mode", "anchor_regime", "stability_margin", "epi_weight"],
        columns="model_family",
        values=[
            "local_maxima_density",
            "peak_targeted_reachability_fraction",
            "random_walk_lag1_autocorrelation",
            "high_order_fraction",
            "sign_epistasis_fraction",
            "reciprocal_sign_epistasis_fraction",
        ],
        aggfunc="mean",
    )
    wide.columns = ["__".join(col) for col in wide.columns]
    wide = wide.reset_index()
    delta_rows = []
    metrics = [
        "high_order_fraction",
        "sign_epistasis_fraction",
        "reciprocal_sign_epistasis_fraction",
        "local_maxima_density",
        "peak_targeted_reachability_fraction",
        "random_walk_lag1_autocorrelation",
    ]
    for metric in metrics:
        full = f"{metric}__full_public_observable"
        core = f"{metric}__additive_core_only_observable"
        tmp = wide[wide["epi_weight"] > 0].copy()
        tmp["delta"] = tmp[full] - tmp[core]
        grouped = (
            tmp.groupby("epi_weight", dropna=False)["delta"]
            .agg(
                mean_delta="mean",
                median_delta="median",
                mean_abs_delta=lambda x: x.abs().mean(),
                max_abs_delta=lambda x: x.abs().max(),
            )
            .reset_index()
        )
        grouped.insert(0, "metric", metric)
        delta_rows.append(grouped)
    full_core_delta = pd.concat(delta_rows, ignore_index=True)
    save_table(full_core_delta, "run-141_full_core_deltas_by_epi_weight")

    metrics_json = {
        "run": "RUN-141",
        "experiment": "EXP-057",
        "mode": summary["mode"],
        "scheduler_exit_code": int((LOG_DIR / "labproj-exit-code.txt").read_text().strip()),
        "stderr_bytes": (LOG_DIR / "stderr.txt").stat().st_size,
        "stdout_lines": len((LOG_DIR / "stdout.txt").read_text().splitlines()),
        "n_parameter_sets": summary["n_parameter_sets"],
        "n_computed_parameter_sets": summary["n_computed_parameter_sets"],
        "n_resumed_parameter_sets": summary["n_resumed_parameter_sets"],
        "n_generator_metric_rows": summary["n_generator_metric_rows"],
        "n_templates": summary["n_templates"],
        "n_support_graphlets": int(support["graphlet_id"].nunique()),
        "support_aware_fraction": float(support["support_aware"].mean()),
        "n_empirical_template_descriptor_rows": summary["n_empirical_template_descriptor_rows"],
        "n_empirical_template_descriptors_computed": summary["n_empirical_template_descriptors_computed"],
        "monotonicity_checks": summary["n_monotonicity_checks"],
        "monotonicity_direction_pass": summary["n_monotonicity_direction_pass"],
        "monotonicity_pass_fraction": summary["monotonicity_pass_fraction"],
        "positive_core_gnk_gate_rows": summary["n_positive_core_gnk_gate_rows"],
        "positive_core_gnk_gate_pass": summary["n_positive_core_gnk_gate_pass"],
        "core_gnk_gate_pass_fraction": summary["core_gnk_gate_pass_fraction"],
        "hard_gate_pass": False,
        "downstream_overlap_joint_heldout_allowed": False,
        "stop_rule": summary["stop_rule"],
        "dvc_outputs_pointer": str(
            RUN_DIR.joinpath("outputs.dvc").relative_to(PROJECT_ROOT)
        ),
        "dvc_logs_pointer": str(RUN_DIR.joinpath("logs.dvc").relative_to(PROJECT_ROOT)),
        "remote_output_dir": "/tank/data/labproj/runtime/project-adapt-env/runs/RUN-141/outputs",
        "checkpoint_root": summary["checkpoint_root"],
        "scratch_design_source": preflight["scratch_design_source"],
        "computed_checkpoint_count": len(checkpoint.get("computed_parameter_sets", [])),
        "resumed_checkpoint_count": len(checkpoint.get("resumed_parameter_sets", [])),
    }
    (TABLE_DIR / "metrics.json").write_text(json.dumps(metrics_json, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
