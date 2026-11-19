from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd


RUN_ID = "RUN-129"
ANALYSIS_ID = "ANA-056"
RUN_SLUG = "run-129"
DATASET = "phot-nonadditive-measurement-operator-certification"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_table(df: pd.DataFrame, csv_path: Path, md_path: Path | None = None) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    if md_path is not None:
        md_path.write_text(df.to_markdown(index=False) + "\n")


def short_error(text: object) -> str:
    if pd.isna(text):
        return ""
    value = str(text)
    if "pairwise_gpmap_preflight_skip" in value:
        return "pairwise_gpmap_preflight_skip: sequence_length=118 exceeds max_pairwise_sequence_length=80"
    if "MultilayerPerceptronGPMap" in value and "build()" in value:
        return "blackbox_ge_neural MAVE-NN/Keras build signature error"
    return value[:180]


def bool_series(series: pd.Series) -> pd.Series:
    return series.map(lambda value: False if pd.isna(value) else bool(value))


def classify_failure(rows: pd.DataFrame) -> str:
    errors = " ".join(str(x) for x in rows["mavenn_fit_error"].dropna().tolist())
    skipped = bool_series(rows["mavenn_fit_skipped_preflight"])
    if skipped.all():
        return "preflight_skipped_pairwise_length"
    if "MultilayerPerceptronGPMap" in errors and "build()" in errors:
        return "mavenn_keras_blackbox_build_api_error"
    if bool_series(rows["fit_failed"]).all():
        return "fit_failed_other"
    return "fit_completed"


def main() -> None:
    analysis_dir = Path(__file__).resolve().parents[1]
    project_root = analysis_dir.parents[1]
    output_dir = project_root / "data" / "processed" / DATASET / RUN_ID / "outputs"
    log_dir = project_root / "data" / "processed" / DATASET / RUN_ID / "logs"
    source_dir = analysis_dir / "source_tables"
    tables_dir = analysis_dir / "tables"
    source_dir.mkdir(exist_ok=True)
    tables_dir.mkdir(exist_ok=True)

    expected_outputs = {
        "preflight.json",
        "graphlet_catalog.csv",
        "operator_profile_catalog.csv",
        "injected_terms.csv",
        "roundtrip_graphlet_metrics.csv",
        "roundtrip_transfer.csv",
        "mavenn_metrics.csv",
        "operator_certification_summary.csv",
        "summary.json",
    }
    missing = sorted(name for name in expected_outputs if not (output_dir / name).exists())
    if missing:
        raise FileNotFoundError(f"Missing expected {RUN_ID} outputs: {missing}")

    audit_rows: list[dict[str, object]] = []
    for root, role in [(output_dir, "run_output"), (log_dir, "run_log")]:
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(project_root)
            audit_rows.append(
                {
                    "run_id": RUN_ID,
                    "role": role,
                    "source_path": str(rel),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
            if role == "run_output":
                destination = source_dir / f"{RUN_SLUG}_{path.name}"
                shutil.copy2(path, destination)

    audit = pd.DataFrame(audit_rows)
    write_table(
        audit,
        tables_dir / f"{RUN_SLUG}_output_collection_audit.csv",
        tables_dir / f"{RUN_SLUG}_output_collection_audit.md",
    )

    summary = json.loads((output_dir / "summary.json").read_text())
    preflight = json.loads((output_dir / "preflight.json").read_text())
    certification = pd.read_csv(output_dir / "operator_certification_summary.csv")
    mavenn = pd.read_csv(output_dir / "mavenn_metrics.csv")
    transfer = pd.read_csv(output_dir / "roundtrip_transfer.csv")
    graphlets = pd.read_csv(output_dir / "graphlet_catalog.csv")
    injected = pd.read_csv(output_dir / "injected_terms.csv")

    profile_failures = (
        mavenn.groupby(["phase", "operator_profile", "gpmap_type"], dropna=False)
        .apply(classify_failure, include_groups=False)
        .reset_index(name="failure_mode")
    )
    profile_failure_counts = (
        mavenn.assign(
            fit_failed_bool=bool_series(mavenn["fit_failed"]),
            preflight_skip_bool=bool_series(mavenn["mavenn_fit_skipped_preflight"]),
        )
        .groupby(["phase", "operator_profile", "gpmap_type"], dropna=False)
        .agg(
            n_fit_attempt_rows=("spectrum_profile", "size"),
            n_failed_fit_attempt_rows=("fit_failed_bool", "sum"),
            n_preflight_skips=("preflight_skip_bool", "sum"),
            spectra=("spectrum_profile", lambda values: ",".join(sorted(map(str, values)))),
            error_examples=("mavenn_fit_error", lambda values: " | ".join(sorted({short_error(x) for x in values if not pd.isna(x)}))),
        )
        .reset_index()
    )
    operator_status = (
        certification.merge(profile_failures, on=["phase", "operator_profile", "gpmap_type"], how="left")
        .merge(profile_failure_counts, on=["phase", "operator_profile", "gpmap_type"], how="left")
    )
    operator_status = operator_status[
        [
            "phase",
            "operator_profile",
            "gpmap_type",
            "operator_certified",
            "n_required_metric_profile_checks",
            "n_ratio_checks_reported",
            "n_ratio_checks_passed",
            "n_fit_failures",
            "n_preflight_skips",
            "failure_mode",
            "error_examples",
        ]
    ]
    write_table(
        operator_status,
        tables_dir / f"{RUN_SLUG}_operator_certification_status.csv",
        tables_dir / f"{RUN_SLUG}_operator_certification_status.md",
    )

    fit_failure_modes = mavenn.copy()
    fit_failure_modes["error_short"] = fit_failure_modes["mavenn_fit_error"].map(short_error)
    fit_failure_modes = fit_failure_modes[
        [
            "phase",
            "operator_profile",
            "spectrum_profile",
            "gpmap_type",
            "fit_failed",
            "mavenn_fit_skipped_preflight",
            "mavenn_estimated_pairwise_parameter_bytes",
            "error_short",
        ]
    ]
    write_table(
        fit_failure_modes,
        tables_dir / f"{RUN_SLUG}_fit_failure_modes.csv",
        tables_dir / f"{RUN_SLUG}_fit_failure_modes.md",
    )

    key_transfer = transfer[bool_series(transfer["certification_metric"])].copy()
    key_transfer = key_transfer[
        [
            "phase",
            "operator_profile",
            "spectrum_profile",
            "metric",
            "injected_mean",
            "recovered_mean",
            "recovered_to_injected_ratio",
            "ratio_within_tolerance",
            "fit_failed",
        ]
    ].sort_values(["operator_profile", "spectrum_profile", "metric"])
    write_table(
        key_transfer,
        tables_dir / f"{RUN_SLUG}_certification_metric_matrix.csv",
        tables_dir / f"{RUN_SLUG}_certification_metric_matrix.md",
    )

    metric_pass_summary = (
        key_transfer.assign(
            ratio_within_tolerance_bool=bool_series(key_transfer["ratio_within_tolerance"]),
            fit_failed_bool=bool_series(key_transfer["fit_failed"]),
        )
        .groupby(["phase", "operator_profile", "metric"], dropna=False)
        .agg(
            n_spectra=("spectrum_profile", "nunique"),
            n_pass=("ratio_within_tolerance_bool", "sum"),
            n_failed_fits=("fit_failed_bool", "sum"),
            n_reported_ratios=("recovered_to_injected_ratio", lambda values: values.notna().sum()),
        )
        .reset_index()
    )
    write_table(
        metric_pass_summary,
        tables_dir / f"{RUN_SLUG}_metric_pass_summary.csv",
        tables_dir / f"{RUN_SLUG}_metric_pass_summary.md",
    )

    graphlet_summary = (
        graphlets.groupby(["graphlet_dimension", "graphlet_split"], dropna=False)
        .agg(n_graphlets=("graphlet_id", "nunique"))
        .reset_index()
        .sort_values(["graphlet_dimension", "graphlet_split"])
    )
    write_table(
        graphlet_summary,
        tables_dir / f"{RUN_SLUG}_graphlet_catalog_summary.csv",
        tables_dir / f"{RUN_SLUG}_graphlet_catalog_summary.md",
    )

    injected_summary = (
        injected.groupby(["spectrum_profile", "order"], dropna=False)
        .agg(
            n_terms=("term_id", "size"),
            mean_abs_coefficient=("coefficient", lambda values: values.abs().mean()),
            max_abs_coefficient=("coefficient", lambda values: values.abs().max()),
        )
        .reset_index()
        .sort_values(["spectrum_profile", "order"])
    )
    write_table(
        injected_summary,
        tables_dir / f"{RUN_SLUG}_injected_spectrum_summary.csv",
        tables_dir / f"{RUN_SLUG}_injected_spectrum_summary.md",
    )

    logs_exit_code = None
    exit_path = log_dir / "labproj-exit-code.txt"
    if exit_path.exists():
        logs_exit_code = exit_path.read_text().strip()

    blackbox_api_failures = int(
        mavenn["mavenn_fit_error"].fillna("").str.contains("MultilayerPerceptronGPMap", regex=False).sum()
    )
    metrics = {
        "analysis_id": ANALYSIS_ID,
        "run_id": RUN_ID,
        "experiment_id": "EXP-050",
        "source_run_id": preflight.get("source_run_id"),
        "mode": preflight.get("mode"),
        "assay_id": preflight.get("assay_id"),
        "n_assay_rows": int(preflight.get("n_assay_rows")),
        "n_graphlets": int(preflight.get("n_graphlets")),
        "n_operator_profiles": int(len(summary.get("operator_profiles", []))),
        "n_known_spectrum_profiles": int(len(summary.get("known_spectrum_profiles", []))),
        "n_injected_terms": int(summary.get("n_injected_terms")),
        "n_fit_attempt_rows": int(len(mavenn)),
        "n_fit_failures": int(bool_series(mavenn["fit_failed"]).sum()),
        "n_preflight_skips": int(bool_series(mavenn["mavenn_fit_skipped_preflight"]).sum()),
        "n_blackbox_api_failures": blackbox_api_failures,
        "n_operator_profiles_certified": int(bool_series(certification["operator_certified"]).sum()),
        "certified_operator_profiles": sorted(
            certification.loc[bool_series(certification["operator_certified"]), "operator_profile"].tolist()
        ),
        "n_required_metric_profile_checks_total": int(certification["n_required_metric_profile_checks"].sum()),
        "n_ratio_checks_reported_total": int(certification["n_ratio_checks_reported"].sum()),
        "n_ratio_checks_passed_total": int(certification["n_ratio_checks_passed"].sum()),
        "n_certification_metric_transfer_rows": int(len(key_transfer)),
        "n_certification_metric_rows_with_recovered_values": int(key_transfer["recovered_mean"].notna().sum()),
        "log_exit_code": logs_exit_code,
        "output_bytes_collected": int(audit.loc[audit["role"] == "run_output", "bytes"].sum()),
        "log_bytes_collected": int(audit.loc[audit["role"] == "run_log", "bytes"].sum()),
        "source_output_dir": str(output_dir.relative_to(project_root)),
        "source_log_dir": str(log_dir.relative_to(project_root)),
        "remote_runtime_output_dir": "/tank/data/labproj/runtime/project-adapt-env/runs/RUN-129/outputs",
    }
    (tables_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
