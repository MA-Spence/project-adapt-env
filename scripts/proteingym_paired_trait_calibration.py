#!/usr/bin/env python3
"""Run paired ProteinGym state/function calibration diagnostics for HYP-007."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from project_adapt_env._compat import ensure_external_paths  # noqa: E402
from project_adapt_env.proteingym_panel import (  # noqa: E402
    PreparedAssay,
    calibration_options_from_config,
    load_config,
    merge_dict,
    prepare_proteingym_panel,
)
from project_adapt_env.utils import (  # noqa: E402
    ProgressTracker,
    atomic_write_dataframe_csv,
    atomic_write_json,
    build_logger,
    flatten_mapping,
    package_version,
    to_builtin,
)

ensure_external_paths()

from adaptenv import LandscapeConfig  # noqa: E402
from adaptenv.calibration import (  # noqa: E402
    LandscapeCalibrator,
    _replace_config,
    summarize_empirical_landscapes,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def make_stage_reporter(progress: ProgressTracker):
    stage_started_at: dict[str, float] = {}

    def report(
        *,
        stage: str,
        completed: int,
        total: int,
        message: str,
        details: dict[str, Any] | None = None,
    ) -> None:
        if stage not in stage_started_at:
            stage_started_at[stage] = time.time()
        progress.write(
            stage=stage,
            completed=completed,
            total=total,
            message=message,
            details=details,
            stage_started_at=stage_started_at[stage],
        )

    return report


def _role_map(config_roles: dict[str, Any]) -> dict[str, str]:
    role_by_assay: dict[str, str] = {}
    for role_name, assay_ids in config_roles.items():
        for assay_id in assay_ids:
            role_by_assay[str(assay_id)] = str(role_name)
    return role_by_assay


def _assay_subset(
    *,
    assay_ids: list[str],
    assays_by_id: dict[str, PreparedAssay],
    source: str,
) -> tuple[list[Any], list[Path], list[str]]:
    landscapes: list[Any] = []
    alignments: list[Path] = []
    wildtypes: list[str] = []
    for assay_id in assay_ids:
        assay = assays_by_id.get(assay_id)
        if assay is None:
            raise ValueError(f"Configured assay id not present in prepared panel: {assay_id}")
        if source == "raw":
            landscape = assay.raw_landscape
        elif source == "latent":
            landscape = assay.latent_landscape
        else:
            raise ValueError(f"Unsupported source: {source}")
        landscapes.append(landscape)
        alignments.append(assay.alignment_path)
        wildtypes.append(assay.wildtype_sequence)
    return landscapes, alignments, wildtypes


def _option_config(
    defaults: dict[str, Any],
    overrides: dict[str, Any] | None,
) -> dict[str, Any]:
    if overrides:
        return merge_dict(defaults, overrides)
    return dict(defaults)


def _prefixed(prefix: str, value: dict[str, Any]) -> dict[str, Any]:
    return flatten_mapping(prefix, to_builtin(value))


def _paired_readout_rows(
    *,
    assays_by_id: dict[str, PreparedAssay],
    state_assays: list[str],
    function_assays: list[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for state_id in state_assays:
        state = assays_by_id[state_id]
        state_frame = pd.read_csv(state.assay_csv_path)
        state_frame = state_frame[["mutated_sequence", "mutant", "DMS_score"]].rename(
            columns={
                "mutant": "state_mutant",
                "DMS_score": "state_score",
            }
        )
        for function_id in function_assays:
            function = assays_by_id[function_id]
            function_frame = pd.read_csv(function.assay_csv_path)
            function_frame = function_frame[
                ["mutated_sequence", "mutant", "DMS_score"]
            ].rename(
                columns={
                    "mutant": "function_mutant",
                    "DMS_score": "function_score",
                }
            )
            merged = state_frame.merge(function_frame, on="mutated_sequence", how="inner")
            spearman = None
            if len(merged) >= 2:
                spearman = float(
                    spearmanr(merged["state_score"], merged["function_score"]).statistic
                )
            rows.append(
                {
                    "state_assay": state_id,
                    "function_assay": function_id,
                    "matched_variants": int(len(merged)),
                    "state_function_spearman": spearman,
                    "state_score_min": float(merged["state_score"].min()) if len(merged) else None,
                    "state_score_max": float(merged["state_score"].max()) if len(merged) else None,
                    "function_score_min": float(merged["function_score"].min()) if len(merged) else None,
                    "function_score_max": float(merged["function_score"].max()) if len(merged) else None,
                }
            )
    return rows


def run_branch(
    *,
    branch_name: str,
    branch_cfg: dict[str, Any],
    config: dict[str, Any],
    assays_by_id: dict[str, PreparedAssay],
) -> tuple[dict[str, Any], dict[str, Any]]:
    source = str(branch_cfg.get("source") or "raw")
    trait_roles = dict(config["trait_roles"])
    state_assays = [
        str(item)
        for item in branch_cfg.get("state_assays", trait_roles.get("state", []))
    ]
    function_assays = [
        str(item)
        for item in branch_cfg.get("function_assays", trait_roles.get("function", []))
    ]
    use_state = bool(branch_cfg.get("use_state", True))
    use_function = bool(branch_cfg.get("use_function", True))

    base_config_dict = merge_dict(
        dict(config["adapt_env"]["base_config"]),
        dict(branch_cfg.get("base_config", {})),
    )
    base_config = LandscapeConfig(**base_config_dict)
    state_options = calibration_options_from_config(
        _option_config(
            dict(config["calibration_defaults"]["state"]),
            branch_cfg.get("state_calibration"),
        )
    )
    function_options = calibration_options_from_config(
        _option_config(
            dict(config["calibration_defaults"]["function"]),
            branch_cfg.get("function_calibration"),
        )
    )

    fitted_parameters: dict[str, Any] = {}
    objective_terms: dict[str, Any] = {}
    state_summary = None
    function_summary = None

    if use_state:
        state_landscapes, state_alignments, state_wildtypes = _assay_subset(
            assay_ids=state_assays,
            assays_by_id=assays_by_id,
            source=source,
        )
        state_summary = summarize_empirical_landscapes(
            state_landscapes,
            kind="unfolding",
            alignment_profiles=state_alignments,
            wildtypes=state_wildtypes,
            options=state_options,
        )
        state_calibrator = LandscapeCalibrator(
            base_config=base_config,
            options=state_options,
        )
        state_updates, state_objectives = state_calibrator.calibrate_unfolding(
            state_summary
        )
        fitted_parameters.update(state_updates)
        objective_terms.update({f"state_{key}": value for key, value in state_objectives.items()})
    else:
        state_updates = {}

    if use_function:
        function_landscapes, function_alignments, function_wildtypes = _assay_subset(
            assay_ids=function_assays,
            assays_by_id=assays_by_id,
            source=source,
        )
        function_summary = summarize_empirical_landscapes(
            function_landscapes,
            kind="functional",
            alignment_profiles=function_alignments,
            wildtypes=function_wildtypes,
            options=function_options,
        )
        function_calibrator = LandscapeCalibrator(
            base_config=base_config,
            options=function_options,
        )
        function_updates, function_objectives = function_calibrator.calibrate_functional(
            function_summary,
            stability_updates=state_updates,
        )
        fitted_parameters.update(function_updates)
        objective_terms.update(
            {f"function_{key}": value for key, value in function_objectives.items()}
        )

    final_config = _replace_config(base_config, **fitted_parameters)
    state_validation = (
        LandscapeCalibrator(base_config=base_config, options=state_options)._validate_fit(
            final_config,
            unfolding_summary=state_summary,
            functional_summary=None,
        )
        if state_summary is not None
        else {}
    )
    function_validation = (
        LandscapeCalibrator(base_config=base_config, options=function_options)._validate_fit(
            final_config,
            unfolding_summary=None,
            functional_summary=function_summary,
        )
        if function_summary is not None
        else {}
    )
    validation = {
        **{f"state_{key}": value for key, value in state_validation.items()},
        **{f"function_{key}": value for key, value in function_validation.items()},
    }

    row = {
        "branch": branch_name,
        "source": source,
        "use_state": use_state,
        "use_function": use_function,
        "state_assays": ",".join(state_assays if use_state else []),
        "function_assays": ",".join(function_assays if use_function else []),
        **_prefixed("validation", validation),
        **_prefixed("fitted", fitted_parameters),
        **_prefixed("objective", objective_terms),
    }
    summary = {
        "source": source,
        "use_state": use_state,
        "use_function": use_function,
        "state_assays": state_assays if use_state else [],
        "function_assays": function_assays if use_function else [],
        "state_synthetic_readout_mode": state_options.synthetic_readout_mode,
        "function_synthetic_readout_mode": function_options.synthetic_readout_mode,
        "fitted_parameters": to_builtin(fitted_parameters),
        "validation": to_builtin(validation),
        "objective_terms": to_builtin(objective_terms),
    }
    return row, summary


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(Path(args.config))
    if args.quick and "quick" in config:
        config = merge_dict(config, config["quick"])

    experiment_id = str(config.get("experiment_id") or config["experiment"]["id"])
    logger = build_logger(
        "project_adapt_env.proteingym_paired_trait",
        output_dir / "progress.log",
    )
    progress = ProgressTracker(
        progress_path=output_dir / "progress.json",
        logger=logger,
        run_label=experiment_id,
    )
    report = make_stage_reporter(progress)

    panel = prepare_proteingym_panel(
        project_root=project_root,
        proteingym_cfg=dict(config["proteingym"]),
        panel_cfg=dict(config["panel"]),
        mmseqs_cfg=dict(config["mmseqs"]),
        mavenn_cfg=dict(config["mavenn"]),
        calibration_max_mutation_count=int(config["calibration_max_mutation_count"]),
        checkpoint_dir=checkpoint_dir,
        progress_callback=lambda event: report(
            stage="panel_preparation",
            completed=int(event["completed"]),
            total=int(event["total"]),
            message=f"Prepared {event['dms_id']}",
            details=event,
        )
        if event.get("event") == "prepared_assay"
        else None,
    )
    assays_by_id = {assay.dms_id: assay for assay in panel.assays}
    role_by_assay = _role_map(dict(config["trait_roles"]))
    missing_roles = sorted(set(role_by_assay) - set(assays_by_id))
    if missing_roles:
        raise ValueError(
            "Trait roles reference ids that were not prepared: "
            + ", ".join(missing_roles)
        )

    branch_rows: list[dict[str, Any]] = []
    branch_summary: dict[str, Any] = {}
    branches_cfg = dict(config["branches"])
    for branch_index, (branch_name, branch_cfg) in enumerate(branches_cfg.items(), start=1):
        row, summary = run_branch(
            branch_name=branch_name,
            branch_cfg=dict(branch_cfg),
            config=config,
            assays_by_id=assays_by_id,
        )
        branch_rows.append(row)
        branch_summary[branch_name] = summary
        report(
            stage="branch_fits",
            completed=branch_index,
            total=len(branches_cfg),
            message=f"Calibrated {branch_name}",
            details={
                "branch": branch_name,
                "use_state": summary["use_state"],
                "use_function": summary["use_function"],
            },
        )

    panel_frame = panel.panel_df.copy()
    panel_frame["trait_role"] = panel_frame["DMS_id"].map(role_by_assay).fillna("unassigned")
    panel_csv_path = output_dir / "selected_panel.csv"
    atomic_write_dataframe_csv(panel_csv_path, panel_frame, index=False)

    mavenn_frame = panel.mavenn_metrics_frame()
    mavenn_frame["trait_role"] = mavenn_frame["dms_id"].map(role_by_assay).fillna("unassigned")
    mavenn_metrics_path = output_dir / "mavenn_assay_metrics.csv"
    atomic_write_dataframe_csv(mavenn_metrics_path, mavenn_frame, index=False)

    branch_validation_path = output_dir / "branch_validations.csv"
    atomic_write_dataframe_csv(branch_validation_path, pd.DataFrame(branch_rows), index=False)

    readout_rows = _paired_readout_rows(
        assays_by_id=assays_by_id,
        state_assays=[str(item) for item in config["trait_roles"].get("state", [])],
        function_assays=[str(item) for item in config["trait_roles"].get("function", [])],
    )
    readout_correlation_path = output_dir / "paired_readout_correlations.csv"
    atomic_write_dataframe_csv(
        readout_correlation_path,
        pd.DataFrame(readout_rows),
        index=False,
    )

    summary_payload = {
        "experiment_id": config["experiment"]["id"],
        "proteingym_version": config["proteingym"]["version"],
        "proteingym_package_version": package_version("proteingym"),
        "mavenn_package_version": package_version("mavenn"),
        "panel_csv_path": str(panel_csv_path),
        "mavenn_assay_metrics_csv_path": str(mavenn_metrics_path),
        "branch_validation_csv_path": str(branch_validation_path),
        "paired_readout_correlation_csv_path": str(readout_correlation_path),
        "trait_roles": role_by_assay,
        "panel": panel.summary_payload(),
        "paired_readout_correlations": readout_rows,
        "branches": branch_summary,
    }
    summary_path = output_dir / "summary.json"
    atomic_write_json(summary_path, summary_payload)

    report(
        stage="completed",
        completed=1,
        total=1,
        message="Paired trait calibration finished",
        details={"summary_json": str(summary_path)},
    )


if __name__ == "__main__":
    main()
