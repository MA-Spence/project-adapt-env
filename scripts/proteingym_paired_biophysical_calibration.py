#!/usr/bin/env python3
"""Run paired ProteinGym abundance-binding calibration for HYP-007."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd


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
from adaptenv.calibration import calibrate_synthetic_landscape  # noqa: E402


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


def branch_row(
    *,
    branch: str,
    source: str,
    use_abundance: bool,
    use_binding: bool,
    abundance_assays: list[str],
    binding_assays: list[str],
    fitted_parameters: dict[str, Any],
    validation: dict[str, Any],
    objective_terms: dict[str, Any],
) -> dict[str, Any]:
    return {
        "branch": branch,
        "source": source,
        "use_abundance": use_abundance,
        "use_binding": use_binding,
        "abundance_assays": ",".join(abundance_assays),
        "binding_assays": ",".join(binding_assays),
        **flatten_mapping("validation", to_builtin(validation)),
        **flatten_mapping("fitted", to_builtin(fitted_parameters)),
        **flatten_mapping("objective", to_builtin(objective_terms)),
    }


def _assay_role_map(config_roles: dict[str, Any]) -> dict[str, str]:
    role_map: dict[str, str] = {}
    for role_name, assay_ids in config_roles.items():
        for assay_id in assay_ids:
            role_map[str(assay_id)] = str(role_name)
    return role_map


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
            raise ValueError(f"Unsupported branch source: {source}")
        landscapes.append(landscape)
        alignments.append(assay.alignment_path)
        wildtypes.append(assay.wildtype_sequence)
    return landscapes, alignments, wildtypes


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

    experiment_id = str(config.get("experiment_id") or "EXP-009")
    logger = build_logger("project_adapt_env.proteingym_paired_biophysical", output_dir / "progress.log")
    progress = ProgressTracker(
        progress_path=output_dir / "progress.json",
        logger=logger,
        run_label=experiment_id,
    )
    report = make_stage_reporter(progress)

    proteingym_cfg = dict(config["proteingym"])
    panel_cfg = dict(config["panel"])
    mmseqs_cfg = dict(config["mmseqs"])
    mavenn_cfg = dict(config["mavenn"])
    base_config_dict = dict(config["adapt_env"]["base_config"])
    branches_cfg = dict(config["branches"])
    assay_roles_cfg = dict(config["assay_roles"])
    role_map = _assay_role_map(assay_roles_cfg)

    panel = prepare_proteingym_panel(
        project_root=project_root,
        proteingym_cfg=proteingym_cfg,
        panel_cfg=panel_cfg,
        mmseqs_cfg=mmseqs_cfg,
        mavenn_cfg=mavenn_cfg,
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

    missing_roles = sorted(set(role_map) - set(assays_by_id))
    if missing_roles:
        raise ValueError(
            "Assay roles reference ids that were not prepared: "
            + ", ".join(missing_roles)
        )

    branch_rows: list[dict[str, Any]] = []
    branch_summary: dict[str, Any] = {}
    branch_names = list(branches_cfg)

    for branch_index, (branch_name, branch_cfg) in enumerate(branches_cfg.items(), start=1):
        source = str(branch_cfg.get("source") or "raw")
        use_abundance = bool(branch_cfg.get("use_abundance", True))
        use_binding = bool(branch_cfg.get("use_binding", True))
        abundance_assays = [
            str(item)
            for item in branch_cfg.get("abundance_assays", assay_roles_cfg.get("abundance", []))
        ]
        binding_assays = [
            str(item)
            for item in branch_cfg.get("binding_assays", assay_roles_cfg.get("binding", []))
        ]
        options = calibration_options_from_config(dict(branch_cfg["calibration"]))

        unfolding_landscapes = None
        unfolding_alignments = None
        unfolding_wildtypes = None
        if use_abundance:
            (
                unfolding_landscapes,
                unfolding_alignments,
                unfolding_wildtypes,
            ) = _assay_subset(
                assay_ids=abundance_assays,
                assays_by_id=assays_by_id,
                source=source,
            )

        functional_landscapes = None
        functional_alignments = None
        functional_wildtypes = None
        if use_binding:
            (
                functional_landscapes,
                functional_alignments,
                functional_wildtypes,
            ) = _assay_subset(
                assay_ids=binding_assays,
                assays_by_id=assays_by_id,
                source=source,
            )

        result = calibrate_synthetic_landscape(
            unfolding_landscapes=unfolding_landscapes,
            functional_landscapes=functional_landscapes,
            unfolding_alignment_profiles=unfolding_alignments,
            functional_alignment_profiles=functional_alignments,
            unfolding_wildtypes=unfolding_wildtypes,
            functional_wildtypes=functional_wildtypes,
            base_config=LandscapeConfig(**base_config_dict),
            options=options,
        )

        row = branch_row(
            branch=branch_name,
            source=source,
            use_abundance=use_abundance,
            use_binding=use_binding,
            abundance_assays=abundance_assays if use_abundance else [],
            binding_assays=binding_assays if use_binding else [],
            fitted_parameters=result.fitted_parameters,
            validation=result.validation,
            objective_terms=result.objective_terms,
        )
        branch_rows.append(row)
        branch_summary[branch_name] = {
            "source": source,
            "use_abundance": use_abundance,
            "use_binding": use_binding,
            "abundance_assays": abundance_assays if use_abundance else [],
            "binding_assays": binding_assays if use_binding else [],
            "fitted_parameters": to_builtin(result.fitted_parameters),
            "validation": to_builtin(result.validation),
            "objective_terms": to_builtin(result.objective_terms),
        }
        report(
            stage="branch_fits",
            completed=branch_index,
            total=len(branch_names),
            message=f"Calibrated {branch_name}",
            details={
                "branch": branch_name,
                "source": source,
                "use_abundance": use_abundance,
                "use_binding": use_binding,
            },
        )

    panel_frame = panel.panel_df.copy()
    panel_frame["assay_role"] = panel_frame["DMS_id"].map(role_map).fillna("unassigned")
    panel_csv_path = output_dir / "selected_panel.csv"
    atomic_write_dataframe_csv(panel_csv_path, panel_frame, index=False)

    mavenn_frame = panel.mavenn_metrics_frame()
    mavenn_frame["assay_role"] = mavenn_frame["dms_id"].map(role_map).fillna("unassigned")
    mavenn_metrics_path = output_dir / "mavenn_assay_metrics.csv"
    atomic_write_dataframe_csv(mavenn_metrics_path, mavenn_frame, index=False)

    branch_validation_path = output_dir / "branch_validations.csv"
    atomic_write_dataframe_csv(branch_validation_path, pd.DataFrame(branch_rows), index=False)

    summary_payload = {
        "experiment_id": config["experiment"]["id"],
        "proteingym_version": proteingym_cfg["version"],
        "proteingym_package_version": package_version("proteingym"),
        "mavenn_package_version": package_version("mavenn"),
        "panel_csv_path": str(panel_csv_path),
        "mavenn_assay_metrics_csv_path": str(mavenn_metrics_path),
        "branch_validation_csv_path": str(branch_validation_path),
        "assay_roles": role_map,
        "panel": panel.summary_payload(),
        "branches": branch_summary,
    }
    summary_path = output_dir / "summary.json"
    atomic_write_json(summary_path, summary_payload)

    report(
        stage="completed",
        completed=1,
        total=1,
        message="Paired abundance-binding calibration finished",
        details={"summary_json": str(summary_path)},
    )


if __name__ == "__main__":
    main()
