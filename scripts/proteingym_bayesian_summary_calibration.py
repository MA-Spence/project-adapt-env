#!/usr/bin/env python3
"""Run shared SMC-ABC summary calibration for HYP-001."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from project_adapt_env._compat import ensure_external_paths  # noqa: E402
from project_adapt_env.adapt_env_bayes import (  # noqa: E402
    AdaptEnvSMCABCProblem,
    build_empirical_target,
    build_parameter_specs,
    posterior_mean_updates,
    restore_summary_target,
    run_synthetic_truth_recovery,
    serialize_summary_target,
    validate_updates,
)
from project_adapt_env.proteingym_panel import (  # noqa: E402
    calibration_options_from_config,
    installed_mavenn_version,
    installed_proteingym_version,
    load_config,
    merge_dict,
    prepare_proteingym_panel,
)
from project_adapt_env.smc_abc import (  # noqa: E402
    SMCABCBackendConfig,
    SMCABCCheckpointConfig,
    SMCABCConfig,
    run_smc_abc,
)
from project_adapt_env.utils import (  # noqa: E402
    ProgressTracker,
    atomic_write_dataframe_csv,
    atomic_write_json,
    build_logger,
    detect_worker_count,
    flatten_mapping,
    to_builtin,
)

ensure_external_paths()

from adaptenv import LandscapeConfig  # noqa: E402
from adaptenv.calibration import (  # noqa: E402
    calibrate_synthetic_landscape,
    summarize_empirical_landscapes,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def branch_row(
    *,
    branch: str,
    fit_mode: str,
    source: str,
    fitted_parameters: dict[str, Any],
    validation: dict[str, Any],
    objective_terms: dict[str, Any],
) -> dict[str, Any]:
    return {
        "branch": branch,
        "fit_mode": fit_mode,
        "source": source,
        **flatten_mapping("validation", to_builtin(validation)),
        **flatten_mapping("fitted", to_builtin(fitted_parameters)),
        **flatten_mapping("objective", to_builtin(objective_terms)),
    }


def weighted_parameter_summary(particles: list[Any], parameter_names: list[str]) -> list[dict[str, Any]]:
    weights = np.asarray([particle.weight for particle in particles], dtype=np.float64)
    total = float(np.sum(weights))
    if total <= 0 or not np.isfinite(total):
        weights = np.full(len(particles), 1.0 / max(len(particles), 1))
    else:
        weights = weights / total
    rows = []
    for name in parameter_names:
        values = np.asarray([particle.parameters[name] for particle in particles], dtype=np.float64)
        mean = float(np.sum(weights * values))
        rows.append(
            {
                "parameter": name,
                "weighted_mean": mean,
                "weighted_std": float(np.sqrt(np.sum(weights * (values - mean) ** 2))),
                "q05": float(np.quantile(values, 0.05)),
                "q50": float(np.quantile(values, 0.50)),
                "q95": float(np.quantile(values, 0.95)),
            }
        )
    return rows


def build_backend_config(backend_cfg: dict[str, Any] | None) -> SMCABCBackendConfig:
    backend_cfg = dict(backend_cfg or {})
    max_workers = int(backend_cfg.get("max_workers") or detect_worker_count(default=1))
    kind = str(backend_cfg.get("kind") or ("ray" if max_workers > 1 else "serial"))
    batch_size = int(backend_cfg.get("batch_size") or max(1, min(max_workers, 8)))
    return SMCABCBackendConfig(
        kind=kind,
        max_workers=max(max_workers, 1),
        batch_size=max(batch_size, 1),
    )


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


def smc_event_message(event: dict[str, Any]) -> str:
    best = event.get("best_distance")
    median = event.get("median_distance")
    best_text = f"{float(best):.4f}" if best is not None else "n/a"
    median_text = f"{float(median):.4f}" if median is not None else "n/a"
    if event.get("event") == "round_complete":
        epsilon = event.get("epsilon")
        epsilon_text = f"{float(epsilon):.4f}" if epsilon is not None else "n/a"
        return f"Round complete; epsilon={epsilon_text}, best={best_text}, median={median_text}"
    return f"Evaluating particles; best={best_text}, median={median_text}"


def restore_resume_cache(
    *,
    output_dir: Path,
    checkpoint_dir: Path,
    resume_cfg: dict[str, Any] | None,
    logger: Any,
) -> None:
    resume_cfg = dict(resume_cfg or {})
    source_run_id = str(resume_cfg.get("source_run_id") or "").strip()
    if not source_run_id:
        return

    runs_dir = output_dir.parent.parent
    source_output_dir = runs_dir / source_run_id / "outputs"
    source_checkpoint_dir = source_output_dir / "checkpoints"
    if not source_checkpoint_dir.is_dir():
        logger.warning(
            "Resume source run %s did not provide checkpoints at %s.",
            source_run_id,
            source_checkpoint_dir,
        )
        return

    copied_files = 0
    copied_bytes = 0
    for source_path in sorted(source_checkpoint_dir.rglob("*")):
        if not source_path.is_file():
            continue
        relative = source_path.relative_to(source_checkpoint_dir)
        target_path = checkpoint_dir / relative
        if target_path.exists():
            continue
        target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, target_path)
        copied_files += 1
        copied_bytes += source_path.stat().st_size

    synthetic_truth_csv = source_output_dir / "synthetic_truth_recovery.csv"
    if synthetic_truth_csv.is_file():
        target_truth_csv = output_dir / "synthetic_truth_recovery.csv"
        if not target_truth_csv.exists():
            shutil.copy2(synthetic_truth_csv, target_truth_csv)
            copied_files += 1
            copied_bytes += synthetic_truth_csv.stat().st_size

    logger.info(
        "Restored resume cache from %s into %s (%d files, %.2f MiB).",
        source_output_dir,
        output_dir,
        copied_files,
        copied_bytes / (1024.0 * 1024.0),
    )


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
    experiment_id = str(config.get("experiment_id") or "EXP-005")

    logger = build_logger("project_adapt_env.proteingym_smc_abc", output_dir / "progress.log")
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
    deterministic_branches = dict(config["deterministic_branches"])
    abc_cfg = dict(config["abc"])
    recovery_cfg = dict(config["synthetic_truth_recovery"])
    resume_cfg = dict(config.get("resume") or {})

    abc_backend = build_backend_config(abc_cfg.get("backend"))
    recovery_backend = build_backend_config(recovery_cfg.get("backend", abc_cfg.get("backend", {})))
    logger.info(
        "Starting %s with backend kind=%s workers=%d batch_size=%d.",
        experiment_id,
        abc_backend.kind,
        abc_backend.max_workers,
        abc_backend.batch_size,
    )
    restore_resume_cache(
        output_dir=output_dir,
        checkpoint_dir=checkpoint_dir,
        resume_cfg=resume_cfg,
        logger=logger,
    )

    mmseqs_cfg["cache_dir"] = str(project_root / mmseqs_cfg["cache_dir"])
    report(
        stage="panel_preparation",
        completed=0,
        total=max(int(panel_cfg["panel_size"]), 1),
        message="Starting ProteinGym panel preparation.",
    )
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
            message=f"Prepared assay {event['dms_id']}.",
            details=event,
        ),
    )

    panel_csv_path = output_dir / "selected_panel.csv"
    atomic_write_dataframe_csv(panel_csv_path, panel.panel_df, index=False)

    mavenn_metrics_path = output_dir / "mavenn_assay_metrics.csv"
    atomic_write_dataframe_csv(mavenn_metrics_path, panel.mavenn_metrics_frame(), index=False)

    shared_branch_rows: list[dict[str, Any]] = []
    branch_summary: dict[str, Any] = {}
    base_config = LandscapeConfig(**base_config_dict)

    total_branches = max(len(deterministic_branches), 1)
    for branch_index, (branch_name, branch_cfg) in enumerate(deterministic_branches.items(), start=1):
        source = str(branch_cfg["source"])
        if source != "raw":
            raise ValueError(f"Unsupported deterministic branch source: {source}")
        options = calibration_options_from_config(branch_cfg["calibration"])
        logger.info("Running deterministic control branch %s.", branch_name)
        result = calibrate_synthetic_landscape(
            functional_landscapes=panel.raw_landscapes,
            functional_alignment_profiles=panel.alignment_paths,
            functional_wildtypes=panel.wildtypes,
            base_config=base_config,
            options=options,
        )
        shared_branch_rows.append(
            branch_row(
                branch=branch_name,
                fit_mode="shared",
                source=source,
                fitted_parameters=result.fitted_parameters,
                validation=result.validation,
                objective_terms=result.objective_terms,
            )
        )
        branch_summary[branch_name] = {
            "fit_mode": "shared",
            "source": source,
            "fitted_parameters": to_builtin(result.fitted_parameters),
            "validation": to_builtin(result.validation),
            "objective_terms": to_builtin(result.objective_terms),
        }
        report(
            stage="deterministic_controls",
            completed=branch_index,
            total=total_branches,
            message=f"Completed deterministic branch {branch_name}.",
            details={"branch": branch_name},
        )

    abc_options = calibration_options_from_config(abc_cfg["calibration"])
    selected_features = tuple(str(item) for item in abc_cfg.get("selected_features", [])) or None
    if selected_features is None:
        raise ValueError("abc.selected_features must be configured explicitly")

    target_cache_path = checkpoint_dir / "empirical_target.json"
    if target_cache_path.is_file():
        empirical_collection = summarize_empirical_landscapes(
            panel.raw_landscapes,
            kind="functional",
            alignment_profiles=panel.alignment_paths,
            wildtypes=panel.wildtypes,
            options=abc_options,
        )
        target_payload = json.loads(target_cache_path.read_text(encoding="utf-8"))
        target = restore_summary_target(
            payload=target_payload,
            empirical_collection=empirical_collection,
        )
        report(
            stage="empirical_target_bootstrap",
            completed=int(abc_cfg["bootstrap_replicates"]),
            total=int(abc_cfg["bootstrap_replicates"]),
            message="Restored cached empirical summary target.",
        )
    else:
        report(
            stage="empirical_target_bootstrap",
            completed=0,
            total=int(abc_cfg["bootstrap_replicates"]),
            message="Building empirical summary target.",
        )
        target = build_empirical_target(
            landscapes=panel.raw_landscapes,
            alignment_paths=panel.alignment_paths,
            wildtypes=panel.wildtypes,
            options=abc_options,
            selected_features=selected_features,
            bootstrap_replicates=int(abc_cfg["bootstrap_replicates"]),
            covariance_shrinkage=float(abc_cfg["covariance_shrinkage"]),
            covariance_ridge=float(abc_cfg["covariance_ridge"]),
            progress_callback=lambda event: report(
                stage="empirical_target_bootstrap",
                completed=int(event["completed"]),
                total=int(event["total"]),
                message="Bootstrapping empirical summary target.",
                details=event,
            ),
        )
        atomic_write_json(target_cache_path, serialize_summary_target(target))

    target_features_path = output_dir / "target_features.csv"
    bootstrap_sd = np.std(target.bootstrap_vectors, axis=0, ddof=1)
    atomic_write_dataframe_csv(
        target_features_path,
        pd.DataFrame(
            {
                "feature_name": target.feature_names,
                "observed_value": target.observed_vector,
                "bootstrap_sd": bootstrap_sd,
            }
        ),
        index=False,
    )

    specs = build_parameter_specs(list(abc_cfg["parameters"]))
    problem = AdaptEnvSMCABCProblem(
        target=target,
        base_config=base_config,
        options=abc_options,
        specs=specs,
        selected_features=selected_features,
        replicates_per_particle=int(abc_cfg["replicates_per_particle"]),
    )
    smc_result = run_smc_abc(
        specs=specs,
        config=SMCABCConfig(**abc_cfg["smc"]),
        simulate=problem.simulate,
        backend=abc_backend,
        checkpoint=SMCABCCheckpointConfig(
            path=str((checkpoint_dir / "smc_empirical.json").resolve()),
            resume=True,
        ),
        progress_callback=lambda event: report(
            stage=f"smc_empirical_round_{int(event['round_index'])}",
            completed=int(event["completed_attempts"]),
            total=int(event["total_attempts"]),
            message=smc_event_message(event),
            details=event,
        ),
        run_label="smc_empirical",
        logger=logger,
    )

    posterior_particles_path = output_dir / "posterior_particles.csv"
    atomic_write_dataframe_csv(
        posterior_particles_path,
        pd.DataFrame(
            [
                {
                    "distance": particle.distance,
                    "weight": particle.weight,
                    **{name: particle.parameters[name] for name in sorted(particle.parameters)},
                }
                for particle in smc_result.particles
            ]
        ),
        index=False,
    )

    posterior_rounds_path = output_dir / "posterior_rounds.csv"
    atomic_write_dataframe_csv(
        posterior_rounds_path,
        pd.DataFrame([to_builtin(summary.__dict__) for summary in smc_result.round_summaries]),
        index=False,
    )

    posterior_parameter_summary_path = output_dir / "posterior_parameter_summary.csv"
    atomic_write_dataframe_csv(
        posterior_parameter_summary_path,
        pd.DataFrame(
            weighted_parameter_summary(
                smc_result.particles,
                parameter_names=[spec.name for spec in specs],
            )
        ),
        index=False,
    )

    best_updates = dict(smc_result.best_particle.parameters)
    posterior_mean = posterior_mean_updates(specs=specs, particles=smc_result.particles)
    _, best_validation = validate_updates(
        empirical_collection=target.empirical_collection,
        base_config=base_config,
        options=abc_options,
        updates=best_updates,
    )
    _, mean_validation = validate_updates(
        empirical_collection=target.empirical_collection,
        base_config=base_config,
        options=abc_options,
        updates=posterior_mean,
    )

    shared_branch_rows.extend(
        [
            branch_row(
                branch="smc_abc_best_raw",
                fit_mode="shared",
                source="raw",
                fitted_parameters=best_updates,
                validation=best_validation,
                objective_terms={"distance": smc_result.best_particle.distance},
            ),
            branch_row(
                branch="smc_abc_posterior_mean_raw",
                fit_mode="shared",
                source="raw",
                fitted_parameters=posterior_mean,
                validation=mean_validation,
                objective_terms={
                    "distance": float(
                        np.average(
                            np.asarray([particle.distance for particle in smc_result.particles], dtype=np.float64),
                            weights=np.asarray([particle.weight for particle in smc_result.particles], dtype=np.float64),
                        )
                    )
                },
            ),
        ]
    )
    branch_summary["smc_abc_best_raw"] = {
        "fit_mode": "shared",
        "source": "raw",
        "fitted_parameters": to_builtin(best_updates),
        "validation": to_builtin(best_validation),
        "objective_terms": {"distance": float(smc_result.best_particle.distance)},
    }
    branch_summary["smc_abc_posterior_mean_raw"] = {
        "fit_mode": "shared",
        "source": "raw",
        "fitted_parameters": to_builtin(posterior_mean),
        "validation": to_builtin(mean_validation),
        "objective_terms": {
            "distance": float(
                np.average(
                    np.asarray([particle.distance for particle in smc_result.particles], dtype=np.float64),
                    weights=np.asarray([particle.weight for particle in smc_result.particles], dtype=np.float64),
                )
            )
        },
    }

    branch_validations_path = output_dir / "branch_validations.csv"
    atomic_write_dataframe_csv(branch_validations_path, pd.DataFrame(shared_branch_rows), index=False)

    synthetic_truth_path = output_dir / "synthetic_truth_recovery.csv"
    recovery_rows = run_synthetic_truth_recovery(
        target=target,
        base_config=base_config,
        options=abc_options,
        specs=specs,
        smc_config=SMCABCConfig(**recovery_cfg["smc"]),
        backend=recovery_backend,
        truths=list(recovery_cfg["truths"]),
        replicates_per_particle=int(recovery_cfg["replicates_per_particle"]),
        selected_features=selected_features,
        checkpoint_dir=checkpoint_dir / "synthetic_truth",
        partial_csv_path=synthetic_truth_path,
        progress_callback=lambda label, event: report(
            stage=f"synthetic_truth_{label}_round_{int(event['round_index'])}",
            completed=int(event["completed_attempts"]),
            total=int(event["total_attempts"]),
            message=smc_event_message(event),
            details={"label": label, **event},
        ),
        logger=logger,
    )
    atomic_write_dataframe_csv(synthetic_truth_path, pd.DataFrame(recovery_rows), index=False)

    summary_payload = {
        "experiment_id": config["experiment"]["id"],
        "proteingym_version": proteingym_cfg["version"],
        "proteingym_package_version": installed_proteingym_version(),
        "mavenn_package_version": installed_mavenn_version(),
        "backend": {
            "abc": to_builtin(abc_backend.__dict__),
            "synthetic_truth_recovery": to_builtin(recovery_backend.__dict__),
        },
        "checkpoint_dir": str(checkpoint_dir),
        "log_path": str(output_dir / "progress.log"),
        "progress_path": str(output_dir / "progress.json"),
        "panel_csv_path": str(panel_csv_path),
        "mavenn_assay_metrics_path": str(mavenn_metrics_path),
        "branch_validation_csv_path": str(branch_validations_path),
        "posterior_particles_csv_path": str(posterior_particles_path),
        "posterior_rounds_csv_path": str(posterior_rounds_path),
        "posterior_parameter_summary_csv_path": str(posterior_parameter_summary_path),
        "target_features_csv_path": str(target_features_path),
        "synthetic_truth_recovery_csv_path": str(synthetic_truth_path),
        "panel": panel.summary_payload(),
        "target_features": {
            "feature_names": target.feature_names,
            "observed_vector": target.observed_vector.tolist(),
        },
        "branches": branch_summary,
        "smc_abc": {
            "selected_features": list(selected_features),
            "replicates_per_particle": int(abc_cfg["replicates_per_particle"]),
            "particle_count": int(abc_cfg["smc"]["n_particles"]),
            "round_count": int(abc_cfg["smc"]["n_rounds"]),
            "best_particle": {
                "parameters": to_builtin(best_updates),
                "distance": float(smc_result.best_particle.distance),
            },
            "posterior_mean": {
                "parameters": to_builtin(posterior_mean),
            },
            "round_summaries": [to_builtin(summary.__dict__) for summary in smc_result.round_summaries],
        },
        "synthetic_truth_recovery": recovery_rows,
    }
    summary_path = output_dir / "summary.json"
    atomic_write_json(summary_path, summary_payload)
    report(
        stage="run_complete",
        completed=1,
        total=1,
        message=f"Run completed; summary written to {summary_path}.",
        details={"summary_path": str(summary_path)},
    )
    print(summary_path)


if __name__ == "__main__":
    main()
