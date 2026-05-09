#!/usr/bin/env python3
"""Run shared SMC-ABC summary calibration for HYP-001."""

from __future__ import annotations

import argparse
import json
import sys
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
    run_synthetic_truth_recovery,
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
from project_adapt_env.smc_abc import SMCABCConfig, run_smc_abc  # noqa: E402
from project_adapt_env.utils import flatten_mapping, to_builtin  # noqa: E402

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
        rows.append(
            {
                "parameter": name,
                "weighted_mean": float(np.sum(weights * values)),
                "weighted_std": float(np.sqrt(np.sum(weights * (values - np.sum(weights * values)) ** 2))),
                "q05": float(np.quantile(values, 0.05)),
                "q50": float(np.quantile(values, 0.50)),
                "q95": float(np.quantile(values, 0.95)),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(Path(args.config))
    if args.quick and "quick" in config:
        config = merge_dict(config, config["quick"])

    proteingym_cfg = dict(config["proteingym"])
    panel_cfg = dict(config["panel"])
    mmseqs_cfg = dict(config["mmseqs"])
    mavenn_cfg = dict(config["mavenn"])
    base_config_dict = dict(config["adapt_env"]["base_config"])
    deterministic_branches = dict(config["deterministic_branches"])
    abc_cfg = dict(config["abc"])
    recovery_cfg = dict(config["synthetic_truth_recovery"])

    mmseqs_cfg["cache_dir"] = str(project_root / mmseqs_cfg["cache_dir"])
    panel = prepare_proteingym_panel(
        project_root=project_root,
        proteingym_cfg=proteingym_cfg,
        panel_cfg=panel_cfg,
        mmseqs_cfg=mmseqs_cfg,
        mavenn_cfg=mavenn_cfg,
        calibration_max_mutation_count=int(config["calibration_max_mutation_count"]),
    )

    panel_csv_path = output_dir / "selected_panel.csv"
    panel.panel_df.to_csv(panel_csv_path, index=False)

    mavenn_metrics_path = output_dir / "mavenn_assay_metrics.csv"
    panel.mavenn_metrics_frame().to_csv(mavenn_metrics_path, index=False)

    shared_branch_rows: list[dict[str, Any]] = []
    branch_summary: dict[str, Any] = {}
    base_config = LandscapeConfig(**base_config_dict)

    for branch_name, branch_cfg in deterministic_branches.items():
        source = str(branch_cfg["source"])
        if source != "raw":
            raise ValueError(f"Unsupported deterministic branch source: {source}")
        options = calibration_options_from_config(branch_cfg["calibration"])
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

    abc_options = calibration_options_from_config(abc_cfg["calibration"])
    selected_features = tuple(str(item) for item in abc_cfg.get("selected_features", [])) or None
    if selected_features is None:
        raise ValueError("abc.selected_features must be configured explicitly")

    target = build_empirical_target(
        landscapes=panel.raw_landscapes,
        alignment_paths=panel.alignment_paths,
        wildtypes=panel.wildtypes,
        options=abc_options,
        selected_features=selected_features,
        bootstrap_replicates=int(abc_cfg["bootstrap_replicates"]),
        covariance_shrinkage=float(abc_cfg["covariance_shrinkage"]),
        covariance_ridge=float(abc_cfg["covariance_ridge"]),
    )
    target_features_path = output_dir / "target_features.csv"
    bootstrap_sd = np.std(target.bootstrap_vectors, axis=0, ddof=1)
    pd.DataFrame(
        {
            "feature_name": target.feature_names,
            "observed_value": target.observed_vector,
            "bootstrap_sd": bootstrap_sd,
        }
    ).to_csv(target_features_path, index=False)

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
    )

    posterior_particles_path = output_dir / "posterior_particles.csv"
    pd.DataFrame(
        [
            {
                "distance": particle.distance,
                "weight": particle.weight,
                **{name: particle.parameters[name] for name in sorted(particle.parameters)},
            }
            for particle in smc_result.particles
        ]
    ).to_csv(posterior_particles_path, index=False)

    posterior_rounds_path = output_dir / "posterior_rounds.csv"
    pd.DataFrame([to_builtin(summary.__dict__) for summary in smc_result.round_summaries]).to_csv(
        posterior_rounds_path,
        index=False,
    )

    posterior_parameter_summary_path = output_dir / "posterior_parameter_summary.csv"
    pd.DataFrame(
        weighted_parameter_summary(
            smc_result.particles,
            parameter_names=[spec.name for spec in specs],
        )
    ).to_csv(posterior_parameter_summary_path, index=False)

    best_updates = dict(smc_result.best_particle.parameters)
    posterior_mean = posterior_mean_updates(specs=specs, particles=smc_result.particles)
    best_config, best_validation = validate_updates(
        empirical_collection=target.empirical_collection,
        base_config=base_config,
        options=abc_options,
        updates=best_updates,
    )
    mean_config, mean_validation = validate_updates(
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
                objective_terms={"distance": float(np.average(
                    np.asarray([particle.distance for particle in smc_result.particles], dtype=np.float64),
                    weights=np.asarray([particle.weight for particle in smc_result.particles], dtype=np.float64),
                ))},
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
    pd.DataFrame(shared_branch_rows).to_csv(branch_validations_path, index=False)

    recovery_rows = run_synthetic_truth_recovery(
        target=target,
        base_config=base_config,
        options=abc_options,
        specs=specs,
        smc_config=SMCABCConfig(**recovery_cfg["smc"]),
        truths=list(recovery_cfg["truths"]),
        replicates_per_particle=int(recovery_cfg["replicates_per_particle"]),
        selected_features=selected_features,
    )
    synthetic_truth_path = output_dir / "synthetic_truth_recovery.csv"
    pd.DataFrame(recovery_rows).to_csv(synthetic_truth_path, index=False)

    summary_payload = {
        "experiment_id": config["experiment"]["id"],
        "proteingym_version": proteingym_cfg["version"],
        "proteingym_package_version": installed_proteingym_version(),
        "mavenn_package_version": installed_mavenn_version(),
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
    summary_path.write_text(
        json.dumps(to_builtin(summary_payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(summary_path)


if __name__ == "__main__":
    main()
