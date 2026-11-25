#!/usr/bin/env python3
"""Joint continuous-r Optuna calibration using RUN-180's direct fitness loss."""
from __future__ import annotations

import argparse
import copy
import importlib.metadata
from pathlib import Path
import shutil

import yaml

import pairwise_calibration_exp065 as repeat
import optuna_pairwise_calibration_exp065 as rollback
import optuna_joint_readout_exp065 as fitting

ROOT, EXP = rollback.ROOT, rollback.EXP
PROTOCOL = "trpb_run180_joint_r_readout_optuna_v1"
HISTORICAL_RATIOS = [0, 0.5, 1, 2, 4, 8]


def historical_projection(config):
    """Validate all unmodified design fields through the existing rollback audit."""
    if config.get("protocol") != PROTOCOL or "interaction_ratios" in config:
        raise ValueError("Joint calibration must not contain a fitting grid")
    if config.get("ratio_bounds") != [0.0, 8.0]:
        raise ValueError("Joint r bounds must be [0, 8]")
    if config.get("optimizer") != dict(sampler="tpe", trials=11712,
                                       n_startup_trials=20, quick_trials=8):
        raise ValueError("Unexpected joint optimization budget or sampler")
    base = copy.deepcopy(config)
    base.pop("ratio_bounds")
    base["interaction_ratios"] = HISTORICAL_RATIOS.copy()
    base["protocol"] = rollback.PROTOCOL
    base["optimizer"]["trials"] = 1952
    return base


def load_baseline(directory, config):
    return rollback.load_baseline(directory, historical_projection(config))


def audit_calibration(ca, cg, target, banks, origin, scale, floor, config,
                      baseline_root, warm_starts):
    """Allow measured cross-platform roundoff, retaining the historical loss audit.

    Linux replay differed from the macOS pair field by at most 1.34e-15;
    additive fields were identical. These tolerances cover numerical roundoff
    without changing fields, normalization, or the objective replay tolerance.
    """
    return rollback.audit_calibration(
        ca, cg, target, banks, origin, scale, floor, config, baseline_root,
        warm_starts, latent_rtol=1e-13, latent_atol=1e-14)


def preflight(config, data, baseline_root, matrix_path, metadata_path):
    report = rollback.preflight(historical_projection(config), data, baseline_root,
                                matrix_path, metadata_path, audit_function=audit_calibration)
    report.update(protocol=PROTOCOL, ratio_bounds=config["ratio_bounds"],
                  historical_ratios_role="Warm starts only; no fitting grid")
    return report


def run(config_path, output, data, matrix_path, metadata_path, source_manifest,
        historical_manifest, profile_control_root, baseline_root, workers=4, quick=False):
    config = yaml.safe_load(config_path.read_text())
    warm_starts, _ = load_baseline(baseline_root, config)
    if importlib.metadata.version("optuna") != "4.9.0":
        raise ValueError("Use the recorded Optuna 4.9.0 environment")

    def validate_design(current, historical):
        repeat.validate_config(historical_projection(current), historical)

    def validate_profile(current, profile):
        base = historical_projection(current)
        for key in rollback.EXTENSIONS:
            base.pop(key)
        repeat.validate_profile_config(base, profile)

    def fit(ca, cg, target, banks, origin, scale, floor, current, original, fits, quick=False):
        audit = audit_calibration(ca, cg, target, banks, origin, scale,
            floor, current, baseline_root, warm_starts)
        audit.update(protocol=PROTOCOL, ratio_bounds=current["ratio_bounds"],
                     historical_ratios_role="Warm starts only; no fitting grid")
        repeat.sweep.dump(output / "run180_replay.json", audit)
        shutil.copyfile(baseline_root / "selected_parameters.json", output / "warm_start_parameters.json")
        shutil.copyfile(baseline_root / "manifest.json", output / "warm_start_manifest.json")
        repeat.sweep.dump(output / "optimizer_environment.json", dict(
            packages={name: importlib.metadata.version(name) for name in
                      ["optuna", "SQLAlchemy", "alembic"]}))
        return fitting.fit_settings(ca, cg, target, banks, origin, scale, floor,
            current, original, fits, quick, optimizer=current["optimizer"], warm_starts=warm_starts)

    repeat.run(config_path, output, data, matrix_path, metadata_path,
        source_manifest, historical_manifest, profile_control_root, workers, quick,
        fit_settings_fn=fit, profile_config_validator=validate_profile,
        config_validator=validate_design, output_mode="RUN180_JOINT_R_OPTUNA",
        intervention_metadata=dict(
            changed_component="Continuous joint Optuna optimization of r and four readout coordinates",
            baseline_run="RUN-180", optimizer=config["optimizer"],
            warm_start=config["warm_start"], ratio_bounds=config["ratio_bounds"],
            endpoint_policy="All six historical (r, readout) fits are candidate starts per observation mode",
            coefficient_distribution="Gaussian (interaction_concentration=None)",
            sigmoid_position="reference anchored (global_intercept=None)"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXP / "joint_readout.yaml")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=ROOT / "data/processed/interaction-strength-exp065-inputs")
    parser.add_argument("--prior-matrix", type=Path, default=ROOT / "data/processed/trpb-profile-additive-inputs/additive_matrix.npy")
    parser.add_argument("--prior-metadata", type=Path, default=ROOT / "data/processed/trpb-profile-additive-inputs/prior_manifest.json")
    parser.add_argument("--source-manifest", type=Path, default=EXP / "joint_readout_source.json")
    parser.add_argument("--historical-source-manifest", type=Path, default=EXP / "interaction_strength_source.json")
    parser.add_argument("--profile-control-root", type=Path, default=ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178/outputs")
    parser.add_argument("--baseline-root", type=Path, default=ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/outputs")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        if args.output_dir.exists() and any(args.output_dir.iterdir()):
            parser.error("Use a new or empty output directory")
        report = preflight(yaml.safe_load(args.config.read_text()), args.data_root,
                           args.baseline_root, args.prior_matrix, args.prior_metadata)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        repeat.sweep.dump(args.output_dir / "run180_replay.json", report)
        print("PASS: all 18 historical warm-start objectives replayed; joint r is continuous [0,8]")
    else:
        kwargs = vars(args)
        kwargs.pop("preflight")
        for old, new in [("config", "config_path"), ("output_dir", "output"),
                         ("data_root", "data"), ("prior_matrix", "matrix_path"),
                         ("prior_metadata", "metadata_path"),
                         ("historical_source_manifest", "historical_manifest")]:
            kwargs[new] = kwargs.pop(old)
        run(**kwargs)


if __name__ == "__main__":
    main()
