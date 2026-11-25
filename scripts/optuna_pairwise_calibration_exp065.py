#!/usr/bin/env python3
"""RUN-180 calibration and confirmation design with Optuna readout optimization.

Historical RUNs are immutable. This adapter reuses the existing graph, sampling,
profile and observation implementation and changes only the fitting procedure.
"""
from __future__ import annotations

import argparse
import copy
import importlib.metadata
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd
import yaml

import pairwise_calibration_exp065 as repeat
import optuna_readout_exp065 as fitting

ROOT = repeat.ROOT
EXP = ROOT / "experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration"
PROTOCOL = "trpb_run180_readout_optuna_v1"
EXTENSIONS = {"optimizer", "warm_start"}


def validate_config(config: dict, baseline: dict):
    """Permit the declared optimizer change, never a silent data/model change."""
    if config.get("protocol") != PROTOCOL:
        raise ValueError("Unexpected Optuna readout protocol")
    ignored = EXTENSIONS | {"protocol", "interpretation"}
    left = {k: v for k, v in config.items() if k not in ignored}
    right = {k: v for k, v in baseline.items() if k not in ignored}
    if left != right:
        keys = sorted(k for k in left.keys() | right.keys() if left.get(k) != right.get(k))
        raise ValueError(f"RUN-180 design changed: {keys}")
    if config["optimizer"] != dict(sampler="tpe", trials=1952,
                                   n_startup_trials=20, quick_trials=8):
        raise ValueError("Unexpected optimizer configuration; amend the protocol explicitly")
    if config["warm_start"]["run"] != "RUN-180":
        raise ValueError("Warm starts must be RUN-180 calibration settings")


def load_baseline(directory: Path, config: dict):
    manifest_path = directory / "manifest.json"
    if repeat.sweep.sha(manifest_path) != config["warm_start"]["manifest_sha256"]:
        raise ValueError("RUN-180 manifest checksum mismatch")
    manifest = json.loads(manifest_path.read_text())
    if manifest["mode"] != "MATCHED_PAIRWISE_REPEAT" or manifest["settings"] != 18:
        raise ValueError("Expected the complete RUN-180 baseline")
    # Authenticate all archived outputs, without using confirmation values to fit.
    for relative, expected in manifest["files"].items():
        path = directory / relative
        if (not path.resolve().is_relative_to(directory.resolve())
                or not path.is_file() or repeat.sweep.sha(path) != expected):
            raise ValueError(f"RUN-180 output checksum mismatch: {relative}")
    baseline = yaml.safe_load((directory / "config.yaml").read_text())
    validate_config(config, baseline)
    return json.loads((directory / "selected_parameters.json").read_text()), manifest


def audit_calibration(ca, cg, target, banks, origin, scale, floor, config,
                      baseline_root, warm_starts, *, latent_rtol=0.0, latent_atol=0.0):
    """Re-evaluate every prior fit on the actual calibration arrays before TPE."""
    if (not np.isfinite([latent_rtol, latent_atol]).all()
            or latent_rtol < 0 or latent_atol < 0):
        raise ValueError("Latent replay tolerances must be finite and nonnegative")
    strict = latent_rtol == 0 and latent_atol == 0
    differences = {}
    with np.load(baseline_root / "latent_calibration.npz") as old:
        for name, actual in [("a", ca), ("g", cg)]:
            if strict:
                np.testing.assert_array_equal(actual, old[name])
            else:
                np.testing.assert_allclose(actual, old[name],
                                           rtol=latent_rtol, atol=latent_atol)
            differences[name] = float(np.max(np.abs(actual - old[name])))
    rows = []
    for setting in warm_starts:
        for key, value in [("calibration_origin", origin), ("calibration_scale", scale),
                           ("near_floor_raw", floor)]:
            if setting[key] != value:
                raise ValueError(f"RUN-180 calibration units changed: {key}")
        objective = fitting.make_objective(
            setting["ratio"], ca, cg, target, banks[setting["observation"]],
            origin, scale, floor, config["fit"])
        actual = objective(setting["parameters"])
        for key, value in actual.items():
            np.testing.assert_allclose(value, setting["objective"][key], rtol=1e-12, atol=1e-14)
        rows.append(dict(ratio=setting["ratio"], observation=setting["observation"],
                         recorded=setting["objective"], replayed=actual))
    return dict(protocol=PROTOCOL, baseline_run="RUN-180", settings=len(rows),
                calibration_latents="array_equal" if strict else "allclose",
                latent_tolerance=dict(rtol=float(latent_rtol), atol=float(latent_atol)),
                maximum_absolute_latent_difference=differences, objective_replays=rows,
                evidence_scope="Engineering replay; not a new calibration or adequacy result")


def preflight(config: dict, data: Path, baseline_root: Path, matrix_path: Path,
              metadata_path: Path, *, audit_function=None):
    """Bounded replay only: no optimization, confirmation or graph computation."""
    warm_starts, manifest = load_baseline(baseline_root, config)
    inputs = json.loads((data / "manifest.json").read_text())
    if repeat.sweep.sha(data / "manifest.json") != manifest["input_manifest_sha256"]:
        raise ValueError("RUN-180 empirical input manifest changed")
    for relative, digest in inputs["files"].items():
        if repeat.sweep.sha(data / relative) != digest:
            raise ValueError(f"Empirical input checksum mismatch: {relative}")
    matrix, _ = repeat.load_prior(matrix_path, metadata_path, config)
    sub = data / "TrpB"
    index = pd.read_csv(sub / "latent_calibration_index.csv")
    with np.load(sub / "latent_calibration.npz") as latent:
        genotypes = latent["genotypes"]
    ca, cg = np.empty(genotypes.shape[:2]), np.empty(genotypes.shape[:2])
    for seed in index.seed.unique():
        selected = np.flatnonzero(index.seed.eq(seed))
        model = repeat.build_model(int(seed), config["model"], matrix)
        assert model.config.interaction_concentration is None
        assert model.config.global_intercept is None
        a, g = model.latent_components(genotypes[selected].reshape(-1, 4))
        ca[selected], cg[selected] = a.reshape(-1, 16), g.reshape(-1, 16)
    origin, scale, floor = (warm_starts[0][k] for k in
                           ["calibration_origin", "calibration_scale", "near_floor_raw"])
    with np.load(sub / "empirical_panel.npz") as empirical:
        split = pd.read_csv(sub / "empirical_panel_index.csv").split.eq("calibration")
        target = (empirical["fitness"][split] - origin) / scale
    original = yaml.safe_load((data / "config.yaml").read_text())
    with tempfile.TemporaryDirectory(prefix="readout-noise-") as temp:
        noise, _ = repeat.sweep.calibrate_noise(
            pd.read_csv(sub / "donor_assignments.csv.gz"),
            set(json.loads((sub / "excluded_noise_donors.json").read_text())),
            original, Path(temp))
    with np.load(sub / "fit_observation_uniforms.npz") as values:
        uniforms = values["uniforms"]
    banks = {"none": None}
    for mode in config["observation_modes"]["TrpB"][1:]:
        banks[mode] = repeat.sweep.observation_bank(noise, uniforms, mode == "shared_stop_sensitivity")
    return (audit_function or audit_calibration)(
        ca, cg, target, banks, origin, scale, floor, config, baseline_root, warm_starts)


def run(config_path: Path, output: Path, data: Path, matrix_path: Path,
        metadata_path: Path, source_manifest: Path, historical_manifest: Path,
        profile_control_root: Path, baseline_root: Path, workers=4, quick=False):
    config = yaml.safe_load(config_path.read_text())
    warm_starts, _ = load_baseline(baseline_root, config)
    if importlib.metadata.version("optuna") != "4.9.0":
        raise ValueError("Use the recorded Optuna 4.9.0 environment")

    def validate_profile(current, profile):
        base = copy.deepcopy(current)
        for key in EXTENSIONS:
            base.pop(key)
        repeat.validate_profile_config(base, profile)

    def fit(ca, cg, target, banks, origin, scale, floor, current, original, fits, quick=False):
        audit = audit_calibration(ca, cg, target, banks, origin, scale, floor,
                                  current, baseline_root, warm_starts)
        repeat.sweep.dump(output / "run180_replay.json", audit)
        shutil.copyfile(baseline_root / "selected_parameters.json", output / "warm_start_parameters.json")
        shutil.copyfile(baseline_root / "manifest.json", output / "warm_start_manifest.json")
        repeat.sweep.dump(output / "optimizer_environment.json", dict(
            packages={name: importlib.metadata.version(name) for name in
                      ["optuna", "SQLAlchemy", "alembic"]}))
        return fitting.fit_settings(ca, cg, target, banks, origin, scale, floor,
                                    current, original, fits, quick,
                                    optimizer=current["optimizer"], warm_starts=warm_starts)

    repeat.run(config_path, output, data, matrix_path, metadata_path,
               source_manifest, historical_manifest, profile_control_root, workers, quick,
               fit_settings_fn=fit, profile_config_validator=validate_profile,
               output_mode="RUN180_PROTOCOL_OPTUNA",
               intervention_metadata=dict(
                   changed_component="readout optimizer: differential evolution to Optuna TPE",
                   baseline_run="RUN-180", optimizer=config["optimizer"],
                   warm_start=config["warm_start"],
                   endpoint_policy="Optimize all r values; preserve previous same-r and r0/r4 fits as candidates",
                   coefficient_distribution="Gaussian (interaction_concentration=None)",
                   sigmoid_position="reference anchored (global_intercept=None)"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXP / "optuna_readout.yaml")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, default=ROOT / "data/processed/interaction-strength-exp065-inputs")
    parser.add_argument("--prior-matrix", type=Path, default=ROOT / "data/processed/trpb-profile-additive-inputs/additive_matrix.npy")
    parser.add_argument("--prior-metadata", type=Path, default=ROOT / "data/processed/trpb-profile-additive-inputs/prior_manifest.json")
    parser.add_argument("--source-manifest", type=Path, default=EXP / "optuna_readout_source.json")
    parser.add_argument("--historical-source-manifest", type=Path, default=EXP / "interaction_strength_source.json")
    parser.add_argument("--profile-control-root", type=Path, default=ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178/outputs")
    parser.add_argument("--baseline-root", type=Path, default=ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/outputs")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--preflight", action="store_true", help="Replay old calibration only, without fitting")
    args = parser.parse_args()
    if args.preflight:
        if args.output_dir.exists() and any(args.output_dir.iterdir()):
            parser.error("Use a new or empty output directory")
        config = yaml.safe_load(args.config.read_text())
        report = preflight(config, args.data_root, args.baseline_root,
                           args.prior_matrix, args.prior_metadata)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        repeat.sweep.dump(args.output_dir / "run180_replay.json", report)
        print("PASS: all 18 RUN-180 objectives and calibration fields replayed")
    else:
        kwargs = vars(args)
        kwargs.pop("preflight")
        kwargs["config_path"] = kwargs.pop("config")
        kwargs["output"] = kwargs.pop("output_dir")
        kwargs["data"] = kwargs.pop("data_root")
        kwargs["matrix_path"] = kwargs.pop("prior_matrix")
        kwargs["metadata_path"] = kwargs.pop("prior_metadata")
        kwargs["historical_manifest"] = kwargs.pop("historical_source_manifest")
        run(**kwargs)


if __name__ == "__main__":
    main()
