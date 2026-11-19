#!/usr/bin/env python3
"""Guarded current-model ProteinGym ABC refit wrapper for EXP-029.

This script intentionally delegates the actual SMC-ABC implementation to
``proteingym_bayesian_summary_calibration.py`` so the calibration pathway is the
same one used by the earlier ProteinGym experiments. Its job is to enforce the
current model semantics before that handoff.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]

EXPECTED_ASSAYS = {
    "PHOT_CHLRE_Chen_2023",
    "SPTN1_CHICK_Tsuboyama_2023_1TUD",
}

PROHIBITED_PARAMETER_NAMES = {
    "n_functional_dims",
    "functional_sigma_base",
    "functional_sigma_anisotropy",
    "functional_weight",
    "peak_distance_from_consensus",
    "noise_amplitude",
    "noise_heterogeneity",
    "noise_taper_mode",
    "noise_taper_exponent",
    "strict_peak_lock",
    "peak_lock_margin",
    "n_higher_order_epistatic_terms",
    "epistasis_order_max",
    "epistasis_order_decay",
    "readout_stability_midpoint",
    "readout_stability_slope",
    "readout_function_exponent",
    "readout_function_midpoint",
    "readout_function_slope",
}

LEGACY_BASE_DEFAULTS: dict[str, Any] = {
    "n_functional_dims": 4,
    "functional_sigma_base": 20.0,
    "functional_sigma_anisotropy": 0.5,
    "functional_weight": 1.0,
    "peak_distance_from_consensus": 0,
    "noise_amplitude": 0.0,
    "noise_taper_mode": "uniform",
    "noise_taper_exponent": 1.0,
    "strict_peak_lock": False,
    "n_higher_order_epistatic_terms": 0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def merge_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_dict(dict(merged[key]), value)
        else:
            merged[key] = value
    return merged


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if quick and isinstance(config.get("quick"), dict):
        config = merge_dict(config, config["quick"])
    return config


def as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def assert_current_semantics(config: dict[str, Any]) -> dict[str, Any]:
    assay_ids = {str(item) for item in config.get("panel", {}).get("assay_ids", [])}
    if assay_ids != EXPECTED_ASSAYS:
        raise ValueError(
            "EXP-029 RUN must use exactly PHOT_CHLRE_Chen_2023 and "
            f"SPTN1_CHICK_Tsuboyama_2023_1TUD; got {sorted(assay_ids)}"
        )

    base = dict(config.get("adapt_env", {}).get("base_config", {}))
    combine_mode = str(base.get("observed_fitness_combine_mode", "")).strip().lower()
    if combine_mode != "stability_margin":
        raise ValueError(
            "Current-model ABC refit must use observed_fitness_combine_mode=stability_margin"
        )
    target = str(base.get("generic_epistasis_target", "")).strip().lower()
    if target != "stability":
        raise ValueError("Current epistasis must be routed into stability")

    for key, default in LEGACY_BASE_DEFAULTS.items():
        if key not in base:
            continue
        value = base[key]
        if isinstance(default, bool):
            if bool(value) != default:
                raise ValueError(f"Legacy base_config key {key} is active: {value!r}")
        elif isinstance(default, str):
            if str(value).strip().lower() != default:
                raise ValueError(f"Legacy base_config key {key} is active: {value!r}")
        else:
            numeric = as_float(value)
            if numeric is None or abs(numeric - float(default)) > 1e-12:
                raise ValueError(f"Legacy base_config key {key} is active: {value!r}")

    if base.get("latent_trait_blocks"):
        raise ValueError("Current-model ABC refit must not configure latent_trait_blocks")
    if base.get("observed_fitness_terms"):
        raise ValueError("Current-model ABC refit must not configure observed_fitness_terms")

    fitted_names = {
        str(item.get("name"))
        for item in config.get("abc", {}).get("parameters", [])
        if isinstance(item, dict)
    }
    prohibited_fit = sorted(fitted_names & PROHIBITED_PARAMETER_NAMES)
    if prohibited_fit:
        raise ValueError(f"Legacy parameters are prohibited in EXP-029 ABC fit: {prohibited_fit}")

    calibration = dict(config.get("abc", {}).get("calibration", {}))
    if str(calibration.get("synthetic_readout_mode", "")).strip().lower() != "stability_margin":
        raise ValueError("abc.calibration.synthetic_readout_mode must be stability_margin")
    if str(calibration.get("empirical_pairwise_target", "")).strip().lower() != "stability":
        raise ValueError("abc.calibration.empirical_pairwise_target must be stability")

    if config.get("deterministic_branches"):
        raise ValueError(
            "EXP-029 current-model refit must use the ABC pathway only; "
            "deterministic legacy branches must remain disabled"
        )

    return {
        "assay_ids": sorted(assay_ids),
        "fitted_parameter_names": sorted(fitted_names),
        "observed_fitness_combine_mode": combine_mode,
        "generic_epistasis_target": target,
        "legacy_parameter_guard": "passed",
    }


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(config_path, quick=args.quick)
    preflight = assert_current_semantics(config)
    (output_dir / "current_model_preflight.json").write_text(
        json.dumps(preflight, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    command = [
        sys.executable,
        "-u",
        str(project_root / "scripts" / "proteingym_bayesian_summary_calibration.py"),
        "--config",
        str(config_path),
        "--output-dir",
        str(output_dir),
        "--project-root",
        str(project_root),
    ]
    if args.quick:
        command.append("--quick")
    subprocess.run(command, check=True)

    summary_path = output_dir / "summary.json"
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        refit_summary = {
            "experiment_id": summary.get("experiment_id"),
            "preflight": preflight,
            "panel": summary.get("panel"),
            "target_features": summary.get("target_features"),
            "smc_abc": summary.get("smc_abc"),
            "source_summary_json": str(summary_path),
        }
        (output_dir / "empirical_refit_summary.json").write_text(
            json.dumps(refit_summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    mavenn_metrics = output_dir / "mavenn_assay_metrics.csv"
    if mavenn_metrics.is_file():
        shutil.copy2(mavenn_metrics, output_dir / "assay_refit_metrics.csv")


if __name__ == "__main__":
    main()
