#!/usr/bin/env python3
"""Fit a synthetic landscape ensemble to inferred reference-free architecture.

This is a new calibration workflow. It does not rewrite or replay EXP-065's
historical fitness-distribution objective. See docs/calibration/README.md.
"""
from __future__ import annotations

import argparse
from dataclasses import fields
import hashlib
import io
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "external" / "Adapt-Env"
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_BUNDLE_KEYS = {"sequences", "fitness", "additive_matrix", "folds"}


def _api():
    # Select this checkout, not an unrelated editable installation. Delay imports
    # so --help can describe the workflow without the calibration extra installed.
    sys.path.insert(0, str(PACKAGE))
    from adaptenv.calibration import CalibrationConfig, CalibrationData, calibrate
    from adaptenv.inference import RFAConfig
    return CalibrationConfig, CalibrationData, RFAConfig, calibrate


def _reject_constant(value: str):
    raise ValueError(f"JSON configuration contains non-finite constant {value}")


def _instantiate(cls, values: dict[str, Any], label: str):
    if not isinstance(values, dict):
        raise ValueError(f"{label} must be a JSON object")
    allowed = {field.name for field in fields(cls) if field.init}
    unknown = set(values) - allowed
    if unknown:
        raise ValueError(f"Unknown {label} fields: {', '.join(sorted(unknown))}")
    # JSON has arrays rather than tuples. Restore immutable sequence fields.
    converted = {
        key: tuple(value) if isinstance(value, list) and (
            key.endswith("_bounds") or key.endswith("_seeds")
            or key == "penalty_grid"
        ) else value
        for key, value in values.items()
    }
    return cls(**converted)


def load_config(path: Path, config_cls, inference_cls):
    payload = path.read_bytes()
    raw = json.loads(payload, parse_constant=_reject_constant)
    if not isinstance(raw, dict):
        raise ValueError("Configuration must be a JSON object")
    values = dict(raw)
    if "inference" in values:
        values["inference"] = _instantiate(inference_cls, values["inference"], "inference")
    config = _instantiate(config_cls, values, "calibration")
    return config, _file_provenance(path, payload)


def _file_provenance(path: Path, payload: bytes) -> dict:
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(payload).hexdigest(),
            "bytes": len(payload)}


def load_bundle(path: Path, data_cls, *, require_additive: bool):
    """Read and hash the same bytes; never allow pickled arrays in input files."""
    payload = path.read_bytes()
    with np.load(io.BytesIO(payload), allow_pickle=False) as archive:
        unknown = set(archive.files) - _BUNDLE_KEYS
        if unknown:
            raise ValueError(f"Unknown NPZ arrays in {path}: {', '.join(sorted(unknown))}")
        required = {"sequences", "fitness"}
        if require_additive:
            required.add("additive_matrix")
        missing = required - set(archive.files)
        if missing:
            raise ValueError(f"Missing NPZ arrays in {path}: {', '.join(sorted(missing))}")
        arrays = {key: archive[key] for key in archive.files}
    data = data_cls(sequences=arrays["sequences"], fitness=arrays["fitness"],
                    folds=arrays.get("folds"))
    matrix = arrays.get("additive_matrix")
    if matrix is not None:
        if (matrix.shape != (data.sequences.shape[1], 20)
                or matrix.dtype.kind not in "fiu" or not np.isfinite(matrix).all()):
            raise ValueError("additive_matrix must be finite numeric L x 20 in AA_ORDER")
        if not np.any(np.ptp(matrix, axis=1) > 0):
            raise ValueError("additive_matrix must contain nonzero within-site variation")
        matrix = np.asarray(matrix, dtype=np.float64)
    provenance = _file_provenance(path, payload)
    provenance.update(sequences=int(data.sequences.shape[0]),
                      sites=int(data.sequences.shape[1]), alphabet=ALPHABET)
    return data, matrix, provenance


def _guard_output(path: Path):
    destination = path.resolve()
    archive = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration"
    for run in ("RUN-175", "RUN-178", "RUN-180"):
        if destination.is_relative_to((archive / run).resolve()):
            raise ValueError(f"Output must not modify the archived {run} record")
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise ValueError("Output must be a new or empty directory; resuming is not supported")


def _validate_panels(data, confirmation, config):
    scale = float(np.nanstd(data.fitness))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Calibration fitness must have nonzero finite variance")
    if config.noise_mode == "replicate_residual":
        if not np.any(np.isfinite(data.fitness).sum(axis=1) >= 2):
            raise ValueError("replicate_residual requires genotypes with >=2 measured replicates")
    if confirmation is not None:
        if confirmation.sequences.shape[1] != data.sequences.shape[1]:
            raise ValueError("Calibration and confirmation must use the same sites")
        if set(map(tuple, data.sequences)) & set(map(tuple, confirmation.sequences)):
            raise ValueError("Confirmation genotypes must be disjoint from calibration")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input", type=Path, required=True,
                        help="Calibration NPZ: sequences, fitness, additive_matrix, optional folds")
    result.add_argument("--config", type=Path, required=True,
                        help="JSON configuration; unknown fields are rejected")
    result.add_argument("--confirmation", type=Path,
                        help="Optional independently held-out NPZ: sequences, fitness, optional folds")
    result.add_argument("--output", type=Path, help="New or empty output directory")
    result.add_argument("--validate-only", action="store_true",
                        help="Check input and configuration structure without fitting or writing outputs")
    return result


def main(argv: list[str] | None = None) -> int:
    arguments = parser()
    args = arguments.parse_args(argv)
    if not args.validate_only and args.output is None:
        arguments.error("--output is required unless --validate-only is used")
    try:
        config_cls, data_cls, inference_cls, fit = _api()
        config, config_provenance = load_config(args.config, config_cls, inference_cls)
        data, matrix, input_provenance = load_bundle(args.input, data_cls, require_additive=True)
        confirmation = None
        confirmation_provenance = None
        if args.confirmation is not None:
            confirmation, confirmation_matrix, confirmation_provenance = load_bundle(
                args.confirmation, data_cls, require_additive=False)
            if confirmation.sequences.shape[1] != data.sequences.shape[1]:
                raise ValueError("Calibration and confirmation must use the same sites")
            if (confirmation_matrix is not None
                    and not np.array_equal(matrix, confirmation_matrix)):
                raise ValueError("Confirmation additive_matrix must match the calibration profile exactly")
        _validate_panels(data, confirmation, config)
        if args.output is not None:
            _guard_output(args.output)
        if args.validate_only:
            print(json.dumps({"status": "structurally_valid", "calibration": input_provenance,
                              "confirmation": confirmation_provenance,
                              "config": config_provenance,
                              "note": "No fitting, scientific acceptance, or output creation performed."},
                             indent=2))
            return 0
        fit(data, matrix, config, output_dir=args.output, confirmation_data=confirmation)
        provenance = {"calibration": input_provenance, "confirmation": confirmation_provenance,
                      "config": config_provenance, "script_sha256": hashlib.sha256(
                          Path(__file__).read_bytes()).hexdigest(),
                      "argv": argv if argv is not None else sys.argv[1:]}
        (args.output / "cli_inputs.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print(f"Calibration complete. Artifacts: {args.output.resolve()}")
        return 0
    except (ValueError, TypeError, OSError, ImportError) as error:
        arguments.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
