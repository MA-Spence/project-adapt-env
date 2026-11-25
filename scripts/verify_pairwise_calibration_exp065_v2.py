#!/usr/bin/env python3
"""EXP-065 validator v2: distinguish donor gzip provenance from donor identity.

The frozen v1 validator is unchanged and hash-pinned below. Its sole amended
comparison is noise_model.json/membership_sha256: each value must authenticate
its own compressed donor file, whose decompressed bytes must match the control
exactly. Every other noise-model field and every other v1 check are retained.
This amendment addresses gzip mtime headers, not changed scientific inputs.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import time

V1_SHA256 = "cc4779d56639f34935543213a7305f5bfe1e352c1ba6e40bcf3b946b52f09ecd"
V1_PATH = Path(__file__).with_name("verify_pairwise_calibration_exp065.py")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_v1():
    if sha(V1_PATH) != V1_SHA256:
        raise ValueError("The frozen v1 validator changed; v2 requires its pinned implementation")
    spec = importlib.util.spec_from_file_location("_exp065_pairwise_frozen_v1", V1_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_membership_amendment(output, control):
    """Check raw provenance and exact decompressed bytes before excluding one key."""
    paths = [Path(root)/"noise" for root in (output, control)]
    models, entries = [], []
    for directory in paths:
        model = json.loads((directory/"noise_model.json").read_text())
        compressed = directory/"donor_assignments.csv.gz"
        digest = sha(compressed)
        if model.get("membership_sha256") != digest:
            raise ValueError(f"membership_sha256 does not match its own compressed donor file: {compressed}")
        models.append(model)
        entries.append(dict(file=str(compressed.resolve()), compressed_sha256=digest))
    other_fields = [{k: v for k, v in model.items() if k != "membership_sha256"}
                    for model in models]
    if other_fields[0] != other_fields[1]:
        raise ValueError("Noise model fields other than membership_sha256 changed")
    digest = hashlib.sha256()
    byte_count = 0
    with gzip.open(paths[0]/"donor_assignments.csv.gz", "rb") as actual, \
            gzip.open(paths[1]/"donor_assignments.csv.gz", "rb") as expected:
        while True:
            left, right = actual.read(1024 * 1024), expected.read(1024 * 1024)
            if left != right:
                raise ValueError("Decompressed donor membership bytes differ")
            if not left:
                break
            digest.update(left)
            byte_count += len(left)
    return dict(candidate=entries[0], control=entries[1],
                decompressed_bytes_equal=True, decompressed_bytes=byte_count,
                decompressed_sha256=digest.hexdigest(),
                other_noise_model_fields_equal=True,
                permitted_difference="membership_sha256 only, after authenticating each gzip file")


def validate_run(output, control, report):
    baseline = load_v1()
    report["frozen_v1"] = dict(path=str(V1_PATH.resolve()), sha256=V1_SHA256)
    report["membership_amendment"] = verify_membership_amendment(output, control)
    noise_paths = {(Path(root)/"noise/noise_model.json").resolve()
                   for root in (output, control)}
    original_read_json = baseline.read_json

    def read_with_verified_membership(path):
        value = original_read_json(path)
        if Path(path).resolve() in noise_paths:
            # Only the two files already authenticated and compared above are
            # normalized. Raw output manifests still hash their original bytes.
            return {k: v for k, v in value.items() if k != "membership_sha256"}
        return value

    baseline.read_json = read_with_verified_membership
    try:
        baseline.validate_run(output, control, report)
    finally:
        baseline.read_json = original_read_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--control-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    report = dict(schema="exp065_pairwise_engineering_validation_v2", passed=False,
                  scope="Engineering reproducibility only; no statistical adequacy or superiority claim",
                  amendment="Authenticate membership_sha256 against each own gzip file; require exact decompressed donor bytes and all remaining noise-model fields to match. All other frozen-v1 checks retained.",
                  output_dir=str(args.output_dir.resolve()), control_dir=str(args.control_dir.resolve()),
                  tolerances=dict(rtol=1e-10, atol=1e-12, additive_atol=1e-14),
                  validator_sha256=sha(Path(__file__)))
    try:
        validate_run(args.output_dir, args.control_dir, report)
        report["passed"] = True
    except (KeyError, ValueError, OSError, TypeError, EOFError) as error:
        report["error"] = f"{type(error).__name__}: {error}"
    report["elapsed_seconds"] = time.perf_counter()-started
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(passed=report["passed"], report=str(args.report),
                         elapsed_seconds=report["elapsed_seconds"], error=report.get("error"))))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
