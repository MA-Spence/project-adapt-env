#!/usr/bin/env python3
"""Run the frozen EXP-065 profile repeat through labproj's local transport.

Default: create a nonexecuting local preview and a verified prelaunch snapshot.
--execute creates a new live RUN; it never resubmits the Slurm preview. Scientific
source/configuration are unchanged. DVC tracking and remote backup are separate.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import yaml


ROOT = Path(__file__).resolve().parents[2]
EXP = Path(__file__).resolve().parent
BUNDLES = {
    "interaction-strength-exp065-inputs": "manifest.json",
    "trpb-profile-additive-inputs": "prior_manifest.json",
}
COLLECTION_NOTE = (
    "Manual verified local collection: labproj scheduler=none executes the live "
    "workspace without source staging and captures streams in execution_result. "
    "Its collect path mapper requires client_data_storage_root_template, absent "
    "from the local target. Sources and complete inputs were snapshotted before "
    "submission; captured streams, outputs and RUN sidecars are preserved here. "
    "DVC tracking and remote backup are not performed by this helper."
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def rooted(root: Path, relative: str) -> Path:
    path = root / relative
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes its recorded root: {relative}")
    return path


def verify_files(directory: Path, files: dict) -> None:
    for relative, expected in files.items():
        if isinstance(expected, dict):
            expected = expected["sha256"]
        path = rooted(directory, relative)
        if not path.is_file() or sha(path) != expected:
            raise ValueError(f"Frozen checksum mismatch: {path}")


def verify_frozen() -> dict:
    source = json.loads((EXP / "alignment_prior_source.json").read_text())
    if len(source["files"]) != 42:
        raise ValueError("Expected the frozen 42-file scientific source manifest")
    verify_files(ROOT, source["files"])
    for name, manifest_name in BUNDLES.items():
        directory = ROOT / "data/processed" / name
        manifest = json.loads((directory / manifest_name).read_text())
        verify_files(directory, manifest["files"])
    config = yaml.safe_load((EXP / "alignment_prior.yaml").read_text())
    prior = config["additive_prior"]
    directory = ROOT / prior["bundle"]
    for file_key, hash_key in [("matrix_file", "matrix_sha256"),
                               ("metadata_file", "metadata_sha256")]:
        if sha(directory / prior[file_key]) != prior[hash_key]:
            raise ValueError(f"Configured profile checksum mismatch: {file_key}")
    return source


def command_record(command: list[str], env: dict, *, required=True) -> dict:
    result = subprocess.run(command, cwd=ROOT, env=env, text=True,
                            capture_output=True, timeout=120, check=False)
    record = dict(command=command, returncode=result.returncode,
                  stdout=result.stdout, stderr=result.stderr)
    if required and result.returncode:
        raise RuntimeError(f"Environment inspection failed: {command!r}\n{result.stderr}")
    return record


def capture_environment(python: Path, env: dict, destination: Path) -> dict:
    probe = (
        "import importlib.metadata,json,platform,sys; "
        "print(json.dumps(dict(python=platform.python_version(), "
        "packages={n:importlib.metadata.version(n) for n in "
        "['numpy','pandas','scipy','pyarrow','PyYAML']}, "
        "executable=sys.executable,platform=platform.platform(), "
        "machine=platform.machine(),python_build=platform.python_build())))"
    )
    actual = json.loads(command_record([str(python), "-c", probe], env)["stdout"])
    required = yaml.safe_load((EXP / "alignment_prior.yaml").read_text())["required_worker_environment"]
    if {key: actual[key] for key in ("python", "packages")} != required:
        raise ValueError(f"Local Python/numerical environment does not match RUN-175: {actual}")
    compiler = shutil.which("clang++", path=env.get("PATH"))
    if compiler is None:
        raise ValueError("Local graph kernel requires clang++ on PATH")
    actual["compiler"] = command_record([compiler, "--version"], env)
    actual["compiler"]["executable"] = compiler
    xcrun = shutil.which("xcrun", path=env.get("PATH"))
    if xcrun is None:
        raise ValueError("Local macOS graph kernel requires xcrun")
    actual["sdk"] = command_record([xcrun, "--show-sdk-path"], env)
    actual["nominal_stack"] = "STACK-000; provenance label only, not activated locally"
    actual["actual_environment"] = "Explicit project-local Python environment; native macOS"
    actual["runtime_controls"] = {key: env[key] for key in (
        "PYTHON", "LABPROJ_QUICK", "ADAPTENV_PROTOCOL", "TMPDIR",
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    destination.mkdir(parents=True, exist_ok=True)
    freeze = command_record([str(python), "-m", "pip", "freeze", "--all"], env)
    numpy_config = command_record([str(python), "-c", "import numpy; numpy.show_config()"], env)
    (destination / "pip-freeze.txt").write_text(freeze["stdout"])
    (destination / "numpy-config.txt").write_text(numpy_config["stdout"])
    dump(destination / "inspection_commands.json", dict(pip_freeze=freeze, numpy_config=numpy_config))
    dump(destination / "actual_environment.json", actual)
    return actual


def snapshot_sources(snapshot: Path, source: dict) -> None:
    metadata = yaml.safe_load((EXP / "metadata.yaml").read_text())
    files = set(source["files"]) | set(metadata["inputs"]["declared_files"])
    files.update(str(path.relative_to(ROOT)) for path in (
        EXP / "metadata.yaml", EXP / "alignment_prior_local.md", Path(__file__),
        EXP / "alignment_prior_local_validation.json", EXP / "alignment_prior_local_tests.log",
        EXP / "alignment_prior_source.json", EXP / "alignment_prior_validation.json",
        ROOT / "labproject.yaml", ROOT / "AGENTS.md"))
    for relative in sorted(files):
        source_path = rooted(ROOT, relative)
        if not source_path.is_file():
            raise ValueError(f"Prelaunch source/context file is missing: {relative}")
        target = snapshot / "source" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, target)
    for name in BUNDLES:
        shutil.copytree(ROOT / "data/processed" / name, snapshot / "inputs" / name)
    verify_files(snapshot / "source", source["files"])
    for name, manifest_name in BUNDLES.items():
        directory = snapshot / "inputs" / name
        verify_files(directory, json.loads((directory / manifest_name).read_text())["files"])
    captured = {str(path.relative_to(snapshot)): sha(path)
                for category in ("source", "inputs", "environment")
                for path in sorted((snapshot / category).rglob("*")) if path.is_file()}
    dump(snapshot / "prelaunch_manifest.json", dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        source_manifest_sha256=sha(EXP / "alignment_prior_source.json"),
        frozen_source_count=42, declared_file_count=len(metadata["inputs"]["declared_files"]),
        files=captured))


def run_records() -> set[Path]:
    return {path for path in (EXP / "runs").glob("RUN-*.yaml")
            if re.fullmatch(r"RUN-\d+\.yaml", path.name)}


def validate_outputs(output: Path) -> dict:
    result = dict(passed=False, errors=[])
    try:
        manifest = json.loads((output / "manifest.json").read_text())
        verify_files(output, manifest["files"])
        for key, expected in dict(mode="MATCHED_PRIOR_REPEAT", settings=18, seed_audits=24,
                                  graph_rows=2163, latent_graph_rows=144,
                                  cube_summary_rows=2163).items():
            if manifest.get(key) != expected:
                raise ValueError(f"Unexpected {key}: {manifest.get(key)!r}; expected {expected!r}")
        for key, path in {
            "config_sha256": EXP / "alignment_prior.yaml",
            "source_manifest_sha256": EXP / "alignment_prior_source.json",
            "input_manifest_sha256": ROOT / "data/processed/interaction-strength-exp065-inputs/manifest.json",
        }.items():
            if manifest[key] != sha(path):
                raise ValueError(f"Output provenance mismatch: {key}")
        required = yaml.safe_load((EXP / "alignment_prior.yaml").read_text())["required_worker_environment"]
        if json.loads((output / "worker_environment.json").read_text()) != required:
            raise ValueError("Output worker environment does not match frozen requirements")
        result.update(passed=True, files_verified=len(manifest["files"]),
                      settings=18, seed_audits=24, manifest_sha256=sha(output / "manifest.json"))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        result["errors"].append(str(exc))
    return result


def collect(run_path: Path, snapshot: Path, cli_returncode: int) -> dict:
    run = yaml.safe_load(run_path.read_text())
    run_id = run["id"]
    output = rooted(ROOT, run["paths"]["output_dir"])
    execution = run.get("execution_result", {})
    for stream in ("stdout", "stderr"):
        path = rooted(ROOT, run["logs"][stream])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(execution.get(stream, ""))
    validation = validate_outputs(output)
    success = (cli_returncode == 0 and run.get("status") == "completed"
               and execution.get("returncode") == 0 and validation["passed"])
    bundle = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration" / run_id
    bundle.mkdir(parents=True, exist_ok=False)
    for category in ("source", "inputs", "environment"):
        shutil.copytree(snapshot / category, bundle / category)
    shutil.copy2(snapshot / "prelaunch_manifest.json", bundle / "prelaunch_manifest.json")
    verify_files(bundle, json.loads((bundle / "prelaunch_manifest.json").read_text())["files"])
    if output.exists():
        shutil.copytree(output, bundle / "outputs")
    if validation["passed"]:
        verify_files(bundle / "outputs", json.loads((bundle / "outputs/manifest.json").read_text())["files"])
    log_dir = rooted(ROOT, run["paths"]["log_dir"])
    shutil.copytree(log_dir, bundle / "logs")
    for name in ("labproj.stdout.txt", "labproj.stderr.txt"):
        shutil.copy2(snapshot / name, bundle / "logs" / name)
    collection = dict(
        run=run_id, target="local", status="verified_complete" if success else "failure_evidence_preserved",
        explanation=COLLECTION_NOTE, bundle=str(bundle.relative_to(ROOT)),
        prelaunch_snapshot=str(snapshot.relative_to(ROOT)),
        cli_returncode=cli_returncode, runner_returncode=execution.get("returncode"),
        output_validation=validation,
        actual_environment=str((bundle / "environment/actual_environment.json").relative_to(ROOT)),
        stream_provenance="Extracted labproj execution_result strings; labproj strips trailing whitespace",
        dvc_tracking="pending", remote_backup="pending_Dutch_unavailable")
    run["local_collection"] = collection
    run_path.write_text(yaml.safe_dump(run, sort_keys=False))
    records = bundle / "run_record"
    records.mkdir()
    for path in (EXP / "runs").glob(run_id + ".*"):
        if path.is_file():
            shutil.copy2(path, records / path.name)
    dump(bundle / "local_collection.json", collection)
    copied = {str(path.relative_to(bundle)): sha(path)
              for path in sorted(bundle.rglob("*")) if path.is_file()}
    dump(bundle / "collection_manifest.json", dict(files=copied))
    verify_files(bundle, copied)
    return collection


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--dry-run", action="store_true")
    action.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    source = verify_frozen()
    # Keep the venv path: resolving its Python symlink would bypass the venv.
    python = ROOT / ".labproj/envs/trpb-py31115/bin/python"
    if not python.is_file():
        raise ValueError(f"Prepared local interpreter is missing: {python}")
    labproj = shutil.which("labproj")
    if labproj is None:
        raise ValueError("labproj is not on PATH")
    preparations = ROOT / ".labproj/local-preparations"
    preparations.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ-")
    snapshot = Path(tempfile.mkdtemp(prefix=stamp, dir=preparations))
    tmpdir = snapshot / "tmp"
    tmpdir.mkdir()
    env = os.environ.copy()
    for key in list(env):
        if key.startswith("ADAPTENV_") or key in {
            "PYTHONPATH", "PYTHONHOME", "LABPROJ_OUTPUT_DIR", "LABPROJ_LOG_DIR",
            "LABPROJ_WORKDIR", "LABPROJ_RUN_ID", "LABPROJ_TARGET", "LABPROJ_QUICK",
        }:
            env.pop(key)
    env.update(PYTHON=str(python), PYTHONDONTWRITEBYTECODE="1", LABPROJ_QUICK="0",
               ADAPTENV_PROTOCOL="alignment-prior", TMPDIR=str(tmpdir),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    capture_environment(python, env, snapshot / "environment")
    snapshot_sources(snapshot, source)
    verify_frozen()
    print(f"Verified prelaunch snapshot: {snapshot}", flush=True)
    before = run_records()
    project_path = ROOT / "labproject.yaml"
    original = project_path.read_bytes()
    project = yaml.safe_load(original)
    project["environment"]["overlay"]["enabled"] = False
    command = [labproj, "submit", "EXP-065", "--target", "local",
               "--execute" if args.execute else "--dry-run"]
    announced = set()
    try:
        project_path.write_text(yaml.safe_dump(project, sort_keys=False))
        process = subprocess.Popen(command, cwd=ROOT, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
            for path in sorted(run_records() - before - announced):
                print(f"Local {'live RUN' if args.execute else 'preview'}: {path.stem}", flush=True)
                announced.add(path)
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                continue
        (snapshot / "labproj.stdout.txt").write_text(stdout)
        (snapshot / "labproj.stderr.txt").write_text(stderr)
    finally:
        project_path.write_bytes(original)
        print("Restored project overlay settings.", flush=True)
    created = run_records() - before
    if len(created) != 1:
        raise RuntimeError(f"Expected exactly one new local RUN, found {len(created)}; logs: {snapshot}")
    run_path = created.pop()
    run = yaml.safe_load(run_path.read_text())
    if run.get("target") != "local" or run.get("experiment") != "EXP-065":
        raise ValueError(f"Unexpected newly created RUN: {run_path}")
    if not args.execute:
        dump(snapshot / "preview.json", dict(run=run["id"], target="local", dry_run=True,
                                             cli_returncode=process.returncode,
                                             run_record=str(run_path.relative_to(ROOT))))
        print(f"Nonexecuting local preview: {run['id']}; snapshot retained at {snapshot}", flush=True)
        return 0 if process.returncode == 0 and run.get("dry_run") else 1
    collection = collect(run_path, snapshot, process.returncode)
    print(json.dumps(collection, indent=2), flush=True)
    return 0 if collection["status"] == "verified_complete" else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"Local transport failed: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
