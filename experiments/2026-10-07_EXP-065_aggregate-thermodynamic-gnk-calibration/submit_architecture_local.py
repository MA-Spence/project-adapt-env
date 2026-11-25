#!/usr/bin/env python3
"""Snapshot and submit EXP-065 architecture calibration through local labproj.

Default/--dry-run registers a nonexecuting preview. --execute launches a detached
supervisor and returns its PID, snapshot and durable log location. --reconcile
SNAPSHOT collects an interrupted supervisor's terminal RUN without rerunning it.
Scientific code and inputs execute from the verified snapshot. No DVC or remote
backup action is performed by this helper.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

import yaml

ROOT = Path(__file__).resolve().parents[2]
EXP_REL = Path(__file__).resolve().parent.relative_to(ROOT)
EXP = ROOT / EXP_REL
INPUT_NAME = "trpb-architecture-calibration-inputs"
METHOD_MANIFEST = Path("docs/model_release/0.7.0rc1/manifest.json")
ENV_LOCK = Path("environment/calibration-py311.lock.txt")
SOURCE_MANIFEST = "architecture_source.json"
CONFIG_NAME = "architecture_calibration.json"
COLLECTION_ROOT = Path("data/processed/aggregate-thermodynamic-gnk-calibration")
METHOD = "thornton_auxiliary_optuna_indirect_inference_v1"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Atomic status records remain readable while a detached supervisor updates.
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def load(path):
    return json.loads(Path(path).read_text())


def rooted(root, relative):
    path = Path(root) / relative
    if not path.resolve().is_relative_to(Path(root).resolve()):
        raise ValueError(f"Path escapes its recorded root: {relative}")
    return path


def verify_files(root, files):
    if not isinstance(files, dict) or not files:
        raise ValueError("A nonempty files mapping is required")
    for relative, record in files.items():
        expected = record["sha256"] if isinstance(record, dict) else record
        path = rooted(root, relative)
        if not path.is_file() or sha(path) != expected:
            raise ValueError(f"Frozen checksum mismatch: {path}")


def verify_frozen():
    source = load(EXP / SOURCE_MANIFEST)
    if source.get("file_count") != len(source["files"]):
        raise ValueError("Source manifest file_count must match its files mapping")
    verify_files(ROOT, source["files"])
    verify_files(ROOT, load(ROOT / METHOD_MANIFEST)["files"])
    inputs = ROOT / "data/processed" / INPUT_NAME
    verify_files(inputs, load(inputs / "manifest.json")["files"])
    for name in ("calibration.npz", "confirmation.npz"):
        if name not in load(inputs / "manifest.json")["files"]:
            raise ValueError(f"Input manifest must explicitly hash {name}")
    return source


def runtime_environment(snapshot):
    env = os.environ.copy()
    for key in list(env):
        if key.startswith("ADAPTENV_") or key in {
            "PYTHONPATH", "PYTHONHOME", "LABPROJ_OUTPUT_DIR", "LABPROJ_LOG_DIR",
            "LABPROJ_WORKDIR", "LABPROJ_RUN_ID", "LABPROJ_TARGET", "LABPROJ_QUICK",
        }:
            env.pop(key)
    python = ROOT / ".labproj/envs/indirect-py311/bin/python"
    if not python.is_file():
        raise ValueError(f"Prepared local interpreter is missing: {python}")
    tmp = snapshot / "tmp"
    (tmp / "dvc-site-cache").mkdir(parents=True, exist_ok=True)
    env.update(PYTHON=str(python), PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1",
               LABPROJ_QUICK="0", ADAPTENV_PROTOCOL="architecture",
               ADAPTENV_ARCHITECTURE_SNAPSHOT=str(snapshot), TMPDIR=str(tmp),
               DVC_SITE_CACHE_DIR=str(tmp / "dvc-site-cache"),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return env


def inspect_command(command, env):
    process = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                             text=True, timeout=120, check=False)
    result = dict(command=command, returncode=process.returncode,
                  stdout=process.stdout, stderr=process.stderr)
    if process.returncode:
        raise RuntimeError(f"Prelaunch inspection failed: {result}")
    return result


def capture_environment(snapshot, env):
    pinned = {}
    for line in (ROOT / ENV_LOCK).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            name, version = line.split("==", 1)
            pinned[name] = version
    probe = (
        "import importlib.metadata as m,json,platform,sys; "
        "print(json.dumps(dict(python=platform.python_version(),executable=sys.executable,"
        "platform=platform.platform(),python_build=platform.python_build(),"
        "packages={d.metadata['Name']:d.version for d in m.distributions()})))"
    )
    actual = json.loads(inspect_command([env["PYTHON"], "-c", probe], env)["stdout"])
    normalized = {re.sub(r"[-_.]+", "-", k).lower(): v for k, v in actual["packages"].items()}
    mismatches = {name: dict(expected=version, actual=normalized.get(re.sub(r"[-_.]+", "-", name).lower()))
                  for name, version in pinned.items()
                  if normalized.get(re.sub(r"[-_.]+", "-", name).lower()) != version}
    if actual["python"] != "3.11.15" or mismatches:
        raise ValueError(f"Frozen calibration environment mismatch: Python {actual['python']}; {mismatches}")
    actual["environment_lock_sha256"] = sha(ROOT / ENV_LOCK)
    actual["runtime_controls"] = {key: env[key] for key in (
        "PYTHON", "LABPROJ_QUICK", "ADAPTENV_PROTOCOL", "ADAPTENV_ARCHITECTURE_SNAPSHOT",
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "TMPDIR", "DVC_SITE_CACHE_DIR")}
    actual["nominal_stack"] = "STACK-000 is a record label; the captured native environment is used"
    directory = snapshot / "environment"
    directory.mkdir()
    dump(directory / "actual_environment.json", actual)
    config = inspect_command([env["PYTHON"], "-c", "import numpy; numpy.show_config()"], env)
    (directory / "numpy-config.txt").write_text(config["stdout"])
    (directory / "installed-packages.txt").write_text("\n".join(f"{k}=={v}" for k, v in sorted(actual["packages"].items())) + "\n")
    shutil.copy2(ROOT / ENV_LOCK, directory / ENV_LOCK.name)


def snapshot_sources(snapshot, source):
    method = load(ROOT / METHOD_MANIFEST)
    files = set(source["files"]) | set(method["files"])
    files.update(map(str, [METHOD_MANIFEST, EXP_REL / SOURCE_MANIFEST, EXP_REL / CONFIG_NAME,
                           EXP_REL / "metadata.yaml", EXP_REL / Path(__file__).name,
                           Path("labproject.yaml"), Path("AGENTS.md")]))
    for relative in sorted(files):
        origin = rooted(ROOT, relative)
        target = snapshot / "source" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(origin, target)
    shutil.copytree(ROOT / "data/processed" / INPUT_NAME, snapshot / "inputs" / INPUT_NAME)
    verify_files(snapshot / "source", source["files"])
    verify_files(snapshot / "source", method["files"])
    input_root = snapshot / "inputs" / INPUT_NAME
    verify_files(input_root, load(input_root / "manifest.json")["files"])
    captured = {str(p.relative_to(snapshot)): sha(p) for category in ("source", "inputs", "environment")
                for p in sorted((snapshot / category).rglob("*")) if p.is_file()}
    dump(snapshot / "prelaunch_manifest.json", dict(created_utc=datetime.now(timezone.utc).isoformat(),
         source_manifest_sha256=sha(EXP / SOURCE_MANIFEST),
         method_manifest_sha256=sha(ROOT / METHOD_MANIFEST),
         input_manifest_sha256=sha(input_root / "manifest.json"), files=captured))


def restore_overlay(snapshot):
    """Restore only our enabled flag if another process edited project metadata."""
    if (snapshot / "overlay-restored.json").exists():
        return
    original_path, patched_path = snapshot / "overlay-original.yaml", snapshot / "overlay-patched.yaml"
    if not original_path.exists() or not patched_path.exists():
        return
    path = ROOT / "labproject.yaml"
    original, patched, current = original_path.read_bytes(), patched_path.read_bytes(), path.read_bytes()
    if current == patched:
        path.write_bytes(original)
    elif current != original:
        state, before = yaml.safe_load(current), yaml.safe_load(original)
        overlay = state.setdefault("environment", {}).setdefault("overlay", {})
        if overlay.get("enabled") is False:
            prior = before.get("environment", {}).get("overlay", {})
            if "enabled" in prior:
                overlay["enabled"] = prior["enabled"]
            else:
                overlay.pop("enabled", None)
            path.write_text(yaml.safe_dump(state, sort_keys=False))
    dump(snapshot / "overlay-restored.json", dict(restored_utc=datetime.now(timezone.utc).isoformat(),
         note="Only the temporary overlay enabled setting was restored; unrelated current metadata preserved"))


def disable_overlay(snapshot):
    path = ROOT / "labproject.yaml"
    original = path.read_bytes()
    project = yaml.safe_load(original)
    project.setdefault("environment", {}).setdefault("overlay", {})["enabled"] = False
    patched = yaml.safe_dump(project, sort_keys=False).encode()
    (snapshot / "overlay-original.yaml").write_bytes(original)
    (snapshot / "overlay-patched.yaml").write_bytes(patched)
    path.write_bytes(patched)


def run_records():
    return {p for p in (EXP / "runs").glob("RUN-*.yaml") if re.fullmatch(r"RUN-\d+\.yaml", p.name)}


def protocol_run_metadata(run, snapshot):
    """The experiment's other protocols retain their own expected outputs."""
    metadata = yaml.safe_load((snapshot / "source" / EXP_REL / "metadata.yaml").read_text())
    expected = metadata.get("architecture_followup", {}).get("expected_outputs")
    if not expected:
        raise ValueError("architecture_followup.expected_outputs is required in frozen metadata")
    run.setdefault("outputs", {})["expected"] = expected
    run["architecture_protocol"] = dict(
        source_manifest_sha256=load(snapshot / "prelaunch_manifest.json")["source_manifest_sha256"],
        prelaunch_snapshot=str(snapshot.relative_to(ROOT)),
        method=METHOD, protocol="architecture")
    return run


def collection_run_metadata(run, collection):
    run["local_collection"] = collection
    run["verification_status"] = "passed" if collection["output_validation"]["passed"] else "failed"
    if not collection["output_validation"]["passed"]:
        run["status"] = "failed"
        run["failure_reason"] = dict(stage="output_verification",
                                     errors=collection["output_validation"]["errors"],
                                     runner_returncode=collection["runner_returncode"])
        if run.get("attempts"):
            run["attempts"][-1]["status"] = "failed"
    return run


def validate_outputs(output, snapshot):
    report = dict(passed=False, errors=[])
    try:
        required = ("result.json", "study.sqlite3", "trials.json", "provenance.json", "cli_inputs.json",
                    "config.json", "empirical_architecture.json", "empirical_coefficients.npz",
                    "confirmation_predictions.npz", "inputs.npz", "confirmation_inputs.npz")
        for name in required:
            if not (output / name).is_file() or (output / name).stat().st_size == 0:
                raise ValueError(f"Missing nonempty output: {name}")
        config = load(snapshot / "source" / EXP_REL / CONFIG_NAME)
        result, trials = load(output / "result.json"), load(output / "trials.json")
        if result.get("method") != METHOD:
            raise ValueError("Unexpected calibration method")
        if result.get("trials") != config["trials"] or len(trials) != config["trials"]:
            raise ValueError("Trial budget not completed")
        best = next((t for t in trials if t["number"] == result["best_trial"]), None)
        if best is None or best["state"] != "COMPLETE" or not result.get("completed_trials", 0):
            raise ValueError("Selected best trial did not complete")
        if not result.get("selected_inference_converged"):
            raise ValueError("Selected/empirical/bootstrap/confirmation auxiliary inference did not converge")
        if result.get("confirmation_scope") != "disjoint_genotypes_and_fresh_seeds":
            raise ValueError("Expected disjoint-genotype/fresh-seed confirmation")
        for key in ("calibration_loss", "confirmation_loss"):
            if not isinstance(result.get(key), (float, int)) or not math.isfinite(result[key]):
                raise ValueError(f"Invalid finite objective: {key}")
        for key in ("calibration_seeds", "confirmation_seeds"):
            if result.get(key) != config[key]:
                raise ValueError(f"Seed list differs from protocol: {key}")
        cli = load(output / "cli_inputs.json")
        for key, source_path in {
            "config": snapshot / "source" / EXP_REL / CONFIG_NAME,
            "calibration": snapshot / "inputs" / INPUT_NAME / "calibration.npz",
            "confirmation": snapshot / "inputs" / INPUT_NAME / "confirmation.npz",
        }.items():
            if cli[key]["sha256"] != sha(source_path):
                raise ValueError(f"CLI consumed different frozen {key} bytes")
        if cli["script_sha256"] != sha(snapshot / "source/scripts/calibrate_architecture.py"):
            raise ValueError("CLI script does not match frozen source")
        provenance = load(output / "provenance.json")
        for name, expected in provenance["source_sha256"].items():
            if sha(rooted(snapshot / "source/external/Adapt-Env", name)) != expected:
                raise ValueError(f"Method source provenance differs: {name}")
        report.update(passed=True, trials=len(trials), completed_trials=result["completed_trials"],
                      failed_trials=sum(t["state"] == "FAIL" for t in trials),
                      result_sha256=sha(output / "result.json"))
    except (OSError, ValueError, KeyError, TypeError, StopIteration) as error:
        report["errors"].append(str(error))
    return report


def collect(run_path, snapshot, cli_returncode):
    run = protocol_run_metadata(yaml.safe_load(run_path.read_text()), snapshot)
    run_id = run["id"]
    if run.get("target") != "local" or run.get("experiment") != "EXP-065" or run.get("dry_run"):
        raise ValueError("Only live local EXP-065 runs can be collected")
    bundle = ROOT / COLLECTION_ROOT / run_id
    if bundle.exists():
        verify_files(bundle, load(bundle / "collection_manifest.json")["files"])
        run = collection_run_metadata(run, load(bundle / "local_collection.json"))
        run_path.write_text(yaml.safe_dump(run, sort_keys=False))
        return load(bundle / "local_collection.json")
    verify_files(snapshot, load(snapshot / "prelaunch_manifest.json")["files"])
    output = rooted(ROOT, run["paths"]["output_dir"])
    validation = validate_outputs(output, snapshot)
    execution = run.get("execution_result", {})
    success = (cli_returncode == 0 and run.get("status") == "completed"
               and execution.get("returncode") == 0 and validation["passed"])
    bundle.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{run_id}-collecting-", dir=bundle.parent))
    for category in ("source", "inputs", "environment"):
        shutil.copytree(snapshot / category, staging / category)
    shutil.copy2(snapshot / "prelaunch_manifest.json", staging / "prelaunch_manifest.json")
    verify_files(staging, load(staging / "prelaunch_manifest.json")["files"])
    if output.exists():
        shutil.copytree(output, staging / "outputs")
    logs = staging / "logs"
    logs.mkdir()
    for stream in ("stdout", "stderr"):
        (logs / f"runner-{stream}.txt").write_text(execution.get(stream, ""))
        logpath = rooted(ROOT, run["logs"][stream])
        logpath.parent.mkdir(parents=True, exist_ok=True)
        if not logpath.exists():
            logpath.write_text(execution.get(stream, ""))
    for name in ("supervisor.log", "labproj.stdout.txt", "labproj.stderr.txt", "supervisor_failure.json"):
        if (snapshot / name).exists():
            shutil.copy2(snapshot / name, logs / name)
    for name in ("launch.json", "supervisor_state.json", "submission.json", "overlay-restored.json"):
        if (snapshot / name).exists():
            shutil.copy2(snapshot / name, staging / name)
    collection = dict(run=run_id, target="local", status="verified_complete" if success else "failure_evidence_preserved",
        bundle=str(bundle.relative_to(ROOT)), prelaunch_snapshot=str(snapshot.relative_to(ROOT)),
        cli_returncode=cli_returncode, runner_returncode=execution.get("returncode"), output_validation=validation,
        source_execution="Frozen scientific CLI, package, configuration and input paths; recorded live run.sh provides initial handoff and labproj registration",
        collection_method="Manual verified collection; local labproj collector lacks storage path mapping",
        stream_provenance="Runner streams from labproj execution_result (trailing whitespace stripped); supervisor log through collection time",
        dvc_tracking="pending_manual_tracking", remote_backup="not_attempted")
    run = collection_run_metadata(run, collection)
    records = staging / "run_record"
    records.mkdir()
    for path in (EXP / "runs").glob(run_id + ".*"):
        if path.is_file():
            shutil.copy2(path, records / path.name)
    (records / run_path.name).write_text(yaml.safe_dump(run, sort_keys=False))
    dump(staging / "local_collection.json", collection)
    files = {str(p.relative_to(staging)): sha(p) for p in sorted(staging.rglob("*")) if p.is_file()}
    dump(staging / "collection_manifest.json", dict(files=files))
    verify_files(staging, files)
    staging.rename(bundle)
    # Preserve unrelated metadata that might have been added while copying data.
    current = protocol_run_metadata(yaml.safe_load(run_path.read_text()), snapshot)
    current = collection_run_metadata(current, collection)
    run_path.write_text(yaml.safe_dump(current, sort_keys=False))
    return collection


def supervise(snapshot, execute):
    env = runtime_environment(snapshot)
    verify_files(snapshot, load(snapshot / "prelaunch_manifest.json")["files"])
    launch = load(snapshot / "launch.json")
    before = run_records()
    command = [launch["labproj"], "submit", "EXP-065", "--target", "local",
               "--execute" if execute else "--dry-run"]
    state = dict(pid=os.getpid(), status="submitting", snapshot=str(snapshot), execute=execute)
    dump(snapshot / "supervisor_state.json", state)
    disable_overlay(snapshot)
    try:
        with (snapshot / "labproj.stdout.txt").open("w") as stdout, (snapshot / "labproj.stderr.txt").open("w") as stderr:
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stdout, stderr=stderr)
            state["labproj_pid"] = process.pid
            last_heartbeat = 0.0
            while process.poll() is None:
                created = run_records() - before
                if len(created) == 1 and "run_record" not in state:
                    run_path = next(iter(created))
                    state.update(status="running" if execute else "preview", run=run_path.stem,
                                 run_record=str(run_path.relative_to(ROOT)))
                    dump(snapshot / "supervisor_state.json", state)
                    # labproj renders its scheduler script before writing the RUN;
                    # its overlay decision is now frozen in that script.
                    restore_overlay(snapshot)
                    print(f"Registered local {'live RUN' if execute else 'preview'}: {run_path.stem}", flush=True)
                if time.monotonic() - last_heartbeat >= 5:
                    state["heartbeat_utc"] = datetime.now(timezone.utc).isoformat()
                    if "run_record" in state:
                        try:
                            current = yaml.safe_load(rooted(ROOT, state["run_record"]).read_text())
                            output = rooted(ROOT, current["paths"]["output_dir"])
                            state["output_files"] = sorted(p.name for p in output.glob("*") if p.is_file())
                        except (OSError, KeyError, TypeError, yaml.YAMLError):
                            pass  # labproj may be writing its terminal record.
                    dump(snapshot / "supervisor_state.json", state)
                    last_heartbeat = time.monotonic()
                time.sleep(0.5)
    finally:
        restore_overlay(snapshot)
    created = run_records() - before
    if len(created) != 1:
        raise RuntimeError(f"Expected one new RUN; found {len(created)}; inspect {snapshot}")
    run_path = next(iter(created))
    run = yaml.safe_load(run_path.read_text())
    state.update(run=run["id"], run_record=str(run_path.relative_to(ROOT)), cli_returncode=process.returncode)
    dump(snapshot / "supervisor_state.json", state)
    dump(snapshot / "submission.json", dict(command=command, returncode=process.returncode, run_record=state["run_record"]))
    if execute:
        state["collection"] = collect(run_path, snapshot, process.returncode)
        success = state["collection"]["status"] == "verified_complete"
    else:
        success = process.returncode == 0 and run.get("dry_run") is True and run.get("target") == "local"
        run_path.write_text(yaml.safe_dump(protocol_run_metadata(run, snapshot), sort_keys=False))
    state["status"] = "completed" if success else "failed"
    dump(snapshot / "supervisor_state.json", state)
    print(json.dumps(state, indent=2), flush=True)
    return 0 if success else 1


def reconcile(snapshot):
    state = load(snapshot / "supervisor_state.json")
    for key in ("pid", "labproj_pid"):
        try:
            os.kill(int(state[key]), 0)
        except (ProcessLookupError, KeyError):
            continue
        raise RuntimeError(f"Supervisor or labproj PID {state[key]} is still running; do not collect mutable outputs")
    restore_overlay(snapshot)
    run_path = rooted(ROOT, state["run_record"])
    run = yaml.safe_load(run_path.read_text())
    if run.get("status") not in ("completed", "failed"):
        raise RuntimeError("RUN has no terminal execution record; inspect surviving processes before reconciliation")
    return collect(run_path, snapshot, state.get("cli_returncode", run.get("execution_result", {}).get("returncode")))


def main(argv=None):
    global ROOT, EXP
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--dry-run", action="store_true")
    action.add_argument("--execute", action="store_true")
    action.add_argument("--reconcile", type=Path)
    action.add_argument("--supervise", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--project-root", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.project_root is not None:
        ROOT = args.project_root.resolve()
        EXP = ROOT / EXP_REL
    if args.reconcile is not None:
        collection = reconcile(args.reconcile.resolve())
        print(json.dumps(collection, indent=2))
        return 0 if collection["status"] == "verified_complete" else 1
    if args.supervise is not None:
        snapshot = args.supervise.resolve()
        try:
            return supervise(snapshot, execute=True)
        except Exception as error:
            restore_overlay(snapshot)
            dump(snapshot / "supervisor_failure.json", dict(error=type(error).__name__, message=str(error)))
            state_path = snapshot / "supervisor_state.json"
            state = load(state_path) if state_path.exists() else {}
            state.update(status="supervisor_failed", error=type(error).__name__, message=str(error))
            dump(state_path, state)
            print(f"Detached local supervisor failed: {error}", file=sys.stderr, flush=True)
            raise
    source = verify_frozen()
    labproj = shutil.which("labproj")
    if labproj is None:
        raise ValueError("labproj is not on PATH")
    preparations = ROOT / ".labproj/local-preparations"
    preparations.mkdir(parents=True, exist_ok=True)
    snapshot = Path(tempfile.mkdtemp(prefix=datetime.now(timezone.utc).strftime("architecture-%Y%m%dT%H%M%SZ-"), dir=preparations))
    env = runtime_environment(snapshot)
    capture_environment(snapshot, env)
    snapshot_sources(snapshot, source)
    command = [env["PYTHON"], str(snapshot / "source/scripts/calibrate_architecture.py"),
               "--input", str(snapshot / "inputs" / INPUT_NAME / "calibration.npz"),
               "--confirmation", str(snapshot / "inputs" / INPUT_NAME / "confirmation.npz"),
               "--config", str(snapshot / "source" / EXP_REL / CONFIG_NAME), "--validate-only"]
    dump(snapshot / "input_validation.json", inspect_command(command, env))
    verify_frozen()
    dump(snapshot / "launch.json", dict(project_root=str(ROOT), labproj=labproj, execute=bool(args.execute),
         created_utc=datetime.now(timezone.utc).isoformat(), snapshot=str(snapshot)))
    if not args.execute:
        return supervise(snapshot, execute=False)
    command = [env["PYTHON"], str(snapshot / "source" / EXP_REL / Path(__file__).name),
               "--supervise", str(snapshot), "--project-root", str(ROOT)]
    caffeinate = shutil.which("caffeinate")
    if caffeinate:
        command = [caffeinate, "-i", *command]
    with (snapshot / "supervisor.log").open("ab") as logfile:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                                   stdout=logfile, stderr=subprocess.STDOUT, start_new_session=True, close_fds=True)
    dump(snapshot / "detached_launch.json", dict(pid=process.pid, command=command, snapshot=str(snapshot)))
    announcement = dict(status="detached_supervisor_started", pid=process.pid, snapshot=str(snapshot),
                        log=str(snapshot / "supervisor.log"), state=str(snapshot / "supervisor_state.json"))
    # Registration usually takes seconds; keep a hard bound if labproj's initial
    # inventory is slow. The supervisor state/heartbeat remains the durable view.
    deadline = time.monotonic() + 60
    while time.monotonic() < deadline and process.poll() is None:
        if (snapshot / "supervisor_state.json").exists():
            state = load(snapshot / "supervisor_state.json")
            if state.get("run"):
                announcement.update(run=state["run"], run_record=state["run_record"])
                break
        if (snapshot / "supervisor_failure.json").exists():
            raise RuntimeError(f"Detached supervisor failed; inspect {snapshot / 'supervisor_failure.json'}")
        time.sleep(0.5)
    if process.poll() not in (None, 0):
        raise RuntimeError(f"Detached supervisor exited {process.returncode}; inspect {snapshot / 'supervisor.log'}")
    print(json.dumps(announcement, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"Local architecture transport failed: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
