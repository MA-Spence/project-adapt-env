#!/usr/bin/env python3
"""Verify a frozen gamma-concentration calibration and submit it through labproj.

The default is a nonexecuting preview. Only --execute submits the job. Source,
inputs and the environment lock are verified before staging and on the worker;
the dedicated Linux environment must be prepared separately.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
EXP_REL = Path(__file__).resolve().parent.relative_to(ROOT)
EXP = ROOT / EXP_REL
BUNDLE_REL = Path("data/processed/trpb-gamma-concentration-slurm-runtime")
PYTHON = "/tank/data/labproj/runtime/project-adapt-env/envs/joint-readout-py311/bin/python"
RUNTIME_NAME = "gamma_concentration_runtime.json"
TRANSPORT_NAME = "gamma_concentration_slurm_transport.json"
PROTOCOL = "gamma-concentration-slurm"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def runtime_pin():
    record = json.loads((EXP / RUNTIME_NAME).read_text())
    digest = record.get("manifest_sha256")
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError("The gamma-concentration runtime requires a SHA-256 manifest pin")
    if record.get("runtime_bundle", BUNDLE_REL.as_posix()) != BUNDLE_REL.as_posix():
        raise ValueError("The runtime pin names a different bundle")
    return record


def verify_bundle(bundle):
    """Reject altered, missing, undeclared or redirected frozen payloads."""
    bundle = Path(bundle).resolve()
    manifest_path = bundle / "prelaunch_manifest.json"
    if manifest_path.is_symlink() or sha(manifest_path) != runtime_pin()["manifest_sha256"]:
        raise ValueError("The gamma-concentration prelaunch manifest changed")
    manifest = json.loads(manifest_path.read_text())
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError("A nonempty prelaunch files mapping is required")
    for relative, expected in files.items():
        key = Path(relative)
        path = bundle / key
        if (key.is_absolute() or ".." in key.parts or key.as_posix() != relative
                or relative == "prelaunch_manifest.json"
                or not path.resolve().is_relative_to(bundle)
                or any((bundle / Path(*key.parts[:i])).is_symlink()
                       for i in range(1, len(key.parts) + 1))):
            raise ValueError(f"Invalid frozen path: {relative}")
        if (not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected)
                or not path.is_file() or sha(path) != expected):
            raise ValueError(f"Frozen checksum mismatch: {relative}")
    actual = {path.relative_to(bundle).as_posix() for path in bundle.rglob("*")
              if path.is_file() or path.is_symlink()}
    if actual != set(files) | {"prelaunch_manifest.json"}:
        raise ValueError("Frozen runtime has missing or undeclared payloads")
    return manifest


def worker_preflight(bundle, log_dir):
    """Verify source/input bytes and record the actual pinned Linux runtime."""
    bundle, log_dir = Path(bundle), Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    manifest = verify_bundle(bundle)
    lock = bundle / "source/environment/readout-py311.lock.txt"
    pinned = {}
    for raw in lock.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            match = re.fullmatch(r"([A-Za-z0-9_.-]+)==([^\s]+)", line)
            if match is None:
                raise ValueError("Worker environment lock must use exact package pins")
            name, version = match.groups()
            if name in pinned:
                raise ValueError("Worker environment lock has duplicate package pins")
            pinned[name] = version
    if not pinned:
        raise ValueError("Worker environment lock is empty")
    packages = {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()}
    normalize = lambda name: re.sub(r"[-_.]+", "-", name).lower()
    normalized = {normalize(name): version for name, version in packages.items()}
    mismatches = {name: {"expected": version, "actual": normalized.get(normalize(name))}
                  for name, version in pinned.items()
                  if normalized.get(normalize(name)) != version}
    actual = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "executable": sys.executable,
        "platform": platform.platform(), "machine": platform.machine(),
        "python_build": platform.python_build(), "packages": packages,
        "environment_lock_sha256": sha(lock), "package_mismatches": mismatches,
        "prelaunch_manifest_sha256": sha(bundle / "prelaunch_manifest.json"),
        "verified_payload_count": len(manifest["files"]),
        "protocol": PROTOCOL, "historical_control": "RUN-188",
        "execution_mode": "joint Gamma alpha/r/readout Optuna TPE fit; controls reused",
        "platform_change": "RUN-186 frozen Linux inputs; local macOS smoke; no cross-platform bitwise identity claimed",
        "runtime_controls": {name: os.environ.get(name) for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "LABPROJ_RUN_ID", "SLURM_JOB_ID", "LABPROJ_OUTPUT_DIR")},
    }
    dump(log_dir / "gamma_concentration_linux_environment.json", actual)
    if actual["python"] != "3.11.15" or platform.system() != "Linux" or mismatches:
        raise ValueError(f"Dedicated Linux environment differs from the frozen pins: {mismatches}; Python={actual['python']}")
    import numpy
    stream = io.StringIO()
    with redirect_stdout(stream):
        numpy.show_config()
    (log_dir / "gamma_concentration_linux_numpy_config.txt").write_text(stream.getvalue())
    print(f"Verified {len(manifest['files'])} frozen payloads and {len(pinned)} pinned packages", flush=True)


def staging_filter(root, files, trees):
    """Keep the self-contained runtime and the small labproj transport record."""
    excludes = []
    allowed = set(files) | set(trees)

    def inspect(folder):
        for entry in os.scandir(folder):
            path = Path(entry.path)
            relative = path.relative_to(root).as_posix()
            if relative in allowed:
                continue
            if entry.is_dir(follow_symlinks=False) and any(
                    item.startswith(relative + "/") for item in allowed):
                inspect(path)
            else:
                excludes.append(relative)
    inspect(root)
    return "# Frozen gamma-concentration calibration and Slurm transport only.\n" + "\n".join(sorted(excludes)) + "\n"


def staged_dispatcher(original):
    first, separator, rest = original.partition(b"\n")
    if not first.startswith(b"#!") or not separator:
        raise ValueError("Experiment run.sh must begin with a shebang")
    if b'== "gamma-concentration-slurm"' not in original:
        raise ValueError("Experiment run.sh lacks its gamma-concentration-slurm handoff")
    return first + b"\nexport ADAPTENV_PROTOCOL=gamma-concentration-slurm\n" + rest


def submit(execute=False):
    manifest = verify_bundle(ROOT / BUNDLE_REL)
    pin = runtime_pin()
    run_script = EXP / "run.sh"
    paths = [ROOT / ".labproj-stageignore", ROOT / "labproject.yaml", ROOT / ".env", run_script]
    originals = {path: path.read_bytes() for path in paths}
    project = yaml.safe_load(originals[ROOT / "labproject.yaml"])
    project["environment"]["overlay"]["enabled"] = False
    transport = re.sub(rb"(?m)^LAB_SSH_HOST=.*$", b"LAB_SSH_HOST=dutch", originals[ROOT / ".env"])
    if not re.search(rb"(?m)^LAB_SSH_HOST=dutch\s*$", transport):
        raise ValueError("Cannot select the recorded Dutch SSH transport")
    dispatcher = staged_dispatcher(originals[run_script])
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "historical_control": "RUN-188", "target": "lab-slurm", "ssh_transport": "dutch",
        "protocol": PROTOCOL, "dedicated_python": PYTHON,
        "runtime_bundle": BUNDLE_REL.as_posix(),
        "prelaunch_manifest_sha256": pin["manifest_sha256"],
        "verified_payload_count": len(manifest["files"]),
        "scientific_change": "Joint Gamma alpha on [0.05,100] log-scale, r on [0,8] and readout; frozen shared-stop law and RUN-180 loss",
        "execution_mode": "fresh calibration, not SQLite resume",
        "original_run_sh_sha256": hashlib.sha256(originals[run_script]).hexdigest(),
        "staged_run_sh_sha256": hashlib.sha256(dispatcher).hexdigest(),
        "transport_files_sha256": {name: sha(EXP / name) for name in (
            "submit_gamma_concentration_slurm.py", "gamma_concentration_slurm_worker.sh", RUNTIME_NAME)},
    }
    dump(EXP / TRANSPORT_NAME, record)
    files = {str(EXP_REL / name) for name in (
        "run.sh", "metadata.yaml", "submit_gamma_concentration_slurm.py",
        "gamma_concentration_slurm_worker.sh", "gamma_concentration.md", "gamma_concentration.yaml",
        "gamma_concentration_source.json", "gamma_concentration_validation.json",
        "gamma_concentration_slurm_storage.json", RUNTIME_NAME, TRANSPORT_NAME)} | {
        "labproject.yaml", ".labproj-stageignore", ".gitignore", ".gitmodules",
        "AGENTS.md", "AIM.md", "PROJECT_STATE.md", "docs/hypotheses.md",
        ".agents/skills/labproj/SKILL.md", str(BUNDLE_REL) + ".dvc"}
    trees = {BUNDLE_REL.as_posix(), "environment", "registry", str(EXP_REL / "runs"),
             str(EXP_REL / "gamma_concentration_remote_environment")}
    filter_text = staging_filter(ROOT, files, trees)
    (EXP / "gamma_concentration_slurm.staging-filter.txt").write_text(filter_text)
    before = set((EXP / "runs").glob("RUN-*.yaml"))
    try:
        (ROOT / ".labproj-stageignore").write_text(filter_text)
        (ROOT / "labproject.yaml").write_text(yaml.safe_dump(project, sort_keys=False))
        (ROOT / ".env").write_bytes(transport)
        run_script.write_bytes(dispatcher)
        command = ["labproj", "submit", "EXP-065", "--target", "lab-slurm",
                   "--execute" if execute else "--dry-run"]
        print("Frozen gamma-concentration runtime; Dutch transport; " +
              ("queue submission" if execute else "nonexecuting preview"), flush=True)
        return subprocess.run(command, cwd=ROOT, timeout=600, check=False).returncode
    finally:
        for path, content in originals.items():
            path.write_bytes(content)
        for path in set((EXP / "runs").glob("RUN-*.yaml")) - before:
            if not path.name.endswith(".manifest.yaml"):
                dump(path.with_suffix(".gamma-concentration-slurm.json"), record)
        print("Restored dispatcher, overlay, SSH transport and staging filter.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--dry-run", action="store_true")
    action.add_argument("--execute", action="store_true")
    action.add_argument("--verify-only", action="store_true")
    action.add_argument("--worker-preflight", type=Path, metavar="BUNDLE")
    parser.add_argument("--log-dir", type=Path)
    args = parser.parse_args()
    if args.worker_preflight:
        if args.log_dir is None:
            parser.error("--worker-preflight requires --log-dir")
        worker_preflight(args.worker_preflight, args.log_dir)
        return 0
    if args.verify_only:
        print(f"Verified {len(verify_bundle(ROOT / BUNDLE_REL)['files'])} frozen gamma-concentration payloads")
        return 0
    return submit(args.execute)


if __name__ == "__main__":
    sys.exit(main())
