#!/usr/bin/env python3
"""Stage RUN-182's unchanged method/inputs and submit EXP-065 to Slurm.

The default is a nonexecuting labproj preview. Only --execute submits a job.
The dedicated Linux environment is prepared separately; neither this helper nor
the worker installs packages. The scheduler starts a fresh identical-budget fit,
not a continuation of the interrupted local Optuna database.
"""
from __future__ import annotations

import argparse
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
from contextlib import redirect_stdout

import yaml

ROOT = Path(__file__).resolve().parents[2]
EXP_REL = Path(__file__).resolve().parent.relative_to(ROOT)
EXP = ROOT / EXP_REL
BUNDLE_REL = Path("data/processed/trpb-architecture-slurm-runtime")
FROZEN_MANIFEST_SHA256 = "4e00fd05143ff27e905ca99050400283cd5c1acec94d55858eee29beb0a02004"
PYTHON = "/tank/data/labproj/runtime/project-adapt-env/envs/architecture-0.7.0rc1-py311/bin/python"
TRANSPORT_NAME = "architecture_slurm_transport.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def verify_bundle(bundle):
    bundle = Path(bundle).resolve()
    manifest_path = bundle / "prelaunch_manifest.json"
    if sha(manifest_path) != FROZEN_MANIFEST_SHA256:
        raise ValueError("The RUN-182 prelaunch manifest changed")
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest.get("files"), dict) or not manifest["files"]:
        raise ValueError("A nonempty prelaunch files mapping is required")
    for relative, expected in manifest["files"].items():
        path = bundle / relative
        if not path.resolve().is_relative_to(bundle) or path.is_symlink():
            raise ValueError(f"Invalid frozen path: {relative}")
        if not path.is_file() or sha(path) != expected:
            raise ValueError(f"Frozen checksum mismatch: {relative}")
    return manifest


def worker_preflight(bundle, log_dir):
    """Verify frozen bytes and actual Linux runtime without executing a fit."""
    bundle, log_dir = Path(bundle), Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    manifest = verify_bundle(bundle)
    lock = bundle / "source/environment/calibration-py311.lock.txt"
    pinned = {}
    for line in lock.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            name, version = line.split("==", 1)
            pinned[name] = version
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
        "source_run": "RUN-182", "execution_mode": "fresh identical-budget Linux restart",
        "platform_change": "macOS arm64 to Linux x86_64; BLAS/runtime recorded, no numerical identity claimed",
        "runtime_controls": {name: os.environ.get(name) for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "LABPROJ_RUN_ID", "SLURM_JOB_ID", "LABPROJ_OUTPUT_DIR")},
    }
    dump(log_dir / "architecture_linux_environment.json", actual)
    if actual["python"] != "3.11.15" or platform.system() != "Linux" or mismatches:
        raise ValueError(f"Dedicated Linux environment differs from the frozen pins: {mismatches}; Python={actual['python']}")
    import numpy
    stream = io.StringIO()
    with redirect_stdout(stream):
        numpy.show_config()
    (log_dir / "architecture_linux_numpy_config.txt").write_text(stream.getvalue())
    print(f"Verified {len(manifest['files'])} RUN-182 payloads and {len(pinned)} pinned packages", flush=True)


def staging_filter(root, files, trees):
    """Keep only self-contained runtime and small labproj transport records."""
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
    return "# Frozen RUN-182 architecture calibration and Slurm transport only.\n" + "\n".join(sorted(excludes)) + "\n"


def staged_dispatcher(original):
    first, separator, rest = original.partition(b"\n")
    if not first.startswith(b"#!") or not separator:
        raise ValueError("Experiment run.sh must begin with a shebang")
    if b'== "architecture-slurm"' not in original:
        raise ValueError("Experiment run.sh lacks its architecture-slurm handoff")
    return first + b"\nexport ADAPTENV_PROTOCOL=architecture-slurm\n" + rest


def submit(execute=False):
    manifest = verify_bundle(ROOT / BUNDLE_REL)
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
        "source_run": "RUN-182", "target": "lab-slurm", "ssh_transport": "dutch",
        "protocol": "architecture-slurm", "dedicated_python": PYTHON,
        "runtime_bundle": BUNDLE_REL.as_posix(),
        "prelaunch_manifest_sha256": FROZEN_MANIFEST_SHA256,
        "verified_payload_count": len(manifest["files"]),
        "scientific_bytes": "Unchanged original RUN-182 prelaunch snapshot",
        "execution_mode": "fresh identical-budget restart; not SQLite resume",
        "original_run_sh_sha256": hashlib.sha256(originals[run_script]).hexdigest(),
        "staged_run_sh_sha256": hashlib.sha256(dispatcher).hexdigest(),
        "transport_files_sha256": {name: sha(EXP / name) for name in (
            "submit_architecture_slurm.py", "architecture_slurm_worker.sh")},
    }
    dump(EXP / TRANSPORT_NAME, record)
    files = {str(EXP_REL / name) for name in (
        "run.sh", "metadata.yaml", "submit_architecture_slurm.py",
        "architecture_slurm_worker.sh", "architecture_slurm.md",
        "architecture_slurm_storage.json", TRANSPORT_NAME)} | {
        "labproject.yaml", ".labproj-stageignore", ".gitignore", ".gitmodules",
        "AGENTS.md", "AIM.md", "PROJECT_STATE.md", "docs/hypotheses.md",
        ".agents/skills/labproj/SKILL.md", str(BUNDLE_REL) + ".dvc"}
    trees = {BUNDLE_REL.as_posix(), "environment", "registry", str(EXP_REL / "runs"),
             str(EXP_REL / "architecture_remote_environment")}
    filter_text = staging_filter(ROOT, files, trees)
    (EXP / "architecture_slurm.staging-filter.txt").write_text(filter_text)
    before = set((EXP / "runs").glob("RUN-*.yaml"))
    try:
        (ROOT / ".labproj-stageignore").write_text(filter_text)
        (ROOT / "labproject.yaml").write_text(yaml.safe_dump(project, sort_keys=False))
        (ROOT / ".env").write_bytes(transport)
        run_script.write_bytes(dispatcher)
        command = ["labproj", "submit", "EXP-065", "--target", "lab-slurm",
                   "--execute" if execute else "--dry-run"]
        print("Frozen RUN-182 runtime; Dutch transport; " + ("queue submission" if execute else "nonexecuting preview"), flush=True)
        result = subprocess.run(command, cwd=ROOT, timeout=600, check=False)
        return result.returncode
    finally:
        # Rendering/staging/submission is bounded; restore all temporary settings.
        for path, content in originals.items():
            path.write_bytes(content)
        # Persist the exact submitted transport record beside every new RUN.
        for path in set((EXP / "runs").glob("RUN-*.yaml")) - before:
            if path.name.endswith(".manifest.yaml"):
                continue
            dump(path.with_suffix(".architecture-slurm.json"), record)
        print("Restored experiment dispatcher, project overlay, SSH transport and staging filter.", flush=True)


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
        print(f"Verified {len(verify_bundle(ROOT / BUNDLE_REL)['files'])} frozen RUN-182 payloads")
        return 0
    return submit(args.execute)


if __name__ == "__main__":
    sys.exit(main())
