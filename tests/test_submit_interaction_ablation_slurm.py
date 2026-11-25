"""Frozen Slurm transport integrity and restoration; never submits real jobs."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/submit_interaction_ablation_slurm.py"
SPEC = importlib.util.spec_from_file_location("interaction_ablation_slurm_transport", SCRIPT)
transport = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(transport)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) if isinstance(value, dict) else value)
    return path


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    root = tmp_path / "project"
    exp_rel = Path("experiments/EXP-065")
    exp = root / exp_rel
    monkeypatch.setattr(transport, "ROOT", root)
    monkeypatch.setattr(transport, "EXP_REL", exp_rel)
    monkeypatch.setattr(transport, "EXP", exp)
    bundle = root / transport.BUNDLE_REL
    payloads = [write(bundle / "source/scripts/model.py", "frozen source\n"),
                write(bundle / "source/data/processed/input/data.npz", "frozen observations"),
                write(bundle / "source/environment/readout-py311.lock.txt", "numpy==2.4.4\npytest==8.4.2\n")]
    manifest = {"files": {path.relative_to(bundle).as_posix(): transport.sha(path)
                          for path in payloads}}
    manifest_path = write(bundle / "prelaunch_manifest.json", manifest)
    write(exp / transport.RUNTIME_NAME, {"manifest_sha256": transport.sha(manifest_path),
                                        "runtime_bundle": transport.BUNDLE_REL.as_posix()})
    write(root / "labproject.yaml", yaml.safe_dump({"environment": {"overlay": {"enabled": True}}}))
    write(root / ".env", "LAB_SSH_HOST=raspi-1\nPRIVATE_VALUE=do-not-print\n")
    write(root / ".labproj-stageignore", "original exclusion\n")
    write(exp / "run.sh", '#!/usr/bin/env bash\nif [[ "$ADAPTENV_PROTOCOL" == "interaction-ablation-slurm" ]]; then\n  exit 0\nfi\n')
    write(exp / "submit_interaction_ablation_slurm.py", "helper fixture\n")
    write(exp / "interaction_ablation_slurm_worker.sh", "worker fixture\n")
    write(exp / "runs/RUN-180.yaml", "historical record\n")
    write(root / "unrelated/large.bin", "must not stage\n")
    return root, exp, bundle, manifest


def repin(prepared):
    _, exp, bundle, _ = prepared
    write(exp / transport.RUNTIME_NAME, {"manifest_sha256": transport.sha(bundle / "prelaunch_manifest.json")})


def test_integrity_verifies_all_source_inputs_and_environment(prepared):
    _, _, bundle, manifest = prepared
    assert transport.verify_bundle(bundle) == manifest
    (bundle / "source/data/processed/input/data.npz").write_text("changed observations")
    with pytest.raises(ValueError, match="checksum"):
        transport.verify_bundle(bundle)


def test_changed_manifest_cannot_redefine_frozen_inputs_without_repinning(prepared):
    _, _, bundle, manifest = prepared
    manifest["files"].pop("source/data/processed/input/data.npz")
    write(bundle / "prelaunch_manifest.json", manifest)
    with pytest.raises(ValueError, match="manifest changed"):
        transport.verify_bundle(bundle)


def test_undeclared_code_cannot_enter_frozen_import_tree(prepared):
    _, _, bundle, _ = prepared
    write(bundle / "source/scripts/sitecustomize.py", "extra injected source")
    with pytest.raises(ValueError, match="undeclared"):
        transport.verify_bundle(bundle)


@pytest.mark.parametrize("relative", ["../outside.txt", "/tmp/outside.txt", "source/../outside.txt"])
def test_escaping_manifest_paths_are_rejected_even_if_pinned(prepared, relative):
    _, _, bundle, manifest = prepared
    manifest["files"][relative] = "0" * 64
    write(bundle / "prelaunch_manifest.json", manifest)
    repin(prepared)
    with pytest.raises(ValueError, match="Invalid frozen path"):
        transport.verify_bundle(bundle)


def test_symlink_to_valid_interior_payload_is_not_frozen_source(prepared):
    _, _, bundle, _ = prepared
    original = bundle / "source/scripts/model.py"
    other = write(bundle / "source/scripts/other.py", original.read_text())
    original.unlink()
    original.symlink_to(other)
    with pytest.raises(ValueError, match="Invalid frozen path"):
        transport.verify_bundle(bundle)


def test_wrong_bundle_pin_rejected(prepared):
    _, exp, bundle, _ = prepared
    record = json.loads((exp / transport.RUNTIME_NAME).read_text())
    record["runtime_bundle"] = "data/processed/wrong"
    write(exp / transport.RUNTIME_NAME, record)
    with pytest.raises(ValueError, match="different bundle"):
        transport.verify_bundle(bundle)


@pytest.mark.parametrize("execute,fail", [(False, False), (True, False), (True, True)])
def test_submit_is_transactional_and_records_exact_dispatcher(prepared, monkeypatch, execute, fail):
    root, exp, _, _ = prepared
    paths = [root / ".labproj-stageignore", root / "labproject.yaml", root / ".env", exp / "run.sh"]
    original = {p: p.read_bytes() for p in paths}

    def launch(command, **kwargs):
        assert command == ["labproj", "submit", "EXP-065", "--target", "lab-slurm",
                           "--execute" if execute else "--dry-run"]
        assert kwargs["cwd"] == root
        assert not yaml.safe_load((root / "labproject.yaml").read_text())["environment"]["overlay"]["enabled"]
        assert "LAB_SSH_HOST=dutch\n" in (root / ".env").read_text()
        dispatcher = (exp / "run.sh").read_bytes()
        assert dispatcher.splitlines()[1] == b"export ADAPTENV_PROTOCOL=interaction-ablation-slurm"
        record = json.loads((exp / transport.TRANSPORT_NAME).read_text())
        assert record["staged_run_sh_sha256"] == transport.sha(exp / "run.sh")
        filter_text = (root / ".labproj-stageignore").read_text()
        assert "unrelated\n" in filter_text
        assert ".env\n" in filter_text
        assert transport.BUNDLE_REL.as_posix() + "\n" not in filter_text
        write(exp / "runs/RUN-999.yaml", "new run record")
        write(exp / "runs/RUN-999.manifest.yaml", "manifest")
        if fail:
            raise OSError("simulated launch failure")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(transport.subprocess, "run", launch)
    if fail:
        with pytest.raises(OSError, match="simulated"):
            transport.submit(execute)
    else:
        assert transport.submit(execute) == 0
    for path, content in original.items():
        assert path.read_bytes() == content
    assert (exp / "runs/RUN-999.interaction-ablation-slurm.json").is_file()
    assert not (exp / "runs/RUN-999.manifest.interaction-ablation-slurm.json").exists()
    assert not (exp / "runs/RUN-180.interaction-ablation-slurm.json").exists()


def test_missing_dispatcher_handoff_fails_before_mutating_transport(prepared):
    root, exp, _, _ = prepared
    (exp / "run.sh").write_text("#!/bin/bash\nexit 0\n")
    originals = {p: p.read_bytes() for p in (root / ".env", root / "labproject.yaml", root / ".labproj-stageignore")}
    with pytest.raises(ValueError, match="handoff"):
        transport.submit()
    assert all(p.read_bytes() == content for p, content in originals.items())


def fake_linux(monkeypatch, compiler, *, pytest_version="8.4.2"):
    write(compiler, "Zig compiler fixture")
    compiler.chmod(0o755)
    monkeypatch.setattr(transport.platform, "python_version", lambda: "3.11.15")
    monkeypatch.setattr(transport.platform, "system", lambda: "Linux")
    monkeypatch.setattr(transport.platform, "platform", lambda: "Linux-fixture")
    monkeypatch.setattr(transport.importlib.metadata, "distributions", lambda: [
        SimpleNamespace(metadata={"Name": "numpy"}, version="2.4.4"),
        SimpleNamespace(metadata={"Name": "pytest"}, version=pytest_version)])
    monkeypatch.setattr(transport.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(stdout="clang fixture via Zig\n"))


def test_worker_records_environment_before_rejecting_mismatched_pin(prepared, monkeypatch, tmp_path):
    _, _, bundle, _ = prepared
    fake_linux(monkeypatch, tmp_path / "zig", pytest_version="8.4.1")
    log_dir = tmp_path / "logs"
    with pytest.raises(ValueError, match="frozen pins"):
        transport.worker_preflight(bundle, log_dir)
    actual = json.loads((log_dir / "interaction_ablation_linux_environment.json").read_text())
    assert actual["package_mismatches"] == {"pytest": {"expected": "8.4.2", "actual": "8.4.1"}}
    assert actual["verified_payload_count"] == 3


def test_valid_worker_environment_has_reproducibility_record(prepared, monkeypatch, tmp_path):
    _, _, bundle, _ = prepared
    fake_linux(monkeypatch, tmp_path / "zig")
    log_dir = tmp_path / "logs"
    transport.worker_preflight(bundle, log_dir)
    actual = json.loads((log_dir / "interaction_ablation_linux_environment.json").read_text())
    assert actual["package_mismatches"] == {}
    assert (log_dir / "interaction_ablation_linux_numpy_config.txt").is_file()


def test_default_cli_does_not_execute(monkeypatch):
    seen = []
    monkeypatch.setattr(transport.sys, "argv", [str(SCRIPT)])
    monkeypatch.setattr(transport, "submit", lambda execute: seen.append(execute) or 0)
    assert transport.main() == 0
    assert seen == [False]
