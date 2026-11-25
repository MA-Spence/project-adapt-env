"""Frozen local-transport provenance; no real RUN is created or executed."""
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/submit_architecture_local.py"
SPEC = importlib.util.spec_from_file_location("architecture_local_transport", SCRIPT)
transport = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(transport)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) if isinstance(value, (dict, list)) else value)
    return path


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    root = tmp_path / "project"
    exp_relative = Path("experiments/EXP-065")
    exp = root / exp_relative
    monkeypatch.setattr(transport, "ROOT", root)
    monkeypatch.setattr(transport, "EXP_REL", exp_relative)
    monkeypatch.setattr(transport, "EXP", exp)
    helper = exp / SCRIPT.name
    monkeypatch.setattr(transport, "__file__", str(helper))
    paths = [write(helper, "transport fixture"),
             write(root / "scripts/calibrate_architecture.py", "frozen CLI"),
             write(root / "external/Adapt-Env/adaptenv/calibration.py", "frozen method"),
             write(root / "external/Adapt-Env/pyproject.toml", "package metadata"),
             write(root / "external/Adapt-Env/LICENSE", "license"),
             write(root / "environment/calibration-py311.lock.txt", "numpy==2.4.4\n")]
    files = {str(p.relative_to(root)): transport.sha(p) for p in paths}
    source = dict(file_count=len(files), files=files)
    write(exp / transport.SOURCE_MANIFEST, source)
    write(root / transport.METHOD_MANIFEST, {"files": files})
    inputs = root / "data/processed" / transport.INPUT_NAME
    calibration = write(inputs / "calibration.npz", "calibration bytes")
    confirmation = write(inputs / "confirmation.npz", "confirmation bytes")
    write(inputs / "manifest.json", {"files": {
        p.name: {"sha256": transport.sha(p), "bytes": p.stat().st_size}
        for p in (calibration, confirmation)}})
    config = dict(trials=3, calibration_seeds=[1, 2], confirmation_seeds=[3, 4])
    write(exp / transport.CONFIG_NAME, config)
    write(exp / "metadata.yaml", yaml.safe_dump({"architecture_followup": {
        "expected_outputs": ["result.json", "study.sqlite3", "trials.json"]}}))
    write(root / "labproject.yaml", yaml.safe_dump({"environment": {"overlay": {"enabled": True}}, "project": {"title": "original"}}))
    write(root / "AGENTS.md", "science record")
    write(root / ".labproj/envs/indirect-py311/bin/python", "interpreter fixture")
    snapshot = root / ".labproj/local-preparations/fixture"
    write(snapshot / "environment/actual_environment.json", {"python": "3.11.15"})
    return root, exp, source, config, snapshot


def snapshot(prepared):
    root, exp, source, config, directory = prepared
    transport.snapshot_sources(directory, source)
    return directory


def outputs(prepared):
    root, exp, _, config, _ = prepared
    output = exp / "outputs/RUN-999"
    for name in ("result.json", "study.sqlite3", "trials.json", "provenance.json", "cli_inputs.json",
                 "config.json", "empirical_architecture.json", "empirical_coefficients.npz",
                 "confirmation_predictions.npz", "inputs.npz", "confirmation_inputs.npz"):
        write(output / name, "placeholder")
    result = dict(method=transport.METHOD, trials=3, completed_trials=2, best_trial=2,
                  selected_inference_converged=True, confirmation_scope="disjoint_genotypes_and_fresh_seeds",
                  calibration_loss=1.5, confirmation_loss=1.8,
                  calibration_seeds=config["calibration_seeds"], confirmation_seeds=config["confirmation_seeds"])
    write(output / "result.json", result)
    write(output / "trials.json", [{"number": 0, "state": "COMPLETE"}, {"number": 1, "state": "FAIL"}, {"number": 2, "state": "COMPLETE"}])
    inputs = root / "data/processed" / transport.INPUT_NAME
    write(output / "cli_inputs.json", dict(config={"sha256": transport.sha(exp / transport.CONFIG_NAME)},
          calibration={"sha256": transport.sha(inputs / "calibration.npz")},
          confirmation={"sha256": transport.sha(inputs / "confirmation.npz")},
          script_sha256=transport.sha(root / "scripts/calibrate_architecture.py")))
    write(output / "provenance.json", {"source_sha256": {
        name: transport.sha(root / "external/Adapt-Env" / name)
        for name in ("adaptenv/calibration.py", "pyproject.toml", "LICENSE")}})
    return output, result


def run_record(prepared, *, returncode=0):
    root, exp, _, _, _ = prepared
    run = dict(id="RUN-999", target="local", experiment="EXP-065", dry_run=False,
               status="completed" if returncode == 0 else "failed",
               paths={"output_dir": str((exp / "outputs/RUN-999").relative_to(root))},
               logs={"stdout": str((exp / "logs/RUN-999/stdout.txt").relative_to(root)),
                     "stderr": str((exp / "logs/RUN-999/stderr.txt").relative_to(root))},
               outputs={"expected": ["historical-output.csv"]},
               execution_result={"returncode": returncode, "stdout": "scientific output", "stderr": "solver messages"})
    return write(exp / "runs/RUN-999.yaml", yaml.safe_dump(run))


def test_snapshot_verifies_source_method_and_input_manifests(prepared):
    root, _, source, _, _ = prepared
    assert transport.verify_frozen() == source
    directory = snapshot(prepared)
    manifest = transport.load(directory / "prelaunch_manifest.json")
    assert "source/docs/model_release/0.7.0rc1/manifest.json" in manifest["files"]
    assert "inputs/trpb-architecture-calibration-inputs/confirmation.npz" in manifest["files"]
    transport.verify_files(directory, manifest["files"])
    (root / "data/processed" / transport.INPUT_NAME / "confirmation.npz").write_text("changed")
    with pytest.raises(ValueError, match="checksum"):
        transport.verify_frozen()


def test_runtime_environment_uses_snapshot_and_scrubs_external_overrides(prepared, monkeypatch):
    root, _, _, _, directory = prepared
    monkeypatch.setenv("ADAPTENV_DATA_ROOT", "/live/override")
    monkeypatch.setenv("PYTHONPATH", "/unexpected/site-packages")
    env = transport.runtime_environment(directory)
    assert env["ADAPTENV_ARCHITECTURE_SNAPSHOT"] == str(directory)
    assert env["ADAPTENV_PROTOCOL"] == "architecture"
    assert env["PYTHON"] == str(root / ".labproj/envs/indirect-py311/bin/python")
    assert "ADAPTENV_DATA_ROOT" not in env and "PYTHONPATH" not in env
    assert env["OMP_NUM_THREADS"] == env["OPENBLAS_NUM_THREADS"] == env["MKL_NUM_THREADS"] == "1"


def test_output_validation_uses_package_relative_provenance_and_allows_failed_trials(prepared):
    directory = snapshot(prepared)
    output, _ = outputs(prepared)
    validation = transport.validate_outputs(output, directory)
    assert validation["passed"], validation
    assert validation["completed_trials"] == 2 and validation["failed_trials"] == 1
    provenance = transport.load(output / "provenance.json")
    provenance["source_sha256"]["adaptenv/calibration.py"] = "changed"
    write(output / "provenance.json", provenance)
    assert not transport.validate_outputs(output, directory)["passed"]


@pytest.mark.parametrize("key,value", [("trials", 2), ("selected_inference_converged", False),
                                      ("method", "other"), ("confirmation_loss", float("nan")),
                                      ("calibration_seeds", [10]),
                                      ("confirmation_scope", "fresh_seeds_on_calibration_genotypes_only")])
def test_partial_or_unconverged_results_never_validate(prepared, key, value):
    directory = snapshot(prepared)
    output, result = outputs(prepared)
    result[key] = value
    write(output / "result.json", result)
    assert not transport.validate_outputs(output, directory)["passed"]


def test_overlay_restoration_preserves_unrelated_changes_and_does_not_repeat(prepared):
    root, _, _, _, directory = prepared
    original = (root / "labproject.yaml").read_bytes()
    transport.disable_overlay(directory)
    current = yaml.safe_load((root / "labproject.yaml").read_text())
    current["project"]["title"] = "concurrent edit"
    (root / "labproject.yaml").write_text(yaml.safe_dump(current))
    transport.restore_overlay(directory)
    current = yaml.safe_load((root / "labproject.yaml").read_text())
    assert current["project"]["title"] == "concurrent edit"
    assert current["environment"]["overlay"]["enabled"] is True
    current["environment"]["overlay"]["enabled"] = False
    (root / "labproject.yaml").write_text(yaml.safe_dump(current))
    transport.restore_overlay(directory)
    assert yaml.safe_load((root / "labproject.yaml").read_text())["environment"]["overlay"]["enabled"] is False
    assert original != (root / "labproject.yaml").read_bytes()


def test_failed_launch_restores_overlay_exactly(prepared, monkeypatch):
    root, _, _, _, directory = prepared
    snapshot(prepared)
    write(directory / "launch.json", {"labproj": "/fake/labproj"})
    original = (root / "labproject.yaml").read_bytes()

    def fail(command, **kwargs):
        assert kwargs["env"]["ADAPTENV_ARCHITECTURE_SNAPSHOT"] == str(directory)
        assert not yaml.safe_load((root / "labproject.yaml").read_text())["environment"]["overlay"]["enabled"]
        raise OSError("fixture launch failure")

    monkeypatch.setattr(transport.subprocess, "Popen", fail)
    with pytest.raises(OSError, match="fixture launch"):
        transport.supervise(directory, execute=True)
    assert (root / "labproject.yaml").read_bytes() == original


@pytest.mark.parametrize("returncode", [0, 1])
def test_collection_preserves_success_failure_and_is_idempotent(prepared, returncode):
    root, _, _, _, _ = prepared
    directory = snapshot(prepared)
    output, _ = outputs(prepared)
    path = run_record(prepared, returncode=returncode)
    if returncode:
        (output / "result.json").unlink()
        write(output / "failure.json", {"stage": "confirmation", "message": "nonconverged"})
    report = transport.collect(path, directory, returncode)
    assert report["status"] == ("verified_complete" if returncode == 0 else "failure_evidence_preserved")
    bundle = root / report["bundle"]
    transport.verify_files(bundle, transport.load(bundle / "collection_manifest.json")["files"])
    assert (bundle / "logs/runner-stderr.txt").read_text() == "solver messages"
    assert yaml.safe_load(path.read_text())["outputs"]["expected"] == ["result.json", "study.sqlite3", "trials.json"]
    assert report["dvc_tracking"] == "pending_manual_tracking"
    assert transport.collect(path, directory, returncode) == report
    if returncode:
        assert (bundle / "outputs/failure.json").is_file()


def test_manifest_paths_must_stay_in_frozen_root(tmp_path):
    with pytest.raises(ValueError, match="escapes"):
        transport.verify_files(tmp_path, {"../elsewhere": "unused"})


def test_zero_exit_with_invalid_artifacts_is_explicit_verification_failure(prepared):
    directory = snapshot(prepared)
    output, result = outputs(prepared)
    result["selected_inference_converged"] = False
    write(output / "result.json", result)
    path = run_record(prepared, returncode=0)
    report = transport.collect(path, directory, 0)
    run = yaml.safe_load(path.read_text())
    assert report["status"] == "failure_evidence_preserved"
    assert run["status"] == "failed" and run["verification_status"] == "failed"
    assert run["failure_reason"]["stage"] == "output_verification"
    assert run["execution_result"]["returncode"] == 0
