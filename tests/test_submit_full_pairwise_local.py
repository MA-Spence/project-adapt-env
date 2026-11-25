"""Transport provenance controls; these do not run or judge a calibration."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / (
    "experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/"
    "submit_full_pairwise_local.py"
)
SPEC = importlib.util.spec_from_file_location("pairwise_local_transport", SCRIPT)
transport = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(transport)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(value, dict):
        value = json.dumps(value)
    path.write_text(value)
    return path


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    root = tmp_path / "project"
    exp = root / "experiments/EXP-065"
    monkeypatch.setattr(transport, "ROOT", root)
    monkeypatch.setattr(transport, "EXP", exp)
    helper = exp / "submit_full_pairwise_local.py"
    monkeypatch.setattr(transport, "__file__", str(helper))
    source_file = write(root / "scripts/model.py", "frozen source\n")
    source = {
        "file_count": 1,
        "files": {"scripts/model.py": transport.sha(source_file)},
    }
    write(exp / "full_pairwise_source.json", source)
    for name, manifest in transport.BUNDLES.items():
        bundle = root / "data/processed" / name
        payload = write(bundle / "payload.bin", "immutable input " + name)
        write(bundle / manifest, {"files": {"payload.bin": transport.sha(payload)}})
    prior = root / "data/processed/trpb-profile-additive-inputs"
    matrix = write(prior / "additive_matrix.npy", "frozen matrix fixture")
    config = {
        "protocol": "full_pairwise_test",
        "additive_prior": {
            "bundle": str(prior.relative_to(root)),
            "matrix_file": matrix.name,
            "matrix_sha256": transport.sha(matrix),
            "metadata_file": "prior_manifest.json",
            "metadata_sha256": transport.sha(prior / "prior_manifest.json"),
        },
        "required_worker_environment": {
            "python": "3.11.15", "packages": {"numpy": "2.4.4"},
        },
    }
    write(exp / "full_pairwise.yaml", yaml.safe_dump(config))
    write(exp / "metadata.yaml", yaml.safe_dump({"inputs": {"declared_files": []}}))
    for relative in (
        "labproject.yaml", "AGENTS.md",
        "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178.dvc",
        "data/processed/interaction-strength-exp065-inputs.dvc",
        "data/processed/trpb-profile-additive-inputs.dvc",
    ):
        write(root / relative, "context\n")
    for name in (helper.name, "full_pairwise.md", "full_pairwise_local_validation.json"):
        write(exp / name, "context\n")
    return root, exp, config, source


def output_fixture(prepared):
    root, exp, config, _ = prepared
    output = root / "outputs/RUN-fixture"
    worker = write(output / "worker_environment.json", config["required_worker_environment"])
    manifest = {
        "mode": "MATCHED_PAIRWISE_REPEAT",
        "public_model": "thermodynamic_pairwise_v3",
        "settings": 18, "seed_audits": 24, "graph_rows": 2163,
        "latent_graph_rows": 144, "cube_summary_rows": 2163,
        "config_sha256": transport.sha(exp / "full_pairwise.yaml"),
        "source_manifest_sha256": transport.sha(exp / "full_pairwise_source.json"),
        "input_manifest_sha256": transport.sha(
            root / "data/processed/interaction-strength-exp065-inputs/manifest.json"
        ),
        "files": {worker.name: transport.sha(worker)},
    }
    write(output / "manifest.json", manifest)
    return output, manifest


def test_frozen_inputs_include_immutable_profile_control(prepared):
    root, exp, _, source = prepared
    assert transport.verify_frozen() == source
    control = root / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178/outputs/payload.bin"
    control.write_text("changed recorded control")
    with pytest.raises(ValueError, match="Frozen checksum mismatch"):
        transport.verify_frozen()


@pytest.mark.parametrize("count", [0, 2, None])
def test_manifest_count_is_consistent_not_hardcoded(prepared, count):
    _, exp, _, source = prepared
    source["file_count"] = count
    write(exp / "full_pairwise_source.json", source)
    with pytest.raises(ValueError, match="file_count"):
        transport.verify_frozen()


def test_snapshot_preserves_nested_control_sources_context_and_hashes(prepared, tmp_path):
    _, _, _, source = prepared
    snapshot = tmp_path / "snapshot"
    write(snapshot / "environment/actual_environment.json", {"python": "3.11.15"})
    transport.snapshot_sources(snapshot, source)
    manifest = json.loads((snapshot / "prelaunch_manifest.json").read_text())
    assert manifest["frozen_source_count"] == 1
    assert "inputs/aggregate-thermodynamic-gnk-calibration/RUN-178/outputs/payload.bin" in manifest["files"]
    assert "source/data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178.dvc" in manifest["files"]
    assert "source/experiments/EXP-065/full_pairwise.md" in manifest["files"]
    transport.verify_files(snapshot, manifest["files"])


@pytest.mark.parametrize("key,value", [
    ("mode", "SMOKE_ONLY"),
    ("public_model", "thermodynamic_gnk_v2"),
    ("settings", 17),
    ("seed_audits", 23),
    ("config_sha256", "changed"),
])
def test_completion_rejects_wrong_model_partial_run_or_provenance(prepared, key, value):
    output, manifest = output_fixture(prepared)
    assert transport.validate_outputs(output)["passed"]
    manifest[key] = value
    write(output / "manifest.json", manifest)
    result = transport.validate_outputs(output)
    assert not result["passed"]
    assert result["errors"]


def test_completion_rejects_changed_actual_environment_even_with_valid_hash(prepared):
    output, manifest = output_fixture(prepared)
    worker = write(output / "worker_environment.json", {"python": "3.13"})
    manifest["files"][worker.name] = transport.sha(worker)
    write(output / "manifest.json", manifest)
    result = transport.validate_outputs(output)
    assert not result["passed"]
    assert "environment" in result["errors"][0]


def test_manifest_paths_cannot_escape_frozen_directory(tmp_path):
    with pytest.raises(ValueError, match="escapes"):
        transport.verify_files(tmp_path, {"../elsewhere": "unused"})


def test_launch_failure_restores_overlay_and_uses_isolated_runtime(prepared, monkeypatch):
    root, _, _, _ = prepared
    project_path = write(root / "labproject.yaml", {
        "environment": {"overlay": {"enabled": True}},
    })
    original = project_path.read_bytes()
    python = write(root / ".labproj/envs/trpb-py31115/bin/python", "interpreter fixture")
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--execute"])
    monkeypatch.setattr(transport.shutil, "which", lambda name: "/fake/bin/labproj")
    monkeypatch.setattr(transport, "capture_environment", lambda *args: None)
    monkeypatch.setattr(transport, "snapshot_sources", lambda *args: None)
    monkeypatch.setenv("ADAPTENV_GRAPH_COMPILER", "/unexpected/compiler")
    monkeypatch.setenv("PYTHONPATH", "/unexpected/packages")

    def fail_launch(command, **kwargs):
        assert command == ["/fake/bin/labproj", "submit", "EXP-065", "--target", "local", "--execute"]
        assert not yaml.safe_load(project_path.read_text())["environment"]["overlay"]["enabled"]
        env = kwargs["env"]
        assert env["PYTHON"] == str(python)
        assert env["ADAPTENV_PROTOCOL"] == "full-pairwise"
        assert env["LABPROJ_QUICK"] == "0"
        assert env["OMP_NUM_THREADS"] == env["OPENBLAS_NUM_THREADS"] == env["MKL_NUM_THREADS"] == "1"
        assert "ADAPTENV_GRAPH_COMPILER" not in env
        assert "PYTHONPATH" not in env
        assert Path(env["DVC_SITE_CACHE_DIR"]).is_relative_to(root / ".labproj/local-preparations")
        assert Path(env["DVC_SITE_CACHE_DIR"]).is_dir()
        raise OSError("simulated subprocess launch failure")

    monkeypatch.setattr(transport.subprocess, "Popen", fail_launch)
    with pytest.raises(OSError, match="simulated"):
        transport.main()
    assert project_path.read_bytes() == original
    assert not transport.run_records()
