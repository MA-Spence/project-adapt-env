"""Input and provenance protections for the new architecture-fitting CLI."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import calibrate_architecture as cli


def make_inputs(tmp_path, **changes):
    values = {
        "sequences": np.array([[i // 20, i % 20] for i in range(40)], dtype=np.int64),
        "fitness": np.linspace(0.1, 0.9, 40),
        "additive_matrix": np.stack([np.arange(20), np.arange(20)[::-1]]).astype(float),
    }
    values.update(changes)
    path = tmp_path / "calibration.npz"
    np.savez(path, **values)
    return path


def make_config(tmp_path, **changes):
    values = {"trials": 2, "startup_trials": 1, "calibration_seeds": [101, 102],
              "confirmation_seeds": [201, 202], "bootstrap_replicates": 2,
              "inference": {"n_folds": 2}}
    values.update(changes)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(values))
    return path


def test_validation_is_read_only_and_hashes_input(tmp_path, capsys):
    source = make_inputs(tmp_path)
    config = make_config(tmp_path)
    output = tmp_path / "not-created"
    assert cli.main(["--input", str(source), "--config", str(config),
                     "--output", str(output), "--validate-only"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "structurally_valid"
    assert result["calibration"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert result["calibration"]["alphabet"] == "ACDEFGHIKLMNPQRSTVWY"
    assert not output.exists()


@pytest.mark.parametrize("change, expected", [
    ({"trail_count": 5}, "Unknown calibration fields"),
    ({"inference": {"fold_count": 2}}, "Unknown inference fields"),
    ({"r_bounds": [float("nan"), 8]}, "non-finite"),
])
def test_configuration_typos_and_nonfinite_values_fail(tmp_path, change, expected):
    config_cls, _, inference_cls, _ = cli._api()
    config = make_config(tmp_path, **change)
    with pytest.raises(ValueError, match=expected):
        cli.load_config(config, config_cls, inference_cls)


def test_bundle_rejects_pickles_unknown_fields_and_invalid_profile(tmp_path):
    _, data_cls, _, _ = cli._api()
    with pytest.raises(ValueError, match="Object arrays"):
        cli.load_bundle(make_inputs(tmp_path, fitness=np.array([{}], dtype=object)),
                        data_cls, require_additive=True)
    with pytest.raises(ValueError, match="Unknown NPZ arrays"):
        cli.load_bundle(make_inputs(tmp_path, private_metadata=np.array([1])),
                        data_cls, require_additive=True)
    with pytest.raises(ValueError, match="nonzero within-site variation"):
        cli.load_bundle(make_inputs(tmp_path, additive_matrix=np.ones((2, 20))),
                        data_cls, require_additive=True)
    with pytest.raises(ValueError, match="L x 20"):
        cli.load_bundle(make_inputs(tmp_path, additive_matrix=np.ones((2, 19))),
                        data_cls, require_additive=True)


def test_conflicting_confirmation_profile_fails_before_fit(tmp_path, monkeypatch, capsys):
    source = make_inputs(tmp_path)
    config = make_config(tmp_path)
    confirmation = tmp_path / "confirmation.npz"
    np.savez(confirmation, sequences=np.array([[2, i] for i in range(20)]),
             fitness=np.linspace(0, 1, 20), additive_matrix=np.arange(40).reshape(2, 20))
    config_cls, data_cls, inference_cls, _ = cli._api()

    def fail_if_called(*args, **kwargs):
        pytest.fail("Fitting must not start after input validation fails")

    monkeypatch.setattr(cli, "_api", lambda: (config_cls, data_cls, inference_cls, fail_if_called))
    with pytest.raises(SystemExit) as exit_status:
        cli.main(["--input", str(source), "--config", str(config),
                  "--confirmation", str(confirmation), "--output", str(tmp_path / "output")])
    assert exit_status.value.code == 2
    assert "must match the calibration profile exactly" in capsys.readouterr().err


@pytest.mark.parametrize("case, expected", [
    ("overlap", "must be disjoint"),
    ("constant", "nonzero finite variance"),
    ("noise_without_replicates", ">=2 measured replicates"),
])
def test_structural_validation_rejects_leakage_and_invalid_observation_design(tmp_path, case, expected):
    config_cls, data_cls, inference_cls, _ = cli._api()
    source = make_inputs(tmp_path, **({"fitness": np.ones(40)} if case == "constant" else {}))
    config = make_config(tmp_path, **({"noise_mode": "replicate_residual"}
                                     if case == "noise_without_replicates" else {}))
    resolved, _ = cli.load_config(config, config_cls, inference_cls)
    data, _, _ = cli.load_bundle(source, data_cls, require_additive=True)
    with pytest.raises(ValueError, match=expected):
        cli._validate_panels(data, data if case == "overlap" else None, resolved)


def test_nonempty_output_and_historical_archives_are_never_mutated(tmp_path, monkeypatch):
    output = tmp_path / "old"
    output.mkdir()
    evidence = output / "result.json"
    evidence.write_text("preserve")
    with pytest.raises(ValueError, match="new or empty"):
        cli._guard_output(output)
    assert evidence.read_text() == "preserve"
    monkeypatch.setattr(cli, "ROOT", tmp_path)
    for run in ("RUN-175", "RUN-178", "RUN-180"):
        with pytest.raises(ValueError, match="archived"):
            cli._guard_output(tmp_path / "data/processed/aggregate-thermodynamic-gnk-calibration"
                              / run / "new-output")


def test_dispatch_preserves_config_and_writes_exact_input_provenance(tmp_path, monkeypatch):
    source = make_inputs(tmp_path)
    config = make_config(tmp_path)
    output = tmp_path / "output"
    config_cls, data_cls, inference_cls, _ = cli._api()
    captured = {}

    def fake_fit(data, additive, resolved, *, output_dir, confirmation_data):
        captured.update(data=data, additive=additive, config=resolved, confirmation=confirmation_data)
        output_dir.mkdir()
        (output_dir / "result.json").write_text("{}\n")
        return {}

    monkeypatch.setattr(cli, "_api", lambda: (config_cls, data_cls, inference_cls, fake_fit))
    argv = ["--input", str(source), "--config", str(config), "--output", str(output)]
    assert cli.main(argv) == 0
    assert captured["config"].calibration_seeds == (101, 102)
    assert captured["config"].inference.n_folds == 2
    assert captured["confirmation"] is None
    assert captured["additive"].shape == (2, 20)
    provenance = json.loads((output / "cli_inputs.json").read_text())
    assert provenance["argv"] == argv
    assert provenance["calibration"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert provenance["config"]["sha256"] == hashlib.sha256(config.read_bytes()).hexdigest()
