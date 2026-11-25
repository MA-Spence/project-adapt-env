"""Continuous-r protocol safeguards and replay of the historical fitting inputs."""
import copy
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import optuna_joint_pairwise_calibration_exp065 as protocol

BASELINE = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/outputs"
DATA = ROOT / "data/processed/interaction-strength-exp065-inputs"
PRIOR = ROOT / "data/processed/trpb-profile-additive-inputs"


@pytest.fixture
def config():
    return yaml.safe_load((protocol.EXP / "joint_readout.yaml").read_text())


def test_projection_preserves_recorded_design_without_putting_grid_in_joint_config(config):
    before = copy.deepcopy(config)
    projected = protocol.historical_projection(config)
    baseline = yaml.safe_load((BASELINE / "config.yaml").read_text())
    protocol.rollback.validate_config(projected, baseline)
    assert projected["interaction_ratios"] == [0, 0.5, 1, 2, 4, 8]
    assert projected["optimizer"]["trials"] == 1952
    assert config == before
    assert "interaction_ratios" not in config
    assert config["ratio_bounds"] == [0.0, 8.0]


@pytest.mark.parametrize("change, expected", [
    ({"interaction_ratios": [0, 0.5, 1, 2, 4, 8]}, "must not contain a fitting grid"),
    ({"ratio_bounds": [0.0, 4.0]}, "bounds must be"),
    ({"ratio_bounds": [0.001, 8.0]}, "bounds must be"),
    ({"protocol": "trpb_run180_readout_optuna_v1"}, "must not contain a fitting grid"),
    ({"optimizer": {"sampler": "tpe", "trials": 1952,
                    "n_startup_trials": 20, "quick_trials": 8}}, "optimization budget"),
])
def test_load_baseline_rejects_changed_joint_search(config, change, expected):
    config.update(change)
    with pytest.raises(ValueError, match=expected):
        protocol.load_baseline(BASELINE, config)


@pytest.mark.parametrize("section, key, value", [
    ("model", "interaction_concentration", 1.0),
    ("model", "global_intercept", 0.0),
    ("model", "interaction_backend", "ranked_gnk"),
    ("observation_modes", "TrpB", ["none"]),
    ("observation_modes", "TrpB", ["shared_stop_sensitivity", "independent_difference", "none"]),
    ("fit", "quantiles", 3),
    ("fit", "scale_bounds", [0.1, 10.0]),
    ("sampling", "synthetic_cubes_per_seed", 4),
])
def test_load_baseline_rejects_changes_beyond_joint_optimization(config, section, key, value):
    config[section][key] = value
    with pytest.raises(ValueError, match="RUN-180 design changed"):
        protocol.load_baseline(BASELINE, config)


def test_authenticated_warm_starts_plan_one_joint_study_per_observation(config, monkeypatch, tmp_path):
    settings, manifest = protocol.load_baseline(BASELINE, config)
    assert len(settings) == 18
    assert manifest["settings"] == 18
    planned = []

    def record_study(mode, ca, cg, target, bank, origin, scale, floor, fit, seed,
                     warm_starts, optimizer, output, quick=False, *, ratio_bounds):
        starts = [item for item in warm_starts if item["observation"] == mode]
        planned.append((mode, seed, optimizer["trials"], ratio_bounds, starts))
        return dict(ratio=0.731, observation=mode, objective={"total": 0.1}), None

    monkeypatch.setattr(protocol.fitting, "fit_readout", record_study)
    modes = config["observation_modes"]["TrpB"]
    first = settings[0]
    selected = protocol.fitting.fit_settings(
        None, None, None, {mode: None for mode in modes},
        first["calibration_origin"], first["calibration_scale"], first["near_floor_raw"],
        config, {"fit": config["fit"]}, tmp_path,
        optimizer=config["optimizer"], warm_starts=settings,
    )
    assert [item[0] for item in planned] == ["none", "independent_difference", "shared_stop_sensitivity"]
    assert len({item[1] for item in planned}) == 3
    assert sum(item[2] for item in planned) == 35136
    assert all(item[3] == [0.0, 8.0] and len(item[4]) == 6 for item in planned)
    assert len(selected) == 3
    assert all(item["ratio"] == 0.731 for item in selected)


def test_actual_preflight_replays_all_historical_objectives_without_optimization(config, monkeypatch):
    def unexpected_optimization(*args, **kwargs):
        raise AssertionError("Preflight must not optimize or evaluate confirmation diagnostics")

    monkeypatch.setattr(protocol.fitting, "fit_settings", unexpected_optimization)
    monkeypatch.setattr(protocol.repeat, "run", unexpected_optimization)
    report = protocol.preflight(config, DATA, BASELINE,
                                PRIOR / "additive_matrix.npy", PRIOR / "prior_manifest.json")
    assert report["protocol"] == protocol.PROTOCOL
    assert report["settings"] == 18
    assert report["calibration_latents"] == "allclose"
    assert report["latent_tolerance"] == {"rtol": 1e-13, "atol": 1e-14}
    assert set(report["maximum_absolute_latent_difference"]) == {"a", "g"}
    assert all(0.0 <= value < 1e-12
               for value in report["maximum_absolute_latent_difference"].values())
    assert report["ratio_bounds"] == [0.0, 8.0]
    assert report["historical_ratios_role"] == "Warm starts only; no fitting grid"
    rows = report["objective_replays"]
    assert len(rows) == 18
    assert {(row["observation"], row["ratio"]) for row in rows} == {
        (mode, ratio) for mode in config["observation_modes"]["TrpB"]
        for ratio in protocol.HISTORICAL_RATIOS
    }
    for row in rows:
        assert row["replayed"] == pytest.approx(row["recorded"], rel=1e-12, abs=1e-14)


def shifted_pair_field(monkeypatch, shift):
    original = protocol.repeat.build_model

    def build(*args, **kwargs):
        model = original(*args, **kwargs)

        def latent_components(genotypes):
            additive, pair = model.latent_components(genotypes)
            return additive, pair + shift

        return SimpleNamespace(config=model.config, latent_components=latent_components)

    monkeypatch.setattr(protocol.repeat, "build_model", build)


def test_joint_preflight_accepts_roundoff_and_replays_objective(config, monkeypatch):
    shifted_pair_field(monkeypatch, 1e-15)
    report = protocol.preflight(config, DATA, BASELINE,
                                PRIOR / "additive_matrix.npy", PRIOR / "prior_manifest.json")
    assert report["settings"] == 18
    differences = report["maximum_absolute_latent_difference"]
    assert differences["a"] == 0.0
    assert 0.0 < differences["g"] < 5e-15
    for row in report["objective_replays"]:
        assert row["replayed"] == pytest.approx(row["recorded"], rel=1e-12, abs=1e-14)


def test_joint_preflight_rejects_nonroundoff_pair_drift(config, monkeypatch):
    shifted_pair_field(monkeypatch, 1e-8)
    with pytest.raises(AssertionError):
        protocol.preflight(config, DATA, BASELINE,
                           PRIOR / "additive_matrix.npy", PRIOR / "prior_manifest.json")


def test_original_rollback_preflight_retains_exact_latent_replay(config, monkeypatch):
    shifted_pair_field(monkeypatch, 1e-15)
    with pytest.raises(AssertionError):
        protocol.rollback.preflight(protocol.historical_projection(config), DATA, BASELINE,
                                    PRIOR / "additive_matrix.npy", PRIOR / "prior_manifest.json")
