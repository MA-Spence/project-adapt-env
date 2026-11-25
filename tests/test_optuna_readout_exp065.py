"""Engineering checks for the restored calibration protocol, not adequacy gates."""
from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import optuna_readout_exp065 as current

FIT = dict(seed=660501, maxiter=60, popsize=8, quantiles=101,
           offset_bounds=[-2., 2.], scale_bounds=[.001, 10.],
           position_bounds=[-14., 8.], steepness_bounds=[.05, 8.])


@pytest.fixture
def problem():
    rng = np.random.default_rng(194)
    a, g = rng.normal(size=(2, 3, 16))
    target = rng.normal(.2, .4, size=(2, 16))
    errors = rng.normal(0, .02, size=(4, 4, 2, 3, 16))
    banks = dict(none=None, independent_difference=(errors, [-.2, .1, .4]),
                 shared_stop_sensitivity=(errors * 2, [-.2, .1, .4]))
    origin, scale, floor = -.03, .2, .001
    starts = []
    for mode, bank in banks.items():
        for ratio in (0., .5, 4.):
            p = np.array([-.1, np.log(2.), -.7, np.log(.7 + .01 * ratio)])
            loss = current.make_objective(ratio, a, g, target, bank,
                                          origin, scale, floor, FIT)(p)
            starts.append(dict(ratio=ratio, observation=mode, parameters=p.tolist(),
                protein="TrpB", dataset_id=current.sweep.IDS["TrpB"],
                additive_prior="empirical_msa_log_profile", objective=loss,
                calibration_origin=origin, calibration_scale=scale, near_floor_raw=floor,
                bounds=current.parameter_bounds(FIT), offset=origin + scale * p[0],
                scale=scale * np.exp(p[1]), position=p[2], steepness=np.exp(p[3]),
                optimizer_converged=False, optimizer_message="fixture"))
    return dict(a=a, g=g, target=target, banks=banks, origin=origin,
                scale=scale, floor=floor, starts=starts)


def do_fit(problem, output, **changes):
    optimizer = dict(trials=10, n_startup_trials=3, quick_trials=8, sampler="tpe")
    args = dict(ratio=.5, mode="none", ca=problem["a"], cg=problem["g"],
                target=problem["target"], bank=None, origin=problem["origin"],
                scale=problem["scale"], floor=problem["floor"], fit=FIT, seed=618,
                warm_starts=problem["starts"], optimizer=optimizer, output=output)
    args.update(changes)
    return current.fit_readout(**args)


@pytest.mark.parametrize("mode", current.MODES)
@pytest.mark.parametrize("ratio", [0., .5, 4.])
def test_loss_matches_historical_function_including_frozen_noise(problem, monkeypatch, mode, ratio):
    p = next(row for row in problem["starts"] if row["ratio"] == ratio
             and row["observation"] == mode)["parameters"]
    expected = current.make_objective(ratio, problem["a"], problem["g"], problem["target"],
        problem["banks"][mode], problem["origin"], problem["scale"], problem["floor"], FIT)(p)
    # Exercise the unchanged historical objective without running its DE search.
    def one_evaluation(objective, *args, **kwargs):
        objective(p)
        return SimpleNamespace(success=False, message="one engineering evaluation")
    monkeypatch.setattr(current.sweep, "differential_evolution", one_evaluation)
    anchors = {r: dict(parameters=p, objective=expected, optimizer_converged=False,
                       optimizer_message="engineering fixture") for r in (0, 4)}
    old, trace = current.sweep.fit_readout(ratio, mode, problem["a"], problem["g"],
        problem["target"], problem["banks"][mode], problem["origin"], problem["scale"],
        problem["floor"], FIT, 618, anchors)
    assert old["objective"] == expected
    assert trace.total.eq(expected["total"]).all()


def test_all_18_recorded_run180_objectives_replay(tmp_path):
    archive = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/outputs"
    inputs = ROOT / "data/processed/interaction-strength-exp065-inputs"
    if not (archive / "selected_parameters.json").is_file():
        pytest.skip("RUN-180 DVC payload is not present; restore for scientific replay")
    settings = json.loads((archive / "selected_parameters.json").read_text())
    assert len(settings) == 18
    config = yaml.safe_load((archive / "config.yaml").read_text())
    original = yaml.safe_load((inputs / "config.yaml").read_text())
    sub = inputs / "TrpB"
    noise, _ = current.sweep.calibrate_noise(
        pd.read_csv(sub / "donor_assignments.csv.gz"),
        set(json.loads((sub / "excluded_noise_donors.json").read_text())), original, tmp_path)
    with np.load(sub / "fit_observation_uniforms.npz") as file:
        uniforms = file["uniforms"]
    banks = {"none": None}
    for mode in current.MODES[1:]:
        banks[mode] = current.sweep.observation_bank(noise, uniforms,
                                                   mode == "shared_stop_sensitivity")
    with np.load(archive / "latent_calibration.npz") as file:
        a, g = file["a"], file["g"]
    index = pd.read_csv(sub / "empirical_panel_index.csv")
    with np.load(sub / "empirical_panel.npz") as file:
        calibration = file["fitness"][index.split.eq("calibration")]
    for row in settings:
        origin, scale, floor = (row[key] for key in
                               ("calibration_origin", "calibration_scale", "near_floor_raw"))
        actual = current.make_objective(row["ratio"], a, g, (calibration-origin)/scale,
            banks[row["observation"]], origin, scale, floor, config["fit"])(row["parameters"])
        np.testing.assert_allclose(list(actual.values()), list(row["objective"].values()),
                                   rtol=1e-12, atol=1e-12)
        current._validate_warm_start(row, row["observation"],
            current.parameter_bounds(config["fit"]), origin, scale, floor)


def test_seeded_tpe_is_reproducible_cannot_degrade_and_persists_history(problem, tmp_path):
    optuna = pytest.importorskip("optuna")
    first, trace = do_fit(problem, tmp_path / "first")
    second, repeated = do_fit(problem, tmp_path / "second")
    pd.testing.assert_frame_equal(trace, repeated)
    assert first == second
    assert len(trace) == first["evaluations"] == 10
    assert trace.iloc[0].candidate == "RUN180_same_ratio"
    assert trace.iloc[1].candidate == "RUN180_r0_endpoint_candidate"
    assert trace.iloc[2].candidate == "RUN180_r4_endpoint_candidate"
    assert first["objective"]["total"] <= first["warm_start_objective"]
    assert first["optimizer_completed"] and not first["optimizer_converged"]
    db = tmp_path / "first" / first["sqlite_history"]
    study = optuna.load_study(study_name="TrpB_none_r0.5", storage=f"sqlite:///{db}")
    assert len(study.trials) == 10
    assert all(trial.state == optuna.trial.TrialState.COMPLETE for trial in study.trials)
    np.testing.assert_array_equal([trial.value for trial in study.trials], trace.total)
    assert study.user_attrs["budget_includes_warm_starts"]
    csv = pd.read_csv(tmp_path / "first" / first["trial_history"])
    pd.testing.assert_frame_equal(csv, trace)
    with pytest.raises(ValueError, match="already exists"):
        do_fit(problem, tmp_path / "first")


@pytest.mark.parametrize("field,new_value,error", [
    ("calibration_scale", .3, "calibration_scale"),
    ("near_floor_raw", .003, "near_floor_raw"),
    ("position", 2., "position"),
    ("parameters", [100., 0., 0., 0.], "outside"),
    ("bounds", [[-2., 2.]] * 4, "bounds"),
    ("additive_prior", "synthetic", "profile-conditioned"),
])
def test_stale_or_changed_warm_start_rejected(problem, tmp_path, field, new_value, error):
    changed = deepcopy(problem["starts"])
    row = next(row for row in changed if row["ratio"] == .5 and row["observation"] == "none")
    row[field] = new_value
    with pytest.raises(ValueError, match=error):
        do_fit(problem, tmp_path, warm_starts=changed)
    assert not list(tmp_path.iterdir())


def test_changed_objective_replay_rejected_before_search(problem, tmp_path):
    changed = deepcopy(problem["starts"])
    changed[1]["objective"]["quantile"] += .01
    with pytest.raises(ValueError, match="objective replay mismatch"):
        do_fit(problem, tmp_path, warm_starts=changed)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("change,error", [
    ({"target": np.ones((2, 16))}, "Unidentifiable"),
    ({"ca": np.full((3, 16), np.nan)}, "finite"),
    ({"scale": 0}, "units"),
    ({"mode": "none", "bank": "bad"}, "bank disagree"),
    ({"optimizer": dict(trials=2, n_startup_trials=1, quick_trials=2, sampler="tpe")}, "budget"),
    ({"optimizer": dict(trials=10, n_startup_trials=3, quick_trials=8, sampler="random")}, "explicit TPE"),
])
def test_invalid_inputs_rejected(problem, tmp_path, change, error):
    with pytest.raises(ValueError, match=error):
        do_fit(problem, tmp_path, **change)
    assert not list(tmp_path.iterdir())


def test_missing_or_duplicate_warm_start_is_rejected(problem, tmp_path):
    for rows in [problem["starts"][1:], problem["starts"] + [problem["starts"][1]]]:
        with pytest.raises(ValueError, match="exactly one"):
            do_fit(problem, tmp_path, warm_starts=rows)


def test_fit_hook_rejects_changed_objective_and_wrong_grid(problem, tmp_path):
    config = dict(fit=deepcopy(FIT), observation_modes={"TrpB": list(current.MODES)},
                  interaction_ratios=[0., .5, 4.])
    args = (problem["a"], problem["g"], problem["target"], problem["banks"],
            problem["origin"], problem["scale"], problem["floor"])
    changed = deepcopy(config)
    changed["fit"]["quantiles"] = 3
    with pytest.raises(ValueError, match="objective and bounds"):
        current.fit_settings(*args, changed, config, tmp_path,
                             optimizer=current.DEFAULT_OPTIMIZER, warm_starts=problem["starts"])
    with pytest.raises(ValueError, match="exact RUN-180"):
        current.fit_settings(*args, config, config, tmp_path,
                             optimizer=current.DEFAULT_OPTIMIZER, warm_starts=problem["starts"][:-1])
