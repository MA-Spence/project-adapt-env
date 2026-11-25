"""Joint search engineering checks; no biological adequacy gate is implied."""
from copy import deepcopy
import inspect
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import optuna_joint_readout_exp065 as current
import optuna_readout_exp065 as restored

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
        for ratio in current.HISTORICAL_WARM_RATIOS:
            p = np.array([-.1, np.log(2.), -.7, np.log(.7 + .01 * ratio)])
            loss = restored.make_objective(ratio, a, g, target, bank,
                                            origin, scale, floor, FIT)(p)
            starts.append(dict(ratio=ratio, observation=mode, parameters=p.tolist(),
                protein="TrpB", dataset_id=restored.sweep.IDS["TrpB"],
                additive_prior="empirical_msa_log_profile", objective=loss,
                calibration_origin=origin, calibration_scale=scale, near_floor_raw=floor,
                bounds=restored.parameter_bounds(FIT), offset=origin + scale * p[0],
                scale=scale * np.exp(p[1]), position=p[2], steepness=np.exp(p[3])))
    return dict(a=a, g=g, target=target, banks=banks, origin=origin,
                scale=scale, floor=floor, starts=starts)


def do_fit(problem, output, **changes):
    optimizer = dict(trials=16, n_startup_trials=6, quick_trials=8, sampler="tpe")
    args = dict(mode="none", ca=problem["a"], cg=problem["g"],
        target=problem["target"], bank=None, origin=problem["origin"],
        scale=problem["scale"], floor=problem["floor"], fit=FIT, seed=618,
        warm_starts=problem["starts"], optimizer=optimizer, output=output,
        ratio_bounds=[0., 8.])
    args.update(changes)
    return current.fit_readout(**args)


@pytest.mark.parametrize("mode", current.MODES)
@pytest.mark.parametrize("ratio", [0., .371235, 3.12345678, 8.])
def test_loss_exactly_matches_rolled_back_objective(problem, mode, ratio):
    args = (problem["a"], problem["g"], problem["target"], problem["banks"][mode],
            problem["origin"], problem["scale"], problem["floor"], FIT)
    parameters = problem["starts"][0]["parameters"]
    assert current.make_objective(*args)(ratio, parameters) == restored.make_objective(
        ratio, *args)(parameters)


def test_adaptive_joint_sampling_persists_five_coordinates_and_protects_best_start(problem, tmp_path, capsys):
    optuna = pytest.importorskip("optuna")
    result, trace = do_fit(problem, tmp_path)
    progress = capsys.readouterr().out
    assert progress.count("FIT_PROGRESS") == 1
    assert "FIT_PROGRESS none 16/16 best_loss=" in progress and "best_r=" in progress
    assert list(trace.r.iloc[:6]) == list(current.HISTORICAL_WARM_RATIOS)
    adaptive = trace.iloc[6:]
    assert adaptive.candidate.eq("optuna_tpe").all()
    assert not adaptive.r.isin(current.HISTORICAL_WARM_RATIOS).any()
    assert adaptive.r.between(0., 8.).all()
    assert adaptive[["r", *current.PARAMETER_NAMES]].nunique().gt(1).all()
    historical = [row["objective"]["total"] for row in problem["starts"]
                  if row["observation"] == "none"]
    assert result["warm_start_objective"] == min(historical)
    assert result["objective"]["total"] <= min(historical)
    assert len(result["parameters"]) == 4
    assert result["interaction_concentration"] is None
    assert result["optimized_coordinates"] == ["r", "p0", "p1", "p2", "p3"]
    assert result["evaluations"] == len(trace) == 16
    assert result["optimizer_completed"] and not result["optimizer_converged"]
    database = tmp_path / result["sqlite_history"]
    study = optuna.load_study(study_name="TrpB_none_joint_r", storage=f"sqlite:///{database}")
    assert study.user_attrs["ratio_bounds"] == [0., 8.]
    assert study.user_attrs["interaction_concentration"] is None
    assert study.user_attrs["budget_includes_warm_starts"]
    for trial, row in zip(study.trials, trace.to_dict("records")):
        assert trial.state == optuna.trial.TrialState.COMPLETE
        assert set(trial.params) == {"r", *current.PARAMETER_NAMES}
        assert trial.distributions["r"].low == 0.
        assert trial.distributions["r"].high == 8.
        assert trial.distributions["r"].step is None
        assert trial.value == row["total"]
        for name in trial.params:
            assert trial.params[name] == row[name]
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / result["trial_history"]), trace)
    with pytest.raises(ValueError, match="already exists"):
        do_fit(problem, tmp_path)


def test_serial_tpe_is_deterministic(problem, tmp_path):
    pytest.importorskip("optuna")
    first, trace = do_fit(problem, tmp_path / "first")
    second, repeated = do_fit(problem, tmp_path / "second")
    assert first == second
    pd.testing.assert_frame_equal(trace, repeated)


@pytest.mark.parametrize("field,value,error", [
    ("calibration_scale", .3, "calibration_scale"),
    ("near_floor_raw", .003, "near_floor_raw"),
    ("position", 2., "position"),
    ("parameters", [100., 0., 0., 0.], "outside"),
    ("additive_prior", "synthetic", "profile-conditioned"),
])
def test_stale_warm_start_is_rejected(problem, tmp_path, field, value, error):
    changed = deepcopy(problem["starts"])
    changed[2][field] = value
    with pytest.raises(ValueError, match=error):
        do_fit(problem, tmp_path, warm_starts=changed)
    assert not list(tmp_path.iterdir())


def test_every_old_start_must_replay_and_be_present(problem, tmp_path):
    changed = deepcopy(problem["starts"])
    changed[4]["objective"]["quantile"] += .01
    with pytest.raises(ValueError, match="objective replay mismatch"):
        do_fit(problem, tmp_path, warm_starts=changed)
    for rows in [problem["starts"][1:], problem["starts"] + [problem["starts"][1]]]:
        with pytest.raises(ValueError, match="exactly one"):
            do_fit(problem, tmp_path, warm_starts=rows)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("bounds", [[0., 4.], [-1., 8.], [0., 8., 9.], [0., np.nan]])
def test_changed_r_domain_is_rejected(problem, tmp_path, bounds):
    with pytest.raises(ValueError, match="ratio_bounds"):
        do_fit(problem, tmp_path, ratio_bounds=bounds)
    assert not list(tmp_path.iterdir())


def test_budget_must_cover_six_candidates(problem, tmp_path):
    optimizer = dict(trials=5, n_startup_trials=2, quick_trials=5, sampler="tpe")
    with pytest.raises(ValueError, match="budget"):
        do_fit(problem, tmp_path, optimizer=optimizer)
    assert not list(tmp_path.iterdir())


def test_fit_hook_selects_three_settings_without_confirmation_input(problem, tmp_path):
    pytest.importorskip("optuna")
    config = dict(fit=deepcopy(FIT), observation_modes={"TrpB": list(current.MODES)},
                  ratio_bounds=[0., 8.])
    args = (problem["a"], problem["g"], problem["target"], problem["banks"],
            problem["origin"], problem["scale"], problem["floor"])
    selected = current.fit_settings(*args, config, config, tmp_path, quick=True,
        optimizer=current.DEFAULT_OPTIMIZER, warm_starts=problem["starts"])
    assert len(selected) == 3
    assert [row["observation"] for row in selected] == list(current.MODES)
    assert all(row["evaluations"] == 8 for row in selected)
    assert all(row["protein"] == "TrpB" for row in selected)
    assert all(len(row["parameters"]) == 4 for row in selected)
    assert all(0 <= row["ratio"] <= 8 for row in selected)
    assert all(row["interaction_concentration"] is None for row in selected)
    assert len(list(tmp_path.glob("*.sqlite3"))) == 3
    for function in (current.fit_settings, current.fit_readout, current.make_objective):
        assert not any("confirm" in name or "heldout" in name
                       for name in inspect.signature(function).parameters)


@pytest.mark.parametrize("amendment,error", [
    ({"interaction_ratios": [0., .5, 1., 2., 4., 8.]}, "search grid"),
    ({"ratio_bounds": [0., 4.]}, "ratio_bounds"),
    ({"observation_modes": {"TrpB": ["none"]}}, "all three"),
    ({"fit": dict(FIT, quantiles=3)}, "objective and bounds"),
])
def test_fit_hook_rejects_protocol_drift(problem, tmp_path, amendment, error):
    original = dict(fit=deepcopy(FIT), observation_modes={"TrpB": list(current.MODES)},
                    ratio_bounds=[0., 8.])
    changed = dict(original, **amendment)
    args = (problem["a"], problem["g"], problem["target"], problem["banks"],
            problem["origin"], problem["scale"], problem["floor"])
    with pytest.raises(ValueError, match=error):
        current.fit_settings(*args, changed, original, tmp_path,
            optimizer=current.DEFAULT_OPTIMIZER, warm_starts=problem["starts"])
    assert not list(tmp_path.iterdir())


def test_out_of_range_trial_and_extra_bank_are_rejected(problem, tmp_path):
    args = (problem["a"], problem["g"], problem["target"], problem["banks"]["none"],
            problem["origin"], problem["scale"], problem["floor"], FIT)
    objective = current.make_objective(*args)
    for ratio in (-.1, 8.1, np.inf, np.nan):
        with pytest.raises(ValueError, match="within"):
            objective(ratio, problem["starts"][0]["parameters"])
    config = dict(fit=FIT, observation_modes={"TrpB": list(current.MODES)}, ratio_bounds=[0., 8.])
    banks = dict(problem["banks"], confirmation=np.full((3, 16), np.nan))
    with pytest.raises(ValueError, match="banks"):
        current.fit_settings(problem["a"], problem["g"], problem["target"], banks,
            problem["origin"], problem["scale"], problem["floor"], config, config, tmp_path,
            optimizer=current.DEFAULT_OPTIMIZER, warm_starts=problem["starts"])
