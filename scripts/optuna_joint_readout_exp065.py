#!/usr/bin/env python3
"""Joint continuous-r Optuna search with the unchanged RUN-180 readout loss.

Each observation mode has one five-dimensional study: r in [0, 8] and the four
historical readout coordinates. The six recorded RUN-180 settings are checked
and queued as warm starts. They are candidates, not a fixed-r search grid.
Gamma concentration is not fitted by this protocol.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import optuna_readout_exp065 as restored

DEFAULT_OPTIMIZER = dict(trials=11712, n_startup_trials=20, quick_trials=8, sampler="tpe")
PARAMETER_NAMES = restored.PARAMETER_NAMES
MODES = restored.MODES
HISTORICAL_WARM_RATIOS = (0., .5, 1., 2., 4., 8.)
RATIO_BOUNDS = (0., 8.)


def validate_ratio_bounds(value):
    """Keep the prospectively requested closed continuous domain explicit."""
    bounds = np.asarray(value, dtype=float)
    if bounds.shape != (2,) or not np.array_equal(bounds, RATIO_BOUNDS):
        raise ValueError("Joint calibration requires continuous ratio_bounds [0, 8]")
    return bounds.tolist()


def make_objective(ca, cg, target, bank, origin, scale, floor, fit):
    """Evaluate precisely the restored loss at an arbitrary nonnegative r.

    This delegates to the historical-compatible implementation rather than
    introducing a second definition of quantiles, local SD or floor mass.
    Only calibration arrays and the frozen observation bank are accepted.
    """
    # Validate calibration inputs before constructing a study or writing files.
    restored.make_objective(0., ca, cg, target, bank, origin, scale, floor, fit)

    def parts(ratio, parameters):
        if not np.isfinite(ratio) or not RATIO_BOUNDS[0] <= ratio <= RATIO_BOUNDS[1]:
            raise ValueError("Interaction ratio must be within [0, 8]")
        return restored.make_objective(
            ratio, ca, cg, target, bank, origin, scale, floor, fit)(parameters)

    return parts


def fit_readout(mode, ca, cg, target, bank, origin, scale, floor, fit, seed,
                warm_starts, optimizer, output: Path, quick=False, *, ratio_bounds):
    """Jointly optimize r and readout, preserving every trial and start source."""
    import optuna

    if mode not in MODES or (mode == "none") != (bank is None):
        raise ValueError("Observation mode and frozen bank disagree")
    ratio_bounds = validate_ratio_bounds(ratio_bounds)
    count = restored._optimizer_config(optimizer, quick)
    bounds = restored.parameter_bounds(fit)
    parts = make_objective(ca, cg, target, bank, origin, scale, floor, fit)
    candidates = []
    mode_starts = [row for row in warm_starts if row.get("observation") == mode]
    if (len(mode_starts) != len(HISTORICAL_WARM_RATIOS)
            or {row.get("ratio") for row in mode_starts} != set(HISTORICAL_WARM_RATIOS)):
        raise ValueError(f"Require exactly one RUN-180 warm start for each historical r in {mode}")
    for source_ratio in HISTORICAL_WARM_RATIOS:
        setting = next(row for row in mode_starts if row["ratio"] == source_ratio)
        p = restored._validate_warm_start(setting, mode, bounds, origin, scale, floor)
        loss = parts(source_ratio, p)
        for key, value in loss.items():
            if not np.isclose(value, setting["objective"][key], rtol=1e-9, atol=1e-10):
                raise ValueError(f"RUN-180 objective replay mismatch: {mode}, r={source_ratio:g}, {key}")
        candidates.append(dict(ratio=source_ratio, parameters=p.tolist(),
            candidate=f"RUN180_r{source_ratio:g}_warm_start", source_ratio=source_ratio,
            replay_objective=loss))
    if count < len(candidates):
        raise ValueError("Trial budget must include every historical warm-start candidate")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    stem = f"TrpB_{mode}_joint_r"
    sqlite_path = output / f"{stem}.sqlite3"
    trace_path = output / f"{stem}.csv.gz"
    start_path = output / f"{stem}_warm_starts.json"
    if any(path.exists() for path in (sqlite_path, trace_path, start_path)):
        raise ValueError("Optuna output already exists; historical search records must not be overwritten")
    restored.sweep.dump(start_path, candidates)
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        study_name=stem, direction="minimize", storage=f"sqlite:///{sqlite_path.resolve()}",
        load_if_exists=False,
        sampler=optuna.samplers.TPESampler(seed=int(seed),
            n_startup_trials=int(optimizer["n_startup_trials"])),
    )
    study.set_user_attr("protocol", "RUN180_direct_readout_joint_continuous_r_optuna_tpe")
    study.set_user_attr("warm_start_source_run", "RUN-180")
    study.set_user_attr("reference_coordinate", "seed_reference_anchored")
    study.set_user_attr("optimizer_seed", int(seed))
    study.set_user_attr("budget_includes_warm_starts", True)
    study.set_user_attr("planned_trials", count)
    study.set_user_attr("optimizer_config", optimizer)
    study.set_user_attr("ratio_bounds", ratio_bounds)
    study.set_user_attr("coordinates", ["r", *PARAMETER_NAMES])
    study.set_user_attr("interaction_concentration", None)
    for candidate in candidates:
        study.enqueue_trial(dict(r=candidate["ratio"],
                                 **dict(zip(PARAMETER_NAMES, candidate["parameters"]))),
            user_attrs=dict(candidate=candidate["candidate"],
                            warm_start_source_ratio=candidate["source_ratio"]))
    trace = []

    def objective(trial):
        ratio = trial.suggest_float("r", *ratio_bounds)
        p = [trial.suggest_float(name, lo, hi)
             for name, (lo, hi) in zip(PARAMETER_NAMES, bounds)]
        losses = parts(ratio, p)
        label = trial.user_attrs.get("candidate", "optuna_tpe")
        trial.set_user_attr("candidate", label)
        trial.set_user_attr("loss_components", losses)
        trace.append(dict(evaluation=trial.number, candidate=label, r=ratio,
                          **dict(zip(PARAMETER_NAMES, p)), **losses))
        return losses["total"]

    def report_progress(current_study, finished_trial):
        completed = finished_trial.number + 1
        if (finished_trial.state == optuna.trial.TrialState.COMPLETE
                and (completed % 200 == 0 or completed == count)):
            best = current_study.best_trial
            print("FIT_PROGRESS", mode, f"{completed}/{count}",
                  f"best_loss={best.value:.10g}", f"best_r={best.params['r']:.10g}",
                  flush=True)

    try:
        study.optimize(objective, n_trials=count, n_jobs=1, show_progress_bar=False,
                       callbacks=[report_progress])
    finally:
        pd.DataFrame(trace).to_csv(trace_path, index=False)
    if len(study.trials) != count or any(
            trial.state != optuna.trial.TrialState.COMPLETE for trial in study.trials):
        raise ValueError("Incomplete Optuna calibration; the fixed trial budget did not finish")
    best = study.best_trial
    ratio = float(best.params["r"])
    p = np.array([best.params[name] for name in PARAMETER_NAMES])
    losses = parts(ratio, p)
    previous_best = min(row["replay_objective"]["total"] for row in candidates)
    if losses["total"] > previous_best:
        raise AssertionError("Selected objective is worse than the best re-evaluated RUN-180 setting")
    setting = dict(
        ratio=ratio, observation=mode, parameters=p.tolist(),
        offset=origin + scale * float(p[0]), scale=scale * float(np.exp(p[1])),
        position=float(p[2]), steepness=float(np.exp(p[3])), objective=losses,
        evaluations=len(trace), optimizer_converged=False, optimizer_completed=True,
        optimizer_message="Fixed Optuna trial budget completed; convergence or global optimality not established.",
        selected_source=best.user_attrs["candidate"], bounds=bounds, ratio_bounds=ratio_bounds,
        relative_boundary_distance=[float(min(v-lo, hi-v)/(hi-lo))
                                    for v, (lo, hi) in zip(p, bounds)],
        ratio_relative_boundary_distance=float(min(ratio-ratio_bounds[0],
            ratio_bounds[1]-ratio)/(ratio_bounds[1]-ratio_bounds[0])),
        calibration_origin=origin, calibration_scale=scale, near_floor_raw=floor,
        optimizer="optuna_tpe", optimizer_seed=int(seed), selected_trial=best.number,
        optimized_coordinates=["r", *PARAMETER_NAMES], interaction_concentration=None,
        warm_start_source_run="RUN-180", warm_start_candidates=len(candidates),
        warm_start_objective=previous_best,
        sqlite_history=sqlite_path.name, trial_history=trace_path.name,
    )
    return setting, pd.DataFrame(trace)


def fit_settings(ca, cg, target, banks, origin, scale, floor, config, original_config,
                 output: Path, quick=False, *, optimizer, warm_starts) -> list:
    """Drop-in fitting hook selecting one continuous-r setting per observation mode."""
    restored._optimizer_config(optimizer, quick)
    if config["fit"] != original_config["fit"]:
        raise ValueError("Joint fit must retain RUN-180's readout objective and bounds")
    if "interaction_ratios" in config:
        raise ValueError("Joint calibration must not contain an interaction_ratios search grid")
    ratio_bounds = validate_ratio_bounds(config["ratio_bounds"])
    modes = config["observation_modes"]["TrpB"]
    if list(modes) != list(MODES):
        raise ValueError("Joint calibration must retain all three RUN-180 observation modes in order")
    if set(banks) != set(modes):
        raise ValueError("Frozen observation banks do not match the protocol modes")
    expected = {(mode, ratio) for mode in modes for ratio in HISTORICAL_WARM_RATIOS}
    actual = [(row.get("observation"), row.get("ratio")) for row in warm_starts]
    if set(actual) != expected or len(actual) != len(expected):
        raise ValueError("Warm starts must cover the exact RUN-180 observation/ratio grid")
    selected = []
    for mode_index, mode in enumerate(modes):
        seed = (restored.sweep.seed_for(restored.sweep.IDS["TrpB"]) + config["fit"]["seed"]
                + 1000 + mode_index * 100)
        setting, _ = fit_readout(mode, ca, cg, target, banks[mode], origin, scale, floor,
            config["fit"], seed, warm_starts, optimizer, output, quick,
            ratio_bounds=ratio_bounds)
        setting.update(protein="TrpB", dataset_id=restored.sweep.IDS["TrpB"],
            additive_prior="empirical_msa_log_profile",
            calibration_lineage="RUN180_direct_objective_joint_continuous_r_Optuna_TPE_warm_start")
        selected.append(setting)
        print("FIT", mode, setting["ratio"], setting["objective"]["total"], flush=True)
    return selected
