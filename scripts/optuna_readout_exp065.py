#!/usr/bin/env python3
"""RUN-180's calibration-only readout objective with a serial Optuna TPE search.

The model, fixed interaction-ratio grid, calibration units, observation banks,
101 quantiles, mean local SD and near-floor target retain their RUN-180 meaning.
Only the numerical search changes. Historical settings are checked and evaluated
as candidates, never treated as new scientific results or confirmation targets.
"""
from __future__ import annotations

from numbers import Integral
from pathlib import Path

import numpy as np
import pandas as pd

import interaction_strength_exp065 as sweep

DEFAULT_OPTIMIZER = dict(trials=1952, n_startup_trials=20, quick_trials=8, sampler="tpe")
PARAMETER_NAMES = ("p0", "p1", "p2", "p3")
MODES = ("none", "independent_difference", "shared_stop_sensitivity")


def parameter_bounds(fit: dict) -> list:
    """Four historical optimizer coordinates; positive scales are log coordinates."""
    if (not isinstance(fit.get("quantiles"), Integral)
            or isinstance(fit["quantiles"], bool) or fit["quantiles"] < 2):
        raise ValueError("quantiles must be an integer of at least two")
    bounds = []
    for key in ("offset_bounds", "scale_bounds", "position_bounds", "steepness_bounds"):
        value = np.asarray(fit[key], dtype=float)
        if value.shape != (2,) or not np.isfinite(value).all() or value[0] >= value[1]:
            raise ValueError(f"Invalid {key}")
        if key in ("scale_bounds", "steepness_bounds"):
            if value[0] <= 0:
                raise ValueError(f"{key} must be positive")
            value = np.log(value)
        bounds.append(value.tolist())
    return bounds


def make_objective(ratio, ca, cg, target, bank, origin, scale, floor, fit):
    """Return the exact historical loss components for a readout parameter vector.

    ``target`` is the calibration-cube fitness in the frozen affine units.
    ``ca``/``cg`` are reference-anchored latent calibration fields. Observation
    errors are evaluated in raw assay units with the supplied fixed bank.
    No confirmation observations, topology or inferred coefficients are accepted.
    """
    parameter_bounds(fit)
    ca, cg, target = (np.asarray(value, dtype=float) for value in (ca, cg, target))
    if (ca.ndim != 2 or ca.shape != cg.shape or not ca.size
            or target.ndim != 2 or not target.size or target.shape[-1] < 2
            or ca.shape[-1] != target.shape[-1]):
        raise ValueError("Expected compatible nonempty cube arrays for ca, cg and target")
    if not all(np.isfinite(value).all() for value in (ca, cg, target)):
        raise ValueError("Calibration arrays must be finite")
    if not np.isfinite([ratio, origin, scale, floor]).all() or scale <= 0 or ratio < 0:
        raise ValueError("Invalid ratio or calibration units")
    if bank is not None:
        if not isinstance(bank, (tuple, list)) or len(bank) != 2:
            raise ValueError("Invalid frozen observation bank")
        arrays, breaks = (np.asarray(value) for value in bank)
        if (arrays.ndim != 5 or arrays.shape[1] != 4 or arrays.shape[2] < 1
                or arrays.shape[-2:] != ca.shape or breaks.ndim != 1
                or arrays.shape[0] != len(breaks) + 1
                or not np.isfinite(arrays).all() or not np.isfinite(breaks).all()
                or np.any(np.diff(breaks) <= 0)):
            raise ValueError("Invalid frozen observation bank shape or values")
    probs = np.linspace(0, 1, fit["quantiles"])
    quantiles = np.quantile(target, probs)
    variance = float(target.var())
    mean_sd = float(target.std(axis=-1).mean())
    if variance <= 1e-24 or mean_sd <= 1e-12:
        raise ValueError("Unidentifiable calibration target: variance or local SD is zero")
    floor_affine = (floor - origin) / scale
    mass = float((target <= floor_affine).mean())

    def parts(parameters):
        p = np.asarray(parameters, dtype=float)
        if p.shape != (4,) or not np.isfinite(p).all():
            raise ValueError("Readout parameters must be four finite values")
        # Arithmetic order intentionally matches interaction_strength.fit_readout.
        raw = origin + scale * sweep.predict(p, ca, cg, ratio)
        if bank is not None:
            r1, r2, _, _ = sweep.observed(raw, bank)
            y = ((r1 + r2) / 2 - origin) / scale
        else:
            y = (raw - origin) / scale
        quantile = float(np.mean((np.quantile(y, probs) - quantiles) ** 2) / variance)
        sd = float((y.std(axis=-1).mean() / mean_sd - 1) ** 2)
        floor_loss = float(((y <= floor_affine).mean() - mass) ** 2)
        result = dict(quantile=quantile, mean_sd=sd, floor_cdf=floor_loss,
                      total=quantile + sd + floor_loss)
        if not np.isfinite(list(result.values())).all():
            raise ValueError("Readout objective is not finite")
        return result

    return parts


def _optimizer_config(optimizer, quick):
    if set(optimizer) != set(DEFAULT_OPTIMIZER) or optimizer["sampler"] != "tpe":
        raise ValueError("Expected explicit TPE trials, n_startup_trials and quick_trials")
    for key in ("trials", "n_startup_trials", "quick_trials"):
        if (not isinstance(optimizer[key], Integral) or isinstance(optimizer[key], bool)
                or optimizer[key] < 1):
            raise ValueError(f"{key} must be a positive integer")
    if optimizer["quick_trials"] > optimizer["trials"]:
        raise ValueError("quick_trials cannot exceed trials")
    return int(optimizer["quick_trials"] if quick else optimizer["trials"])


def _validate_warm_start(setting, mode, bounds, origin, scale, floor):
    """Reject stale units, reference coordinates or out-of-domain candidates."""
    if setting.get("observation") != mode:
        raise ValueError("Warm-start observation mode differs")
    if (setting.get("protein") != "TrpB"
            or setting.get("dataset_id") != sweep.IDS["TrpB"]
            or setting.get("additive_prior") != "empirical_msa_log_profile"):
        raise ValueError("Warm start must describe the profile-conditioned TrpB model")
    if not np.isfinite(setting.get("ratio", np.nan)) or setting["ratio"] < 0:
        raise ValueError("Invalid warm-start ratio")
    p = np.asarray(setting["parameters"], dtype=float)
    if p.shape != (4,) or not np.isfinite(p).all():
        raise ValueError("Warm start must have four finite readout coordinates")
    for name, expected in (("calibration_origin", origin), ("calibration_scale", scale),
                           ("near_floor_raw", floor)):
        if not np.isclose(setting[name], expected, rtol=1e-12, atol=1e-14):
            raise ValueError(f"Warm-start {name} differs from the frozen calibration units")
    if (np.shape(setting["bounds"]) != (4, 2)
            or not np.allclose(setting["bounds"], bounds, rtol=0, atol=1e-14)):
        raise ValueError("Warm-start bounds differ from the historical search domain")
    if any(not lo <= value <= hi for value, (lo, hi) in zip(p, bounds)):
        raise ValueError("Warm-start parameters are outside the search bounds")
    for name, expected in (("offset", origin + scale * float(p[0])),
                           ("scale", scale * float(np.exp(p[1]))),
                           ("position", float(p[2])),
                           ("steepness", float(np.exp(p[3])))):
        if not np.isclose(setting[name], expected, rtol=1e-12, atol=1e-14):
            raise ValueError(f"Warm-start {name} does not match its historical coordinates")
    if not np.isfinite([setting["objective"][key]
                        for key in ("quantile", "mean_sd", "floor_cdf", "total")]).all():
        raise ValueError("Warm-start objective must be finite")
    return p


def fit_readout(ratio, mode, ca, cg, target, bank, origin, scale, floor, fit, seed,
                warm_starts, optimizer, output: Path, quick=False):
    """Optimize one fixed r/observation setting and retain its complete trial record."""
    import optuna

    if mode not in MODES or (mode == "none") != (bank is None):
        raise ValueError("Observation mode and frozen bank disagree")
    count = _optimizer_config(optimizer, quick)
    bounds = parameter_bounds(fit)
    parts = make_objective(ratio, ca, cg, target, bank, origin, scale, floor, fit)
    # Same-r RUN-180 candidate first, then its r=0/r=4 endpoint candidates.
    # Cross-r candidates are scored at this r, not by their old objective value.
    candidates = []
    for source_ratio in dict.fromkeys((float(ratio), 0.0, 4.0)):
        matches = [row for row in warm_starts if row.get("observation") == mode
                   and row.get("ratio") == source_ratio]
        if len(matches) != 1:
            raise ValueError(f"Require exactly one RUN-180 warm start for {mode}, r={source_ratio:g}")
        setting = matches[0]
        p = _validate_warm_start(setting, mode, bounds, origin, scale, floor)
        label = ("RUN180_same_ratio" if source_ratio == ratio
                 else f"RUN180_r{source_ratio:g}_endpoint_candidate")
        loss = parts(p)
        if source_ratio == ratio:
            for key, value in loss.items():
                if not np.isclose(value, setting["objective"][key], rtol=1e-9, atol=1e-10):
                    raise ValueError(f"RUN-180 objective replay mismatch: {mode}, r={ratio:g}, {key}")
        candidate = dict(parameters=p.tolist(), candidate=label,
                         source_ratio=source_ratio, replay_objective=loss)
        if not any(np.array_equal(p, row["parameters"]) for row in candidates):
            candidates.append(candidate)
    if count < len(candidates):
        raise ValueError("Trial budget must include every distinct warm-start candidate")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    stem = f"TrpB_{mode}_r{ratio:g}"
    sqlite_path = output / f"{stem}.sqlite3"
    trace_path = output / f"{stem}.csv.gz"
    start_path = output / f"{stem}_warm_starts.json"
    if any(path.exists() for path in (sqlite_path, trace_path, start_path)):
        raise ValueError("Optuna output already exists; historical search records must not be overwritten")
    sweep.dump(start_path, candidates)
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        study_name=stem, direction="minimize", storage=f"sqlite:///{sqlite_path.resolve()}",
        load_if_exists=False,
        sampler=optuna.samplers.TPESampler(seed=int(seed),
            n_startup_trials=int(optimizer["n_startup_trials"])),
    )
    study.set_user_attr("protocol", "RUN180_direct_readout_with_optuna_tpe")
    study.set_user_attr("warm_start_source_run", "RUN-180")
    study.set_user_attr("reference_coordinate", "seed_reference_anchored")
    study.set_user_attr("optimizer_seed", int(seed))
    study.set_user_attr("budget_includes_warm_starts", True)
    study.set_user_attr("planned_trials", count)
    study.set_user_attr("optimizer_config", optimizer)
    for candidate in candidates:
        study.enqueue_trial(dict(zip(PARAMETER_NAMES, candidate["parameters"])),
                            user_attrs=dict(candidate=candidate["candidate"],
                                            warm_start_source_ratio=candidate["source_ratio"]))
    trace = []

    def objective(trial):
        p = [trial.suggest_float(name, lo, hi)
             for name, (lo, hi) in zip(PARAMETER_NAMES, bounds)]
        losses = parts(p)
        label = trial.user_attrs.get("candidate", "optuna_tpe")
        trial.set_user_attr("candidate", label)
        trial.set_user_attr("loss_components", losses)
        trace.append(dict(evaluation=trial.number, candidate=label,
                          **dict(zip(PARAMETER_NAMES, p)), **losses))
        return losses["total"]

    try:
        study.optimize(objective, n_trials=count, n_jobs=1, show_progress_bar=False)
    finally:
        pd.DataFrame(trace).to_csv(trace_path, index=False)
    if len(study.trials) != count or any(
            trial.state != optuna.trial.TrialState.COMPLETE for trial in study.trials):
        raise ValueError("Incomplete Optuna calibration; the fixed trial budget did not finish")
    best = study.best_trial
    p = np.array([best.params[name] for name in PARAMETER_NAMES])
    losses = parts(p)
    if losses["total"] > candidates[0]["replay_objective"]["total"]:
        raise AssertionError("Selected objective is worse than the re-evaluated RUN-180 setting")
    setting = dict(
        ratio=float(ratio), observation=mode, parameters=p.tolist(),
        offset=origin + scale * float(p[0]), scale=scale * float(np.exp(p[1])),
        position=float(p[2]), steepness=float(np.exp(p[3])), objective=losses,
        evaluations=len(trace), optimizer_converged=False, optimizer_completed=True,
        optimizer_message="Fixed Optuna trial budget completed; convergence or global optimality not established.",
        selected_source=best.user_attrs["candidate"], bounds=bounds,
        relative_boundary_distance=[float(min(v-lo, hi-v)/(hi-lo))
                                    for v, (lo, hi) in zip(p, bounds)],
        calibration_origin=origin, calibration_scale=scale, near_floor_raw=floor,
        optimizer="optuna_tpe", optimizer_seed=int(seed), selected_trial=best.number,
        warm_start_source_run="RUN-180", warm_start_candidates=len(candidates),
        warm_start_objective=candidates[0]["replay_objective"]["total"],
        sqlite_history=sqlite_path.name, trial_history=trace_path.name,
    )
    return setting, pd.DataFrame(trace)


def fit_settings(ca, cg, target, banks, origin, scale, floor, config, original_config,
                 output: Path, quick=False, *, optimizer, warm_starts) -> list:
    """Replacement fitting hook for the unchanged full-pairwise EXP-065 pipeline."""
    _optimizer_config(optimizer, quick)
    if config["fit"] != original_config["fit"]:
        raise ValueError("Rollback must retain RUN-180's readout objective and bounds")
    modes = config["observation_modes"]["TrpB"]
    ratios = config["interaction_ratios"]
    if len(set(modes)) != len(modes) or len(set(ratios)) != len(ratios):
        raise ValueError("Duplicate observation modes or ratios")
    if set(banks) != set(modes):
        raise ValueError("Frozen observation banks do not match the protocol modes")
    expected = {(mode, float(ratio)) for mode in modes for ratio in ratios}
    actual = [(row.get("observation"), row.get("ratio")) for row in warm_starts]
    if set(actual) != expected or len(actual) != len(expected):
        raise ValueError("Warm starts must cover the exact RUN-180 observation/ratio grid")
    selected = []
    for mode_index, mode in enumerate(modes):
        for ratio_index, ratio in enumerate(ratios):
            seed = (sweep.seed_for(sweep.IDS["TrpB"]) + config["fit"]["seed"]
                    + 1000 + mode_index * 100 + ratio_index)
            setting, _ = fit_readout(
                ratio, mode, ca, cg, target, banks[mode], origin, scale, floor,
                config["fit"], seed, warm_starts, optimizer, output, quick,
            )
            setting.update(protein="TrpB", dataset_id=sweep.IDS["TrpB"],
                           additive_prior="empirical_msa_log_profile",
                           calibration_lineage="RUN180_direct_objective_Optuna_TPE_warm_start")
            selected.append(setting)
            print("FIT", mode, ratio, setting["objective"]["total"], flush=True)
    return selected
