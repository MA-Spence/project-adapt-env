"""Frozen-output analysis for ANA-071. No fitting and no live AdaptEnv imports."""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit

ANALYSIS_NAME = "2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy"
COLORS = {"calibration": "#2478A6", "confirmation": "#D06837", "empirical": "#252525", "bootstrap": "#8962A9"}
LABELS = [
    "log(1 + additive latent variance)", "log(1 + pair latent variance)", "Pair share of latent variance",
    "Additive |effect|: median / RMS", "Additive |effect|: 90th / RMS", "Additive |effect|: 99th / RMS",
    "Pair |effect|: median / RMS", "Pair |effect|: 90th / RMS", "Pair |effect|: 99th / RMS",
    "Additive: top 10% energy share", "Pair: top 10% energy share",
    "Fitness: 10th percentile", "Fitness: median", "Fitness: 90th percentile",
    "Inferred sigmoid lower limit", "Inferred sigmoid upper limit",
    "Inferred occupancy below 0.05", "Inferred occupancy above 0.95",
]
PARAMETERS = ["r", "concentration", "steepness", "position", "lower", "span"]
PARAMETER_LABELS = {
    "r": "Interaction strength r", "concentration": "Effect concentration alpha", "steepness": "Sigmoid steepness",
    "position": "Sigmoid position", "lower": "Output lower bound (standardized)", "span": "Output span (standardized)",
}


def find_root(start: Path | None = None) -> Path:
    here = (start or Path.cwd()).resolve()
    for path in (here, *here.parents):
        if (path / "AIM.md").exists() and (path / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-184/outputs/result.json").exists():
            return path
    raise FileNotFoundError("Run from this repository after restoring the RUN-184 DVC bundle")


def read_json(path):
    return json.loads(Path(path).read_text())


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def r2(actual, predicted):
    actual, predicted = np.asarray(actual), np.asarray(predicted)
    finite = np.isfinite(actual) & np.isfinite(predicted)
    actual, predicted = actual[finite], predicted[finite]
    return float(1 - np.square(actual - predicted).sum() / np.square(actual - actual.mean()).sum())


def save_figure(fig, directory, name):
    for extension in ("png", "pdf"):
        fig.savefig(directory / f"{name}.{extension}", dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def analyze(root: Path | None = None):
    root = root or find_root()
    outputs = root / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-184/outputs"
    analysis = root / "analyses" / ANALYSIS_NAME
    table_dir, figure_dir = analysis / "tables", analysis / "figures"
    table_dir.mkdir(exist_ok=True)
    figure_dir.mkdir(exist_ok=True)
    result, empirical, trials, config = [read_json(outputs / name) for name in
                                        ("result.json", "empirical_architecture.json", "trials.json", "config.json")]
    best = next(trial for trial in trials if trial["number"] == result["best_trial"])
    names = empirical["summary_names"]
    scales = np.asarray(empirical["fixed_scales"])
    bootstraps = np.asarray(empirical["bootstrap_summaries"])
    np.testing.assert_allclose(scales, np.sqrt(bootstraps.var(axis=0, ddof=1) + config["scale_floor"] ** 2), rtol=1e-12)
    targets = {"calibration": np.asarray(empirical["summary"]), "confirmation": np.asarray(result["confirmation_architecture"]["target"])}
    simulated = {"calibration": np.asarray(best["diagnostics"]["summaries_by_seed"]),
                 "confirmation": np.asarray(result["confirmation_architecture"]["simulated"])}
    residuals, contributions, rows = {}, {}, []
    origin, scale = result["phenotype_origin"], result["phenotype_scale"]
    for panel in targets:
        residuals[panel] = (simulated[panel] - targets[panel]) / scales
        contributions[panel] = np.square(residuals[panel].mean(axis=0)) / len(names)
        np.testing.assert_allclose(contributions[panel].sum(), result[f"{panel}_loss"], rtol=1e-12)
        for j, name in enumerate(names):
            for seed, value, residual in zip(result[f"{panel}_seeds"], simulated[panel][:, j], residuals[panel][:, j]):
                is_fitness = name.startswith("fitness_") or name in ("fitted_lower", "fitted_upper")
                rows.append(dict(panel=panel, seed=seed, summary=name, label=LABELS[j], empirical=targets[panel][j],
                                 synthetic=value, fixed_scale=scales[j], standardized_residual=residual,
                                 ensemble_loss_contribution=contributions[panel][j],
                                 empirical_native=(origin + scale * targets[panel][j]) if is_fitness else targets[panel][j],
                                 synthetic_native=(origin + scale * value) if is_fitness else value,
                                 units="native fitness" if is_fitness else "summary units"))
    architecture = pd.DataFrame(rows)
    architecture.to_csv(table_dir / "architecture_summary_by_seed.csv", index=False)
    aggregate = architecture.groupby(["panel", "summary", "label"], sort=False).agg(
        empirical=("empirical", "first"), synthetic_mean=("synthetic", "mean"), synthetic_min=("synthetic", "min"),
        synthetic_max=("synthetic", "max"), empirical_native=("empirical_native", "first"),
        synthetic_native_mean=("synthetic_native", "mean"), fixed_scale=("fixed_scale", "first"),
        standardized_mean_residual=("standardized_residual", "mean"), loss_contribution=("ensemble_loss_contribution", "first")).reset_index()
    aggregate.to_csv(table_dir / "architecture_summary_aggregate.csv", index=False)
    bootstrap_rows = []
    for b, (values, fit) in enumerate(zip(bootstraps, empirical["bootstrap_inference"])):
        for j, name in enumerate(names):
            bootstrap_rows.append(dict(bootstrap=b, summary=name, value=values[j], full_data_target=targets["calibration"][j],
                                       fixed_scale=scales[j], offset_from_full_target_in_fixed_scales=(values[j] - targets["calibration"][j]) / scales[j],
                                       selected_penalty=fit["selected_penalty"], selected_cv_mean_r2=fit["convergence"]["selected_cv_mean_r2"]))
    pd.DataFrame(bootstrap_rows).to_csv(table_dir / "architecture_bootstrap.csv", index=False)

    # Read archived empirical coefficients once, and predict the disjoint panel.
    # These coefficients were fitted to empirical calibration observations; this
    # is a different model from the randomized, summary-calibrated generator.
    with np.load(outputs / "empirical_coefficients.npz", allow_pickle=False) as data:
        coefficients = {key: data[key] for key in data.files}
    with np.load(outputs / "confirmation_inputs.npz", allow_pickle=False) as data:
        sequences, observed = data["sequences"], data["fitness"]
    latent = np.full(len(sequences), float(coefficients["intercept"]))
    for site in range(sequences.shape[1]):
        latent += coefficients["additive"][site, sequences[:, site]]
    for table, (i, j) in zip(coefficients["pairwise"], coefficients["pairs"]):
        latent += table[sequences[:, i], sequences[:, j]]
    prediction = origin + scale * (float(coefficients["lower"]) +
                                  (float(coefficients["upper"]) - float(coefficients["lower"])) * expit(latent))
    heldout = dict(model="empirical calibration coefficient model; not synthetic generator", genotypes=len(sequences),
                   individual_replicate_r2=r2(observed, np.broadcast_to(prediction[:, None], observed.shape)),
                   genotype_mean_r2=r2(observed.mean(axis=1), prediction),
                   genotype_mean_rmse=float(np.sqrt(np.square(observed.mean(axis=1) - prediction).mean())),
                   fitting_scope="Empirical calibration coefficients; no refitting or confirmation selection")
    pd.DataFrame([heldout]).to_csv(table_dir / "architecture_empirical_model_heldout_prediction.csv", index=False)

    fit_rows, cv_rows = [], []
    fits = [("calibration", "empirical", -1, empirical["inference"]),
            ("confirmation", "empirical", -1, result["confirmation_architecture"]["empirical_inference"])]
    fits.extend(("calibration", "synthetic", seed, fit) for seed, fit in zip(result["calibration_seeds"], best["diagnostics"]["auxiliary_inference"]))
    fits.extend(("confirmation", "synthetic", seed, fit) for seed, fit in zip(result["confirmation_seeds"], result["confirmation_architecture"]["synthetic_inference"]))
    for panel, kind, seed, fit in fits:
        convergence = fit["convergence"]
        fit_rows.append(dict(panel=panel, kind=kind, seed=seed, selected_penalty=fit["selected_penalty"],
                             selected_cv_mean_r2=convergence["selected_cv_mean_r2"], converged=fit["converged"],
                             iterations=convergence["iterations"], lower=fit["lower"], upper=fit["upper"]))
        for row in fit["cv_records"]:
            cv_rows.append(dict(panel=panel, kind=kind, seed=seed, selected_penalty=fit["selected_penalty"], **row))
    fits_frame = pd.DataFrame(fit_rows)
    fits_frame.to_csv(table_dir / "architecture_inference_fits.csv", index=False)
    pd.DataFrame(cv_rows).to_csv(table_dir / "architecture_cv_folds.csv", index=False)

    trial_rows = []
    for trial in trials:
        row = dict(number=trial["number"], state=trial["state"], loss=trial["value"], **trial["parameters"])
        diagnostic = trial["diagnostics"]
        row["elapsed_seconds"] = diagnostic.get("elapsed_seconds")
        if "failure" in diagnostic:
            failure_fit = diagnostic["failure"].get("inference", {})
            failure = failure_fit.get("convergence", {})
            row.update(failure_message=diagnostic["failure"].get("message"), optimizer_message=failure.get("message"),
                       full_fit_converged=failure.get("converged"), selected_cv_converged=failure.get("selected_cv_converged"),
                       full_fit_iterations=failure.get("iterations"), selected_penalty=failure_fit.get("selected_penalty"))
        trial_rows.append(row)
    search = pd.DataFrame(trial_rows)
    search.to_csv(table_dir / "architecture_search_trials.csv", index=False)
    parameter_rows = []
    for parameter in PARAMETERS:
        low, high = config[f"{parameter}_bounds"]
        value = result["best_parameters_standardized"][parameter]
        log_search = parameter in ("concentration", "steepness", "span")
        position = (np.log(value / low) / np.log(high / low)) if log_search else ((value - low) / (high - low))
        parameter_rows.append(dict(parameter=parameter, selected_standardized=value, selected_native=result["best_parameters"][parameter],
                                   lower_search_bound=low, upper_search_bound=high, search_is_log=log_search,
                                   position_in_search_coordinate=position,
                                   nearest_boundary_fraction_linear=(min(value-low, high-value)/(high-low)),
                                   near_best_trial_count=len(result["near_best_region"]["trials"]),
                                   attempted_min=float(search[parameter].min()), attempted_max=float(search[parameter].max()),
                                   completed_min=float(search.loc[search.state == "COMPLETE", parameter].min()),
                                   completed_max=float(search.loc[search.state == "COMPLETE", parameter].max())))
    parameters = pd.DataFrame(parameter_rows)
    parameters.to_csv(table_dir / "architecture_selected_parameters.csv", index=False)

    summary = {
        "source": "Frozen RUN-184; no auxiliary or generator refitting",
        "losses_recomputed": {panel: float(contributions[panel].sum()) for panel in targets},
        "loss_interpretation": "Mean squared fixed-scale summary distance; not R2, likelihood, or a calibrated goodness-of-fit statistic",
        "attempts": len(trials), "successful": int((search.state == "COMPLETE").sum()), "failed": int((search.state == "FAIL").sum()),
        "failure_optimizer_messages": dict(Counter(search.loc[search.state == "FAIL", "optimizer_message"])),
        "best_trial": result["best_trial"], "near_best_trials": result["near_best_region"]["trials"],
        "selected_parameters": parameter_rows, "empirical_model_heldout_prediction": heldout,
        "selected_inference_fits": fit_rows,
        "bootstrap": {
            "count": len(bootstraps), "all_converged": all(fit["converged"] for fit in empirical["bootstrap_inference"]),
            "penalty_counts": dict(Counter(str(fit["selected_penalty"]) for fit in empirical["bootstrap_inference"])),
            "additive_log1p_variance_target": float(targets["calibration"][0]),
            "additive_log1p_variance_min_median_max": np.quantile(bootstraps[:, 0], [0, .5, 1]).tolist(),
            "pair_share_target": float(targets["calibration"][2]),
            "pair_share_min_median_max": np.quantile(bootstraps[:, 2], [0, .5, 1]).tolist(),
            "largest_mean_offset_in_fixed_scales": float(np.max(np.abs((bootstraps.mean(axis=0) - targets["calibration"]) / scales))),
            "interpretation": "Bootstrap architectures shift substantially from the original fit despite convergence. Causes are not identified; spread around shifted fits is not assurance of target uncertainty coverage.",
        },
        "gamma": {
            "alpha": result["best_parameters"]["concentration"],
            "raw_variance_mixture_cv": float(1 / np.sqrt(result["best_parameters"]["concentration"])),
            "raw_coefficient_excess_kurtosis": float(3 / result["best_parameters"]["concentration"]),
            "latent_generator_pair_fraction": float(result["best_parameters"]["r"] ** 2 / (2 + result["best_parameters"]["r"] ** 2)),
            "interpretation": "Near-Gaussian raw coefficient distribution. Formulas precede centering and realized RMS normalization. This run does not establish a benefit of finite-alpha concentration.",
        },
        "architecture_by_panel": {},
    }
    for panel in targets:
        order = np.argsort(contributions[panel])[::-1]
        summary["architecture_by_panel"][panel] = dict(
            empirical=targets[panel].tolist(), synthetic_mean=simulated[panel].mean(axis=0).tolist(),
            standardized_mean_residual=residuals[panel].mean(axis=0).tolist(),
            count_mean_residual_magnitude_above_one_fixed_scale=int((np.abs(residuals[panel].mean(axis=0)) > 1).sum()),
            contribution_ranking=[dict(summary=names[j], contribution=float(contributions[panel][j]),
                                       fraction_of_loss=float(contributions[panel][j] / contributions[panel].sum())) for j in order])
    summary["summary_names"] = names
    save_json(table_dir / "architecture_findings.json", summary)
    return dict(root=root, analysis=analysis, figures=figure_dir, tables=table_dir, result=result, empirical=empirical, config=config,
                targets=targets, simulated=simulated, residuals=residuals, contributions=contributions, scales=scales,
                bootstrap=bootstraps, search=search, parameters=parameters, fits=fits_frame, heldout=heldout,
                observed=observed, empirical_prediction=prediction, summary=summary)


def set_style():
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 180})


def plot_architecture(data):
    fig, axes = plt.subplots(1, 2, figsize=(16, 10.5), gridspec_kw={"width_ratios": [1.5, 1]}, sharey=True)
    y = np.arange(len(LABELS))
    for panel, offset in (("calibration", -.15), ("confirmation", .15)):
        residual = data["residuals"][panel]
        for k in range(residual.shape[0]):
            axes[0].scatter(residual[k], y + offset + (k - 1.5) * .035, color=COLORS[panel], s=16, alpha=.65)
        axes[0].scatter(residual.mean(axis=0), y + offset, color=COLORS[panel], marker="D", s=28,
                        edgecolor="white", linewidth=.4, label=panel.capitalize())
        axes[1].barh(y + offset, data["contributions"][panel], height=.28, color=COLORS[panel], alpha=.9)
    axes[0].axvline(0, color="#222222", linewidth=1)
    for x in (-1, 1):
        axes[0].axvline(x, color="#aaaaaa", linewidth=.7, linestyle=":")
    axes[0].set_yticks(y, LABELS)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("(Synthetic − empirical) / frozen summary scale")
    axes[0].set_title("A  Most coefficient-shape and fitness targets remain mismatched", loc="left", pad=15)
    axes[0].legend(loc="lower left", frameon=False)
    axes[1].set_xlabel("Contribution to total loss (sums to 24.03 / 27.35)")
    axes[1].set_title("B  Fitness upper tail is the largest loss component", loc="left", pad=15)
    for ax in axes:
        for division in (2.5, 8.5, 10.5, 13.5, 15.5):
            ax.axhline(division, color="#dddddd", linewidth=.7)
        ax.grid(axis="x", color="#eeeeee", zorder=0)
    fig.suptitle("RUN-184: agreement on inferred architecture is incomplete", fontsize=17, x=.04, ha="left", y=.99)
    fig.text(.04, .025,
             "Dots show all four generator seeds; diamonds show their means. No confidence intervals are implied.\n"
             "Each panel uses its own empirical target and the same calibration-only weights. Fixed-scale units are not validated standard errors;\n"
             "loss is squared ensemble-mean discrepancy / 18. Uniform-product latent variance is not observed fitness epistasis; limits / occupancy are auxiliary-inferred.",
             ha="left", va="bottom", fontsize=10)
    fig.subplots_adjust(left=.265, right=.98, top=.91, bottom=.14, wspace=.12)
    save_figure(fig, data["figures"], "03_architecture")


def plot_inference(data):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
    observed, prediction = data["observed"].mean(axis=1), data["empirical_prediction"]
    ax = axes[0, 0]
    ax.hexbin(observed, prediction, gridsize=45, mincnt=1, bins="log", cmap="Blues", linewidths=0)
    limits = (min(observed.min(), prediction.min()), max(observed.max(), prediction.max()))
    ax.plot(limits, limits, color="#777777", linestyle=":", linewidth=1)
    ax.set(xlabel="Observed confirmation fitness (genotype mean)", ylabel="Prediction from empirical calibration coefficients")
    ax.set_title("A  Empirical coefficient fit predicts held-out genotypes", loc="left", pad=12)
    ax.text(.04, .96, f"Held-out R² = {data['heldout']['genotype_mean_r2']:.3f} (genotype means)\n"
            f"R² = {data['heldout']['individual_replicate_r2']:.3f} (individual replicates)\n8,192 disjoint genotypes; no refit", transform=ax.transAxes, va="top", fontsize=10)
    ax = axes[0, 1]
    groups = [("Empirical\ncalibration", "calibration", "empirical"), ("Empirical\nconfirmation", "confirmation", "empirical"),
              ("Synthetic\ncalibration", "calibration", "synthetic"), ("Synthetic\nconfirmation", "confirmation", "synthetic")]
    for i, (_, panel, kind) in enumerate(groups):
        frame = data["fits"].query("panel == @panel and kind == @kind")
        values = frame.selected_cv_mean_r2.to_numpy()
        ax.scatter(i + np.linspace(-.1, .1, len(values)) if len(values) > 1 else [i], values,
                   color=COLORS[panel], s=45, marker="o" if kind == "empirical" else "D")
        ax.text(i, values.min() - .045, "λ = " + ", ".join(f"{x:g}" for x in frame.selected_penalty.unique()), ha="center", fontsize=9)
    ax.set_xticks(range(4), [row[0] for row in groups])
    ax.set(ylim=(.56, 1.055), ylabel="Mean 10-fold R² at the selected penalty")
    ax.set_title("B  Auxiliary inference has different prediction tasks", loc="left", pad=12)
    ax.axhline(1, color="#cccccc", linewidth=.8)
    ax.text(.02, .03, "Synthetic observations have no measurement error.\nFour synthetic dots = four seeds; these are CV scores,\nnot the held-out prediction test in A.", transform=ax.transAxes, fontsize=9)
    ax = axes[1, 0]
    boot = data["bootstrap"]
    ax.scatter(boot[:, 0], boot[:, 2], color=COLORS["bootstrap"], s=35, alpha=.7, label="20 genotype bootstraps")
    for panel, marker in (("calibration", "D"), ("confirmation", "s")):
        values = data["targets"][panel]
        ax.scatter(values[0], values[2], color=COLORS[panel], marker=marker, s=70, label=f"Original {panel} fit")
    ax.set(xlabel="log(1 + inferred additive latent variance)", ylabel="Inferred pair share of latent variance", ylim=(-.06, .9))
    ax.set_title("C  Every bootstrap fit shifts to a different architecture", loc="left", pad=12)
    ax.legend(frameon=False, loc="center right", fontsize=9)
    ax = axes[1, 1]
    indices = [0, 2, 4, 5, 9, 14]
    offsets = (boot - data["targets"]["calibration"]) / data["scales"]
    for i, j in enumerate(indices):
        ax.scatter(offsets[:, j], np.full(len(boot), i) + np.linspace(-.13, .13, len(boot)),
                   color=COLORS["bootstrap"], s=15, alpha=.65)
    ax.axvline(0, color="#444444", linewidth=1)
    ax.set_yticks(range(len(indices)), [LABELS[i] for i in indices], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("(Bootstrap summary − original target) / frozen scale")
    ax.set_title("D  Weight scales do not represent the shift from the target", loc="left", pad=12)
    fig.suptitle("Predictive empirical inference does not validate the synthetic ensemble", fontsize=17, x=.04, ha="left", y=.995)
    fig.text(.04, .015, "All selected fits passed the declared convergence checks. This does not prove global optima or coefficient identifiability.\n"
             "Bootstrap scales use variance around the bootstrap mean plus a fixed floor. The shift in C–D is not included as variance around\n"
             "the original target; its cause is not identified here. Thus convergence and a smaller weighted loss are insufficient graduation criteria.", fontsize=10)
    fig.subplots_adjust(left=.07, right=.98, top=.92, bottom=.14, hspace=.5, wspace=.66)
    save_figure(fig, data["figures"], "04_inference_and_weighting")


def plot_search(data):
    fig = plt.figure(figsize=(15, 12))
    grid = fig.add_gridspec(3, 3, height_ratios=[1.1, 1, 1], left=.065, right=.975, top=.92, bottom=.15, hspace=.6, wspace=.27)
    ax = fig.add_subplot(grid[0, :])
    complete = data["search"].query("state == 'COMPLETE'")
    failed = data["search"].query("state == 'FAIL'")
    ax.scatter(complete.number, complete.loss, color=COLORS["calibration"], s=32, label="Completed trial")
    ax.step(complete.number, np.minimum.accumulate(complete.loss), color="#222222", where="post", label="Best completed loss so far")
    ax.scatter(failed.number, np.ones(len(failed)) * 1.025, transform=ax.get_xaxis_transform(), marker="x", color="#B44E4E", s=22, clip_on=False,
               label="Failed trial (no loss assigned)")
    ax.scatter([data["result"]["best_trial"]], [data["result"]["calibration_loss"]], marker="*", s=160, color=COLORS["confirmation"], zorder=5)
    ax.set(yscale="log", xlabel="Optuna trial index (zero-based)", ylabel="Calibration loss", xlim=(-1, 100), ylim=(19, 900))
    ax.set_title("A  Only 37 of 100 attempts produced valid objective values", loc="left", pad=25)
    ax.legend(frameon=False, ncol=3, loc="upper right", fontsize=9)
    ax.text(.69, .16, "Best: trial 63, loss 24.03\nOnly one trial within best + 1", transform=ax.transAxes, fontsize=11)
    for i, parameter in enumerate(PARAMETERS):
        ax = fig.add_subplot(grid[1 + i // 3, i % 3])
        ax.scatter(complete[parameter], complete.loss, color=COLORS["calibration"], alpha=.75, s=22)
        ax.scatter(failed[parameter], np.ones(len(failed)) * 1.025, transform=ax.get_xaxis_transform(), marker="x", color="#B44E4E", s=15,
                   alpha=.65, clip_on=False)
        bestvalue = data["result"]["best_parameters_standardized"][parameter]
        ax.scatter([bestvalue], [data["result"]["calibration_loss"]], marker="*", s=140, color=COLORS["confirmation"], zorder=4)
        low, high = data["config"][f"{parameter}_bounds"]
        if parameter in ("concentration", "steepness", "span"):
            ax.set_xscale("log")
        ax.set(xlim=(low, high), yscale="log", ylim=(19, 900), xlabel=PARAMETER_LABELS[parameter])
        ax.set_title(f"{chr(66+i)}  Selected {bestvalue:.4g}; bounds [{low:g}, {high:g}]", loc="left", fontsize=11, pad=14)
        if i % 3 == 0:
            ax.set_ylabel("Calibration loss")
        ax.grid(axis="y", color="#eeeeee")
    alpha = data["result"]["best_parameters"]["concentration"]
    fig.suptitle("Finite search selects weak, nearly Gaussian interactions; it does not establish an optimum", fontsize=16, x=.04, ha="left", y=.995)
    fig.text(.04, .025, "Blue points are valid trials; red crosses mark numerical failures above the axes (not assigned objective values). Orange stars mark the selected trial.\n"
             "Six panels are projections of the same joint search, not controlled one-parameter sweeps. Bounds are the frozen search bounds.\n"
             f"Selected alpha = {alpha:.2f} is near the Gaussian limit: raw Gamma-variance CV = {1/np.sqrt(alpha):.3f}, raw coefficient excess kurtosis = {3/alpha:.3f}.\n"
             "These formulas precede table centering / RMS normalization. No Gaussian ablation was run. One near-best trial is not a confidence region.", fontsize=10)
    save_figure(fig, data["figures"], "05_search_and_parameters")


def run(root: Path | None = None):
    set_style()
    data = analyze(root)
    plot_architecture(data)
    plot_inference(data)
    plot_search(data)
    return data


if __name__ == "__main__":
    outcome = run()
    print(json.dumps({"losses": outcome["summary"]["losses_recomputed"], "heldout_empirical_model": outcome["heldout"],
                      "figures": str(outcome["figures"]), "tables": str(outcome["tables"])}, indent=2))
