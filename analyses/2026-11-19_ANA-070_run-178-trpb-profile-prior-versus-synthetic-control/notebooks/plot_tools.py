"""Figures for the paired RUN-178 profile / RUN-175 synthetic TrpB comparison.

Inputs are the verified, wide graph/cube merge produced by data_tools. Figures
show every specified setting; no selection or inferential intervals occur here.
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, PercentFormatter
import numpy as np
import pandas as pd


RATIOS = [0, .5, 1, 2, 4, 8]
EPSILONS = [0, .001, .01]
PRIORS = ["Synthetic", "Profile"]
COLORS = {"Synthetic": "#0072B2", "Profile": "#D55E00"}
EMPIRICAL = "#202020"
MODES = ["none", "independent_difference", "shared_stop_sensitivity"]
MODE_LABELS = {
    "none": "No added noise · intrinsic map",
    "independent_difference": "Independent difference · observed map",
    "shared_stop_sensitivity": "Shared-stop sensitivity · observed map",
}
METRICS = {
    "uniform_steps": ("Expected uphill substitutions", "Accepted substitutions", False),
    "uniform_hit_top": ("Upper-1% walk success", "Probability", True),
    "epistatic_fraction": ("Equal-cube epistatic fraction", "Order ≥2 WHT power fraction", True),
    "ho_fraction": ("Equal-cube high-order fraction", "Order ≥3 WHT power fraction", True),
    "pooled_epistatic_fraction": ("Pooled epistatic fraction", "Order ≥2 WHT power fraction", True),
    "pooled_ho_fraction": ("Pooled high-order fraction", "Order ≥3 WHT power fraction", True),
    "sign_epistasis_fraction": ("Sign epistasis", "Sign-epistatic face fraction", True),
    "slope_median": ("Diminishing returns", "Median fixed-direction slope", False),
    "sd_fitness": ("Confirmation local spread", "Mean local fitness SD", False),
    "power_availability": ("WHT availability", "Fraction with resolvable power", True),
    "slope_availability": ("DR availability", "Fraction of estimable slopes", True),
}


def _style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "figure.titlesize": 16, "legend.fontsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .7, "pdf.fonttype": 42,
        "savefig.dpi": 170, "legend.frameon": False,
    })


def _save(fig, folder, name, caption):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fig.text(.025, .018, caption, va="bottom", fontsize=9, color="#454545", linespacing=1.45)
    png = folder / (name + ".png")
    fig.savefig(png, dpi=170, facecolor="white")
    fig.savefig(folder / (name + ".pdf"), facecolor="white")
    plt.close(fig)
    return png


def _legend(fig, y, *, paired=False):
    handles = [Line2D([], [], color=COLORS[p], marker="o", lw=2, ms=5,
                      label=p + (" prior · RUN-175" if p == "Synthetic" else " prior · RUN-178"))
               for p in PRIORS]
    handles.append(Line2D([], [], color=EMPIRICAL, lw=1.5, ls="--", label="Measured TrpB"))
    if paired:
        handles.append(Line2D([], [], color="#929292", lw=1, label="Same generator seed"))
    fig.legend(handles=handles, ncol=len(handles), loc="upper center", bbox_to_anchor=(.5, y))


def _empirical(data, epsilon=.001):
    selected = data.loc[data["layer"].eq("empirical")
                        & data["prior"].eq("Profile")
                        & np.isclose(data["relative_resolution"], epsilon)]
    if len(selected) != 1:
        raise ValueError(f"Expected one Profile-copy empirical row at epsilon={epsilon}; got {len(selected)}")
    return selected.iloc[0]


def _models(data, *, observation="shared_stop_sensitivity", cohort="fresh", epsilon=.001, ratio=None):
    layer = "intrinsic" if observation == "none" else "observed"
    mask = (data["observation"].eq(observation) & data["cohort"].eq(cohort)
            & data["layer"].eq(layer) & np.isclose(data["relative_resolution"], epsilon))
    if ratio is not None:
        mask &= np.isclose(data["ratio"], ratio)
    return data.loc[mask]


def _axis(ax, metric, *, title=None, letter=None):
    label, ylabel, percentage = METRICS[metric]
    ax.set_title((letter + "  " if letter else "") + (title or label), loc="left", pad=11)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#e5e5e5", lw=.6, zorder=0)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_locator(MaxNLocator(5))
    if percentage:
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        lo, hi = ax.get_ylim()
        ax.set_ylim(max(-.015, lo), min(1.035, hi))
    if metric in ["power_availability", "slope_availability"]:
        ax.set_ylim(-.025, 1.04)


def _ratio_axis(ax, *, xlabel=True):
    ax.set_xticks(range(len(RATIOS)), [f"{r:g}" for r in RATIOS])
    ax.set_xlim(-.18, len(RATIOS)-.82)
    if xlabel:
        ax.set_xlabel("Interaction strength r")
    ax.axvline(1, color="#b6b6b6", lw=.8, ls=":", zorder=0)


def _band(ax, subset, metric, prior, *, shade=True):
    z = subset.loc[subset["prior"].eq(prior)]
    if z.duplicated(["ratio", "seed"]).any():
        raise ValueError("Model subset must have one row per ratio and generator seed")
    grouped = z.groupby("ratio")[metric]
    stats = pd.DataFrame({"median": grouped.median(), "lo": grouped.quantile(.1),
                          "hi": grouped.quantile(.9)}).reindex(RATIOS)
    x = np.arange(len(RATIOS))
    if shade:
        ax.fill_between(x, stats["lo"].to_numpy(), stats["hi"].to_numpy(),
                        color=COLORS[prior], alpha=.14, lw=0, zorder=1)
    ax.plot(x, stats["median"].to_numpy(), color=COLORS[prior], marker="o", ms=4,
            lw=1.7, zorder=3, label=prior)


def _reference_line(ax, empirical, metric):
    value = float(empirical[metric])
    if np.isfinite(value):
        ax.axhline(value, color=EMPIRICAL, lw=1.35, ls="--", zorder=2)


def _value(value, metric):
    return f"{value:.1%}" if METRICS[metric][2] else f"{value:.3g}"


def paired_reference(data, figures):
    """Paired seed-level values at the inherited illustrative r=0.5 regime."""
    _style()
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 9.2))
    fig.subplots_adjust(left=.078, right=.985, bottom=.16, top=.79, hspace=.65, wspace=.30)
    fig.suptitle("TrpB: the profile prior changes agreement unevenly across paired seeds", y=.975)
    fig.text(.5, .932, "Inherited r = 0.5 · shared-stop observed maps · 8 fresh seeds · ε/range = 0.001",
             ha="center", fontsize=11)
    _legend(fig, .91, paired=True)
    subset = _models(data, ratio=.5)
    empirical = _empirical(data)
    metrics = ["uniform_steps", "uniform_hit_top", "epistatic_fraction",
               "pooled_ho_fraction", "sign_epistasis_fraction", "slope_median"]
    for i, (ax, metric) in enumerate(zip(axes.flat, metrics)):
        paired = subset.pivot(index="seed", columns="prior", values=metric).reindex(columns=PRIORS)
        paired = paired.replace([np.inf, -np.inf], np.nan).dropna()
        reference = float(empirical[metric])
        error_delta = abs(paired["Profile"]-reference) - abs(paired["Synthetic"]-reference)
        closer = (error_delta < -1e-10).sum()
        jitter = np.linspace(-.07, .07, len(paired))
        for (_, row), shift in zip(paired.iterrows(), jitter):
            ax.plot([shift, 1+shift], row[PRIORS].to_numpy(), color="#a8a8a8", lw=.85, alpha=.8, zorder=1)
        for x, prior in enumerate(PRIORS):
            values = paired[prior].to_numpy()
            ax.scatter(x+jitter, values, color=COLORS[prior], edgecolors="white", linewidths=.4, s=32, zorder=4)
            median = np.median(values) if len(values) else np.nan
            ax.plot([x-.19, x+.19], [median, median], color="white", lw=7, zorder=5)
            ax.plot([x-.19, x+.19], [median, median], color=COLORS[prior], lw=3.6, zorder=6)
        _reference_line(ax, empirical, metric)
        _axis(ax, metric, letter=chr(97+i))
        ax.set_title(f"{chr(97+i)}  {METRICS[metric][0]}\nProfile closer to empirical: {closer}/{len(paired)} seeds",
                     loc="left", fontsize=11, linespacing=1.55, pad=10)
        ax.set_xticks([0, 1], ["Synthetic", "Profile"])
        ax.set_xlim(-.4, 1.4)
        medians = paired.median()
        ax.text(.5, -.23, "Medians: " + _value(medians["Synthetic"], metric) + " → "
                + _value(medians["Profile"], metric) + "  |  empirical " + _value(reference, metric),
                transform=ax.transAxes, ha="center", fontsize=9)
    return _save(fig, figures, "01_paired_inherited_regime",
        "Dots and links pair the same generator seeds; thick colored bars mark medians. “Closer” compares each seed’s absolute error against measured TrpB.\n"
        "An improved ensemble median need not improve most paired realizations. Walks use uniform uphill moves; the upper 1% is defined separately for each map.\n"
        "WHT and DR use confirmation cubes. The inherited regime was highlighted in RUN-175; it was not selected again using the profile outcomes.")


def sweep(data, figures):
    """All strengths under the observation sensitivity used by the prior claim."""
    _style()
    fig, axes = plt.subplots(3, 3, figsize=(14.2, 12.2))
    fig.subplots_adjust(left=.078, right=.985, bottom=.12, top=.86, hspace=.49, wspace=.30)
    fig.suptitle("TrpB: joint search and landscape statistics across every tested strength", y=.98)
    fig.text(.5, .947, "Shared-stop observed maps · 8 fresh seeds · ε/range = 0.001", ha="center", fontsize=11)
    _legend(fig, .927)
    subset, empirical = _models(data), _empirical(data)
    metrics = ["uniform_steps", "uniform_hit_top", "epistatic_fraction", "ho_fraction",
               "pooled_epistatic_fraction", "pooled_ho_fraction", "sign_epistasis_fraction",
               "slope_median", "sd_fitness"]
    for i, (ax, metric) in enumerate(zip(axes.flat, metrics)):
        for prior in PRIORS:
            _band(ax, subset, metric, prior)
        _reference_line(ax, empirical, metric)
        _ratio_axis(ax)
        _axis(ax, metric, letter=chr(97+i))
    return _save(fig, figures, "02_interaction_sweep",
        "Lines are medians; bands span the 10th–90th percentiles of generator realizations, not confidence intervals. All six strengths are retained; spacing is categorical.\n"
        "The dotted guide marks inherited r = 0.5, not a new optimum. Equal-cube fractions weight cubes equally; pooled fractions use summed WHT power.\n"
        "WHT, sign epistasis, DR and walk outcomes were not fitted. Local SD was fitted on calibration cubes; panel i measures its transfer to confirmation cubes.")


def observation(data, figures):
    """Compare readouts fitted under each fixed observation assumption."""
    _style()
    fig, axes = plt.subplots(3, 3, figsize=(14.4, 11.7))
    fig.subplots_adjust(left=.078, right=.985, bottom=.115, top=.855, hspace=.47, wspace=.30)
    fig.suptitle("TrpB: observation assumptions change the prior comparison", y=.979)
    fig.text(.5, .945, "All strengths · 8 fresh seeds · ε/range = 0.001", ha="center", fontsize=11)
    _legend(fig, .925)
    empirical = _empirical(data)
    metrics = ["uniform_steps", "uniform_hit_top", "epistatic_fraction"]
    row_labels = ["Uniform-uphill walk length\n(expected substitutions)",
                  "Upper-1% walk success\n(probability)",
                  "Equal-cube epistatic fraction\n(WHT order ≥2)"]
    for row, metric in enumerate(metrics):
        for col, mode in enumerate(MODES):
            ax = axes[row, col]
            subset = _models(data, observation=mode)
            for prior in PRIORS:
                _band(ax, subset, metric, prior)
            _reference_line(ax, empirical, metric)
            _ratio_axis(ax)
            _axis(ax, metric, letter=chr(97+row*3+col), title=MODE_LABELS[mode])
            ax.set_ylabel(row_labels[row])
        if METRICS[metric][2]:
            common_limits = (-.025, 1.04)
        else:
            common_limits = (min(ax.get_ylim()[0] for ax in axes[row]), max(ax.get_ylim()[1] for ax in axes[row]))
        for ax in axes[row]:
            ax.set_ylim(common_limits)
    return _save(fig, figures, "03_observation_models",
        "Readouts are fitted under each observation assumption; noise donor models stay frozen. No added noise uses intrinsic maps; other columns use observed maps.\n"
        "Lines and bands show medians and 10–90% generator variation, not confidence intervals. Dotted vertical guides mark inherited r = 0.5.\n"
        "Shared-stop is a specified measurement sensitivity, not an identified assay-error mechanism. Agreement after adding it does not identify intrinsic biological epistasis.")


def robustness(data, figures):
    """Resolution and seed-cohort sensitivity at the unchanged r=0.5 setting."""
    _style()
    fig, axes = plt.subplots(2, 5, figsize=(18.2, 8.4))
    fig.subplots_adjust(left=.065, right=.988, bottom=.17, top=.80, hspace=.57, wspace=.48)
    fig.suptitle("TrpB: resolution and seed cohort bound the inherited-regime comparison", y=.974)
    fig.text(.5, .93, "r = 0.5 · shared-stop observed maps · all three recorded resolutions", ha="center", fontsize=11)
    _legend(fig, .905)
    metrics = ["uniform_steps", "uniform_hit_top", "epistatic_fraction", "pooled_ho_fraction", "slope_median"]
    for row, (cohort, n) in enumerate([("fresh", 8), ("reused", 16)]):
        for col, metric in enumerate(metrics):
            ax = axes[row, col]
            for prior, offset in [("Synthetic", -.065), ("Profile", .065)]:
                medians, low, high = [], [], []
                for epsilon in EPSILONS:
                    z = _models(data, cohort=cohort, epsilon=epsilon, ratio=.5)
                    values = z.loc[z["prior"].eq(prior), metric].replace([np.inf, -np.inf], np.nan).dropna()
                    medians.append(values.median())
                    low.append(values.quantile(.1))
                    high.append(values.quantile(.9))
                medians = np.asarray(medians)
                ax.errorbar(np.arange(3)+offset, medians,
                            yerr=np.vstack([medians-np.asarray(low), np.asarray(high)-medians]),
                            color=COLORS[prior], marker="o", ms=4, lw=1.4,
                            elinewidth=1, capsize=3, zorder=3)
            empirical = [float(_empirical(data, epsilon)[metric]) for epsilon in EPSILONS]
            ax.plot(range(3), empirical, color=EMPIRICAL, marker="D", ms=3, ls="--", lw=1.2, zorder=2)
            _axis(ax, metric, letter=chr(97+row*5+col))
            ax.set_xticks(range(3), ["0", "0.001", "0.01"])
            ax.set_xlim(-.3, 2.3)
            ax.set_xlabel("ε / calibration range")
            if col == 0:
                ax.text(-.34, .5, f"{cohort.capitalize()} seeds (n = {n})", rotation=90,
                        transform=ax.transAxes, ha="center", va="center", fontsize=11, fontweight="bold")
    for col in range(len(metrics)):
        limits = (min(ax.get_ylim()[0] for ax in axes[:, col]), max(ax.get_ylim()[1] for ax in axes[:, col]))
        for ax in axes[:, col]:
            ax.set_ylim(limits)
    return _save(fig, figures, "04_sensitivity",
        "Colored markers are cohort medians; error bars span 10–90% generator variation, not confidence intervals. Black diamonds use the matching empirical resolution.\n"
        "Fresh and reused seeds are shown separately; neither cohort is selected for agreement. The middle resolution is primary; all three values use the same frozen range.\n"
        "Expected walks, upper-1% success and descriptor values can change when small effects become unresolved. WHT/DR are retained only when estimable.")


def calibration(data, fits, figures):
    """Separate fitted calibration loss from confirmation spread and availability."""
    _style()
    fig, axes = plt.subplots(4, 3, figsize=(14.4, 14.3))
    fig.subplots_adjust(left=.078, right=.985, bottom=.10, top=.88, hspace=.49, wspace=.30)
    fig.suptitle("TrpB: calibration quality and the availability of unfitted descriptors", y=.983)
    fig.text(.5, .955, "All strengths and observation assumptions · confirmation uses 8 fresh seeds, ε/range = 0.001",
             ha="center", fontsize=11)
    _legend(fig, .935)
    empirical = _empirical(data)
    for col, mode in enumerate(MODES):
        ax = axes[0, col]
        for prior in PRIORS:
            z = fits.loc[fits["prior"].eq(prior) & fits["observation"].eq(mode)]
            if z.duplicated("ratio").any():
                raise ValueError("Expected one fitted readout per prior/observation/ratio")
            values = z.set_index("ratio").reindex(RATIOS)["objective_total"]
            ax.plot(range(6), values, color=COLORS[prior], marker="o", ms=4, lw=1.7)
        _ratio_axis(ax)
        ax.set_title(f"{chr(97+col)}  {MODE_LABELS[mode]}", loc="left", pad=11)
        ax.set_ylabel("Calibration objective (fitted)")
        ax.grid(axis="y", color="#e5e5e5", lw=.6)
        ax.set_ylim(bottom=0)
        for row, metric in enumerate(["sd_fitness", "power_availability", "slope_availability"], start=1):
            ax = axes[row, col]
            subset = _models(data, observation=mode)
            for prior in PRIORS:
                _band(ax, subset, metric, prior)
            _reference_line(ax, empirical, metric)
            _ratio_axis(ax)
            _axis(ax, metric, letter=chr(97+row*3+col), title=METRICS[metric][0])
    spread_limits = (min(ax.get_ylim()[0] for ax in axes[1]), max(ax.get_ylim()[1] for ax in axes[1]))
    for ax in axes[1]:
        ax.set_ylim(spread_limits)
    return _save(fig, figures, "05_calibration_and_availability",
        "Top row: finite-budget fitted loss combines fitness quantiles, mean local SD and floor mass. Lower is better on those targets; it is not an overall realism score.\n"
        "RUN-178 refits all readouts, including r = 0 and r = 4 endpoints; RUN-175 controls are reused. Confirmation local SD tests transfer of a fitted descriptor.\n"
        "Lower rows: medians and 10–90% generator ranges, not confidence intervals. Availability reports estimable WHT cubes and DR slopes; unavailable values are never replaced by zero.")
