"""Reproducible RUN-164 audits, descriptive uncertainty, and scientific figures.

No parameters are fitted here. The primary distance is independently recomputed
from the frozen RUN tables. Bootstrap intervals are conditional on the selected
parameters and the observed cube/seed panel; they do not include refitting.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import math
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

ANALYSIS = Path(__file__).resolve().parents[1]
ROOT = ANALYSIS.parents[1]
BUNDLE = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-164"
FIGURES = ANALYSIS / "figures"
TABLES = ANALYSIS / "tables"
DATASETS = ["gb1_wu2016_four_site_graphfla", "trpb_johnston2024_four_site_graphfla", "cr9114_h1_phillips2021_graphfla"]
NAMES = dict(zip(DATASETS, ["GB1", "TrpB", "CR9114-H1"]))
FAMILIES = ["additive", "bare_gnk", "composite"]
LABELS = {"empirical": "Measured", "additive": "Additive + sigmoid", "bare_gnk": "Bare GNK", "composite": "Composite"}
COLORS = {"empirical": "#222222", "additive": "#0072B2", "bare_gnk": "#D55E00", "composite": "#009E73"}
FEATURES = ["fitness_sd", "walsh_epistatic_fraction", "walsh_high_order_fraction",
            "beneficial_diminishing_returns_slope", "reference_beneficial_fraction",
            "neutral_edge_fraction", "sign_epistasis_fraction", "reciprocal_sign_epistasis_fraction",
            "local_maxima_density", "accessible_path_fraction"]
FEATURE_LABELS = ["Fitness SD", "WHT variance at orders ≥2", "WHT variance at orders ≥3",
                  "Beneficial-mutation slope", "Anchor beneficial fraction", "Neutral-edge fraction",
                  "Sign epistasis", "Reciprocal sign epistasis", "Strict-peak density", "Accessible-path fraction"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                    "axes.spines.right": False, "axes.titlesize": 12, "axes.labelsize": 10,
                    "figure.facecolor": "white", "savefig.facecolor": "white", "pdf.fonttype": 42})


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump_json(path, payload):
    def finite(x):
        if isinstance(x, dict): return {str(k): finite(v) for k, v in x.items()}
        if isinstance(x, (list, tuple, np.ndarray)): return [finite(v) for v in x]
        if isinstance(x, (float, np.floating)): return float(x) if np.isfinite(x) else None
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, np.bool_): return bool(x)
        return x
    Path(path).write_text(json.dumps(finite(payload), indent=2, allow_nan=False) + "\n")


def independent_cube_metrics(values, resolution):
    """Direct orthogonal contrast matrix and fixed-mutation regressions (4 sites)."""
    f = np.asarray(values, float)
    assert f.shape == (16,) and np.isfinite(f).all()
    masks = np.arange(16)
    h = np.array([[(-1) ** int(int(i & j).bit_count()) for j in masks] for i in masks])
    c = h @ f / 16
    energy = np.array([sum(c[i] ** 2 for i in range(1, 16) if i.bit_count() == k) for k in range(1, 5)])
    total = energy.sum()
    resolved = total > max(1e-18, resolution ** 2)
    slopes = []
    for bit in range(4):
        low = masks[(masks & (1 << bit)) == 0]
        x, y = f[low], f[low | (1 << bit)] - f[low]
        if np.var(x) > max(1e-18, resolution ** 2) and y.mean() > resolution:
            slopes.append(np.mean((x - x.mean()) * (y - y.mean())) / np.var(x))
    edges = [f[n ^ (1 << b)] - f[n] for n in masks for b in range(4) if not n & (1 << b)]
    peaks = sum(all(f[n] - f[n ^ (1 << b)] > resolution for b in range(4)) for n in masks)
    return {"total_variance": total, "walsh_epistatic_fraction": energy[1:].sum() / total if resolved else np.nan,
            "walsh_high_order_fraction": energy[2:].sum() / total if resolved else np.nan,
            "beneficial_diminishing_returns_slope": np.mean(slopes) if slopes else np.nan,
            "neutral_edge_fraction": np.mean(np.abs(edges) <= resolution), "local_maxima_density": peaks / 16}


def embedding(reference, rows):
    """The registered training-only scaling and explicit missingness coordinates."""
    ref, x = reference[FEATURES].to_numpy(float), rows[FEATURES].to_numpy(float)
    center = np.array([np.mean(c[np.isfinite(c)]) if np.isfinite(c).any() else 0 for c in ref.T])
    scale = np.array([max(np.std(c[np.isfinite(c)]), .05) if np.isfinite(c).any() else .05 for c in ref.T])
    mask = np.isfinite(x)
    z = np.zeros_like(x)
    np.divide(x - center, scale, out=z, where=mask)
    return np.concatenate([z, mask.astype(float)], axis=1) / math.sqrt(len(FEATURES)), scale


def distance(x, y):
    return max(0., float(2 * cdist(x, y).mean() - cdist(x, x).mean() - cdist(y, y).mean()))


def load_and_audit():
    FIGURES.mkdir(exist_ok=True); TABLES.mkdir(exist_ok=True)
    if not (BUNDLE / "outputs/manifest.json").exists():
        raise FileNotFoundError("Restore data/processed/aggregate-thermodynamic-gnk-calibration/RUN-164.dvc with DVC before running this notebook.")
    output = BUNDLE / "outputs"
    ctx = {name: pd.read_csv(output / file) for name, file in {
        "empirical": "empirical_cubes.csv", "generated": "confirmation.csv", "calibration": "calibration.csv",
        "resolution": "resolution_sensitivity.csv", "scaling": "scaling.csv"}.items()}
    for name, file in {"manifest": "manifest.json", "summary": "summary.json", "selected": "selected_parameters.json"}.items():
        ctx[name] = json.loads((output / file).read_text())
    manifest, summary = ctx["manifest"], ctx["summary"]
    assert manifest["quick"] is False and manifest["synthetic_smoke"] is False
    assert summary["status"] == "DESCRIPTIVE_CONFIRMATION_COMPLETE"
    assert summary["features"] == FEATURES
    cfg = manifest["config"]
    assert cfg["sampling"]["calibration_cubes"] == cfg["sampling"]["confirmation_cubes"] == 24
    assert len(cfg["fit"]["calibration_seeds"]) == 4 and len(cfg["fit"]["confirmation_seeds"]) == 8
    assert not set(cfg["fit"]["calibration_seeds"]) & set(cfg["fit"]["confirmation_seeds"])
    verified = []
    for line in (BUNDLE / "remote-output-sha256.txt").read_text().splitlines():
        expected, relative = line.split(maxsplit=1)
        relative = relative.strip().lstrip("*")
        assert sha256(BUNDLE / relative) == expected, relative
        verified.append(relative)
    assert len(verified) == 11 and (BUNDLE / "logs/labproj-exit-code.txt").read_text().strip() == "0"
    for relative, expected in manifest["source_sha256"].items():
        assert sha256(BUNDLE / "source" / relative) == expected, relative
    input_checks = []
    for path, expected in manifest["input_sha256"].items():
        local = ROOT / "data/processed/ensemble-landscape-validation-datasets/datasets" / Path(path).parent.name / Path(path).name
        if local.exists():
            assert sha256(local) == expected, str(local)
        input_checks.append({"dataset": Path(path).parent.name, "sha256": expected, "local_full_input_verified": local.exists()})
    audit_rows = []
    for dataset in DATASETS:
        rows = ctx["empirical"].query("dataset == @dataset")
        assert len(rows) == 48 and rows.groupby("split").size().to_dict() == {"calibration": 24, "confirmation": 24}
        seen = set()
        for _, row in rows.iterrows():
            nodes = [tuple(v) for v in json.loads(row.node_keys)]
            assert len(nodes) == len(set(nodes)) == 16
            assert not seen.intersection(nodes)
            seen.update(nodes)
            raw = np.array(json.loads(row.raw_values))
            values = (raw - row.raw_offset) / row.raw_scale
            recomputed = independent_cube_metrics(values, cfg["measurement"]["resolution"])
            for k, v in recomputed.items():
                np.testing.assert_allclose(v, row[k], rtol=1e-9, atol=1e-10, equal_nan=True, err_msg=f"{dataset} cube {row.cube} {k}")
            audit_rows.append({"dataset": dataset, "cube": row.cube, "split": row.split,
                               "spectral_resolved": bool(row.spectral_resolved), "metrics_reproduced": True})
        for family in FAMILIES:
            g = ctx["generated"].query("dataset == @dataset and family == @family")
            assert len(g) == 64 and (g.groupby("seed").size() == 8).all()
            assert set(g.seed) == set(cfg["fit"]["confirmation_seeds"])
            selected = next(x for x in ctx["selected"] if x["dataset"] == dataset and x["family"] == family)
            candidates = ctx["calibration"].query("dataset == @dataset and family == @family")
            assert len(candidates) == {"additive": 16, "bare_gnk": 24, "composite": 160}[family]
            np.testing.assert_allclose(candidates.loss.min(), selected["calibration_loss"], atol=1e-12)
            assert any(json.loads(p) == selected["parameters"] and np.isclose(loss, candidates.loss.min())
                       for p, loss in zip(candidates.parameters, candidates.loss))
    audit = {"run": "RUN-164", "remote_files_verified": verified, "frozen_source_files_verified": list(manifest["source_sha256"]),
             "full_inputs": input_checks, "empirical_cubes_recomputed": len(audit_rows), "unique_genotypes_per_assay": 768,
             "calibration_and_confirmation_genotypes_disjoint": True, "generator_seed_panels_disjoint": True,
             "notebook_parameters_refitted": False, "bootstrap_replicates": 2000, "bootstrap_seed": 63001}
    dump_json(TABLES / "audit.json", audit)
    ctx["audit"] = audit
    return ctx


def derive_comparisons(ctx, bootstrap_replicates=2000):
    rng = np.random.default_rng(63001)
    rows, contrasts, descriptors = [], [], []
    ctx["bootstrap"] = {}
    for dataset in DATASETS:
        train = ctx["empirical"].query("dataset == @dataset and split == 'calibration'")
        test = ctx["empirical"].query("dataset == @dataset and split == 'confirmation'")
        x, scale = embedding(train, test)
        wx = rng.multinomial(len(x), np.full(len(x), 1 / len(x)), size=bootstrap_replicates) / len(x)
        # Paired seed-cluster bootstrap: retain the eight cubes within a seed.
        ws = rng.multinomial(8, np.full(8, 1 / 8), size=bootstrap_replicates) / 8
        wy = np.repeat(ws / 8, 8, axis=1)
        dxx = cdist(x, x)
        xterm = np.einsum("bi,ij,bj->b", wx, dxx, wx, optimize=True)
        family_results = {}
        for family in FAMILIES:
            g = ctx["generated"].query("dataset == @dataset and family == @family").sort_values(["seed", "cube"])
            y, _ = embedding(train, g)
            point = distance(x, y)
            reported = next(o for o in ctx["summary"]["outcomes"] if o["dataset"] == dataset and o["family"] == family)
            np.testing.assert_allclose(point, reported["confirmation_joint_loss"], rtol=1e-10, atol=1e-12)
            boot = np.maximum(0, 2 * np.einsum("bi,ij,bj->b", wx, cdist(x, y), wy, optimize=True)
                              - xterm - np.einsum("bi,ij,bj->b", wy, cdist(y, y), wy, optimize=True))
            lower, upper = np.quantile(boot, [.025, .975])
            family_results[family] = (point, boot)
            rows.append({"dataset": dataset, "assay": NAMES[dataset], "family": family, "joint_distance": point,
                         "bootstrap_low": lower, "bootstrap_high": upper, "calibration_distance": reported["calibration_loss"],
                         "empirical_cubes": 24, "generated_cubes": 64, "generator_seeds": 8,
                         "parameters": json.dumps(reported["parameters"], sort_keys=True)})
            for k, s in zip(FEATURES, scale):
                empirical_mean, generated_mean = test[k].mean(), g[k].mean()
                np.testing.assert_allclose(generated_mean - empirical_mean, reported["descriptor_mean_error"][k], rtol=1e-8, atol=1e-9)
                descriptors.append({"dataset": dataset, "assay": NAMES[dataset], "family": family, "feature": k,
                    "empirical_mean": empirical_mean, "generated_mean": generated_mean,
                    "empirical_median": test[k].median(), "generated_median": g[k].median(),
                    "mean_error": generated_mean - empirical_mean, "standardized_mean_error": (generated_mean - empirical_mean) / s,
                    "empirical_resolved_n": int(test[k].notna().sum()), "generated_resolved_n": int(g[k].notna().sum()),
                    "empirical_resolved_fraction": test[k].notna().mean(), "generated_resolved_fraction": g[k].notna().mean()})
        for baseline in ["additive", "bare_gnk"]:
            point = family_results["composite"][0] - family_results[baseline][0]
            boot = family_results["composite"][1] - family_results[baseline][1]
            low, high = np.quantile(boot, [.025, .975])
            contrasts.append({"dataset": dataset, "assay": NAMES[dataset], "comparison": f"composite - {baseline}",
                              "difference": point, "bootstrap_low": low, "bootstrap_high": high})
        ctx["bootstrap"][dataset] = family_results
    ctx["comparison"] = pd.DataFrame(rows); ctx["contrasts"] = pd.DataFrame(contrasts); ctx["descriptors"] = pd.DataFrame(descriptors)
    ctx["comparison"].to_csv(TABLES / "aggregate_summary.csv", index=False)
    ctx["contrasts"].to_csv(TABLES / "model_contrasts.csv", index=False)
    ctx["descriptors"].to_csv(TABLES / "descriptor_summary.csv", index=False)
    return ctx["comparison"]


def replay_scaling(ctx):
    """Replay the frozen scaling fixtures, auditing RMS without squaring underflow."""
    sys.path.insert(0, str(BUNDLE / "source/external/Adapt-Env"))
    from adaptenv import ThermodynamicConfig, ThermodynamicLandscape
    assert Path(sys.modules[ThermodynamicLandscape.__module__].__file__).resolve().is_relative_to((BUNDLE / "source").resolve())
    rows = []
    def stable_rms(v):
        scale = np.max(np.abs(v))
        return float(scale * np.sqrt(np.mean((v / scale) ** 2))) if scale else 0.
    for _, stored in ctx["scaling"].iterrows():
        length, seed, n = int(stored.length), int(stored.seed), int(stored.n_mutations)
        model = ThermodynamicLandscape(ThermodynamicConfig(length=length, seed=seed))
        # The stored fingerprint hashes floating-point field coefficients. Keep
        # cross-platform bitwise identity separate from the numeric checks below.
        fingerprint_matches = model.fingerprint == stored.fingerprint
        rng = np.random.default_rng(seed + 4711)
        src = rng.integers(0, 20, (n, length)); sites = rng.integers(0, length, n)
        dst = src.copy(); dst[np.arange(n), sites] = (dst[np.arange(n), sites] + rng.integers(1, 20, n)) % 20
        a = model.evaluate_batch_components(src); b = model.evaluate_batch_components(dst)
        delta = b["fitness"] - a["fitness"]
        original = float(np.sqrt(np.mean(delta ** 2)))
        for key, actual in [("additive_mutation_rms", stable_rms(b["additive"] - a["additive"])),
                            ("interaction_mutation_rms", stable_rms(b["interaction"] - a["interaction"]))]:
            np.testing.assert_allclose(actual, stored[key], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(original, stored.observable_mutation_rms, rtol=1e-10, atol=0)
        rows.append({**stored.to_dict(), "stable_observable_rms": stable_rms(delta),
                     "replay_fingerprint": model.fingerprint, "bitwise_fingerprint_matches": fingerprint_matches,
                     "max_absolute_observable_delta": float(np.max(np.abs(delta))),
                     "median_logit": float(np.median(a["logit"])), "zero_fitness_fraction": float(np.mean(a["fitness"] == 0)),
                     "zero_delta_fraction": float(np.mean(delta == 0)),
                     "scaling_reference_rank": model.config.reference_rank, "scaling_beta": model.config.steepness,
                     "scaling_ratio": model.config.interaction_ratio, "scaling_p0": model.config.reference_occupancy})
    ctx["scaling_replay"] = pd.DataFrame(rows)
    ctx["scaling_replay"].to_csv(TABLES / "scaling_replay.csv", index=False)
    return ctx["scaling_replay"]


def save(fig, name):
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return FIGURES / f"{name}.png"


def legend_handles(include_empirical=False):
    return [Line2D([], [], color=COLORS[k], lw=2.5, label=LABELS[k])
            for k in (["empirical"] if include_empirical else []) + FAMILIES]


def figure_joint(ctx):
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7), gridspec_kw={"height_ratios": [1.35, 1]}, layout="constrained")
    for col, dataset in enumerate(DATASETS):
        rows = ctx["comparison"].query("dataset == @dataset").set_index("family").loc[FAMILIES]
        ax = axes[0, col]
        for j, (family, r) in enumerate(rows.iterrows()):
            ax.bar(j, r.joint_distance, color=COLORS[family], alpha=.85, width=.6)
            ax.vlines(j, r.bootstrap_low, r.bootstrap_high, color="#252525", lw=1.4)
            ax.plot([j - .07, j + .07], [r.bootstrap_low] * 2, color="#252525", lw=1.3)
            ax.plot([j - .07, j + .07], [r.bootstrap_high] * 2, color="#252525", lw=1.3)
            ax.text(j, r.joint_distance + .035, f"{r.joint_distance:.3f}", ha="center", va="bottom", fontsize=9,
                    bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none", "pad": 1})
        ax.set(xticks=range(3), xticklabels=["Additive", "GNK", "Composite"], title=NAMES[dataset], ylim=(0, 3.7))
        ax.grid(axis="y", alpha=.15); ax.set_axisbelow(True)
        if col == 0: ax.set_ylabel("Held-out joint distance\n(lower is closer)")
        ax = axes[1, col]
        cs = ctx["contrasts"].query("dataset == @dataset")
        ax.axvline(0, color="#777777", lw=1)
        for j, (_, r) in enumerate(cs.iterrows()):
            ax.hlines(j, r.bootstrap_low, r.bootstrap_high, color=COLORS["composite"], lw=2)
            ax.scatter(r.difference, j, color=COLORS["composite"], s=45, zorder=4)
        ax.set(yticks=[0, 1], yticklabels=["vs additive", "vs GNK"], ylim=(-.6, 1.6), xlabel="Composite − baseline distance")
        ax.invert_yaxis(); ax.grid(axis="x", alpha=.15)
    fig.suptitle("Held-out point distances favour composite for GB1 and bare GNK for TrpB / CR9114-H1", fontsize=13)
    fig.supxlabel("Bars/dots: frozen held-out estimates. Whiskers: 95% conditional bootstrap intervals; no refitting. Negative contrasts favour composite.", fontsize=9)
    return save(fig, "01_heldout_comparison")


def figure_mismatch(ctx):
    fig, axes = plt.subplots(1, 3, figsize=(13, 6.2), layout="constrained", sharey=True)
    for ax, dataset in zip(axes, DATASETS):
        rows = ctx["descriptors"].query("dataset == @dataset")
        z = rows.pivot(index="feature", columns="family", values="standardized_mean_error").loc[FEATURES, FAMILIES].to_numpy()
        im = ax.imshow(z, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
        for i in range(len(FEATURES)):
            for j in range(3):
                v = z[i, j]
                ax.text(j, i, f"{v:+.1f}", ha="center", va="center", color="white" if abs(v) > 2.4 else "#222", fontsize=9)
        ax.set(title=NAMES[dataset], xticks=range(3), xticklabels=["Additive", "GNK", "Composite"], yticks=range(10))
    axes[0].set_yticklabels(FEATURE_LABELS)
    fig.colorbar(im, ax=axes, shrink=.7, label="(Generated mean − measured mean) / calibration scale", extend="both")
    fig.suptitle("Lower overall distance still leaves specific aggregate mismatches", fontsize=14)
    fig.supxlabel("Means use available values only. Scale = max(calibration SD, 0.05); colour clipped at ±4, annotations show actual errors.\nMissingness is part of the joint distance and is shown explicitly in Figure 3.", fontsize=9)
    return save(fig, "02_descriptor_mismatches")


def figure_distributions(ctx):
    metrics = ["walsh_epistatic_fraction", "walsh_high_order_fraction", "beneficial_diminishing_returns_slope"]
    labels = ["Fraction of WHT variance at orders ≥2", "Fraction of WHT variance at orders ≥3", "Beneficial-mutation slope (symlog)"]
    fig, axes = plt.subplots(3, 3, figsize=(12.5, 10), layout="constrained")
    for col, dataset in enumerate(DATASETS):
        measured = ctx["empirical"].query("dataset == @dataset and split == 'confirmation'")
        for row, (metric, label) in enumerate(zip(metrics, labels)):
            ax = axes[row, col]; counts = []
            for family in ["empirical"] + FAMILIES:
                frame = measured if family == "empirical" else ctx["generated"].query("dataset == @dataset and family == @family")
                v = np.sort(frame[metric].dropna().to_numpy())
                counts.append(f"{len(v)}/{len(frame)}")
                if len(v):
                    ax.step(np.r_[v[0], v], np.r_[0, np.arange(1, len(v) + 1) / len(v)], where="post",
                            color=COLORS[family], lw=2.3 if family == "empirical" else 1.7)
            ax.text(.98, .04, "Available E/A/G/C\n" + " · ".join(counts), transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=7.6, bbox={"facecolor": "white", "alpha": .9, "edgecolor": "none"})
            ax.set(xlabel=label, ylim=(0, 1.03)); ax.grid(alpha=.13)
            if row < 2: ax.set_xlim(0, 1)
            else:
                ax.set_xscale("symlog", linthresh=.15)
                ax.axvline(0, color="#777", lw=.8, ls=":")
            if row == 0: ax.set_title(NAMES[dataset])
            if col == 0: ax.set_ylabel("Fraction of available cubes")
    fig.legend(handles=legend_handles(True), loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("WHT and diminishing-returns distributions reveal what the joint score hides", fontsize=14)
    return save(fig, "03_wht_and_diminishing_returns")


def figure_resolution(ctx):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4), sharey=True, layout="constrained")
    resolutions = [1e-12, 1e-3, 1e-2]
    for ax, dataset in zip(axes, DATASETS):
        ax.axvspan(.82, 1.18, color="#999", alpha=.12)
        for family in FAMILIES:
            rows = ctx["resolution"].query("dataset == @dataset and family == @family").set_index("resolution")
            ax.plot(range(3), [rows.loc[t, "confirmation_joint_loss"] for t in resolutions], "o-", color=COLORS[family], lw=2)
        ax.set(title=NAMES[dataset], xticks=range(3), xticklabels=["10⁻¹²", "10⁻³\n(primary)", "10⁻²"],
               xlabel="Resolution / calibration fitness range", ylim=(0, 3.1))
        ax.grid(axis="y", alpha=.15)
    axes[0].set_ylabel("Held-out joint distance")
    fig.legend(handles=legend_handles(), loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Model rankings depend on the measurement resolution", fontsize=14)
    return save(fig, "04_resolution_sensitivity")


def figure_scaling(ctx):
    df = ctx["scaling_replay"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    ax = axes[0, 0]
    for key, label, colour in [("build_seconds", "Build", "#0072B2"), ("paired_query_seconds", "512 paired latent queries", "#D55E00")]:
        ax.scatter(df.length, df[key], color=colour, alpha=.55, s=25)
        means = df.groupby("length")[key].mean(); ax.plot(means.index, means, "o-", color=colour, label=label)
    ax.set(xscale="log", yscale="log", ylabel="Seconds", title="A  Cost on the Slurm CPU node")
    ax.legend(frameon=False, fontsize=9)
    ax = axes[0, 1]
    for key, label, colour in [("additive_mutation_rms", "Additive", "#0072B2"), ("interaction_mutation_rms", "GNK", "#D55E00")]:
        means = df.groupby("length")[key].mean()
        ax.scatter(df.length, df[key], color=colour, alpha=.6, s=25)
        ax.plot(means.index, means, "o-", color=colour, label=label)
    ax.axhline(1, color="#555", ls="--", lw=1)
    ax.set(xscale="log", ylim=(.88, 1.12), ylabel="Single-mutation RMS", title="B  Latent normalization remains near one")
    ax.legend(frameon=False, fontsize=9)
    ax = axes[1, 0]
    for _, r in df.iterrows():
        if r.stable_observable_rms > 0:
            ax.scatter(r.length, np.log10(r.stable_observable_rms), color="#009E73", s=35)
        else:
            ax.scatter(r.length, -330, marker="x", color="#009E73", s=45)
    ax.axhline(-320, color="#aaa", lw=.8, ls=":")
    ax.set(xscale="log", ylim=(-343, 0), ylabel="log₁₀ observable mutation RMS", title="C  Observable effects collapse on random backgrounds",
           yticks=[0, -50, -100, -150, -200, -250, -300, -330], yticklabels=["0", "−50", "−100", "−150", "−200", "−250", "−300", "0*"])
    ax.text(.02, .06, "*Crosses in bottom strip: exactly zero in float64.\nStable RMS avoids underflow from squaring tiny effects.",
            transform=ax.transAxes, fontsize=8, bbox={"facecolor": "white", "alpha": .9, "edgecolor": "none"})
    ax = axes[1, 1]
    ax.scatter(df.length, df.median_logit, color="#222", s=30)
    means = df.groupby("length").median_logit.mean(); ax.plot(means.index, means, "o-", color="#222")
    ax.axhline(math.log(.99 / .01), color="#777", ls="--", label="Reference logit")
    ax.set(xscale="log", ylabel="Median logit on sampled backgrounds", title="D  Backgrounds move far below the reference")
    ax.legend(frameon=False, fontsize=9)
    for ax in axes.flat:
        ax.set_xticks([32, 128, 512], ["32", "128", "512"]); ax.set_xlabel("Sequence length"); ax.grid(alpha=.12)
        ax.xaxis.set_minor_locator(NullLocator())
    fig.suptitle("Scalable evaluation and latent effects do not establish scalable observable variation", fontsize=14)
    fig.supxlabel("Scaling fixtures use default reference rank 1, β=4, r=0.25, p₀=0.99; these differ from the fitted four-site models. Three seeds per length.", fontsize=9)
    return save(fig, "05_scaling_and_saturation")


def write_metrics(ctx):
    metrics = {
        "analysis": "ANA-063", "run": "RUN-164", "status": "descriptive_analysis_complete",
        "audit": ctx["audit"], "heldout": ctx["comparison"].to_dict("records"),
        "paired_conditional_bootstrap": ctx["contrasts"].to_dict("records"),
        "scaling_replay": ctx["scaling_replay"].to_dict("records"),
        "limits": ["24 confirmation cubes per assay and 8 synthetic seeds; intervals are conditional, not refitted",
                   "Cubes are genotype-disjoint but may share sites and alleles; intervals do not establish population-wide coverage",
                   "GB1 has unresolved spectra and only six available beneficial-slope cube summaries at primary resolution",
                   "Resolution sensitivity changes model rankings; no noise model or replicate correction is available",
                   "Four-active-site fits and separate default-parameter scaling do not prove large-L empirical recovery",
                   "Scaling replay checks recorded RMS numerically at rtol 1e-10; bitwise coefficient fingerprints may differ across platforms",
                   "No empirical residue interaction identities, physical free energies, or whole-protein WHT inferred"],
        "analysis_python": sys.version,
        "analysis_packages": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy", "matplotlib", "nbformat", "nbclient", "ipykernel"]},
        "analysis_source_sha256": sha256(__file__),
    }
    dump_json(TABLES / "metrics.json", metrics)
    return metrics
