#!/usr/bin/env python3
"""EXP-065 follow-up: paired readouts, conditional statistics, embedded cubes.

All scientific comparisons are descriptive. Assertions protect data separation
and pairing, not a scientific acceptance threshold. The original RUN-164 runner
and its outputs are not changed by this follow-up.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from aggregate_thermodynamic_calibration import (
    FEATURES, JointDistance, ThermodynamicConfig, ThermodynamicLandscape,
    clean, describe, digest, measured_cubes, sigmoid, write_json,
)
import numpy as np
import pandas as pd
import scipy
import yaml

DIM = 4
BITS = np.array([[(n >> i) & 1 for i in range(DIM)] for n in range(16)], bool)


def stable_rms(values):
    x = np.asarray(values, float)
    scale = np.max(np.abs(x))
    return float(scale * np.sqrt(np.mean((x / scale) ** 2))) if scale else 0.0


def affine_calibration(raw_cubes, empirical_cubes):
    """Match mean within-cube SD and overall mean, using calibration cubes only."""
    raw = np.asarray(raw_cubes, float)
    empirical = np.asarray(empirical_cubes, float)
    spread = float(raw.std(axis=1).mean())
    if spread <= 1e-14:
        raise ValueError("readout has no numerically usable calibration spread")
    scale = float(empirical.std(axis=1).mean()) / spread
    if scale <= 0:
        raise ValueError("empirical calibration panel has no within-cube variation")
    return {"scale": scale, "offset": float(empirical.mean() - scale * raw.mean())}


def raw_readout(a, g, parameter):
    z = np.asarray(a) + parameter["ratio"] * np.asarray(g)
    if parameter["readout"] == "linear":
        return z
    mu = math.log(parameter["p0"] / (1 - parameter["p0"]))
    return sigmoid(mu + parameter["beta"] * z)


def observe(a, g, parameter):
    return parameter["offset"] + parameter["scale"] * raw_readout(a, g, parameter)


def fresh_empirical(path, excluded, *, count, seed):
    """Draw disjoint candidates, then reject any cube containing prior RUN nodes."""
    candidates, audits = measured_cubes(path, dimension=DIM, count=count * 3,
                                        seed=seed, missing_fitness="exclude")
    retained, records = [], []
    rejected = 0
    for cube, audit in zip(candidates, audits):
        if any(tuple(n) in excluded for n in audit["node_keys"]):
            rejected += 1
            continue
        retained.append(cube)
        records.append(audit)
        if len(retained) == count:
            break
    if len(retained) != count:
        raise ValueError("Insufficient new cubes after excluding prior RUN genotypes")
    return np.asarray(retained), records, rejected


def cube_nodes(anchor, alternatives, sites):
    sites = np.asarray(sites, int)
    nodes = np.repeat(np.asarray(anchor)[None, :], 16, axis=0)
    nodes[:, sites] = np.where(BITS, np.asarray(alternatives)[sites], np.asarray(anchor)[sites])
    return nodes


def shared_support_count(model, sites):
    selected = set(map(int, sites))
    supports = {term.support for term in model.field.terms} if model.field else set()
    return sum(len(selected.intersection(s)) >= 2 for s in supports)


def make_model(cfg, length, seed, reference_rank):
    return ThermodynamicLandscape(ThermodynamicConfig(
        length=length, seed=seed, reference_rank=reference_rank,
        representation_rank=cfg["model"]["representation_rank"],
        order_weights=tuple(cfg["model"]["order_weights"])))


def active_panels(cfg, seeds, split):
    rows = []
    for seed in seeds:
        model = make_model(cfg, DIM, seed, cfg["model"]["reference_rank"])
        rng = np.random.default_rng(seed + 87111)
        for cube in range(cfg["sampling"]["synthetic_cubes_per_seed"]):
            anchor = rng.integers(0, 20, DIM)
            alt = (anchor + rng.integers(1, 20, DIM)) % 20
            a, g = model.latent_components(cube_nodes(anchor, alt, np.arange(DIM)))
            rows.append({"phase": "active_four_site", "split": split, "length": DIM,
                         "seed": seed, "cube": cube, "library": 0, "site_policy": "all_four",
                         "ensemble": "uniform", "reference_rank": cfg["model"]["reference_rank"],
                         "sites": json.dumps(list(range(DIM))), "anchor": json.dumps(anchor.tolist()),
                         "alternatives": json.dumps(alt.tolist()), "hamming_distance": int(np.sum(anchor != model.reference)),
                         "shared_supports": shared_support_count(model, range(DIM)),
                         "fingerprint": model.fingerprint, "a": a, "g": g})
    return rows


def embedded_panels(cfg):
    """Libraries fix four mutable sites; backgrounds vary at controlled distances.

    Uniform-site libraries and libraries occupying a sampled four-site GNK
    support are interpretation controls, not inferred protein contact maps.
    """
    for length, rank, seed in itertools.product(cfg["embedding"]["lengths"],
            cfg["embedding"]["reference_ranks"], cfg["embedding"]["seeds"]):
        start = time.perf_counter()
        model = make_model(cfg, length, seed, rank)
        elapsed = time.perf_counter() - start
        supports = sorted({t.support for t in model.field.terms if len(t.support) == DIM})
        rng = np.random.default_rng(seed + 19001 + length)
        for library in range(cfg["embedding"]["libraries_per_seed"]):
            policies = {"uniform_sites": sorted(rng.choice(length, DIM, replace=False).tolist()),
                        "shared_support": list(supports[int(rng.integers(len(supports)))])}
            for distance in [*cfg["embedding"]["distances"], None]:
                for cube in range(cfg["embedding"]["cubes_per_library"]):
                    anchor = model.reference.copy()
                    if distance is None:
                        anchor = rng.integers(0, 20, length)
                    elif distance:
                        changed = rng.choice(length, distance, replace=False)
                        anchor[changed] = (anchor[changed] + rng.integers(1, 20, distance)) % 20
                    alt = (anchor + rng.integers(1, 20, length)) % 20
                    # Exactly the same background and alternatives across site policies.
                    for policy, sites in policies.items():
                        nodes = cube_nodes(anchor, alt, sites)
                        a, g = model.latent_components(nodes)
                        yield {"phase": "embedded", "split": "transfer", "length": length,
                               "seed": seed, "cube": cube, "library": library, "site_policy": policy,
                               "ensemble": "uniform" if distance is None else f"distance_{distance}",
                               "reference_rank": rank, "sites": json.dumps(sites),
                               "anchor": json.dumps(anchor.tolist()), "alternatives": json.dumps(alt.tolist()),
                               "hamming_distance": int(np.sum(anchor != model.reference)),
                               "shared_supports": shared_support_count(model, sites),
                               "fingerprint": model.fingerprint, "build_seconds": elapsed, "a": a, "g": g}
        print(f"Embedded panel: L={length}, reference rank={rank}, seed={seed}", flush=True)


def metric_rows(cubes, resolution):
    return [describe(cube, DIM, resolution) for cube in cubes]


def fit_pairs(cfg, train_panels, empirical):
    a = np.asarray([r["a"] for r in train_panels]); g = np.asarray([r["g"] for r in train_panels])
    resolution = cfg["measurement"]["resolution"]
    train_metrics = metric_rows(empirical, resolution)
    distance = JointDistance(train_metrics)
    fits, candidates = [], []
    for ratio in cfg["fit"]["interaction_ratios"]:
        for readout in ("linear", "sigmoid"):
            grid = [(None, None)] if readout == "linear" else itertools.product(
                cfg["fit"]["steepness"], cfg["fit"]["reference_occupancy"])
            scored = []
            for beta, p0 in grid:
                parameter = {"ratio": ratio, "readout": readout, "beta": beta, "p0": p0}
                raw = raw_readout(a, g, parameter)
                try:
                    parameter.update(affine_calibration(raw, empirical))
                    loss = distance(train_metrics, metric_rows(observe(a, g, parameter), resolution))
                except ValueError as exc:
                    candidates.append({**parameter, "loss": math.nan, "excluded_reason": str(exc)})
                    continue
                candidates.append({**parameter, "loss": loss, "excluded_reason": ""})
                scored.append({**parameter, "calibration_loss": loss})
            if not scored:
                raise ValueError(f"No usable calibration candidates for {readout}, ratio {ratio}")
            fits.append(min(scored, key=lambda p: p["calibration_loss"]))
    return fits, candidates, distance


def add_observation(rows, values_rows, edge_rows, base, values, parameter, resolutions):
    meta = {k: v for k, v in base.items() if k not in ("a", "g", "anchor", "alternatives")}
    meta.update({k: parameter[k] for k in ("readout", "ratio", "beta", "p0", "offset", "scale")})
    primary = resolutions[0]
    low = np.concatenate([np.array([i for i in range(16) if not i & (1 << bit)]) for bit in range(DIM)])
    bits = np.repeat(np.arange(DIM), 8)
    high = low | (1 << bits)
    delta = values[high] - values[low]
    for resolution in resolutions:
        rows.append({**meta, "resolution": resolution, **describe(values, DIM, resolution),
                     "edge_effect_rms_stable": stable_rms(delta),
                     "zero_edge_fraction": float(np.mean(delta == 0)),
                     "fitness_median": float(np.median(values))})
    values_rows.append({**meta, "values": json.dumps(values.tolist())})
    for bit, source, target, effect in zip(bits, low, high, delta):
        edge_rows.append({**meta, "mutation": int(bit), "source_node": int(source),
                          "target_node": int(target), "background_fitness": float(values[source]),
                          "effect": float(effect), "resolution_for_strata": primary})


def conditional_table(frame, empirical_calibration):
    boundaries = np.unique(np.quantile([r["fitness_mean"] for r in empirical_calibration], [.25, .5, .75]))
    edges = np.r_[-np.inf, boundaries, np.inf]
    frame = frame.copy()
    frame["background_bin"] = pd.cut(frame.fitness_mean, edges, labels=False, include_lowest=True)
    groups = ["phase", "split", "readout", "ratio", "p0", "length", "reference_rank",
              "site_policy", "ensemble", "resolution", "background_bin"]
    metrics = ["fitness_mean", "fitness_sd", "walsh_epistatic_fraction", "walsh_high_order_fraction",
               "beneficial_diminishing_returns_slope", "sign_epistasis_fraction", "neutral_edge_fraction",
               *[f"walsh_energy_order_{k}" for k in range(1, 5)]]
    result = []
    for keys, panel in frame.groupby(groups, dropna=False):
        record = dict(zip(groups, keys)); record["n_cubes"] = len(panel)
        for metric in metrics:
            available = panel[metric].dropna()
            record[f"{metric}_mean"] = available.mean() if len(available) else math.nan
            record[f"{metric}_median"] = available.median() if len(available) else math.nan
            record[f"{metric}_available"] = len(available)
        result.append(record)
    return pd.DataFrame(result), boundaries


def figures(frame, comparisons, selected_pair, output, resolution, quick):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
    output = output / "figures"; output.mkdir(exist_ok=True)
    prefix = "SMOKE ONLY — " if quick else ""
    colors = {"linear": "#0072B2", "sigmoid": "#D55E00", "empirical": "#222222"}
    primary = frame[(frame.resolution == resolution) & (frame.split == "confirmation")]
    empirical = primary[primary.phase == "empirical"]
    active = primary[primary.phase == "active_four_site"]
    comparison = pd.DataFrame(comparisons)
    comparison = comparison[comparison.resolution == resolution]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for readout in ["linear", "sigmoid"]:
        part = comparison[comparison.readout == readout].sort_values("ratio")
        axes[0].plot(part.ratio, part.confirmation_joint_distance, "o-", color=colors[readout], label=readout)
        for ax, key in zip(axes[1:], ["fitness_sd", "walsh_high_order_fraction"]):
            group = active[active.readout == readout].groupby("ratio")[key].mean()
            ax.plot(group.index, group.values, "o-", color=colors[readout])
    for ax, key in zip(axes[1:], ["fitness_sd", "walsh_high_order_fraction"]):
        ax.axhline(empirical[key].mean(), color="black", linestyle="--", label="measured mean")
    for ax, label in zip(axes, ["Joint descriptor distance", "Mean within-cube fitness SD", "Mean WHT fraction at orders ≥3"]):
        ax.set(xlabel="Latent interaction/additive ratio", ylabel=label); ax.grid(alpha=.15)
    axes[0].legend(); axes[1].legend()
    fig.suptitle(prefix + "TrpB: paired readouts after calibration of observable amplitude")
    for ext in ["png", "pdf"]: fig.savefig(output / f"01_paired_readouts.{ext}", dpi=160)
    plt.close(fig)
    ratio = selected_pair[0]["ratio"]
    chosen = active[active.ratio == ratio]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    specs = [("high_order_power", "Absolute WHT power at orders ≥3"),
             ("walsh_high_order_fraction", "WHT variance fraction at orders ≥3"),
             ("beneficial_diminishing_returns_slope", "Beneficial-mutation slope")]
    for row, panel in enumerate([chosen, primary[primary.phase == "operating_position"]]):
        group_col = "readout" if row == 0 else "p0"
        palette = colors if row == 0 else {.1:"#009E73", .5:"#CC79A7", .9:"#E69F00"}
        for ax, (key, label) in zip(axes[row], specs):
            for name, group in panel.groupby(group_col):
                values = group.walsh_energy_order_3 + group.walsh_energy_order_4 if key == "high_order_power" else group[key]
                ax.scatter(group.fitness_mean, values, s=17, alpha=.6, color=palette[name], label=str(name))
            values = empirical.walsh_energy_order_3 + empirical.walsh_energy_order_4 if key == "high_order_power" else empirical[key]
            ax.scatter(empirical.fitness_mean, values, color="black", marker="x", label="measured", s=28)
            ax.set(xlabel="Mean cube fitness", ylabel=label); ax.grid(alpha=.15)
            if key == "high_order_power": ax.set_yscale("symlog", linthresh=1e-7)
            if key == "beneficial_diminishing_returns_slope": ax.set_yscale("symlog", linthresh=.2)
        axes[row, 0].legend(title="Readout" if row == 0 else "Operating occupancy")
    fig.suptitle(prefix + f"Conditional relationships on held-out cubes (paired ratio r={ratio:g})")
    for ext in ["png", "pdf"]: fig.savefig(output / f"02_conditional_relationships.{ext}", dpi=160)
    plt.close(fig)
    embedded = frame[(frame.phase == "embedded") & (frame.reference_rank == .5) & (frame.resolution == resolution)]
    ensemble_colors = {"distance_0":"#0072B2", "distance_4":"#009E73", "distance_16":"#D55E00", "uniform":"#CC79A7"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for row, readout in enumerate(["linear", "sigmoid"]):
        panel = embedded[embedded.readout == readout]
        for ax, key in zip(axes[row], ["edge_effect_rms_stable", "walsh_high_order_fraction"]):
            for (ensemble, policy), group in panel.groupby(["ensemble", "site_policy"]):
                means = group.groupby("length")[key].mean()
                ax.plot(means.index, means.values, marker="o", linestyle="-" if policy == "uniform_sites" else "--", color=ensemble_colors[ensemble])
            ax.set(xlabel="Sequence length", ylabel=f"{readout}: " + ("stable edge RMS" if key == "edge_effect_rms_stable" else "WHT fraction at orders ≥3"))
            ax.set_xscale("log", base=2); ax.set_xticks(sorted(panel.length.unique()), labels=sorted(panel.length.unique())); ax.grid(alpha=.15)
            if key == "edge_effect_rms_stable": ax.set_yscale("symlog", linthresh=1e-6)
    handles = [Line2D([0], [0], color=c, label=e.replace("_", " ")) for e,c in ensemble_colors.items()]
    handles += [Line2D([0],[0],color="black",linestyle=s,label=l) for s,l in [("-","uniform sites"),("--","shared support")]]
    axes[0,0].legend(handles=handles, fontsize=8)
    fig.suptitle(prefix + "Frozen readouts in longer landscapes: reference rank 0.5\nMeans across seeds/libraries; rank-1 stress controls retained in tables")
    for ext in ["png", "pdf"]: fig.savefig(output / f"03_embedded_neighbourhoods.{ext}", dpi=160)
    plt.close(fig)


def run(config, output, *, data_root=None, quick=False):
    cfg = json.loads(json.dumps(config))
    if quick:
        cfg["sampling"].update(seed=659501, calibration_cubes=4, confirmation_cubes=4, synthetic_cubes_per_seed=2)
        cfg["fit"].update(calibration_seeds=[659601, 659602], confirmation_seeds=[659701, 659702],
                          interaction_ratios=[0.0, 1.0], steepness=[.125, .5], reference_occupancy=[.5])
        cfg["embedding"].update(lengths=[16, 32], seeds=[659801], libraries_per_seed=1,
                                cubes_per_library=1, distances=[0, 4])
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    data_root = Path(data_root or ROOT / cfg["data"]["root"])
    path = data_root / "datasets" / cfg["data"]["dataset"] / "variants.parquet"
    if not path.exists():
        path = path.with_suffix(".csv")
    prior_path = ROOT / cfg["data"]["prior_nodes"]
    prior = json.loads(prior_path.read_text())
    excluded = {tuple(n) for n in prior["node_keys"]}
    if not quick:
        excluded.update(tuple(n) for n in prior.get("engineering_node_keys", []))
    ntrain, ntest = cfg["sampling"]["calibration_cubes"], cfg["sampling"]["confirmation_cubes"]
    train_seeds, test_seeds = cfg["fit"]["calibration_seeds"], cfg["fit"]["confirmation_seeds"]
    assert not set(train_seeds) & set(test_seeds)
    raw, audit, rejected = fresh_empirical(path, excluded, count=ntrain + ntest, seed=cfg["sampling"]["seed"])
    all_nodes = [tuple(n) for row in audit for n in row["node_keys"]]
    assert len(all_nodes) == len(set(all_nodes)) and not excluded.intersection(all_nodes)
    offset = float(raw[:ntrain].min()); scale = float(np.ptp(raw[:ntrain]))
    if scale <= 0:
        raise ValueError("Cannot define scale from a flat empirical calibration panel")
    empirical = (raw - offset) / scale
    resolutions = list(dict.fromkeys([cfg["measurement"]["resolution"], *cfg["measurement"]["sensitivity_resolutions"]]))
    rows, value_rows, edge_rows, latent_rows = [], [], [], []
    empirical_rows = []
    empirical_parameter = {"readout": "empirical", "ratio": None, "beta": None, "p0": None, "offset": offset, "scale": scale}
    for index, (values, record) in enumerate(zip(empirical, audit)):
        base = {"phase": "empirical", "split": "calibration" if index < ntrain else "confirmation",
                "length": DIM, "seed": cfg["sampling"]["seed"], "cube": index, "library": 0,
                "site_policy": "measured_four_sites", "ensemble": "observed_complete",
                "reference_rank": None, "cube_id": f"empirical_{index}"}
        add_observation(rows, value_rows, edge_rows, base, values, empirical_parameter, resolutions)
        empirical_rows.append({**base, "node_keys": json.dumps(record["node_keys"]),
                               "raw_values": json.dumps(raw[index].tolist()), "raw_offset": offset,
                               "raw_scale": scale, "sampling_audit": json.dumps(record)})
    print(f"Sampled {ntrain}+{ntest} fresh empirical cubes; excluded {len(excluded)} prior/exploratory genotypes", flush=True)
    train_panels = active_panels(cfg, train_seeds, "calibration")
    test_panels = active_panels(cfg, test_seeds, "confirmation")
    for i, base in enumerate(train_panels + test_panels):
        base["cube_id"] = f"active_{i}"
    fits, candidates, distance = fit_pairs(cfg, train_panels, empirical[:ntrain])
    print(f"Calibrated {len(candidates)} candidates; {len(fits)} paired readout/ratio arms", flush=True)
    for base in train_panels + test_panels:
        latent_rows.append({**{k: v for k, v in base.items() if k not in ("a", "g")},
                            "additive": json.dumps(base["a"].tolist()), "gnk": json.dumps(base["g"].tolist())})
        for parameter in fits:
            add_observation(rows, value_rows, edge_rows, base, observe(base["a"], base["g"], parameter), parameter, resolutions)
    chosen_sigmoid = min((p for p in fits if p["readout"] == "sigmoid"), key=lambda p: p["calibration_loss"])
    chosen_linear = next(p for p in fits if p["readout"] == "linear" and p["ratio"] == chosen_sigmoid["ratio"])
    for base in test_panels:
        for p0 in cfg["conditional"]["operating_occupancies"]:
            parameter = {**chosen_sigmoid, "p0": p0}
            shifted = {**base, "phase": "operating_position"}
            add_observation(rows, value_rows, edge_rows, shifted, observe(base["a"], base["g"], parameter), parameter, resolutions)
    for i, base in enumerate(embedded_panels(cfg)):
        base["cube_id"] = f"embedded_{i}"
        latent_rows.append({**{k: v for k, v in base.items() if k not in ("a", "g")},
                            "additive": json.dumps(base["a"].tolist()), "gnk": json.dumps(base["g"].tolist())})
        for parameter in (chosen_linear, chosen_sigmoid):
            values = observe(base["a"], base["g"], parameter)
            mu = 0.0 if parameter["readout"] == "linear" else math.log(parameter["p0"] / (1 - parameter["p0"]))
            latent = base["a"] + parameter["ratio"] * base["g"]
            enriched = {**base, "latent_mean": float(latent.mean()),
                        "logit_mean": float(np.mean(mu + parameter["beta"] * latent)) if parameter["beta"] is not None else math.nan}
            add_observation(rows, value_rows, edge_rows, enriched, values, parameter, resolutions)
    frame = pd.DataFrame(rows)
    outcomes = []
    for parameter in fits:
        for resolution in resolutions:
            measured = frame[(frame.phase == "empirical") & (frame.split == "confirmation") & (frame.resolution == resolution)]
            generated = frame[(frame.phase == "active_four_site") & (frame.split == "confirmation") &
                (frame.readout == parameter["readout"]) & (frame.ratio == parameter["ratio"]) & (frame.resolution == resolution)]
            outcomes.append({**parameter, "resolution": resolution,
                             "confirmation_joint_distance": distance(measured.to_dict("records"), generated.to_dict("records")),
                             "n_measured": len(measured), "n_generated": len(generated)})
    primary = frame[frame.resolution == resolutions[0]]
    conditioning, boundaries = conditional_table(frame, primary[(primary.phase == "empirical") & (primary.split == "calibration")].to_dict("records"))
    pd.DataFrame(empirical_rows).to_csv(output / "empirical_cubes.csv", index=False)
    pd.DataFrame(candidates).to_csv(output / "calibration.csv", index=False)
    frame.to_csv(output / "cube_statistics.csv", index=False)
    pd.DataFrame(value_rows).to_csv(output / "cube_values.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(output / "mutation_effects.csv", index=False)
    pd.DataFrame(latent_rows).to_csv(output / "latent_cubes.csv", index=False)
    pd.DataFrame(outcomes).to_csv(output / "readout_comparison.csv", index=False)
    conditioning.to_csv(output / "conditional_summary.csv", index=False)
    figures(frame, outcomes, [chosen_linear, chosen_sigmoid], output, resolutions[0], quick)
    write_json(output / "selected_parameters.json", {"within_ratio_pairs": fits,
        "embedded_pair": [chosen_linear, chosen_sigmoid], "conditional_bin_boundaries": boundaries,
        "affine_calibration": "calibration-only mean within-cube SD and overall mean; no clipping; frozen for transfer"})
    source_files = ["scripts/paired_readout_exp065.py", "scripts/aggregate_thermodynamic_calibration.py",
        "src/project_adapt_env/dense_descriptive.py", "src/project_adapt_env/summary_statistics.py",
        "external/Adapt-Env/adaptenv/thermodynamic.py", "external/Adapt-Env/adaptenv/extensive_gnk.py",
        "external/Adapt-Env/adaptenv/ranked_categorical_gnk.py", config["protocol_path"],
        str(Path(config["protocol_path"]).with_suffix(".yaml"))]
    write_json(output / "manifest.json", {"experiment": "EXP-065", "protocol": "paired_readout_followup_v1",
        "quick": quick, "config": cfg, "python": platform.python_version(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "input_sha256": {str(path): digest(path), str(prior_path): digest(prior_path)},
        "source_sha256": {p: digest(ROOT / p) for p in source_files},
        "prior_run": prior["run"], "prior_run_nodes_excluded": len(prior["node_keys"]),
        "known_nodes_excluded": len(excluded), "engineering_nodes_excluded": 0 if quick else len(prior.get("engineering_node_keys", [])),
        "overlapping_candidate_cubes_rejected": rejected,
        "fresh_genotypes": len(all_nodes), "genotype_and_seed_splits_disjoint": True,
        "run_id": __import__("os").environ.get("LABPROJ_RUN_ID"),
        "noise_information": "Prepared assay table contains no independent replicate errors; no noise correction or invented noise model",
        "metric_features": FEATURES})
    write_json(output / "summary.json", {"status": "SMOKE_ONLY" if quick else "DESCRIPTIVE_FOLLOWUP_COMPLETE",
        "empirical_cubes": ntrain + ntest, "selected_readout_ratio_arms": len(fits),
        "latent_cubes": len(latent_rows), "observed_cube_variants": len(value_rows),
        "statistic_rows": len(frame), "paired_confirmation": outcomes,
        "scope": "Aggregate TrpB comparison, operating-position diagnostic and frozen-parameter embedded transfer; no empirical interaction identities or success gates"})
    print(f"Wrote all follow-up outputs to {output}", flush=True)
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    run(yaml.safe_load(args.config.read_text()), args.output_dir, data_root=args.data_root, quick=args.quick)
