"""Independent, observed-panel graph calculations for ANA-071.

This analysis module does not import the fitted model or its diagnostic function.
Only genotypes already present in the archived confirmation panel are joined.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PanelGraph:
    sequences: np.ndarray
    edges: np.ndarray  # lower allele index -> higher allele index at one site
    edge_sites: np.ndarray
    faces: np.ndarray  # base, change_i, change_j, change_i_and_j; i < j
    degrees: np.ndarray


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def panel_graph(sequences, alphabet_size=20):
    """Join sparse observed genotypes using integer addresses, without scoring holes.

    The lookup is only an address-to-observed-row index (160,000 integers for
    the four-site panel). A value of -1 marks an absent genotype. This is not
    landscape completion: no fitness value is assigned to an absent genotype.
    Edge enumeration is vectorized over observed rows; each complete face is
    emitted once in increasing site and allele order.
    """
    x = np.asarray(sequences, dtype=np.int64)
    assert x.ndim == 2 and np.all((x >= 0) & (x < alphabet_size))
    n, length = x.shape
    powers = alphabet_size ** np.arange(length - 1, -1, -1, dtype=np.int64)
    codes = x @ powers
    assert len(np.unique(codes)) == n
    address = np.full(alphabet_size ** length, -1, dtype=np.int64)
    address[codes] = np.arange(n)
    edge_blocks, site_blocks = [], []
    for site in range(length):
        for allele in range(alphabet_size):
            left = np.flatnonzero(x[:, site] < allele)
            target = codes[left] + (allele - x[left, site]) * powers[site]
            right = address[target]
            keep = right >= 0
            edge_blocks.append(np.column_stack((left[keep], right[keep])))
            site_blocks.append(np.full(int(keep.sum()), site, dtype=np.int64))
    edges = np.concatenate(edge_blocks)
    edge_sites = np.concatenate(site_blocks)
    face_blocks = []
    for i, j in combinations(range(length), 2):
        e = edges[edge_sites == i]
        for allele in range(alphabet_size):
            left = e[x[e[:, 0], j] < allele]
            delta = (allele - x[left[:, 0], j]) * powers[j]
            single = address[codes[left[:, 0]] + delta]
            double = address[codes[left[:, 1]] + delta]
            keep = (single >= 0) & (double >= 0)
            face_blocks.append(np.column_stack((left[keep], single[keep], double[keep])))
    faces = np.concatenate(face_blocks)
    assert len(edges) == len(np.unique(np.sort(edges, axis=1), axis=0))
    assert len(faces) == len(np.unique(np.sort(faces, axis=1), axis=0))
    assert np.all(np.sum(x[edges[:, 0]] != x[edges[:, 1]], axis=1) == 1)
    assert np.all(np.sum(x[faces[:, 0]] != x[faces[:, 3]], axis=1) == 2)
    degrees = np.bincount(edges.ravel(), minlength=n)
    return PanelGraph(x, edges, edge_sites, faces, degrees)


def regression_slope(background, effect):
    b, d = np.asarray(background), np.asarray(effect)
    centered = b - b.mean()
    denominator = np.dot(centered, centered)
    return float(np.dot(centered, d - d.mean()) / denominator) if denominator > 0 else np.nan


def face_states(graph, values, epsilon=0.0):
    a, b, c, d = np.asarray(values)[graph.faces].T
    effects = np.column_stack((b - a, d - c, c - a, d - b))
    signs = np.where(effects > epsilon, 1, np.where(effects < -epsilon, -1, 0))
    flip1, flip2 = signs[:, 0] * signs[:, 1] == -1, signs[:, 2] * signs[:, 3] == -1
    return flip1 | flip2, flip1 & flip2, np.all(np.abs(effects) > epsilon, axis=1)


def graph_metrics(graph, values, epsilon=0.0):
    """Exact expectation under uniform start and uniform improving-neighbor choice."""
    y = np.asarray(values, dtype=float)
    assert y.shape == (len(graph.sequences),) and np.all(np.isfinite(y))
    u, v = graph.edges.T
    delta = y[v] - y[u]
    src = np.concatenate((u[delta > epsilon], v[delta < -epsilon]))
    dst = np.concatenate((v[delta > epsilon], u[delta < -epsilon]))
    # Compressed outgoing adjacency; strict improvement makes this graph acyclic.
    order = np.argsort(src, kind="stable")
    neighbors = dst[order]
    counts = np.bincount(src, minlength=len(y))
    offsets = np.r_[0, np.cumsum(counts)]
    steps = np.zeros(len(y))
    for node in np.argsort(y)[::-1]:
        if counts[node]:
            steps[node] = 1.0 + steps[neighbors[offsets[node]:offsets[node + 1]]].mean()
    signed, reciprocal, resolved = face_states(graph, y, epsilon)
    return {
        "nodes": len(y), "edges": len(u), "complete_faces": len(graph.faces),
        "epsilon": float(epsilon),
        "sign_epistasis_fraction": float(signed.mean()),
        "reciprocal_sign_epistasis_fraction": float(reciprocal.mean()),
        "all_four_effects_resolved_fraction": float(resolved.mean()),
        "sign_epistasis_among_fully_resolved_faces": float(signed[resolved].mean()) if resolved.any() else np.nan,
        "reciprocal_sign_epistasis_among_fully_resolved_faces": float(reciprocal[resolved].mean()) if resolved.any() else np.nan,
        "fixed_orientation_diminishing_returns_slope": regression_slope(y[u], delta),
        "uniform_uphill_steps": float(steps.mean()),
        "panel_local_peak_fraction": float(np.mean(counts == 0)),
        "panel_local_peak_fraction_nonisolated": float(np.mean(counts[graph.degrees > 0] == 0)),
        "improving_oriented_edges": int(np.count_nonzero(np.abs(delta) > epsilon)),
        "neutral_edge_fraction": float(np.mean(np.abs(delta) <= epsilon)),
    }


QUANTILES = [0.0, .001, .01, .1, .5, .9, .95, .99, .999, 1.0]


def distribution_metrics(values):
    y = np.asarray(values, dtype=float).ravel()
    assert np.all(np.isfinite(y))
    return {
        "n": int(len(y)), "mean": float(y.mean()), "sd": float(y.std(ddof=0)),
        "minimum": float(y.min()), "maximum": float(y.max()),
        **{f"q{q:g}": float(value) for q, value in zip(QUANTILES, np.quantile(y, QUANTILES))},
    }


def dfe_metrics(graph, values):
    y = np.asarray(values)
    delta = y[graph.edges[:, 1]] - y[graph.edges[:, 0]]
    out = {f"signed_{key}": value for key, value in distribution_metrics(delta).items()}
    out.update({f"absolute_{key}": value for key, value in distribution_metrics(np.abs(delta)).items()})
    out["rms"] = float(np.sqrt(np.mean(delta ** 2)))
    return out


def mutation_slopes(graph, values, minimum_backgrounds=10):
    """Descriptive regressions for fixed mutation identity, with an explicit cutoff.

    Both response and predictor contain background fitness. Even within one
    mutation identity, negative slopes can include mathematical coupling and
    measurement noise. These slopes are not an independent mechanistic test.
    """
    y = np.asarray(values)
    u, v = graph.edges.T
    a = graph.sequences[u, graph.edge_sites]
    b = graph.sequences[v, graph.edge_sites]
    codes = graph.edge_sites * 400 + a * 20 + b
    rows = []
    for code in np.unique(codes):
        select = codes == code
        if select.sum() < minimum_backgrounds:
            continue
        rows.append({
            "site_index": int(code // 400), "allele_from": int((code % 400) // 20),
            "allele_to": int(code % 20), "backgrounds": int(select.sum()),
            "slope": regression_slope(y[u[select]], y[v[select]] - y[u[select]]),
        })
    return rows


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_landscape_analysis(project, analysis):
    project, analysis = Path(project), Path(analysis)
    source = project / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-184/outputs"
    tables = analysis / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    result = json.loads((source / "result.json").read_text())
    with np.load(source / "confirmation_inputs.npz") as data:
        sequences, replicates = data["sequences"].copy(), data["fitness"].copy()
    with np.load(source / "confirmation_predictions.npz") as predictions:
        np.testing.assert_array_equal(sequences, predictions["sequences"])
        seeds, synthetic = predictions["seeds"].copy(), predictions["values"].copy()
    with np.load(source / "inputs.npz") as data:
        calibration_x, calibration_y = data["sequences"].copy(), data["fitness"].copy()
    assert replicates.shape == (8192, 2) and synthetic.shape == (4, 8192, 2)
    assert len(set(map(tuple, sequences)) & set(map(tuple, calibration_x))) == 0
    np.testing.assert_array_equal(synthetic[:, :, 0], synthetic[:, :, 1])
    assert seeds.tolist() == result["confirmation_seeds"]
    graph = panel_graph(sequences)
    assert (len(graph.sequences), len(graph.edges), len(graph.faces)) == (8192, 15940, 566)
    epsilon = .001 * float(np.ptp(calibration_y.mean(axis=1)))
    series = {
        "measured_mean": replicates.mean(axis=1),
        "measured_rep1": replicates[:, 0],
        "measured_rep2": replicates[:, 1],
        **{f"seed_{seed}": synthetic[i].mean(axis=1) for i, seed in enumerate(seeds)},
    }
    kind = lambda key: "synthetic" if key.startswith("seed_") else "measured"
    lower = float(result["best_parameters"]["lower"])
    upper = lower + float(result["best_parameters"]["span"])
    fitness_rows = []
    for name, values in {**series, "measured_pooled_replicates": replicates.ravel()}.items():
        fitness_rows.append({"series": name, "kind": kind(name), **distribution_metrics(values),
            "fraction_below_model_lower_bound": float(np.mean(values < lower)),
            "fraction_within_model_bounds": float(np.mean((values >= lower) & (values <= upper))),
            "fraction_above_model_upper_bound": float(np.mean(values > upper)),
        })
    fitness = pd.DataFrame(fitness_rows)
    graph_rows, mutation_rows = [], []
    for name, values in series.items():
        for threshold_name, threshold in [("zero", 0.0), ("calibration_range_0.001", epsilon)]:
            graph_rows.append({"series": name, "kind": kind(name), "threshold": threshold_name,
                               **graph_metrics(graph, values, threshold)})
        mutation_rows.extend({"series": name, "kind": kind(name), **row} for row in mutation_slopes(graph, values))
    metrics = pd.DataFrame(graph_rows)
    mutations = pd.DataFrame(mutation_rows)
    dfe = pd.DataFrame([{"series": name, "kind": kind(name), **dfe_metrics(graph, values)} for name, values in series.items()])
    comparisons = []
    for metric in ["sign_epistasis_fraction", "reciprocal_sign_epistasis_fraction",
                   "fixed_orientation_diminishing_returns_slope", "uniform_uphill_steps",
                   "panel_local_peak_fraction"]:
        for threshold in ["zero", "calibration_range_0.001"]:
            rows = metrics.loc[metrics.threshold == threshold].set_index("series")
            target = float(rows.loc["measured_mean", metric])
            sims = rows.loc[[f"seed_{s}" for s in seeds], metric].to_numpy()
            comparisons.append({"metric": metric, "threshold": threshold,
                "measured_mean": target, "measured_rep1": float(rows.loc["measured_rep1", metric]),
                "measured_rep2": float(rows.loc["measured_rep2", metric]),
                "synthetic_min": float(sims.min()), "synthetic_median": float(np.median(sims)),
                "synthetic_max": float(sims.max()), "synthetic_median_minus_measured": float(np.median(sims) - target),
                "measured_inside_four_seed_range": bool(sims.min() <= target <= sims.max()),
            })
    comparison = pd.DataFrame(comparisons)
    verification = []
    recorded_metrics = ["nodes", "edges", "complete_faces", "sign_epistasis_fraction",
        "reciprocal_sign_epistasis_fraction", "fixed_orientation_diminishing_returns_slope",
        "uniform_uphill_steps", "panel_local_peak_fraction"]
    for name, recorded in [("measured_mean", result["unfitted_diagnostics"]["empirical"])] + [
        (f"seed_{seed}", row) for seed, row in zip(seeds, result["unfitted_diagnostics"]["synthetic"])]:
        recomputed = metrics.loc[(metrics.series == name) & (metrics.threshold == "zero")].iloc[0]
        for metric in recorded_metrics:
            error = abs(float(recomputed[metric]) - float(recorded[metric]))
            assert error < 1e-12, (name, metric, error)
            verification.append({"series": name, "metric": metric, "recorded": recorded[metric],
                                 "recomputed": float(recomputed[metric]), "absolute_error": error, "passed": True})
    u, v = graph.edges.T
    effects1, effects2 = replicates[v, 0] - replicates[u, 0], replicates[v, 1] - replicates[u, 1]
    noise = {
        "replicate_pearson_r": float(np.corrcoef(replicates.T)[0, 1]),
        "replicate_rmse": float(np.sqrt(np.mean((replicates[:, 1] - replicates[:, 0]) ** 2))),
        "replicate2_minus_replicate1_mean": float(np.mean(replicates[:, 1] - replicates[:, 0])),
        "edge_effect_pearson_r": float(np.corrcoef(effects1, effects2)[0, 1]),
        "edge_effect_rmse": float(np.sqrt(np.mean((effects1 - effects2) ** 2))),
        "edge_effect_sign_disagreement_zero": float(np.mean(effects1 * effects2 < 0)),
        "cross_replicate_pooled_slope_background1_effect2": regression_slope(replicates[u, 0], effects2),
        "cross_replicate_pooled_slope_background2_effect1": regression_slope(replicates[u, 1], effects1),
        "interpretation": "Two observed replicates expose sensitivity; they do not identify intrinsic fitness, independent errors or an assay noise model. Cross-replicate slopes are descriptive and not noise-corrected estimates.",
    }
    for threshold_name, threshold in [("zero", 0.0), ("calibration_range_0.001", epsilon)]:
        s1, r1, _ = face_states(graph, replicates[:, 0], threshold)
        s2, r2, _ = face_states(graph, replicates[:, 1], threshold)
        resolved = (np.abs(effects1) > threshold) & (np.abs(effects2) > threshold)
        noise[f"edge_resolved_in_both_fraction_{threshold_name}"] = float(resolved.mean())
        noise[f"edge_sign_disagreement_among_both_resolved_{threshold_name}"] = float(np.mean(effects1[resolved] * effects2[resolved] < 0))
        noise[f"face_sign_epistasis_class_disagreement_{threshold_name}"] = float(np.mean(s1 != s2))
        noise[f"face_reciprocal_class_disagreement_{threshold_name}"] = float(np.mean(r1 != r2))
    summary = {
        "scope": "Disjoint 8,192-genotype confirmation panel; four new generator seeds; same previously studied four-site TrpB dataset.",
        "native_fitness_units": True,
        "graph": {"nodes": len(sequences), "edges": len(graph.edges), "complete_faces": len(graph.faces),
            "fraction_of_q20_four_site_space": len(sequences) / 20**4,
            "mean_observed_degree": float(graph.degrees.mean()), "possible_degree": 76,
            "isolated_nodes": int(np.sum(graph.degrees == 0)),
            "isolated_node_fraction": float(np.mean(graph.degrees == 0)),
            "face_denominator": "All 566 observed complete two-site faces, including faces with unresolved effects; fully resolved conditional values are separately tabulated.",
            "interpretation": "Walk lengths and local peak fractions describe this sparse observed induced graph only, not whole-landscape searchability."},
        "positive_epsilon": epsilon,
        "epsilon_rule": "0.001 times the calibration panel replicate-mean fitness range; no confirmation-derived threshold.",
        "synthetic_noise_mode": "none; archived synthetic replicates are identical",
        "sigmoid_bounds": {"lower": lower, "upper": upper, "span": upper - lower},
        "fitness": fitness.to_dict(orient="records"),
        "unfitted_comparison": comparison.to_dict(orient="records"),
        "dfe": dfe.to_dict(orient="records"),
        "replicate_sensitivity": noise,
        "mutation_slopes": {name: {"mutations": int(len(rows)), "median": float(rows.slope.median()),
            "q10": float(rows.slope.quantile(.1)), "q90": float(rows.slope.quantile(.9))}
            for name, rows in mutations.groupby("series")},
        "diminishing_returns_caveat": "Fixed-direction slopes regress y_mutant - y_background on y_background. The shared subtraction creates mathematical coupling; allele-index orientation is arbitrary. Negative slopes, pooled or within mutation, alone do not establish biological diminishing returns.",
        "seed_interpretation": "All four confirmation seeds are shown; their min/max or median describes this ensemble, not confidence intervals or biological replication.",
        "source_sha256": {name: sha256(source / name) for name in ["confirmation_inputs.npz", "confirmation_predictions.npz", "inputs.npz", "result.json"]},
    }
    for label, frame in [("fitness", fitness), ("graph", metrics), ("dfe", dfe), ("mutations", mutations), ("comparison", comparison)]:
        frame.to_csv(tables / f"landscape_{label}.csv", index=False)
    save_json(tables / "landscape_summary.json", summary)
    save_json(tables / "landscape_verification.json", {"passed": True,
        "checks": len(verification), "maximum_absolute_error": max(row["absolute_error"] for row in verification),
        "description": "Independent vectorized graph construction and exact dynamic programming versus archived zero-threshold unfitted diagnostics; no model fit is rerun.",
        "rows": verification})
    return {"source": source, "result": result, "graph": graph, "series": series, "seeds": seeds,
            "replicates": replicates, "fitness": fitness, "metrics": metrics, "dfe": dfe,
            "mutations": mutations, "comparison": comparison, "summary": summary, "epsilon": epsilon}
