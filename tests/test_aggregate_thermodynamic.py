from __future__ import annotations

import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from adaptenv import ThermodynamicConfig, ThermodynamicLandscape
from adaptenv.thermodynamic import sigmoid
from project_adapt_env.dense_descriptive import exact_summary, spectral_metrics
from project_adapt_env.summary_statistics import mutation_summary, population_edge_effects

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("aggregate_runner", ROOT / "scripts/aggregate_thermodynamic_calibration.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def test_flat_cube_is_one_plateau_not_sixteen_peaks():
    result = exact_summary(np.zeros(16), dim=4)
    assert result["n_local_maxima"] == 0
    assert result["n_weak_local_maxima"] == 16
    assert result["n_terminal_plateaus"] == result["basin_count"] == 1
    assert result["neutral_edge_fraction"] == 1
    assert result["largest_plateau_fraction"] == 1
    assert np.isnan(result["walsh_epistatic_fraction"])
    assert result["total_variance"] == 0


def test_resolution_applies_to_signs_peaks_accessibility_and_spectra():
    values = np.arange(16) * 1e-7
    unresolved = exact_summary(values, dim=4, resolution=1e-4)
    assert unresolved["n_local_maxima"] == 0
    assert unresolved["accessible_path_fraction"] == 0
    assert not unresolved["spectral_resolved"]
    resolved = exact_summary(values, dim=4, resolution=1e-12)
    assert resolved["n_local_maxima"] == 1
    assert resolved["accessible_path_fraction"] == 1


def test_walsh_reports_pairwise_and_high_order_separately():
    x = np.array(list(itertools.product([-1, 1], repeat=3)))
    values = 2 * x[:, 0] + x[:, 0] * x[:, 1] + 3 * np.prod(x, axis=1)
    result = spectral_metrics(values, dim=3)
    assert result["total_variance"] == pytest.approx(14)
    assert result["walsh_epistatic_fraction"] == pytest.approx(10/14)
    assert result["walsh_high_order_fraction"] == pytest.approx(9/14)
    assert sum(result[f"walsh_energy_order_{k}"] for k in range(1, 4)) == pytest.approx(14)


def test_reference_dfe_does_not_invent_reverse_mutations():
    states, values = [("A",), ("B",), ("C",)], [10., 5., 1.]
    effects, weights = population_edge_effects(states, values, background_weights=[1., 0., 0.])
    np.testing.assert_array_equal(effects, [-5., -9.])
    assert np.dot(weights, effects) == -7
    assert not np.any(effects > 0)


def test_same_mutation_background_slopes_and_replicate_correction():
    bits = np.array([[(n >> i) & 1 for i in range(4)] for n in range(16)])
    additive = bits @ np.array([0.1, 0.2, 0.3, 0.4])
    assert mutation_summary(additive, dim=4)["diminishing_returns_slope"] == pytest.approx(0, abs=1e-12)
    concave = additive - 0.25 * additive**2
    truth = mutation_summary(concave, dim=4)
    assert truth["diminishing_returns_slope"] < 0
    reps = concave[None, :] + np.random.default_rng(22).normal(0, 0.03, (80, 16))
    corrected = mutation_summary(reps.mean(axis=0), dim=4, replicates=reps)
    assert corrected["diminishing_returns_error_corrected"]
    assert corrected["diminishing_returns_slope"] == pytest.approx(truth["diminishing_returns_slope"], abs=0.025)


@pytest.mark.parametrize("length", [1, 2, 3, 4, 5, 16])
def test_generator_small_lengths_and_reference_occupancy(length):
    model = ThermodynamicLandscape(ThermodynamicConfig(length=length, seed=7))
    assert model.evaluate(model.reference) == pytest.approx(model.config.reference_occupancy)
    assert model.summary()["global_optimum"] == "not_certified"
    assert model.summary()["n_supports"] <= sum(__import__("math").comb(length, k) for k in (2, 3, 4) if k <= length)


def test_determinism_batching_and_incremental_mutations():
    cfg = ThermodynamicConfig(length=20, seed=28, batch_size=3, reference_rank=0.5)
    model = ThermodynamicLandscape(cfg)
    replay = ThermodynamicLandscape(cfg)
    rng = np.random.default_rng(80)
    seq = rng.integers(0, 20, (15, 20))
    values = model.evaluate_batch(seq)
    np.testing.assert_array_equal(values, replay.evaluate_batch(seq))
    np.testing.assert_array_equal(values[::-1], model.evaluate_batch(seq[::-1]))
    sites, alleles = np.repeat(np.arange(20), 20), np.tile(np.arange(20), 20)
    mutants = np.repeat(seq[:1], len(sites), axis=0)
    mutants[np.arange(len(sites)), sites] = alleles
    np.testing.assert_allclose(model.mutation_effects(seq[0], sites, alleles), model.evaluate_batch(mutants)-values[0], atol=3e-14)
    assert model.fingerprint == replay.fingerprint


def test_exact_uniform_mutation_normalization_on_two_sites():
    model = ThermodynamicLandscape(ThermodynamicConfig(length=2, seed=11))
    sequences = np.array(list(itertools.product(range(20), repeat=2)))
    a, g = model.latent_components(sequences)
    a, g = a.reshape(20, 20), g.reshape(20, 20)
    def energy(field):
        e0 = (field[:, None, :] - field[None, :, :])**2
        e1 = (field[:, :, None] - field[:, None, :])**2
        return (e0.sum() + e1.sum()) / (2*20*20*19)
    assert energy(a) == pytest.approx(1, abs=1e-12)
    assert energy(g) == pytest.approx(1, abs=1e-12)


def test_empirical_stability_pool_is_reproducible_and_explicit():
    model = ThermodynamicLandscape(ThermodynamicConfig(length=4), stability_effect_samples=[-3., -1., 0.1, 0.3])
    assert model.summary()["core_prior"] == "sampled_stabilizing_positive_effects"
    assert model.evaluate(model.reference) == pytest.approx(model.config.reference_occupancy)
    with pytest.raises(ValueError):
        model.evaluate([0.5, 0, 0, 0])
    with pytest.raises(ValueError):
        ThermodynamicConfig(steepness=float("nan"))


def test_measured_sampling_has_disjoint_nodes_and_no_imputation(tmp_path):
    states = list(itertools.product("ABCDEF", repeat=3))
    frame = pd.DataFrame(states, columns=["pos1", "pos2", "pos3"])
    frame["fitness"] = np.arange(len(frame))
    path = tmp_path / "variants.csv"; frame.to_csv(path, index=False)
    cubes, audit = runner.measured_cubes(path, dimension=2, count=12, seed=14)
    keys = [tuple(n) for a in audit for n in a["node_keys"]]
    assert len(keys) == len(set(keys)) == 48
    assert all(len(c) == 4 for c in cubes)
    np.testing.assert_array_equal(cubes, runner.measured_cubes(path, dimension=2, count=12, seed=14)[0])


def test_parquet_and_csv_sample_identical_measured_cubes(tmp_path):
    pytest.importorskip("pyarrow")
    states = list(itertools.product("ABCD", repeat=3))
    frame = pd.DataFrame(states, columns=["pos1", "pos2", "pos3"])
    frame["fitness"] = np.arange(len(frame), dtype=float) / 4
    csv, parquet = tmp_path / "variants.csv", tmp_path / "variants.parquet"
    frame.to_csv(csv, index=False); frame.to_parquet(parquet, index=False)
    a, aa = runner.measured_cubes(csv, dimension=2, count=6, seed=15)
    b, bb = runner.measured_cubes(parquet, dimension=2, count=6, seed=15)
    np.testing.assert_array_equal(a, b)
    assert aa == bb


def test_nonfinite_measurements_require_declared_exclusion(tmp_path):
    states = list(itertools.product("ABCD", repeat=2))
    frame = pd.DataFrame(states, columns=["pos1", "pos2"])
    frame["fitness"] = np.arange(len(frame), dtype=float)
    frame.loc[0, "fitness"] = np.nan
    path = tmp_path / "variants.csv"; frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="nonfinite"):
        runner.measured_cubes(path, dimension=2, count=2, seed=6)
    cubes, audit = runner.measured_cubes(path, dimension=2, count=2, seed=6, missing_fitness="exclude")
    assert all(np.isfinite(c).all() for c in cubes)
    assert all(a["excluded_nonfinite_rows"] == 1 for a in audit)
    assert not any(tuple(n) == ("A", "A") for a in audit for n in a["node_keys"])


def test_beneficial_and_deleterious_slopes_do_not_cancel():
    # All logits are on the concave branch; positive and negative mutations
    # have opposite slope signs even though they share the same nonlinear map.
    nodes = np.arange(16)
    logits = 5 + sum(effect * ((nodes >> bit) & 1)
                     for bit, effect in enumerate([0.7, -0.7, 0.5, -0.5]))
    result = mutation_summary(sigmoid(logits), dim=4)
    assert result["beneficial_slope_n_mutations"] == 2
    assert result["deleterious_slope_n_mutations"] == 2
    assert result["beneficial_diminishing_returns_slope"] < 0
    assert result["deleterious_background_slope"] > 0


def test_joint_distance_detects_changed_covariance_with_identical_marginals():
    def row(x, y):
        return dict(zip(runner.FEATURES, [x, y] + [0.]*(len(runner.FEATURES)-2)))
    data = [row(-1, -1), row(1, 1)]
    mismatched = [row(-1, 1), row(1, -1)]
    distance = runner.JointDistance(data)
    assert distance(data, data) == pytest.approx(0)
    assert distance(data, mismatched) > 0.1
    unresolved = [dict(data[0], walsh_epistatic_fraction=np.nan), data[1]]
    assert distance(data, unresolved) > 0


def test_smoke_writes_strict_manifest_and_never_claims_science(tmp_path):
    import yaml
    path = next((ROOT / "experiments").glob("*EXP-065*/config.yaml"))
    result = runner.run(yaml.safe_load(path.read_text()), tmp_path, smoke=True)
    assert result["status"] == "SMOKE_ONLY"
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["input_sha256"] == {}
    assert manifest["source_sha256"]
    production = yaml.safe_load(path.read_text())
    assert not set(manifest["config"]["fit"]["confirmation_seeds"]) & set(production["fit"]["confirmation_seeds"])
    assert set(o["family"] for o in result["outcomes"]) == {"additive", "bare_gnk", "composite"}
    assert all(np.isfinite(o["confirmation_joint_loss"]) for o in result["outcomes"])
