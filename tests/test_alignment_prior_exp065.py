"""Controls for the matched prior intervention, not empirical acceptance tests."""
import copy
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import alignment_prior_exp065 as repeat

EXPERIMENT = ROOT / "experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration"
INPUTS = ROOT / "data/processed/interaction-strength-exp065-inputs"


def test_prior_intervention_preserves_reference_gnk_and_public_readout():
    config = dict(reference_rank=.5, representation_rank=8,
                  order_weights=(.35, .45, .20))
    frequencies = np.random.default_rng(911).dirichlet(np.ones(20), size=4)
    matrix = np.log(20 * frequencies)
    nodes = np.random.default_rng(913).integers(20, size=(128, 4))
    for seed in [57117694, 57117794]:
        control = repeat.sweep.ThermodynamicLandscape(
            repeat.sweep.ThermodynamicConfig(length=4, seed=seed, **config)
        )
        changed = repeat.build_model(seed, config, matrix)
        np.testing.assert_array_equal(changed.reference, control.reference)
        _, old_g = control.latent_components(nodes)
        new_a, new_g = changed.latent_components(nodes)
        np.testing.assert_allclose(new_g, old_g, atol=1e-12, rtol=1e-12)
        expected_a = (matrix[np.arange(4), nodes].sum(axis=1)
                      - matrix[np.arange(4), changed.reference].sum()) / changed.sigma_a
        np.testing.assert_allclose(new_a, expected_a, atol=1e-12, rtol=1e-12)
        parameters = np.array([-.1, np.log(1.2), -.7, np.log(.6)])
        for ratio in [0, .5, 1, 2, 4, 8]:
            public = repeat.build_model(
                seed, config, matrix, interaction_ratio=ratio,
                reference_occupancy=float(repeat.expit(parameters[2])),
                steepness=float(np.exp(parameters[3])), fitness_min=parameters[0],
                fitness_max=parameters[0] + np.exp(parameters[1]),
            )
            expected = repeat.sweep.predict(parameters, new_a, new_g, ratio)
            np.testing.assert_allclose(public.evaluate_batch(nodes), expected,
                                       atol=1e-12, rtol=1e-12)


def test_panel_uses_recorded_seed_and_geometry_and_replays_frozen_g():
    index = pd.read_csv(INPUTS / "TrpB/latent_calibration_index.csv").iloc[:16].reset_index(drop=True)
    with np.load(INPUTS / "TrpB/latent_calibration.npz") as values:
        frozen = {key: values[key][:16].copy() for key in values.files}
    config = yaml.safe_load((EXPERIMENT / "interaction_strength.yaml").read_text())["model"]
    matrix = np.log(20 * np.random.default_rng(71).dirichlet(np.ones(20), size=4))
    new_index, additive, interactions, audit = repeat.regenerate_panel(index, frozen, config, matrix)
    pd.testing.assert_frame_equal(new_index.drop(columns="fingerprint"), index.drop(columns="fingerprint"))
    np.testing.assert_allclose(interactions, frozen["g"], rtol=1e-10, atol=1e-10)
    assert not np.allclose(additive, frozen["a"])
    assert len(audit) == 2
    assert all(item["maximum_normalized_g_difference"] < 1e-10 for item in audit)
    corrupted = {key: value.copy() for key, value in frozen.items()}
    corrupted["g"][0, 0] += .1
    with pytest.raises(AssertionError):
        repeat.regenerate_panel(index, corrupted, config, matrix)


def test_endpoint_refits_use_original_optimizer_seeds_and_feed_intermediate_candidates(tmp_path, monkeypatch):
    config = yaml.safe_load((EXPERIMENT / "interaction_strength.yaml").read_text())
    original = yaml.safe_load((INPUTS / "config.yaml").read_text())
    calls, sweep_calls = [], []
    calibration = np.arange(32).reshape(2, 16) / 32
    banks = {mode: object() for mode in config["observation_modes"]["TrpB"]}

    def endpoint(arm, a, g, target, bank, origin, scale, fitting, seed, quick, floor):
        calls.append((arm, seed, bank))
        assert a is calibration and g is calibration and target is calibration
        return dict(parameters=[float(seed), 0., 0., 0.], evaluations=1952), pd.DataFrame({"value": [seed]})

    def intermediate(ratio, mode, a, g, target, bank, origin, scale, floor,
                     fitting, seed, anchors, quick):
        assert a is calibration and g is calibration and target is calibration
        assert set(anchors) == {0, 4}
        for r in [0, 4]:
            expected = (repeat.production.seed_for(repeat.DATASET)
                        + original["fit"]["seed"]
                        + config["observation_modes"]["TrpB"].index(mode)*10
                        + original["model"]["arms"].index(f"r{r}_sigmoid"))
            assert anchors[r]["parameters"][0] == expected
        sweep_calls.append((ratio, mode, seed))
        label = "frozen_RUN172_anchor" if ratio in [0, 4] else "RUN172_r0_candidate"
        return dict(selected_source=label, objective={"total": 0.}), pd.DataFrame({"candidate": [label]})

    monkeypatch.setattr(repeat.production, "fit_map", endpoint)
    monkeypatch.setattr(repeat.sweep, "fit_readout", intermediate)
    results = repeat.fit_settings(calibration, calibration, calibration, banks,
                                  0., 1., .05, config, original, tmp_path)
    assert len(calls) == 6  # Two refitted endpoints per observation mode.
    assert len(sweep_calls) == len(results) == 18
    for result in results:
        assert not result["selected_source"].startswith("RUN172")
    for mode_index, mode in enumerate(config["observation_modes"]["TrpB"]):
        for ratio_index, ratio in enumerate(config["interaction_ratios"]):
            expected = (repeat.production.seed_for(repeat.DATASET) + config["fit"]["seed"]
                        + 1000 + mode_index*100 + ratio_index)
            assert (ratio, mode, expected) in sweep_calls


def test_frozen_kernels_and_matched_design_cannot_silently_drift(tmp_path):
    historical = EXPERIMENT / "interaction_strength_source.json"
    metadata = json.loads(historical.read_text())
    current = tmp_path / "source.json"
    current.write_text(json.dumps({"files": {
        name: metadata["files"][name] for name in repeat.REUSED_SOURCE_FILES
    }}))
    assert set(repeat.verify_sources(current, historical)) == set(repeat.REUSED_SOURCE_FILES)
    original = yaml.safe_load((EXPERIMENT / "interaction_strength.yaml").read_text())
    changed = copy.deepcopy(original)
    changed["proteins"] = ["TrpB"]
    repeat.validate_config(changed, original)
    changed["sampling"]["synthetic_cubes_per_seed"] = 7
    with pytest.raises(ValueError, match="sampling"):
        repeat.validate_config(changed, original)


def test_prior_requires_frozen_files_and_canonical_residue_order(tmp_path):
    matrix_file, metadata_file = tmp_path / "matrix.npy", tmp_path / "metadata.json"
    matrix = np.arange(80, dtype=float).reshape(4, 20)
    np.save(matrix_file, matrix)
    metadata = {"alphabet": list(repeat.AA_ORDER),
                "assay_positions_1based": repeat.ASSAY_POSITIONS,
                "matrix_rows_correspond_to_assay_positions_1based": repeat.ASSAY_POSITIONS}
    metadata_file.write_text(json.dumps(metadata))
    config = {"additive_prior": {"matrix_sha256": repeat.sweep.sha(matrix_file),
                                 "metadata_sha256": repeat.sweep.sha(metadata_file)}}
    loaded, _ = repeat.load_prior(matrix_file, metadata_file, config)
    np.testing.assert_array_equal(loaded, matrix)
    metadata["alphabet"] = list(repeat.AA_ORDER)[::-1]
    metadata_file.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="checksum"):
        repeat.load_prior(matrix_file, metadata_file, config)
    config["additive_prior"]["metadata_sha256"] = repeat.sweep.sha(metadata_file)
    with pytest.raises(ValueError, match="AA_ORDER"):
        repeat.load_prior(matrix_file, metadata_file, config)
    metadata["alphabet"] = list(repeat.AA_ORDER)
    metadata["matrix_rows_correspond_to_assay_positions_1based"] = [184, 183, 227, 228]
    metadata_file.write_text(json.dumps(metadata))
    config["additive_prior"]["metadata_sha256"] = repeat.sweep.sha(metadata_file)
    with pytest.raises(ValueError, match="ordered TrpB positions"):
        repeat.load_prior(matrix_file, metadata_file, config)
