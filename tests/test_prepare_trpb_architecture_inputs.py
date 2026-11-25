"""Frozen biological input identity, split independence and grouped-fold checks."""
import csv
import json
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_trpb_architecture_inputs as preparation
from adaptenv.calibration import CalibrationConfig, CalibrationData, _folds
from adaptenv.inference import RFAConfig

BUNDLE = ROOT / "data/processed/trpb-architecture-calibration-inputs"


def test_selection_is_geometry_only_disjoint_and_repeatable():
    x = np.indices((4, 4, 4, 4)).reshape(4, -1).T
    settings = dict(n_calibration=80, n_confirmation=70, selection_seed=11,
                    n_folds=5, inference_seed=12)
    first = preparation.choose_panels(x, **settings)
    second = preparation.choose_panels(x, **settings)
    expected = np.random.default_rng(11).permutation(len(x))
    np.testing.assert_array_equal(first["calibration"][0], expected[:80])
    np.testing.assert_array_equal(first["confirmation"][0], expected[80:150])
    assert not set(first["calibration"][0]) & set(first["confirmation"][0])
    for name in first:
        for a, b in zip(first[name], second[name]):
            np.testing.assert_array_equal(a, b)


def test_existing_input_bundle_is_never_overwritten(tmp_path):
    record = tmp_path / "manifest.json"
    record.write_text("preserve")
    with pytest.raises(ValueError, match="not overwritten"):
        preparation.prepare(tmp_path)
    assert record.read_text() == "preserve"


@pytest.mark.skipif(not (BUNDLE / "manifest.json").exists(), reason="DVC input bundle not materialized")
def test_frozen_panels_keep_actual_replicates_profile_mask_and_fold_rules():
    manifest = json.loads((BUNDLE / "manifest.json").read_text())
    for relative, expected in manifest["files"].items():
        assert preparation.sha(BUNDLE / relative) == expected
    assert manifest["selection"]["uses_fitness_values"] is False
    assert manifest["latent_coefficient_degrees_of_freedom"] == 2243
    matrix = np.load(BUNDLE / "source/profile/additive_matrix.npy")
    panels = {}
    for name in ("calibration", "confirmation"):
        with np.load(BUNDLE / f"{name}.npz", allow_pickle=False) as archive:
            panel = {key: archive[key] for key in archive.files}
        assert panel["sequences"].shape == (8192, 4)
        assert panel["fitness"].shape == (8192, 2)
        assert np.isfinite(panel["fitness"]).all()
        assert len(np.unique(panel["sequences"], axis=0)) == 8192
        assert not np.array_equal(panel["fitness"][:, 0], panel["fitness"][:, 1])
        np.testing.assert_array_equal(panel["additive_matrix"], matrix)
        data = CalibrationData(panel["sequences"], panel["fitness"])
        config = CalibrationConfig(inference=RFAConfig(n_folds=10, seed=20261122))
        np.testing.assert_array_equal(panel["folds"], _folds(data, config))
        assert min(np.sum(panel["folds"] != fold) for fold in np.unique(panel["folds"])) > 2243
        for record in preparation.coverage(panel["sequences"], panel["folds"]):
            assert record["observed_cells"] == 400
            assert record["minimum_cell_count"] >= 8
            assert record["minimum_training_cell_count"] >= 5
        panels[name] = panel
    assert not set(map(tuple, panels["calibration"]["sequences"])) & set(map(tuple, panels["confirmation"]["sequences"]))

    # Verify against the original author table, not merely the derived source CSV.
    selected = {}
    for panel in panels.values():
        for sequence, values in zip(panel["sequences"], panel["fitness"]):
            name = "".join(preparation.AA_ORDER[int(a)] for a in sequence)
            selected[name] = values
    seen = set()
    author_table = BUNDLE / "source/replicates/source/four-site_summary_AA_data.csv"
    with author_table.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["AAs"] in selected:
                np.testing.assert_allclose(selected[row["AAs"]],
                                           [float(row["mu_1-bg/max"]), float(row["mu_2-bg/max"])],
                                           rtol=0, atol=1e-14)
                seen.add(row["AAs"])
    assert seen == set(selected)
