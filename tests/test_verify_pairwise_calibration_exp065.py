"""Negative controls for the independent saved-output validator."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pairwise_validator", ROOT/"scripts/verify_pairwise_calibration_exp065.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)
CONTROL = ROOT/"data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178/outputs"


def test_array_comparison_detects_shift_and_nan_mask_change():
    baseline = np.array([0., 1., np.nan])
    assert validator.check_array(baseline, baseline, "fixture")["array_equal"]
    for changed in [np.array([.01, 1., np.nan]), np.array([0., 1., 0.])]:
        with pytest.raises(validator.ValidationError, match="Numeric mismatch"):
            validator.check_array(changed, baseline, "fixture")


def test_summary_rejects_duplicate_and_missing_design_cells():
    frame = pd.read_csv(CONTROL/"latent_graph_summary.csv")
    audit = pd.read_csv(CONTROL/"seed_audit.csv")
    seeds = list(audit[["seed", "cohort"]].itertuples(index=False, name=None))
    assert validator.coverage(frame, seeds, latent=True)["rows"] == 144
    with pytest.raises(validator.ValidationError, match="coverage"):
        validator.coverage(frame.iloc[:-1], seeds, latent=True)
    with pytest.raises(validator.ValidationError, match="Duplicate"):
        validator.coverage(pd.concat([frame, frame.iloc[:1]]), seeds, latent=True)


def test_configuration_rejects_scientific_drift():
    old = yaml.safe_load((CONTROL/"config.yaml").read_text())
    new = deepcopy(old)
    new.update(model=deepcopy(validator.MODEL), profile_control_run="RUN-178")
    validator.validate_config(new, old)
    new["fit"]["maxiter"] += 1
    with pytest.raises(validator.ValidationError, match="Undeclared configuration"):
        validator.validate_config(new, old)


def test_manifest_rejects_corruption_and_root_escape(tmp_path):
    payload = tmp_path/"payload.txt"
    payload.write_text("frozen")
    files = {payload.name: validator.sha(payload)}
    validator.check_files(tmp_path, files)
    payload.write_text("changed")
    with pytest.raises(validator.ValidationError, match="Checksum"):
        validator.check_files(tmp_path, files)
    with pytest.raises(validator.ValidationError, match="escapes"):
        validator.check_files(tmp_path, {"../outside": "0"*64})


def test_production_validation_rejects_control_as_candidate():
    with pytest.raises(validator.ValidationError, match="MATCHED_PAIRWISE_REPEAT"):
        validator.validate_run(CONTROL, CONTROL, {})


def test_missing_old_source_is_explicit_and_not_a_current_source_fallback(tmp_path):
    result = validator.optional_archival_sources({"files": {"old.py": "0"*64}}, [tmp_path])
    assert result["available"] is False and result["required"] is False
    with pytest.raises(validator.ValidationError, match="No complete payload"):
        validator.manifested_root({"files": {"new.py": "0"*64}}, [tmp_path])


def test_r0_numeric_comparison_retains_objectives_but_ignores_provenance():
    old = {"parameters": [1., 2.], "objective": {"total": .25}, "selected_source": "old"}
    new = deepcopy(old)
    new["selected_source"] = "new"
    validator.compare_numeric_json(new, old, "r0")
    new["objective"]["total"] = .3
    with pytest.raises(validator.ValidationError, match="Numeric mismatch"):
        validator.compare_numeric_json(new, old, "r0")
