"""Protocol-isolation checks, not empirical fitness-model acceptance gates."""
import ast
import copy
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import pairwise_calibration_exp065 as repeat

EXP = ROOT / "experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration"
INPUTS = ROOT / "data/processed/interaction-strength-exp065-inputs"
CONTROL = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178/outputs"


def new_config():
    config = yaml.safe_load((CONTROL / "config.yaml").read_text())
    config.update(protocol="trpb_full_pairwise_repeat_v1", interpretation="Changed interaction ensemble",
                  model=copy.deepcopy(repeat.EXPECTED_MODEL), profile_control_run="RUN-178",
                  profile_control_bundle=str(CONTROL.relative_to(ROOT)),
                  profile_control_manifest_sha256=repeat.sweep.sha(CONTROL / "manifest.json"))
    return config


def calibration_inputs():
    index = pd.read_csv(INPUTS / "TrpB/latent_calibration_index.csv")
    control_index = pd.read_csv(CONTROL / "latent_calibration_index.csv")
    with np.load(INPUTS / "TrpB/latent_calibration.npz") as values:
        original = {key: values[key] for key in values.files}
    with np.load(CONTROL / "latent_calibration.npz") as values:
        control = {key: values[key] for key in values.files}
    audits = json.loads((CONTROL / "calibration_prior_audit.json").read_text())
    matrix = np.load(CONTROL / "additive_matrix.npy")
    return index, original, matrix, control, control_index, audits


def test_only_declared_model_and_metadata_changes_are_allowed():
    config = new_config()
    historical = yaml.safe_load((EXP / "interaction_strength.yaml").read_text())
    profile = yaml.safe_load((CONTROL / "config.yaml").read_text())
    repeat.validate_config(config, historical)
    repeat.validate_profile_config(config, profile)
    for key in ["interaction_ratios", "seed_design", "fit", "sampling", "measurement"]:
        changed = copy.deepcopy(config)
        changed[key] = None
        with pytest.raises(ValueError, match=key):
            repeat.validate_config(changed, historical)
    for key, value in [("representation_rank", 8), ("target_degrees", [2.0]),
                       ("orders", [2, 3]), ("reference_rank", 0.6),
                       ("additive_prior", "synthetic")]:
        changed = copy.deepcopy(config)
        changed["model"][key] = value
        with pytest.raises(ValueError, match="full_pairwise"):
            repeat.validate_config(changed, historical)
    for key in ["additive_prior", "required_worker_environment", "endpoint_calibration_run",
                "source_run", "observation_modes", "input_bundle"]:
        changed = copy.deepcopy(config)
        changed[key] = None
        with pytest.raises(ValueError, match="RUN-178"):
            repeat.validate_profile_config(changed, profile)
    changed = copy.deepcopy(config)
    changed["undeclared_setting"] = 1
    with pytest.raises(ValueError, match="undeclared_setting"):
        repeat.validate_profile_config(changed, profile)


def test_actual_control_manifest_covers_all_frozen_outputs():
    manifest = repeat.verify_profile_control(CONTROL, new_config())
    assert len(manifest["files"]) == 99
    assert manifest["input_manifest_sha256"] == repeat.sweep.sha(INPUTS / "manifest.json")


def test_control_manifest_and_every_payload_are_checked(tmp_path):
    files = {}
    for i in range(99):
        path = tmp_path / f"payload{i}"
        path.write_text(str(i))
        files[path.name] = repeat.sweep.sha(path)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(dict(mode="MATCHED_PRIOR_REPEAT",
                                      public_model="thermodynamic_gnk_v2", files=files)))
    config = dict(profile_control_run="RUN-178",
                  profile_control_manifest_sha256=repeat.sweep.sha(manifest))
    repeat.verify_profile_control(tmp_path, config)
    (tmp_path / "payload98").write_text("changed")
    with pytest.raises(ValueError, match="payload98"):
        repeat.verify_profile_control(tmp_path, config)
    config["profile_control_manifest_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="manifest checksum"):
        repeat.verify_profile_control(tmp_path, config)


def test_calibration_keeps_profile_reference_geometry_and_changes_interactions():
    index, original, matrix, control, control_index, audits = calibration_inputs()
    updated, a, g, result = repeat.regenerate_panel(
        index, original, repeat.EXPECTED_MODEL, matrix, control, control_index, audits)
    pd.testing.assert_frame_equal(updated.drop(columns="fingerprint"),
                                  control_index.drop(columns="fingerprint"))
    np.testing.assert_array_equal(a, control["a"])
    assert np.max(np.abs(g-control["g"])) > 0.1
    assert len(result) == 8
    assert all(row["n_pair_supports"] == 6 for row in result)
    assert all(row["maximum_normalized_a_difference"] == 0 for row in result)
    np.testing.assert_allclose([row["normalized_interaction_mutation_energy"] for row in result], 1.)
    # r=0 is an exact negative control before the unchanged fitting algorithm.
    for setting in json.loads((CONTROL / "selected_parameters.json").read_text()):
        if setting["ratio"] == 0:
            np.testing.assert_array_equal(
                repeat.sweep.predict(setting["parameters"], a, g, 0),
                repeat.sweep.predict(setting["parameters"], control["a"], control["g"], 0))


@pytest.mark.parametrize("corruption", ["additive", "reference", "geometry", "index"])
def test_calibration_rejects_unintended_changes(corruption):
    index, original, matrix, control, control_index, audits = calibration_inputs()
    if corruption == "additive":
        control["a"][0, 0] += 0.01
    elif corruption == "reference":
        audits[0]["reference"][0] = (audits[0]["reference"][0] + 1) % 20
    elif corruption == "geometry":
        control["genotypes"][0, 0, 0] = (control["genotypes"][0, 0, 0] + 1) % 20
    else:
        control_index.loc[0, "cube"] += 100
    with pytest.raises(AssertionError):
        repeat.regenerate_panel(index, original, repeat.EXPECTED_MODEL, matrix,
                                control, control_index, audits)


def test_scientific_kernel_exception_is_only_known_whitespace(tmp_path, monkeypatch):
    historical = EXP / "interaction_strength_source.json"
    old_graph = (ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration"
                 / "RUN-178/source/scripts/landscape_graph.py")
    for relative in repeat.REUSED_SOURCE_FILES:
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, destination)
    archived = tmp_path / old_graph.relative_to(ROOT)
    archived.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(old_graph, archived)
    monkeypatch.setattr(repeat, "ROOT", tmp_path)
    source = tmp_path / "source.json"

    def freeze():
        source.write_text(json.dumps({"files": {
            relative: repeat.sweep.sha(tmp_path / relative)
            for relative in repeat.REUSED_SOURCE_FILES}}))

    freeze()
    audit = repeat.verify_sources(source, historical)
    assert audit["scripts/landscape_graph.py"]["equivalence"] == "trailing_whitespace_only_and_ast_identical"
    graph = tmp_path / "scripts/landscape_graph.py"
    graph.write_text(graph.read_text() + "\nINJECTED_SETTING = 1\n")
    freeze()
    with pytest.raises(ValueError, match="semantically"):
        repeat.verify_sources(source, historical)


def test_fitting_algorithm_ast_is_identical_except_provenance_labels():
    labels = {
        "RUN172_endpoint_prior_only_refit": "RUN172_endpoint_profile_pairwise_refit",
        "matched_alignment_endpoint_refit": "matched_pairwise_endpoint_refit",
        "matched_alignment_r0_endpoint_candidate": "matched_pairwise_r0_endpoint_candidate",
        "matched_alignment_r4_endpoint_candidate": "matched_pairwise_r4_endpoint_candidate",
    }

    class CanonicalLabels(ast.NodeTransformer):
        def visit_Constant(self, node):
            if isinstance(node.value, str) and node.value in labels:
                node.value = labels[node.value]
            return node

    def fitting_ast(path):
        tree = ast.parse(path.read_text())
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                        and node.name == "fit_settings")
        return ast.dump(CanonicalLabels().visit(function))

    assert fitting_ast(ROOT / "scripts/pairwise_calibration_exp065.py") == fitting_ast(
        ROOT / "scripts/alignment_prior_exp065.py")
