from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "observable_ruggedness_coverage_probe_exp064.py"
SPEC = importlib.util.spec_from_file_location("observable_ruggedness_coverage_probe_exp064", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def load_exp064_config() -> dict:
    path = (
        PROJECT_ROOT
        / "experiments"
        / "2026-07-24_EXP-064_observable-ruggedness-coverage-probe"
        / "config.yaml"
    )
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_design_key_excludes_amplitude_and_operating_point_is_explicit() -> None:
    first = runner.design_key(
        template="gb1",
        support_mode="structured",
        representation="binary_r8",
        beta=2.0,
        p0=0.5,
        seed=6401,
    )
    repeated = runner.design_key(
        template="gb1",
        support_mode="structured",
        representation="binary_r8",
        beta=2.0,
        p0=0.5,
        seed=6401,
    )
    changed = runner.design_key(
        template="gb1",
        support_mode="structured",
        representation="binary_r8",
        beta=4.0,
        p0=0.5,
        seed=6401,
    )
    assert first == repeated
    assert first != changed


def test_observable_reference_occupancy_and_amplitude_pairing() -> None:
    additive = np.asarray([3.0, 2.0])
    gnk = np.asarray([4.0, -2.0])
    at_zero = runner.observable_fitness(
        additive,
        gnk,
        additive_reference=3.0,
        sigma_a=2.0,
        sigma_g=1.5,
        beta=4.0,
        p0=0.7,
        amplitude=0.0,
        output_scale=10.0,
    )
    changed_gnk = runner.observable_fitness(
        additive,
        100.0 * gnk,
        additive_reference=3.0,
        sigma_a=2.0,
        sigma_g=1.5,
        beta=4.0,
        p0=0.7,
        amplitude=0.0,
        output_scale=10.0,
    )
    assert at_zero[0] == pytest.approx(7.0)
    np.testing.assert_array_equal(at_zero, changed_gnk)


def test_frozen_stage4_descriptor_path_recovers_additive_control() -> None:
    values = np.asarray([node.bit_count() for node in range(16)], dtype=float)
    descriptors = runner.descriptors_for_values(values, dimension=4)
    assert descriptors["local_maxima_density"] == pytest.approx(1.0 / 16.0)
    assert descriptors["accessible_path_fraction"] == pytest.approx(1.0)
    assert descriptors["walsh_high_order_fraction"] == pytest.approx(0.0)
    assert descriptors["sign_epistasis_fraction"] == pytest.approx(0.0)


def test_empirical_band_ingestion_reproduces_published_res050_tables() -> None:
    config = load_exp064_config()
    config["_runtime"] = {"quick": False}
    bands, checks, rugged = runner.resolve_empirical_bands(config, PROJECT_ROOT)
    assert len(bands) == 10 * 8
    assert len(rugged) == 3
    assert all(bool(row["pass"]) for row in checks)
    assert {row["template"] for row in bands} == set(
        config["data"]["tier1_dataset_ids"]
    )


def test_routed_sequences_preserve_reference_and_mutation_count() -> None:
    reference = np.arange(12, dtype=np.int64) % 20
    sequences = np.repeat(reference[None, :], 3, axis=0)
    sequences[1, 2] = (sequences[1, 2] + 3) % 20
    sequences[2, [1, 7]] = (sequences[2, [1, 7]] + np.asarray([2, 5])) % 20
    permutation = np.asarray([4, 9, 0, 7, 2, 11, 5, 3, 10, 8, 1, 6])
    routed = runner._route_sequences(sequences, reference, permutation)
    np.testing.assert_array_equal(routed[0], reference)
    np.testing.assert_array_equal(
        np.sum(sequences != reference, axis=1),
        np.sum(routed != reference, axis=1),
    )


def test_quick_grid_retains_required_primary_and_specificity_controls() -> None:
    config = load_exp064_config()
    quick = config["quick"]["generator"]
    assert {row["id"] for row in quick["representations"]} >= {
        "binary_r8",
        "binary_r1",
        "dense_q20_r8",
    }
    assert set(quick["support_modes"]) == {"structured", "shuffled", "random"}
    assert 0.0 in quick["amplitudes"]
    assert any(float(value) > 0.0 for value in quick["amplitudes"])
