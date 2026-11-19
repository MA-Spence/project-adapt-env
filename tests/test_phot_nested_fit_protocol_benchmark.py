from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from project_adapt_env.smc_abc import ParameterSpec


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from phot_nested_fit_protocol_benchmark import (  # noqa: E402
    apply_affine_readout,
    build_equivalence_rows,
    epistasis_property_metrics,
    fit_positive_affine_readout,
    narrow_specs_around,
    polynomial_feature_matrix,
    prefixed_metrics,
    raw_to_unit_vector,
    unit_vector_to_raw,
)


def test_raw_unit_roundtrip_respects_float_and_integer_specs() -> None:
    specs = [
        ParameterSpec("x", "float", low=-2.0, high=2.0),
        ParameterSpec("k", "int", low=0, high=3),
    ]
    raw = {"x": 1.0, "k": 2}

    unit = raw_to_unit_vector(raw, specs)
    recovered = unit_vector_to_raw(unit, specs)

    assert recovered["x"] == pytest.approx(1.0)
    assert recovered["k"] == 2


def test_polynomial_feature_matrix_contains_bias_linear_square_and_interactions() -> None:
    x = np.asarray([[0.5, 0.25]], dtype=np.float64)

    design = polynomial_feature_matrix(x)

    assert design.shape == (1, 6)
    assert design[0].tolist() == pytest.approx([1.0, 0.5, 0.25, 0.25, 0.0625, 0.125])


def test_narrow_specs_around_keeps_bounds_and_integer_support() -> None:
    specs = [
        ParameterSpec("x", "float", low=0.0, high=10.0),
        ParameterSpec("k", "int", low=0, high=3),
    ]

    narrowed = narrow_specs_around(specs, {"x": 0.25, "k": 0}, width_fraction=0.2)

    assert narrowed[0].low == pytest.approx(0.0)
    assert 0.25 <= narrowed[0].high <= 2.0
    assert narrowed[1].low == 0.0
    assert narrowed[1].high >= 1.0


def test_epistasis_property_metrics_detect_pairwise_magnitude_and_diminishing_returns() -> None:
    # Negative interaction: both singles improve over background, but the double
    # is worse than the additive expectation.
    values = np.asarray([0.0, 1.0, 1.0, 1.25], dtype=np.float64)

    metrics = epistasis_property_metrics(values, dim=2)

    assert metrics["magnitude_epistasis_mean"] == pytest.approx(0.75)
    assert metrics["signed_pairwise_epistasis_mean"] == pytest.approx(-0.75)
    assert metrics["diminishing_returns_epistasis_fraction"] == pytest.approx(1.0)


def test_positive_affine_readout_rescales_pairwise_epistasis_magnitude() -> None:
    empirical = np.asarray([0.0, 1.0, 1.0, 1.25], dtype=np.float64)
    overlarge_model = empirical * 100.0

    readout = fit_positive_affine_readout([overlarge_model], [empirical])
    calibrated = apply_affine_readout(overlarge_model, readout)

    assert readout["readout_slope"] == pytest.approx(0.01)
    assert calibrated.tolist() == pytest.approx(empirical.tolist())
    assert epistasis_property_metrics(calibrated, dim=2)["magnitude_epistasis_mean"] == pytest.approx(0.75)


def test_prefixed_epistasis_metrics_are_available_for_raw_latent_reporting() -> None:
    metrics = epistasis_property_metrics(np.asarray([0.0, 2.0, 2.0, 3.0], dtype=np.float64), dim=2)
    raw = prefixed_metrics(metrics, "raw_")

    assert raw["raw_magnitude_epistasis_mean"] == pytest.approx(1.0)
    assert raw["raw_signed_pairwise_epistasis_mean"] == pytest.approx(-1.0)


def test_build_equivalence_rows_checks_family_and_best_any_abc_references() -> None:
    protocol = [
        {
            "graphlet_split": "test",
            "fit_method": "deterministic_summary_loss",
            "candidate_family": "empirical_module",
            "candidate_model": "cheap",
            "parameter_set_id": "cheap-1",
            "distance_to_empirical_primary": 0.51,
            "distance_component__accessible_path_fraction": 1.2,
        }
    ]
    comparator = pd.DataFrame(
        [
            {
                "graphlet_split": "test",
                "candidate_family": "empirical_module",
                "parameter_set_id": "abc-emp",
                "distance_to_empirical_primary": 0.50,
            },
            {
                "graphlet_split": "test",
                "candidate_family": "random_module",
                "parameter_set_id": "abc-best",
                "distance_to_empirical_primary": 0.49,
            },
        ]
    )
    config = {
        "equivalence": {
            "absolute_distance_margin": 0.05,
            "relative_distance_margin": 0.10,
            "max_component_distance_margin": 2.0,
        }
    }

    rows = build_equivalence_rows(protocol, comparator, config)

    assert {row["reference"] for row in rows} == {"same_family_full_abc", "best_any_full_abc"}
    assert all(row["operationally_equivalent"] for row in rows)
