from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from project_adapt_env.smc_abc import ParameterSpec


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from phot_trait_space_bounded_gate_fast_fit import (  # noqa: E402
    apply_observation_readout,
    candidate_updates,
    center_raw_from_updates,
    fit_observation_readout,
    raw_latent_objective,
    run_softmin_capability_probe,
    softmin_latent_values,
    stage_configs,
)


def test_center_raw_from_updates_uses_fitted_values_and_geometric_bound_default() -> None:
    specs = [
        ParameterSpec("epistasis_strength", "float", low=0.0, high=0.1),
        ParameterSpec("n_epistatic_terms", "int", low=0, high=500),
        ParameterSpec("trait_space_epistasis_bound", "float", low=0.01, high=1.0),
    ]

    raw = center_raw_from_updates(
        {"epistasis_strength": 0.042, "n_epistatic_terms": 336},
        specs,
    )

    assert raw["epistasis_strength"] == pytest.approx(0.042)
    assert raw["n_epistatic_terms"] == 336
    assert raw["trait_space_epistasis_bound"] == pytest.approx(0.1)


def test_candidate_updates_force_bounded_trait_space_stability_route() -> None:
    config = {"bounded_trait_space_fast_fit": {"bound_mode": "clip"}}
    updates = candidate_updates(
        base_updates={"generic_epistasis_target": "score", "epistasis_strength": 0.03},
        raw={"trait_space_epistasis_bound": 0.25},
        config=config,
    )

    assert updates["generic_epistasis_target"] == "stability"
    assert updates["trait_space_epistasis_bound"] == pytest.approx(0.25)
    assert updates["trait_space_epistasis_bound_mode"] == "clip"
    assert updates["epistasis_strength"] == pytest.approx(0.03)


def test_candidate_updates_add_global_readout_and_ratio_strength() -> None:
    config = {
        "bounded_trait_space_fast_fit": {
            "bound_mode": "tanh",
            "global_readout": {
                "enabled": True,
                "source": "stability_margin",
                "transform": "sigmoid",
                "combine_mode": "sum",
                "fitness_output_scale": 1.0,
            },
            "specific_epistasis_variance_ratio": {
                "strength_scale": 0.04,
                "reference_n_terms": 220,
                "max_strength": 0.06,
                "empirical_pairwise_multiplier": 0.5,
            },
        }
    }
    updates = candidate_updates(
        base_updates={"stability_scale": 1.0, "n_epistatic_terms": 220},
        raw={
            "global_readout_midpoint": -1.0,
            "global_readout_slope": 2.5,
            "global_readout_weight": 1.2,
            "specific_epistasis_variance_ratio": 0.25,
            "trait_space_epistasis_bound": 0.1,
        },
        config=config,
    )

    assert updates["observed_fitness_combine_mode"] == "sum"
    assert updates["observed_fitness_terms"] == [
        {
            "source": "stability_margin",
            "transform": "sigmoid",
            "midpoint": pytest.approx(-1.0),
            "slope": pytest.approx(2.5),
            "weight": pytest.approx(1.2),
            "exponent": pytest.approx(1.0),
        }
    ]
    assert updates["fitness_output_scale"] == pytest.approx(1.0)
    assert updates["epistasis_strength"] == pytest.approx(0.01)
    assert updates["empirical_pairwise_strength"] == pytest.approx(0.005)
    assert "specific_epistasis_variance_ratio" not in updates


def test_candidate_updates_softmin_latent_pops_generator_parameters_and_disables_j() -> None:
    config = {
        "bounded_trait_space_fast_fit": {
            "bound_mode": "tanh",
            "softmin_latent": {"enabled": True, "disable_specific_epistasis": True},
        }
    }

    updates = candidate_updates(
        base_updates={"n_epistatic_terms": 220, "epistasis_strength": 0.03, "empirical_pairwise_strength": 0.02},
        raw={
            "trait_space_epistasis_bound": 0.1,
            "softmin_tau": 0.5,
            "softmin_n_traits": 5,
            "softmin_auxiliary_scale": 1.25,
            "softmin_auxiliary_bias": -0.2,
            "softmin_auxiliary_correlation": 0.4,
            "softmin_auxiliary_order": 3,
            "softmin_output_scale": 1.1,
            "softmin_structured_residual_scale": 0.25,
        },
        config=config,
    )

    assert updates["n_epistatic_terms"] == 0
    assert updates["epistasis_strength"] == pytest.approx(0.0)
    assert updates["empirical_pairwise_strength"] == pytest.approx(0.0)
    assert "softmin_tau" not in updates
    assert "softmin_n_traits" not in updates
    assert "softmin_auxiliary_scale" not in updates
    assert "softmin_structured_residual_scale" not in updates


def test_candidate_updates_sparse_walsh_field_replaces_pairwise_j_path() -> None:
    config = {
        "bounded_trait_space_fast_fit": {
            "bound_mode": "tanh",
            "sparse_walsh_field": {
                "enabled": True,
                "terms": 128,
                "min_order": 2,
                "max_order": 4,
                "order_weights": [0.2, 0.3, 0.5],
                "cluster_positions": [[0, 1, 2], [2, 3, 4]],
                "seed": 20260617,
                "epi_weight": 0.025,
            },
        }
    }

    updates = candidate_updates(
        base_updates={
            "n_epistatic_terms": 220,
            "epistasis_strength": 0.03,
            "empirical_pairwise_strength": 0.02,
        },
        raw={
            "sparse_walsh_epi_weight": 0.04,
            "sparse_walsh_terms": 256,
            "sparse_walsh_order2_weight": 0.1,
            "sparse_walsh_order3_weight": 0.2,
            "sparse_walsh_order4_weight": 0.7,
        },
        config=config,
    )

    assert updates["n_epistatic_terms"] == 0
    assert updates["n_higher_order_epistatic_terms"] == 0
    assert updates["empirical_pairwise_strength"] == pytest.approx(0.0)
    assert updates["generic_epistasis_target"] == "stability"
    assert updates["epistasis_strength"] == pytest.approx(0.04)
    assert updates["sparse_walsh_terms"] == 256
    assert updates["sparse_walsh_min_order"] == 2
    assert updates["sparse_walsh_max_order"] == 4
    assert updates["sparse_walsh_order_weights"] == pytest.approx([0.1, 0.2, 0.7])
    assert updates["sparse_walsh_cluster_positions"] == [[0, 1, 2], [2, 3, 4]]
    assert "sparse_walsh_epi_weight" not in updates


def test_softmin_latent_values_can_anchor_output_scale() -> None:
    graphlet = {"graphlet_id": "g0", "graphlet_dimension": 3}
    base_values = np.asarray([0.0, 1.0, 0.5, 1.5, -0.2, 0.8, 0.3, 1.3])
    config = {
        "bounded_trait_space_fast_fit": {
            "softmin_latent": {
                "enabled": True,
                "seed": 123,
                "normalize_output": True,
                "n_traits": 3,
                "auxiliary_order": 2,
                "output_center": "zero",
            }
        }
    }

    values = softmin_latent_values(
        base_values=base_values,
        graphlet=graphlet,
        raw={"softmin_output_scale": 0.27},
        config=config,
    )

    assert np.mean(values) == pytest.approx(0.0, abs=1e-12)
    assert np.std(values) == pytest.approx(0.27)


def test_stage_configs_merge_cumulatively() -> None:
    config = {
        "bounded_trait_space_fast_fit": {"random_candidates": 1},
        "cumulative_fix_stages": [
            {"label": "A", "fixes": "A", "bounded_trait_space_fast_fit": {"random_candidates": 2}},
            {"label": "A+B", "fixes": "A+B", "bounded_trait_space_fast_fit": {"local_rounds": 1}},
        ],
    }

    stages = stage_configs(config)

    assert stages[0]["bounded_trait_space_fast_fit"]["random_candidates"] == 2
    assert "local_rounds" not in stages[0]["bounded_trait_space_fast_fit"]
    assert stages[1]["bounded_trait_space_fast_fit"]["random_candidates"] == 2
    assert stages[1]["bounded_trait_space_fast_fit"]["local_rounds"] == 1
    assert stages[1]["_fix_stage"]["label"] == "A+B"


def test_monotone_piecewise_readout_is_applied() -> None:
    graphlets = [
        {
            "graphlet_id": "g0",
            "graphlet_split": "calibration",
            "empirical_values": {"latent_phi": np.asarray([0.0, 0.0, 1.0, 1.0])},
        }
    ]
    config = {
        "observation_calibration": {
            "enabled": True,
            "method": "monotone_piecewise",
            "max_knots": 4,
            "min_bin_size": 1,
        },
        "graphlet_evaluation": {"primary_score_column": "latent_phi"},
    }

    readout = fit_observation_readout(
        values_by_graphlet={"g0": np.asarray([0.0, 1.0, 2.0, 3.0])},
        graphlets=graphlets,
        config=config,
    )
    transformed = apply_observation_readout(np.asarray([0.0, 1.0, 2.0, 3.0]), readout)

    assert readout["readout_method"] == "monotone_piecewise"
    assert np.all(np.diff(transformed) >= -1e-12)
    assert transformed[0] == pytest.approx(0.0)
    assert transformed[-1] == pytest.approx(1.0)


def test_monotone_sigmoid_readout_is_restricted_two_parameter_map() -> None:
    graphlets = [
        {
            "graphlet_id": "g0",
            "graphlet_split": "calibration",
            "empirical_values": {"latent_phi": np.asarray([-1.0, -0.5, 0.5, 1.0])},
        }
    ]
    config = {
        "observation_calibration": {
            "enabled": True,
            "method": "monotone_sigmoid",
            "n_midpoint_grid": 5,
            "slope_multipliers": [0.5, 1.0, 2.0],
            "y_quantile_low": 0.0,
            "y_quantile_high": 1.0,
        },
        "graphlet_evaluation": {"primary_score_column": "latent_phi"},
    }

    readout = fit_observation_readout(
        values_by_graphlet={"g0": np.asarray([-2.0, -1.0, 1.0, 2.0])},
        graphlets=graphlets,
        config=config,
    )
    transformed = apply_observation_readout(np.asarray([-2.0, -1.0, 1.0, 2.0]), readout)

    assert readout["readout_method"] == "monotone_sigmoid"
    assert "readout_x_knots" not in readout
    assert np.all(np.diff(transformed) >= -1e-12)
    assert np.min(transformed) >= -1.0 - 1e-12
    assert np.max(transformed) <= 1.0 + 1e-12


def test_raw_latent_objective_scores_raw_pathology_terms() -> None:
    target_rows = [
        {"single_mutant_effect_sd": 2.0, "n_local_maxima": 4.0, "accessible_path_fraction": 0.25},
        {"single_mutant_effect_sd": 2.0, "n_local_maxima": 6.0, "accessible_path_fraction": 0.35},
    ]
    active_rows = [
        {"raw_single_mutant_effect_sd": 1.0, "raw_n_local_maxima": 1.0, "raw_accessible_path_fraction": 0.02},
        {"raw_single_mutant_effect_sd": 1.0, "raw_n_local_maxima": 1.0, "raw_accessible_path_fraction": 0.02},
    ]
    config = {
        "bounded_trait_space_fast_fit": {
            "raw_latent_objective": {
                "enabled": True,
                "objective_weight": 0.5,
                "terms": [
                    "single_mutant_effect_sd",
                    "n_local_maxima",
                    "accessible_path_fraction",
                ],
                "metric_weights": {"accessible_path_fraction": 4.0},
            }
        }
    }

    penalty, extras = raw_latent_objective(
        active_rows=active_rows,
        target_rows=pytest.importorskip("pandas").DataFrame(target_rows),
        config=config,
    )

    assert penalty > 0.0
    assert extras["raw_latent_objective_penalty"] == pytest.approx(penalty)
    assert extras["raw_latent_distance_component__accessible_path_fraction"] > 0.0


def test_model_free_softmin_capability_probe_emits_k_trait_rows() -> None:
    graphlets = [
        {
            "graphlet_id": "g0",
            "graphlet_split": "calibration",
            "split": "calibration",
            "graphlet_dimension": 2,
            "mutation_features": "0:A,1:B",
            "empirical_values": {"latent_phi": np.asarray([0.0, 0.7, 0.5, 1.0])},
        }
    ]
    config = {
        "graphlet_evaluation": {"primary_score_column": "latent_phi"},
        "analysis": {"random_walk": {"n_walks_per_graphlet": 2, "steps": 4, "lag_max": 2}},
        "bounded_trait_space_fast_fit": {
            "softmin_capability_probe": {
                "enabled": True,
                "k_values": [2, 3],
                "tau_values": [0.25],
                "random_candidates": 2,
                "seed": 123,
            }
        },
        "_fix_stage": {"label": "A+F+B", "fixes": "A+F+B", "index": 3},
    }

    rows = run_softmin_capability_probe(config=config, graphlets=graphlets)

    assert {row["probe_k_traits"] for row in rows} == {2, 3}
    assert all(row["fix_stage"] == "A+F+B" for row in rows)
    assert all(np.isfinite(float(row["probe_fit_rmse"])) for row in rows)
