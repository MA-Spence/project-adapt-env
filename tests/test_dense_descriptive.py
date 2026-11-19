from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from project_adapt_env.dense_descriptive import (
    CompleteGraphlet,
    component_failure_table,
    enumerate_complete_graphlets,
    exact_random_walk_metrics,
    exact_summary,
    model_distance_table,
    node_masks,
    values_for_graphlet,
)


def _score_map(features: tuple[tuple[int, str], ...], values: np.ndarray) -> dict[frozenset[tuple[int, str]], dict[str, float]]:
    graphlet = CompleteGraphlet("toy", len(features), features)
    return {
        mask: {"score": float(value)}
        for mask, value in zip(node_masks(graphlet.features), np.asarray(values, dtype=np.float64))
    }


def test_node_masks_use_binary_hypercube_order() -> None:
    features = ((2, "A"), (5, "C"), (7, "D"))

    assert node_masks(features) == [
        frozenset(),
        frozenset({(2, "A")}),
        frozenset({(5, "C")}),
        frozenset({(2, "A"), (5, "C")}),
        frozenset({(7, "D")}),
        frozenset({(2, "A"), (7, "D")}),
        frozenset({(5, "C"), (7, "D")}),
        frozenset({(2, "A"), (5, "C"), (7, "D")}),
    ]


def test_enumerate_complete_graphlets_is_deterministic_and_excludes_incomplete() -> None:
    full_features = ((0, "A"), (1, "B"), (2, "C"))
    score_map = _score_map(full_features, np.arange(8, dtype=np.float64))
    incomplete_features = ((0, "A"), (1, "B"), (3, "D"))
    incomplete = _score_map(incomplete_features, np.arange(8, dtype=np.float64))
    missing = frozenset({(0, "A"), (1, "B"), (3, "D")})
    incomplete.pop(missing)
    score_map.update(incomplete)

    first = enumerate_complete_graphlets(
        score_map,
        assay_id="toy_assay",
        score_columns=["score"],
        dimensions=[3],
        max_graphlets_per_dimension=10,
        top_feature_pool=4,
    )
    second = enumerate_complete_graphlets(
        score_map,
        assay_id="toy_assay",
        score_columns=["score"],
        dimensions=[3],
        max_graphlets_per_dimension=10,
        top_feature_pool=4,
    )

    assert first == second
    assert [graphlet.features for graphlet in first] == [full_features]


def test_values_for_graphlet_follow_binary_node_order() -> None:
    features = ((0, "A"), (1, "B"), (2, "C"))
    score_map = _score_map(features, np.arange(8, dtype=np.float64))
    graphlet = CompleteGraphlet("toy", 3, features)

    assert values_for_graphlet(score_map, graphlet, score_column="score").tolist() == list(range(8))


def test_exact_summary_of_monotone_additive_cube_has_no_high_order_or_sign_epistasis() -> None:
    values = np.asarray(
        [float(((node >> 0) & 1) + 2 * ((node >> 1) & 1) + 3 * ((node >> 2) & 1)) for node in range(8)],
        dtype=np.float64,
    )

    summary = exact_summary(values, dim=3)

    assert summary["walsh_high_order_fraction"] == pytest.approx(0.0)
    assert summary["sign_epistasis_fraction"] == pytest.approx(0.0)
    assert summary["reciprocal_sign_epistasis_fraction"] == pytest.approx(0.0)
    assert summary["n_local_maxima"] == pytest.approx(1.0)
    assert summary["n_local_minima"] == pytest.approx(1.0)
    assert summary["accessible_path_fraction"] == pytest.approx(1.0)


def test_exact_random_walk_lag1_is_deterministic() -> None:
    values = np.asarray([0.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64)

    first = exact_random_walk_metrics(values, dim=3)
    second = exact_random_walk_metrics(values, dim=3)

    assert first == second
    assert first["random_walk_step_sd"] > 0.0


def test_failure_table_flags_known_wrong_high_order_summary() -> None:
    additive = np.asarray(
        [float(((node >> 0) & 1) + ((node >> 1) & 1) + ((node >> 2) & 1)) for node in range(8)],
        dtype=np.float64,
    )
    parity = np.asarray([1.0 if node.bit_count() % 2 == 0 else -1.0 for node in range(8)], dtype=np.float64)
    rows = []
    for graphlet_id in ["g0", "g1", "g2"]:
        rows.append(
            {
                "graphlet_id": graphlet_id,
                "candidate_model": "empirical_observed__score",
                "score_column": "score",
                **exact_summary(additive, dim=3),
            }
        )
        rows.append(
            {
                "graphlet_id": graphlet_id,
                "candidate_model": "wrong_high_order_model",
                "score_column": "score",
                **exact_summary(parity, dim=3),
            }
        )
    failures = component_failure_table(
        pd.DataFrame(rows),
        target_candidate_model="empirical_observed__score",
        metric_columns=["walsh_high_order_fraction", "accessible_path_fraction", "sign_epistasis_fraction"],
    )
    distances = model_distance_table(failures)

    top_failure = failures[failures["candidate_model"] == "wrong_high_order_model"].iloc[0]
    assert top_failure["metric"] == "walsh_high_order_fraction"
    assert top_failure["abs_scaled_error"] > 1.0
    assert distances.loc[distances["candidate_model"] == "wrong_high_order_model", "n_components"].iloc[0] == 3
