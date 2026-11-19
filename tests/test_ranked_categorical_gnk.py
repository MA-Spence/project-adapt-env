from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from project_adapt_env.extensive_gnk import ExtensiveGNKField  # noqa: E402
from project_adapt_env.ranked_categorical_gnk import (  # noqa: E402
    RankedCategoricalGNKField,
)


def _inputs(length: int = 24) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(7800 + int(length))
    matrix = rng.normal(size=(length, 20))
    reference = rng.integers(0, 20, size=length, dtype=np.int64)
    return matrix, reference


def _ranked(
    *,
    contrast_family: str = "binary",
    representation_rank: int = 1,
    length: int = 24,
) -> RankedCategoricalGNKField:
    matrix, reference = _inputs(length)
    return RankedCategoricalGNKField(
        additive_matrix=matrix,
        reference_sequence=reference,
        seed=9127,
        orders=(2, 3, 4),
        order_weights=(0.35, 0.45, 0.20),
        target_degrees=(4.0, 5.0, 6.0),
        contrast_family=contrast_family,
        representation_rank=representation_rank,
    )


def test_binary_rank_one_exactly_replays_extensive_field() -> None:
    matrix, reference = _inputs()
    legacy = ExtensiveGNKField(
        additive_matrix=matrix,
        reference_sequence=reference,
        seed=9127,
        orders=(2, 3, 4),
        order_weights=(0.35, 0.45, 0.20),
        target_degrees=(4.0, 5.0, 6.0),
    )
    ranked = _ranked()
    sequences = np.random.default_rng(411).integers(
        0, 20, size=(96, ranked.length), dtype=np.int64
    )

    assert ranked.term_counts == legacy.term_counts
    assert ranked.fingerprint == legacy.fingerprint
    assert ranked.reference_raw_score == legacy.reference_raw_score
    np.testing.assert_array_equal(
        ranked.score_batch(sequences, reference_center=False),
        legacy.score_batch(sequences, reference_center=False),
    )
    np.testing.assert_array_equal(
        ranked.score_batch(sequences),
        legacy.score_batch(sequences),
    )
    for order in ranked.orders:
        np.testing.assert_array_equal(
            ranked.score_batch_order(sequences, order),
            legacy.score_batch_order(sequences, order),
        )
        np.testing.assert_array_equal(
            ranked.support_degrees(order), legacy.support_degrees(order)
        )
        assert ranked.support_pairs(order) == legacy.support_pairs(order)


@pytest.mark.parametrize("contrast_family", ["binary", "dense_q20"])
@pytest.mark.parametrize("representation_rank", [1, 8])
def test_ranked_field_is_deterministic_centered_and_order_resolved(
    contrast_family: str,
    representation_rank: int,
) -> None:
    first = _ranked(
        contrast_family=contrast_family,
        representation_rank=representation_rank,
    )
    second = _ranked(
        contrast_family=contrast_family,
        representation_rank=representation_rank,
    )
    sequences = np.random.default_rng(991).integers(
        0, 20, size=(64, first.length), dtype=np.int64
    )

    assert first.fingerprint == second.fingerprint
    np.testing.assert_array_equal(
        first.score_batch(sequences), second.score_batch(sequences)
    )
    assert first.score_batch(first.reference_sequence[None, :])[0] == 0.0
    np.testing.assert_array_equal(
        first.centered_score_batch(sequences), first.score_batch(sequences)
    )
    by_order = first.score_batch_by_order(sequences)
    np.testing.assert_allclose(
        sum(by_order.values()),
        first.score_batch(sequences),
        rtol=2e-13,
        atol=2e-13,
    )
    for order in first.orders:
        np.testing.assert_array_equal(
            first.score_batch_order(sequences, order), by_order[order]
        )
        np.testing.assert_allclose(
            first.score_batch_ablate_order(sequences, order),
            sum(
                values
                for active_order, values in by_order.items()
                if active_order != order
            ),
            rtol=0.0,
            atol=0.0,
        )


@pytest.mark.parametrize("contrast_family", ["binary", "dense_q20"])
def test_rank_eight_uses_exact_same_support_operator_normalization(
    contrast_family: str,
) -> None:
    field = _ranked(
        contrast_family=contrast_family,
        representation_rank=8,
    )
    total = 0.0
    variance = 0.0
    for order in field.orders:
        realized = 0.0
        for support_id in field.support_ids(order):
            coefficients = field.support_coefficients(support_id)
            operator_gram = field.support_operator_gram(support_id)
            field_gram = field.support_field_gram(support_id)
            assert operator_gram.shape == (8, 8)
            assert field_gram.shape == (8, 8)
            np.testing.assert_allclose(
                operator_gram, operator_gram.T, rtol=0.0, atol=1e-15
            )
            np.testing.assert_allclose(
                field_gram, field_gram.T, rtol=0.0, atol=1e-15
            )
            realized += float(
                coefficients @ operator_gram @ coefficients
            )
            variance += float(coefficients @ field_gram @ coefficients)
        target = (
            field.order_weights[order] * field.additive_mutation_energy
        )
        assert field.expected_order_mutation_energy[order] == pytest.approx(
            target, rel=2e-13, abs=2e-13
        )
        assert realized == pytest.approx(target, rel=2e-13, abs=2e-13)
        total += realized

    assert total == pytest.approx(
        field.additive_mutation_energy, rel=2e-13, abs=2e-13
    )
    assert field.expected_field_variance() == pytest.approx(
        variance, rel=2e-13, abs=2e-13
    )


@pytest.mark.parametrize("contrast_family", ["binary", "dense_q20"])
def test_representation_rank_does_not_change_support_topology(
    contrast_family: str,
) -> None:
    rank_one = _ranked(
        contrast_family=contrast_family, representation_rank=1
    )
    rank_eight = _ranked(
        contrast_family=contrast_family, representation_rank=8
    )

    assert rank_one.term_counts == rank_eight.term_counts
    assert rank_one.support_fingerprint == rank_eight.support_fingerprint
    for order in rank_one.orders:
        assert rank_one.support_ids(order) == rank_eight.support_ids(order)
        assert rank_one.supports(order) == rank_eight.supports(order)
        assert rank_one.support_pairs(order) == rank_eight.support_pairs(order)
        np.testing.assert_array_equal(
            rank_one.support_degrees(order),
            rank_eight.support_degrees(order),
        )
        assert rank_eight.component_counts[order] == (
            8 * rank_eight.term_counts[order]
        )
        assert int(np.sum(rank_eight.support_degrees(order))) == (
            rank_eight.term_counts[order] * order
        )


def test_binary_rank_eight_normalization_and_variance_include_cross_terms() -> None:
    field = _ranked(
        contrast_family="binary", representation_rank=8
    )
    diagonal_operator_energy = 0.0
    diagonal_field_variance = 0.0
    for support_id in field.support_ids():
        coefficients = field.support_coefficients(support_id)
        diagonal_operator_energy += float(
            np.sum(
                coefficients
                * coefficients
                * np.diag(field.support_operator_gram(support_id))
            )
        )
        diagonal_field_variance += float(
            np.sum(
                coefficients
                * coefficients
                * np.diag(field.support_field_gram(support_id))
            )
        )

    # Non-orthogonal binary components make the diagonal-only shortcut wrong.
    assert abs(
        diagonal_operator_energy - field.additive_mutation_energy
    ) > 1e-4 * field.additive_mutation_energy
    assert abs(
        diagonal_field_variance - field.expected_field_variance()
    ) > 1e-4 * field.expected_field_variance()


def test_dense_q20_components_are_centered_unit_variance_and_not_binary() -> None:
    field = _ranked(
        contrast_family="dense_q20", representation_rank=8
    )

    for term in field.terms[:24]:
        for contrast in term.contrasts:
            assert float(np.mean(contrast)) == pytest.approx(0.0, abs=1e-15)
            assert float(np.mean(contrast * contrast)) == pytest.approx(
                1.0, rel=5e-15, abs=5e-15
            )
            assert len(np.unique(contrast)) == 20
            assert float(np.max(np.abs(contrast))) < 4.0

    first_support = field.support_ids()[0]
    components = [
        term for term in field.terms if term.term_id == first_support
    ]
    assert [term.component_id for term in components] == list(range(8))
    for support_index in range(len(components[0].support)):
        contrast_matrix = np.asarray(
            [
                term.contrasts[support_index]
                for term in components
            ]
        )
        np.testing.assert_allclose(
            contrast_matrix @ contrast_matrix.T / 20.0,
            np.eye(8),
            rtol=0.0,
            atol=4e-15,
        )


def test_base_coefficient_stream_and_dense_rank_are_paired_across_arms() -> None:
    binary = _ranked(
        contrast_family="binary", representation_rank=8
    )
    dense = _ranked(
        contrast_family="dense_q20", representation_rank=8
    )
    dense_rank_one = _ranked(
        contrast_family="dense_q20", representation_rank=1
    )

    assert [
        (term.term_id, term.component_id, term.base_coefficient)
        for term in binary.terms
    ] == [
        (term.term_id, term.component_id, term.base_coefficient)
        for term in dense.terms
    ]
    for rank_one_term in dense_rank_one.terms:
        rank_eight_term = next(
            term
            for term in dense.terms
            if term.term_id == rank_one_term.term_id
            and term.component_id == 0
        )
        assert rank_one_term.base_coefficient == rank_eight_term.base_coefficient
        for rank_one_contrast, rank_eight_contrast in zip(
            rank_one_term.contrasts, rank_eight_term.contrasts
        ):
            np.testing.assert_array_equal(
                rank_one_contrast, rank_eight_contrast
            )

    for support_id in dense.support_ids()[:20]:
        eigenvalues = np.linalg.eigvalsh(
            dense.support_operator_gram(support_id)
        )
        assert int(np.sum(eigenvalues > 1e-12)) == 8


def test_ranked_field_validates_experiment_axis_and_uniform_q20_operator() -> None:
    matrix, reference = _inputs()
    with pytest.raises(ValueError, match="representation_rank"):
        RankedCategoricalGNKField(
            additive_matrix=matrix,
            reference_sequence=reference,
            seed=1,
            representation_rank=4,
        )
    with pytest.raises(ValueError, match="contrast_family"):
        RankedCategoricalGNKField(
            additive_matrix=matrix,
            reference_sequence=reference,
            seed=1,
            contrast_family="unknown",
        )
    nonuniform = np.full((matrix.shape[0], 20), 1.0 / 20.0)
    nonuniform[:, 0] += 0.01
    nonuniform[:, 1] -= 0.01
    with pytest.raises(ValueError, match="uniform q20"):
        RankedCategoricalGNKField(
            additive_matrix=matrix,
            reference_sequence=reference,
            seed=1,
            background_frequencies=nonuniform,
        )
    with pytest.raises(ValueError, match="q20"):
        RankedCategoricalGNKField(
            additive_matrix=matrix[:, :19],
            reference_sequence=np.minimum(reference, 18),
            seed=1,
        )


@pytest.mark.parametrize("contrast_family", ["binary", "dense_q20"])
def test_rank_eight_empirical_mutation_rms_matches_additive_scale(
    contrast_family: str,
) -> None:
    field = _ranked(
        contrast_family=contrast_family,
        representation_rank=8,
        length=18,
    )
    rng = np.random.default_rng(481)
    n_samples = 8_000
    source = rng.integers(
        0, 20, size=(n_samples, field.length), dtype=np.int64
    )
    sites = rng.integers(0, field.length, size=n_samples, dtype=np.int64)
    rows = np.arange(n_samples, dtype=np.int64)
    source_alleles = source[rows, sites]
    alternatives = rng.integers(0, 19, size=n_samples, dtype=np.int64)
    alternatives += alternatives >= source_alleles
    mutated = source.copy()
    mutated[rows, sites] = alternatives

    gnk_delta = field.score_batch(mutated) - field.score_batch(source)
    additive_delta = (
        field.additive_score_batch(mutated)
        - field.additive_score_batch(source)
    )
    ratio = math.sqrt(
        float(np.mean(gnk_delta * gnk_delta))
        / float(np.mean(additive_delta * additive_delta))
    )
    assert ratio == pytest.approx(1.0, rel=0.08)
