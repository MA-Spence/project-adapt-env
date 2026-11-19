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

from project_adapt_env.extensive_gnk import (  # noqa: E402
    ExtensiveGNKField,
    balanced_categorical_contrast,
    uniform_alternative_substitution_energy,
)


def _field(
    length: int,
    *,
    alphabet_size: int = 20,
    data_seed: int = 1701,
    field_seed: int = 2601,
) -> ExtensiveGNKField:
    rng = np.random.default_rng(data_seed + length)
    additive_matrix = rng.normal(size=(length, alphabet_size))
    reference = rng.integers(0, alphabet_size, size=length, dtype=np.int64)
    return ExtensiveGNKField(
        additive_matrix=additive_matrix,
        reference_sequence=reference,
        seed=field_seed,
        orders=(2, 3, 4),
        order_weights=(0.35, 0.45, 0.20),
        target_degrees=(4.0, 5.0, 6.0),
    )


def _sample_single_substitutions(
    *,
    length: int,
    alphabet_size: int,
    n_samples: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    source = rng.integers(
        0,
        alphabet_size,
        size=(n_samples, length),
        dtype=np.int64,
    )
    sites = rng.integers(0, length, size=n_samples, dtype=np.int64)
    rows = np.arange(n_samples, dtype=np.int64)
    source_alleles = source[rows, sites]
    alternatives = rng.integers(
        0,
        alphabet_size - 1,
        size=n_samples,
        dtype=np.int64,
    )
    alternatives += alternatives >= source_alleles
    mutated = source.copy()
    mutated[rows, sites] = alternatives
    return source, mutated


def test_field_fingerprint_and_scores_are_deterministic() -> None:
    first = _field(24)
    second = _field(24)
    sequences = np.random.default_rng(81).integers(
        0,
        first.alphabet_size,
        size=(64, first.length),
        dtype=np.int64,
    )

    assert first.fingerprint == second.fingerprint
    assert len(first.fingerprint) == 64
    np.testing.assert_array_equal(
        first.score_batch(sequences),
        second.score_batch(sequences),
    )
    assert _field(24, field_seed=2602).fingerprint != first.fingerprint


def test_term_counts_and_order_specific_coverage_are_balanced() -> None:
    field = _field(24)

    for order in field.orders:
        expected_count = math.ceil(
            field.length * field.target_degrees[order] / order
        )
        order_terms = [term for term in field.terms if term.order == order]
        supports = [term.support for term in order_terms]
        degrees = field.support_degrees(order)

        assert field.term_counts[order] == expected_count
        assert len(order_terms) == expected_count
        assert len(set(supports)) == expected_count
        assert all(len(support) == order for support in supports)
        assert int(np.sum(degrees)) == expected_count * order
        assert int(np.min(degrees)) >= 1
        assert int(np.max(degrees) - np.min(degrees)) <= 1


def test_q20_contrast_is_bounded_balanced_and_operator_normalized() -> None:
    alphabet_size = 20
    frequencies = np.full(alphabet_size, 1.0 / alphabet_size)
    contrast = balanced_categorical_contrast(
        alphabet_size,
        seed=17,
        term_id=9,
        site=4,
        source_frequencies=frequencies,
    )

    assert float(np.mean(contrast)) == pytest.approx(0.0, abs=1e-15)
    assert float(np.mean(contrast * contrast)) == pytest.approx(1.0, abs=1e-15)
    assert float(np.max(np.abs(contrast))) <= 1.0
    assert set(np.unique(contrast)) == {-1.0, 1.0}
    assert uniform_alternative_substitution_energy(
        contrast,
        frequencies,
    ) == pytest.approx(2.0 * alphabet_size / (alphabet_size - 1), rel=1e-14)


def test_reference_centering_and_mutation_ratio_scaling_are_exact() -> None:
    field = _field(24)
    rng = np.random.default_rng(91)
    sequences = rng.integers(
        0,
        field.alphabet_size,
        size=(48, field.length),
        dtype=np.int64,
    )
    reference = field.reference_sequence[None, :]

    assert field.score_batch(reference, mutation_rms_ratio=0.73)[0] == 0.0
    assert field.score_batch(
        reference,
        mutation_rms_ratio=0.73,
        reference_center=False,
    )[0] == pytest.approx(0.73 * field.reference_raw_score)

    unit_scores = field.score_batch(sequences, mutation_rms_ratio=1.0)
    scaled_scores = field.score_batch(sequences, mutation_rms_ratio=0.37)
    np.testing.assert_allclose(scaled_scores, 0.37 * unit_scores, rtol=0.0, atol=1e-13)
    np.testing.assert_array_equal(
        field.score_batch(sequences, mutation_rms_ratio=0.0),
        np.zeros(sequences.shape[0]),
    )

    by_order = field.score_batch_by_order(sequences)
    np.testing.assert_allclose(
        sum(by_order.values()),
        unit_scores,
        rtol=1e-13,
        atol=1e-13,
    )


def test_order_energy_allocation_exactly_matches_additive_mutation_energy() -> None:
    field = _field(24)

    realized_total = 0.0
    for order in field.orders:
        explicit = sum(
            term.coefficient**2 * term.operator_energy
            for term in field.terms
            if term.order == order
        )
        target = field.order_weights[order] * field.additive_mutation_energy

        assert field.expected_order_mutation_energy[order] == pytest.approx(
            target,
            rel=1e-13,
            abs=1e-13,
        )
        assert explicit == pytest.approx(target, rel=1e-13, abs=1e-13)
        realized_total += explicit

    assert realized_total == pytest.approx(
        field.additive_mutation_energy,
        rel=1e-13,
        abs=1e-13,
    )


def test_empirical_single_mutation_ratio_is_stable_across_lengths() -> None:
    target_ratio = 0.65
    observed_ratios: list[float] = []

    for length in (24, 48):
        field = _field(length)
        source, mutated = _sample_single_substitutions(
            length=length,
            alphabet_size=field.alphabet_size,
            n_samples=8_000,
            seed=7300 + length,
        )
        gnk_effects = field.score_batch(
            mutated,
            mutation_rms_ratio=target_ratio,
        ) - field.score_batch(source, mutation_rms_ratio=target_ratio)
        additive_effects = field.additive_score_batch(
            mutated
        ) - field.additive_score_batch(source)
        observed_ratio = math.sqrt(
            float(np.mean(gnk_effects * gnk_effects))
            / float(np.mean(additive_effects * additive_effects))
        )
        observed_ratios.append(observed_ratio)

        assert observed_ratio == pytest.approx(target_ratio, rel=0.10)

    assert abs(observed_ratios[0] - observed_ratios[1]) <= 0.05
