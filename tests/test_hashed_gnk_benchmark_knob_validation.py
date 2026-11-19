from __future__ import annotations

import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from hashed_gnk_benchmark_knob_validation import (  # noqa: E402
    build_terms,
    evaluate_parameter_set,
    stable_gaussian,
    summarize_core_guardrails,
)


def _config() -> dict:
    return {
        "generator": {
            "seed": 20260618,
            "sequence_length": 4,
            "support_modes": ["structured_block"],
            "seeds": [11],
            "epi_weight_grid": [0.0],
            "n_terms": 12,
            "order_weights": {2: 0.5, 3: 0.5},
            "core_slope": 0.8,
            "stability_offset": 1.6,
            "readout_scale": 0.5,
            "transition_half_width": 0.5,
            "low_plateau_threshold": 0.05,
            "high_plateau_threshold": 0.95,
        },
        "expectations": {
            "n_local_maxima": "increase",
            "accessible_path_fraction": "decrease",
        },
    }


def test_stable_gaussian_is_repeatable() -> None:
    first = stable_gaussian("term", 3, (1, 2), seed=17)
    second = stable_gaussian("term", 3, (1, 2), seed=17)

    assert first == pytest.approx(second)


def test_build_terms_is_deterministic_and_respects_order_weights() -> None:
    first = build_terms(
        length=6,
        seed=17,
        support_mode="structured_block",
        order_weights={2: 0.25, 3: 0.75},
        total_terms=8,
    )
    second = build_terms(
        length=6,
        seed=17,
        support_mode="structured_block",
        order_weights={2: 0.25, 3: 0.75},
        total_terms=8,
    )

    assert first == second
    assert {term.order for term in first} == {2, 3}
    assert all(len(term.support) == term.order for term in first)
    assert all(0 <= pos < 6 for term in first for pos in term.support)


def test_epi_zero_keeps_additive_latent_core_and_emits_guardrails() -> None:
    config = _config()
    parameter = {
        "phase": "WS-B",
        "parameter_set_id": "structured_block__ew0.000__seed11",
        "support_mode": "structured_block",
        "seed": 11,
        "sequence_length": 4,
        "n_terms": 12,
        "epi_weight": 0.0,
        "core_slope": 0.8,
        "stability_offset": 1.6,
        "order_weights": {2: 0.5, 3: 0.5},
    }

    metrics, supports, hamming, epistasis = evaluate_parameter_set(config, parameter)
    latent = next(row for row in metrics if row["score_layer"] == "latent_margin")
    guardrails = summarize_core_guardrails(metrics, hamming, epistasis)

    assert len(metrics) == 2
    assert supports
    assert hamming
    assert epistasis
    assert latent["walsh_high_order_fraction"] == pytest.approx(0.0)
    assert guardrails
    assert guardrails[0]["beneficial_fraction_declines_toward_peak"]
