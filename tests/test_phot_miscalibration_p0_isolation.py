from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from phot_miscalibration_p0_isolation import (  # noqa: E402
    allocate_sparse_counts,
    bootstrap_metric_bands,
    normalize_assay_scale,
    raw_term_values_for_masks,
    support_for_term,
)


def test_support_for_term_is_deterministic_and_exact_order() -> None:
    clusters = [[0, 2, 4, 6], [1, 3, 5, 7]]

    first = support_for_term(term_id=11, order=3, seed=20260617, clusters=clusters, sequence_length=10)
    second = support_for_term(term_id=11, order=3, seed=20260617, clusters=clusters, sequence_length=10)

    assert first == second
    assert len(first) == 3
    assert len(set(first)) == 3
    assert all(0 <= pos < 10 for pos in first)


def test_allocate_sparse_counts_drops_zero_weight_orders() -> None:
    assert allocate_sparse_counts({2: 0.0, 3: 1.0}, n_terms=5) == {3: 5}


def test_bootstrap_metric_bands_reports_grouped_uncertainty() -> None:
    rows = [
        {
            "phase": "P0.1",
            "candidate_model": "target",
            "parameter_set_id": "cached",
            "score_column": "latent_phi",
            "walsh_high_order_fraction": value,
        }
        for value in [0.2, 0.4, 0.6, 0.8]
    ]

    bands = bootstrap_metric_bands(
        rows,
        metrics=["walsh_high_order_fraction"],
        iterations=50,
        seed=123,
    )

    assert len(bands) == 1
    band = bands[0]
    assert band["phase"] == "P0.1"
    assert band["candidate_model"] == "target"
    assert band["n_graphlets"] == 4
    assert band["graphlet_mean"] == pytest.approx(0.5)
    assert band["bootstrap_mean_q05"] <= band["bootstrap_mean_q50"] <= band["bootstrap_mean_q95"]


def test_synthetic_roundtrip_normalizes_wildtype_on_assay_scale() -> None:
    terms = [
        {
            "support": ((0, "A"), (1, "B")),
            "coefficient": 2.0,
        }
    ]
    masks = [
        frozenset(),
        frozenset({(0, "A")}),
        frozenset({(1, "B")}),
        frozenset({(0, "A"), (1, "B")}),
    ]

    raw = raw_term_values_for_masks(masks, terms=terms, additive_weight=0.5)
    normalized = normalize_assay_scale(raw)
    wt_raw = raw_term_values_for_masks([frozenset()], terms=terms, additive_weight=0.5)[0]
    wt_on_assay_scale = (wt_raw - float(np.mean(raw))) / float(np.std(raw))

    assert float(np.mean(normalized)) == pytest.approx(0.0)
    assert float(np.std(normalized)) == pytest.approx(1.0)
    assert wt_on_assay_scale == pytest.approx(normalized[0])
    assert wt_on_assay_scale != pytest.approx(0.0)
