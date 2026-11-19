from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from phot_nonadditive_operator_certification import pairwise_parameter_bytes, preflight_skip_reason, ratio_pass  # noqa: E402


def test_ratio_pass_applies_preregistered_recovered_to_injected_window() -> None:
    ratio, passed = ratio_pass(
        injected=2.0,
        recovered=2.6,
        ratio_min=0.7,
        ratio_max=1.4,
        near_zero_tolerance=1e-6,
    )

    assert ratio == pytest.approx(1.3)
    assert passed


def test_ratio_pass_uses_absolute_tolerance_for_near_zero_injected_metric() -> None:
    ratio, passed = ratio_pass(
        injected=1e-12,
        recovered=5e-7,
        ratio_min=0.7,
        ratio_max=1.4,
        near_zero_tolerance=1e-6,
    )

    assert math.isnan(ratio)
    assert passed


def test_ratio_pass_rejects_nonfinite_recovery() -> None:
    ratio, passed = ratio_pass(
        injected=1.0,
        recovered=math.nan,
        ratio_min=0.7,
        ratio_max=1.4,
        near_zero_tolerance=1e-6,
    )

    assert math.isnan(ratio)
    assert not passed


def test_pairwise_preflight_skips_infeasible_full_length_profile() -> None:
    cfg = {"alphabet": "protein", "gpmap_type": "pairwise"}
    reason = preflight_skip_reason(
        {
            "operator_preflight": {
                "enabled": True,
                "max_pairwise_sequence_length": 80,
                "max_pairwise_parameter_bytes": 10_000,
            }
        },
        wildtype="A" * 100,
        cfg=cfg,
    )

    assert pairwise_parameter_bytes("A" * 100, cfg) > 10_000
    assert "max_pairwise_sequence_length" in reason
