from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (
    PROJECT_ROOT / "scripts",
    PROJECT_ROOT / "src",
    PROJECT_ROOT / "external" / "Adapt-Env",
):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import categorical_representation_rank_rescue_exp063 as runner  # noqa: E402


CONFIG_PATH = PROJECT_ROOT / (
    "experiments/2026-07-23_EXP-063_paired-categorical-code-and-"
    "representation-rank-gnk-texture-rescue/config.yaml"
)


def _flags(**overrides: bool) -> dict[str, bool]:
    values = {
        "technical_pass": True,
        "edge_pass": True,
        "actual_texture_pass": True,
        "midpoint_control_pass": True,
        "directionality_pass": True,
        "length_stability_pass": True,
        "package_pass_before_paired_superiority": True,
    }
    values.update(overrides)
    return values


def test_quick_configuration_retains_four_arms_and_two_length_extremes() -> None:
    config = runner.load_config(CONFIG_PATH, quick=True)
    runner.validate_config(config)
    references, _ = runner.load_references(PROJECT_ROOT, config)

    assert [row["id"] for row in config["candidate"]["arms"]] == list(
        runner.ARM_COMPLEXITY
    )
    assert config["intervention"]["dimensionless_amplitudes"] == [0.0, 1.0]
    assert config["acceptance"]["expected_length_bins"] == 2
    assert len(references) == 2
    assert {row["length_bin"] for row in references} == {
        "short_47_52",
        "L448",
    }


def test_decision_prefers_single_factor_rescue_and_guards_incumbent_replay() -> None:
    failed_incumbent = _flags(
        actual_texture_pass=False,
        package_pass_before_paired_superiority=False,
    )
    arm_flags = {
        "binary_r1": failed_incumbent,
        "dense_q20_r1": _flags(),
        "binary_r8": _flags(),
        "dense_q20_r8": _flags(),
    }

    assert runner.classify_decision(
        quick=False,
        pairing_pass=True,
        core_pass=True,
        arm_flags=arm_flags,
        superiority={
            "dense_q20_r1": True,
            "binary_r8": False,
            "dense_q20_r8": True,
        },
    ) == "DENSE_Q20_RANK1_SUFFICIENT"

    arm_flags["binary_r1"] = _flags()
    assert runner.classify_decision(
        quick=False,
        pairing_pass=True,
        core_pass=True,
        arm_flags=arm_flags,
        superiority={
            "dense_q20_r1": True,
            "binary_r8": True,
            "dense_q20_r8": True,
        },
    ) == "INCUMBENT_PREVALENCE_FAILURE_NOT_REPLICATED"


def test_decision_reports_the_failed_rescue_stage() -> None:
    arm_flags = {
        "binary_r1": _flags(
            actual_texture_pass=False,
            package_pass_before_paired_superiority=False,
        ),
        "dense_q20_r1": _flags(
            actual_texture_pass=False,
            directionality_pass=False,
            length_stability_pass=False,
            package_pass_before_paired_superiority=False,
        ),
        "binary_r8": _flags(
            actual_texture_pass=False,
            directionality_pass=False,
            length_stability_pass=False,
            package_pass_before_paired_superiority=False,
        ),
        "dense_q20_r8": _flags(
            actual_texture_pass=False,
            directionality_pass=False,
            length_stability_pass=False,
            package_pass_before_paired_superiority=False,
        ),
    }

    assert runner.classify_decision(
        quick=False,
        pairing_pass=True,
        core_pass=True,
        arm_flags=arm_flags,
        superiority={arm: False for arm in runner.RESCUE_ARMS},
    ) == "NO_INTERVENTION_RESCUES_TEXTURE_PREVALENCE"


def test_strict_json_replaces_nonfinite_diagnostics_with_null(
    tmp_path: Path,
) -> None:
    path = tmp_path / "strict.json"
    runner.strict_json_atomic(
        path,
        {
            "finite": np.float64(1.25),
            "nan": float("nan"),
            "positive_infinity": float("inf"),
        },
    )

    raw = path.read_text(encoding="utf-8")
    parsed = json.loads(raw)
    assert "NaN" not in raw
    assert "Infinity" not in raw
    assert parsed == {
        "finite": 1.25,
        "nan": None,
        "positive_infinity": None,
    }


def test_uniform_control_excludes_every_candidate_support_pair() -> None:
    config = runner.load_config(CONFIG_PATH, quick=True)
    rng = np.random.default_rng(63063)
    matrix = rng.normal(size=(24, 20))
    reference = rng.integers(0, 20, size=24, dtype=np.int64)
    field = runner.build_field(
        matrix=matrix,
        reference=reference,
        field_seed=6321,
        arm="binary_r1",
        config=config,
    )
    square = runner.generate_squares(
        field,
        reference,
        core_seed=6311,
        gnk_seed=6321,
        stratum="uniform_unsupported_site_pair",
        background="reference",
        count=64,
    )
    supported = {
        tuple(pair)
        for order in field.orders
        for pair in field.support_pairs(order)
    }

    assert all(tuple(pair) not in supported for pair in square["sites"])


def test_exp063_applies_registered_rank_correlation_tolerance_locally() -> None:
    config = runner.load_config(CONFIG_PATH, quick=True)
    helper = runner.helper_config(config)
    registered = config["acceptance"]["dose_response"]

    assert helper["acceptance"]["dose_response"]["spearman_min"] == (
        registered["spearman_min"] - registered["rho_comparison_atol"]
    )
