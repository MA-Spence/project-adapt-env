from __future__ import annotations

import copy
import math
from pathlib import Path
import sys

import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT / "scripts", PROJECT_ROOT / "src", PROJECT_ROOT / "external" / "Adapt-Env"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import extensive_gnk_sigmoid_transmission_exp061 as runner  # noqa: E402
from project_adapt_env.extensive_gnk import ExtensiveGNKField  # noqa: E402


CONFIG_PATH = (
    PROJECT_ROOT
    / "experiments/2026-07-23_EXP-061_operator-normalized-gnk-sigmoid-transmission-paired-residual-validation/config.yaml"
)


def test_stable_residual_preserves_plateau_signal_lost_by_direct_subtraction() -> None:
    direct = runner.direct_sigmoid_residual(160.0, 0.01, output_scale=10.0, steepness=20.0)
    stable = runner.stable_sigmoid_residual(160.0, 0.01, output_scale=10.0, steepness=20.0)
    assert float(direct) == 0.0
    assert float(stable) > 0.0
    assert math.isfinite(float(stable))


def test_direct_and_stable_residual_agree_in_resolvable_dynamic_range() -> None:
    logits = np.asarray([-4.0, -1.0, 0.0, 1.0, 4.0])
    shifts = np.asarray([-0.03, 0.02, 0.01, -0.01, 0.04])
    direct = runner.direct_sigmoid_residual(logits, shifts, output_scale=10.0, steepness=20.0)
    stable = runner.stable_sigmoid_residual(logits, shifts, output_scale=10.0, steepness=20.0)
    np.testing.assert_allclose(stable, direct, rtol=1e-13, atol=1e-15)


def test_normalized_secant_gain_has_correct_zero_shift_limit() -> None:
    occupancies = np.asarray([0.1, 0.5, 0.9])
    logits = np.log(occupancies / (1.0 - occupancies))
    shifts = np.zeros(3)
    residual = runner.stable_sigmoid_residual(logits, shifts, output_scale=10.0, steepness=20.0)
    gain = runner.normalized_secant_gain(
        logits,
        shifts,
        residual,
        output_scale=10.0,
        steepness=20.0,
        floor=1e-10,
    )
    np.testing.assert_allclose(gain, 4.0 * occupancies * (1.0 - occupancies), atol=1e-15)


def test_complete_reference_single_mutant_panel_is_exact_and_deterministic() -> None:
    reference = np.arange(7, dtype=np.int64)
    first_source, first_target = runner.complete_single_mutants(reference)
    second_source, second_target = runner.complete_single_mutants(reference)
    assert first_source.shape == (19 * reference.size, reference.size)
    np.testing.assert_array_equal(first_source, second_source)
    np.testing.assert_array_equal(first_target, second_target)
    np.testing.assert_array_equal(
        np.sum(first_source != first_target, axis=1),
        np.ones(19 * reference.size, dtype=np.int64),
    )
    for site in range(reference.size):
        observed = set(first_target[19 * site : 19 * (site + 1), site].tolist())
        assert observed == set(range(20)) - {int(reference[site])}


def test_additive_surrogate_has_zero_latent_plaquette_epistasis() -> None:
    rng = np.random.default_rng(41)
    additive = rng.normal(size=(8, 20))
    reference = np.arange(8, dtype=np.int64)
    surrogate = runner.additive_surrogate_matrix(additive, seed=50041)
    base = reference.copy()
    first = base.copy(); first[2] = (first[2] + 1) % 20
    second = base.copy(); second[5] = (second[5] + 3) % 20
    both = first.copy(); both[5] = second[5]
    values = runner.score_additive_field(
        surrogate,
        reference,
        np.stack([base, first, second, both]),
    )
    epsilon = values[3] - values[1] - values[2] + values[0]
    assert abs(float(epsilon)) <= 1e-12


def test_candidate_additivized_contrast_is_exact_double_endpoint_secant() -> None:
    a = np.asarray([0.2, -0.1, 0.4, 0.3])
    g = np.asarray([0.1, -0.2, 0.5, 1.1])
    g_ab_null = g[1] + g[2] - g[0]
    theta = -0.05
    z = 20.0 * (a - theta)
    r = runner.stable_sigmoid_residual(z, 0.03 * g, output_scale=10.0, steepness=20.0)
    r_ab_null = runner.stable_sigmoid_residual(
        z[3], 0.03 * g_ab_null, output_scale=10.0, steepness=20.0,
    )
    epsilon_actual = r[3] - r[1] - r[2] + r[0]
    epsilon_null = r_ab_null - r[1] - r[2] + r[0]
    tau = epsilon_actual - epsilon_null
    np.testing.assert_allclose(tau, r[3] - r_ab_null, rtol=0.0, atol=1e-15)
    assert np.sign(float(tau)) == np.sign(float(g[3] - g_ab_null))


def test_double_endpoint_secant_retains_tiny_texture_on_extreme_plateau() -> None:
    ab_null_logit = 160.0
    weight = 0.03
    latent_epsilon = 1.0e-8
    direct_tau = runner.direct_sigmoid_residual(
        ab_null_logit,
        weight * latent_epsilon,
        output_scale=10.0,
        steepness=20.0,
    )
    stable_tau = runner.stable_sigmoid_residual(
        ab_null_logit,
        weight * latent_epsilon,
        output_scale=10.0,
        steepness=20.0,
    )
    assert float(direct_tau) == 0.0
    assert float(stable_tau) > 0.0
    assert math.isfinite(float(stable_tau))


def test_dimensionless_weight_conversion_uses_per_field_operator_rms() -> None:
    config = runner.load_config(CONFIG_PATH, quick=False)
    rows = runner.weight_grid(config, candidate_mutation_rms=2.0)
    a_half = next(row for row in rows if row["weight_label"] == "a_0.5")
    assert math.isclose(float(a_half["weight"]), 0.5 / (20.0 * 2.0))
    assert next(row for row in rows if row["weight_label"] == "fixed_0.6")["weight"] == 0.6


def test_square_generation_is_repeatable_and_treatment_independent() -> None:
    rng = np.random.default_rng(9)
    additive = rng.normal(size=(12, 20))
    reference = np.arange(12, dtype=np.int64) % 20
    field = ExtensiveGNKField(
        additive_matrix=additive,
        reference_sequence=reference,
        seed=30009,
        target_degrees=(2, 2, 2),
    )
    kwargs = dict(
        generator_seed=2601,
        stratum="candidate_supported_order_3",
        background="uniform_q20",
        count=24,
    )
    first = runner.generate_squares(field, reference, **kwargs)
    second = runner.generate_squares(field, reference, **kwargs)
    assert first["square_hash"] == second["square_hash"]
    for key in ("base", "first", "second", "both", "sites"):
        np.testing.assert_array_equal(first[key], second[key])
    assert np.all(np.sum(first["base"] != first["both"], axis=1) == 2)


def test_source_config_resolution_and_quick_filter_are_hash_stable() -> None:
    full = runner.load_config(CONFIG_PATH, quick=False)
    runner.validate_config(full)
    references, source_hash = runner.load_references(PROJECT_ROOT, full)
    assert [row["length"] for row in references] == [48, 93, 118, 220, 238, 448]
    assert len(source_hash) == 64

    quick = runner.load_config(CONFIG_PATH, quick=True)
    runner.validate_config(quick)
    quick_references, quick_source_hash = runner.load_references(PROJECT_ROOT, quick)
    assert [row["reference_id"] for row in quick_references] == ["sr43c_48", "phot_118", "gfp_238"]
    assert quick_source_hash == source_hash
    assert quick["experiment"]["disable_scientific_decision"] is True


def test_decision_branches_do_not_let_stable_diagnostic_rescue_deployment() -> None:
    config = runner.load_config(CONFIG_PATH, quick=False)
    passed, decision = runner.classify_decision(
        {
            "technical": True,
            "controlled_edges": True,
            "frozen_texture": True,
            "actual_texture": True,
            "public": False,
        },
        config,
    )
    assert passed is False
    assert decision == "TRANSITION_REGION_GNK_TRANSMISSION_VALIDATED_PUBLIC_POINT_SATURATED"

    passed, decision = runner.classify_decision(
        {
            "technical": True,
            "controlled_edges": True,
            "frozen_texture": True,
            "actual_texture": True,
            "public": True,
        },
        config,
    )
    assert passed is True
    assert decision == "CANDIDATE_GNK_MARGIN8_NOISELESS_PUBLIC_READOUT_TRANSMISSION_VALIDATED"


def test_config_rejects_nonexact_topology_gate() -> None:
    config = runner.load_config(CONFIG_PATH, quick=False)
    invalid = copy.deepcopy(config)
    invalid["acceptance"]["topology_non_tie_agreement_min"] = 0.999
    try:
        runner.validate_config(invalid)
    except ValueError as error:
        assert "topology agreement" in str(error)
    else:
        raise AssertionError("nonexact topology gate should be rejected")


def test_expected_output_contract_is_complete() -> None:
    config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    expected = set(config["outputs"]["expected"])
    assert expected == {
        "preflight.json",
        "design_catalog.csv",
        "occupancy_calibration.csv",
        "paired_residual_metrics.csv",
        "mutation_transmission_metrics.csv",
        "epistasis_transmission_metrics.csv",
        "topology_audit.csv",
        "scale_law_fits.csv",
        "acceptance_checks.csv",
        "summary.json",
    }
