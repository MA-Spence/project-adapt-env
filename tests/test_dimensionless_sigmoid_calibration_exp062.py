from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (
    PROJECT_ROOT / "scripts",
    PROJECT_ROOT / "src",
    PROJECT_ROOT / "external" / "Adapt-Env",
):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import dimensionless_sigmoid_calibration_exp062 as runner  # noqa: E402


def test_dimensionless_fitness_matches_raw_parameterization() -> None:
    additive = np.asarray([-3.5, -1.0, 1.25, 5.0], dtype=np.float64)
    gnk = np.asarray([0.4, -1.2, 0.0, 1.7], dtype=np.float64)
    additive_reference = 1.25
    sigma_a = 2.75
    sigma_g = 1.6
    beta = 3.25
    mu = math.log(0.95 / 0.05)
    amplitude = 0.6
    output_scale = 10.0

    dimensionless = runner.dimensionless_fitness(
        additive,
        gnk,
        additive_reference_score=additive_reference,
        sigma_a=sigma_a,
        sigma_g=sigma_g,
        beta=beta,
        mu=mu,
        amplitude=amplitude,
        output_scale=output_scale,
    )
    raw = runner.raw_parameterization(
        additive_reference_score=additive_reference,
        sigma_a=sigma_a,
        sigma_g=sigma_g,
        beta=beta,
        mu=mu,
        amplitude=amplitude,
    )

    assert raw["steepness"] == pytest.approx(beta / sigma_a)
    assert raw["threshold"] == pytest.approx(
        additive_reference - mu / raw["steepness"]
    )
    assert raw["gnk_weight"] == pytest.approx(
        amplitude / (raw["steepness"] * sigma_g)
    )
    raw_logit = raw["steepness"] * (
        additive + raw["gnk_weight"] * gnk - raw["threshold"]
    )
    raw_fitness = output_scale / (1.0 + np.exp(-raw_logit))
    np.testing.assert_allclose(dimensionless, raw_fitness, rtol=2e-15, atol=2e-15)


def test_global_scale_reports_realized_beta_and_exact_scale_recovers_target() -> None:
    beta = 4.0
    calibration_sigma_a = 2.5
    heldout_task_sigma_a = 3.0

    globally_realized = runner.realized_additive_logit_rms(
        beta=beta,
        task_sigma_a=heldout_task_sigma_a,
        scale_sigma_a=calibration_sigma_a,
    )
    exactly_realized = runner.realized_additive_logit_rms(
        beta=beta,
        task_sigma_a=heldout_task_sigma_a,
        scale_sigma_a=heldout_task_sigma_a,
    )

    assert globally_realized == pytest.approx(4.8)
    assert exactly_realized == pytest.approx(beta)
    assert globally_realized != pytest.approx(exactly_realized)


def test_core_calibration_selection_is_core_only_and_order_independent() -> None:
    mu_95 = math.log(0.95 / 0.05)
    mu_99 = math.log(0.99 / 0.01)
    rows = [
        {
            "beta": 3.0,
            "mu": mu_99,
            "core_pass": True,
            "gnk_texture_pass": True,
        },
        {
            "beta": 4.0,
            "mu": mu_95,
            "core_pass": True,
            "gnk_texture_pass": False,
        },
        {
            "beta": 4.0,
            "mu": mu_99,
            "core_pass": True,
            "gnk_texture_pass": False,
        },
        {
            "beta": 8.0,
            "mu": mu_99,
            "core_pass": False,
            "gnk_texture_pass": True,
        },
    ]
    beta_order = [0.5, 1.0, 1.5, 2.0, 2.25, 2.5, 3.0, 4.0, 8.0]
    mu_order = [mu_95, mu_99]

    selected = runner.select_core_calibration(
        rows,
        beta_order=beta_order,
        mu_order=mu_order,
    )
    reversed_selected = runner.select_core_calibration(
        list(reversed(rows)),
        beta_order=beta_order,
        mu_order=mu_order,
    )
    gnk_flipped = [
        {**row, "gnk_texture_pass": not row["gnk_texture_pass"]} for row in rows
    ]
    selected_after_gnk_flip = runner.select_core_calibration(
        gnk_flipped,
        beta_order=beta_order,
        mu_order=mu_order,
    )

    assert selected is not None
    assert selected["beta"] == 4.0
    assert selected["mu"] == mu_99
    assert (reversed_selected["beta"], reversed_selected["mu"]) == (
        selected["beta"],
        selected["mu"],
    )
    assert (
        selected_after_gnk_flip["beta"],
        selected_after_gnk_flip["mu"],
    ) == (selected["beta"], selected["mu"])

    no_qualifying_rows = [{**row, "core_pass": False} for row in rows]
    assert runner.select_core_calibration(
        no_qualifying_rows,
        beta_order=beta_order,
        mu_order=mu_order,
    ) is None


def test_zero_amplitude_paired_residual_is_exactly_zero() -> None:
    additive = np.asarray([-100.0, -2.0, 0.0, 3.0, 100.0], dtype=np.float64)
    gnk = np.asarray([10.0, -3.0, 4.0, -8.0, 2.0], dtype=np.float64)
    residual = runner.paired_residual(
        additive,
        gnk,
        additive_reference_score=0.0,
        sigma_a=2.75,
        sigma_g=1.5,
        beta=4.0,
        mu=math.log(0.99 / 0.01),
        amplitude=0.0,
        output_scale=10.0,
    )

    np.testing.assert_array_equal(residual, np.zeros_like(additive))


def test_stratified_permutation_preserves_blocks_and_replays() -> None:
    blocks = [
        (np.asarray([-2.0, -1.0, 1.0, 2.0]), np.asarray([-1.5, -0.8, 0.7, 1.8])),
        (np.asarray([-20.0, -10.0, 10.0, 20.0]), np.asarray([-9.0, -4.0, 6.0, 11.0])),
    ]
    first = runner.stratified_permutation_projection(blocks, replicates=99, seed=62062)
    replay = runner.stratified_permutation_projection(blocks, replicates=99, seed=62062)
    denominator = sum(float(np.dot(x, x)) for x, _ in blocks)
    expected = sum(float(np.dot(x, y)) for x, y in blocks) / denominator

    assert first == replay
    assert first["observed_projection"] == pytest.approx(expected)
    assert 0.0 < first["one_sided_pvalue"] <= 1.0


def test_candidate_additivized_texture_is_the_double_endpoint_secant() -> None:
    additive = [
        np.asarray([0.0, 0.2]),
        np.asarray([-0.5, 0.1]),
        np.asarray([-0.3, -0.2]),
        np.asarray([-0.9, -0.4]),
    ]
    gnk = [
        np.asarray([0.0, 0.0]),
        np.asarray([0.4, -0.2]),
        np.asarray([-0.1, 0.3]),
        np.asarray([0.8, 0.7]),
    ]
    values = runner.texture_arrays(
        additive_nodes=additive,
        gnk_nodes=gnk,
        core_mode="actual_additive_core",
        reference_score=0.0,
        scale_sigma_a=2.5,
        sigma_g=1.7,
        beta=4.0,
        mu=math.log(99.0),
        amplitude=1.0,
        output_scale=10.0,
    )
    additivized = gnk[1] + gnk[2] - gnk[0]
    expected = runner.paired_residual(
        additive[3],
        gnk[3],
        additive_reference_score=0.0,
        sigma_a=2.5,
        sigma_g=1.7,
        beta=4.0,
        mu=math.log(99.0),
        amplitude=1.0,
        output_scale=10.0,
    ) - runner.paired_residual(
        additive[3],
        additivized,
        additive_reference_score=0.0,
        sigma_a=2.5,
        sigma_g=1.7,
        beta=4.0,
        mu=math.log(99.0),
        amplitude=1.0,
        output_scale=10.0,
    )

    np.testing.assert_allclose(values["tau"], expected, rtol=2e-15, atol=2e-15)
    np.testing.assert_allclose(
        values["latent_epsilon"], gnk[3] - gnk[1] - gnk[2] + gnk[0]
    )


def test_full_reference_panel_composition_is_enforced() -> None:
    config_path = PROJECT_ROOT / (
        "experiments/2026-07-23_EXP-062_dimensionless-additive-core-sigmoid-"
        "calibration-and-held-out-gnk-dose-response-validation/config.yaml"
    )
    config = runner.load_config(config_path, quick=False)
    references, _ = runner.load_references(PROJECT_ROOT, config)
    counts = (
        pd.DataFrame(references)
        .groupby(["length_bin", "split"])
        .size()
        .unstack(fill_value=0)
    )

    assert len(counts) == 6
    assert (counts["calibration"] == 2).all()
    assert (counts["confirmation"] == 2).all()


def test_directionality_requires_every_length_bin_not_only_overall_proteins() -> None:
    doses = [0.25, 0.5, 1.0, 2.0]
    rows = []
    for bin_index in range(6):
        for protein_index in range(2):
            passing = bin_index < 5
            values = [0.1, 0.2, 0.3, 0.4] if passing else [0.4, 0.3, 0.2, 0.1]
            for dose, value in zip(doses, values):
                rows.append({
                    "scale_mode": runner.PRIMARY_GLOBAL,
                    "additive_core_mode": "actual_additive_core",
                    "field_arm": "extensive_balanced_bounded",
                    "pooled_supported_reference": True,
                    "amplitude": dose,
                    "reference_id": f"p{bin_index}_{protein_index}",
                    "length_bin": f"b{bin_index}",
                    "length": 50 + 50 * bin_index,
                    "core_seed": 1,
                    "gnk_seed": 2,
                    "direct_tau_y_rms": value,
                    "direct_tau_practical_resolution_fraction": value,
                })
    config = {
        "intervention": {"directionality_amplitudes": doses},
        "acceptance": {
            "bootstrap_replicates": 50,
            "bootstrap_seed": 62062,
            "dose_response": {
                "spearman_min": 0.8,
                "minimum_protein_pass_fraction": 0.8,
                "minimum_protein_pass_fraction_per_length_bin": 0.8,
                "minimum_cell_pass_fraction_per_length_bin": 0.8,
                "rms_growth_floor": 0.02,
                "practical_fraction_growth_floor": 0.05,
            },
        },
    }
    passed, _, summary = runner.dose_response_gate(
        pd.DataFrame(rows), mode=runner.PRIMARY_GLOBAL, config=config
    )

    assert summary["protein_pass_fraction"] > 0.8
    assert summary["minimum_protein_pass_fraction_per_length_bin"] == 0.0
    assert summary["minimum_cell_pass_fraction_per_length_bin"] == 0.0
    assert not passed


def test_transmission_requires_edge_and_texture_to_pass_in_the_same_cells() -> None:
    edge_rows = []
    texture_rows = []
    for core_seed in range(5):
        ids = {
            "reference_id": "heldout_protein",
            "length_bin": "L100_149",
            "core_seed": core_seed,
            "gnk_seed": 6221,
        }
        edge_rows.append({
            **ids,
            "scale_mode": runner.PRIMARY_GLOBAL,
            "sample_id": "reference_shell_h1_complete",
            "amplitude": 1.0,
            "direct_practical_resolution_fraction": 0.5 if core_seed < 4 else 0.0,
            "direct_correlation": 0.9 if core_seed < 4 else 0.0,
        })
        common_texture = {
            **ids,
            "scale_mode": runner.PRIMARY_GLOBAL,
            "field_arm": "extensive_balanced_bounded",
            "pooled_supported_reference": True,
            "amplitude": 1.0,
            "direct_tau_projection": 0.5,
            "direct_tau_correlation": 0.9,
            "permutation_pass": core_seed > 0,
            "direct_tau_practical_resolution_fraction": 0.5 if core_seed > 0 else 0.0,
        }
        texture_rows.extend([
            {
                **common_texture,
                "additive_core_mode": "actual_additive_core",
                "direct_tau_y_rms": 1.0,
            },
            {
                **common_texture,
                "additive_core_mode": "frozen_additive_core",
                "direct_tau_y_rms": 1.0,
            },
        ])
    config = {
        "acceptance": {
            "bootstrap_replicates": 50,
            "bootstrap_seed": 62062,
            "transmission": {
                "edge_practical_fraction_min": 0.1,
                "texture_practical_fraction_min": 0.1,
                "correlation_min": 0.5,
                "actual_to_frozen_rms_ratio_min": 0.5,
                "minimum_cell_pass_fraction_per_length_bin": 0.8,
                "projection_bootstrap_lower_min": 0.0,
            },
        },
    }

    passed, values = runner.transmission_gate(
        pd.DataFrame(edge_rows),
        pd.DataFrame(texture_rows),
        mode=runner.PRIMARY_GLOBAL,
        core_mode="actual_additive_core",
        config=config,
    )

    assert values["edge_minimum_length_bin_pass_fraction"] == pytest.approx(0.8)
    assert values["texture_minimum_length_bin_pass_fraction"] == pytest.approx(0.8)
    assert values["joint_minimum_length_bin_pass_fraction"] == pytest.approx(0.6)
    assert not passed


def test_acceptance_cleanly_stops_after_heldout_core_non_generalization() -> None:
    config_path = PROJECT_ROOT / (
        "experiments/2026-07-23_EXP-062_dimensionless-additive-core-sigmoid-"
        "calibration-and-held-out-gnk-dose-response-validation/config.yaml"
    )
    config = runner.load_config(config_path, quick=True)
    core_rows = pd.DataFrame([
        {
            "panel": "confirmation",
            "scale_mode": runner.PRIMARY_GLOBAL,
            "length_bin": "L100_149",
            "core_pass": False,
        },
        {
            "panel": "confirmation",
            "scale_mode": runner.EXACT_BETA,
            "length_bin": "L100_149",
            "core_pass": True,
        },
    ])
    checks, flags, dose_rows, fit_rows = runner.acceptance_checks(
        references=[
            {"split": "calibration", "sequence_sha256": "calibration-sequence"},
            {"split": "confirmation", "sequence_sha256": "confirmation-sequence"},
        ],
        scale_rows=[{"finite": True}],
        selection={"beta": 4.0, "mu": math.log(99.0), "reference_occupancy": 0.99},
        core_rows=core_rows,
        technical_rows=pd.DataFrame(),
        edge=pd.DataFrame(),
        texture=pd.DataFrame(),
        config=config,
    )

    assert flags["technical_pass"]
    assert not flags["confirmation_core_pass"]
    assert not flags["global_transmission_pass"]
    assert dose_rows == []
    assert all(row["stage_skipped_core_non_generalization"] for row in fit_rows)
    transmission = next(
        row for row in checks if row["check_id"] == "global_actual_core_gnk_transmission"
    )
    assert "stage_skipped_core_non_generalization" in str(transmission["value"])


@pytest.mark.parametrize(
    ("flags", "expected"),
    [
        (
            {
                "quick": True,
                "technical_pass": False,
                "selection_exists": False,
                "confirmation_core_pass": False,
                "global_transmission_pass": False,
                "global_directionality_pass": False,
                "global_length_stability_pass": False,
                "exact_full_pass": False,
                "frozen_texture_pass": False,
            },
            "QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION",
        ),
        (
            {
                "quick": False,
                "technical_pass": False,
                "selection_exists": False,
                "confirmation_core_pass": False,
                "global_transmission_pass": True,
                "global_directionality_pass": True,
                "global_length_stability_pass": True,
                "exact_full_pass": True,
                "frozen_texture_pass": True,
            },
            "EXECUTION_OR_PAIRING_INVALID",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": False,
                "confirmation_core_pass": True,
                "global_transmission_pass": True,
                "global_directionality_pass": True,
                "global_length_stability_pass": True,
                "exact_full_pass": True,
                "frozen_texture_pass": True,
            },
            "NO_CORE_OPERATING_POINT_QUALIFIES",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": False,
                "global_transmission_pass": True,
                "global_directionality_pass": True,
                "global_length_stability_pass": True,
                "exact_full_pass": True,
                "frozen_texture_pass": True,
            },
            "CORE_CALIBRATION_DOES_NOT_GENERALIZE",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": True,
                "global_transmission_pass": True,
                "global_directionality_pass": False,
                "global_length_stability_pass": True,
                "exact_full_pass": True,
                "frozen_texture_pass": True,
            },
            "CALIBRATED_GNK_TRANSMISSION_NOT_DIRECTIONAL",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": True,
                "global_transmission_pass": True,
                "global_directionality_pass": True,
                "global_length_stability_pass": False,
                "exact_full_pass": True,
                "frozen_texture_pass": True,
            },
            "CALIBRATED_GNK_TRANSMISSION_NOT_LENGTH_STABLE",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": True,
                "global_transmission_pass": True,
                "global_directionality_pass": True,
                "global_length_stability_pass": True,
                "exact_full_pass": False,
                "frozen_texture_pass": True,
            },
            "GLOBAL_OPERATOR_CALIBRATED_GNK_TRANSMISSION_VALIDATED_NOISELESS",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": True,
                "global_transmission_pass": False,
                "global_directionality_pass": False,
                "global_length_stability_pass": False,
                "exact_full_pass": True,
                "frozen_texture_pass": True,
            },
            "PER_LANDSCAPE_CORE_NORMALIZATION_REQUIRED",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": True,
                "global_transmission_pass": True,
                "global_directionality_pass": True,
                "global_length_stability_pass": True,
                "exact_full_pass": True,
                "frozen_texture_pass": False,
            },
            "CALIBRATED_SIGMOID_TRANSMISSION_FAILED",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": True,
                "global_transmission_pass": False,
                "global_directionality_pass": False,
                "global_length_stability_pass": False,
                "exact_full_pass": False,
                "frozen_texture_pass": True,
            },
            "GNK_TEXTURE_REMAINS_FROZEN_CORE_ONLY",
        ),
        (
            {
                "quick": False,
                "technical_pass": True,
                "selection_exists": True,
                "confirmation_core_pass": True,
                "global_transmission_pass": False,
                "global_directionality_pass": False,
                "global_length_stability_pass": False,
                "exact_full_pass": False,
                "frozen_texture_pass": False,
            },
            "CALIBRATED_SIGMOID_TRANSMISSION_FAILED",
        ),
    ],
)
def test_terminal_decision_precedence(flags: dict[str, bool], expected: str) -> None:
    assert runner.classify_decision(**flags) == expected
