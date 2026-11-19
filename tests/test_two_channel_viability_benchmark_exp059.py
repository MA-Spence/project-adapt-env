from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import two_channel_viability_benchmark_exp059 as runner  # noqa: E402


def minimal_config() -> dict:
    return {
        "generator": {
            "sequence_length": 48,
            "alphabet_size": 20,
            "weights": [0.0, 1.0],
            "working_weight": 1.0,
            "support_modes": ["structured", "shuffled", "random"],
            "architectures": [runner.ARCH_PRODUCT, runner.ARCH_MIXED],
            "direct_negative_control": runner.ARCH_CLOSED,
            "functional_sigma_base": 2.0,
            "gnk_terms": 32,
        },
        "calibration": {
            "seeds": [701],
            "dimension": 2,
            "alphabet_size": 20,
            "graphlets_per_seed": 1,
            "operating_points": [
                {"id": "m14_s2", "stability_margin": 14.0, "threshold_steepness": 2.0},
                {"id": "m22_s2", "stability_margin": 22.0, "threshold_steepness": 2.0},
            ],
        },
        "confirmation": {"seeds": [1701], "graphlets_per_template": 1},
        "probe": {
            "dead_sequence_samples": 16,
            "eval_chunk_size": 512,
            "site_seed_offset": 51000,
            "allele_seed_offset": 61000,
            "allele_selection": "stability_ranked_including_anchor",
        },
        "descriptors": {
            "fixed_fitness_tolerance": 1e-6,
            "single_step_effect_samples": 32,
            "pairwise_epistasis_samples": 32,
            "basin_starts": 8,
            "exact_path_starts": 4,
        },
        "gates": {
            "preflight": {
                "viable_gate_threshold": 0.90,
                "viable_graphlet_min_occupancy": 0.0,
                "dead_gate_threshold": 0.10,
                "dead_sample_min_fraction": 0.0,
                "effective_capacity_min_sd": 0.0,
                "effective_capacity_interior_low": 0.0,
                "effective_capacity_interior_high": 1.0,
                "effective_capacity_min_interior_fraction": 0.0,
                "raw_gnk_score_min_sd": 0.0,
                "product_reconstruction_atol": 1e-8,
                "closed_control_target_effective_margin": -2.0,
                "closed_control_gate_threshold": 0.10,
                "closed_control_min_low_gate_fraction": 0.90,
            },
            "g1_viability": {
                "stability_titration": {
                    "achieved_margins": [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0],
                    "transition_window": [-0.5, 0.5],
                    "saturated_window": [2.0, 4.0],
                }
            },
        },
        "_runtime": {"quick": True},
    }


def fake_design(*, passed: bool = True) -> runner.GraphletDesign:
    return runner.GraphletDesign(
        arm="calibration",
        operating_point_id="m22_s2",
        seed=701,
        dataset_id="embedded_q20_d2",
        graphlet_index=0,
        anchor_regime="viable_designed_peak",
        sampling_group_id="calibration__seed701__embedded_q20_d2__g000",
        graphlet_id="calibration__m22_s2__seed701__embedded_q20_d2__g000",
        anchor_sequence=(0,) * 48,
        sites=(0, 1),
        choices=((0, 1), (0, 1)),
        shape=(2, 2),
        target_margin=22.0,
        actual_margin=22.0,
        gate_occupancy=1.0,
        effective_capacity_sd=0.2,
        effective_capacity_interior_fraction=0.8,
        dead_gate_fraction=1.0,
        product_reconstruction_error=0.0,
        raw_gnk_score_sd=1.0,
        closed_control_low_gate_fraction=1.0,
        preflight_pass=passed,
        rejection_reason="" if passed else "test_rejection",
    )


def test_default_stage_order_and_operating_point_order_are_preserved() -> None:
    config = minimal_config()
    assert runner.configured_stage_order({"experiment": {"stage_order": ["calibration", "confirmation"]}}) == [
        "calibration",
        "confirmation",
    ]
    points = runner.operating_points(config, "calibration", PROJECT_ROOT)
    assert [point.operating_point_id for point in points] == ["m14_s2", "m22_s2"]


def test_full_dimension_is_four_but_explicit_quick_can_reduce() -> None:
    config = minimal_config()
    with pytest.raises(ValueError, match="four-site"):
        runner.validate_global_config(config, "calibration", quick=False)
    runner.validate_global_config(config, "calibration", quick=True)


def test_sampling_design_excludes_operating_point() -> None:
    config = minimal_config()
    points = runner.operating_points(config, "calibration", PROJECT_ROOT)
    first, _ = runner.prepare_designs_for_operating_point(
        config, arm="calibration", point=points[0], seeds=[701], project_root=PROJECT_ROOT
    )
    second, _ = runner.prepare_designs_for_operating_point(
        config, arm="calibration", point=points[1], seeds=[701], project_root=PROJECT_ROOT
    )
    left, right = first[0], second[0]
    assert left.sampling_group_id == right.sampling_group_id
    assert left.sites == right.sites
    assert left.choices == right.choices
    assert left.anchor_sequence == right.anchor_sequence
    assert left.operating_point_id != right.operating_point_id


def test_stability_ranked_choices_include_designed_peak_allele(monkeypatch: pytest.MonkeyPatch) -> None:
    config = minimal_config()
    monkeypatch.setattr(
        runner,
        "_template_specs",
        lambda *_args, **_kwargs: [
            {
                "dataset_id": "toy",
                "dimension": 2,
                "alphabet_cardinalities": [3, 4],
                "max_alphabet_cardinality": 4,
                "n_expected_nodes": 12,
            }
        ],
    )
    point = runner.OperatingPoint("m22_s2", 22.0, 2.0)
    designs, _ = runner.prepare_designs_for_operating_point(
        config, arm="confirmation", point=point, seeds=[1701], project_root=PROJECT_ROOT
    )
    design = designs[0]
    env = runner.build_environment(
        config,
        seed=1701,
        operating_point={"stability_margin": 0.0, "threshold_steepness": 1.0},
        routing_mode="structured",
        gnk_weight=0.0,
    )
    peak = np.asarray(env._designed_peak)
    for site, choices in zip(design.sites, design.choices):
        assert int(peak[site]) == choices[0]
        expected_alternatives = sorted(
            [allele for allele in range(20) if allele != int(peak[site])],
            key=lambda allele: -float(env.stability_matrix[site, allele]),
        )
        assert set(choices[1:]) == set(expected_alternatives[: len(choices) - 1])


def test_product_derivation_matches_direct_adaptenv_evaluation() -> None:
    config = minimal_config()
    point = {"stability_margin": 18.0, "threshold_steepness": 4.0}
    rng = np.random.default_rng(23)
    sequences = rng.integers(0, 20, size=(24, 48), dtype=np.int64)
    base_env = runner.build_environment(
        config, seed=701, operating_point=point, routing_mode="structured", gnk_weight=0.0
    )
    texture_env = runner.build_environment(
        config, seed=701, operating_point=point, routing_mode="structured", gnk_weight=1.0
    )
    direct_env = runner.build_environment(
        config, seed=701, operating_point=point, routing_mode="structured", gnk_weight=0.6
    )
    keys = ("effective_stability_margin_used", "stability_gate", "functional_capacity")
    base = runner.evaluate_components(base_env, sequences, keys=keys, chunk_size=32)
    texture = runner.evaluate_components(
        texture_env, sequences, keys=("gnk_epistatic_score",), chunk_size=32
    )["gnk_epistatic_score"]
    direct = runner.evaluate_components(
        direct_env,
        sequences,
        keys=("fitness", "effective_functional_capacity", "stability_gate"),
        chunk_size=32,
    )
    derived = runner.derive_observable(
        runner.ARCH_PRODUCT,
        weight=0.6,
        base_margin=base["effective_stability_margin_used"],
        base_gate=base["stability_gate"],
        base_function_capacity=base["functional_capacity"],
        raw_gnk_score=texture,
        threshold_steepness=4.0,
        fitness_output_scale=10.0,
    )
    np.testing.assert_allclose(derived["fitness"], direct["fitness"], rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(
        derived["effective_functional_capacity"], direct["effective_functional_capacity"], rtol=0.0, atol=1e-12
    )


def test_metric_plan_and_qary_metrics_are_paired_and_exact() -> None:
    config = minimal_config()
    first = runner.make_metric_plan(config, shape=(3, 3), sampling_group_id="shared")
    second = runner.make_metric_plan(config, shape=(3, 3), sampling_group_id="shared")
    assert first.plan_hash == second.plan_hash
    assert first.plan_hash != runner.make_metric_plan(config, shape=(3, 3), sampling_group_id="other").plan_hash

    coordinates = np.indices((3, 3)).reshape(2, -1).T
    additive = coordinates[:, 0] + 2.0 * coordinates[:, 1]
    interacting = additive + coordinates[:, 0] * coordinates[:, 1]
    assert runner.exact_interaction_fraction(additive, (3, 3)) == pytest.approx(0.0, abs=1e-14)
    assert runner.exact_interaction_fraction(interacting, (3, 3)) > 0.0


def test_fixed_tolerance_local_maxima_and_exact_shortest_paths() -> None:
    values = np.asarray([0.0, 1.0, 1.0, 2.0])
    viable = np.ones(4, dtype=bool)
    count, density = runner.masked_strict_local_maxima(values, (2, 2), viable, tolerance=1e-6)
    assert (count, density) == pytest.approx((1, 0.25))
    audit = runner.exact_shortest_path_audit(values, (2, 2), viable, [0], tolerance=1e-6)
    assert audit["accessible_shortest_path_fraction"] == pytest.approx(1.0)


def test_sign_test_counts_ties_nonpositive_and_bootstrap_uses_base_seed() -> None:
    assert runner.exact_one_sided_sign_p([1.0] * 7 + [0.0]) == pytest.approx(9.0 / 256.0)
    first = runner.fixed_seed_bootstrap_lower(
        [0.1, 0.2, 0.3], replicates=100, confidence=0.95, seed_key="x", base_seed=59059
    )
    second = runner.fixed_seed_bootstrap_lower(
        [0.1, 0.2, 0.3], replicates=100, confidence=0.95, seed_key="x", base_seed=59059
    )
    assert first == second


def test_component_preflight_records_closed_control_and_raw_gnk() -> None:
    config = minimal_config()
    point = runner.OperatingPoint("m22_s2", 22.0, 2.0)
    designs, _ = runner.prepare_designs_for_operating_point(
        config, arm="calibration", point=point, seeds=[701], project_root=PROJECT_ROOT
    )
    checked = runner.component_preflight(config, designs, point=point)[0]
    assert np.isfinite(checked.raw_gnk_score_sd)
    assert checked.closed_control_low_gate_fraction >= 0.90
    assert "closed_control" not in checked.rejection_reason


def test_titration_honours_achieved_margins() -> None:
    config = minimal_config()
    result = runner._titration_summary(
        config,
        architecture=runner.ARCH_PRODUCT,
        weight=1.0,
        point=runner.OperatingPoint("m22_s2", 22.0, 2.0),
        capacity_at_background=0.5,
        raw_gnk_at_background=0.5,
    )
    assert np.isfinite(result["titration_transition_slope"])
    assert np.isfinite(result["titration_saturated_slope"])


def test_quick_orchestration_emits_declared_outputs_on_normal_stop(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = minimal_config()
    config["experiment"] = {"stage_order": ["calibration", "confirmation"]}
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    output_dir = tmp_path / "outputs"
    design = replace(fake_design(), preflight_pass=False, rejection_reason="test_rejection")
    monkeypatch.setattr(runner, "prepare_designs_for_operating_point", lambda *_a, **_k: ([design], []))
    monkeypatch.setattr(runner, "component_preflight", lambda _c, designs, **_k: list(designs))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "runner",
            "--config",
            str(config_path),
            "--output-dir",
            str(output_dir),
            "--project-root",
            str(PROJECT_ROOT),
            "--quick",
        ],
    )
    assert runner.main() == 0
    expected = {
        "preflight.json",
        "design_catalog.csv",
        "component_metrics.csv",
        "descriptor_metrics.csv",
        "calibration_decision.csv",
        "frozen_operating_point.json",
        "confirmation_seed_metrics.csv",
        "viability_checks.csv",
        "gnk_transmission_checks.csv",
        "topography_checks.csv",
        "routing_equivalence.csv",
        "mechanism_decision.csv",
        "summary.json",
    }
    assert expected <= {path.name for path in output_dir.iterdir()}
    frozen = json.loads((output_dir / "frozen_operating_point.json").read_text(encoding="utf-8"))
    assert frozen["selected_operating_point"] is None
    summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["decision"] == "NO_COMPONENT_QUALIFIED_POINT"
