from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = PROJECT_ROOT / "scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import extensive_gnk_latent_scaling_exp060 as runner  # noqa: E402


def _valid_config() -> dict:
    acceptance = {
        "bootstrap_replicates": 20,
        "bootstrap_seed": 60060,
        "reference_center_atol": 1e-12,
        "support_max_degree_range": 1,
        "support_max_cv": 0.25,
        "support_max_duplicate_fraction": 0.01,
        "contrast_moment_atol": 1e-12,
        "local_ratio_median_band": [0.8, 1.2],
        "local_ratio_seed_band": [0.6, 1.4],
        "local_ratio_min_seed_fraction": 0.8,
        "local_slope_abs_max": 0.2,
        "local_slope_ci_abs_max": 0.3,
        "min_slope_improvement": 0.2,
        "field_variance_slope_band": [0.7, 1.3],
        "field_variance_per_length_cv_max": 0.3,
        "fractional_distance_slope_band": [0.3, 0.7],
        "distance_ratio_slope_abs_max": 0.2,
        "fixed_distance_slope_abs_max": 0.2,
        "weight_linearity_relative_error_max": 0.02,
        "weight_linearity_r2_min": 0.999,
        "order_variance_fraction_atol": 0.15,
        "curvature_slope_abs_max": 0.2,
        "q20_median_pearson_kurtosis_max": 6.0,
        "legacy_local_slope_expected_band": [-0.8, -0.2],
    }
    return {
        "experiment_id": "EXP-060",
        "experiment": {"disable_scientific_decision": False},
        "data": {},
        "references": {
            "synthetic": [
                {"id": "short", "length": 12, "seed": 101},
                {"id": "medium", "length": 16, "seed": 102},
                {"id": "long", "length": 20, "seed": 103},
            ],
            "datasets": [],
        },
        "replication": {"seeds": [2601, 2602]},
        "legacy_control": {
            "terms": 24,
            "min_order": 2,
            "max_order": 4,
            "order_weights": [0.35, 0.45, 0.20],
        },
        "candidate": {
            "orders": [2, 3, 4],
            "order_weights": [0.35, 0.45, 0.20],
            "target_degrees": [2.0, 2.0, 2.0],
        },
        "intervention": {"weights": [0.0, 1.0]},
        "probe": {
            "regimes": ["reference_shell", "uniform_background"],
            "hamming_distances": [1, 2, 4, 8, "0.1L", "0.5L"],
            "pairs_per_distance": 8,
            "eval_chunk_size": 32,
            "supported_pair_curvature_samples": 16,
        },
        "acceptance": acceptance,
        "outputs": {},
        "decision_rule": {
            "pass": "LATENT_PASS",
            "technical_fail": "TECHNICAL_FAIL",
            "local_scale_fail": "LOCAL_SCALE_FAIL",
        },
        "quick": {
            "experiment": {"disable_scientific_decision": True},
            "references": {
                "synthetic": [
                    {"id": "quick_short", "length": 12, "seed": 201},
                    {"id": "quick_long", "length": 20, "seed": 202},
                ]
            },
            "replication": {"seeds": [2601]},
            "probe": {
                "pairs_per_distance": 4,
                "supported_pair_curvature_samples": 4,
            },
        },
    }


def test_validate_config_accepts_complete_design_and_rejects_missing_requirements() -> None:
    config = _valid_config()
    runner.validate_config(config)

    missing_regime = copy.deepcopy(config)
    missing_regime["probe"]["regimes"] = ["reference_shell"]
    with pytest.raises(ValueError, match="uniform_background"):
        runner.validate_config(missing_regime)

    missing_gate = copy.deepcopy(config)
    del missing_gate["acceptance"]["curvature_slope_abs_max"]
    with pytest.raises(ValueError, match="curvature_slope_abs_max"):
        runner.validate_config(missing_gate)


def test_inline_references_are_hash_checked_without_runtime_dataset_access(
    tmp_path: Path,
) -> None:
    config = _valid_config()
    sequences = ("A" * 12, "C" * 16, "D" * 20)
    config["references"] = {
        "synthetic": [],
        "datasets": [],
        "inline": [
            {
                "id": f"inline_{len(sequence)}",
                "dataset_id": f"source_{len(sequence)}",
                "expected_length": len(sequence),
                "source": f"missing/source_{len(sequence)}.csv",
                "source_field": "target_sequence",
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            }
            for sequence in sequences
        ],
    }

    runner.validate_config(config)
    catalog = runner.load_reference_catalog(tmp_path, config)
    assert [row["length"] for row in catalog] == [12, 16, 20]
    assert [row["sequence"] for row in catalog] == list(sequences)

    config["references"]["inline"][0]["sequence_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="hash mismatch"):
        runner.load_reference_catalog(tmp_path, config)


def test_quick_main_is_classified_as_smoke_only_without_scientific_pass(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(_valid_config(), sort_keys=False), encoding="utf-8")
    output_dir = tmp_path / "outputs"
    monkeypatch.setattr(
        runner,
        "parse_args",
        lambda: SimpleNamespace(
            config=str(config_path),
            output_dir=str(output_dir),
            project_root=str(PROJECT_ROOT),
            quick=True,
        ),
    )
    monkeypatch.setenv("LABPROJ_WORKERS", "1")

    def fake_run_task(payload: dict) -> dict:
        reference = payload["reference"]
        return {
            "task": {
                "reference_id": reference["reference_id"],
                "length": reference["length"],
                "seed": payload["seed"],
                "task_signature": payload["task_signature"],
            },
            "pairing_rows": [],
            "support_rows": [],
            "contrast_rows": [],
            "metric_rows": [],
            "order_rows": [],
        }

    monkeypatch.setattr(runner, "run_task", fake_run_task)
    monkeypatch.setattr(runner, "make_scale_fits", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(
        runner,
        "acceptance_checks",
        lambda **_kwargs: [
            {
                "check_id": "smoke_execution",
                "passed": True,
                "required": True,
                "value": 1,
                "criterion": "technical smoke completed",
            }
        ],
    )

    runner.main()

    summary = json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["reference_count"] == 2
    assert summary["scientific_decision_enabled"] is False
    assert summary["scientific_pass"] is None
    assert summary["decision"] == "QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION"
    assert summary["next_stage_blocked_until_latent_pass"] is True


@pytest.mark.parametrize("regime", ["reference_shell", "uniform_background"])
def test_pair_generation_is_exact_repeatable_and_treatment_independent(regime: str) -> None:
    reference = np.arange(12, dtype=np.int64) % 20
    kwargs = {
        "reference": reference,
        "seed": 2601,
        "regime": regime,
        "distance_label": "fixed_4",
        "hamming_distance": 4,
        "n_pairs": 64,
        "alphabet_size": 20,
    }
    first_source, first_target = runner.generate_pairs(**kwargs)
    second_source, second_target = runner.generate_pairs(**kwargs)

    np.testing.assert_array_equal(first_source, second_source)
    np.testing.assert_array_equal(first_target, second_target)
    np.testing.assert_array_equal(
        np.sum(first_source != first_target, axis=1),
        np.full(64, 4),
    )
    assert np.all(first_source[first_source != first_target] != first_target[first_source != first_target])
    if regime == "reference_shell":
        np.testing.assert_array_equal(first_source, np.repeat(reference[None, :], 64, axis=0))

    first_hash = runner.json_hash(
        {"sources": runner.array_hash(first_source), "targets": runner.array_hash(first_target)}
    )
    second_hash = runner.json_hash(
        {"sources": runner.array_hash(second_source), "targets": runner.array_hash(second_target)}
    )
    assert first_hash == second_hash


def test_supported_curvature_uses_the_requested_order_for_pairs_and_scores() -> None:
    class OrderSpyField:
        length = 6
        alphabet_size = 2

        def __init__(self) -> None:
            self.support_calls: list[int] = []
            self.score_calls: list[int] = []

        def support_pairs(self, order: int) -> tuple[tuple[int, int], ...]:
            self.support_calls.append(order)
            return ((0, 1),) if order == 2 else ((2, 3),)

        def score_batch_order(self, sequences: np.ndarray, order: int) -> np.ndarray:
            self.score_calls.append(order)
            left, right, scale = (0, 1, 1.0) if order == 2 else (2, 3, 2.0)
            return scale * sequences[:, left].astype(float) * sequences[:, right].astype(float)

    order_two = OrderSpyField()
    result_two = runner.sample_supported_curvature(
        order_two,
        order=2,
        seed=81,
        n_samples=64,
        chunk_size=16,
    )
    order_three = OrderSpyField()
    result_three = runner.sample_supported_curvature(
        order_three,
        order=3,
        seed=81,
        n_samples=64,
        chunk_size=16,
    )

    assert order_two.support_calls == [2]
    assert set(order_two.score_calls) == {2}
    assert order_three.support_calls == [3]
    assert set(order_three.score_calls) == {3}
    assert result_two["n"] == result_three["n"] == 64
    assert result_two["curvature_rms"] == pytest.approx(1.0)
    assert result_three["curvature_rms"] == pytest.approx(2.0)


def test_checkpoint_is_reused_only_when_task_signature_matches(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    checkpoint = tmp_path / "checkpoint.json"
    cached = {
        "task": {"task_signature": "matching-signature"},
        "pairing_rows": [],
        "support_rows": [],
        "contrast_rows": [],
        "metric_rows": [],
        "order_rows": [],
    }
    checkpoint.write_text(json.dumps(cached), encoding="utf-8")
    base_payload = {
        "config": {"probe": {"eval_chunk_size": 8}},
        "reference": {"reference_id": "toy", "sequence": "AC"},
        "seed": 17,
        "checkpoint_path": str(checkpoint),
    }

    matching = {**base_payload, "task_signature": "matching-signature"}
    assert runner.run_task(matching) == cached

    monkeypatch.setattr(
        runner,
        "build_additive_core",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("recompute-stale")),
    )
    stale = {**base_payload, "task_signature": "new-signature"}
    with pytest.raises(RuntimeError, match="recompute-stale"):
        runner.run_task(stale)
