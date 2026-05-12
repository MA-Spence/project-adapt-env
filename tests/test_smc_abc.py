from __future__ import annotations

import importlib.util
import numpy as np
import pytest

import project_adapt_env.smc_abc as smc_abc_module
from project_adapt_env.adapt_env_bayes import build_parameter_specs, posterior_mean_updates
from project_adapt_env.smc_abc import (
    SMCABCBackendConfig,
    SMCABCCheckpointConfig,
    SMCABCConfig,
    SimulationResult,
    run_smc_abc,
)


def test_build_parameter_specs_round_trips_types() -> None:
    specs = build_parameter_specs(
        [
            {"name": "x", "kind": "float", "low": -2.0, "high": 2.0},
            {"name": "k", "kind": "int", "low": 1, "high": 4},
        ]
    )
    assert [spec.name for spec in specs] == ["x", "k"]
    assert specs[0].kind == "float"
    assert specs[1].kind == "int"


def test_smc_abc_recovers_simple_simulator() -> None:
    specs = build_parameter_specs(
        [
            {"name": "x", "kind": "float", "low": -4.0, "high": 4.0},
            {"name": "k", "kind": "int", "low": 0, "high": 6},
        ]
    )
    target = np.asarray([1.25, 3.0], dtype=np.float64)

    def simulate(parameters: dict[str, float | int]) -> SimulationResult:
        vector = np.asarray(
            [
                float(parameters["x"]),
                float(parameters["k"]),
            ],
            dtype=np.float64,
        )
        distance = float(np.linalg.norm(vector - target))
        return SimulationResult(distance=distance, summary_vector=vector)

    result = run_smc_abc(
        specs=specs,
        config=SMCABCConfig(
            n_particles=24,
            n_rounds=4,
            initial_pool_multiplier=5,
            proposal_pool_multiplier=6,
            epsilon_quantile=0.6,
            seed=17,
        ),
        simulate=simulate,
    )

    best = result.best_particle
    assert abs(float(best.parameters["x"]) - target[0]) < 0.5
    assert abs(int(best.parameters["k"]) - int(target[1])) <= 1

    posterior_mean = posterior_mean_updates(specs=specs, particles=result.particles)
    assert abs(float(posterior_mean["x"]) - target[0]) < 0.75
    assert abs(int(posterior_mean["k"]) - int(target[1])) <= 1


def test_smc_abc_resumes_from_checkpoint_after_interruption(tmp_path) -> None:
    specs = build_parameter_specs(
        [
            {"name": "x", "kind": "float", "low": -3.0, "high": 3.0},
            {"name": "k", "kind": "int", "low": 0, "high": 5},
        ]
    )
    target = np.asarray([0.75, 2.0], dtype=np.float64)
    checkpoint_path = tmp_path / "smc_resume.json"
    call_counter = {"count": 0}

    def flaky_simulate(parameters: dict[str, float | int]) -> SimulationResult:
        call_counter["count"] += 1
        if call_counter["count"] == 5:
            raise RuntimeError("synthetic interruption")
        vector = np.asarray(
            [float(parameters["x"]), float(parameters["k"])],
            dtype=np.float64,
        )
        return SimulationResult(
            distance=float(np.linalg.norm(vector - target)),
            summary_vector=vector,
        )

    with pytest.raises(RuntimeError, match="synthetic interruption"):
        run_smc_abc(
            specs=specs,
            config=SMCABCConfig(
                n_particles=12,
                n_rounds=3,
                initial_pool_multiplier=3,
                proposal_pool_multiplier=4,
                epsilon_quantile=0.6,
                seed=11,
            ),
            simulate=flaky_simulate,
            backend=SMCABCBackendConfig(kind="serial", max_workers=1, batch_size=2),
            checkpoint=SMCABCCheckpointConfig(path=str(checkpoint_path), resume=True),
        )

    assert checkpoint_path.is_file()

    def stable_simulate(parameters: dict[str, float | int]) -> SimulationResult:
        vector = np.asarray(
            [float(parameters["x"]), float(parameters["k"])],
            dtype=np.float64,
        )
        return SimulationResult(
            distance=float(np.linalg.norm(vector - target)),
            summary_vector=vector,
        )

    resumed = run_smc_abc(
        specs=specs,
        config=SMCABCConfig(
            n_particles=12,
            n_rounds=3,
            initial_pool_multiplier=3,
            proposal_pool_multiplier=4,
            epsilon_quantile=0.6,
            seed=11,
        ),
        simulate=stable_simulate,
        backend=SMCABCBackendConfig(kind="serial", max_workers=1, batch_size=2),
        checkpoint=SMCABCCheckpointConfig(path=str(checkpoint_path), resume=True),
    )

    assert len(resumed.round_summaries) == 3
    assert abs(float(resumed.best_particle.parameters["x"]) - target[0]) < 0.75
    assert abs(int(resumed.best_particle.parameters["k"]) - int(target[1])) <= 1


def test_smc_abc_retries_recoverable_batch_failure_without_resampling(monkeypatch) -> None:
    specs = build_parameter_specs(
        [
            {"name": "x", "kind": "float", "low": -2.0, "high": 2.0},
            {"name": "k", "kind": "int", "low": 0, "high": 4},
        ]
    )
    target = np.asarray([0.5, 2.0], dtype=np.float64)

    class RecoverableExecutor:
        def __init__(self) -> None:
            self.calls = 0
            self.recoveries = 0
            self.seen_batches: list[list[dict[str, float | int]]] = []

        def evaluate_batch(
            self,
            parameters_batch: list[dict[str, float | int]],
        ) -> list[SimulationResult]:
            self.calls += 1
            copied_batch = [dict(parameters) for parameters in parameters_batch]
            self.seen_batches.append(copied_batch)
            if self.calls == 1:
                raise RuntimeError("keepalive watchdog timeout")
            results = []
            for parameters in parameters_batch:
                vector = np.asarray(
                    [float(parameters["x"]), float(parameters["k"])],
                    dtype=np.float64,
                )
                results.append(
                    SimulationResult(
                        distance=float(np.linalg.norm(vector - target)),
                        summary_vector=vector,
                    )
                )
            return results

        def max_batch_retries(self) -> int:
            return 1

        def recover_from_batch_error(self, error: Exception) -> bool:
            if "keepalive watchdog timeout" not in str(error):
                return False
            self.recoveries += 1
            return True

        def close(self) -> None:
            return None

    executor = RecoverableExecutor()
    monkeypatch.setattr(smc_abc_module, "_build_executor", lambda **_: executor)

    result = run_smc_abc(
        specs=specs,
        config=SMCABCConfig(
            n_particles=8,
            n_rounds=2,
            initial_pool_multiplier=3,
            proposal_pool_multiplier=3,
            epsilon_quantile=0.6,
            seed=23,
        ),
        simulate=lambda parameters: SimulationResult(
            distance=float(np.linalg.norm(np.asarray([float(parameters["x"]), float(parameters["k"])]) - target)),
            summary_vector=np.asarray([float(parameters["x"]), float(parameters["k"])], dtype=np.float64),
        ),
        backend=SMCABCBackendConfig(kind="serial", max_workers=1, batch_size=2),
    )

    assert executor.recoveries == 1
    assert executor.seen_batches[0] == executor.seen_batches[1]
    assert len(result.round_summaries) == 2
    assert abs(float(result.best_particle.parameters["x"]) - target[0]) < 1.0
    assert abs(int(result.best_particle.parameters["k"]) - int(target[1])) <= 1


@pytest.mark.skipif(importlib.util.find_spec("ray") is None, reason="ray is not installed")
def test_smc_abc_ray_backend_imports_project_package_without_pythonpath(monkeypatch) -> None:
    monkeypatch.delenv("PYTHONPATH", raising=False)
    specs = build_parameter_specs(
        [
            {"name": "x", "kind": "float", "low": -2.0, "high": 2.0},
            {"name": "k", "kind": "int", "low": 0, "high": 4},
        ]
    )
    target = np.asarray([0.5, 2.0], dtype=np.float64)

    def simulate(parameters: dict[str, float | int]) -> SimulationResult:
        vector = np.asarray(
            [float(parameters["x"]), float(parameters["k"])],
            dtype=np.float64,
        )
        return SimulationResult(
            distance=float(np.linalg.norm(vector - target)),
            summary_vector=vector,
        )

    result = run_smc_abc(
        specs=specs,
        config=SMCABCConfig(
            n_particles=8,
            n_rounds=2,
            initial_pool_multiplier=3,
            proposal_pool_multiplier=3,
            epsilon_quantile=0.6,
            seed=19,
        ),
        simulate=simulate,
        backend=SMCABCBackendConfig(kind="ray", max_workers=2, batch_size=2),
        run_label="ray_import_test",
    )

    assert abs(float(result.best_particle.parameters["x"]) - target[0]) < 1.0
    assert abs(int(result.best_particle.parameters["k"]) - int(target[1])) <= 1
