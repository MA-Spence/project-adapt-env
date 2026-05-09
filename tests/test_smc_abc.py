from __future__ import annotations

import numpy as np

from project_adapt_env.adapt_env_bayes import build_parameter_specs, posterior_mean_updates
from project_adapt_env.smc_abc import SMCABCConfig, SimulationResult, run_smc_abc


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
