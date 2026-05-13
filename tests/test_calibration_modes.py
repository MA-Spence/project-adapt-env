from __future__ import annotations

import numpy as np
import pytest

from project_adapt_env._compat import ensure_external_paths

ensure_external_paths()

from project_adapt_env.adapt_env_bayes import (
    AdaptEnvSMCABCProblem,
    build_validation_objective_target,
)
from project_adapt_env.proteingym_panel import EmpiricalSequence
from project_adapt_env.smc_abc import ParameterSpec

from adaptenv import FitnessLandscapeEnv, LandscapeConfig
from adaptenv.calibration import (
    CalibrationOptions,
    _replace_config,
    _synthetic_readout_from_components,
)


class _EnumeratedEmpiricalLandscape:
    def __init__(self, *, name: str, wildtype: str, fitness_by_sequence: dict[str, float]) -> None:
        self.name = name
        self.wildtype = EmpiricalSequence.from_string(
            wildtype,
            name="wt",
            id="wt",
        )
        self.wildtype_sequence = self.wildtype
        self.sequences = [self.wildtype]
        self._fitness_by_sequence = dict(fitness_by_sequence)
        for sequence in sorted(fitness_by_sequence):
            if sequence == wildtype:
                continue
            self.sequences.append(EmpiricalSequence.from_string(sequence))

    def get_fitness(self, sequence: object) -> float:
        if hasattr(sequence, "to_str"):
            key = str(sequence.to_str())
        else:
            key = "".join(str(symbol) for symbol in np.asarray(sequence, dtype=object))
        return float(self._fitness_by_sequence[key])


def _pairwise_override_matrix() -> np.ndarray:
    matrix = np.zeros((20, 20), dtype=np.float64)
    matrix[1, 1] = 2.0
    return matrix


def test_empirical_pairwise_target_stability_changes_effective_stability() -> None:
    config = LandscapeConfig(
        L=4,
        seed=7,
        stability_margin=0.5,
        threshold_steepness=5.0,
        noise_amplitude=0.0,
        n_epistatic_terms=1,
        n_higher_order_epistatic_terms=0,
    )
    reference = np.asarray([0, 0, 0, 0], dtype=np.int64)
    override = {
        "positions": np.asarray([[0, 1]], dtype=np.int64),
        "matrices": _pairwise_override_matrix()[None, :, :],
    }
    mutant = np.asarray([1, 1, 0, 0], dtype=np.int64)

    score_landscape = FitnessLandscapeEnv(
        config,
        reference_sequence=reference,
        pairwise_epistasis_override=override,
        _empirical_pairwise_target="score",
    )
    stability_landscape = FitnessLandscapeEnv(
        config,
        reference_sequence=reference,
        pairwise_epistasis_override=override,
        _empirical_pairwise_target="stability",
    )

    score_components = score_landscape.evaluate(mutant, return_components=True)
    stability_components = stability_landscape.evaluate(mutant, return_components=True)

    assert stability_components["empirical_pairwise_score"] == pytest.approx(
        score_components["empirical_pairwise_score"]
    )
    assert score_components["effective_stability"] == pytest.approx(
        score_components["stability"]
    )
    assert stability_components["effective_stability"] == pytest.approx(
        stability_components["stability"]
        + stability_components["empirical_pairwise_score"]
    )
    assert stability_components["effective_stability_margin_used"] == pytest.approx(
        stability_components["stability_margin_used"]
        + stability_components["empirical_pairwise_score"]
    )


def test_empirical_pairwise_target_function_changes_effective_functional_fitness() -> None:
    config = LandscapeConfig(
        L=4,
        seed=13,
        stability_margin=0.5,
        threshold_steepness=5.0,
        noise_amplitude=0.0,
        n_epistatic_terms=1,
        n_higher_order_epistatic_terms=0,
    )
    reference = np.asarray([0, 0, 0, 0], dtype=np.int64)
    override = {
        "positions": np.asarray([[0, 1]], dtype=np.int64),
        "matrices": _pairwise_override_matrix()[None, :, :],
    }
    mutant = np.asarray([1, 1, 0, 0], dtype=np.int64)

    score_landscape = FitnessLandscapeEnv(
        config,
        reference_sequence=reference,
        pairwise_epistasis_override=override,
        _empirical_pairwise_target="score",
    )
    function_landscape = FitnessLandscapeEnv(
        config,
        reference_sequence=reference,
        pairwise_epistasis_override=override,
        _empirical_pairwise_target="function",
    )

    score_components = score_landscape.evaluate(mutant, return_components=True)
    function_components = function_landscape.evaluate(mutant, return_components=True)

    assert function_components["empirical_pairwise_score"] == pytest.approx(
        score_components["empirical_pairwise_score"]
    )
    assert score_components["effective_functional_fitness"] == pytest.approx(
        score_components["functional_fitness"]
    )
    assert score_components["effective_functional_capacity"] == pytest.approx(
        score_components["functional_capacity"]
    )
    assert function_components["effective_functional_fitness"] > function_components[
        "functional_fitness"
    ]
    assert function_components["effective_functional_capacity"] > function_components[
        "functional_capacity"
    ]


def test_synthetic_readout_mode_selects_requested_latent() -> None:
    config = LandscapeConfig(
        L=4,
        seed=11,
        stability_margin=0.5,
        threshold_steepness=5.0,
        noise_amplitude=0.0,
        n_epistatic_terms=1,
        n_higher_order_epistatic_terms=0,
    )
    reference = np.asarray([0, 0, 0, 0], dtype=np.int64)
    override = {
        "positions": np.asarray([[0, 1]], dtype=np.int64),
        "matrices": _pairwise_override_matrix()[None, :, :],
    }
    mutant = np.asarray([[1, 1, 0, 0]], dtype=np.int64)

    landscape = FitnessLandscapeEnv(
        config,
        reference_sequence=reference,
        pairwise_epistasis_override=override,
        _empirical_pairwise_target="stability",
    )
    components = landscape.evaluate_batch_components(mutant)

    fitness_values = _synthetic_readout_from_components(
        components,
        options=CalibrationOptions(synthetic_readout_mode="fitness"),
    )
    stability_values = _synthetic_readout_from_components(
        components,
        options=CalibrationOptions(synthetic_readout_mode="stability"),
    )
    margin_values = _synthetic_readout_from_components(
        components,
        options=CalibrationOptions(synthetic_readout_mode="stability_margin"),
    )
    stability_function_values = _synthetic_readout_from_components(
        components,
        options=CalibrationOptions(synthetic_readout_mode="stability_function"),
    )
    stability_binding_values = _synthetic_readout_from_components(
        components,
        options=CalibrationOptions(synthetic_readout_mode="stability_binding"),
    )

    assert fitness_values[0] == pytest.approx(float(components["fitness"][0]))
    assert stability_values[0] == pytest.approx(
        float(components["effective_stability"][0])
    )
    assert margin_values[0] == pytest.approx(
        float(components["effective_stability_margin_used"][0])
    )
    assert stability_function_values[0] == pytest.approx(
        float(components["biophysical_function_readout"][0])
    )
    assert stability_binding_values[0] == pytest.approx(
        float(components["biophysical_binding_readout"][0])
    )


def test_biophysical_binding_readout_increases_with_stability_and_function() -> None:
    config = LandscapeConfig(
        L=4,
        seed=19,
        stability_margin=0.5,
        threshold_steepness=5.0,
        noise_amplitude=0.0,
        n_epistatic_terms=0,
        n_higher_order_epistatic_terms=0,
        readout_stability_midpoint=0.0,
        readout_stability_slope=1.5,
        readout_function_midpoint=0.5,
        readout_function_slope=8.0,
    )
    reference = np.asarray([0, 0, 0, 0], dtype=np.int64)
    landscape = FitnessLandscapeEnv(config, reference_sequence=reference)

    low = landscape._biophysical_binding_readout_batch(
        effective_stability_margins=np.asarray([-2.0], dtype=np.float64),
        effective_functional_capacities=np.asarray([0.2], dtype=np.float64),
    )[0]
    high_stability = landscape._biophysical_binding_readout_batch(
        effective_stability_margins=np.asarray([2.0], dtype=np.float64),
        effective_functional_capacities=np.asarray([0.2], dtype=np.float64),
    )[0]
    high_function = landscape._biophysical_binding_readout_batch(
        effective_stability_margins=np.asarray([-2.0], dtype=np.float64),
        effective_functional_capacities=np.asarray([0.8], dtype=np.float64),
    )[0]
    high_both = landscape._biophysical_binding_readout_batch(
        effective_stability_margins=np.asarray([2.0], dtype=np.float64),
        effective_functional_capacities=np.asarray([0.8], dtype=np.float64),
    )[0]

    assert high_stability > low
    assert high_function > low
    assert high_both > high_stability
    assert high_both > high_function


def test_replace_config_maps_legacy_functional_updates_to_primary_trait_block() -> None:
    config = LandscapeConfig(
        L=4,
        latent_trait_blocks=(
            {
                "name": "readout",
                "dims": 1,
                "sigma_base": 5.0,
                "sigma_anisotropy": 0.0,
                "weight": 1.0,
            },
        ),
        observed_fitness_combine_mode="product",
        observed_fitness_terms=(
            {"source": "stability_gate"},
            {"source": "trait:readout:capacity"},
        ),
        seed=23,
    )

    updated = _replace_config(
        config,
        n_functional_dims=3,
        functional_sigma_base=9.0,
    )

    assert updated.n_functional_dims == 3
    assert updated.functional_sigma_base == pytest.approx(9.0)
    assert updated.latent_trait_blocks[0].dims == 3
    assert updated.latent_trait_blocks[0].sigma_base == pytest.approx(9.0)


def test_validation_objective_problem_prefers_matched_primary_trait_sigma() -> None:
    base_config = LandscapeConfig(
        L=4,
        use_blosum_kernel=False,
        stability_margin=1.5,
        threshold_steepness=6.0,
        latent_trait_blocks=(
            {
                "name": "readout",
                "dims": 1,
                "sigma_base": 6.0,
                "sigma_anisotropy": 0.0,
                "weight": 1.0,
            },
        ),
        observed_fitness_combine_mode="product",
        observed_fitness_terms=(
            {"source": "stability_gate"},
            {"source": "trait:readout:capacity"},
        ),
        n_epistatic_terms=1,
        epistasis_strength=0.0,
        empirical_pairwise_strength=0.025,
        n_higher_order_epistatic_terms=0,
        peak_distance_from_consensus=2,
        noise_amplitude=0.0,
        seed=41,
    )
    synthetic_empirical = FitnessLandscapeEnv(base_config, reference_sequence="AAAA")
    reference = synthetic_empirical.seq_to_str(synthetic_empirical.reference)
    variants: list[str] = [reference]
    for pos in range(synthetic_empirical.L):
        mutant = synthetic_empirical.reference.copy()
        mutant[pos] = (int(mutant[pos]) + 1) % 20
        variants.append(synthetic_empirical.seq_to_str(mutant))
    for left in range(synthetic_empirical.L):
        for right in range(left + 1, synthetic_empirical.L):
            mutant = synthetic_empirical.reference.copy()
            mutant[left] = (int(mutant[left]) + 1) % 20
            mutant[right] = (int(mutant[right]) + 1) % 20
            variants.append(synthetic_empirical.seq_to_str(mutant))
    fitness_by_sequence = {
        sequence: float(synthetic_empirical.evaluate(sequence))
        for sequence in variants
    }
    empirical = _EnumeratedEmpiricalLandscape(
        name="toy_readout",
        wildtype=reference,
        fitness_by_sequence=fitness_by_sequence,
    )
    options = CalibrationOptions(
        synthetic_seed=77,
        synthetic_readout_mode="fitness",
        empirical_pairwise_target="trait:readout",
        max_observed_sequences=256,
        max_double_mutants=256,
    )
    target = build_validation_objective_target(
        landscapes=[empirical],
        alignment_paths=[None],
        wildtypes=[reference],
        options=options,
    )
    problem = AdaptEnvSMCABCProblem(
        target=target,
        base_config=base_config,
        options=options,
        specs=[
            ParameterSpec(
                name="functional_sigma_base",
                kind="float",
                low=2.0,
                high=12.0,
            )
        ],
        replicates_per_particle=1,
        distance_mode="validation_objective",
    )

    matched = problem.simulate({"functional_sigma_base": 6.0}).distance
    mismatched = problem.simulate({"functional_sigma_base": 2.5}).distance

    assert matched < mismatched
