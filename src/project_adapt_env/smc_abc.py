"""Sequential Monte Carlo ABC for simulator-based calibration."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np


@dataclass(frozen=True)
class ParameterSpec:
    """One bounded calibration parameter."""

    name: str
    kind: str
    low: float
    high: float
    initial_scale_fraction: float = 0.15
    min_scale_fraction: float = 0.05

    def sample_prior(self, rng: np.random.RandomState) -> float | int:
        if self.kind == "float":
            return float(rng.uniform(self.low, self.high))
        if self.kind == "int":
            return int(rng.randint(int(self.low), int(self.high) + 1))
        raise ValueError(f"Unsupported parameter kind: {self.kind}")

    def in_support(self, value: float | int) -> bool:
        if self.kind == "int":
            return int(self.low) <= int(value) <= int(self.high)
        return self.low <= float(value) <= self.high

    def midpoint(self) -> float:
        return 0.5 * (self.low + self.high)

    def support_width(self) -> float:
        return float(self.high - self.low)


@dataclass
class SimulationResult:
    """Result from one simulator evaluation."""

    distance: float
    summary_vector: np.ndarray
    extras: dict[str, Any] = field(default_factory=dict)


@dataclass
class Particle:
    """One weighted particle."""

    parameters: dict[str, float | int]
    distance: float
    weight: float
    summary_vector: np.ndarray
    extras: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class SMCABCConfig:
    """Settings for the adaptive SMC-ABC loop."""

    n_particles: int = 48
    n_rounds: int = 4
    initial_pool_multiplier: int = 6
    proposal_pool_multiplier: int = 8
    epsilon_quantile: float = 0.6
    seed: int = 0


@dataclass
class SMCRoundSummary:
    """One adaptive population update."""

    round_index: int
    epsilon: float
    attempts: int
    best_distance: float
    median_distance: float


@dataclass
class SMCABCResult:
    """Posterior approximation from SMC-ABC."""

    particles: list[Particle]
    round_summaries: list[SMCRoundSummary]

    @property
    def best_particle(self) -> Particle:
        return min(self.particles, key=lambda particle: particle.distance)


def _reflect_into_bounds(value: float, low: float, high: float) -> float:
    if not np.isfinite(value):
        return 0.5 * (low + high)
    if low == high:
        return low
    width = high - low
    while value < low or value > high:
        if value < low:
            value = low + (low - value)
        if value > high:
            value = high - (value - high)
    if value < low:
        value = low
    if value > high:
        value = high
    return float(value if width > 0 else low)


def _parameter_value_array(particles: list[Particle], name: str) -> np.ndarray:
    return np.asarray([particle.parameters[name] for particle in particles], dtype=np.float64)


def _normalized_weights(particles: list[Particle]) -> np.ndarray:
    weights = np.asarray([particle.weight for particle in particles], dtype=np.float64)
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0:
        return np.full(len(particles), 1.0 / max(len(particles), 1))
    return weights / total


def _proposal_scales(specs: list[ParameterSpec], particles: list[Particle]) -> dict[str, float]:
    weights = _normalized_weights(particles)
    scales: dict[str, float] = {}
    for spec in specs:
        values = _parameter_value_array(particles, spec.name)
        mean = float(np.sum(weights * values))
        variance = float(np.sum(weights * (values - mean) ** 2))
        scale = math.sqrt(max(2.0 * variance, 0.0))
        minimum = max(spec.support_width() * spec.min_scale_fraction, 1e-6)
        if spec.kind == "int":
            minimum = max(minimum, 1.0)
        scales[spec.name] = float(max(scale, minimum))
    return scales


def _sample_prior(specs: list[ParameterSpec], rng: np.random.RandomState) -> dict[str, float | int]:
    return {spec.name: spec.sample_prior(rng) for spec in specs}


def _propose_parameters(
    *,
    specs: list[ParameterSpec],
    ancestor: Particle,
    scales: dict[str, float],
    rng: np.random.RandomState,
) -> dict[str, float | int]:
    proposal: dict[str, float | int] = {}
    for spec in specs:
        center = float(ancestor.parameters[spec.name])
        scale = float(scales[spec.name])
        raw = rng.normal(center, scale)
        reflected = _reflect_into_bounds(raw, spec.low, spec.high)
        if spec.kind == "int":
            proposal[spec.name] = int(np.clip(round(reflected), int(spec.low), int(spec.high)))
        else:
            proposal[spec.name] = float(reflected)
    return proposal


def _gaussian_logpdf(x: float, mean: float, scale: float) -> float:
    var = max(scale**2, 1e-12)
    return float(-0.5 * (math.log(2.0 * math.pi * var) + ((x - mean) ** 2) / var))


def _kernel_logpdf(
    specs: list[ParameterSpec],
    parameters: dict[str, float | int],
    ancestor: dict[str, float | int],
    scales: dict[str, float],
) -> float:
    return float(
        sum(
            _gaussian_logpdf(
                float(parameters[spec.name]),
                float(ancestor[spec.name]),
                float(scales[spec.name]),
            )
            for spec in specs
        )
    )


def _logsumexp(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -np.inf
    anchor = float(np.max(finite))
    return float(anchor + math.log(np.sum(np.exp(finite - anchor))))


def _mixture_logpdf(
    specs: list[ParameterSpec],
    parameters: dict[str, float | int],
    particles: list[Particle],
    scales: dict[str, float],
) -> float:
    weights = _normalized_weights(particles)
    terms = []
    for weight, particle in zip(weights, particles):
        logpdf = _kernel_logpdf(specs, parameters, particle.parameters, scales)
        terms.append(math.log(max(float(weight), 1e-300)) + logpdf)
    return _logsumexp(np.asarray(terms, dtype=np.float64))


def _prior_logpdf(specs: list[ParameterSpec], parameters: dict[str, float | int]) -> float:
    for spec in specs:
        if not spec.in_support(parameters[spec.name]):
            return -np.inf
    return 0.0


def _reweight_particles(
    *,
    specs: list[ParameterSpec],
    particles: list[Particle],
    previous_particles: list[Particle],
    scales: dict[str, float],
) -> None:
    if not previous_particles:
        uniform = 1.0 / max(len(particles), 1)
        for particle in particles:
            particle.weight = uniform
        return
    raw_weights = []
    for particle in particles:
        numerator = _prior_logpdf(specs, particle.parameters)
        denominator = _mixture_logpdf(specs, particle.parameters, previous_particles, scales)
        raw_weights.append(math.exp(numerator - denominator) if np.isfinite(denominator) else 0.0)
    total = float(np.sum(raw_weights))
    if not np.isfinite(total) or total <= 0:
        uniform = 1.0 / max(len(particles), 1)
        for particle in particles:
            particle.weight = uniform
        return
    for particle, weight in zip(particles, raw_weights):
        particle.weight = float(weight / total)


def run_smc_abc(
    *,
    specs: list[ParameterSpec],
    config: SMCABCConfig,
    simulate: Callable[[dict[str, float | int]], SimulationResult],
) -> SMCABCResult:
    """Run an adaptive SMC-ABC calibration loop."""

    if config.n_particles < 2:
        raise ValueError("SMCABCConfig.n_particles must be at least 2")
    rng = np.random.RandomState(config.seed)
    round_summaries: list[SMCRoundSummary] = []

    initial_pool_size = max(config.initial_pool_multiplier * config.n_particles, config.n_particles)
    initial_particles: list[Particle] = []
    for _ in range(initial_pool_size):
        parameters = _sample_prior(specs, rng)
        simulation = simulate(parameters)
        initial_particles.append(
            Particle(
                parameters=parameters,
                distance=float(simulation.distance),
                weight=1.0,
                summary_vector=np.asarray(simulation.summary_vector, dtype=np.float64),
                extras=dict(simulation.extras),
            )
        )
    initial_particles.sort(key=lambda particle: particle.distance)
    particles = initial_particles[: config.n_particles]
    uniform = 1.0 / float(config.n_particles)
    for particle in particles:
        particle.weight = uniform
    distances = np.asarray([particle.distance for particle in particles], dtype=np.float64)
    round_summaries.append(
        SMCRoundSummary(
            round_index=0,
            epsilon=float(np.max(distances)),
            attempts=initial_pool_size,
            best_distance=float(np.min(distances)),
            median_distance=float(np.median(distances)),
        )
    )

    for round_index in range(1, config.n_rounds):
        previous_particles = particles
        previous_weights = _normalized_weights(previous_particles)
        proposal_scales = _proposal_scales(specs, previous_particles)
        proposal_pool_size = max(config.proposal_pool_multiplier * config.n_particles, config.n_particles)
        proposal_particles: list[Particle] = []
        for _ in range(proposal_pool_size):
            ancestor_index = int(rng.choice(len(previous_particles), p=previous_weights))
            ancestor = previous_particles[ancestor_index]
            parameters = _propose_parameters(
                specs=specs,
                ancestor=ancestor,
                scales=proposal_scales,
                rng=rng,
            )
            simulation = simulate(parameters)
            proposal_particles.append(
                Particle(
                    parameters=parameters,
                    distance=float(simulation.distance),
                    weight=1.0,
                    summary_vector=np.asarray(simulation.summary_vector, dtype=np.float64),
                    extras=dict(simulation.extras),
                )
            )
        proposal_particles.sort(key=lambda particle: particle.distance)
        epsilon = float(
            np.quantile(
                np.asarray([particle.distance for particle in proposal_particles], dtype=np.float64),
                config.epsilon_quantile,
            )
        )
        accepted = [particle for particle in proposal_particles if particle.distance <= epsilon]
        if len(accepted) < config.n_particles:
            accepted = proposal_particles[: config.n_particles]
            epsilon = float(accepted[-1].distance)
        else:
            accepted = accepted[: config.n_particles]
        _reweight_particles(
            specs=specs,
            particles=accepted,
            previous_particles=previous_particles,
            scales=proposal_scales,
        )
        particles = accepted
        distances = np.asarray([particle.distance for particle in particles], dtype=np.float64)
        round_summaries.append(
            SMCRoundSummary(
                round_index=round_index,
                epsilon=epsilon,
                attempts=proposal_pool_size,
                best_distance=float(np.min(distances)),
                median_distance=float(np.median(distances)),
            )
        )

    return SMCABCResult(particles=particles, round_summaries=round_summaries)
