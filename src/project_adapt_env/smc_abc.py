"""Sequential Monte Carlo ABC for simulator-based calibration."""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .utils import atomic_write_json, to_builtin


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


@dataclass(frozen=True)
class SMCABCBackendConfig:
    """Parallel execution settings for simulator evaluations."""

    kind: str = "serial"
    max_workers: int = 1
    batch_size: int = 1


@dataclass(frozen=True)
class SMCABCCheckpointConfig:
    """Checkpoint configuration for resumable SMC-ABC runs."""

    path: str | None = None
    resume: bool = True


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


class _SerialSimulationExecutor:
    def __init__(self, simulate: Callable[[dict[str, float | int]], SimulationResult]) -> None:
        self._simulate = simulate

    def evaluate_batch(
        self,
        parameters_batch: list[dict[str, float | int]],
    ) -> list[SimulationResult]:
        return [self._simulate(parameters) for parameters in parameters_batch]

    def max_batch_retries(self) -> int:
        return 0

    def recover_from_batch_error(self, error: Exception) -> bool:
        del error
        return False

    def close(self) -> None:
        return None


class _RaySimulationExecutor:
    def __init__(
        self,
        *,
        simulate: Callable[[dict[str, float | int]], SimulationResult],
        max_workers: int,
    ) -> None:
        try:
            import ray
        except ModuleNotFoundError as exc:  # pragma: no cover - runtime guard
            raise RuntimeError(
                "ray is required for SMCABCBackendConfig(kind='ray'). "
                "Ensure the experiment overlay includes ray."
            ) from exc

        self._ray = ray
        self._owns_runtime = not ray.is_initialized()
        self._simulate = simulate
        self._worker_count = max(int(max_workers), 1)
        src_root = Path(__file__).resolve().parents[1]
        pythonpath_entries = [str(src_root)]
        existing_pythonpath = os.environ.get("PYTHONPATH", "").strip()
        if existing_pythonpath:
            pythonpath_entries.append(existing_pythonpath)
        runtime_env = {
            "env_vars": {
                "PYTHONPATH": os.pathsep.join(pythonpath_entries),
            }
        }
        if self._owns_runtime:
            ray.init(
                num_cpus=max(int(max_workers), 1),
                include_dashboard=False,
                ignore_reinit_error=True,
                log_to_driver=False,
                runtime_env=runtime_env,
            )

        @ray.remote(num_cpus=1)
        class _SimulatorActor:
            def __init__(self, fn: Callable[[dict[str, float | int]], SimulationResult]) -> None:
                self._simulate = fn

            def simulate_one(self, parameters: dict[str, float | int]) -> SimulationResult:
                return self._simulate(parameters)

        self._actor_class = _SimulatorActor
        self._actors: list[Any] = []
        self._spawn_actors()

    def evaluate_batch(
        self,
        parameters_batch: list[dict[str, float | int]],
    ) -> list[SimulationResult]:
        if not parameters_batch:
            return []
        actor_count = min(len(self._actors), len(parameters_batch))
        refs = []
        for index, parameters in enumerate(parameters_batch):
            actor = self._actors[index % actor_count]
            refs.append(actor.simulate_one.remote(parameters))
        return list(self._ray.get(refs))

    def max_batch_retries(self) -> int:
        return 2

    def recover_from_batch_error(self, error: Exception) -> bool:
        if not self._is_recoverable_error(error):
            return False
        self._kill_actors()
        self._spawn_actors()
        return True

    def _is_recoverable_error(self, error: Exception) -> bool:
        recoverable_names = (
            "ActorUnavailableError",
            "RayActorError",
            "WorkerCrashedError",
            "ObjectLostError",
        )
        recoverable_types = tuple(
            error_type
            for name in recoverable_names
            if (error_type := getattr(self._ray.exceptions, name, None)) is not None
        )
        if recoverable_types and isinstance(error, recoverable_types):
            return True
        message = str(error)
        return "keepalive watchdog timeout" in message or "actor is temporarily unavailable" in message

    def _spawn_actors(self) -> None:
        self._actors = [
            self._actor_class.remote(self._simulate) for _ in range(self._worker_count)
        ]

    def _kill_actors(self) -> None:
        for actor in self._actors:
            try:
                self._ray.kill(actor, no_restart=True)
            except Exception:  # pragma: no cover - best-effort shutdown
                pass
        self._actors = []

    def close(self) -> None:
        self._kill_actors()
        if self._owns_runtime:
            self._ray.shutdown()


def _build_executor(
    *,
    simulate: Callable[[dict[str, float | int]], SimulationResult],
    backend: SMCABCBackendConfig,
) -> _SerialSimulationExecutor | _RaySimulationExecutor:
    kind = str(backend.kind).strip().lower()
    if kind == "serial":
        return _SerialSimulationExecutor(simulate)
    if kind == "ray":
        return _RaySimulationExecutor(simulate=simulate, max_workers=max(int(backend.max_workers), 1))
    raise ValueError(f"Unsupported SMC backend kind: {backend.kind}")


def _reflect_into_bounds(value: float, low: float, high: float) -> float:
    if not np.isfinite(value):
        return 0.5 * (low + high)
    if low == high:
        return low
    while value < low or value > high:
        if value < low:
            value = low + (low - value)
        if value > high:
            value = high - (value - high)
    return float(np.clip(value, low, high))


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


def _serialize_rng_state(state: tuple[Any, ...]) -> dict[str, Any]:
    return {
        "name": str(state[0]),
        "keys": np.asarray(state[1], dtype=np.uint32).tolist(),
        "pos": int(state[2]),
        "has_gauss": int(state[3]),
        "cached_gaussian": float(state[4]),
    }


def _deserialize_rng_state(payload: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(payload["name"]),
        np.asarray(payload["keys"], dtype=np.uint32),
        int(payload["pos"]),
        int(payload["has_gauss"]),
        float(payload["cached_gaussian"]),
    )


def _particle_to_payload(particle: Particle) -> dict[str, Any]:
    return {
        "parameters": to_builtin(particle.parameters),
        "distance": float(particle.distance),
        "weight": float(particle.weight),
        "summary_vector": np.asarray(particle.summary_vector, dtype=np.float64).tolist(),
        "extras": to_builtin(particle.extras),
    }


def _particle_from_payload(payload: dict[str, Any]) -> Particle:
    return Particle(
        parameters=dict(payload["parameters"]),
        distance=float(payload["distance"]),
        weight=float(payload["weight"]),
        summary_vector=np.asarray(payload["summary_vector"], dtype=np.float64),
        extras=dict(payload.get("extras", {})),
    )


def _round_summary_to_payload(summary: SMCRoundSummary) -> dict[str, Any]:
    return {
        "round_index": int(summary.round_index),
        "epsilon": float(summary.epsilon),
        "attempts": int(summary.attempts),
        "best_distance": float(summary.best_distance),
        "median_distance": float(summary.median_distance),
    }


def _round_summary_from_payload(payload: dict[str, Any]) -> SMCRoundSummary:
    return SMCRoundSummary(
        round_index=int(payload["round_index"]),
        epsilon=float(payload["epsilon"]),
        attempts=int(payload["attempts"]),
        best_distance=float(payload["best_distance"]),
        median_distance=float(payload["median_distance"]),
    )


def _checkpoint_payload(
    *,
    phase: str,
    current_round_index: int,
    particles: list[Particle],
    pool_particles: list[Particle],
    round_summaries: list[SMCRoundSummary],
    rng: np.random.RandomState,
) -> dict[str, Any]:
    return {
        "phase": phase,
        "current_round_index": int(current_round_index),
        "particles": [_particle_to_payload(particle) for particle in particles],
        "pool_particles": [_particle_to_payload(particle) for particle in pool_particles],
        "round_summaries": [_round_summary_to_payload(summary) for summary in round_summaries],
        "rng_state": _serialize_rng_state(rng.get_state()),
    }


def _save_checkpoint(path: Path | None, payload: dict[str, Any]) -> None:
    if path is None:
        return
    atomic_write_json(path, payload)


def _load_checkpoint(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    import json

    return json.loads(path.read_text(encoding="utf-8"))


def _emit_progress(
    progress_callback: Callable[[dict[str, Any]], None] | None,
    *,
    run_label: str,
    event: str,
    phase: str,
    round_index: int,
    completed_attempts: int,
    total_attempts: int,
    pool_particles: list[Particle],
    extra: dict[str, Any] | None = None,
) -> None:
    if progress_callback is None:
        return
    distances = np.asarray([particle.distance for particle in pool_particles], dtype=np.float64)
    payload = {
        "run_label": run_label,
        "event": event,
        "phase": phase,
        "round_index": int(round_index),
        "completed_attempts": int(completed_attempts),
        "total_attempts": int(total_attempts),
        "best_distance": float(np.min(distances)) if distances.size else None,
        "median_distance": float(np.median(distances)) if distances.size else None,
    }
    if extra:
        payload.update(to_builtin(extra))
    progress_callback(payload)


def _evaluate_pool(
    *,
    specs: list[ParameterSpec],
    total_attempts: int,
    pool_particles: list[Particle],
    executor: _SerialSimulationExecutor | _RaySimulationExecutor,
    rng: np.random.RandomState,
    backend: SMCABCBackendConfig,
    parameter_factory: Callable[[np.random.RandomState], dict[str, float | int]],
    progress_callback: Callable[[dict[str, Any]], None] | None,
    checkpoint_path: Path | None,
    checkpoint_context: dict[str, Any],
    run_label: str,
    phase: str,
    round_index: int,
    logger: Any | None = None,
) -> list[Particle]:
    del specs  # kept for interface symmetry with proposal path
    batch_size = max(int(backend.batch_size), 1)
    while len(pool_particles) < total_attempts:
        remaining = total_attempts - len(pool_particles)
        current_batch_size = min(batch_size, remaining)
        parameters_batch = [parameter_factory(rng) for _ in range(current_batch_size)]
        batch_retry_limit = max(int(executor.max_batch_retries()), 0)
        batch_retry_count = 0
        while True:
            try:
                simulations = executor.evaluate_batch(parameters_batch)
                break
            except Exception as exc:
                if batch_retry_count >= batch_retry_limit or not executor.recover_from_batch_error(exc):
                    raise
                batch_retry_count += 1
                if logger is not None:
                    logger.warning(
                        (
                            "Recovered %s batch failure at phase=%s round=%d "
                            "(completed=%d/%d, retry=%d/%d)."
                        ),
                        run_label,
                        phase,
                        round_index,
                        len(pool_particles),
                        total_attempts,
                        batch_retry_count,
                        batch_retry_limit,
                        exc_info=exc,
                    )
        for parameters, simulation in zip(parameters_batch, simulations):
            pool_particles.append(
                Particle(
                    parameters=parameters,
                    distance=float(simulation.distance),
                    weight=1.0,
                    summary_vector=np.asarray(simulation.summary_vector, dtype=np.float64),
                    extras=dict(simulation.extras),
                )
            )
        _save_checkpoint(
            checkpoint_path,
            _checkpoint_payload(
                phase=phase,
                current_round_index=round_index,
                particles=list(checkpoint_context["particles"]),
                pool_particles=pool_particles,
                round_summaries=list(checkpoint_context["round_summaries"]),
                rng=rng,
            ),
        )
        _emit_progress(
            progress_callback,
            run_label=run_label,
            event="batch_complete",
            phase=phase,
            round_index=round_index,
            completed_attempts=len(pool_particles),
            total_attempts=total_attempts,
            pool_particles=pool_particles,
        )
    return pool_particles


def run_smc_abc(
    *,
    specs: list[ParameterSpec],
    config: SMCABCConfig,
    simulate: Callable[[dict[str, float | int]], SimulationResult],
    backend: SMCABCBackendConfig | None = None,
    checkpoint: SMCABCCheckpointConfig | None = None,
    progress_callback: Callable[[dict[str, Any]], None] | None = None,
    run_label: str = "smc_abc",
    logger: Any | None = None,
) -> SMCABCResult:
    """Run an adaptive, optionally parallel and resumable SMC-ABC loop."""

    if config.n_particles < 2:
        raise ValueError("SMCABCConfig.n_particles must be at least 2")

    backend = backend or SMCABCBackendConfig()
    checkpoint_path = Path(checkpoint.path).resolve() if checkpoint and checkpoint.path else None
    checkpoint_payload = (
        _load_checkpoint(checkpoint_path) if checkpoint and checkpoint.resume else None
    )

    rng = np.random.RandomState(config.seed)
    round_summaries: list[SMCRoundSummary] = []
    particles: list[Particle] = []
    phase = "initial_pool"
    current_round_index = 0
    pool_particles: list[Particle] = []

    if checkpoint_payload is not None:
        phase = str(checkpoint_payload["phase"])
        current_round_index = int(checkpoint_payload["current_round_index"])
        particles = [_particle_from_payload(item) for item in checkpoint_payload.get("particles", [])]
        pool_particles = [
            _particle_from_payload(item) for item in checkpoint_payload.get("pool_particles", [])
        ]
        round_summaries = [
            _round_summary_from_payload(item)
            for item in checkpoint_payload.get("round_summaries", [])
        ]
        rng.set_state(_deserialize_rng_state(dict(checkpoint_payload["rng_state"])))
        if logger is not None:
            logger.info(
                "Resuming %s from checkpoint %s (phase=%s round=%d pool=%d).",
                run_label,
                checkpoint_path,
                phase,
                current_round_index,
                len(pool_particles),
            )
        if phase == "completed":
            return SMCABCResult(particles=particles, round_summaries=round_summaries)

    initial_pool_size = max(config.initial_pool_multiplier * config.n_particles, config.n_particles)
    proposal_pool_size = max(config.proposal_pool_multiplier * config.n_particles, config.n_particles)

    executor = _build_executor(simulate=simulate, backend=backend)
    try:
        if phase == "initial_pool":
            context = {"particles": [], "round_summaries": round_summaries}
            pool_particles = _evaluate_pool(
                specs=specs,
                total_attempts=initial_pool_size,
                pool_particles=pool_particles,
                executor=executor,
                rng=rng,
                backend=backend,
                parameter_factory=lambda state: _sample_prior(specs, state),
                progress_callback=progress_callback,
                checkpoint_path=checkpoint_path,
                checkpoint_context=context,
                run_label=run_label,
                phase="initial_pool",
                round_index=0,
                logger=logger,
            )
            pool_particles.sort(key=lambda particle: particle.distance)
            particles = pool_particles[: config.n_particles]
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
            pool_particles = []
            phase = "proposal_pool" if config.n_rounds > 1 else "completed"
            current_round_index = 1
            _save_checkpoint(
                checkpoint_path,
                _checkpoint_payload(
                    phase=phase,
                    current_round_index=current_round_index,
                    particles=particles,
                    pool_particles=pool_particles,
                    round_summaries=round_summaries,
                    rng=rng,
                ),
            )
            _emit_progress(
                progress_callback,
                run_label=run_label,
                event="round_complete",
                phase="initial_pool",
                round_index=0,
                completed_attempts=initial_pool_size,
                total_attempts=initial_pool_size,
                pool_particles=particles,
                extra={"epsilon": float(np.max(distances))},
            )

        while phase == "proposal_pool" and current_round_index < config.n_rounds:
            previous_particles = particles
            previous_weights = _normalized_weights(previous_particles)
            proposal_scales = _proposal_scales(specs, previous_particles)
            context = {
                "particles": previous_particles,
                "round_summaries": round_summaries,
            }
            pool_particles = _evaluate_pool(
                specs=specs,
                total_attempts=proposal_pool_size,
                pool_particles=pool_particles,
                executor=executor,
                rng=rng,
                backend=backend,
                parameter_factory=lambda state, prev=previous_particles, weights=previous_weights, scales=proposal_scales: _propose_parameters(
                    specs=specs,
                    ancestor=prev[int(state.choice(len(prev), p=weights))],
                    scales=scales,
                    rng=state,
                ),
                progress_callback=progress_callback,
                checkpoint_path=checkpoint_path,
                checkpoint_context=context,
                run_label=run_label,
                phase="proposal_pool",
                round_index=current_round_index,
                logger=logger,
            )
            pool_particles.sort(key=lambda particle: particle.distance)
            epsilon = float(
                np.quantile(
                    np.asarray([particle.distance for particle in pool_particles], dtype=np.float64),
                    config.epsilon_quantile,
                )
            )
            accepted = [particle for particle in pool_particles if particle.distance <= epsilon]
            if len(accepted) < config.n_particles:
                accepted = pool_particles[: config.n_particles]
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
                    round_index=current_round_index,
                    epsilon=epsilon,
                    attempts=proposal_pool_size,
                    best_distance=float(np.min(distances)),
                    median_distance=float(np.median(distances)),
                )
            )
            pool_particles = []
            current_round_index += 1
            phase = "proposal_pool" if current_round_index < config.n_rounds else "completed"
            _save_checkpoint(
                checkpoint_path,
                _checkpoint_payload(
                    phase=phase,
                    current_round_index=current_round_index,
                    particles=particles,
                    pool_particles=pool_particles,
                    round_summaries=round_summaries,
                    rng=rng,
                ),
            )
            _emit_progress(
                progress_callback,
                run_label=run_label,
                event="round_complete",
                phase="proposal_pool",
                round_index=current_round_index - 1,
                completed_attempts=proposal_pool_size,
                total_attempts=proposal_pool_size,
                pool_particles=particles,
                extra={"epsilon": epsilon},
            )

        _save_checkpoint(
            checkpoint_path,
            _checkpoint_payload(
                phase="completed",
                current_round_index=current_round_index,
                particles=particles,
                pool_particles=[],
                round_summaries=round_summaries,
                rng=rng,
            ),
        )
        return SMCABCResult(particles=particles, round_summaries=round_summaries)
    finally:
        executor.close()
