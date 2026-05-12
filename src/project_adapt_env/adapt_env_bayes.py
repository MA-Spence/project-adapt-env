"""Project-local Bayesian calibration for Adapt-Env summary matching."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np

from ._compat import ensure_external_paths
from .smc_abc import (
    ParameterSpec,
    SMCABCBackendConfig,
    SMCABCCheckpointConfig,
    SMCABCConfig,
    SMCABCResult,
    SimulationResult,
    run_smc_abc,
)
from .utils import atomic_write_dataframe_csv, atomic_write_json, to_builtin, weighted_quantile

ensure_external_paths()

from adaptenv import LandscapeConfig  # noqa: E402
from adaptenv.calibration import (  # noqa: E402
    CalibrationOptions,
    LandscapeCalibrator,
    _aggregate_synthetic_stats,
    _build_synthetic_landscape_for_summary,
    _summarize_effect_distribution,
    _synthetic_summary_for_landscape,
    summarize_empirical_landscapes,
)


DEFAULT_FEATURES = (
    "mean",
    "variance",
    "skewness",
    "fraction_neutral",
    "fraction_beneficial",
    "fraction_deleterious",
    "fraction_lethal",
    "conservation_correlation",
    "epistasis_variance",
)


@dataclass
class SummaryTarget:
    """Observed summary vector plus bootstrap covariance for ABC distance."""

    feature_names: list[str]
    observed_vector: np.ndarray
    covariance: np.ndarray
    inverse_covariance: np.ndarray
    bootstrap_vectors: np.ndarray
    empirical_collection: Any


def build_parameter_specs(config: list[dict[str, Any]]) -> list[ParameterSpec]:
    return [
        ParameterSpec(
            name=str(item["name"]),
            kind=str(item["kind"]),
            low=float(item["low"]),
            high=float(item["high"]),
            initial_scale_fraction=float(item.get("initial_scale_fraction", 0.15)),
            min_scale_fraction=float(item.get("min_scale_fraction", 0.05)),
        )
        for item in config
    ]


def _coerce_float(value: float | None) -> float:
    if value is None or not np.isfinite(value):
        return 0.0
    return float(value)


def _feature_row_from_effect_summary(
    *,
    assay_name: str,
    effect_summary: dict[str, Any],
    conservation_correlation: float | None,
    epistasis_variance: float | None,
    selected_features: tuple[str, ...],
) -> dict[str, float]:
    row: dict[str, float] = {}
    for feature in selected_features:
        if feature == "conservation_correlation":
            value = conservation_correlation
        elif feature == "epistasis_variance":
            value = epistasis_variance
        else:
            value = effect_summary.get(feature)
        row[f"{assay_name}__{feature}"] = _coerce_float(value)
    return row


def _flatten_feature_rows(
    rows: list[dict[str, float]],
    *,
    feature_names: list[str] | None = None,
) -> tuple[list[str], np.ndarray]:
    merged: dict[str, float] = {}
    for row in rows:
        merged.update(row)
    ordered_names = list(feature_names) if feature_names is not None else sorted(merged)
    vector = np.asarray([merged[name] for name in ordered_names], dtype=np.float64)
    return ordered_names, vector


def _bootstrap_empirical_vector(
    summaries: list[Any],
    *,
    options: CalibrationOptions,
    selected_features: tuple[str, ...],
    rng: np.random.RandomState,
    feature_names: list[str],
) -> np.ndarray:
    rows: list[dict[str, float]] = []
    for summary in summaries:
        observed = np.asarray(summary.observed_mask, dtype=bool)
        effects = np.asarray(summary.single_mutant_effects[observed], dtype=np.float64)
        fitnesses = np.asarray(summary.single_mutant_fitnesses[observed], dtype=np.float64)
        if effects.size:
            idx = rng.choice(effects.size, size=effects.size, replace=True)
            effect_summary = _summarize_effect_distribution(
                effects[idx],
                mutant_fitnesses=fitnesses[idx],
                wildtype_fitness=summary.wildtype_fitness,
                options=options,
            )
        else:
            effect_summary = summary.effect_summary

        conservation_correlation = summary.conservation_sensitivity_correlation
        if summary.conservation is not None and summary.position_sensitivity is not None:
            x = np.asarray(summary.conservation, dtype=np.float64)
            y = np.asarray(summary.position_sensitivity, dtype=np.float64)
            valid = np.isfinite(x) & np.isfinite(y)
            if np.sum(valid) >= 2:
                valid_idx = np.flatnonzero(valid)
                sample_idx = rng.choice(valid_idx, size=valid_idx.size, replace=True)
                x_sample = x[sample_idx]
                y_sample = y[sample_idx]
                if np.std(x_sample) > 0 and np.std(y_sample) > 0:
                    conservation_correlation = float(np.corrcoef(x_sample, y_sample)[0, 1])

        epistasis_variance = summary.epistasis_variance
        pairwise = summary.pairwise_epistasis
        if pairwise is not None:
            pairwise = np.asarray(pairwise, dtype=np.float64)
            pairwise = pairwise[np.isfinite(pairwise)]
            if pairwise.size:
                pair_idx = rng.choice(pairwise.size, size=pairwise.size, replace=True)
                epistasis_variance = float(np.var(pairwise[pair_idx]))

        rows.append(
            _feature_row_from_effect_summary(
                assay_name=str(summary.name),
                effect_summary=effect_summary,
                conservation_correlation=conservation_correlation,
                epistasis_variance=epistasis_variance,
                selected_features=selected_features,
            )
        )
    _, vector = _flatten_feature_rows(rows, feature_names=feature_names)
    return vector


def build_empirical_target(
    *,
    landscapes: list[Any],
    alignment_paths: list[Any],
    wildtypes: list[Any],
    options: CalibrationOptions,
    selected_features: tuple[str, ...] = DEFAULT_FEATURES,
    bootstrap_replicates: int = 256,
    covariance_shrinkage: float = 0.25,
    covariance_ridge: float = 1e-6,
    progress_callback: Any | None = None,
) -> SummaryTarget:
    collection = summarize_empirical_landscapes(
        landscapes,
        kind="functional",
        alignment_profiles=alignment_paths,
        wildtypes=wildtypes,
        options=options,
    )
    observed_rows = [
        _feature_row_from_effect_summary(
            assay_name=str(summary.name),
            effect_summary=summary.effect_summary,
            conservation_correlation=summary.conservation_sensitivity_correlation,
            epistasis_variance=summary.epistasis_variance,
            selected_features=selected_features,
        )
        for summary in collection.per_landscape
    ]
    feature_names, observed_vector = _flatten_feature_rows(observed_rows)

    rng = np.random.RandomState(options.synthetic_seed + 404)
    bootstrap_vectors_list = []
    total_replicates = int(bootstrap_replicates)
    for replicate_index in range(total_replicates):
        bootstrap_vectors_list.append(
            _bootstrap_empirical_vector(
                list(collection.per_landscape),
                options=options,
                selected_features=selected_features,
                rng=rng,
                feature_names=feature_names,
            )
        )
        if progress_callback is not None and (
            (replicate_index + 1) == total_replicates or (replicate_index + 1) % 16 == 0
        ):
            progress_callback(
                {
                    "event": "bootstrap_progress",
                    "completed": replicate_index + 1,
                    "total": total_replicates,
                }
            )
    bootstrap_vectors = np.stack(bootstrap_vectors_list, axis=0)
    covariance = np.cov(bootstrap_vectors, rowvar=False)
    covariance = np.asarray(covariance, dtype=np.float64)
    if covariance.ndim == 0:
        covariance = covariance.reshape(1, 1)
    diagonal = np.diag(np.diag(covariance))
    shrinkage = float(np.clip(covariance_shrinkage, 0.0, 1.0))
    covariance = (1.0 - shrinkage) * covariance + shrinkage * diagonal
    ridge_scale = max(float(np.nanmedian(np.diag(diagonal))), 1.0)
    covariance = covariance + np.eye(covariance.shape[0], dtype=np.float64) * float(
        covariance_ridge * ridge_scale
    )
    inverse_covariance = np.linalg.pinv(covariance)

    return SummaryTarget(
        feature_names=feature_names,
        observed_vector=observed_vector,
        covariance=covariance,
        inverse_covariance=inverse_covariance,
        bootstrap_vectors=bootstrap_vectors,
        empirical_collection=collection,
    )


def serialize_summary_target(target: SummaryTarget) -> dict[str, Any]:
    return {
        "feature_names": list(target.feature_names),
        "observed_vector": np.asarray(target.observed_vector, dtype=np.float64).tolist(),
        "covariance": np.asarray(target.covariance, dtype=np.float64).tolist(),
        "inverse_covariance": np.asarray(target.inverse_covariance, dtype=np.float64).tolist(),
        "bootstrap_vectors": np.asarray(target.bootstrap_vectors, dtype=np.float64).tolist(),
    }


def restore_summary_target(
    *,
    payload: dict[str, Any],
    empirical_collection: Any,
) -> SummaryTarget:
    return SummaryTarget(
        feature_names=[str(item) for item in payload["feature_names"]],
        observed_vector=np.asarray(payload["observed_vector"], dtype=np.float64),
        covariance=np.asarray(payload["covariance"], dtype=np.float64),
        inverse_covariance=np.asarray(payload["inverse_covariance"], dtype=np.float64),
        bootstrap_vectors=np.asarray(payload["bootstrap_vectors"], dtype=np.float64),
        empirical_collection=empirical_collection,
    )


def _simulate_vector_once(
    *,
    summaries: list[Any],
    base_config: LandscapeConfig,
    options: CalibrationOptions,
    updates: dict[str, Any],
    selected_features: tuple[str, ...],
) -> np.ndarray:
    rows: list[dict[str, float]] = []
    for idx, summary in enumerate(summaries):
        landscape = _build_synthetic_landscape_for_summary(
            summary,
            base_config=base_config,
            options=options,
            landscape_index=idx,
            updates=updates,
        )
        synthetic = _synthetic_summary_for_landscape(
            landscape,
            summary,
            options=options,
        )
        rows.append(
            _feature_row_from_effect_summary(
                assay_name=str(summary.name),
                effect_summary=synthetic["effect_summary"],
                conservation_correlation=synthetic["conservation_correlation"],
                epistasis_variance=synthetic["epistasis_variance"],
                selected_features=selected_features,
            )
        )
    _, vector = _flatten_feature_rows(rows)
    return vector


def posterior_mean_updates(
    *,
    specs: list[ParameterSpec],
    particles: list[Any],
) -> dict[str, float | int]:
    weights = np.asarray([particle.weight for particle in particles], dtype=np.float64)
    total = float(np.sum(weights))
    if total <= 0 or not np.isfinite(total):
        weights = np.full(len(particles), 1.0 / max(len(particles), 1))
    else:
        weights = weights / total
    updates: dict[str, float | int] = {}
    for spec in specs:
        values = np.asarray([particle.parameters[spec.name] for particle in particles], dtype=np.float64)
        mean = float(np.sum(weights * values))
        if spec.kind == "int":
            updates[spec.name] = int(np.clip(round(mean), int(spec.low), int(spec.high)))
        else:
            updates[spec.name] = float(np.clip(mean, spec.low, spec.high))
    return updates


class AdaptEnvSMCABCProblem:
    """Bridge Adapt-Env summaries into the generic SMC-ABC loop."""

    def __init__(
        self,
        *,
        target: SummaryTarget,
        base_config: LandscapeConfig,
        options: CalibrationOptions,
        specs: list[ParameterSpec],
        selected_features: tuple[str, ...] = DEFAULT_FEATURES,
        replicates_per_particle: int = 1,
    ) -> None:
        self.target = target
        self.base_config = base_config
        self.options = options
        self.specs = specs
        self.selected_features = selected_features
        self.replicates_per_particle = max(int(replicates_per_particle), 1)
        self._cache: dict[tuple[tuple[str, float | int], ...], SimulationResult] = {}

    def _cache_key(self, parameters: dict[str, float | int]) -> tuple[tuple[str, float | int], ...]:
        normalized = []
        for spec in self.specs:
            value = parameters[spec.name]
            if spec.kind == "float":
                normalized.append((spec.name, round(float(value), 8)))
            else:
                normalized.append((spec.name, int(value)))
        return tuple(normalized)

    def simulate(self, parameters: dict[str, float | int]) -> SimulationResult:
        cache_key = self._cache_key(parameters)
        if cache_key in self._cache:
            return self._cache[cache_key]

        vectors = []
        for replicate_index in range(self.replicates_per_particle):
            replicate_options = replace(
                self.options,
                synthetic_seed=int(self.options.synthetic_seed) + 1009 * replicate_index,
            )
            vector = _simulate_vector_once(
                summaries=list(self.target.empirical_collection.per_landscape),
                base_config=self.base_config,
                options=replicate_options,
                updates=parameters,
                selected_features=self.selected_features,
            )
            vectors.append(vector)
        summary_vector = np.mean(np.stack(vectors, axis=0), axis=0)
        delta = summary_vector - self.target.observed_vector
        distance = float(
            math.sqrt(
                max(
                    float(delta.T @ self.target.inverse_covariance @ delta) / max(delta.size, 1),
                    0.0,
                )
            )
        )
        result = SimulationResult(
            distance=distance,
            summary_vector=summary_vector,
            extras={
                "delta": delta,
                "replicate_count": self.replicates_per_particle,
            },
        )
        self._cache[cache_key] = result
        return result


def validate_updates(
    *,
    empirical_collection: Any,
    base_config: LandscapeConfig,
    options: CalibrationOptions,
    updates: dict[str, Any],
) -> tuple[LandscapeConfig, dict[str, float | None]]:
    config = replace(base_config, **updates)
    calibrator = LandscapeCalibrator(base_config=base_config, options=options)
    validation = calibrator._validate_fit(
        config,
        unfolding_summary=None,
        functional_summary=empirical_collection,
    )
    return config, validation


def run_synthetic_truth_recovery(
    *,
    target: SummaryTarget,
    base_config: LandscapeConfig,
    options: CalibrationOptions,
    specs: list[ParameterSpec],
    smc_config: SMCABCConfig,
    backend: SMCABCBackendConfig | None,
    truths: list[dict[str, Any]],
    replicates_per_particle: int,
    selected_features: tuple[str, ...] = DEFAULT_FEATURES,
    checkpoint_dir: Path | None = None,
    partial_csv_path: Path | None = None,
    progress_callback: Any | None = None,
    logger: Any | None = None,
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    completed_by_label: dict[str, dict[str, Any]] = {}
    if checkpoint_dir is not None:
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
    if partial_csv_path is not None and partial_csv_path.is_file():
        try:
            import pandas as pd

            existing_frame = pd.read_csv(partial_csv_path)
            for row in existing_frame.to_dict(orient="records"):
                completed_by_label[str(row["label"])] = row
        except Exception:
            completed_by_label = {}
    for truth_index, truth in enumerate(truths):
        truth_updates = dict(truth["updates"])
        truth_label = str(truth["label"])
        if truth_label in completed_by_label:
            results.append(dict(completed_by_label[truth_label]))
            continue
        truth_options = replace(options, synthetic_seed=int(options.synthetic_seed) + 5000 + truth_index * 97)
        truth_vector = _simulate_vector_once(
            summaries=list(target.empirical_collection.per_landscape),
            base_config=base_config,
            options=truth_options,
            updates=truth_updates,
            selected_features=selected_features,
        )
        truth_target = SummaryTarget(
            feature_names=list(target.feature_names),
            observed_vector=truth_vector,
            covariance=target.covariance,
            inverse_covariance=target.inverse_covariance,
            bootstrap_vectors=target.bootstrap_vectors,
            empirical_collection=target.empirical_collection,
        )
        problem = AdaptEnvSMCABCProblem(
            target=truth_target,
            base_config=base_config,
            options=options,
            specs=specs,
            selected_features=selected_features,
            replicates_per_particle=replicates_per_particle,
        )
        result = run_smc_abc(
            specs=specs,
            config=smc_config,
            simulate=problem.simulate,
            backend=backend,
            checkpoint=SMCABCCheckpointConfig(
                path=str((checkpoint_dir / f"synthetic_truth_{truth_label}.json").resolve())
                if checkpoint_dir is not None
                else None,
                resume=True,
            ),
            progress_callback=(
                (lambda event, label=truth_label: progress_callback(label, event))
                if progress_callback is not None
                else None
            ),
            run_label=f"synthetic_truth::{truth_label}",
            logger=logger,
        )
        weights = np.asarray([particle.weight for particle in result.particles], dtype=np.float64)
        weights = weights / max(float(np.sum(weights)), 1e-12)
        posterior_mean = posterior_mean_updates(specs=specs, particles=result.particles)
        recovery_row: dict[str, Any] = {
            "label": truth_label,
            "best_distance": float(result.best_particle.distance),
        }
        for spec in specs:
            values = np.asarray([particle.parameters[spec.name] for particle in result.particles], dtype=np.float64)
            q05 = weighted_quantile(values, weights, 0.05)
            q95 = weighted_quantile(values, weights, 0.95)
            truth_value = float(truth_updates[spec.name])
            recovery_row[f"{spec.name}__truth"] = truth_updates[spec.name]
            recovery_row[f"{spec.name}__best"] = result.best_particle.parameters[spec.name]
            recovery_row[f"{spec.name}__posterior_mean"] = posterior_mean[spec.name]
            recovery_row[f"{spec.name}__q05"] = q05
            recovery_row[f"{spec.name}__q95"] = q95
            recovery_row[f"{spec.name}__truth_in_q90"] = bool(q05 <= truth_value <= q95)
        results.append(to_builtin(recovery_row))
        if partial_csv_path is not None:
            import pandas as pd

            atomic_write_dataframe_csv(partial_csv_path, pd.DataFrame(results), index=False)
        if checkpoint_dir is not None:
            atomic_write_json(
                checkpoint_dir / "synthetic_truth_recovery_progress.json",
                {"completed_labels": [str(row["label"]) for row in results]},
            )
    return results
