# CLAIM-001: Single-Family Empirical-Like Statistics Are Recoverable With Assay-Specific Thermodynamic Calibration

- Status: completed
- Created: 2026-05-14
- Supersedes: HYP-001
- Aims: AIM-001
- Linked results: RES-001, RES-002, RES-003, RES-004, RES-005, RES-006, RES-007, RES-011, RES-013

## Claim

The current evidence supports a reduced version of `HYP-001`: Adapt-Env can be
parameterised to recover empirical-like local landscape statistics for a single,
isolated protein family or assay system at a time. The evidence does not support
the stronger original claim that one shared synthetic parameter regime can cover
heterogeneous empirical systems without system-specific calibration.

## Interpretation

The central shift is from universal realism to conditional realism. The model is
scientifically useful when the target family, readout class, and calibration
objective are explicit. Under that framing, the landscape is a controlled
synthetic analogue of one empirical system, not a universal generator of all DMS
landscapes.

This is still a useful result for benchmarking ML-based protein engineering
tasks. A benchmark does not need to reconstruct an assay variant-by-variant, but
it does need defensible local statistics, a clear calibration target, and
transparent limits on external validity. The current results support that
reduced use case: synthetic screening on a calibrated family-specific landscape,
with the fitted system reported as part of the benchmark definition.

## Evidence

- `RES-001`, `RES-003`, `RES-004`, and `RES-005` consistently weaken the
  original shared-regime version of `HYP-001`. Across mixed or homogeneous
  panels, shared calibration failed to recover the empirical summary-statistic
  envelope strongly enough to justify a universal parameter regime.
- `RES-002` shows that the uncalibrated model family can span broad single- and
  double-mutant DFE regimes. This supports model capacity, but not empirical
  adequacy by itself.
- `RES-006` narrows the failure mode: fitting a single assay improves over the
  shared panel but still does not fully rescue the original hypothesis. This is
  evidence that pooling heterogeneous systems was a real confounder, while also
  showing that assay-specific fitting alone is not automatically sufficient.
- `RES-007` provides the strongest positive result for the reduced claim:
  stability-targeted single-assay calibration materially improves recovery on
  SPTN1.
- `RES-011` and `RES-013` independently support the same direction on PHOT by
  showing that stability-readout controls outperform more complex alternatives,
  with `RES-013` giving the cleaner validation-objective comparison.

## Implications

`HYP-001` is superseded rather than simply supported. The project should not
claim universal empirical realism across protein families or assay classes. It
can claim a more defensible system-by-system calibration result: for a defined
family and readout, the model can generate empirical-like local statistics when
the calibration target aligns with the dominant biophysical signal.

This makes the model acceptable as a synthetic benchmark environment only under
that reduced interpretation. Benchmark reports should specify the calibrated
family, readout, objective, and fitted parameter regime. The resulting benchmark
tests ML search under a plausible calibrated landscape, not under a validated
universal biological prior.

## Boundaries

This claim is about local empirical-like statistics, not exact reconstruction of
an empirical assay. It does not certify global optima, universal transfer across
families, or faithful recovery of all higher-order biological mechanisms.
