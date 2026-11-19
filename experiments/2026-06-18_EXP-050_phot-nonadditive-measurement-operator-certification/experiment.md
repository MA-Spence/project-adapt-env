# EXP-050: PHOT non-additive measurement operator certification

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-006, HYP-009
- Previous experiments: EXP-049
- Motivating results: RES-048
- Scratch source: scratch/RES-048_forward_plan_tunable_benchmark_generator.md

## Question

Can non-additive MAVE-NN GP-map operators recover injected high-order and local-ruggedness spectra under the same dense complete-graphlet measurement path used for PHOT?

## Pre-experiment prediction

Pairwise and black-box/GE GP-map profiles will recover injected high-order fraction, sign epistasis, local maxima, and accessibility substantially better than the additive operator, but residual smoothing may remain.

## Rationale

RES-048 showed strong additive-map smoothing and exposed that the prior PHOT target was not measured under a single certified operator. This WS-A experiment re-runs the round-trip certification with non-additive GP-map profiles before any empirical overlap or benchmark-region claim is attempted.

## Experimental design

This is WS-A from `scratch/RES-048_forward_plan_tunable_benchmark_generator.md`.
It reuses the RES-048/RUN-125 round-trip harness, but replaces the additive
MAVE-NN GP-map with non-additive operator profiles.

The run builds the same dense PHOT complete graphlets from frozen RUN-077 source
data, injects three known synthetic spectra, fits each spectrum through two
non-additive MAVE-NN profiles, and compares direct injected statistics with
recovered latent-operator statistics.

Operator profiles:

- `pairwise_low_regularization`: MAVE-NN pairwise GP-map. For full-length
  PHOT this is subject to a preregistered parameter-count preflight; if the
  pairwise tensor is infeasible, the run records `fit_skipped_preflight` rather
  than allowing an OOM kill.
- `blackbox_ge_neural`: MAVE-NN black-box/GE profile.

Injected spectra:

- `pure_order3`
- `order2_order3_mixture`
- `low_order_control`

Pass criterion: for high-order fraction, local sign epistasis, local maxima, and
accessibility, recovered/injected ratios must fall in `[0.7, 1.4]` across all
three injected spectra. If no operator is identity-like, the transfer function is
still recorded and must be applied symmetrically in downstream empirical/model
comparisons.

## Inputs

- Frozen PHOT source output: `RUN-077`
- Scratch plan: `scratch/RES-048_forward_plan_tunable_benchmark_generator.md`
- Round-trip harness: `scripts/phot_nonadditive_operator_certification.py`

## Configuration

- `config.yaml` preregisters graphlet selection, injected spectra, operator
  profiles, and certification tolerances.

## Execution

- Completed run: `RUN-129` on `lab-slurm` as Slurm job `1937`, inferred
  scheduler exit code `0`.
- Failed/canceled execution context: `RUN-126` failed by OOM before the
  corrected preflight behavior; `RUN-128` was canceled before execution.

## Expected outputs

- preflight.json
- graphlet_catalog.csv
- operator_profile_catalog.csv
- injected_terms.csv
- roundtrip_graphlet_metrics.csv
- roundtrip_transfer.csv
- mavenn_metrics.csv
- operator_certification_summary.csv
- summary.json

## Analysis plan

- Completed analysis: `ANA-056`
- Completed result: `RES-056`

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- `RUN-129` produced all expected tables, but no non-additive operator was
  certified. The pairwise profile was preflight-skipped for all three spectra
  at PHOT length `118`, and the black-box/GE profile failed for all three
  spectra with an installed MAVE-NN/Keras `MultilayerPerceptronGPMap build()`
  interface error.
- `RES-056` records this as a negative certification result for the current
  non-additive PHOT measurement route, not as a recovered-landscape comparison
  or a biological PHOT ruggedness measurement.
