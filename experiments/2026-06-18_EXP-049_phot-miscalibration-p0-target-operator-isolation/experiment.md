# EXP-049: PHOT miscalibration P0 target-operator isolation

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Are the PHOT target summaries and sparse-Walsh support sampler measured on a common footing, as tested by P0.1 target uncertainty, P0.2 round-trip measurement calibration, and P0.3 support-containment audit?

## Pre-experiment prediction

The local-energy and high-order empirical target components will be unstable under bootstrap or MAVE-NN refit, injected order-3 spectra will be distorted by the empirical measurement operator, and support containment will quantify whether sparse-Walsh supports are visible to sampled dense graphlets.

## Rationale

RES-047 and scratch/RES-047_miscalibration_root_cause_and_isolation_plan.md argue that the spectrum-calibrated sparse-Walsh forward gate may be failing because target summaries, MAVE-NN latent measurement, and dense graphlet sampling are not on a common footing. This experiment pauses new generator engineering and isolates P0.1, P0.2, and P0.3 in one PHOT-only diagnostic.

## Experimental design

- Single PHOT-only isolation run implementing Phase 0 from
  `scratch/RES-047_miscalibration_root_cause_and_isolation_plan.md`.
- P0.1: bootstrap cached PHOT target summaries over dense complete graphlets,
  audit available assay replicate columns, and refit MAVE-NN `latent_phi` with
  five seeds across two regularisation profiles before recomputing the same
  summaries.
- P0.2: synthesize assay-scale PHOT datasets from known injected order spectra,
  run the identical empirical pipeline (MAVE-NN `latent_phi`, dense complete
  graphlet selection, exact summaries), and report recovered-vs-injected
  transfer functions.
- P0.3: reproduce the EXP-048 sparse-Walsh support-selection rule and report
  the fraction of generated order-k terms whose support is fully contained in
  sampled `dense_descriptive` graphlets.
- This experiment deliberately does not repair or refit a generator; it tests
  whether the empirical target, measurement operator, and graphlet sampler are
  commensurate.

## Inputs

- Frozen PHOT source assay, cached MAVE-NN `latent_phi`, and wildtype metadata
  from `RUN-077`.
- EXP-048 / RUN-123 sparse-Walsh support-selection settings and ANA-047 module
  catalog used only for P0.3 support-containment auditing.
- Advisory scratch note:
  `scratch/RES-047_miscalibration_root_cause_and_isolation_plan.md`.

## Configuration

- `config.yaml` fixes the dense graphlet operator, MAVE-NN refit profiles,
  injected synthetic spectra, bootstrap budget, and sparse-Walsh support audit
  settings.

## Execution

- Use `labproj submit` for RUN generation.
- The runner writes a `preflight.json` before executing P0.1, P0.2, and P0.3.

## Expected outputs

- preflight.json; graphlet_catalog.csv; p0_1_target_bootstrap.csv; p0_1_mavenn_refit_metrics.csv; p0_2_roundtrip_transfer.csv; p0_3_support_containment_summary.csv; summary.json

## Analysis plan

- Create an ANA record after execution.

## Completion criteria

- [ ] Config committed
- [ ] Run script committed
- [ ] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
