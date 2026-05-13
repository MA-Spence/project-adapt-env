# EXP-007: ProteinGym Single-Assay Stability-Readout SMC-ABC Diagnostic for SPTN1_CHICK

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001

## Question

Does fitting the SPTN1_CHICK_Tsuboyama_2023_1TUD stability assay through a stability-targeted synthetic readout, instead of collapsed synthetic fitness, improve empirical recovery under the same single-assay SMC-ABC workflow?

## Pre-experiment prediction

If the main remaining misspecification is that Adapt-Env collapses stability assays into a generic scalar fitness too early, then a stability-targeted readout with pairwise overrides routed through latent stability should improve single-mutant holdout, double-mutant holdout, and epistasis-prediction metrics relative to the RUN-022 controls on the same assay.

## Rationale

RES-006 showed that removing cross-assay pooling helps but does not rescue the SPTN1 single-assay fit. Code review and local diagnostics now indicate that the current calibration path explains this stability assay mainly through a near-binary folding gate applied to collapsed synthetic fitness, while the functional layer is nearly flat in the wild-type neighborhood. This experiment keeps the exact SPTN1 assay and baseline controls from EXP-006, but adds stability-readout calibration branches that expose latent stability to the observation model and move empirical pairwise structure upstream of the nonlinear gate.

## Experimental design

- Template: generic
- Use the exact single-assay `SPTN1_CHICK_Tsuboyama_2023_1TUD` panel from `EXP-006`.
- Keep the current `baseline_shared_raw` and `predictive_richpair_shared_raw` branches as paired controls.
- Add `baseline_shared_stability_readout` and `predictive_richpair_shared_stability_readout` branches that fit the assay through `synthetic_readout_mode: stability_margin` with `empirical_pairwise_target: stability`.
- Run the same SMC-ABC workflow with the stability-readout calibration options so the Bayesian comparison remains paired to `RUN-022`.
- Compare holdout ranking, double-mutant recovery, epistasis prediction, and reference-to-peak behavior across the control and stability-readout branches.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) paired branch comparisons between the existing RUN-022-style controls and new stability-readout controls on SPTN1_CHICK_Tsuboyama_2023_1TUD; 2) assay-specific mavenn diagnostics; 3) SMC-ABC posterior particles and round diagnostics under the stability-readout path; 4) synthetic-truth recovery metrics for the modified calibration path; and 5) a structured summary.json

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
