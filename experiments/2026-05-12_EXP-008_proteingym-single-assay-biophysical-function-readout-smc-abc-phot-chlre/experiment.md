# EXP-008: ProteinGym Single-Assay Biophysical-Function-Readout SMC-ABC Diagnostic for PHOT_CHLRE

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001, HYP-007

## Question

Does replacing collapsed synthetic fitness with an explicit stability-plus-function biophysical readout improve empirical recovery for the multi-mutant `PHOT_CHLRE_Chen_2023` functional assay under the same single-assay SMC-ABC workflow used in the most recent stability round?

## Pre-experiment prediction

If the misspecification exposed by `RUN-024` extends beyond unfolding assays, then a functional assay fit through an explicit stability-plus-function readout, with empirical pairwise overrides routed through the function layer, should improve single-mutant holdout, double-mutant holdout, and epistasis-prediction metrics relative to the paired raw-fitness controls on the same assay.

## Rationale

`RUN-024` showed that the recent single-assay improvement came from exposing the right latent variable to calibration rather than from changing the optimizer alone. `HYP-007` makes the stronger claim that functional assays should not be fit through the same collapsed scalar used for viability-style readouts. This experiment therefore keeps the exact single-assay Bayesian workflow, control branches, and resource profile from `EXP-007`, but moves to one activity assay with substantial multi-mutant coverage and adds calibration branches that treat assay signal as a stability-plus-function biophysical readout instead of collapsed synthetic fitness.

## Experimental design

- Template: generic
- Use a single ProteinGym functional assay with strong multi-mutant coverage: `PHOT_CHLRE_Chen_2023`.
- Keep the current `baseline_shared_raw` and `predictive_richpair_shared_raw` branches as paired controls.
- Add `baseline_shared_biophysical_function_readout` and `predictive_richpair_shared_biophysical_function_readout` branches that fit the assay through `synthetic_readout_mode: stability_function` with `empirical_pairwise_target: function`.
- Fit extra readout parameters that control the stability midpoint, stability slope, and function exponent of the explicit biophysical readout.
- Run the same SMC-ABC workflow with the biophysical-function calibration options so the Bayesian comparison remains directly paired to the recent single-assay stability workflow.
- Compare holdout ranking, double-mutant recovery, epistasis prediction, and reference-to-peak behavior across the raw-fitness controls and the new biophysical-function branches.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) paired branch comparisons between the existing raw-fitness controls and new biophysical-function controls on `PHOT_CHLRE_Chen_2023`; 2) assay-specific mavenn diagnostics; 3) SMC-ABC posterior particles and round diagnostics under the biophysical-function path; 4) synthetic-truth recovery metrics for the modified calibration path; and 5) a structured `summary.json`

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
