# EXP-019: ProteinGym KCNJ2 Surface-Function Paired Latent Trait Diagnostic

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

Does paired surface trafficking plus ion-conduction function on KCNJ2_MOUSE support a two-trait state/function calibration better than a function-only single-trait fit?

## Pre-experiment prediction

If the multi-latent direction is more appropriate for mechanistically separable non-folding-limited assays, then KCNJ2_MOUSE surface abundance plus channel function should benefit from a state/function decomposition relative to a function-only single-trait fit. Because the assays contain single mutants only, this experiment will not test double-mutant epistasis recovery.

## Rationale

KCNJ2_MOUSE separates surface trafficking from potassium-channel function, giving a cleaner biological split than PHOT scalar fluorescence and a different paired-readout system than KRAS abundance/binding. It is therefore an appropriate single-mutant diagnostic of whether measured state and function traits are non-redundant for HYP-007.

## Experimental design

- Prepare the matched ProteinGym assays
  `KCNJ2_MOUSE_Coyote-Maestas_2022_surface` and
  `KCNJ2_MOUSE_Coyote-Maestas_2022_function`.
- Treat surface trafficking as the state/stability-facing readout and
  ion-conduction function as the function-facing readout.
- Fit three deterministic calibration branches with identical data preparation:
  - `function_only_single_trait`
  - `paired_state_function_single_trait`
  - `paired_state_function_explicit_two_trait`
- In the paired branches, first fit the state readout against a stability-margin
  objective, then fit the function readout with the state-derived updates held
  fixed.
- In the explicit two-trait branch, use the built-in stability gate plus one
  named `readout` trait block, with observed fitness formed as
  `stability_gate * trait:readout:capacity`.
- Because this ProteinGym pair is single-mutant-only, do not interpret this
  run as evidence about double-mutant ranking or epistasis recovery.

## Inputs

- ProteinGym v1.3 reference metadata and substitutions parquet.
- `KCNJ2_MOUSE_Coyote-Maestas_2022_surface`
- `KCNJ2_MOUSE_Coyote-Maestas_2022_function`

## Configuration

- `config.yaml` fixes the KCNJ2 assay pair, state/function roles, branch
  definitions, and Slurm-scale MAVE-NN/calibration settings.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- paired readout correlations
- MAVE-NN assay diagnostics
- branch-level state/function calibration metrics
- and a summary.json

## Analysis plan

- Create an ANA record after execution.

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

- Corrected live execution `RUN-053` completed on `lab-slurm` as scheduler job
  `196` with exit code `0`.
- Outputs were collected under
  `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053`
  and DVC-tracked as
  `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053.dvc`.
- The matched KCNJ2 surface/function readouts were weakly correlated
  (`Spearman 0.290` over `6789` variants). The MAVE-NN observation layer was
  weaker for function (`0.389` test Spearman) than for surface trafficking
  (`0.538` test Spearman).
- The paired state/function branches did not improve function recovery over
  the function-only branch. Function-only holdout Spearman was `0.066`, paired
  single-trait was `0.065`, and paired explicit two-trait was `0.021`; function
  KS was `0.460`, `0.502`, and `0.511` respectively.
- Because both assays contain only single mutants, this run does not test
  double-mutant ranking or epistasis recovery.
- `ANA-019` and `RES-019` record that this KCNJ2 paired state/function
  diagnostic weakens the current `HYP-007` implementation without directly
  updating `HYP-001`.
