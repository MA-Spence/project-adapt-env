# EXP-018: ProteinGym OXDA Expression-Activity Paired Latent Trait Diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

Does paired expression plus enzyme activity on OXDA_RHOTO support a two-trait state/function calibration better than a function-only single-trait fit?

## Pre-experiment prediction

If HYP-007 is failing on KRAS/PHOT partly because the tested systems are structurally mismatched, then OXDA_RHOTO expression plus enzymatic activity should show improved activity holdout recovery or readout agreement when expression constrains the state/stability-facing layer and activity constrains a separate function-facing layer. Because the assays contain single mutants only, this experiment will not test double-mutant epistasis recovery.

## Rationale

OXDA_RHOTO provides matched ProteinGym expression and enzyme-activity measurements with weak readout correlation, making it an immediate non-KRAS paired-readout diagnostic for HYP-007. The system is not merely a folding/proteolysis assay, but the absence of multiple mutants limits inference to single-mutant recovery and cross-readout structure.

## Experimental design

- Prepare the matched ProteinGym assays
  `OXDA_RHOTO_Vanella_2023_expression` and
  `OXDA_RHOTO_Vanella_2023_activity`.
- Treat expression as the state/stability-facing readout and enzyme activity
  as the function-facing readout.
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
- `OXDA_RHOTO_Vanella_2023_expression`
- `OXDA_RHOTO_Vanella_2023_activity`

## Configuration

- `config.yaml` fixes the OXDA assay pair, state/function roles, branch
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

- [ ] Config committed
- [ ] Run script committed
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- `RUN-051` completed on `lab-slurm` as scheduler job `194` and produced the
  expected paired-readout and branch-validation outputs.
- `ANA-018` and `RES-018` record the interpretation: paired expression/state
  information did not improve activity recovery over the function-only branch,
  and the fitted reference remained effectively at the peak.
- The experiment remains limited to single-mutant recovery because both OXDA
  ProteinGym assays contain no multiple mutants.
