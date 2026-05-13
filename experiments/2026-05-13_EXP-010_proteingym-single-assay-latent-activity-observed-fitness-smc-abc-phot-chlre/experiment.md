# EXP-010: ProteinGym Single-Assay Latent-Activity-Observed-Fitness SMC-ABC Diagnostic for PHOT_CHLRE

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001, HYP-007

## Question

Does fitting the PHOT_CHLRE multi-mutant activity assay through the new observed-fitness composition path, with scalar fitness formed as a latent activity readout rather than a collapsed post hoc score, improve recovery under the single-assay SMC-ABC workflow?

## Pre-experiment prediction

If RES-008 failed partly because the previous activity model still treated the assay readout as an auxiliary calibration view over a legacy collapsed generator, then promoting activity to the public observed-fitness readout and routing epistasis through the latent function layer should improve holdout ranking, epistasis prediction, and reference-to-peak behavior on PHOT_CHLRE relative to the earlier experiment.

## Rationale

RES-007 showed that calibration improves when the fitted readout matches the latent causal quantity. RES-008 then showed that the first PHOT_CHLRE stability-plus-function implementation did not outperform the raw control, which suggests the remaining mismatch may be structural rather than merely parametric. The new Adapt-Env observed-fitness composition path now allows the scalar benchmark fitness itself to be generated from latent molecular components instead of leaving activity as a special secondary readout layered onto the historical collapsed score. This experiment returns to the same multi-mutant PHOT_CHLRE assay used in EXP-008 and tests whether making activity the actual observed fitness, with epistasis routed through the latent function layer before collapse, improves empirical recovery under otherwise matched single-assay calibration.

## Experimental design

- Template: generic
- Use the same single ProteinGym activity assay with strong multi-mutant coverage
  as `EXP-008`: `PHOT_CHLRE_Chen_2023`.
- Keep the same single-assay SMC-ABC workflow, resource profile, summary
  features, and synthetic-truth recovery structure used in the previous PHOT
  round.
- Replace the legacy collapsed observed fitness with the new generalized
  observed-fitness composition path:
  - `observed_fitness_combine_mode: product`
  - `observed_fitness_terms: [stability_gate, function_capacity]`
- Route both generic and empirical pairwise epistasis through the latent
  function layer before collapse into scalar observed fitness.
- Use two deterministic controls under that generator:
  - `baseline_shared_activity_readout`
  - `predictive_richpair_shared_activity_readout`
- Run one Bayesian shared fit on the same generator with
  `synthetic_readout_mode: fitness`, so calibration targets the public scalar
  activity readout directly.
- Interpret the historical comparison against `EXP-008` / `RES-008`, since the
  current Bayesian runner supports one base generator per experiment.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) PHOT_CHLRE deterministic and Bayesian validation metrics under the
  latent-activity observed-fitness generator; 2) assay-specific `mavenn`
  diagnostics; 3) posterior particles, round summaries, and synthetic-truth
  recovery for the new parameterization; and 4) a structured `summary.json`

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
