# EXP-014: PHOT Explicit Two-Trait Validation-Objective Replicate

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

Does a direct replicate of the current explicit stability-plus-readout trait model reproduce the RES-010 failure pattern under the validation objective?

## Pre-experiment prediction

If RES-010 reflects structural mismatch rather than stochastic SMC noise, then a new run of the current explicit two-trait formulation should again show weak single- and double-mutant recovery, near-peak geometry, and posterior pressure on empirical pairwise strength.

## Rationale

RES-010 is the first clean explicit two-latent PHOT test, but interpretation is stronger if the failure is reproducible under a new experiment entry paired with the stability-objective control and downstream mismatch diagnostics.

## Experimental design

- Template: generic
- Configure the experiment in config.yaml.
- Run via run.sh or labproj submit.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- Validation-objective SMC posterior for the current explicit two-trait model
- deterministic two-trait controls
- PHOT MAVE-NN diagnostics
- posterior summaries
- branch validations
- and summary.json.

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
