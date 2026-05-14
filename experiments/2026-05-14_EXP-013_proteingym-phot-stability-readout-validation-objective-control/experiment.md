# EXP-013: PHOT Stability-Readout Validation-Objective Control

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

Does the simpler stability-targeted PHOT control still beat activity-oriented and two-trait formulations when fitted under the same held-out validation objective used in EXP-012?

## Pre-experiment prediction

If the RES-011 improvement mainly reflects a better causal latent mapping rather than objective mismatch, then the stability-readout control should retain strong single- and double-mutant ranking under validation-objective SMC while exposing whether epistasis prediction and reference-to-peak geometry remain unresolved.

## Rationale

RES-011 showed that stability-targeted semantics improved PHOT_CHLRE ranking and KS, but it used the older bootstrap summary-vector SMC objective. RES-010 used the held-out validation objective for the explicit two-trait model. This experiment removes that objective mismatch by fitting the stability-targeted control with the same validation-objective machinery.

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

- Validation-objective SMC posterior
- deterministic stability-control comparisons
- PHOT assay MAVE-NN diagnostics
- branch validations
- posterior summaries
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
