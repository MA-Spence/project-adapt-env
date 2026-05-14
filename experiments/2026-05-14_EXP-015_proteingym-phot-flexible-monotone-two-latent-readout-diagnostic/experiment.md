# EXP-015: PHOT Flexible Monotone Two-Latent Readout Diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

Does replacing the forced stability-gate times readout-capacity collapse with a flexible monotone two-latent readout surface improve PHOT_CHLRE recovery?

## Pre-experiment prediction

If the structural mismatch is mainly the imposed product/gate readout, then a constrained monotone readout over stability and readout-trait scores should improve held-out ranking and epistasis prediction without requiring a near-peak reference solution.

## Rationale

RES-010 tested the explicit two-trait implementation but forced the public scalar to stability_gate * trait:readout:capacity. RES-011 suggests stability alone is a strong predictive surrogate. This experiment isolates whether the failure is the hard-coded nonlinear collapse rather than the existence of two latent traits.

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

- PHOT panel preparation
- fitted monotone two-latent readout diagnostic metrics
- held-out single/double ranking metrics
- epistasis-prediction metrics
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
