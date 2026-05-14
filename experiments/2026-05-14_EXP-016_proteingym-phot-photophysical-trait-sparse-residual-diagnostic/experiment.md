# EXP-016: PHOT Photophysical-Trait Sparse-Residual Diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

After accounting for stability and a photophysical readout trait, is the remaining PHOT_CHLRE mismatch explained by sparse residual epistasis rather than by the latent mapping itself?

## Pre-experiment prediction

If PHOT fluorescence depends on a brightness or FMN-occupancy trait plus a small set of specific interactions, then adding sparse residual epistasis after the additive latent readout should improve epistasis prediction and double-mutant recovery more than broad extra latent capacity does.

## Rationale

PHOT_CHLRE reports CreiLOV fluorescence rather than a clean activity phenotype. The missing trait may be photophysical brightness or FMN occupancy, and residual pairwise structure may identify specific chromophore-pocket or allosteric interactions rather than generic ruggedness.

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

- PHOT photophysical-trait diagnostic metrics
- residual epistasis fit diagnostics
- held-out single/double ranking
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
