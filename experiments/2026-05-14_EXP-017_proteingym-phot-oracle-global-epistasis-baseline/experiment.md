# EXP-017: PHOT Oracle Global-Epistasis Baseline

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

Can an oracle additive latent score plus flexible monotone global-epistasis nonlinearity outperform the mechanistic PHOT variants, thereby implicating restrictive AdaptEnv latent priors rather than measurement noise or scalar-assay non-identifiability?

## Pre-experiment prediction

If AdaptEnv's latent priors are too restrictive, a MAVE-NN/MoCHI-style additive latent score with flexible monotone global epistasis should beat the mechanistic variants on held-out single mutants, double mutants, and epistasis prediction. If it does not, the dominant limitation is more likely measurement noise or insufficient identifiability from the scalar PHOT assay.

## Rationale

The oracle baseline separates mechanistic misspecification from assay information limits. It gives the data a flexible additive-to-observed map without requiring the AdaptEnv stability/activity priors to be correct.

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

- MAVE-NN/MoCHI-style oracle baseline metrics on PHOT_CHLRE
- including held-out single/double ranking
- epistasis-prediction metrics
- comparison to RES-010/RES-011 metrics
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
