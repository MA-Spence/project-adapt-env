# EXP-006: ProteinGym Single-Assay SMC-ABC Calibration Diagnostic for SPTN1_CHICK

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001

## Question

Does the Bayesian SMC-ABC calibration from EXP-005 recover one empirical assay from the paired six-assay panel when pooling is removed, using SPTN1_CHICK_Tsuboyama_2023_1TUD as a single-dataset diagnostic?

## Pre-experiment prediction

If pooling across assays is a major bottleneck, then fitting the same Bayesian calibration path to only SPTN1_CHICK_Tsuboyama_2023_1TUD should materially improve single-mutant and double-mutant holdout recovery relative to the shared six-assay fit while retaining nonzero epistatic structure.

## Rationale

RES-003 and RES-004 already weakened the claim that pooling alone explains the calibration failure, but both relied on deterministic fitters. RES-005 then showed that the Bayesian fitter can recover structured parameters on model-matched synthetic targets and no longer fails by trivial collapse, yet the shared six-assay empirical fit remained weak. This experiment isolates whether removing cross-assay pooling on the strongest per-assay candidate from the paired panel makes the empirical inverse problem materially easier under the same SMC-ABC path.

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

- 1) a fixed single-assay ProteinGym target for SPTN1_CHICK_Tsuboyama_2023_1TUD; 2) assay-specific mavenn diagnostics; 3) deterministic branch comparison metrics; 4) SMC-ABC posterior particles and round diagnostics; 5) synthetic-truth recovery metrics; and 6) a structured summary.json

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
