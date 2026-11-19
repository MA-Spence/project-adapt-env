# EXP-055: Stage 4D joint multi-descriptor realism gate

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009

## Question

Can one public-model knob setting jointly match multi-peakedness, navigability, order spectrum, sign epistasis, and DFE bands on at least two Tier-1 assays?

## Pre-experiment prediction

Intermediate GNK weight near the RUN-127 joint regime will be the best candidate, but the test can fail if the empirical rugged-yet-navigable joint lies outside the reachable manifold.

## Rationale

This is Step 4D, the hard joint-realism test from the science plan, using TrpB/DHFR/GB1-style Tier-1 anchors and explicit ablation controls.

## Experimental design

This implements Step 4D from
`scratch/science_plan_part1_ensemble_landscape_validation.md`.

The run applies the hard joint-realism gate: a single public-model knob setting
must jointly overlap empirical bands for the registered DFE, multi-peakedness,
navigability/accessibility, order-spectrum, sign-epistasis, and
scale-normalized magnitude-epistasis descriptors on at least two Tier-1 assays.
The headline empirical anchors are GB1, TrpB, DHFR, and CR9114-family binary
landscapes.

Controls:
- additive-core-only controls are expected to fail ruggedness/accessibility;
- bare-GNK controls are expected to fail the core/readout DFE/global block;
- structured, shuffled, and random routing are retained so joint realism is not
  conflated with one routing identity.

## Inputs

- `scratch/science_plan_part1_ensemble_landscape_validation.md`
- `scripts/ensemble_landscape_validation_stage4.py`
- `data/processed/ensemble-landscape-validation-datasets.dvc`
- `external/Adapt-Env/adaptenv/landscape.py`

## Configuration

- `config.yaml` defines the Tier-1 anchors, joint descriptor set, minimum
  fraction-in-band rule, and minimum number of assays hit.

## Execution

- Use `labproj submit EXP-055 --target lab-slurm`.

## Expected outputs

- joint_realism_coverage.csv
- joint_realism_candidates.csv
- posterior_knob_regions.csv
- descriptor_scope_audit.csv
- summary.json

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
