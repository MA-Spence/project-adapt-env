# EXP-056: Stage 4E held-out descriptor generalisation

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009

## Question

Do knob regions calibrated on a subset of Tier-1 landscapes predict held-out landscape descriptors within empirical bands?

## Pre-experiment prediction

The full additive-core plus GNK model will generalise better than core-only or bare-GNK controls, but failures should identify missing mechanisms rather than be forced into a universal-realism claim.

## Rationale

This is Step 4E and turns coverage into out-of-sample descriptor prediction across held-out empirical landscapes.

## Experimental design

This implements Step 4E from
`scratch/science_plan_part1_ensemble_landscape_validation.md`.

The run performs leave-one-Tier-1-landscape-out descriptor prediction. For each
held-out empirical landscape, it selects a public-model knob region using the
remaining Tier-1 landscapes and then reports whether the held-out descriptor
set falls inside the empirical bands. The test is out-of-sample descriptor
generalisation, not universal variant-level reconstruction.

Controls:
- the same additive-core-only, bare-GNK, and routing-control regimes are present
  in the generator grid;
- higher-order/topographic descriptors are never computed from sparse
  ProteinGym-only systems;
- failures are retained as missing-mechanism evidence rather than forced into a
  realism claim.

## Inputs

- `scratch/science_plan_part1_ensemble_landscape_validation.md`
- `scripts/ensemble_landscape_validation_stage4.py`
- `data/processed/ensemble-landscape-validation-datasets.dvc`
- `external/Adapt-Env/adaptenv/landscape.py`

## Configuration

- `config.yaml` defines the held-out Tier-1 panel, overlap metrics, and
  generator grid.

## Execution

- Use `labproj submit EXP-056 --target lab-slurm`.

## Expected outputs

- heldout_predictions.csv
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
