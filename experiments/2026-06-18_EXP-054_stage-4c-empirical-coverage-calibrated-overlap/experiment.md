# EXP-054: Stage 4C empirical coverage and calibrated overlap

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009

## Question

Does the generator reachable manifold cover the empirical descriptor envelope, and what posterior-like knob regions overlap each empirical assay band without point-fitting?

## Pre-experiment prediction

At least one region of full public-model knob space will overlap each core Tier-1 empirical band on the registered descriptor set, while additive-core-only and bare-GNK controls will miss complementary blocks.

## Rationale

This is Step 4C: coverage plus calibrated overlap as bands/posterior regions, not one-assay inverse fitting.

## Experimental design

This implements Step 4C from
`scratch/science_plan_part1_ensemble_landscape_validation.md`.

The run recomputes the empirical descriptor bands and generator manifold through
the same frozen code path, then asks whether empirical bands are covered by the
reachable generator manifold. It reports posterior-like knob regions as
descriptor-distance bands, not point fits to individual assays.

Controls:
- additive-core-only and bare-GNK controls are scored through the same overlap
  tables as the full public model;
- routing controls are kept matched on the global knob grid;
- sparse ProteinGym entries contribute only descriptors allowed by the scope
  audit.

## Inputs

- `scratch/science_plan_part1_ensemble_landscape_validation.md`
- `scripts/ensemble_landscape_validation_stage4.py`
- `data/processed/ensemble-landscape-validation-datasets.dvc`
- `external/Adapt-Env/adaptenv/landscape.py`

## Configuration

- `config.yaml` defines the overlap metric set and generator grid.

## Execution

- Use `labproj submit EXP-054 --target lab-slurm`.

## Expected outputs

- coverage_metrics.csv
- posterior_knob_regions.csv
- control_overlap_failures.csv
- empirical_bands.csv
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
