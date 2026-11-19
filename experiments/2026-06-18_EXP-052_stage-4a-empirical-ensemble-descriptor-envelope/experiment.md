# EXP-052: Stage 4A empirical ensemble descriptor envelope

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009

## Question

What operator-free descriptor envelope is observed across the current Tier-1 combinatorial landscape library, with ProteinGym restricted to order-1/global summaries?

## Pre-experiment prediction

Tier-1 landscapes will show broad but structured bands for DFE, ruggedness, accessibility, sign epistasis, and order spectrum; sparse ProteinGym assays will provide only order-1/global summaries under the scope audit.

## Rationale

This is Step 4A from scratch/science_plan_part1_ensemble_landscape_validation.md and creates the empirical reference distribution before any generator overlap claim.

## Experimental design

This implements Step 4A from
`scratch/science_plan_part1_ensemble_landscape_validation.md`.

The run applies one frozen descriptor pipeline to the current ensemble dataset
bundle under `data/processed/ensemble-landscape-validation-datasets`.
Tier-1 complete or near-complete state graphs are eligible for operator-free
higher-order, ruggedness, accessibility, Walsh/spectral, local-maxima, and
sign/reciprocal-epistasis descriptors. Sparse ProteinGym systems are explicitly
restricted to order-1 or marginal fitness summaries; the run writes
`descriptor_scope_audit.csv` so this restriction is machine-checkable.

Controls:
- scope audit blocks higher-order claims on sparse single-substitution DMS;
- full binary landscapes are measured as full cubes when complete;
- multiallelic complete landscapes are measured through complete sampled binary
  subcubes, never through an inferred MAVE-NN or other measurement operator.

## Inputs

- `scratch/science_plan_part1_ensemble_landscape_validation.md`
- `scripts/ensemble_landscape_validation_stage4.py`
- `data/processed/ensemble-landscape-validation-datasets.dvc`

## Configuration

- `config.yaml` lists the Tier-1 and ProteinGym/Tier-3 datasets, graphlet
  sampling limits, and the descriptor set used downstream by Steps 4C-4E.

## Execution

- Use `labproj submit EXP-052 --target lab-slurm`.

## Expected outputs

- descriptor_metrics.csv
- graphlet_metrics.csv
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
