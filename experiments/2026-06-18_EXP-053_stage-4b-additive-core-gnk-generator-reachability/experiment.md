# EXP-053: Stage 4B additive-core GNK generator reachability

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009

## Question

Does the locked public additive thermodynamic core plus hashed GNK texture produce a protein-scale reachable descriptor manifold with directional knob maps and required ablation failures?

## Pre-experiment prediction

Increasing GNK weight will raise ruggedness, local maxima, sign epistasis, and high-order fraction while lowering accessibility/autocorrelation; additive-core-only and bare-GNK controls will fail complementary descriptor blocks.

## Rationale

This is Step 4B and extends RUN-127 from small enumerable coordinates to the public AdaptEnv additive_thermodynamic_gnk surface with core-removal, GNK-removal, and routing controls.

## Experimental design

This implements Step 4B from
`scratch/science_plan_part1_ensemble_landscape_validation.md`.

The run evaluates the public AdaptEnv `additive_thermodynamic_gnk` protocol on
protein-scale sequence lengths using sampled binary subcubes. It sweeps
`gnk_epistasis_weight`, stability margin, seed, and support routing, then
records the same descriptor family used for the empirical envelope.

Controls:
- additive-core-only observable (`gnk_epistasis_weight = 0`) must retain the
  DFE/global-readout block but fail specific ruggedness/accessibility texture;
- bare GNK texture must fail the core/readout block;
- structured, shuffled, and random GNK routing are compared at matched global
  knob settings.

## Inputs

- `scratch/science_plan_part1_ensemble_landscape_validation.md`
- `scripts/ensemble_landscape_validation_stage4.py`
- `external/Adapt-Env/adaptenv/landscape.py`

## Configuration

- `config.yaml` preregisters the generator grid, monotonicity expectations, and
  ablation controls.

## Execution

- Use `labproj submit EXP-053 --target lab-slurm`.

## Expected outputs

- generator_parameter_grid.csv
- generator_descriptor_metrics.csv
- generator_knob_monotonicity.csv
- core_texture_ablation_controls.csv
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
