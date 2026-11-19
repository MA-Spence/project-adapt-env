# EXP-051: Hashed GNK benchmark generator knob-map validation

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009
- Previous experiments: EXP-048, EXP-049
- Motivating results: RES-047, RES-048
- Scratch source: scratch/RES-048_forward_plan_tunable_benchmark_generator.md

## Question

Do hashed GNK/RMF-style epistatic-texture knobs, added inside the retained stability-threshold core and readout, set local maxima, accessibility, Walsh spectrum, and sign/reciprocal epistasis in the predicted directions?

## Pre-experiment prediction

The core-only condition will retain neutrality and diminishing-returns structure, while increasing hashed GNK epistasis weight will increase local maxima and high-order/sign epistasis and reduce accessibility; structured and shuffled routing should be matched for global statistics but differ in support localization.

## Rationale

RES-047 and RES-048 argue against another PHOT inverse-fit loop. This WS-B experiment replaces the failing residual epistasis module with a lazy hashed GNK texture layer and validates knob-to-statistic maps directly on enumerable model coordinates.

## Experimental design

This is WS-B from `scratch/RES-048_forward_plan_tunable_benchmark_generator.md`.
It implements a direct enumerable prototype of the proposed composite generator:

- retain an additive stability/Fuji backbone;
- retain the saturating readout that creates neutral plateaus and global
  diminishing-returns epistasis;
- replace the failing residual epistasis module with a deterministic hashed GNK
  texture layer whose support is structured, shuffled, or random.

The experiment enumerates binary small-L coordinates and directly measures
knob-to-statistic curves. The main knob is `epi_weight`, which plays the
RMF/GNK interaction-amplitude role. The validation metrics include local maxima,
accessibility, Walsh high-order fraction, sign/reciprocal epistasis, random-walk
autocorrelation, hamming-layer beneficial fractions, and threshold-stratified
pairwise epistasis.

The core guardrail is explicit: `epi_weight=0` must still show the stability
core/readout behaviour rather than collapsing to a bare GNK/RMF model. Nonzero
`epi_weight` then tests the hashed GNK texture and routing controls.

## Inputs

- Scratch plan: `scratch/RES-048_forward_plan_tunable_benchmark_generator.md`
- Runner: `scripts/hashed_gnk_benchmark_knob_validation.py`

## Configuration

- `config.yaml` preregisters support modes, seeds, epistasis-weight grid,
  core/readout parameters, order spectrum, and expected monotonic directions.

## Execution

- Use `labproj submit EXP-051 --target lab-slurm`.

## Expected outputs

- preflight.json
- parameter_grid.csv
- graphlet_metrics.csv
- knob_summary.csv
- monotonicity_checks.csv
- support_catalog.csv
- hamming_layer_metrics.csv
- threshold_epistasis_metrics.csv
- core_readout_guardrail_checks.csv
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
