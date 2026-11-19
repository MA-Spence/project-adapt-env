# EXP-035: PHOT replicated held-out split graphlet diagnostic

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008

## Question

Is the RUN-079 no-module held-out PHOT graphlet improvement stable across replicated calibration/test graphlet splits, or is it a split artifact?

## Pre-experiment prediction

If the no-module coupled-operator signal is real, calibration-selected no-module candidates should beat fitted-current, smooth-core, empirical-module, random-module, and shuffled-module controls on held-out graphlet distance across most replicated splits and should improve secondary WHT/topology/dynamics components rather than one component only. If the advantage is split-specific, the median and win-rate against fitted current should collapse or be matched by random/shuffled controls.

## Rationale

RUN-079 gave the strongest held-out positive signal for HYP-008, but it was observed on one fixed calibration/test graphlet split and did not validate local empirical modules. This experiment keeps the same PHOT base and candidate families and resamples the split many times to test whether the no-module improvement is robust before any new operator is proposed.

## Experimental design

- Use `PHOT_CHLRE_Chen_2023` as the only empirical system.
- Reuse the fixed `RUN-079` graphlet candidate metrics from
  `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079`;
  do not refit the PHOT base, rerun ABC, or introduce a new operator.
- Preserve the `RUN-079` candidate families: fitted current, smooth core,
  no-module coupled operator, empirical-module coupled operator, random-module
  control, shuffled-module control, uncoupled local-module control,
  residual-only diagnostic, and GNK local null.
- Generate many balanced calibration/test graphlet splits stratified by
  graphlet dimension.
- Within each replicated split, select the best candidate inside each family
  using calibration graphlet distance only.
- Evaluate the calibration-selected candidates on held-out graphlets. The
  primary metric is held-out distance to empirical `latent_phi` graphlet
  metrics. Secondary components are SE, RSE, local maxima, accessibility,
  lag-1 random-walk autocorrelation, WHT high-order fraction, and GFT energy
  fractions.
- Treat the original `RUN-079` split as a reference row, not as the only
  decision point.

## Inputs

- `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079.dvc`
- `RUN-079/graphlet_metrics.csv`
- `RUN-079/graphlet_catalog.csv`
- `RUN-079/parameter_grid.csv`

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- preflight.json
- split_assignments.csv
- split_selected_distances.csv
- split_family_summary.csv
- split_family_pairwise_deltas.csv
- component_summary.csv
- selection_counts.csv

## Analysis plan

- Create an ANA record after execution.
- Primary support would require no-module calibration-selected held-out
  graphlet distance to beat fitted current across most random splits and to
  have a negative median delta versus fitted current.
- A weaker capacity-only finding is possible if no-module sometimes wins but
  has low win rate, high variance, or is matched by random/shuffled modules.
- If random/shuffled modules equal or beat empirical-module and no-module
  families, the result should be interpreted as another weakening of the
  current local-module form of `HYP-008`.

## Completion criteria

- [x] Config recorded
- [x] Run script recorded
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- `RUN-089` completed and was collected into
  `data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089`.
  The output directory is DVC-tracked.
- `ANA-035` and `RES-035` record the interpretation. The no-module
  graphlet-distance improvement from `RUN-079` is stable against fitted-current
  and smooth-core controls across replicated splits, but the specificity
  controls fail because the uncoupled-module control is better overall and
  random/shuffled/empirical module families also improve against current.
- The experiment weakens the current `HYP-008` mechanism-specific
  interpretation while preserving a narrow capacity claim for PHOT graphlet
  distance movement.
