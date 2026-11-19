# EXP-027: TEM-1 fixed-scale residual ratio sweep graphlet diagnostic

## Status

- completed: `RUN-069` completed on `lab-slurm` with scheduler exit code `0`
  on 2026-05-15. Outputs are collected locally and DVC-tracked.

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005

## Question

Does scanning fixed residual-to-additive scale ratios reveal a threshold where high-order residual variance becomes exposed as adaptive-walk additive gains diminish, without using graphlet-local DFE normalization?

## Pre-experiment prediction

If the corrected relative-scale mechanism is present, increasing fixed residual scale should monotonically increase residual exposure, WHT high-order fraction, SE/RSE, local maxima, reduced accessibility, and shorter random-walk autocorrelation, with stronger effects in later walk stages where available beneficial effects are smaller. Smooth-core, pairwise, exact-match, placement, and GNK controls distinguish this from direct local scaling or non-procedural ruggedness.

## Rationale

RES-026 showed that the locally DFE-scaled contact residual creates topology but is not the intended fixed residual/additive mechanism. This repeat removes the local-DFE primary calibration setting, scans fixed residual scale values directly, and reduces evaluated graphlets by the scale-grid size so the primary scale-sweep workload remains comparable to the previous single-scale diagnostic.

## Experimental design

- Repeat the `EXP-026` dense graphlet diagnostic with the same TEM-1 reference,
  smooth thermodynamic core, start classes, adaptive-walk sampling scheme,
  graphlet dimensions, contact/random graphlet sources, output-reduction
  settings, Walsh/GFT order summaries, topology metrics, accessibility metrics,
  and random-walk autocorrelation metrics.
- Remove the local-DFE-scaled primary residual setting that confounded
  `RUN-067`. The primary contact residual grid uses only
  `residual_scale_mode = constant_walk_initial_beneficial_effect`.
- Scan fixed residual/additive scale ratios using
  `residual_scale_relative_to_local_dfe_sd = 1.0, 4.0, 16.0`. In this fixed
  mode the multiplier is applied to the walk-initial mean available beneficial
  effect, not to the current graphlet-local DFE scale.
- Reduce evaluated graphlets by the grid size: `EXP-026` used `90` graphlets;
  this experiment uses `30 = 90 / 3` graphlets. This preserves the primary
  fixed-scale sweep workload while allowing the same statistics to be reported
  for each scale value.
- Keep the major controls: smooth core, RES-024 exact-match haplotype control,
  RES-023 DCA pair override, shuffled-contact residual, random-position
  residual, GNK local null, and residual-only diagnostic. The residual controls
  use the same fixed scale grid as the primary model.
- Keep graphlet output compact: pair-row output is disabled and only the top
  Walsh/GFT terms are recorded, while all order-level WHT/GFT summaries remain
  in the metric rows.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- `config.yaml` records the scale-grid design explicitly under
  `scale_grid_design`, including the baseline graphlet budget, grid size,
  fixed scaling mode, scale values, and reduced graphlet budget.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- module_parameter_sweep.csv
- graphlet_catalog.csv
- calibration_graphlet_metrics.csv
- test_graphlet_metrics.csv
- graphlet_walsh_terms.csv
- graphlet_gft_coefficients.csv
- graphlet_accessibility.csv

## Analysis plan

- Create an ANA record after execution using the same review statistics as
  `ANA-026`: sampling audit, diminishing-return proxies by walk stage,
  calibration/test transfer, scale-wise model comparison, start-class
  stratification, walk-stage relative-scale summaries, Walsh/GFT order
  summaries, accessibility, local maxima/basin metrics, SE/RSE, and
  random-walk autocorrelation.
- The central test is directional. For the fixed-scale primary residual,
  increasing scale should increase residual exposure, high-order WHT fraction,
  SE/RSE, local maxima, lower accessibility, and shorter autocorrelation. The
  strongest evidence for the intended mechanism would be a larger effect in
  late graphlets where available beneficial effects are smaller.
- A result should be interpreted as negative for the proposed mechanism if only
  local-scale controls or GNK produce topology, or if scale increases topology
  only by making high-order variance dominant and non-empirical-looking.

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- 2026-05-15: Local quick smoke completed with synthetic TEM-1-like alignment.
  The smoke produced `15` graphlets, `300` graphlet/model metric rows,
  `2400` capped Walsh rows, and `2400` capped GFT rows. Graphlets were balanced
  across five start classes (`3` each) and walk stages (`5` early, `5` middle,
  `5` late), confirming the reduced graphlet cap still samples the intended
  adaptive-walk axis.
- 2026-05-15: `RUN-068` dry-run submission created. `RUN-069` submitted live to
  `lab-slurm` as scheduler job `228`.
- 2026-05-15: `labproj run reconcile RUN-069` updated the run to `completed`
  with scheduler exit code `0`. A collection dry-run was attempted, but local
  collection is blocked until `/Volumes/DutchZFS/labproj` is mounted.
- 2026-05-15: `RUN-069` outputs were subsequently collected locally and
  DVC-tracked at
  `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069.dvc`.
- 2026-05-15: `ANA-027` records the analysis and `RES-027` records the result.
  The fixed residual scale knob is directional, but the late-walk/peak-stage
  unmasking interpretation remains weak.
