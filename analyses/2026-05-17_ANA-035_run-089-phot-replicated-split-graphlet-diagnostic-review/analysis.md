# ANA-035: RUN-089 PHOT replicated split graphlet diagnostic review

## Purpose

Analyze whether `RUN-079`'s no-module held-out PHOT graphlet improvement is
stable across replicated calibration/test graphlet splits and whether it
survives fitted-current, smooth-core, empirical-module, random-module,
shuffled-module, and uncoupled-module controls.

## Linked experiments/runs

- Experiments: `EXP-035`
- Runs: `RUN-089`
- Source data:
  `data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089`
- DVC pointer:
  `data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089.dvc`

## Scientific context

`RUN-079` was the strongest positive PHOT coupled-operator signal for
`HYP-008`: the no-module candidate beat fitted-current and smooth-core controls
on the original held-out graphlet split. The weakness was that this was a
single graphlet split and the empirical-module mechanism did not validate
against no-module/random/shuffled controls. `EXP-035` therefore did not test a
new operator. It resampled the calibration/test graphlet split many times using
the same `PHOT_CHLRE_Chen_2023` base and the same candidate families from
`RUN-079`.

Support for `HYP-008` would require more than a lower graphlet distance against
the fitted current model. The no-module or empirical local-module family would
need a stable held-out advantage and the advantage would need to be specific
relative to random, shuffled, and uncoupled controls. A broad improvement shared
by non-empirical controls is capacity evidence, not evidence for the proposed
local non-diffuse mechanism.

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-089_preflight_and_output_audit.md`
- `tables/run-089_family_summary.md`
- `tables/run-089_pairwise_control_deltas.md`
- `tables/run-089_original_split_reference.md`
- `tables/run-089_component_summary.md`
- `tables/run-089_top_selection_counts.md`
- `tables/run-089_hyp008_decision_table.md`

## Run and data audit

`RUN-089` completed and was collected with `labproj collect`. The collected
output directory contains the planned 8 runtime output files and is tracked by
DVC with pointer
`data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089.dvc`
(`md5 e8e1e0fda6b45f489c0e584ebba888dc.dir`, 8 files, 11,739,787 bytes).

The run used 512 random balanced calibration/test graphlet splits plus the
original `RUN-079` split as a reference row. Each random split used 24 graphlets
stratified by graphlet dimension, 50% calibration graphlets, 50% held-out test
graphlets, 3,624 graphlet metric rows, and 144 candidate parameter rows. Within
each candidate family, selection was by calibration graphlet distance only.

## Key observations

- The no-module improvement is stable against the frozen current-model
  controls. Across 512 random splits, no-module median held-out graphlet
  distance is `0.4839`, versus `0.5232` for fitted current and `0.5358` for
  smooth core. It beats fitted current in `93.8%` of splits and smooth core in
  `97.7%` of splits; median deltas are `-0.0468` and `-0.0577`.
- The no-module result does not beat the single uncoupled-module landscape
  control. The uncoupled-module control has lower median held-out distance
  (`0.4469`), lower mean distance (`0.4907`), is best in `51.2%` of splits, and
  beats no-module in `63.5%` of random splits.
- The improvement is not specific to empirical local modules. Random,
  shuffled, and empirical module families also beat fitted current often:
  random-module median distance `0.4793` with `86.1%` win rate versus current;
  shuffled-module `0.4999` with `82.6%`; empirical-module `0.4943` with
  `84.4%`. Empirical-module performance is not reliably better than shuffled
  module performance.
- The original `RUN-079` split was unusually favorable to no-module. On that
  fixed split, no-module test distance is `0.3137`, far below fitted current
  `0.4016` and smooth core `0.4119`. The replicated median no-module distance
  is `0.4839`, so the original split overstated the size and specificity of
  the no-module advantage.
- The secondary component package is still mismatched. No-module modestly
  improves several components relative to fitted current, including
  accessibility component distance (`0.7599` versus `0.9537`) and WHT
  high-order component distance (`0.5532` versus `0.5938`). But it still
  underestimates empirical accessibility (`0.0525` versus `0.2181`), local
  maxima (`4.08` versus `5.17`), and lag-1 autocorrelation (`0.0607` versus
  `0.2031`), while overestimating WHT high-order fraction (`0.5131` versus
  `0.3583`).
- GNK and residual-only diagnostics remain poor controls for this specific
  target. GNK median held-out distance is `0.7608`; residual-only median
  distance is `2.7185`.

## Interpretation

`RUN-089` answers the registered diagnostic question directly. The no-module
held-out graphlet gain from `RUN-079` is not merely a one-split artifact when
the comparison is fitted current or smooth core. Under that narrow metric, the
candidate family is improving the model on held-out graphlet data.

The result does not support the local non-diffuse mechanism in `HYP-008`.
Specificity controls fail: an uncoupled-module control is better overall, and
random/shuffled/empirical module variants all show broadly similar gains
against the fitted current model. This is evidence that the `RUN-079` candidate
set contains enough capacity or calibration-selection flexibility to move the
PHOT graphlet distance, not evidence that empirical local modules or the
no-module nonlinear coupling capture the proposed mechanism.

The best description is therefore mixed but mostly negative for the hypothesis:
`RUN-089` confirms a stable no-module graphlet-distance improvement against
the fitted current and smooth-core baselines, while weakening `HYP-008` because
the improvement fails the uncoupled, random-module, shuffled-module, and
secondary-component specificity checks.

## Limitations

- The diagnostic uses only `PHOT_CHLRE_Chen_2023`.
- The empirical target remains the processed `latent_phi` graphlet layer.
- The experiment reuses `RUN-079` candidate metrics and does not test a joint
  refit or a new operator.
- Candidate selection is deliberately graphlet-split based; it does not add an
  assay guardrail or an independent biological validation layer.
- The result weakens the current `RUN-079` local-module/no-module
  interpretation, not every possible nonlinear local trait operator.
