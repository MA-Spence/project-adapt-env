# ANA-036: RUN-091 PHOT coupled family discriminating ABC review

## Purpose

Analyze whether `RUN-091`'s ABC refit rescues PHOT coupled-family held-out
graphlet performance while preserving specificity against no-module,
random-module, shuffled-module, and uncoupled-module controls.

## Linked experiments/runs

- Experiments: `EXP-036`
- Runs: `RUN-091`
- Source data:
  `data/processed/phot-coupled-family-discriminating-abc/RUN-091`
- DVC pointer:
  `data/processed/phot-coupled-family-discriminating-abc/RUN-091.dvc`

## Scientific context

`EXP-036` was registered after `RES-035` showed that the earlier no-module
held-out graphlet improvement was stable against fitted-current and smooth-core
baselines but failed specificity against uncoupled, random, and shuffled
controls. The follow-up question was not whether to add a new operator. It was
whether SMC-ABC fitting of the same candidate families could recover a
specific empirical-module advantage on the held-out PHOT graphlet split.

Support for `HYP-008` would require the calibration-selected empirical-module
family to beat fitted current, smooth core, no-module, random-module,
shuffled-module, and uncoupled-module controls on held-out graphlet distance,
while preserving component guardrails and showing non-diffuse module-local
energy. A gain that is matched by random modules or beaten by an uncoupled
landscape is capacity evidence, not support for the proposed local
non-diffuse mechanism.

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-091_preflight_and_collection_audit.md`
- `tables/run-091_selected_family_distances.md`
- `tables/run-091_pairwise_control_deltas.md`
- `tables/run-091_component_guardrails.md`
- `tables/run-091_component_distance_summary.md`
- `tables/run-091_raw_test_metric_means.md`
- `tables/run-091_module_energy_localization.md`
- `tables/run-091_hyp008_decision_table.md`

## Run and data audit

`RUN-091` completed and was collected with `labproj`. The collected output
directory contains 29 files and is tracked by DVC with pointer
`data/processed/phot-coupled-family-discriminating-abc/RUN-091.dvc`
(`md5 bd1334082ba0a381c3bec0afc8f5a1a1.dir`, 29 files, 114,905,657 bytes).

The run used `PHOT_CHLRE_Chen_2023`, the `RUN-077` system-specific fitted PHOT
base, 24 complete graphlets, 12 calibration graphlets, and 12 held-out test
graphlets across dimensions 4, 5, and 6. The independently fitted ABC families
were no-module, empirical-module, random-module, shuffled-module, and
uncoupled-module. The comparison table also retained fitted-current,
smooth-core, fitted-current posterior mean, residual-only, GNK, and empirical
readout-residual diagnostics.

## Key observations

- The empirical-module ABC family improves the fixed held-out graphlet distance
  against the frozen current-model baselines. Its test distance is `0.2553`,
  versus `0.4016` for fitted current and `0.4119` for smooth core.
- It also beats the no-module and shuffled-module controls on this split:
  `0.2553` versus `0.3489` for no-module and `0.3535` for shuffled-module.
  The no-module family still beats fitted current and smooth core, but it is no
  longer the best coupled family after ABC.
- The specificity controls still fail. The random-module family is effectively
  tied and fractionally better (`0.255188` versus `0.255286`, delta
  `+9.72e-05` for empirical minus random), and the uncoupled-module control is
  clearly better (`0.2237`, delta `+0.0316` for empirical minus uncoupled).
- Component guardrails do not explain the failure as a hard constraint
  violation. Empirical-module, random-module, shuffled-module, uncoupled-module,
  no-module, fitted-current, and smooth-core rows all have zero guardrail
  violations on held-out graphlets. The residual-only diagnostic fails with two
  violations and distance `2.6756`; GNK has no guardrail violations but poor
  held-out distance (`0.8230`).
- The empirical-module raw held-out metrics are a partial rescue, not a full
  mechanistic match. It moves local maxima toward the empirical target
  (`4.33` versus `5.17`) and lag-1 autocorrelation (`0.156` versus `0.191`)
  better than the fitted current model, but accessibility remains far too low
  (`0.0183` versus empirical `0.1975`).
- The winning uncoupled control matches several secondary components better
  without validating empirical local coupling: accessibility is `0.1319`,
  lag-1 autocorrelation is `0.1913`, WHT high-order fraction is `0.3721`
  versus empirical `0.3788`, and GFT low/mid/high fractions are closer as a
  package, although local maxima remain too low (`3.08`).
- The empirical module energy is diffuse rather than localized. The selected
  empirical-module candidate has a plaquette effective module count of
  `9.11/12` and normalized entropy `0.937`, close to the empirical readout
  residual's diffuse pattern and not consistent with a concentrated local
  module mechanism.

## Interpretation

`RUN-091` improves the narrow held-out graphlet-distance performance of the
coupled-family model class. ABC fitting is doing something useful relative to
the previous grid search: the empirical-module family now beats fitted current,
smooth core, no-module, and shuffled-module controls on the fixed held-out
graphlet split.

That is not enough for `HYP-008`. The empirical-module result fails the two
most important specificity checks: random modules are indistinguishable under
the primary distance, and a single uncoupled auxiliary landscape is the best
held-out family. The selected empirical-module energy is also diffuse, not
locally concentrated. The result therefore supports a narrow capacity claim for
ABC-fitted coupled families, but weakens the current local non-diffuse
empirical-module mechanism interpretation.

The strongest scientific reading is mixed but negative for mechanism:
`RUN-091` shows that ABC fitting can improve the model on held-out PHOT
graphlet data, but the improvement is not specific to empirical local modules
and is better explained by flexible auxiliary-landscape capacity.

## Limitations

- This is a single fixed PHOT calibration/test graphlet split, not a replicated
  split diagnostic.
- The empirical target is the processed `latent_phi` graphlet layer from the
  current PHOT pipeline.
- The run freezes the `RUN-077` PHOT base and tests coupled-family residual
  fitting, not a fresh current-model refit.
- The comparison tests graphlet generalization and component guardrails, not
  independent assay-level validation.
- The result weakens this operationalization of `HYP-008`, not every possible
  nonlinear multi-trait local operator.
