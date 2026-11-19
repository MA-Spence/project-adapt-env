# RES-036: RUN-091 ABC improves PHOT graphlet distance but fails random and uncoupled specificity

## Summary

`RUN-091` shows that ABC fitting of the coupled-family controls improves PHOT
held-out graphlet distance for the empirical-module family relative to fitted
current, smooth core, no-module, and shuffled-module controls. It still fails
mechanism specificity: the random-module family is effectively tied and the
uncoupled-module control is better. Selected empirical-module energy remains
diffuse. This weakens the current local non-diffuse `HYP-008` mechanism while
preserving a narrow ABC capacity result.

## Generated from

- Analyses: ANA-036
- Experiments: EXP-036
- Runs: RUN-091

## Relevant hypotheses

- Supports: None
- Weakens: HYP-008
- Refutes: None

## Evidence

- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/metrics.json
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_component_distance_summary.md
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_component_guardrails.md
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_hyp008_decision_table.md
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_module_energy_localization.md
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_pairwise_control_deltas.md
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_preflight_and_collection_audit.md
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_raw_test_metric_means.md
- analyses/2026-05-17_ANA-036_run-091-phot-coupled-family-discriminating-abc-review/tables/run-091_selected_family_distances.md
- results/RES-036_run-091-abc-improves-phot-graphlet-distance-fails-random-uncoupled-specificity/metrics.json
- data/processed/phot-coupled-family-discriminating-abc/RUN-091.dvc
- data/processed/phot-coupled-family-discriminating-abc/RUN-091/summary.json
- data/processed/phot-coupled-family-discriminating-abc/RUN-091/selected_family_distances.csv
- data/processed/phot-coupled-family-discriminating-abc/RUN-091/selected_family_pairwise_deltas.csv
- data/processed/phot-coupled-family-discriminating-abc/RUN-091/component_guardrail_summary.csv
- data/processed/phot-coupled-family-discriminating-abc/RUN-091/module_energy_summary.csv

## Interpretation

The ABC refit materially improves the primary fixed held-out graphlet metric.
The empirical-module test distance is `0.2553`, compared with `0.4016` for
fitted current, `0.4119` for smooth core, `0.3489` for no-module, and `0.3535`
for shuffled-module. Under this metric, ABC fitting is better than the earlier
grid search and does improve model performance on held-out PHOT graphlets.

The improvement is not mechanism-specific. Random-module is fractionally better
than empirical-module (`0.255188` versus `0.255286`; empirical minus random
delta `+9.72e-05`), and the uncoupled-module control is clearly better
(`0.2237`; empirical minus uncoupled delta `+0.0316`). Because the best control
does not require empirical local coupling, the empirical-module gain should not
be cited as support for local non-diffuse ruggedness.

The secondary metrics are also mixed. The empirical-module candidate improves
lag-1 autocorrelation (`0.156` versus empirical `0.191`) and local maxima
(`4.33` versus empirical `5.17`) relative to fitted-current controls, but
accessibility remains far too low (`0.0183` versus empirical `0.1975`). The
uncoupled control better matches accessibility, lag-1 autocorrelation, WHT
high-order fraction, and GFT energy fractions as a package, even though its
local maxima remain too low.

The locality diagnostic is negative. The selected empirical-module candidate
has plaquette effective module count `9.11/12` and normalized entropy `0.937`,
which is diffuse rather than concentrated. Component guardrails were not
violated, but passing guardrails is insufficient when random/uncoupled
specificity and module-locality checks fail.

## Effect on hypothesis

`RES-036` weakens `HYP-008` as currently operationalized. It supports the narrow
claim that ABC-fitted coupled-family models can lower PHOT held-out graphlet
distance relative to current and smooth baselines. It does not support the
stronger hypothesis that empirical local modules identify the local
non-diffuse ruggedness mechanism.

## Limitations

- PHOT is the only empirical system tested.
- This is one fixed calibration/test graphlet split, not a replicated split
  stability analysis.
- The target is the processed `latent_phi` graphlet layer.
- The `RUN-077` PHOT base is frozen; this does not test a fresh base refit.
- The result is graphlet-evidence only and does not establish assay-level
  generalization.

## Downstream use

- Cite `RES-036` as evidence that ABC fitting can improve fixed held-out PHOT
  graphlet distance relative to fitted-current and smooth-core controls.
- Do not cite `RUN-091` as support for empirical local-module `HYP-008`
  mechanism without noting the random-module tie, uncoupled-module win, and
  diffuse module-energy diagnostic.
