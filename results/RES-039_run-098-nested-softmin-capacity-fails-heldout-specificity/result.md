# RES-039: RUN-098 nested soft-min capacity fails calibration-selected heldout specificity

## Summary

Nested third-landscape soft-min capacity does not improve calibration-selected held-out PHOT full-objective distance over fitted current, smooth core, or the locked first-stage model; random/shuffled controls also beat empirical-module, so the run weakens the current empirical-module implementation of HYP-008.

## Generated from

- Analyses: ANA-039

## Relevant hypotheses

- Supports: None
- Weakens: HYP-008
- Refutes: None

## Evidence

- analyses/2026-05-18_ANA-039_run-098-phot-nested-softmin-ensemble-capacity-review/metrics.json
- analyses/2026-05-18_ANA-039_run-098-phot-nested-softmin-ensemble-capacity-review/tables/run-098_preflight_and_collection_audit.md
- analyses/2026-05-18_ANA-039_run-098-phot-nested-softmin-ensemble-capacity-review/tables/run-098_selected_test_distances.md
- analyses/2026-05-18_ANA-039_run-098-phot-nested-softmin-ensemble-capacity-review/tables/run-098_localization_components_test.md
- analyses/2026-05-18_ANA-039_run-098-phot-nested-softmin-ensemble-capacity-review/tables/run-098_module_locality_test.md
- analyses/2026-05-18_ANA-039_run-098-phot-nested-softmin-ensemble-capacity-review/tables/run-098_test_oracle_capacity_diagnostic.md
- analyses/2026-05-18_ANA-039_run-098-phot-nested-softmin-ensemble-capacity-review/tables/run-098_hyp008_decision_table.md
- results/RES-039_run-098-nested-softmin-capacity-fails-heldout-specificity/metrics.json
- data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc
- data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098/selected_family_distances.csv
- data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098/selected_family_pairwise_deltas.csv
- data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098/component_guardrail_summary.csv
- data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098/module_energy_summary.csv
- data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098/model_distances.csv

## Interpretation

The nested third-landscape soft-min ensemble does not pass the
calibration-selected held-out test for `HYP-008`. The empirical-module nested
candidate has held-out full-objective distance `0.543180`, worse than the
locked first-stage model (`0.529590`), fitted current (`0.521007`), smooth core
(`0.526819`), no-module (`0.537199`), random-module (`0.482630`), and
shuffled-module (`0.510961`). It only beats uncoupled-module (`0.640141`).

The best calibration-selected nested family is the random-module negative
control (`0.482630`). That indicates generic latent/coupled capacity can still
move the held-out metric, but it does not support a residue-local empirical
module mechanism. The local-ruggedness package is also not empirical-module
specific: Local7 distance is `0.419862` for shuffled-module, `0.476565` for
no-module, `0.525045` for empirical-module, and `0.674567` for fitted
current/smooth core.

The module-locality evidence is negative. The selected empirical-module model
is diffuse across modules on held-out graphlets, with edge effective module
count `11.727783/12`, edge normalized entropy `0.995298`, plaquette effective
module count `10.560660/12`, and plaquette normalized entropy `0.972360`.

The empirical readout residual remains the best overall diagnostic
(`0.211455` full distance; `0.135199` Local7), showing that residual PHOT
structure exists. Test-oracle posterior particles also show capacity
(`0.395321` no-module and `0.409176` empirical-module), but those particles
were not selected by calibration and do not count as held-out support.

## Effect on hypothesis

`RES-039` weakens `HYP-008` as currently operationalized. It specifically
weakens the idea that adding another nested coupled latent landscape will turn
the current empirical-module construction into a held-out, locally
non-diffuse, mechanism-specific explanation of PHOT ruggedness. The remaining
support is only a narrow capacity statement: the posterior contains candidates
and negative-control families that can move the metric, but calibration does
not select an empirical-module mechanism that beats the required controls.

## Limitations

- PHOT is the only empirical system in this run.
- The comparison uses one fixed calibration/test graphlet split.
- Test-oracle posterior particles are reported only as capacity diagnostics.
- The local ruggedness metrics are dense-graphlet summaries, not full
  landscape-wide localization estimates.
- The empirical layer is processed `latent_phi`.

## Downstream use

- Cite `RES-039` as evidence against pursuing nested ensemble capacity alone as
  the next mechanism-specific fix for `HYP-008`.
- Do not cite `RES-039` as support for empirical module localization.
- Treat the test-oracle rows as a calibration/selection failure signal, not as
  held-out validation.
