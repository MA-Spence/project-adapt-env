# RES-038: RUN-096 local ruggedness objective improves controls but not empirical-module specificity

## Summary

Adding local Dirichlet-energy localization terms to the PHOT ABC objective improves fitted-family held-out distance over fitted current/smooth core, but the calibration-selected no-module and uncoupled controls beat the empirical-module family on held-out objective distance, so the result weakens the current empirical-module mechanism.

## Generated from

- Analyses: ANA-038

## Relevant hypotheses

- Supports: None
- Weakens: HYP-008
- Refutes: None

## Evidence

- analyses/2026-05-17_ANA-038_run-096-phot-abc-local-ruggedness-objective-review/metrics.json
- analyses/2026-05-17_ANA-038_run-096-phot-abc-local-ruggedness-objective-review/tables/run-096_hyp008_decision_table.md
- analyses/2026-05-17_ANA-038_run-096-phot-abc-local-ruggedness-objective-review/tables/run-096_localization_components_test.md
- analyses/2026-05-17_ANA-038_run-096-phot-abc-local-ruggedness-objective-review/tables/run-096_preflight_and_collection_audit.md
- analyses/2026-05-17_ANA-038_run-096-phot-abc-local-ruggedness-objective-review/tables/run-096_selected_test_distances.md
- results/RES-038_run-096-local-ruggedness-objective-no-empirical-module-specificity/metrics.json
- data/processed/phot-abc-local-ruggedness-objective/RUN-096.dvc
- data/processed/phot-abc-local-ruggedness-objective/RUN-096/selected_family_distances.csv
- data/processed/phot-abc-local-ruggedness-objective/RUN-096/selected_family_pairwise_deltas.csv
- data/processed/phot-abc-local-ruggedness-objective/RUN-096/component_guardrail_summary.csv
- data/processed/phot-abc-local-ruggedness-objective/RUN-096/module_energy_summary.csv

## Interpretation

Adding local Dirichlet-energy localization metrics to the ABC objective improves
the fitted families over fitted current and smooth core under the expanded
held-out graphlet objective, but not in the mechanism-specific way required by
`HYP-008`. The calibration-selected no-module candidate is best among fitted
families on held-out full-objective distance (`0.277056`), followed by
uncoupled-module (`0.291161`). Empirical-module is worse (`0.365272`), though
it still beats fitted current (`0.521007`) and smooth core (`0.526819`).

The local ruggedness components show the same pattern. The fitted families are
closer than fitted current/smooth core on the Local7 component package, but the
selected empirical-module particle is not best: Local7 distance is `0.462290`
for empirical-module, worse than uncoupled-module (`0.289976`), no-module
(`0.326809`), random-module (`0.439467`), and shuffled-module (`0.452016`).
The empirical readout residual remains best overall (`0.135199` Local7,
`0.211455` full distance).

The module-locality evidence is also negative. The selected empirical-module
candidate remains diffuse across its modules, with held-out edge effective
module count `11.010808/12`, edge normalized entropy `0.982059`, plaquette
effective module count `9.341583/12`, and plaquette normalized entropy
`0.944129`.

A test-oracle empirical-module posterior particle reaches held-out distance
`0.259136`, better than selected no-module, but it was not selected by the
calibration objective. That is useful capacity evidence, not a valid held-out
support result.

## Effect on hypothesis

`RES-038` weakens `HYP-008` as currently operationalized. It supports only a
narrow capacity statement: ABC-fitted auxiliary/coupled landscapes can improve
PHOT held-out graphlet and localization distances relative to fitted
current/smooth core. It does not support the stronger empirical local-module
mechanism, because empirical-module is not calibration-selected best on
held-out data, loses to no-module and uncoupled controls, is not best on the
local ruggedness component package, and remains diffuse over modules.

## Limitations

- PHOT is the only empirical system in this run.
- The comparison uses one fixed calibration/test graphlet split.
- Test-oracle posterior particles are not counted as support.
- The localization statistics are dense-graphlet summaries, not full-landscape
  localization estimates.
- The empirical layer is processed `latent_phi`.

## Downstream use

- Cite `RES-038` as evidence that local-ruggedness-aware ABC improves PHOT
  fitted-family held-out distance over fitted current/smooth core.
- Do not cite `RES-038` as support for the empirical calibration module
  mechanism in `HYP-008`.
- Treat future work as needing a better mechanism or module discovery strategy,
  not just reweighting the current empirical-module objective.
