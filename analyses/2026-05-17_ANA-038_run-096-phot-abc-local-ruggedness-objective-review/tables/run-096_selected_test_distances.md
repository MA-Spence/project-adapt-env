# RUN-096 Calibration-Selected Distances

Lower values are closer to empirical PHOT. `Original9` is the mean distance
over the pre-RUN-096 graphlet metrics. `Local7` is the mean distance over the
added Dirichlet-energy localization terms.

## Held-Out Test

| Candidate family | Full objective | Original9 | Local7 |
|---|---:|---:|---:|
| empirical_readout_residual | 0.211455 | 0.270766 | 0.135199 |
| no_module | 0.277056 | 0.238360 | 0.326809 |
| uncoupled_module | 0.291161 | 0.292083 | 0.289976 |
| empirical_module | 0.365272 | 0.289814 | 0.462290 |
| random_module | 0.366768 | 0.310224 | 0.439467 |
| shuffled_module | 0.379458 | 0.323024 | 0.452016 |
| fitted_current | 0.521007 | 0.401571 | 0.674567 |
| smooth_core | 0.526819 | 0.411904 | 0.674567 |
| fitted_current_posterior_mean | 0.581217 | 0.440819 | 0.761728 |
| gnk_local_null | 0.832663 | 0.822951 | 0.845150 |
| residual_only_diagnostic | 2.675584 | 2.675584 | NA |

## Calibration

| Candidate family | Full objective | Original9 | Local7 |
|---|---:|---:|---:|
| uncoupled_module | 0.386328 | 0.313672 | 0.479744 |
| empirical_module | 0.388435 | 0.321996 | 0.473857 |
| no_module | 0.408368 | 0.430185 | 0.380317 |
| random_module | 0.416311 | 0.470523 | 0.346610 |
| shuffled_module | 0.445200 | 0.446410 | 0.443644 |
| empirical_readout_residual | 0.489077 | 0.542832 | 0.419964 |
| gnk_local_null | 0.851370 | 0.775948 | 0.948341 |
| fitted_current | 0.871060 | 0.731836 | 1.050061 |
| smooth_core | 0.878677 | 0.745378 | 1.050061 |
| fitted_current_posterior_mean | 0.904798 | 0.783857 | 1.060292 |
| residual_only_diagnostic | 2.709940 | 2.709940 | NA |
*** Add File: /Users/matthewspence/Documents/Documents/project-adapt-env/analyses/2026-05-17_ANA-038_run-096-phot-abc-local-ruggedness-objective-review/tables/run-096_localization_components_test.md
# RUN-096 Held-Out Localization Components

Raw held-out graphlet means for empirical PHOT and the calibration-selected
models. Concentration metrics are not all monotone in the same direction; the
distance table should be used for the normalized fit decision.

| Candidate | Edge Moran I | Edge top10 | Edge effective | Site top1 | Site effective | Plaquette top1 | Plaquette effective |
|---|---:|---:|---:|---:|---:|---:|---:|
| empirical_target | -0.003912 | 0.406000 | 0.408275 | 0.416216 | 0.748406 | 0.195004 | 0.812337 |
| empirical_readout_residual | 0.028335 | 0.415638 | 0.403175 | 0.387376 | 0.798293 | 0.195378 | 0.817757 |
| no_module | -0.007756 | 0.386679 | 0.420759 | 0.353390 | 0.814399 | 0.221319 | 0.745926 |
| uncoupled_module | 0.012438 | 0.396168 | 0.410892 | 0.343267 | 0.846321 | 0.210944 | 0.762215 |
| empirical_module | 0.047487 | 0.390309 | 0.422990 | 0.302669 | 0.885943 | 0.218381 | 0.749199 |
| random_module | 0.039033 | 0.395764 | 0.412935 | 0.300254 | 0.901063 | 0.221355 | 0.750644 |
| shuffled_module | 0.039124 | 0.414808 | 0.396689 | 0.320456 | 0.866550 | 0.230263 | 0.729763 |
| fitted_current | -0.005576 | 0.322656 | 0.455256 | 0.298503 | 0.890423 | 0.226926 | 0.730060 |
| smooth_core | -0.005576 | 0.322656 | 0.455256 | 0.298503 | 0.890423 | 0.226926 | 0.730060 |
*** Add File: /Users/matthewspence/Documents/Documents/project-adapt-env/analyses/2026-05-17_ANA-038_run-096-phot-abc-local-ruggedness-objective-review/tables/run-096_hyp008_decision_table.md
# RUN-096 HYP-008 Decision Table

| Criterion | Observation | Decision |
|---|---|---|
| Fit improves over fitted current/smooth core | No-module, uncoupled, empirical-module, random, and shuffled selected families all beat fitted current (`0.521007`) and smooth core (`0.526819`) on held-out full objective distance. | Capacity signal present |
| Local ruggedness terms improve over fitted current/smooth core | Local7 distances for selected fitted families are `0.289976` to `0.462290`, below fitted current/smooth core `0.674567`. | Directional local-metric improvement |
| Empirical-module family wins held-out selection | Empirical-module held-out full distance is `0.365272`, worse than no-module `0.277056` and uncoupled-module `0.291161`. | Fails |
| Empirical-module beats random/shuffled controls | Empirical-module narrowly beats random-module by `0.001496` and shuffled-module by `0.014186`, but this is not enough because no-module and uncoupled controls are better. | Partial, not sufficient |
| Empirical-module local metrics are best | Empirical-module Local7 distance is `0.462290`, worse than uncoupled `0.289976`, no-module `0.326809`, random `0.439467`, and shuffled `0.452016`. | Fails |
| Module energy is non-diffuse in empirical modules | Selected empirical-module test energy is diffuse: edge effective module count `11.010808/12`, edge normalized entropy `0.982059`, plaquette effective module count `9.341583/12`, plaquette normalized entropy `0.944129`. | Fails |
| Calibration-selected evidence only | A test-oracle empirical-module particle has test distance `0.259136`, but it was not calibration-selected. | Capacity only, not support |

Overall: RUN-096 supports the narrower statement that ABC-fitted coupled
families can move PHOT held-out graphlet and localization distances away from
the fitted current/smooth core. It does not support the tested empirical local
module mechanism required by HYP-008.
