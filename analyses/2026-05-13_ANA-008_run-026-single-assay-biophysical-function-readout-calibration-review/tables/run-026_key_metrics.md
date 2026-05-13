# RUN-026 key metrics

## Assay summary

| Quantity | Value |
| --- | --- |
| Assays | `1` |
| DMS assay | `PHOT_CHLRE_Chen_2023` |
| Assay class | `FACS` activity |
| Taxon | `Eukaryote` |
| Sequence length | `118` residues |
| Total measured variants | `167,529` |
| Single mutants | `2,122` |
| Multiple mutants | `165,407` |
| `mavenn` test Spearman | `0.680` |
| `mavenn` test NRMSE | `0.935` |
| Bayesian calibration mode | `synthetic_readout_mode: stability_function` |
| Bayesian pairwise target | `empirical_pairwise_target: function` |

Note: the Bayesian outputs in `summary.json` and `branch_validations.csv` still use the historical `smc_abc_*_raw` labels, but `EXP-008` `config.yaml` confirms that the SMC path itself used the biophysical-function calibration semantics above.

## Within-run branch comparison

| Metric | `baseline_shared_raw` | `predictive_richpair_shared_raw` | `baseline_shared_biophysical_function_readout` | `predictive_richpair_shared_biophysical_function_readout` | `smc_abc_best_raw` | `smc_abc_posterior_mean_raw` |
| --- | --- | --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.066` | `0.338` | `0.216` | `0.349` | `0.293` | `0.270` |
| Single-mutant holdout NRMSE | `0.996` | `0.957` | `0.972` | `0.942` | `0.974` | `0.954` |
| Double-mutant holdout Spearman | `0.319` | `0.342` | `0.081` | `0.551` | `0.562` | `-0.142` |
| Double-mutant holdout NRMSE | `0.989` | `1.066` | `1.067` | `0.825` | `0.857` | `1.235` |
| Epistasis-prediction Spearman | `0.468` | `0.272` | `0.051` | `0.095` | `0.109` | `0.233` |
| Epistasis-prediction KS | `0.205` | `0.290` | `0.858` | `0.778` | `0.812` | `0.619` |
| Functional KS | `0.807` | `0.296` | `0.405` | `0.357` | `0.367` | `0.561` |
| Reference fraction of peak | `0.278` | `0.032` | `1.000` | `1.000` | `0.083` | `0.112` |
| Reference distance to peak | `16.0` | `91.0` | `0.0` | `0.0` | `57.0` | `53.0` |
| Fitted `epistasis_strength` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0648` | `0.0456` |
| Fitted `empirical_pairwise_strength` | `0.0250` | `0.0500` | `0.0250` | `0.0500` | `0.0448` | `0.0373` |

## Key within-run contrasts

| Contrast | Single-mutant holdout Spearman | Double-mutant holdout Spearman | Epistasis-prediction Spearman | Functional KS | Reference fraction of peak |
| --- | --- | --- | --- | --- | --- |
| `baseline_shared_biophysical_function_readout - baseline_shared_raw` | `+0.150` | `-0.238` | `-0.416` | `-0.402` | `+0.722` |
| `predictive_richpair_shared_biophysical_function_readout - predictive_richpair_shared_raw` | `+0.010` | `+0.209` | `-0.177` | `+0.061` | `+0.968` |
| `smc_abc_best_raw - predictive_richpair_shared_raw` | `-0.046` | `+0.220` | `-0.163` | `+0.070` | `+0.051` |

Interpretation: the explicit biophysical-function readout helps some metrics, especially double-mutant holdout, but it does not improve single-mutant ranking, double-mutant ranking, epistasis prediction, and KS together. The deterministic biophysical branches also collapse the reference onto the fitted peak.

## Posterior summary for selected parameters

| Parameter | Weighted mean | q05 | q95 |
| --- | --- | --- | --- |
| `epistasis_strength` | `0.0456` | `0.0052` | `0.0769` |
| `empirical_pairwise_strength` | `0.0373` | `0.0020` | `0.0703` |
| `readout_stability_midpoint` | `1.272` | `-1.162` | `3.716` |
| `readout_stability_slope` | `2.510` | `0.819` | `3.830` |
| `readout_function_exponent` | `1.199` | `0.576` | `1.892` |
| `n_functional_dims` | `3.86` | `3.0` | `5.0` |

## Empirical SMC round summary

| Round | Epsilon | Attempts | Best distance | Median distance |
| --- | --- | --- | --- | --- |
| `0` | `14.067` | `384` | `13.237` | `13.802` |
| `1` | `15.602` | `512` | `11.919` | `13.434` |
| `2` | `15.616` | `512` | `11.938` | `13.424` |
| `3` | `15.777` | `512` | `11.765` | `13.431` |

Note: the empirical SMC round metrics in `posterior_rounds.csv` describe the Bayesian fit used in the current branch labels. They are numerically stable, but the posterior mean still does not give a robust empirical improvement.

## Synthetic-truth recovery summary

| Truth | Best distance | Parameters within posterior q90 | Misses |
| --- | --- | --- | --- |
| `moderate_epistatic` | `0.955` | `13 / 13` | `None` |
| `flatter_low_epistasis` | `4.574` | `13 / 13` | `None` |

Note: the row-level `synthetic_truth_recovery.csv` supports the `13 / 13` counts above. The aggregate `summary.json` synthetic-truth block underreports those counts and should not be treated as the authoritative source for this run.
