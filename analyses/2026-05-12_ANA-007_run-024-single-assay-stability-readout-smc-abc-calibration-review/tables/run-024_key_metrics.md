# RUN-024 key metrics

## Assay summary

| Quantity | Value |
| --- | --- |
| Assays | `1` |
| DMS assay | `SPTN1_CHICK_Tsuboyama_2023_1TUD` |
| Assay class | `cDNA display proteolysis` stability |
| Taxon | `Eukaryote` |
| Sequence length | `60` residues |
| Total measured variants | `3,201` |
| Single mutants | `1,051` |
| Multiple mutants | `2,150` |
| `mavenn` test Spearman | `0.908` |
| `mavenn` test NRMSE | `0.542` |
| Bayesian calibration mode | `synthetic_readout_mode: stability_margin` |
| Bayesian pairwise target | `empirical_pairwise_target: stability` |

Note: the Bayesian outputs in `summary.json` and `branch_validations.csv` still use the historical `smc_abc_*_raw` labels, but `EXP-007` `config.yaml` confirms that the SMC path itself used the stability-targeted calibration semantics above.

## Within-run branch comparison

| Metric | `baseline_shared_raw` | `predictive_richpair_shared_raw` | `baseline_shared_stability_readout` | `predictive_richpair_shared_stability_readout` | `smc_abc_best_raw` | `smc_abc_posterior_mean_raw` |
| --- | --- | --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.340` | `0.328` | `0.501` | `0.534` | `0.556` | `0.554` |
| Single-mutant holdout NRMSE | `0.946` | `0.940` | `0.876` | `0.848` | `0.876` | `0.861` |
| Double-mutant holdout Spearman | `0.504` | `0.281` | `0.572` | `0.504` | `0.483` | `0.486` |
| Double-mutant holdout NRMSE | `0.922` | `1.089` | `0.852` | `1.052` | `0.993` | `0.981` |
| Epistasis-prediction Spearman | `0.097` | `0.427` | `0.520` | `0.517` | `0.540` | `0.523` |
| Epistasis-prediction KS | `0.289` | `0.325` | `0.224` | `0.203` | `0.214` | `0.227` |
| Functional KS | `0.325` | `0.331` | `0.216` | `0.223` | `0.233` | `0.227` |
| Reference fraction of peak | `1.000` | `0.083` | `0.992` | `0.992` | `0.039` | `0.070` |
| Reference distance to peak | `0.0` | `57.0` | `3.0` | `1.0` | `50.0` | `54.0` |
| Fitted `epistasis_strength` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0686` | `0.0382` |
| Fitted `empirical_pairwise_strength` | `0.0000` | `0.0250` | `0.0500` | `0.0250` | `0.0367` | `0.0417` |

## Comparison to earlier records

| Metric | `RUN-019 smc_abc_best_raw` | `RUN-022 smc_abc_best_raw` | `RUN-012 SPTN1 per-assay raw` | `RUN-024 smc_abc_best_raw` |
| --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.195` | `0.278` | `0.406` | `0.556` |
| Double-mutant holdout Spearman | `0.119` | `0.362` | `0.164` | `0.483` |
| Epistasis-prediction Spearman | `0.265` | `0.245` | `0.529` | `0.540` |
| Functional KS | `0.483` | `0.337` | `0.414` | `0.233` |
| Reference fraction of peak | `0.049` | `0.048` | `0.039` | `0.039` |
| Reference distance to peak | `56.5` | `53.0` | `59.0` | `50.0` |
| Fitted `epistasis_strength` | `0.0527` | `0.0497` | `0.0000` | `0.0686` |
| Fitted `empirical_pairwise_strength` | `0.0008` | `0.0036` | `0.0500` | `0.0367` |

## Posterior summary for selected parameters

| Parameter | Weighted mean | q05 | q95 |
| --- | --- | --- | --- |
| `epistasis_strength` | `0.0382` | `0.0042` | `0.0734` |
| `empirical_pairwise_strength` | `0.0417` | `0.0077` | `0.0765` |
| `noise_amplitude` | `0.0052` | `0.0011` | `0.0091` |
| `peak_distance_from_consensus` | `1.89` | `0.0` | `4.0` |
| `n_functional_dims` | `3.71` | `2.0` | `6.0` |

## Empirical SMC round summary

| Round | Epsilon | Attempts | Best distance | Median distance |
| --- | --- | --- | --- | --- |
| `0` | `9.731` | `384` | `8.780` | `9.560` |
| `1` | `10.090` | `512` | `8.494` | `9.209` |
| `2` | `9.878` | `512` | `8.263` | `9.076` |
| `3` | `9.752` | `512` | `8.014` | `9.060` |

## Synthetic-truth recovery summary

| Truth | Best distance | Parameters within posterior q90 | Misses |
| --- | --- | --- | --- |
| `moderate_epistatic` | `0.418` | `10 / 10` | `None` |
| `flatter_low_epistasis` | `0.440` | `10 / 10` | `None` |
