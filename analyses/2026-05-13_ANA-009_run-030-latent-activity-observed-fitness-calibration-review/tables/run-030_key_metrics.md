# RUN-030 key metrics

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
| `mavenn` test Spearman | `0.689` |
| `mavenn` test NRMSE | `0.921` |
| Bayesian calibration mode | `synthetic_readout_mode: fitness` |
| Observed-fitness combine mode | `observed_fitness_combine_mode: product` |
| Observed-fitness terms | `[stability_gate, function_capacity]` |
| Generic epistasis target | `generic_epistasis_target: function` |
| Bayesian pairwise target | `empirical_pairwise_target: function` |

## Within-run branch comparison

| Metric | `baseline_shared_activity_readout` | `predictive_richpair_shared_activity_readout` | `smc_abc_best_raw` | `smc_abc_posterior_mean_raw` |
| --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.216` | `0.349` | `0.294` | `0.094` |
| Single-mutant holdout NRMSE | `0.972` | `0.942` | `0.975` | `1.005` |
| Double-mutant holdout Spearman | `0.081` | `0.551` | `0.363` | `-0.052` |
| Double-mutant holdout NRMSE | `1.067` | `0.825` | `0.971` | `1.142` |
| Epistasis-prediction Spearman | `0.051` | `0.095` | `-0.026` | `-0.200` |
| Epistasis-prediction KS | `0.858` | `0.778` | `0.682` | `0.722` |
| Functional KS | `0.405` | `0.356` | `0.253` | `0.612` |
| Reference fraction of peak | `1.000` | `1.000` | `0.999` | `0.959` |
| Reference distance to peak | `0.0` | `0.0` | `5.0` | `4.0` |
| Fitted `epistasis_strength` | `0.0000` | `0.0000` | `0.0424` | `0.0499` |
| Fitted `empirical_pairwise_strength` | `0.0250` | `0.0500` | `0.0296` | `0.0306` |

Interpretation: the Bayesian fit lowers functional KS, but it does not improve
the epistasis-relevant predictive package and it still leaves the fitted
reference very near the peak.

## Key within-run contrast

| Contrast | Single-mutant holdout Spearman | Double-mutant holdout Spearman | Epistasis-prediction Spearman | Functional KS | Reference fraction of peak |
| --- | --- | --- | --- | --- | --- |
| `smc_abc_best_raw - predictive_richpair_shared_activity_readout` | `-0.054` | `-0.188` | `-0.121` | `-0.104` | `-0.001` |

Interpretation: relative to the strongest deterministic control, the best
Bayesian fit wins only on KS and barely moves the reference away from the peak.

## Cross-run comparison to RUN-026

### Deterministic branch equivalence

| Metric | `RUN-026 predictive_richpair_shared_biophysical_function_readout` | `RUN-030 predictive_richpair_shared_activity_readout` | Delta |
| --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.349` | `0.349` | `+0.000` |
| Double-mutant holdout Spearman | `0.551` | `0.551` | `+0.000` |
| Epistasis-prediction Spearman | `0.095` | `0.095` | `+0.000` |
| Functional KS | `0.357` | `0.356` | `-0.001` |
| Reference fraction of peak | `1.000` | `1.000` | `+0.000` |
| Reference distance to peak | `0.0` | `0.0` | `+0.0` |

Interpretation: the generalized observed-fitness composition path largely
reproduces the earlier deterministic activity-readout behavior on this assay.

### Bayesian best-particle comparison

| Metric | `RUN-026 smc_abc_best_raw` | `RUN-030 smc_abc_best_raw` | Delta |
| --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.293` | `0.294` | `+0.002` |
| Double-mutant holdout Spearman | `0.562` | `0.363` | `-0.199` |
| Epistasis-prediction Spearman | `0.109` | `-0.026` | `-0.135` |
| Functional KS | `0.367` | `0.253` | `-0.114` |
| Reference fraction of peak | `0.083` | `0.999` | `+0.917` |
| Reference distance to peak | `57.0` | `5.0` | `-52.0` |

Interpretation: promoting activity to the public observed-fitness generator
lowers KS but materially worsens double-mutant recovery, epistasis prediction,
and peak geometry relative to the earlier Bayesian activity-readout fit.

## Posterior summary for selected parameters

| Parameter | Weighted mean | q05 | q95 |
| --- | --- | --- | --- |
| `epistasis_strength` | `0.0499` | `0.0151` | `0.0752` |
| `empirical_pairwise_strength` | `0.0306` | `0.0053` | `0.0773` |
| `functional_sigma_base` | `16.676` | `10.362` | `24.419` |
| `n_functional_dims` | `4.02` | `2.0` | `6.0` |
| `stability_margin` | `8.542` | `5.670` | `11.765` |
| `peak_distance_from_consensus` | `1.74` | `1.0` | `4.0` |

## Empirical SMC round summary

| Round | Epsilon | Attempts | Best distance | Median distance |
| --- | --- | --- | --- | --- |
| `0` | `22.838` | `384` | `12.217` | `18.789` |
| `1` | `26.588` | `512` | `12.365` | `14.812` |
| `2` | `26.215` | `512` | `12.106` | `13.382` |
| `3` | `25.621` | `512` | `9.989` | `13.144` |

## Synthetic-truth recovery summary

| Truth | Best distance | Parameters within posterior q90 | Misses |
| --- | --- | --- | --- |
| `moderate_epistatic` | `1.025` | `10 / 10` | `None` |
| `flatter_low_epistasis` | `9.979` | `10 / 10` | `None` |

Interpretation: the inverse problem remains identifiable on matched synthetic
data, so the negative empirical result is not explained by a trivial failure of
the SMC machinery.
