# RUN-038 key metrics

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
| `mavenn` test Spearman | `0.686` |
| `mavenn` test NRMSE | `0.916` |
| Bayesian calibration mode | `synthetic_readout_mode: stability_margin` |
| Bayesian pairwise target | `empirical_pairwise_target: stability` |
| ABC distance mode | bootstrap summary-vector distance |

Note: the Bayesian outputs in `summary.json` and `branch_validations.csv` still
use the historical `smc_abc_*_raw` labels, but `EXP-011` `config.yaml`
confirms that the SMC path itself used the stability-targeted calibration
semantics above.

## Summary target features

| Feature | Observed value | Bootstrap SD |
| --- | --- | --- |
| `PHOT_CHLRE_Chen_2023__conservation_correlation` | `0.547` | `0.059` |
| `PHOT_CHLRE_Chen_2023__epistasis_variance` | `0.036` | `0.004` |
| `PHOT_CHLRE_Chen_2023__fraction_beneficial` | `0.008` | `0.002` |
| `PHOT_CHLRE_Chen_2023__fraction_deleterious` | `0.331` | `0.010` |
| `PHOT_CHLRE_Chen_2023__fraction_lethal` | `0.039` | `0.004` |
| `PHOT_CHLRE_Chen_2023__fraction_neutral` | `0.661` | `0.010` |
| `PHOT_CHLRE_Chen_2023__mean` | `-0.250` | `0.008` |
| `PHOT_CHLRE_Chen_2023__skewness` | `-1.376` | `0.045` |
| `PHOT_CHLRE_Chen_2023__variance` | `0.124` | `0.005` |

## Within-run branch comparison

| Metric | `baseline_shared_raw` | `predictive_richpair_shared_raw` | `baseline_shared_stability_readout` | `predictive_richpair_shared_stability_readout` | `smc_abc_best_raw` | `smc_abc_posterior_mean_raw` |
| --- | --- | --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.066` | `0.338` | `0.504` | `0.505` | `0.505` | `0.504` |
| Single-mutant holdout NRMSE | `0.996` | `0.957` | `0.851` | `0.871` | `0.866` | `0.874` |
| Double-mutant holdout Spearman | `0.315` | `0.338` | `0.538` | `0.701` | `0.721` | `0.694` |
| Double-mutant holdout NRMSE | `1.002` | `1.062` | `0.956` | `0.773` | `0.803` | `0.752` |
| Epistasis-prediction Spearman | `0.435` | `0.259` | `-0.250` | `0.008` | `0.066` | `-0.156` |
| Epistasis-prediction KS | `0.267` | `0.295` | `0.750` | `0.739` | `0.676` | `0.767` |
| Functional KS | `0.785` | `0.291` | `0.176` | `0.178` | `0.166` | `0.182` |
| Reference fraction of peak | `0.279` | `0.032` | `1.000` | `0.987` | `0.867` | `0.337` |
| Reference distance to peak | `16.0` | `91.0` | `0.0` | `3.0` | `4.0` | `21.0` |
| Fitted `epistasis_strength` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0451` | `0.0343` |
| Fitted `empirical_pairwise_strength` | `0.0250` | `0.0500` | `0.0250` | `0.0500` | `0.0790` | `0.0439` |

## Key within-run contrast

| Contrast | Single-mutant holdout Spearman | Double-mutant holdout Spearman | Epistasis-prediction Spearman | Functional KS | Reference fraction of peak |
| --- | --- | --- | --- | --- | --- |
| `smc_abc_best_raw - predictive_richpair_shared_stability_readout` | `-0.000` | `+0.020` | `+0.058` | `-0.012` | `-0.120` |

Interpretation: the Bayesian control fit only modestly improves the strongest
deterministic stability-readout branch. The main gain is double-mutant holdout,
while the peak-geometry problem remains.

## Cross-run comparison to prior PHOT Bayesian fits

| Metric | `RUN-026 smc_abc_best_raw` | `RUN-030 smc_abc_best_raw` | `RUN-036 smc_abc_best_raw` | `RUN-038 smc_abc_best_raw` |
| --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.293` | `0.294` | `0.198` | `0.505` |
| Double-mutant holdout Spearman | `0.562` | `0.363` | `0.251` | `0.721` |
| Epistasis-prediction Spearman | `0.109` | `-0.026` | `0.259` | `0.066` |
| Functional KS | `0.367` | `0.253` | `0.496` | `0.166` |
| Reference fraction of peak | `0.083` | `0.999` | `0.986` | `0.867` |
| Reference distance to peak | `57.0` | `5.0` | `6.0` | `4.0` |

Interpretation: `RUN-038` is the strongest PHOT fit so far on ranking and KS,
but not on epistasis prediction or peak geometry.

## Comparison to the earlier successful stability-assay fit

| Metric | `RUN-024 smc_abc_best_raw` | `RUN-038 smc_abc_best_raw` |
| --- | --- | --- |
| Single-mutant holdout Spearman | `0.556` | `0.505` |
| Double-mutant holdout Spearman | `0.483` | `0.721` |
| Epistasis-prediction Spearman | `0.540` | `0.066` |
| Functional KS | `0.233` | `0.166` |
| Reference fraction of peak | `0.039` | `0.867` |
| Reference distance to peak | `50.0` | `4.0` |

Interpretation: the stability-targeted formulation does not perform as cleanly
on the activity readout as it did on the unfolding assay. It improves double-
mutant ranking and KS on `PHOT_CHLRE`, but the activity fit remains much weaker
on epistasis prediction and peak geometry.

## Posterior summary for selected parameters

| Parameter | Weighted mean | q05 | q95 |
| --- | --- | --- | --- |
| `stability_scale` | `0.949` | `0.541` | `1.466` |
| `stability_margin` | `8.350` | `5.843` | `11.793` |
| `blosum_blend` | `0.558` | `0.239` | `0.785` |
| `stability_conservation_power` | `0.956` | `0.762` | `1.212` |
| `functional_sigma_base` | `19.278` | `10.753` | `28.649` |
| `n_functional_dims` | `4.039` | `3.0` | `6.0` |
| `peak_distance_from_consensus` | `1.894` | `0.0` | `4.0` |
| `epistasis_strength` | `0.0343` | `0.0023` | `0.0743` |
| `empirical_pairwise_strength` | `0.0439` | `0.0037` | `0.0786` |
| `noise_amplitude` | `0.0050` | `0.0005` | `0.0098` |

## Empirical SMC round summary

| Round | Epsilon | Attempts | Best distance | Median distance |
| --- | --- | --- | --- | --- |
| `0` | `11.903` | `384` | `10.648` | `11.556` |
| `1` | `12.990` | `512` | `10.909` | `11.428` |
| `2` | `12.951` | `512` | `10.699` | `11.405` |
| `3` | `12.973` | `512` | `10.959` | `11.430` |

Interpretation: the control converges to a usable posterior and recovers the
empirical target substantially better than the recent PHOT activity-oriented
fits on ranking and KS, but the round summaries do not show a strong monotonic
collapse onto one clean solution.

## Synthetic-truth recovery summary

| Truth | Best distance | Parameters within posterior q90 | Misses |
| --- | --- | --- | --- |
| `moderate_epistatic` | `0.165` | `10 / 10` | `None` |
| `flatter_low_epistasis` | `0.261` | `10 / 10` | `None` |
