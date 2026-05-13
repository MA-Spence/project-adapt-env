# RUN-036 key metrics

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
| `mavenn` test NRMSE | `0.865` |
| Explicit latent traits | built-in `stability` plus named `readout` block |
| Observed-fitness combine mode | `product` |
| Observed-fitness terms | `[stability_gate, trait:readout:capacity]` |
| Generic epistasis target | `score`, with `epistasis_strength = 0` |
| Bayesian pairwise target | `trait:readout` |
| ABC distance mode | `validation_objective` |
| Objective target features | `[objective_core, objective_double, objective_total]` |
| Synthetic-truth recovery | not run in this mode |

## Validation-objective target

| Feature | Observed value | Bootstrap SD |
| --- | --- | --- |
| `objective_core` | `0.0` | `0.0` |
| `objective_double` | `0.0` | `0.0` |
| `objective_total` | `0.0` | `0.0` |

Interpretation: the SMC fit directly minimizes the held-out validation objective
to zero rather than matching an auxiliary summary-statistics vector.

## Within-run branch comparison

| Metric | `baseline_shared_two_trait_readout` | `predictive_richpair_shared_two_trait_readout` | `smc_abc_best_raw` | `smc_abc_posterior_mean_raw` |
| --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.191` | `0.161` | `0.198` | `0.198` |
| Single-mutant holdout NRMSE | `0.983` | `1.000` | `0.981` | `0.982` |
| Double-mutant holdout Spearman | `0.216` | `0.192` | `0.251` | `0.251` |
| Double-mutant holdout NRMSE | `1.073` | `1.099` | `1.119` | `1.118` |
| Epistasis-prediction Spearman | `0.240` | `0.143` | `0.259` | `0.254` |
| Epistasis-prediction KS | `0.750` | `0.795` | `0.483` | `0.500` |
| Functional KS | `0.305` | `0.344` | `0.496` | `0.489` |
| Reference fraction of peak | `1.000` | `1.000` | `0.986` | `0.986` |
| Reference distance to peak | `2.0` | `10.0` | `6.0` | `5.0` |
| Fitted `functional_sigma_base` | `10.0` | `28.0` | `17.591` | `20.996` |
| Fitted `empirical_pairwise_strength` | `0.0000` | `0.0500` | `0.0782` | `0.0768` |
| Fitted `stability_margin` | `8.0` | `8.0` | `8.033` | `9.330` |

Interpretation: the Bayesian fit improves the within-run ranking and epistasis
metrics modestly, but it pays for those gains with a substantially worse
functional KS and it still leaves the fitted reference very near the peak.

## Key within-run contrast

| Contrast | Single-mutant holdout Spearman | Double-mutant holdout Spearman | Epistasis-prediction Spearman | Functional KS | Reference fraction of peak |
| --- | --- | --- | --- | --- | --- |
| `smc_abc_best_raw - baseline_shared_two_trait_readout` | `+0.007` | `+0.035` | `+0.019` | `+0.191` | `-0.014` |

Interpretation: relative to the stronger deterministic predictive branch, the
Bayesian fit gains only small ranking improvements while worsening the main
distributional fit metric.

## Cross-run comparison to prior PHOT Bayesian fits

| Metric | `RUN-026 smc_abc_best_raw` | `RUN-030 smc_abc_best_raw` | `RUN-036 smc_abc_best_raw` |
| --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.293` | `0.294` | `0.198` |
| Double-mutant holdout Spearman | `0.562` | `0.363` | `0.251` |
| Epistasis-prediction Spearman | `0.109` | `-0.026` | `0.259` |
| Functional KS | `0.367` | `0.253` | `0.496` |
| Reference fraction of peak | `0.083` | `0.999` | `0.986` |
| Reference distance to peak | `57.0` | `5.0` | `6.0` |

Interpretation: the cleaner two-trait validation-objective fit improves
epistasis-prediction Spearman relative to both earlier PHOT Bayesian runs, but
it is worse on single-mutant ranking, double-mutant ranking, KS, and it does
not recover the non-pathological reference geometry seen in `RUN-026`.

## Delta versus RUN-030 Bayesian best

| Metric | Delta (`RUN-036 - RUN-030`) |
| --- | --- |
| Single-mutant holdout Spearman | `-0.096` |
| Double-mutant holdout Spearman | `-0.111` |
| Epistasis-prediction Spearman | `+0.284` |
| Functional KS | `+0.243` |
| Reference fraction of peak | `-0.014` |
| Reference distance to peak | `+1.0` |

## Delta versus RUN-026 Bayesian best

| Metric | Delta (`RUN-036 - RUN-026`) |
| --- | --- |
| Single-mutant holdout Spearman | `-0.095` |
| Double-mutant holdout Spearman | `-0.311` |
| Epistasis-prediction Spearman | `+0.150` |
| Functional KS | `+0.130` |
| Reference fraction of peak | `+0.903` |
| Reference distance to peak | `-51.0` |

## Posterior summary for selected parameters

| Parameter | Weighted mean | q05 | q95 |
| --- | --- | --- | --- |
| `stability_scale` | `0.959` | `0.590` | `1.371` |
| `stability_margin` | `9.330` | `5.887` | `11.741` |
| `blosum_blend` | `0.486` | `0.244` | `0.773` |
| `stability_conservation_power` | `0.999` | `0.784` | `1.215` |
| `functional_sigma_base` | `20.996` | `13.087` | `29.283` |
| `empirical_pairwise_strength` | `0.0768` | `0.0725` | `0.0796` |

Interpretation: the posterior drives the only fitted epistatic lever,
`empirical_pairwise_strength`, to the top of its prior range without producing
a strong empirical fit.

## Empirical SMC round summary

| Round | Epsilon | Attempts | Best distance | Median distance |
| --- | --- | --- | --- | --- |
| `0` | `24.383` | `384` | `23.346` | `23.768` |
| `1` | `24.785` | `512` | `23.233` | `23.504` |
| `2` | `23.778` | `512` | `23.257` | `23.362` |
| `3` | `23.742` | `512` | `23.248` | `23.334` |

Interpretation: the posterior improves only weakly across rounds, which is
consistent with a shallow fit landscape under the cleaned-up two-trait model.
