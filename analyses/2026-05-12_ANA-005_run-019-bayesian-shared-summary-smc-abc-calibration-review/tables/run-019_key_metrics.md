# RUN-019 key metrics

## Panel summary

| Quantity | Value |
| --- | --- |
| Assays | `6` |
| Assay class | `cDNA display proteolysis` stability |
| Taxa | `Human 1`, `Prokaryote 1`, `Eukaryote 3`, `Virus 1` |
| Total measured variants | `23,279` |
| Single mutants | `6,424` |
| Multiple mutants | `16,855` |
| Sequence length range | `44` to `72` residues |
| Mean `mavenn` test Spearman | `0.890` |
| Mean `mavenn` test NRMSE | `0.480` |

## Shared branch comparison

| Metric | `baseline_shared_raw` | `predictive_richpair_shared_raw` | `smc_abc_best_raw` | `smc_abc_posterior_mean_raw` |
| --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.152` | `0.243` | `0.195` | `-0.017` |
| Single-mutant holdout NRMSE | `1.068` | `1.117` | `1.092` | `1.191` |
| Double-mutant holdout Spearman | `-0.016` | `0.098` | `0.119` | `-0.011` |
| Double-mutant holdout NRMSE | `1.120` | `1.163` | `1.226` | `1.094` |
| Epistasis-prediction Spearman | `0.472` | `0.094` | `0.265` | `-0.025` |
| Epistasis-prediction KS | `0.373` | `0.367` | `0.310` | `0.388` |
| Functional KS | `0.557` | `0.483` | `0.483` | `0.499` |
| Reference fraction of peak | `0.090` | `0.998` | `0.049` | `0.065` |
| Reference distance to peak | `57.5` | `3.0` | `56.5` | `49.0` |
| Fitted `epistasis_strength` | `0.0000` | `0.0000` | `0.0527` | `0.0429` |
| Fitted `empirical_pairwise_strength` | `0.0250` | `0.0000` | `0.0008` | `0.0050` |
| Fitted `noise_amplitude` | `0.0000` | `0.0000` | `0.0004` | `0.0052` |

## Bayesian fit versus best deterministic control

| Metric | Best deterministic control | `smc_abc_best_raw` | Delta |
| --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.243` | `0.195` | `-0.048` |
| Double-mutant holdout Spearman | `0.098` | `0.119` | `+0.021` |
| Epistasis-prediction Spearman | `0.094` | `0.265` | `+0.171` |
| Functional KS | `0.483` | `0.483` | `-0.001` |
| Reference fraction of peak | `0.998` | `0.049` | `-0.948` |
| Reference distance to peak | `3.0` | `56.5` | `+53.5` |

## Posterior summary for selected parameters

| Parameter | Weighted mean | q05 | q95 |
| --- | --- | --- | --- |
| `epistasis_strength` | `0.0429` | `0.0048` | `0.0754` |
| `empirical_pairwise_strength` | `0.0050` | `0.0005` | `0.0107` |
| `noise_amplitude` | `0.0052` | `0.0005` | `0.0093` |
| `peak_distance_from_consensus` | `2.07` | `0` | `4` |
| `n_functional_dims` | `4.00` | `3` | `6` |

## Synthetic-truth recovery summary

| Truth | Best distance | Parameters within posterior q90 | Misses |
| --- | --- | --- | --- |
| `moderate_epistatic` | `2.503` | `10 / 10` | `None` |
| `flatter_low_epistasis` | `4.206` | `9 / 10` | `noise_amplitude` |
