# RUN-022 key metrics

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
| `mavenn` test Spearman | `0.907` |
| `mavenn` test NRMSE | `0.483` |

## Within-run branch comparison

| Metric | `baseline_shared_raw` | `predictive_richpair_shared_raw` | `smc_abc_best_raw` | `smc_abc_posterior_mean_raw` |
| --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.340` | `0.310` | `0.278` | `0.286` |
| Single-mutant holdout NRMSE | `0.943` | `0.930` | `0.976` | `0.940` |
| Double-mutant holdout Spearman | `0.504` | `0.283` | `0.362` | `0.287` |
| Double-mutant holdout NRMSE | `0.921` | `1.153` | `1.038` | `1.052` |
| Epistasis-prediction Spearman | `0.095` | `0.380` | `0.245` | `0.546` |
| Epistasis-prediction KS | `0.360` | `0.242` | `0.268` | `0.305` |
| Functional KS | `0.355` | `0.301` | `0.337` | `0.347` |
| Reference fraction of peak | `1.000` | `0.042` | `0.048` | `0.048` |
| Reference distance to peak | `0.0` | `57.0` | `53.0` | `56.0` |
| Fitted `epistasis_strength` | `0.0000` | `0.0000` | `0.0497` | `0.0344` |
| Fitted `empirical_pairwise_strength` | `0.0000` | `0.0500` | `0.0036` | `0.0400` |

## Comparison to earlier records

| Metric | `RUN-019 smc_abc_best_raw` | `RUN-012 SPTN1 per-assay raw` | `RUN-022 smc_abc_best_raw` |
| --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.195` | `0.406` | `0.278` |
| Double-mutant holdout Spearman | `0.119` | `0.164` | `0.362` |
| Epistasis-prediction Spearman | `0.265` | `0.529` | `0.245` |
| Functional KS | `0.483` | `0.414` | `0.337` |
| Reference fraction of peak | `0.049` | `0.039` | `0.048` |
| Fitted `epistasis_strength` | `0.0527` | `0.0000` | `0.0497` |
| Fitted `empirical_pairwise_strength` | `0.0008` | `0.0500` | `0.0036` |

## Posterior summary for selected parameters

| Parameter | Weighted mean | q05 | q95 |
| --- | --- | --- | --- |
| `epistasis_strength` | `0.0344` | `0.0025` | `0.0696` |
| `empirical_pairwise_strength` | `0.0400` | `0.0023` | `0.0753` |
| `noise_amplitude` | `0.0031` | `0.0001` | `0.0078` |
| `peak_distance_from_consensus` | `2.11` | `0` | `4` |
| `n_functional_dims` | `3.86` | `2` | `6` |

## Synthetic-truth recovery summary

| Truth | Best distance | Parameters within posterior q90 | Misses |
| --- | --- | --- | --- |
| `moderate_epistatic` | `0.230` | `10 / 10` | `None` |
| `flatter_low_epistasis` | `0.509` | `8 / 10` | `empirical_pairwise_strength`, `noise_amplitude` |
