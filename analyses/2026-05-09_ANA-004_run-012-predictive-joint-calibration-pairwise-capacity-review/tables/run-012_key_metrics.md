# RUN-012 key metrics

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
| Mean `mavenn` test NRMSE | `0.455` |

## Branch comparison

| Metric | `baseline_shared_raw` | `predictive_shared_raw` | `predictive_richpair_shared_raw` | `predictive_richpair_shared_latent` | Per-assay raw mean |
| --- | --- | --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.174` | `0.202` | `0.243` | `0.053` | `0.201` |
| Single-mutant holdout NRMSE | `1.070` | `1.042` | `1.117` | `1.158` | `1.157` |
| Double-mutant holdout Spearman | `-0.035` | `-0.197` | `0.098` | `0.029` | `0.089` |
| Double-mutant holdout NRMSE | `1.118` | `1.346` | `1.163` | `1.203` | `1.391` |
| Epistasis-prediction Spearman | `0.559` | `0.505` | `0.094` | `-0.018` | `0.294` |
| Epistasis-prediction KS | `0.435` | `0.330` | `0.334` | `0.694` | `0.365` |
| Functional KS | `0.539` | `0.539` | `0.463` | `0.488` | `0.570` |
| Reference fraction of peak | `0.091` | `0.095` | `0.998` | `0.900` | `0.196` |
| Fitted `epistasis_strength` | `0.000` | `0.000` | `0.000` | `0.000` | `0.000` in all `6` fits |
| Fitted `empirical_pairwise_strength` | `0.025` | `0.025` | `0.000` | `0.000` | mean `0.038`, range `0.000` to `0.050` |

## Paired comparison to `RUN-010` shared raw baseline

| Metric | `RUN-010 shared_raw` | Best `RUN-012` shared branch | Delta |
| --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.232` | `0.243` | `+0.011` |
| Double-mutant holdout Spearman | `-0.017` | `0.098` | `+0.115` |
| Functional KS | `0.571` | `0.463` | `-0.109` |
| Epistasis-prediction Spearman | `0.500` | `0.094` | `-0.407` |
| Reference fraction of peak | `0.104` | `0.998` | `+0.894` |

## Per-assay raw fit summary

| Quantity | Value |
| --- | --- |
| Assays with non-null single-mutant holdout Spearman | `5 / 6` |
| Assays with single-mutant holdout Spearman `>= 0.2` | `2 / 5` non-null |
| Assays with non-null double-mutant holdout Spearman | `5 / 6` |
| Assays with double-mutant holdout Spearman `>= 0.2` | `1 / 5` non-null |
| Assays with nonzero `empirical_pairwise_strength` | `5 / 6` |
| Assays with nonzero `epistasis_strength` | `0 / 6` |
| Best single-mutant holdout Spearman | `0.406` for `SPTN1_CHICK_Tsuboyama_2023_1TUD` |
| Best double-mutant holdout Spearman | `0.266` for `HECD1_HUMAN_Tsuboyama_2023_3DKM` |
| Best epistasis-prediction Spearman | `0.571` for `POLG_PESV_Tsuboyama_2023_2MXD` |
