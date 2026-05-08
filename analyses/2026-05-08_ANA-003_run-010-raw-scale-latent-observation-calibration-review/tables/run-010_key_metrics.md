# RUN-010 key metrics

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
| Mean `mavenn` test Spearman | `0.894` |
| Mean `mavenn` test NRMSE | `0.491` |

## Branch comparison

| Metric | Shared raw branch | Shared latent branch | Per-assay latent mean |
| --- | --- | --- | --- |
| Single-mutant holdout Spearman | `0.232` | `0.084` | `0.136` |
| Single-mutant holdout NRMSE | `1.075` | `1.127` | `1.219` |
| Double-mutant holdout Spearman | `-0.017` | `0.010` | `0.116` |
| Double-mutant holdout NRMSE | `1.114` | `1.106` | `1.100` |
| Epistasis-prediction Spearman | `0.500` | `-0.016` | `0.236` |
| Epistasis-prediction KS | `0.338` | `0.648` | `0.606` |
| Functional KS | `0.571` | `0.479` | `0.546` |
| Reference fraction of peak | `0.104` | `0.897` | `0.181` |
| Fitted `epistasis_strength` | `0.000` | `0.000` | `0.000` in all `6` fits |
| Fitted `empirical_pairwise_strength` | `0.025` | `0.000` | `0.000` to `0.050` |

## Per-assay latent fit summary

| Quantity | Value |
| --- | --- |
| Assays with holdout Spearman `>= 0.2` | `2 / 5` non-null |
| Assays with double-holdout Spearman `>= 0.2` | `1 / 5` non-null |
| Assays with epistasis-prediction Spearman `>= 0.2` | `4 / 5` non-null |
| Best single-mutant holdout Spearman | `0.333` for `SPTN1_CHICK_Tsuboyama_2023_1TUD` |
| Best double-mutant holdout Spearman | `0.442` for `SPTN1_CHICK_Tsuboyama_2023_1TUD` |
