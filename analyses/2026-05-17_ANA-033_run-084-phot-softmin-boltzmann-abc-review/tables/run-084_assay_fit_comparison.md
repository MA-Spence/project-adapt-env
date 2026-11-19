# RUN-084 assay fit comparison

| Candidate | Parameter set | Objective | Holdout Spearman | Single Spearman | Double holdout Spearman | Epistasis-prediction Spearman |
|---|---|---:|---:|---:|---:|---:|
| soft-min/Boltzmann posterior mean | `abc_posterior_mean` | 0.4570 | 0.7527 | 0.4571 | 0.5359 | -0.0881 |
| fitted-current posterior mean | `PHOT_CHLRE_Chen_2023__posterior_mean` | 0.4631 | 0.7532 | 0.4820 | 0.5980 | -0.0511 |
| soft-min/Boltzmann test-oracle particle | `abc_particle_rank_023` | 0.4752 | 0.7428 | 0.4735 | 0.5683 | -0.1526 |
| smooth core | `PHOT_CHLRE_Chen_2023__smooth_core_A` | 0.4814 | 0.7573 | 0.4782 | 0.6395 | -0.1745 |
| fitted current | `PHOT_CHLRE_Chen_2023__current_best` | 0.4856 | 0.7574 | 0.4782 | 0.6395 | -0.2042 |
| soft-min/Boltzmann, calibration selected | `abc_best_particle` | 0.4860 | 0.7555 | 0.4378 | 0.5331 | -0.1835 |

The posterior mean has the best aggregate assay objective, but the improvement
is small and does not improve the main validation correlations as a package.
Assay improvement alone was preregistered as a guardrail, not support for
`HYP-008`.
