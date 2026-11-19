# RUN-123 Comparator Context

| candidate family | candidate model | parameter set | test distance |
| --- | --- | --- | --- |
| empirical_readout_residual | empirical_global_readout_residual | latent_phi__minus_poly_fitted_current_best | 0.2115 |
| random_module | nested_softmin_random_module | random_module__posterior_rank001 | 0.4826 |
| shuffled_module | nested_softmin_shuffled_module | shuffled_module__posterior_rank002 | 0.5110 |
| fitted_current | fitted_current_best | PHOT_CHLRE_Chen_2023__fitted_current_best | 0.5210 |
| smooth_core | fitted_smooth_core_only | PHOT_CHLRE_Chen_2023__fitted_smooth_core_only | 0.5268 |
| locked_first_stage | nested_locked_first_stage_c | RUN-084_rank001_locked_first_stage | 0.5296 |
| no_module | nested_softmin_no_module | no_module__posterior_rank001 | 0.5372 |
| empirical_module | nested_softmin_empirical_module | empirical_module__posterior_rank007 | 0.5432 |
| fitted_current_posterior_mean | fitted_current_posterior_mean | PHOT_CHLRE_Chen_2023__fitted_current_posterior_mean | 0.5812 |
| uncoupled_module | nested_softmin_uncoupled_module | uncoupled_module__posterior_rank011 | 0.6401 |

RUN-123 generated no sparse-Walsh FIT candidate rows, so these values are context only; no operational equivalence comparison was possible.
