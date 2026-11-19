# RUN-098 selected held-out test distances

Calibration-selected rows only; lower distance is better. `original9_mean` is the original graphlet/topology/spectral component package and `local7_mean` is the local Dirichlet-energy localization package.

| family | model | selected_parameter | full16_distance | original9_mean | local7_mean | guardrail_violations |
| --- | --- | --- | --- | --- | --- | --- |
| empirical_readout_residual | empirical_global_readout_residual | latent_phi__minus_poly_fitted_current_best | 0.211455 | 0.270766 | 0.135199 | 0 |
| random_module | nested_softmin_random_module | random_module__posterior_rank001 | 0.482630 | 0.427748 | 0.553193 | 0 |
| shuffled_module | nested_softmin_shuffled_module | shuffled_module__posterior_rank002 | 0.510961 | 0.581815 | 0.419862 | 0 |
| fitted_current | fitted_current_best | PHOT_CHLRE_Chen_2023__fitted_current_best | 0.521007 | 0.401571 | 0.674567 | 0 |
| smooth_core | fitted_smooth_core_only | PHOT_CHLRE_Chen_2023__fitted_smooth_core_only | 0.526819 | 0.411904 | 0.674567 | 0 |
| locked_first_stage | nested_locked_first_stage_c | RUN-084_rank001_locked_first_stage | 0.529590 | 0.465355 | 0.612177 | 0 |
| no_module | nested_softmin_no_module | no_module__posterior_rank001 | 0.537199 | 0.584359 | 0.476565 | 0 |
| empirical_module | nested_softmin_empirical_module | empirical_module__posterior_rank007 | 0.543180 | 0.557285 | 0.525045 | 0 |
| fitted_current_posterior_mean | fitted_current_posterior_mean | PHOT_CHLRE_Chen_2023__fitted_current_posterior_mean | 0.581217 | 0.440819 | 0.761728 | 0 |
| uncoupled_module | nested_softmin_uncoupled_module | uncoupled_module__posterior_rank011 | 0.640141 | 0.736476 | 0.516281 | 0 |
| gnk_local_null | gnk_local_null | PHOT_CHLRE_Chen_2023__gnk_k2 | 0.809475 | 0.765671 | 0.865795 | 1 |
| residual_only_diagnostic | fitted_residual_only_diagnostic | PHOT_CHLRE_Chen_2023__fitted_residual_only_diagnostic | 2.675584 | 2.675584 |  | 2 |
