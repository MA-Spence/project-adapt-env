# RUN-098 test-oracle capacity diagnostic

Best held-out row per family across all posterior particles. These rows are not used for the primary scientific decision because held-out graphlets were not available to the calibration selector.

| family | model | best_test_parameter | posterior_rank | test_distance | original9_mean | local7_mean |
| --- | --- | --- | --- | --- | --- | --- |
| empirical_readout_residual | empirical_global_readout_residual | latent_phi__minus_poly_fitted_current_best |  | 0.211455 | 0.270766 | 0.135199 |
| no_module | nested_softmin_no_module | no_module__posterior_rank014 | 14 | 0.395321 | 0.399523 | 0.389919 |
| empirical_module | nested_softmin_empirical_module | empirical_module__posterior_rank016 | 16 | 0.409176 | 0.435009 | 0.375961 |
| random_module | nested_softmin_random_module | random_module__posterior_rank012 | 12 | 0.418306 | 0.441020 | 0.389102 |
| shuffled_module | nested_softmin_shuffled_module | shuffled_module__posterior_rank003 | 3 | 0.456758 | 0.474790 | 0.433573 |
| fitted_current | fitted_current_best | PHOT_CHLRE_Chen_2023__fitted_current_best |  | 0.521007 | 0.401571 | 0.674567 |
| smooth_core | fitted_smooth_core_only | PHOT_CHLRE_Chen_2023__fitted_smooth_core_only |  | 0.526819 | 0.411904 | 0.674567 |
| locked_first_stage | nested_locked_first_stage_c | RUN-084_rank001_locked_first_stage |  | 0.529590 | 0.465355 | 0.612177 |
| uncoupled_module | nested_softmin_uncoupled_module | uncoupled_module__posterior_rank002 | 2 | 0.568014 | 0.526309 | 0.621636 |
| fitted_current_posterior_mean | fitted_current_posterior_mean | PHOT_CHLRE_Chen_2023__fitted_current_posterior_mean |  | 0.581217 | 0.440819 | 0.761728 |
| gnk_local_null | gnk_local_null | PHOT_CHLRE_Chen_2023__gnk_k2 |  | 0.809475 | 0.765671 | 0.865795 |
| residual_only_diagnostic | fitted_residual_only_diagnostic | PHOT_CHLRE_Chen_2023__fitted_residual_only_diagnostic |  | 2.675584 | 2.675584 |  |
