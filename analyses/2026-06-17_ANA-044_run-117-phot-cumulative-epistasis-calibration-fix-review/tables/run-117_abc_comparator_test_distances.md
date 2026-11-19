# RUN-098 comparator test distances copied into RUN-117

| candidate_family              | candidate_model                   | parameter_set_id                                      |   distance_to_empirical_primary | selection_split   |
|:------------------------------|:----------------------------------|:------------------------------------------------------|--------------------------------:|:------------------|
| empirical_readout_residual    | empirical_global_readout_residual | latent_phi__minus_poly_fitted_current_best            |                          0.2115 | calibration       |
| random_module                 | nested_softmin_random_module      | random_module__posterior_rank001                      |                          0.4826 | calibration       |
| shuffled_module               | nested_softmin_shuffled_module    | shuffled_module__posterior_rank002                    |                          0.511  | calibration       |
| fitted_current                | fitted_current_best               | PHOT_CHLRE_Chen_2023__fitted_current_best             |                          0.521  | calibration       |
| smooth_core                   | fitted_smooth_core_only           | PHOT_CHLRE_Chen_2023__fitted_smooth_core_only         |                          0.5268 | calibration       |
| locked_first_stage            | nested_locked_first_stage_c       | RUN-084_rank001_locked_first_stage                    |                          0.5296 | calibration       |
| no_module                     | nested_softmin_no_module          | no_module__posterior_rank001                          |                          0.5372 | calibration       |
| empirical_module              | nested_softmin_empirical_module   | empirical_module__posterior_rank007                   |                          0.5432 | calibration       |
| fitted_current_posterior_mean | fitted_current_posterior_mean     | PHOT_CHLRE_Chen_2023__fitted_current_posterior_mean   |                          0.5812 | calibration       |
| uncoupled_module              | nested_softmin_uncoupled_module   | uncoupled_module__posterior_rank011                   |                          0.6401 | calibration       |
| gnk_local_null                | gnk_local_null                    | PHOT_CHLRE_Chen_2023__gnk_k2                          |                          0.8095 | calibration       |
| residual_only_diagnostic      | fitted_residual_only_diagnostic   | PHOT_CHLRE_Chen_2023__fitted_residual_only_diagnostic |                          2.676  | calibration       |
