# RUN-079 Calibration-Selected Held-Out Distances by Module Source

| module_source                   | candidate_model                   | parameter_set_id                                          | base_anchor             |   calibration_distance |   heldout_test_distance |
|:--------------------------------|:----------------------------------|:----------------------------------------------------------|:------------------------|-----------------------:|------------------------:|
| empirical_calibration_plaquette | empirical_global_readout_residual | latent_phi__minus_poly_fitted_current_best                | fitted_current_best     |               0.542832 |                0.270766 |
| none                            | coupled_no_module                 | coupled_no_module__fitted_current_best__grid0003          | fitted_current_best     |               0.627044 |                0.313741 |
| random_module                   | coupled_random_module             | coupled_random_module__fitted_current_best__grid0045      | fitted_current_best     |               0.609824 |                0.354558 |
| shuffled_module                 | coupled_shuffled_module           | coupled_shuffled_module__fitted_current_best__grid0029    | fitted_current_best     |               0.617167 |                0.39455  |
| uncoupled_empirical_module      | uncoupled_local_module            | uncoupled_local_module__fitted_smooth_core_only__grid0135 | fitted_smooth_core_only |               0.521885 |                0.395115 |
| baseline_control                | fitted_current_best               | PHOT_CHLRE_Chen_2023__fitted_current_best                 | fitted_current_best     |               0.731836 |                0.401571 |
| gnk_control                     | gnk_local_null                    | PHOT_CHLRE_Chen_2023__gnk_k2                              |                         |               0.775977 |                0.733265 |
