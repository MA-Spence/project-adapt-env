# RUN-091 Calibration-Selected Family Distances

| test rank | label | calibration distance | test distance | all distance | selected parameter | base anchor | module source | posterior rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | uncoupled module | 0.3785 | 0.2237 | 0.2637 | uncoupled_module__posterior_rank006 | fitted_current_best | uncoupled_empirical_module | 6.0000 |
| 2 | random module | 0.3791 | 0.2552 | 0.2752 | random_module__posterior_rank011 | fitted_smooth_core_only | random_module | 11.0000 |
| 3 | empirical module | 0.4452 | 0.2553 | 0.3173 | empirical_module__posterior_rank001 | fitted_smooth_core_only | empirical_calibration_plaquette | 1.0000 |
| 4 | empirical readout residual | 0.5428 | 0.2708 | 0.3628 | latent_phi__minus_poly_fitted_current_best | fitted_current_best | empirical_calibration_plaquette |  |
| 5 | no module | 0.3886 | 0.3489 | 0.3409 | no_module__posterior_rank001 | fitted_smooth_core_only | none | 1.0000 |
| 6 | shuffled module | 0.4389 | 0.3535 | 0.3744 | shuffled_module__posterior_rank003 | fitted_current_best | shuffled_module | 3.0000 |
| 7 | fitted current | 0.7318 | 0.4016 | 0.4909 | PHOT_CHLRE_Chen_2023__fitted_current_best | fitted_current_best | baseline_control |  |
| 8 | smooth core | 0.7454 | 0.4119 | 0.5026 | PHOT_CHLRE_Chen_2023__fitted_smooth_core_only | fitted_smooth_core_only | baseline_control |  |
| 9 | fitted current posterior mean | 0.7839 | 0.4408 | 0.5527 | PHOT_CHLRE_Chen_2023__fitted_current_posterior_mean | fitted_current_posterior_mean | baseline_control |  |
| 10 | GNK local null | 0.7759 | 0.8230 | 0.7707 | PHOT_CHLRE_Chen_2023__gnk_k2 |  | gnk_control |  |
| 11 | residual-only diagnostic | 2.7099 | 2.6756 | 2.6407 | PHOT_CHLRE_Chen_2023__fitted_residual_only_diagnostic | fitted_residual_only_diagnostic | baseline_control |  |
