# RUN-084 graphlet distance comparison

Distances are to the empirical observed `latent_phi` graphlet package. Lower is
better. Calibration-selected evidence is the registered support test; test-oracle
rows are capacity diagnostics only.

| Split | Candidate | Parameter set | Distance |
|---|---|---|---:|
| calibration | empirical readout residual | `latent_phi__minus_poly_fitted_smooth_core_only` | 0.5428 |
| calibration | soft-min/Boltzmann | `abc_best_particle` | 0.6100 |
| calibration | soft-min/Boltzmann | `abc_posterior_mean` | 0.7179 |
| calibration | fitted current | `PHOT_CHLRE_Chen_2023__fitted_current_best` | 0.7318 |
| calibration | smooth core | `PHOT_CHLRE_Chen_2023__fitted_smooth_core_only` | 0.7454 |
| calibration | fitted-current posterior mean | `PHOT_CHLRE_Chen_2023__fitted_current_posterior_mean` | 0.7839 |
| calibration | GNK local null | `PHOT_CHLRE_Chen_2023__gnk_k2` | 0.8411 |
| test | empirical readout residual | `latent_phi__minus_poly_fitted_smooth_core_only` | 0.2708 |
| test | soft-min/Boltzmann, test oracle | `abc_particle_rank_023` | 0.3964 |
| test | fitted current | `PHOT_CHLRE_Chen_2023__fitted_current_best` | 0.4016 |
| test | soft-min/Boltzmann | `abc_particle_rank_024` | 0.4029 |
| test | smooth core | `PHOT_CHLRE_Chen_2023__fitted_smooth_core_only` | 0.4119 |
| test | fitted-current posterior mean | `PHOT_CHLRE_Chen_2023__fitted_current_posterior_mean` | 0.4408 |
| test | soft-min/Boltzmann, calibration selected | `abc_best_particle` | 0.4434 |
| test | soft-min/Boltzmann posterior mean | `abc_posterior_mean` | 0.4543 |
| test | GNK local null | `PHOT_CHLRE_Chen_2023__gnk_k2` | 0.7605 |

The calibration-selected particle generalizes poorly: it beats controls on the
calibration split but is worse than fitted current and smooth core on the test
split. The one test-oracle soft-min particle is not valid held-out selection
evidence because its calibration distance is `0.7126`.
