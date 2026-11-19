# RUN-084 locality tests

## Graphlet energy concentration on held-out graphlets

| Candidate | Parameter set | Edge energy | Edge effective graphlets | Edge entropy | Plaquette energy | Plaquette effective graphlets | Plaquette entropy |
|---|---|---:|---:|---:|---:|---:|---:|
| empirical readout residual | `latent_phi__minus_poly_fitted_smooth_core_only` | 4.02 | 9.04 | 0.938 | 15.40 | 6.72 | 0.841 |
| soft-min/Boltzmann, calibration selected | `abc_best_particle` | 41,436.92 | 2.02 | 0.443 | 179,907.93 | 1.91 | 0.410 |
| soft-min/Boltzmann, posterior mean | `abc_posterior_mean` | 30,423.82 | 1.82 | 0.398 | 130,586.20 | 1.86 | 0.400 |
| soft-min/Boltzmann, test oracle | `abc_particle_rank_023` | 22,259.98 | 1.34 | 0.187 | 92,458.33 | 1.44 | 0.211 |
| soft-min/Boltzmann, low-innovation near control | `abc_particle_rank_024` | 5,251.48 | 1.27 | 0.153 | 21,525.72 | 1.35 | 0.177 |

## Feature/pair energy concentration on held-out graphlets

| Candidate | Parameter set | Kind | Features/pairs | Total energy | Effective count | Top 25% fraction | Entropy |
|---|---|---|---:|---:|---:|---:|---:|
| empirical readout residual | `latent_phi__minus_poly_fitted_smooth_core_only` | edge | 13 | 4.02 | 7.65 | 0.648 | 0.859 |
| empirical readout residual | `latent_phi__minus_poly_fitted_smooth_core_only` | plaquette | 46 | 15.40 | 18.48 | 0.727 | 0.847 |
| soft-min/Boltzmann, calibration selected | `abc_best_particle` | edge | 13 | 41,436.92 | 5.48 | 0.829 | 0.714 |
| soft-min/Boltzmann, calibration selected | `abc_best_particle` | plaquette | 46 | 179,907.93 | 11.33 | 0.895 | 0.703 |
| soft-min/Boltzmann, test oracle | `abc_particle_rank_023` | edge | 13 | 22,259.98 | 5.44 | 0.794 | 0.687 |
| soft-min/Boltzmann, test oracle | `abc_particle_rank_023` | plaquette | 46 | 92,458.33 | 12.07 | 0.919 | 0.688 |

## Shared-feature covariance on held-out graphlets

| Candidate | Parameter set | Edge-energy similarity Spearman | Plaquette-energy similarity Spearman | Switch-edge similarity Spearman | Switch-plaquette similarity Spearman |
|---|---|---:|---:|---:|---:|
| empirical readout residual | `latent_phi__minus_poly_fitted_smooth_core_only` | -0.1998 | -0.1670 |  |  |
| soft-min/Boltzmann, calibration selected | `abc_best_particle` | 0.0663 | 0.0934 | 0.0741 | 0.0639 |
| soft-min/Boltzmann, posterior mean | `abc_posterior_mean` | 0.0046 | 0.0005 | 0.2208 | 0.2337 |
| soft-min/Boltzmann, test oracle | `abc_particle_rank_023` | -0.2236 | -0.2236 | -0.2283 | -0.2265 |

The soft-min candidates create concentrated state-transition innovation, but
the absolute innovation energies are many orders of magnitude larger than the
empirical residual and shared-feature covariance is weak or negative. This is
not empirical-like local, non-diffuse ruggedness.
