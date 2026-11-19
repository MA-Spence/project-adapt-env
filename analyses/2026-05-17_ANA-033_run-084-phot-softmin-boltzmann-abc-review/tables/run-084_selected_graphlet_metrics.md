# RUN-084 selected held-out graphlet metrics

| Candidate | Parameter set | SE | RSE | Local maxima | Accessibility | RW lag-1 | WHT high-order | GFT low | GFT mid | GFT high |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| empirical observed | `latent_phi` | 0.6146 | 0.2583 | 5.17 | 0.1975 | 0.1905 | 0.3788 | 0.6212 | 0.2236 | 0.1552 |
| empirical readout residual | `latent_phi__minus_poly_fitted_smooth_core_only` | 0.6444 | 0.2448 | 5.25 | 0.1378 | 0.1072 | 0.4925 | 0.5045 | 0.3493 | 0.1462 |
| fitted current | `PHOT_CHLRE_Chen_2023__fitted_current_best` | 0.5681 | 0.2326 | 3.50 | 0.0280 | 0.0948 | 0.4984 | 0.5016 | 0.3140 | 0.1844 |
| smooth core | `PHOT_CHLRE_Chen_2023__fitted_smooth_core_only` | 0.5681 | 0.2326 | 3.50 | 0.0280 | 0.0724 | 0.4984 | 0.5016 | 0.3140 | 0.1844 |
| soft-min/Boltzmann, calibration selected | `abc_best_particle` | 0.5986 | 0.2608 | 4.25 | 0.0037 | 0.0405 | 0.5268 | 0.4370 | 0.3596 | 0.2034 |
| soft-min/Boltzmann, posterior mean | `abc_posterior_mean` | 0.5691 | 0.2413 | 3.50 | 0.0238 | 0.0618 | 0.5098 | 0.4589 | 0.3518 | 0.1892 |
| soft-min/Boltzmann, test oracle | `abc_particle_rank_023` | 0.5694 | 0.2368 | 3.67 | 0.0275 | 0.0868 | 0.4972 | 0.5002 | 0.3197 | 0.1801 |
| GNK local null | `PHOT_CHLRE_Chen_2023__gnk_k2` | 0.4858 | 0.1733 | 3.25 | 0.0036 | 0.1764 | 0.1556 | 0.3864 | 0.5613 | 0.0523 |

The calibration-selected soft-min particle increases local maxima relative to
fitted current, but it destroys accessibility and does not repair the spectral
or random-walk package.
