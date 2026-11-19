# RUN-059 Parameter-Grid Summary

| Parameter set | Ruggedness score | Terms | Strength | Covar bias | Rank | HO terms | Abs epistasis | Epistasis var | Nonzero frac | Strong frac | Sign frac | Random SD | Lag1 autocorr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| additive_control | 0 | 0 | 0 | 0 | 2 | 0 | 2.526e-12 | 1.778e-23 | 0 | 0 | 0 | 235.5 | 0.9993 |
| run057_sparse_baseline | 0.0500 | 320.0 | 0.0500 | 0 | 2 | 0 | 3.926e-03 | 3.459e-03 | 8.391e-03 | 2.894e-04 | 0 | 229.5 | 0.9994 |
| sparse_strong_pairwise | 0.2000 | 320.0 | 0.2000 | 0 | 2 | 0 | 0.0130 | 0.0325 | 8.681e-03 | 1.447e-03 | 0 | 227.7 | 0.9994 |
| dense_moderate_pairwise | 0.3062 | 3000 | 0.1000 | 0 | 2 | 0 | 0.0710 | 0.1153 | 0.0752 | 8.391e-03 | 0 | 235.3 | 0.9994 |
| dense_strong_pairwise | 0.7655 | 3000 | 0.2500 | 0 | 2 | 0 | 0.1948 | 0.7912 | 0.0891 | 9.549e-03 | 2.894e-04 | 227.5 | 0.9993 |
| dense_high_rank_pairwise | 1.5309 | 3000 | 0.2500 | 0 | 8 | 0 | 0.1724 | 0.6547 | 0.0752 | 8.970e-03 | 2.894e-04 | 228.6 | 0.9994 |
| dense_covarying_pairwise | 1.5155 | 3000 | 0.2500 | 0.8000 | 4 | 0 | 0.2809 | 1.7159 | 0.0775 | 3.472e-03 | 0 | 223.7 | 0.9994 |
| dense_covarying_strong_pairwise | 2.4249 | 3000 | 0.4000 | 0.8000 | 4 | 0 | 0.4063 | 4.0045 | 0.0709 | 2.894e-03 | 0 | 231.0 | 0.9993 |
| very_dense_covarying_high_rank_pairwise | 3.6373 | 6000 | 0.3000 | 0.8000 | 8 | 0 | 0.6172 | 4.2792 | 0.1458 | 4.051e-03 | 0 | 240.1 | 0.9993 |
| higher_order_stress | 3.8050 | 6000 | 0.3000 | 0.8000 | 8 | 5000 | 0.6497 | 4.7054 | 0.1510 | 5.208e-03 | 2.894e-04 | 234.7 | 0.9994 |

The epistasis diagnostics increase strongly with the designed ruggedness axis, but lag-1 random-walk autocorrelation remains near 0.9993 across all regimes.
