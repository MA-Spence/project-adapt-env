# RUN-098 selected module-locality diagnostics

Held-out module-energy concentration for selected empirical/random/shuffled/uncoupled module-bearing comparators. Effective module counts close to 12 and entropies close to 1 indicate diffuse energy rather than localized module concentration.

| model | parameter | edge_eff_modules | edge_top10_fraction | edge_entropy | plaquette_eff_modules | plaquette_top10_fraction | plaquette_entropy |
| --- | --- | --- | --- | --- | --- | --- | --- |
| empirical_global_readout_residual | latent_phi__minus_poly_fitted_current_best | 11.547570 | 0.215110 | 0.992100 | 9.721268 | 0.284071 | 0.952094 |
| nested_softmin_empirical_module | empirical_module__posterior_rank007 | 11.727783 | 0.203749 | 0.995298 | 10.560660 | 0.260297 | 0.972360 |
| nested_softmin_random_module | random_module__posterior_rank001 | 9.389596 | 0.303430 | 0.941230 | 5.451778 | 0.506625 | 0.742538 |
| nested_softmin_shuffled_module | shuffled_module__posterior_rank002 | 10.221679 | 0.267982 | 0.963609 | 5.104343 | 0.548253 | 0.742964 |
| nested_softmin_uncoupled_module | uncoupled_module__posterior_rank011 | 11.745080 | 0.203189 | 0.995581 | 10.564001 | 0.261077 | 0.972141 |
