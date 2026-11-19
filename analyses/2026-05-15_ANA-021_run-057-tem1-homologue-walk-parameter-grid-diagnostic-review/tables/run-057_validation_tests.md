# RUN-057 Validation Tests

| Test | Direction | Pass | Spearman / median | Low mean | High mean | Threshold |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| fitness_declines_with_consensus_distance | decrease | true | -0.9609 | 6983 | 2595 | 0.2 |
| beneficial_fraction_increases_with_consensus_distance | increase | true | 0.9531 | 0.000917 | 0.03960 | 0.2 |
| adaptive_total_gain_increases_with_consensus_distance | increase | true | 0.8609 | 576.9 | 3137 | 0.2 |
| adaptive_walks_decelerate | median below max | true | 0.6849 |  |  | 0.75 |
| random_walk_autocorrelation_declines_with_consensus_distance | decrease | true | -0.2015 | 0.9989 | 0.9987 | 0.1 |
| stability_margin_increases_neutral_fraction | increase | false | 0.1927 | 0.3616 | 0.3616 | 0.2 |
| stability_margin_decreases_deleterious_fraction | decrease | false | -0.1376 | 0.6241 | 0.6241 | 0.2 |
| epistasis_strength_increases_abs_epistasis | increase | true | 0.9487 | 2.482e-12 | 0.004305 | 0.2 |
| epistasis_strength_increases_random_step_variance | increase | true | 0.4216 | 224.14 | 227.23 | 0.2 |

The two failures both concern the nominal `stability_margin` sweep. They should
be interpreted together with the propagation check, not as evidence that real
TEM-1-like systems lack a stability-margin effect.
