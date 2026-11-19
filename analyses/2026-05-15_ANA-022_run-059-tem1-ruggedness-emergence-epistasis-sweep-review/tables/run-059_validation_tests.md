# RUN-059 Validation Tests

| Test | Direction | Pass | Spearman / median | Low mean | High mean | Threshold |
| --- | --- | --- | --- | --- | --- | --- |
| fitness_declines_with_consensus_distance | decrease | true | -0.6211 | 8736 | 4415 | 0.2000 |
| beneficial_fraction_increases_with_consensus_distance | increase | true | 0.8076 | 6.840e-04 | 0.0363 | 0.2000 |
| adaptive_total_gain_increases_with_consensus_distance | increase | true | 0.8383 | 518.6 | 3290 | 0.2000 |
| adaptive_walks_decelerate |  | true | 0.7387 |  |  | 0.7500 |
| random_walk_autocorrelation_declines_with_consensus_distance | decrease | true | -0.2860 | 0.9994 | 0.9993 | 0.1000 |
| epistasis_strength_increases_abs_epistasis | increase | true | 0.9355 | 0.0220 | 0.3869 | 0.2000 |
| epistasis_strength_increases_random_step_variance | increase | false | -0.0123 | 232.0 | 230.9 | 0.2000 |
| ruggedness_score_increases_abs_epistasis | increase | true | 0.9636 | 0.0220 | 0.4614 | 0.2000 |
| ruggedness_score_increases_epistasis_variance | increase | true | 0.9636 | 0.0378 | 3.4110 | 0.2000 |
| ruggedness_score_increases_nonzero_epistasis | increase | true | 0.8146 | 0.0231 | 0.1107 | 0.2000 |
| ruggedness_score_increases_sign_epistasis | increase | true | 0.4179 | 0 | 1.447e-04 | 0.2000 |
| ruggedness_score_increases_random_step_variance | increase | false | 0.1152 | 232.0 | 233.6 | 0.2000 |
| ruggedness_score_decreases_random_lag1_autocorrelation | decrease | false | 0.0303 | 0.9993 | 0.9994 | 0.1000 |

The run passed 10/13 preregistered checks. The failures are not in the consensus-distance adaptive-walk axis; they are in the expected random-walk response to epistasis/ruggedness score.
