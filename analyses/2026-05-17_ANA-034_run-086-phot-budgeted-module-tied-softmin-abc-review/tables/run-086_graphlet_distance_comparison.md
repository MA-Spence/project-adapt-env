# RUN-086 Graphlet Distance Comparison

| Candidate | Calibration distance | Held-out test distance |
|---|---:|---:|
| Empirical readout residual diagnostic | 0.5428 | 0.2708 |
| Fitted current best | 0.7318 | 0.4016 |
| Fitted current posterior mean | 0.7839 | 0.4408 |
| Fitted smooth core | 0.7454 | 0.4119 |
| GNK local null | 0.8411 | 0.7605 |
| Budgeted soft-min ABC best particle | 0.7240 | 0.4060 |
| Budgeted soft-min posterior mean | 0.7240 | 0.4061 |
| ABC best, no module | 0.7263 | 0.4144 |
| ABC best, random module | 0.7239 | 0.4081 |
| ABC best, shuffled module | 0.7431 | 0.4029 |
| Best posterior particle on held-out test | n/a | 0.4033 |

The registered support criterion is not met because the calibration-selected
budgeted model does not beat fitted current on held-out graphlets.

