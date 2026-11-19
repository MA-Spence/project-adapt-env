# RUN-084 SMC posterior summary

| Round | Epsilon | Attempts | Best distance | Median distance |
|---:|---:|---:|---:|---:|
| 0 | 0.8258 | 160 | 0.7267 | 0.7717 |
| 1 | 0.8404 | 160 | 0.6494 | 0.7650 |
| 2 | 0.7800 | 160 | 0.6583 | 0.7442 |

| Selected parameter | Value |
|---|---:|
| Selected particle | `abc_best_particle` |
| Total ABC distance | 0.6583 |
| Calibration graphlet distance | 0.6100 |
| Assay guardrail distance | 0.3220 |
| `tau` | 1.0854 |
| `log_tau` | 0.0819 |
| Auxiliary correlation | 0.6708 |
| Auxiliary scale | 0.6741 |
| Auxiliary offset SD | 1.1962 |
| Auxiliary pairwise scale | 0.0464 |
| Auxiliary pairwise density | 0.0400 |
| Mean values-to-`a` correlation | 0.9449 |
| Mean edge-gradient correlation to `a` | 0.9505 |
| Mean innovation-to-`a` SD ratio | 0.1853 |
| Mean auxiliary-weight SD | 0.2923 |
| Total regularization penalty | 0.0000 |

The posterior did not show strong monotone SMC contraction. The best distance
improved from round 0 to round 1, then the final round best was slightly worse
than round 1 while the median distance improved modestly.
