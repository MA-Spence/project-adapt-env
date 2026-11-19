# RUN-057 Parameter Grid Summary

| Parameter set | Stability margin | Epistasis strength | Beneficial | Neutral | Deleterious | Abs epistasis | Negative epistasis | Sign epistasis | Random step SD | Lag-1 autocorr |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| pset-000 | 4.0 | 0.0 | 0.01431 | 0.3616 | 0.6241 | 2.466e-12 | 0 | 0 | 225.3173 | 0.9988 |
| pset-001 | 8.0 | 0.0 | 0.01431 | 0.3616 | 0.6241 | 2.558e-12 | 0 | 0 | 219.8398 | 0.9988 |
| pset-002 | 12.0 | 0.0 | 0.01431 | 0.3616 | 0.6241 | 2.421e-12 | 0 | 0 | 227.2659 | 0.9987 |
| pset-003 | 4.0 | 0.025 | 0.01431 | 0.3616 | 0.6241 | 0.003099 | 0 | 0 | 220.1294 | 0.9987 |
| pset-004 | 8.0 | 0.025 | 0.01429 | 0.3617 | 0.6240 | 0.001203 | 0 | 0 | 226.1671 | 0.9988 |
| pset-005 | 12.0 | 0.025 | 0.01428 | 0.3615 | 0.6242 | 0.001805 | 0.000579 | 0 | 228.0210 | 0.9988 |
| pset-006 | 4.0 | 0.05 | 0.01432 | 0.3616 | 0.6241 | 0.004678 | 0 | 0 | 228.8497 | 0.9988 |
| pset-007 | 8.0 | 0.05 | 0.01427 | 0.3616 | 0.6241 | 0.004092 | 0 | 0 | 227.1034 | 0.9989 |
| pset-008 | 12.0 | 0.05 | 0.01431 | 0.3618 | 0.6239 | 0.004144 | 0 | 0 | 225.7479 | 0.9989 |

## Grouped knob effects

| Stability margin | Neutral | Deleterious | Beneficial | Mean beneficial effect | Random step SD | Lag-1 autocorr |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4.0 | 0.3616 | 0.6241 | 0.01432 | 160.8882 | 224.7654 | 0.9988 |
| 8.0 | 0.3616 | 0.6241 | 0.01429 | 160.9653 | 224.3701 | 0.9988 |
| 12.0 | 0.3616 | 0.6241 | 0.01430 | 161.0981 | 227.0116 | 0.9988 |

| Epistasis strength | Abs epistasis | Negative epistasis | Sign epistasis | Random step SD | Lag-1 autocorr |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.0 | 2.482e-12 | 0 | 0 | 224.1410 | 0.9988 |
| 0.025 | 0.002036 | 0.000193 | 0 | 224.7725 | 0.9988 |
| 0.05 | 0.004305 | 0 | 0 | 227.2336 | 0.9989 |

The `epistasis_strength` knob moved absolute epistasis and random-step variance.
The nominal `stability_margin` knob had negligible effect on the exported DFE
class fractions in this run.
