# RUN-086 Locality Tests

| Quantity | Calibration | Held-out test |
|---|---:|---:|
| ABC best graphlet edge innovation energy | 0.0 | 6.3e-18 |
| ABC best graphlet plaquette innovation energy | 0.0 | 2.0e-17 |
| ABC best edge effective module count | n/a | 1.00 |
| ABC best plaquette effective module count | n/a | 1.00 |
| ABC best feature-edge effective module count | n/a | 3.65 |
| ABC best feature-plaquette effective module count | n/a | 4.80 |
| Shared-feature versus edge-energy similarity Spearman | 0.4334 | 0.3851 |
| Shared-feature versus plaquette-energy similarity Spearman | 0.4334 | 0.3851 |

The energy concentration numbers are not interpretable as biological locality
because the absolute innovation energy is essentially zero. The budget removes
the `RUN-084` over-innovation failure by collapsing the operator onto the base,
not by producing empirical-like local ruggedness.

