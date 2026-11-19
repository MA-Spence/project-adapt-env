# RUN-059 Ruggedness Regime Deltas

| Regime | Abs epistasis | Abs epistasis delta vs additive | Epistasis var | Nonzero frac | Sign frac | Random SD | Lag1 autocorr |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RUN-057 pset, margin 8/strength 0.05 | 4.092e-03 | 4.092e-03 |  |  | 0 | 227.1 | 0.9989 |
| RUN-059 sparse baseline | 3.926e-03 | 3.926e-03 | 3.459e-03 | 8.391e-03 | 0 | 229.5 | 0.9994 |
| RUN-059 higher-order stress | 0.6497 | 0.6497 | 4.7054 | 0.1510 | 2.894e-04 | 234.7 | 0.9994 |

The stress regime has much larger sampled pairwise epistasis than both the additive control and the RUN-057-like sparse baseline, but it does not produce a materially lower lag-1 autocorrelation.
