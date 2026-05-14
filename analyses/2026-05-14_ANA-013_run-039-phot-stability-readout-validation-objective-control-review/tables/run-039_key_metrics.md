# RUN-039 Key Metrics

`RUN-039` is the completed live run for `EXP-013`, the PHOT stability-readout
validation-objective control on `PHOT_CHLRE_Chen_2023`.

## Assay Layer

| Metric | Value |
| --- | ---: |
| Variants loaded | 167,529 |
| Single mutants in reference | 2,122 |
| Multiple mutants in reference | 165,407 |
| MAVE-NN fit variants | 3,276 |
| MAVE-NN test variants | 310 |
| MAVE-NN test Spearman | 0.689 |
| MAVE-NN test NRMSE | 0.966 |

## Branch Validation Metrics

| Branch | Single holdout Spearman | Double holdout Spearman | Epistasis Spearman | Functional KS | Reference fraction of peak | Reference distance to peak |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| RUN-039 baseline raw | 0.066 | 0.442 | 0.519 | 0.846 | 0.284 | 12 |
| RUN-039 predictive raw | 0.338 | 0.415 | 0.397 | 0.329 | 0.079 | 84 |
| RUN-039 baseline stability readout | 0.504 | 0.538 | -0.254 | 0.183 | 1.000 | 0 |
| RUN-039 predictive stability readout | 0.505 | 0.703 | 0.032 | 0.190 | 0.987 | 3 |
| RUN-039 SMC best stability-targeted fit | 0.506 | 0.783 | 0.634 | 0.176 | 0.323 | 45 |
| RUN-039 SMC posterior-mean stability-targeted fit | 0.507 | 0.749 | -0.147 | 0.197 | 0.287 | 22 |
| RUN-036 explicit two-latent validation-objective best | 0.198 | 0.251 | 0.259 | 0.496 | 0.986 | 6 |
| RUN-038 stability-readout summary-vector best | 0.505 | 0.721 | 0.066 | 0.166 | 0.867 | 4 |

## Main Comparisons

| Comparison | Single holdout delta | Double holdout delta | Epistasis delta | Functional KS delta | Reference fraction delta | Reference distance delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| RUN-039 SMC best minus RUN-036 explicit two-latent best | +0.308 | +0.531 | +0.376 | -0.320 | -0.663 | +39 |
| RUN-039 SMC best minus RUN-038 stability summary-vector best | +0.001 | +0.062 | +0.568 | +0.010 | -0.544 | +41 |
| RUN-039 SMC best minus RUN-039 predictive stability deterministic | +0.001 | +0.080 | +0.602 | -0.014 | -0.664 | +42 |

Negative KS deltas are improvements. Negative reference-fraction deltas are
improvements when they move the fitted reference away from an implausible
near-peak solution.
