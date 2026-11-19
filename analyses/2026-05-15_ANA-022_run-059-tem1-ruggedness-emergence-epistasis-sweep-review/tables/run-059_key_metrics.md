# RUN-059 Key Metrics

`RUN-059` is the completed live run for `EXP-022`, the TEM-1 homologue ruggedness-emergence epistasis parameter sweep.

| Metric | Value |
| --- | --- |
| Parameter sets | 10 |
| Alignment sequences | 768 |
| Effective alignment sequences | 604.1 |
| Selected homologue starts | 9 |
| Background neighborhoods | 180 |
| Walk summaries | 540 |
| Walk trajectory rows | 25212 |
| Sampled epistasis pairs | 34560 |
| Validation tests passed | 10/13 |
| Greedy adaptive median deceleration ratio | 0.7387 |
| Fitness vs consensus-distance Spearman | -0.6211 |
| Beneficial fraction vs consensus-distance Spearman | 0.8076 |
| Adaptive total gain vs consensus-distance Spearman | 0.8383 |
| Ruggedness score vs abs epistasis Spearman | 0.9636 |
| Ruggedness score vs epistasis variance Spearman | 0.9636 |
| Ruggedness score vs nonzero epistasis Spearman | 0.8146 |
| Ruggedness score vs random-step SD Spearman | 0.1152 |
| Ruggedness score vs random lag-1 autocorr Spearman | 0.0303 |
| RUN-059 sparse baseline abs epistasis | 3.926e-03 |
| Higher-order stress abs epistasis | 0.6497 |
| Higher-order stress nonzero epistasis fraction | 0.1510 |
| Higher-order stress sign epistasis fraction | 2.894e-04 |
| Random lag-1 autocorrelation range | 0.9993 to 0.9994 |
| Random-step SD range | 223.7 to 240.1 |

## Interpretation

The strong epistasis regimes create clear sampled pairwise epistasis signal, but the random-walk diagnostics remain smooth. This separates local pairwise epistasis from global random-walk ruggedness in the current TEM-1 stability-margin setup.
