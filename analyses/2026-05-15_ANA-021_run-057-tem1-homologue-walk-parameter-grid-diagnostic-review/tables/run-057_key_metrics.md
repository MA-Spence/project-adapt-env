# RUN-057 Key Metrics

`RUN-057` is the completed live run for `EXP-021`, the TEM-1 homologue
adaptive/random walk parameter-grid diagnostic.

| Metric | Value |
| --- | ---: |
| Parameter sets | 9 |
| Alignment sequences | 768 |
| Selected homologue starts | 9 |
| Background neighborhoods | 162 |
| Walk summaries | 486 |
| Walk trajectory rows | 16,355 |
| Sampled epistasis pairs | 15,552 |
| Validation tests passed | 7/9 |
| Greedy adaptive median deceleration ratio | 0.6849 |
| Fitness vs consensus-distance Spearman | -0.9609 |
| Beneficial fraction vs consensus-distance Spearman | 0.9531 |
| Greedy total gain vs consensus-distance Spearman | 0.8609 |
| Random lag-1 autocorrelation vs consensus-distance Spearman | -0.2015 |
| Epistasis strength vs absolute epistasis Spearman | 0.9487 |
| Epistasis strength vs random-step SD Spearman | 0.4216 |
| Stability margin vs neutral fraction Spearman | 0.1927 |
| Stability margin vs deleterious fraction Spearman | -0.1376 |

## Interpretation

The homologue consensus-distance axis and adaptive-walk diagnostics behaved
coherently. The remaining failures are concentrated in the stability-margin
knob tests rather than in the walk design itself.
