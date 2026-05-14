# RUN-041 Key Metrics

## Assay and Observation Layer

| Quantity | Value |
| --- | ---: |
| Assay | `PHOT_CHLRE_Chen_2023` |
| Taxon | `Eukaryote` |
| Selection type | `FACS` |
| Coarse selection type | `Activity` |
| Sequence length | 118 |
| Total variants | 167529 |
| Single mutants | 2122 |
| Multiple mutants | 165407 |
| MAVE-NN test Spearman | 0.657395 |
| MAVE-NN test NRMSE | 0.947080 |

## Diagnostic Configuration

| Quantity | Value |
| --- | ---: |
| Mode | `flexible_monotone_two_latent` |
| Train variants | 15997 |
| Single holdout variants | 404 |
| Double holdout variants | 43 |
| Max train mutation count | 5 |
| Max train variants | 24000 |
| Ridge alpha | 2.0 |
| 2D monotone bins | 16 |
| Max epistasis pairs | 30000 |

## Model Validation Metrics

| Model | Single rho | Single NRMSE | Double rho | Double NRMSE | Epistasis rho | Epistasis KS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `flexible_monotone_stability_activity_surface` | 0.473048 | 1.208746 | 0.806672 | 0.716887 | 0.034560 | 0.813953 |

## Comparisons

| Comparison | Single rho delta | Double rho delta | Epistasis rho delta | Epistasis KS delta |
| --- | ---: | ---: | ---: | ---: |
| Versus `RUN-040` explicit product/gate two-trait best | 0.274967 | 0.554613 | -0.229734 | 0.348044 |
| Versus `RUN-039` stability-targeted validation-objective best | -0.033285 | 0.023878 | -0.599836 | 0.063953 |

Positive Spearman deltas are improvements. Negative KS deltas are improvements.
`RUN-041` does not produce functional KS or reference-to-peak geometry because
it is a score-prediction diagnostic rather than a full landscape posterior.
