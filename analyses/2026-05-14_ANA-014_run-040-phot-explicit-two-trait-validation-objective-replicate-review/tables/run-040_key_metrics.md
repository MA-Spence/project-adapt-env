# RUN-040 Key Metrics

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
| Calibration variants | 2298 |
| MAVE-NN test Spearman | 0.677022 |
| MAVE-NN test NRMSE | 0.878079 |

## Branch Validation Metrics

| Branch | Single holdout rho | Double holdout rho | Epistasis rho | Functional KS | Ref fraction of peak | Ref distance |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `baseline_shared_two_trait_readout` | 0.190798 | 0.216129 | 0.239920 | 0.305783 | 0.999940 | 2 |
| `predictive_richpair_shared_two_trait_readout` | 0.161123 | 0.191804 | 0.142736 | 0.347148 | 0.999981 | 10 |
| `smc_abc_best_raw` | 0.198081 | 0.252059 | 0.264294 | 0.498231 | 0.985958 | 6 |
| `smc_abc_posterior_mean_raw` | 0.198081 | 0.251435 | 0.254912 | 0.482606 | 0.985959 | 5 |

## Best Particle Parameters

| Parameter | Value |
| --- | ---: |
| `stability_scale` | 1.062694 |
| `stability_margin` | 7.862859 |
| `blosum_blend` | 0.760450 |
| `stability_conservation_power` | 1.078672 |
| `functional_sigma_base` | 12.109340 |
| `empirical_pairwise_strength` | 0.079566 |
| Validation-objective distance | 23.228047 |

## Comparisons

| Comparison | Single rho delta | Double rho delta | Epistasis rho delta | KS delta | Ref fraction delta | Ref distance delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Best particle vs `RUN-036` explicit two-trait best | 0.000000 | 0.000624 | 0.005770 | 0.002036 | -0.000000 | 0 |
| Best particle vs `RUN-039` stability-targeted best | -0.308252 | -0.530735 | -0.370101 | 0.321780 | 0.663154 | -39 |
| Best particle vs within-run baseline deterministic | 0.007283 | 0.035930 | 0.024375 | 0.192448 | -0.013982 | 4 |
| Best particle vs within-run rich-pair deterministic | 0.036958 | 0.060255 | 0.121558 | 0.151083 | -0.014023 | -4 |

Positive deltas in Spearman metrics are improvements. Negative KS deltas and
larger positive reference-distance deltas are improvements.
