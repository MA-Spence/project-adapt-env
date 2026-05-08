# RUN-008 Key Metrics

## Experiment size

- Grid settings: `81`
- Seeds per setting: `4`
- Landscapes evaluated: `324`
- Sequence length: `80`
- Sampled double mutants per landscape: `10,000`

## Global summary

- Fraction of settings where doubles are more deleterious than singles: `1.0000`
- Fraction of settings where doubles are more lethal than singles: `0.8333`
- Highest mean epistasis setting: `stability_margin-5__functional_sigma_base-28__n_functional_dims-2__epistasis_strength-0.06`
- Lowest mean epistasis setting: `stability_margin-12__functional_sigma_base-28__n_functional_dims-8__epistasis_strength-0`

## Setting ranges

| Metric | Minimum | Maximum |
| --- | ---: | ---: |
| Single neutral fraction | `0.9020` | `1.0000` |
| Single lethal fraction | `0.0000` | `0.0092` |
| Double lethal fraction | `0.0001` | `0.3285` |
| Double beneficial fraction | `0.0000` | `0.0404` |
| Mean absolute epistasis | `0.0095` | `3.1256` |

## Stability-margin means

| Stability margin | Single neutral | Single lethal | Double neutral | Double lethal | Mean absolute epistasis |
| --- | ---: | ---: | ---: | ---: | ---: |
| `5.0` | `0.9757` | `0.0092` | `0.6107` | `0.3285` | `3.0815` |
| `8.0` | `0.9870` | `0.0000` | `0.9295` | `0.0212` | `0.2557` |
| `12.0` | `0.9870` | `0.0000` | `0.9537` | `0.0001` | `0.0496` |

## Functional-sigma means

| Functional sigma base | Single neutral | Single beneficial | Double beneficial |
| --- | ---: | ---: | ---: |
| `12.0` | `0.9594` | `0.0065` | `0.0128` |
| `20.0` | `0.9942` | `0.0000` | `0.0046` |
| `28.0` | `0.9962` | `0.0000` | `0.0047` |

## Epistasis-strength means

| Epistasis strength | Single neutral | Mean absolute epistasis | Positive epistasis fraction | Negative epistasis fraction |
| --- | ---: | ---: | ---: | ---: |
| `0.00` | `0.9832` | `1.1149` | `0.0004` | `0.1232` |
| `0.03` | `0.9832` | `1.1283` | `0.0043` | `0.1277` |
| `0.06` | `0.9832` | `1.1434` | `0.0122` | `0.1353` |

## Interpretation

- The uncalibrated family spans a broad envelope of local DFE and epistasis regimes, especially for double mutants.
- `stability_margin` is the clearest harshness control in this scan.
- `epistasis_strength` changes double-mutant epistasis much more than it changes the single-mutant DFE.
- This is a model-capacity result, not an empirical realism result.
