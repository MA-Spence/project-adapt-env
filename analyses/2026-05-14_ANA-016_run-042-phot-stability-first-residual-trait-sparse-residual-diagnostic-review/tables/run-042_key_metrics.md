# RUN-042 key metrics

## Assay and diagnostic setup

| Quantity | Value |
| --- | ---: |
| Assay | `PHOT_CHLRE_Chen_2023` |
| Sequence length | 118 |
| Total measured variants | 167529 |
| Single mutants | 2122 |
| Multiple mutants | 165407 |
| MAVE-NN test Spearman | 0.68197723238683 |
| MAVE-NN test NRMSE | 0.8850599616066555 |
| Diagnostic training variants | 15997 |
| Single holdout variants | 404 |
| Double holdout variants | 43 |

## Model validation metrics

| Model | Single Spearman | Double Spearman | Epistasis Spearman | Epistasis KS | Notes |
| --- | ---: | ---: | ---: | ---: | --- |
| `stability_first_monotone` | 0.5143295454154426 | 0.5548338891206218 | -0.19771174313031462 | 0.9069767441860465 | stability proxy plus monotone readout |
| `stability_plus_photophysical_residual_trait` | 0.49721211179556724 | 0.8225429590965381 | -0.0006844543458324198 | 0.8604651162790697 | stability first, then additive residual trait and monotone 2D readout |
| `stability_plus_photophysical_trait_sparse_pair_residual` | 0.49721211179556724 | 0.8225429590965381 | -0.0006844543458324198 | 0.8604651162790697 | same held-out metrics as residual-trait model; 133 fitted pair effects |

## Key comparisons for residual-trait model

| Comparator | Delta single Spearman | Delta double Spearman | Delta epistasis Spearman | Delta epistasis KS |
| --- | ---: | ---: | ---: | ---: |
| Within-run stability-first model | -0.01711743361987539 | 0.26770906997591626 | 0.1970272887844822 | -0.046511627906976716 |
| `RUN-040` explicit product/gate two-trait best | 0.2991312117232027 | 0.5704840565932303 | -0.2649787671434405 | 0.39455602536997886 |
| `RUN-039` stability-targeted validation-objective best | -0.009120792032560776 | 0.039749249568323486 | -0.6350802280636401 | 0.11046511627906974 |
| `RUN-041` flexible monotone two-latent diagnostic | 0.02416427295235979 | 0.015871033718230176 | -0.03524448022072804 | 0.046511627906976716 |

## Interpretation

- Stability-first residualization exposes additional double-mutant ranking
  signal in `PHOT_CHLRE`.
- The additional residual trait does not recover epistasis: Spearman remains
  approximately zero and KS remains poor.
- The sparse pair-residual layer fits nonzero pair effects but has no observed
  held-out validation gain over the residual-trait model.
