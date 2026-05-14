# RUN-043 key metrics

## Assay and diagnostic setup

| Quantity | Value |
| --- | ---: |
| Assay | `PHOT_CHLRE_Chen_2023` |
| Sequence length | 118 |
| Total measured variants | 167529 |
| Single mutants | 2122 |
| Multiple mutants | 165407 |
| MAVE-NN test Spearman | 0.6816352371202583 |
| MAVE-NN test NRMSE | 0.8717570070126405 |
| Diagnostic training variants | 15997 |
| Single holdout variants | 404 |
| Double holdout variants | 43 |

## Model validation metrics

| Model | Single Spearman | Double Spearman | Epistasis Spearman | Epistasis KS | Notes |
| --- | ---: | ---: | ---: | ---: | --- |
| `oracle_additive_linear` | 0.004563680311934349 | 0.8051948051948052 | 0.2519929622734269 | 0.9069767441860465 | additive latent score fit directly to assay score |
| `oracle_additive_monotone_global_epistasis` | 0.004563680311934349 | 0.8029624320695554 | 0.23950468136514647 | 0.8837209302325582 | monotone global-epistasis map over additive score |

## Key comparisons for monotone oracle model

| Comparator | Delta single Spearman | Delta double Spearman | Delta epistasis Spearman | Delta epistasis KS |
| --- | ---: | ---: | ---: | ---: |
| `RUN-039` stability-targeted validation-objective best | -0.5017692235161937 | 0.020168722541340833 | -0.3948910923526612 | 0.13372093023255816 |
| `RUN-040` explicit product/gate two-trait best | -0.19351721976043024 | 0.5509035295662477 | -0.0247896314324616 | 0.41781183932346727 |
| `RUN-041` flexible monotone two-latent diagnostic | -0.4684841585312731 | -0.003709493308752476 | 0.20494465549025084 | 0.06976744186046513 |
| `RUN-042` stability-first residual-trait diagnostic | -0.4926484314836329 | -0.019580527026982653 | 0.24018913571097888 | 0.023255813953488413 |

## Monotone map effect relative to linear oracle

| Metric | Monotone minus linear |
| --- | ---: |
| Single Spearman | 0.0 |
| Single NRMSE | -0.0304278038187189 |
| Double Spearman | -0.0022323731252498202 |
| Double NRMSE | -0.04036687292113661 |
| Epistasis Spearman | -0.012488280908280414 |
| Epistasis KS | -0.023255813953488302 |

## Interpretation

- The oracle additive/global-epistasis baseline does not beat the recent PHOT
  mechanistic variants across the full metric package.
- The baseline can rank held-out doubles, but it essentially fails on
  held-out single mutants and remains weak on epistasis recovery.
- This fails the preregistered condition for concluding that restrictive
  AdaptEnv latent priors alone explain the PHOT structural mismatch.
