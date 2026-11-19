# RUN-086 Assay Fit Comparison

| Rank | Candidate | Objective | Holdout Spearman | Single Spearman | Double holdout Spearman | Epistasis Spearman | KS |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | Fitted current posterior mean | 0.4631 | 0.7532 | 0.4820 | 0.5980 | -0.0511 | 0.0560 |
| 2 | Fitted smooth core | 0.4814 | 0.7573 | 0.4782 | 0.6395 | -0.1745 | 0.0562 |
| 3 | Fitted current best | 0.4856 | 0.7574 | 0.4782 | 0.6395 | -0.2042 | 0.0561 |
| 4 | Best budgeted candidate by assay objective | 0.5350 | 0.7573 | 0.4782 | 0.6395 | -0.0276 | 0.0581 |
| 16 | Budgeted soft-min posterior mean | 0.5558 | 0.7573 | 0.4782 | 0.6395 | -0.1745 | 0.0562 |
| 27 | Budgeted soft-min ABC best particle | 0.5560 | 0.7573 | 0.4782 | 0.6395 | -0.1745 | 0.0562 |

The budgeted operator does not improve the PHOT assay guardrail. Its apparent
ranking metrics are mostly inherited from the base landscape.

