# RUN-086 Selected Held-Out Graphlet Metrics

| Candidate | SE | RSE | Local maxima | Accessibility | Lag-1 autocorr. | WHT high | GFT high |
|---|---:|---:|---:|---:|---:|---:|---:|
| Empirical observed `latent_phi` | 0.6146 | 0.2583 | 5.17 | 0.1975 | 0.1905 | 0.3788 | 0.1552 |
| Fitted current best | 0.5681 | 0.2326 | 3.50 | 0.0280 | 0.0948 | 0.4984 | 0.1844 |
| Fitted smooth core | 0.5681 | 0.2326 | 3.50 | 0.0280 | 0.0724 | 0.4984 | 0.1844 |
| Budgeted soft-min ABC best | 0.5681 | 0.2326 | 3.50 | 0.0280 | 0.0851 | 0.4984 | 0.1844 |
| ABC best, shuffled module | 0.5681 | 0.2326 | 3.50 | 0.0280 | 0.0919 | 0.4984 | 0.1844 |
| Empirical readout residual diagnostic | 0.6444 | 0.2448 | 5.25 | 0.1378 | 0.1072 | 0.4925 | 0.1462 |

The selected budgeted particle inherits the fitted-current/smooth-core graphlet
package rather than recovering the empirical topology and accessibility
pattern.

