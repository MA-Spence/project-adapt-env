# RUN-006 Key Metrics

## Panel

- Assays: 8
- Taxa: 2 Human, 2 Prokaryote, 2 Eukaryote, 2 Virus
- Assay types: 4 FACS, 2 Growth, 1 cDNA display proteolysis, 1 Antibiotics resistance
- Total single mutants: 27,249
- Sequence-length range: 60 to 364 residues

## Fitted shared parameters

- `n_functional_dims`: `4`
- `functional_sigma_base`: `12.0`
- `peak_distance_from_consensus`: `0`
- `epistasis_strength`: `0.0`
- `empirical_pairwise_strength`: `0.0`
- `noise_amplitude`: `0.0`

## Calibration validation

- Train NRMSE: `0.9917`
- Holdout NRMSE: `0.9983`
- Holdout Spearman: `0.1030`
- Functional KS distance: `0.5372`
- Reference fraction of peak: `1.0`

## Mean absolute assay-wise mismatches

| Metric | Mean absolute difference |
| --- | ---: |
| Conservation-sensitivity correlation | `0.2288` |
| Fraction above `+1 SD` | `0.0980` |
| Fraction below `-1 SD` | `0.0455` |
| 5th percentile | `0.4144` |
| 95th percentile | `0.5508` |
| Score skewness | `2.4266` |

## Best and worst assay-level aggregate mismatch

| Assay | Aggregate mismatch sum |
| --- | ---: |
| `PKN1_HUMAN_Tsuboyama_2023_1URF` | `1.693` |
| `OXDA_RHOTO_Vanella_2023_activity` | `7.623` |

## Interpretation

- The fitted shared regime recovers some coarse local statistics, but only after collapsing epistasis, pairwise structure, and observation noise to zero.
- The main failure mode is a compressed synthetic upper tail together with much stronger negative skew than most empirical assays.
- This run therefore weakens `HYP-001` rather than supporting it.
