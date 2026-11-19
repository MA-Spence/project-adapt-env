# ANA-022: RUN-059 TEM-1 ruggedness emergence epistasis sweep review

## Purpose

Assess whether RUN-059 shows emergence of ruggedness under stronger TEM-1 epistasis parameterizations, and determine how the result updates HYP-004 and HYP-005 relative to RUN-057.

## Linked experiments/runs

- Experiments: EXP-022
- Runs: RUN-059

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-059_key_metrics.md`
- `tables/run-059_validation_tests.md`
- `tables/run-059_parameter_grid_summary.md`
- `tables/run-059_ruggedness_regime_deltas.md`
- `tables/run-059_background_summary_by_start_stratum.md`
- `tables/run-059_selected_regimes_by_start_stratum.md`
- `tables/run-059_walk_summary_by_mode_and_stratum.md`

## Inputs

- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/summary.json`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/literature_baselines.json`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/alignment_profile_summary.json`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/homologue_starts.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/background_neighborhoods.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/walk_summaries.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/walk_trajectories.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/epistasis_pairs.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/parameter_grid_summary.csv`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/parameter_grid_summary.csv`
- `experiments/2026-05-15_EXP-022_tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/config.yaml`
- `experiments/2026-05-15_EXP-022_tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/runs/RUN-059.yaml`
- `results/RES-021_run-057-tem1-homologue-walks-support-consensus-axis-dynamics-but-expose-stability-margin-knob-failure/result.md`

## Analysis performed

- Confirmed from run metadata and labproj reconciliation that `RUN-059` completed successfully as the live `EXP-022` execution on `lab-slurm` with scheduler job `206` and exit code `0`.
- Reviewed the experiment question under `HYP-004` and `HYP-005`: whether stronger epistasis density, strength, rank, covariation bias, and higher-order stress settings cause a TEM-1 homologue stability-margin landscape to become measurably rugged.
- Checked that the run repeated the `EXP-021` TEM-1 beta-lactamase setup using the reviewed UniProt `P62593` reference and a 768-sequence MMSeqs/ColabFold alignment profile with `604.1` effective sequences.
- Summarized the run across `10` explicit parameter regimes, `9` selected homologue starts, `180` sampled backgrounds, `540` walk summaries, `25,212` walk-trajectory rows, and `34,560` sampled epistasis pairs.
- Recomputed compact review tables from the exported CSV/JSON files using only the Python standard library.
- Compared the RUN-059 sparse baseline against the analogous RUN-057 `stability_margin = 8`, `epistasis_strength = 0.05` condition.

## Main observations

- The run passed `10/13` validation checks. The consensus-distance adaptive-walk checks passed, and the direct pairwise epistasis checks passed; the failed checks are the expected random-walk response to stronger epistasis.
- The TEM-1 homologue consensus-distance axis remains coherent. Farther-from-consensus starts have lower mean fitness and higher beneficial availability, and greedy adaptive-walk total gain increases with consensus distance.
- Greedy adaptive walks still decelerate, but the median deceleration ratio is close to the preregistered limit (`0.7387` against a maximum of `0.75`).
- Stronger epistasis regimes create clear local pairwise ruggedness in the sampled neighbourhoods. Ruggedness score correlates strongly with mean absolute epistasis (`Spearman 0.9636`), epistasis variance (`0.9636`), and nonzero epistasis fraction (`0.8146`).
- The higher-order stress regime raises mean absolute epistasis from effectively zero in the additive control and `0.00393` in the RUN-057-like sparse baseline to `0.6497`, with epistasis variance `4.7054` and nonzero epistasis fraction `0.1510`.
- Sign epistasis remains nearly absent even in the strongest settings. The higher-order stress regime has sign-epistasis fraction `0.000289`, so the positive validation check for sign epistasis is statistically directional but biologically weak in magnitude.
- Random-walk ruggedness does not emerge in the same way. Ruggedness score has only weak correlation with random-step standard deviation (`Spearman 0.1152`, failed) and does not decrease lag-1 autocorrelation (`Spearman 0.0303`, failed). Lag-1 autocorrelation remains between `0.99928` and `0.99943` across all parameter regimes.
- Dense covariation-enriched regimes substantially shift the background fitness scale through large positive epistatic scores. This is evidence that the knobs can perturb the landscape, but it is also a caveat: the strongest regimes are stress tests of model behaviour rather than calibrated biological TEM-1 parameterizations.

## Interpretation

RUN-059 shows where ruggedness emerges in this implementation: it appears first as measurable pairwise epistasis in sampled local neighbourhoods, not as a strongly corrugated random-walk landscape. Increasing density, rank, covariation targeting, and higher-order stress terms moves epistasis magnitude, variance, and nonzero fraction by large amounts relative to the RUN-057 smooth baseline. However, unconditional random walks remain almost perfectly autocorrelated, and sign epistasis remains extremely rare.

The result therefore separates three concepts that should not be conflated in the record:

- Direct pairwise epistasis is now clearly present under stronger settings.
- Adaptive-walk opportunity still follows the TEM-1 consensus-distance axis.
- Biological-style ruggedness in the sense of random-walk decorrelation, frequent sign epistasis, or strongly constraining peak-proximal interactions remains weak.

## Result records created

- `RES-022`

## Hypothesis updates

- `HYP-005` is supported only for the direct epistasis-knob subcomponent: stronger and denser epistasis terms produce the expected directional increases in sampled epistasis magnitude, variance, and nonzero fraction.
- `HYP-004` is weakened for the richer rugged-landscape interpretation. The adaptive-walk consensus-axis behaviour remains intact, but RUN-059 does not show the random-walk decorrelation or substantial sign/constraining epistasis expected from strongly rugged biological landscapes.
