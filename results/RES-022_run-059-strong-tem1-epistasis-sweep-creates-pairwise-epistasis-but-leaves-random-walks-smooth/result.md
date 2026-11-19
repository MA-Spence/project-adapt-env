# RES-022: RUN-059 strong TEM-1 epistasis sweep creates pairwise epistasis but leaves random walks smooth

## Summary

RUN-059 passes 10/13 validation checks: stronger TEM-1 epistasis regimes increase sampled pairwise epistasis magnitude, variance, and nonzero fraction, while the consensus-distance adaptive-walk pattern from RUN-057 persists. However, random-step variance and lag-1 autocorrelation do not respond directionally to the ruggedness score, lag-1 autocorrelation remains near 0.9993 across regimes, and sign epistasis remains nearly absent even in the strongest settings. The result supports HYP-005 only for the direct epistasis-knob subcomponent and weakens the richer HYP-004 rugged-landscape interpretation.

## Generated from

- Analyses: `ANA-022`

## Relevant hypotheses

- Supports: `HYP-005`
- Weakens: `HYP-004`
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-022_run-059-tem1-ruggedness-emergence-epistasis-sweep-review/tables/run-059_key_metrics.md`
- `analyses/2026-05-15_ANA-022_run-059-tem1-ruggedness-emergence-epistasis-sweep-review/tables/run-059_validation_tests.md`
- `analyses/2026-05-15_ANA-022_run-059-tem1-ruggedness-emergence-epistasis-sweep-review/tables/run-059_parameter_grid_summary.md`
- `analyses/2026-05-15_ANA-022_run-059-tem1-ruggedness-emergence-epistasis-sweep-review/tables/run-059_ruggedness_regime_deltas.md`
- `analyses/2026-05-15_ANA-022_run-059-tem1-ruggedness-emergence-epistasis-sweep-review/tables/run-059_background_summary_by_start_stratum.md`
- `analyses/2026-05-15_ANA-022_run-059-tem1-ruggedness-emergence-epistasis-sweep-review/tables/run-059_selected_regimes_by_start_stratum.md`
- `analyses/2026-05-15_ANA-022_run-059-tem1-ruggedness-emergence-epistasis-sweep-review/tables/run-059_walk_summary_by_mode_and_stratum.md`
- `results/RES-022_run-059-strong-tem1-epistasis-sweep-creates-pairwise-epistasis-but-leaves-random-walks-smooth/metrics.json`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/summary.json`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/literature_baselines.json`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/alignment_profile_summary.json`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/homologue_starts.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/background_neighborhoods.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/walk_summaries.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/walk_trajectories.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/epistasis_pairs.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059/parameter_grid_summary.csv`
- `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059.dvc`
- `experiments/2026-05-15_EXP-022_tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/config.yaml`
- `experiments/2026-05-15_EXP-022_tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/runs/RUN-059.yaml`

## Interpretation

- `RUN-059` cleanly answers the immediate follow-up to `RES-021`: the earlier smoothness was not only because the analysis failed to record epistasis parameters. When epistasis density, strength, rank, covariation bias, and a higher-order stress condition are increased, sampled local pairwise epistasis rises strongly.
- The strongest stress regime increased mean absolute epistasis to `0.6497`, compared with `0.00393` in the RUN-059 sparse baseline and `0.00409` in the analogous RUN-057 condition. Epistasis variance increased to `4.7054`, and nonzero epistasis fraction increased to `0.1510`.
- The consensus-distance adaptive-walk axis remains intact. Farther-from-consensus TEM-1 homologues still have lower starting fitness, higher beneficial availability, larger greedy adaptive gains, and decelerating greedy walks.
- The random-walk landscape remains smooth by the metrics in this experiment. Ruggedness score did not directionally increase random-step standard deviation and did not decrease lag-1 autocorrelation; lag-1 autocorrelation stayed in the narrow range `0.99928` to `0.99943`.
- Sign epistasis is still not a substantive feature of the sampled neighbourhoods. The strongest stress regime has sign-epistasis fraction `0.000289`; this is too small to treat as biologically rich sign-epistatic ruggedness.
- Dense covariation-enriched regimes also shift the background fitness scale through large positive epistatic scores. That confirms the model can be perturbed strongly, but it makes the highest regimes better interpreted as stress tests than as calibrated TEM-1 biological settings.

## Effect on hypothesis

- `HYP-005` is supported only for the direct epistasis-knob subcomponent. The epistasis controls move epistasis magnitude, variance, and nonzero fraction directionally and with large effects.
- `HYP-004` is weakened for the richer rugged-landscape claim. The adaptive-walk and consensus-distance components remain supported, but `RUN-059` does not reproduce strong random-walk decorrelation or frequent constraining/sign epistasis expected from genuinely rugged biological landscapes.

## Limitations

- The result uses sampled local neighbourhoods and sampled mutation pairs, not exhaustive enumeration of the full 286-residue TEM-1 landscape.
- The ruggedness parameter score is a diagnostic ordering axis, not a biological measurement.
- Higher-order exact-match terms are a legacy/stress mechanism in AdaptEnv, not part of the supported public biological core.
- The strong covariation-enriched regimes alter the fitness scale, so their epistasis magnitudes should not be read as calibrated TEM-1 values.
- Consensus distance is an alignment-derived proxy for stability/robustness distance, not a direct empirical TEM-1 fitness or peak-distance measurement.

## Downstream use

- Use `RES-022` as the completed result record for `EXP-022` / `RUN-059`.
- Cite `RUN-059` as evidence that strong epistasis settings can induce measurable sampled pairwise epistasis in the TEM-1 stability-margin model.
- Do not cite `RUN-059` as evidence that the current TEM-1 configuration produces biologically rich random-walk ruggedness or frequent sign epistasis.
