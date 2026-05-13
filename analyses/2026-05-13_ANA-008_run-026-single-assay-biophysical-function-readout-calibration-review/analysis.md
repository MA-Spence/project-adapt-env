# ANA-008: RUN-026 single-assay biophysical-function-readout calibration review

## Purpose

Assess whether the explicit stability-plus-function biophysical readout in `RUN-026` improves empirical recovery on `PHOT_CHLRE_Chen_2023` relative to the paired raw-readout controls, and determine what the result implies for `HYP-007`.

## Linked experiments/runs

- Experiments: `EXP-008`
- Runs: `RUN-026`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: `python3`

## Inputs

- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/summary.json`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/selected_panel.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/branch_validations.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/posterior_particles.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/synthetic_truth_recovery.csv`
- `experiments/2026-05-12_EXP-008_proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-12_EXP-008_proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/runs/RUN-026.yaml`

## Analysis performed

- Reviewed the reconciled run record and confirmed that `RUN-026` completed successfully on `lab-slurm` under scheduler job `95`.
- Rechecked the fixed `PHOT_CHLRE_Chen_2023` panel selection and the assay-specific `mavenn` diagnostics to determine whether the activity assay measurement layer was strong enough to interpret branch differences.
- Verified from `EXP-008` `config.yaml` that the new deterministic branches and the SMC-ABC path used `synthetic_readout_mode: stability_function` with `empirical_pairwise_target: function`, because the exported Bayesian branch labels still retain the historical `smc_abc_*_raw` suffixes.
- Compared the two raw deterministic controls, the two new biophysical-function deterministic controls, and the two Bayesian summaries on single-mutant holdout, double-mutant holdout, epistasis-prediction, KS, and reference-to-peak metrics.
- Compared the best Bayesian biophysical-function fit against the strongest raw-readout control within the same run to test the specific `HYP-007` prediction that the richer readout should improve single-mutant ranking, double-mutant ranking, and epistasis prediction together.
- Inspected the posterior parameter summary to determine whether the new readout parameters were numerically identified or simply absorbed by a diffuse posterior.
- Checked the preregistered synthetic-truth recovery outputs to determine whether any empirical failure was accompanied by a broken inverse problem.
- Cross-checked `summary.json` against the row-level `synthetic_truth_recovery.csv` because the aggregate truth-count fields in the summary artifact differ from the detailed CSV.

## Outputs

- Figures: none
- Tables: `tables/run-026_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-13_ANA-008_run-026-single-assay-biophysical-function-readout-calibration-review/tables/run-026_key_metrics.md`

## Main observations

- `RUN-026` completed successfully and produced the declared durable outputs under `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026`, including the branch comparison table, posterior particles, posterior parameter summary, and synthetic-truth recovery diagnostics.
- The empirical target was one large functional ProteinGym assay, `PHOT_CHLRE_Chen_2023`, with sequence length `118`, `167,529` measured variants total, `2,122` single mutants, and `165,407` multiple mutants. This gives substantial epistasis coverage, so the negative result is not explained by a double-mutant-starved assay.
- The measurement layer was usable but materially weaker than in the recent `SPTN1` stability runs. The assay-specific `mavenn` model reached test Spearman `0.680` and test NRMSE `0.935`, so the observation model is informative but not especially clean. That weakens confidence in very fine metric differences, but not in the larger within-run branch tradeoffs.
- The new biophysical-function deterministic branches did not dominate the paired raw-readout controls. `baseline_shared_biophysical_function_readout` improved single-mutant holdout Spearman from `0.066` to `0.216` and lowered KS from `0.807` to `0.405`, but it simultaneously collapsed double-mutant holdout Spearman from `0.319` to `0.081`, dropped epistasis-prediction Spearman from `0.468` to `0.051`, and moved the reference all the way to the fitted peak (`fraction_of_peak 1.000`, distance `0`).
- The richer predictive deterministic branch showed the same problem in a different form. `predictive_richpair_shared_biophysical_function_readout` slightly improved single-mutant holdout Spearman (`0.338` to `0.349`) and strongly improved double-mutant holdout Spearman (`0.342` to `0.551`), but epistasis-prediction Spearman fell from `0.272` to `0.095`, functional KS worsened from `0.296` to `0.357`, and the fitted reference again collapsed exactly onto the peak.
- The best Bayesian biophysical-function fit was cleaner than the deterministic biophysical branches in one important respect: it avoided the peak-collapse pathology and retained nonzero structured effects, with `functional_reference_fraction_of_peak 0.083`, `functional_reference_distance_to_peak 57`, `epistasis_strength 0.0648`, and `empirical_pairwise_strength 0.0448`.
- Even that best Bayesian fit did not satisfy the `HYP-007` prediction. Relative to the strongest raw-readout control in the same run, `predictive_richpair_shared_raw`, the best Bayesian biophysical-function fit improved double-mutant holdout Spearman from `0.342` to `0.562`, but single-mutant holdout Spearman fell from `0.338` to `0.293`, epistasis-prediction Spearman fell from `0.272` to `0.109`, and functional KS worsened from `0.296` to `0.367`. The richer readout therefore did not improve the full predictive package together.
- The posterior mean was weaker than the best particle rather than reinforcing it. Its single-mutant holdout Spearman was only `0.270`, double-mutant holdout Spearman turned negative (`-0.142`), and functional KS degraded to `0.561`. That indicates the posterior mass did not concentrate on a robustly improved solution.
- The new readout parameters were not numerically degenerate, but they were also not tightly pinned down. The posterior q90 intervals were broad for `readout_stability_midpoint` (`-1.162` to `3.716`), `readout_stability_slope` (`0.819` to `3.830`), and `readout_function_exponent` (`0.576` to `1.892`). This is consistent with a fit that can trade among many assay-readout settings without yielding a clear empirical win.
- The inverse problem itself remained well behaved. The row-level `synthetic_truth_recovery.csv` shows that both preregistered truths fell within the posterior q90 interval for all `13/13` fitted parameters. That means the negative empirical result is not simply a failure of the SMC machinery to recover parameters on matched synthetic data.
- There is a small artifact inconsistency: the aggregate synthetic-truth summary embedded in `summary.json` reports only `10` and `8` parameters within q90 for the two truths, while the detailed `synthetic_truth_recovery.csv` reports `13/13` for both. The row-level CSV is the more defensible source because it exposes every parameter-specific truth check directly.
- The scientific implication is therefore a qualified negative result. On this single activity assay, the explicit stability-plus-function readout does not deliver the predicted simultaneous gain in single-mutant ranking, double-mutant ranking, and epistasis prediction over the current raw-readout controls. The run weakens `HYP-007`, but it does not refute it outright because this is still one assay, the measurement layer is only moderate, and the implemented readout is a minimal two-trait extension rather than a richer orthogonal latent decomposition.

## Result records created

- `RES-008`

## Hypothesis updates

- `HYP-007` is weakened. The preregistered direction-of-effect prediction was not met on `PHOT_CHLRE_Chen_2023`: the new biophysical-function readout improved some metrics but degraded others, and no branch produced a clear across-the-board improvement over the best raw-readout control.
- `HYP-001` is not directly updated by this result. `EXP-008` is a single-assay diagnostic of a modified assay-readout family rather than a shared-regime multi-assay test, so it does not by itself support or weaken the broader shared-landscape claim.
