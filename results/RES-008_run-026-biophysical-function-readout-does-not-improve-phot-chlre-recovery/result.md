# RES-008: RUN-026 biophysical-function readout does not improve PHOT_CHLRE recovery

## Summary

On the single functional ProteinGym assay `PHOT_CHLRE_Chen_2023`, the explicit stability-plus-function biophysical readout did not outperform the paired raw-readout controls across the combined package of single-mutant holdout, double-mutant holdout, epistasis prediction, and KS metrics. The result therefore weakens `HYP-007` rather than supporting it, and it does not directly update `HYP-001`.

Amendment, `2026-05-13`: later code review showed that `EXP-008` did not yet instantiate an explicit fitted latent-trait block for the assay readout. It tested a richer special readout layered onto the legacy single functional block, so it should be interpreted as a negative result for that implementation, not as a clean direct test of the intended explicit two-latent-trait fitting question.

## Generated from

- Analyses: `ANA-008`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-13_ANA-008_run-026-single-assay-biophysical-function-readout-calibration-review/tables/run-026_key_metrics.md`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/summary.json`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/selected_panel.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/branch_validations.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/synthetic_truth_recovery.csv`
- `experiments/2026-05-12_EXP-008_proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-12_EXP-008_proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/runs/RUN-026.yaml`

## Interpretation

- `RUN-026` completed the intended single-assay diagnostic on `PHOT_CHLRE_Chen_2023`, a large functional assay with `2,122` single mutants and `165,407` multiple mutants, so the experiment had enough multi-mutant coverage to test whether the richer readout improved epistasis-relevant recovery.
- Subsequent code review narrowed what this run actually tested. `EXP-008` used `synthetic_readout_mode = stability_function`, but it did not supply explicit `latent_trait_blocks`; the generator therefore still fell back to one legacy FGM-style functional block, and the fitted `functional_sigma_base` / `n_functional_dims` parameters did not represent an explicit second typed latent trait.
- The measurement layer was usable but not especially strong. The assay-specific `mavenn` model reached test Spearman `0.680` and test NRMSE `0.935`, which is adequate for broad branch comparisons but weakens confidence in very small metric differences.
- The deterministic biophysical-function branches did not dominate the paired raw-readout controls. They improved some metrics, especially double-mutant holdout in the richer branch, but both biophysical deterministic fits collapsed the reference exactly onto the fitted peak (`fraction_of_peak 1.000`, distance `0`), which is scientifically pathological.
- The best Bayesian biophysical-function fit avoided that peak-collapse pathology and retained nonzero structured effects, but it still did not produce the predicted across-the-board empirical win. Relative to the strongest raw-readout control in the same run, `predictive_richpair_shared_raw`, the best Bayesian fit improved double-mutant holdout Spearman from `0.342` to `0.562`, but single-mutant holdout Spearman fell from `0.338` to `0.293`, epistasis-prediction Spearman fell from `0.272` to `0.109`, and functional KS worsened from `0.296` to `0.367`.
- The posterior mean was weaker than the best particle, with double-mutant holdout Spearman turning negative (`-0.142`) and functional KS worsening to `0.561`. That indicates the posterior did not concentrate on a robustly better empirical solution.
- The inverse problem itself was not broken. The row-level `synthetic_truth_recovery.csv` shows both preregistered truths within the posterior q90 interval for `13/13` fitted parameters, so the negative empirical result is not explained by a trivial failure of the SMC machinery on matched synthetic data.
- There is one artifact inconsistency that matters for interpretation: the aggregate synthetic-truth block embedded in `summary.json` underreports the q90 recovery counts, while the row-level `synthetic_truth_recovery.csv` reports `13/13` for both truths. The CSV is the authoritative source because it exposes the parameter-level checks directly.

## Effect on hypothesis

- `HYP-007` is weakened, but more narrowly than the original summary implied. The specific implemented direction tested here, a richer stability-plus-function readout layered onto the legacy functional block, did not improve empirical recovery relative to the collapsed/raw-readout controls on an activity assay with substantial multi-mutant coverage.
- `HYP-001` is not directly updated. This experiment tests one modified single-assay readout family on one functional assay, not the broader shared-regime claim across multiple empirical systems.

## Limitations

- The result is still based on one assay, `PHOT_CHLRE_Chen_2023`, so it weakens `HYP-007` as a general modeling direction but does not refute it.
- This was not yet a clean explicit two-latent-trait fit. The later code review found that the run left the generator on the legacy single functional block rather than fitting a named readout trait block directly.
- The assay measurement layer is only moderate by recent project standards (`mavenn` test Spearman `0.680`), so the run is more informative about large branch-level tradeoffs than about fine distinctions among near-tied fits.
- The implemented biophysical-function readout is a minimal two-trait extension rather than a richer latent decomposition of multiple molecular phenotypes.
- The exported Bayesian branch names still retain the historical `smc_abc_*_raw` suffixes, so correct interpretation depends on the experiment config rather than the artifact labels alone.

## Downstream use

- Use `RES-008` as the project record for whether the current explicit stability-plus-function readout improved empirical recovery on the first functional ProteinGym activity assay tested under `EXP-008`.
- Use this result together with `RES-007` to separate two claims that are easy to conflate: stability-targeted readouts improved the recent `SPTN1` stability assay, but the first explicit function-readout extension did not transfer that gain cleanly to a functional activity assay.
