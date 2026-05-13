# RES-009: RUN-030 latent-activity observed-fitness readout does not rescue PHOT_CHLRE recovery

## Summary

On `PHOT_CHLRE_Chen_2023`, promoting activity to the public observed-fitness
readout did not improve the full empirical recovery package relative to either
the within-run activity-readout controls or the earlier `RUN-026`
biophysical-function baseline. The new formulation lowers KS, but it worsens
double-mutant recovery, epistasis prediction, and reference-to-peak geometry.

Amendment, `2026-05-13`: later code review showed that `EXP-010` still did not
fit an explicit second latent readout trait block. It used the new
observed-fitness composition path, but without `latent_trait_blocks` the
generator still fell back to one legacy functional block. This result therefore
does not directly answer the explicit two-latent-trait fitting question.
The same review also found two additional technical confounders: both generic
and empirical epistasis were routed into that same legacy functional layer, and
the Bayesian SMC stage still optimized the older bootstrap summary-vector
distance rather than the later held-out validation objective.

## Generated from

- Analyses: `ANA-009`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-13_ANA-009_run-030-latent-activity-observed-fitness-calibration-review/tables/run-030_key_metrics.md`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/summary.json`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/selected_panel.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/branch_validations.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/synthetic_truth_recovery.csv`
- `results/RES-008_run-026-biophysical-function-readout-does-not-improve-phot-chlre-recovery/metrics.json`
- `experiments/2026-05-13_EXP-010_proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-13_EXP-010_proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/runs/RUN-030.yaml`

## Interpretation

- `RUN-030` completed the intended follow-up on the same single activity assay
  used in `RUN-026`, `PHOT_CHLRE_Chen_2023`, so the comparison is paired at the
  assay level rather than confounded by a different target system.
- Subsequent code review narrowed the causal claim that can be made from this
  run. `EXP-010` changed the public observed-fitness composition path, but it
  did not configure explicit `latent_trait_blocks`, so the fitted
  `functional_sigma_base` / `n_functional_dims` parameters still referred to the
  legacy single functional block rather than to an explicit readout trait.
- The same code review found that `generic_epistasis_target: function` and
  `empirical_pairwise_target: function` were both acting on that same internal
  functional layer. `epistasis_strength` and `empirical_pairwise_strength`
  therefore behaved as partially competing levers rather than as a clean split
  between generic ruggedness and assay-matched pairwise effects.
- The measurement layer was slightly stronger than before, not weaker. The
  assay-specific `mavenn` model reached test Spearman `0.689` and test NRMSE
  `0.921`, versus `0.680` and `0.935` in `RUN-026`. The negative outcome is
  therefore not explained by a degraded observation model.
- The deterministic `EXP-010` activity-readout branches were effectively the
  same as the deterministic biophysical-function branches from `RUN-026`. The
  strongest deterministic branch again reached single-mutant holdout Spearman
  `0.349`, double-mutant holdout Spearman `0.551`, epistasis-prediction
  Spearman `0.095`, and functional KS `0.356`, while collapsing the reference
  exactly onto the peak (`fraction_of_peak 1.000`, distance `0`). This implies
  the generalized observed-fitness composition path mostly reproduces the
  earlier deterministic activity-readout behavior on this assay.
- The best Bayesian fit under the new public-activity semantics did not rescue
  the empirical package. Relative to the strongest within-run deterministic
  branch, it lowered functional KS from `0.356` to `0.253`, but single-mutant
  holdout Spearman fell from `0.349` to `0.294`, double-mutant holdout Spearman
  fell from `0.551` to `0.363`, epistasis-prediction Spearman fell from `0.095`
  to `-0.026`, and the fitted reference still remained extremely near the peak
  (`fraction_of_peak 0.9995`, distance `5`).
- The historical comparison to `RUN-026` is also negative. Relative to the
  earlier Bayesian best fit on the same assay, single-mutant holdout remained
  essentially unchanged (`0.293` to `0.294`) and functional KS improved
  (`0.367` to `0.253`), but double-mutant holdout Spearman fell sharply
  (`0.562` to `0.363`), epistasis-prediction Spearman turned negative
  (`0.109` to `-0.026`), and the reference moved from a non-pathological
  location (`fraction_of_peak 0.083`, distance `57`) back to a near-peak state
  (`0.9995`, distance `5`).
- The posterior mean was weaker than the best particle and remained near the
  peak, so the posterior did not concentrate on a robustly improved empirical
  solution.
- The inverse problem itself was still well behaved. Both preregistered truths
  fell within the posterior q90 interval for `10/10` fitted parameters in the
  row-level `synthetic_truth_recovery.csv`. This is a negative empirical result,
  not a trivial failure of the SMC machinery on matched synthetic data.
- Because `RUN-030` still used the older bootstrap summary-vector SMC distance,
  the Bayesian posterior was optimized against a target that was not identical
  to the held-out predictive metrics later used to judge success. That objective
  mismatch is consistent with the observed pattern: the run improves a coarse
  distributional metric such as KS and still recovers synthetic truths, yet
  fails on the epistasis-relevant holdout metrics that mattered most
  scientifically.

## Effect on hypothesis

- `HYP-007` is weakened again at the implementation level. `EXP-010` was
  designed as a more direct version of the activity-readout idea from
  `RUN-026`, but promoting activity to the public scalar fitness without an
  explicit fitted readout trait block, while also fitting competing epistasis
  levers under the older summary-vector objective, did not improve the full
  predictive package on `PHOT_CHLRE_Chen_2023`.
- `HYP-001` is not directly updated. This remains a one-assay activity
  diagnostic rather than a shared-regime multi-assay test.

## Limitations

- The result is still based on one assay, `PHOT_CHLRE_Chen_2023`, so it weakens
  the current activity-readout implementation more strongly than it falsifies
  the entire multi-latent direction.
- This was not yet an explicit two-latent-trait fit. The later code review
  found that the run changed the observed readout composition but still relied
  on the legacy single functional block internally.
- The run also left both generic and empirical epistasis targeting the same
  internal functional quantity, which makes the fitted epistasis controls
  harder to interpret mechanistically.
- The assay measurement layer is only moderate by recent project standards, so
  the run is more informative about large branch-level tradeoffs than about very
  fine metric differences.
- `EXP-010` changes the public scalar generator but does not yet fit explicit
  additional latent blocks directly inside the Bayesian parameterization; it
  still varies only the existing top-level landscape controls.
- The Bayesian fit used the older bootstrap summary-vector distance rather than
  the later held-out validation objective, so the fitted posterior was not
  directly optimized for the predictive metrics that defined success in the
  scientific interpretation.
- The exported Bayesian branch names still retain the historical
  `smc_abc_*_raw` suffixes, so correct interpretation depends on the experiment
  config rather than the artifact labels alone.

## Downstream use

- Use `RES-009` as the project record for whether the generalized
  observed-fitness composition path, with activity promoted to public scalar
  fitness, improved the `PHOT_CHLRE` activity fit.
- Use `RES-009` together with `RES-008` to separate two activity-assay claims:
  the earlier special stability-plus-function readout did not rescue
  `PHOT_CHLRE`, and the more direct public-activity formulation also does not
  rescue it under the current model family.
