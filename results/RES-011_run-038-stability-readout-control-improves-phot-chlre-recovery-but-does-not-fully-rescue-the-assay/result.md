# RES-011: RUN-038 stability-readout control improves PHOT_CHLRE recovery but does not fully rescue the assay

## Summary

On `PHOT_CHLRE_Chen_2023`, the simpler `EXP-011` stability-readout control
materially improves single-mutant ranking, double-mutant ranking, and KS
relative to the recent PHOT Bayesian activity formulations, but epistasis
prediction remains weak and the best fit still sits too near the peak to count
as a clean rescue.

## Generated from

- Analyses: `ANA-011`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-14_ANA-011_run-038-single-assay-stability-readout-control-calibration-review/tables/run-038_key_metrics.md`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/summary.json`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/selected_panel.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/branch_validations.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/synthetic_truth_recovery.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/target_features.csv`
- `results/RES-007_run-024-stability-targeted-single-assay-calibration-materially-improves-sptn1-recovery/metrics.json`
- `results/RES-008_run-026-biophysical-function-readout-does-not-improve-phot-chlre-recovery/metrics.json`
- `results/RES-009_run-030-latent-activity-observed-fitness-readout-does-not-rescue-phot-chlre-recovery/metrics.json`
- `results/RES-010_run-036-explicit-two-latent-trait-validation-objective-does-not-rescue-phot-chlre-recovery/metrics.json`
- `experiments/2026-05-13_EXP-011_proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-13_EXP-011_proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/runs/RUN-038.yaml`

## Interpretation

- `RUN-038` completed the intended same-assay control on
  `PHOT_CHLRE_Chen_2023`, so it directly answers the unresolved comparison left
  open by `RES-010`: whether the simpler `RUN-024`-style stability-targeted
  formulation can recover the PHOT activity assay better than the more recent
  PHOT activity-oriented implementations.
- The measurement layer remained moderate and comparable to the earlier PHOT
  runs. The assay-specific `mavenn` model reached test Spearman `0.686` and
  test NRMSE `0.916`, so the control result is not explained by a different
  observation-model quality regime.
- The stability-readout semantics mattered immediately. Within the run, both
  deterministic stability-readout controls strongly outperformed the raw-readout
  controls on single-mutant holdout ranking and KS, and the richer predictive
  stability branch also improved double-mutant holdout ranking to `0.701`.
  This means the negative PHOT outcome was not robust to reverting to the
  simpler stability-targeted semantics.
- Those deterministic gains were still not a clean rescue because the
  stability-readout deterministic branches kept the fitted reference at or very
  near the peak and had weak epistasis-prediction Spearman (`-0.250` and
  `0.008`).
- The best Bayesian control fit is the strongest PHOT Bayesian result so far on
  ranking and KS. Relative to the earlier PHOT Bayesian baselines, it improved
  single-mutant holdout Spearman from `0.293/0.294/0.198` in
  `RUN-026/RUN-030/RUN-036` to `0.505`, improved double-mutant holdout
  Spearman from `0.562/0.363/0.251` to `0.721`, and lowered functional KS from
  `0.367/0.253/0.496` to `0.166`.
- Even so, the run does not fully rescue `PHOT_CHLRE`. Epistasis-prediction
  Spearman stayed weak at `0.066`, below both `RUN-026` (`0.109`) and
  `RUN-036` (`0.259`), and the best particle still sat too near the peak
  (`fraction_of_peak 0.867`, distance `4`) to count as a clean mechanistic fit.
- Relative to the earlier successful `SPTN1` stability result from `RUN-024`,
  this activity-readout control still performs less cleanly as a full package.
  `RUN-038` is slightly worse on single-mutant holdout (`0.505` versus
  `0.556`) and dramatically worse on epistasis prediction (`0.066` versus
  `0.540`) and reference geometry (`fraction_of_peak 0.867`, distance `4`
  versus `0.039`, distance `50`), even though it is better on double-mutant
  holdout and functional KS.
- The posterior mean reinforces the main tradeoff rather than removing it. It
  keeps nearly the same ranking quality and improves geometry relative to the
  best particle, but epistasis-prediction Spearman turns negative. The run
  therefore does not define one clean all-metric posterior solution.
- The inverse problem itself remained well behaved. Both preregistered
  synthetic truths fell within the posterior q90 interval for `10/10` fitted
  parameters, so the remaining PHOT weaknesses are not explained by a trivial
  SMC failure on matched synthetic data.
- The most defensible interpretation is therefore narrow. The simpler
  stability-targeted control outperforms the recent PHOT activity-oriented
  Bayesian implementations on the main ranking and KS metrics, but it still
  fails to recover strong epistasis prediction or a comfortably non-pathological
  reference-to-peak geometry. The control improves the assay fit materially
  without fully rescuing it.

## Effect on hypothesis

- `HYP-007` is weakened at the implementation level on
  `PHOT_CHLRE_Chen_2023`. The simpler stability-targeted control outperforms the
  recent PHOT activity-oriented Bayesian formulations on the main ranking and KS
  metrics, so the current richer activity/readout implementations are not yet
  delivering the predicted empirical gain on this assay.
- `HYP-001` is not directly updated. `RUN-038` is still a one-assay activity
  diagnostic rather than a shared-regime multi-assay test, and it does not
  establish the broader summary-statistics claim across assay families.

## Limitations

- The result is still based on one assay, `PHOT_CHLRE_Chen_2023`, so it
  constrains the current PHOT implementations more strongly than it constrains
  the broader space of possible multi-latent formulations.
- The assay measurement layer is only moderate, so the run is more informative
  about large tradeoffs among ranking, KS, epistasis prediction, and geometry
  than about very fine distinctions.
- `EXP-011` uses the older bootstrap summary-vector SMC distance rather than the
  later held-out `validation_objective` mode, so its posterior is not directly
  optimized against the same objective used in `RUN-036`.
- Although the stability-only control clearly improves the PHOT fit relative to
  the recent activity-oriented variants, the best particle still sits too near
  the peak and the epistasis-prediction metric remains weak.
- A good predictive control fit on an activity assay does not by itself prove
  that a stability-only causal chain is mechanistically correct for that assay.

## Downstream use

- Use `RES-011` as the project record for whether the simpler
  stability-targeted control transfers to the `PHOT_CHLRE` activity readout.
- Use `RES-011` together with `RES-007` to answer the direct scientific
  comparison: the stability-targeted formulation does improve the activity assay
  relative to the recent PHOT activity-oriented runs, but it does not perform as
  cleanly on the full package as it did on the earlier unfolding/stability
  assay.
- Use `RES-011` together with `RES-008`, `RES-009`, and `RES-010` to separate
  two PHOT claims: the recent richer PHOT implementations were not the best fit
  on this assay, but the simpler control still does not fully rescue the assay.
