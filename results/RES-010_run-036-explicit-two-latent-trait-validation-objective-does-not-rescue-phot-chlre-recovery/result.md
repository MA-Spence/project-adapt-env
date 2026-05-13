# RES-010: RUN-036 explicit two-latent-trait validation-objective fit does not rescue PHOT_CHLRE recovery

## Summary

On `PHOT_CHLRE_Chen_2023`, the clean `EXP-012` explicit two-latent-trait fit
improves epistasis-prediction Spearman relative to the earlier PHOT Bayesian
runs, but it still does not recover a strong overall predictive package. Single-
and double-mutant holdout ranking remain weak, functional KS worsens, and the
fitted reference stays very near the peak.

## Generated from

- Analyses: `ANA-010`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-13_ANA-010_run-036-explicit-two-latent-trait-validation-objective-calibration-review/tables/run-036_key_metrics.md`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/summary.json`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/selected_panel.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/branch_validations.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/target_features.csv`
- `results/RES-008_run-026-biophysical-function-readout-does-not-improve-phot-chlre-recovery/metrics.json`
- `results/RES-009_run-030-latent-activity-observed-fitness-readout-does-not-rescue-phot-chlre-recovery/metrics.json`
- `experiments/2026-05-13_EXP-012_proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-13_EXP-012_proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/runs/RUN-036.yaml`

## Interpretation

- `RUN-036` completed the intended clean replacement test on the same
  `PHOT_CHLRE_Chen_2023` activity assay used in `RUN-026` and `RUN-030`, so the
  historical comparison is paired at the assay level.
- Unlike the earlier amended PHOT results, this run directly tested the
  intended two-latent formulation. `EXP-012` explicitly instantiated a named
  `readout` trait block, kept stability as the second latent trait, removed the
  competing generic epistasis lever by fixing `epistasis_strength = 0`, routed
  only `empirical_pairwise_strength` into `trait:readout`, and optimized the
  held-out validation objective directly.
- The measurement layer was moderate but usable. The assay-specific `mavenn`
  model reached test Spearman `0.686` and test NRMSE `0.865`, so the result is
  not explained by a failed observation model.
- Within `RUN-036`, the Bayesian fit modestly improved the main predictive
  metrics over the deterministic two-trait controls. Relative to the stronger
  deterministic predictive branch, single-mutant holdout Spearman increased
  from `0.191` to `0.198`, double-mutant holdout Spearman from `0.216` to
  `0.251`, and epistasis-prediction Spearman from `0.240` to `0.259`.
- Those gains did not rescue the broader empirical fit. Functional KS worsened
  from `0.305` to `0.496`, and the fitted reference still remained near the
  peak (`fraction_of_peak 0.986`, distance `6`).
- The cross-run comparison is still negative overall. Relative to `RUN-030`,
  `RUN-036` substantially improved epistasis-prediction Spearman
  (`-0.026 -> 0.259`), but single-mutant holdout Spearman fell
  (`0.294 -> 0.198`), double-mutant holdout Spearman fell (`0.363 -> 0.251`),
  functional KS worsened (`0.253 -> 0.496`), and the reference remained near
  the peak.
- Relative to the earlier `RUN-026` Bayesian best fit, `RUN-036` again
  improved epistasis-prediction Spearman (`0.109 -> 0.259`) but was worse on
  single-mutant ranking (`0.293 -> 0.198`), double-mutant ranking
  (`0.562 -> 0.251`), KS (`0.367 -> 0.496`), and reference geometry
  (`fraction_of_peak 0.083 -> 0.986`, distance `57 -> 6`).
- The posterior behavior points toward structural strain rather than an unused
  epistasis lever. `empirical_pairwise_strength` was the only fitted epistatic
  term, and its posterior mean `0.0768` with `q95 0.0796` sat essentially at
  the top of the configured prior range `0.00-0.08`. The model therefore used
  the remaining pairwise freedom aggressively without finding a strong fit.
- SMC convergence was shallow. The best validation-objective distance improved
  only slightly across rounds (`23.346 -> 23.248`), and the posterior mean
  remained similar to the best particle on the main metrics. This is consistent
  with a model family that can trade metrics against each other but not satisfy
  the full predictive target well on this assay.
- The main limitation is the missing matched synthetic-truth recovery. Because
  `validation_objective` mode currently disables that scaffold, `RUN-036` cannot
  by itself determine how much of the remaining difficulty comes from inverse-
  problem properties of the new objective versus residual structural mismatch in
  the latent causal chain.
- Even with that limitation, this is the most direct PHOT test so far. The run
  removed the earlier latent-block confounders and objective mismatch, yet still
  failed to recover a strong overall predictive package. That makes the
  negative result more directly informative about the current explicit two-
  latent `PHOT_CHLRE` model family than `RES-008` or `RES-009`.

## Effect on hypothesis

- `HYP-007` is weakened more directly on `PHOT_CHLRE_Chen_2023`. Once the
  explicit second trait block, single remaining epistasis lever, and held-out
  validation objective were all implemented as intended, the model still did
  not improve the full recovery package over the prior PHOT baselines.
- `HYP-001` is not directly updated. This remains a one-assay activity
  diagnostic rather than a shared-regime multi-assay test.

## Limitations

- The result is still based on one assay, `PHOT_CHLRE_Chen_2023`, so it weakens
  the current explicit two-trait PHOT implementation more strongly than it
  falsifies the entire broader multi-latent direction.
- The assay measurement layer is only moderate, so the run is more informative
  about large tradeoffs among ranking, KS, and geometry than about very small
  metric differences.
- `validation_objective` mode currently disables matched synthetic-truth
  recovery, so this run lacks the inverse-problem control that earlier PHOT
  runs still provided.
- The paired simpler control from `EXP-011` is not yet available as a completed
  result, so the clean two-trait run can currently be compared only to
  `RUN-026` and `RUN-030`, not yet to the intended same-assay stability-only
  control.
- The exported Bayesian branch labels still retain the historical
  `smc_abc_*_raw` suffixes, so correct interpretation depends on the experiment
  config rather than the artifact labels alone.

## Downstream use

- Use `RES-010` as the project record for whether the clean explicit
  two-latent-trait, validation-objective formulation improved the single-assay
  `PHOT_CHLRE` activity fit.
- Use `RES-010` together with the amended `RES-008` and `RES-009` records to
  separate three PHOT claims: the earlier special stability-plus-function
  formulation did not rescue the assay, the public-activity reformulation also
  did not rescue it, and the later clean explicit two-trait validation-objective
  version still does not rescue the full predictive package.
