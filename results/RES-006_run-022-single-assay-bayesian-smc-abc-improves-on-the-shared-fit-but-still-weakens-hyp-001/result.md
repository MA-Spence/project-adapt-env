# RES-006: RUN-022 single-assay Bayesian SMC-ABC improves on the shared fit but still weakens HYP-001

## Summary

RUN-022 restricted the Bayesian SMC-ABC fit to the single SPTN1_CHICK_Tsuboyama_2023_1TUD assay and improved several metrics relative to the earlier six-assay shared fit, but it still did not recover the assay cleanly enough to support HYP-001, so the result still weakens the hypothesis while narrowing the failure mode.

## Generated from

- Analyses: ANA-006

## Relevant hypotheses

- Supports: None
- Weakens: HYP-001
- Refutes: None

## Evidence

- `analyses/2026-05-12_ANA-006_run-022-single-assay-bayesian-smc-abc-calibration-review/tables/run-022_key_metrics.md`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/summary.json`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/selected_panel.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/branch_validations.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/posterior_particles.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/synthetic_truth_recovery.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/per_assay_branch_fits.csv`
- `results/RES-005_run-019-bayesian-smc-abc-recovers-nonzero-epistasis-but-still-weakens-hyp-001/metrics.json`
- `experiments/2026-05-12_EXP-006_proteingym-single-assay-smc-abc-calibration-sptn1-chick/config.yaml`
- `experiments/2026-05-12_EXP-006_proteingym-single-assay-smc-abc-calibration-sptn1-chick/runs/RUN-022.yaml`

## Interpretation

- `RUN-022` completed the intended single-assay diagnostic after restoring cached outputs from the failed `RUN-021` attempt. The empirical target was one `Tsuboyama 2023` `cDNA display proteolysis` stability assay, `SPTN1_CHICK_Tsuboyama_2023_1TUD`, with `3,201` measured variants including `2,150` multiple mutants.
- The measurement layer remained strong. The assay-specific `mavenn` model reached test Spearman `0.907` and test NRMSE `0.483`, so the downstream fit is not explained by a failure of the raw-scale observation model.
- The Bayesian single-assay fit retained nonzero structured effects. The best empirical particle fit `epistasis_strength = 0.0497` and `empirical_pairwise_strength = 0.0036`, while keeping the fitted reference far from the peak (`fraction_of_peak = 0.048`, distance `53`).
- Relative to the earlier six-assay shared Bayesian fit from `RUN-019`, the best `RUN-022` Bayesian fit improved single-mutant holdout Spearman from `0.195` to `0.278`, improved double-mutant holdout Spearman from `0.119` to `0.362`, and lowered functional KS from `0.483` to `0.337`. That means cross-assay pooling was materially contributing to the earlier failure.
- Those gains still did not rescue the scientific claim. Within `RUN-022`, the naive `baseline_shared_raw` control achieved higher single- and double-mutant holdout Spearman (`0.340` and `0.504`) but only by collapsing to the pathological reference-at-peak solution. The better-behaved `predictive_richpair_shared_raw` control still beat the Bayesian best fit on single-mutant holdout Spearman (`0.310` versus `0.278`), epistasis-prediction Spearman (`0.380` versus `0.245`), and functional KS (`0.301` versus `0.337`), while the Bayesian fit was better only on double-mutant holdout Spearman and on retaining nonzero epistasis.
- Relative to the earlier strongest assay-specific deterministic `SPTN1` fit from `RUN-012`, the Bayesian single-assay fit improved double-mutant holdout Spearman (`0.362` versus `0.164`) and functional KS (`0.337` versus `0.414`) but worsened single-mutant holdout Spearman (`0.278` versus `0.406`) and epistasis-prediction Spearman (`0.245` versus `0.529`). The run therefore does not establish a clean new assay-level ceiling.
- The preregistered single-assay synthetic-truth recovery remained much stronger than the empirical reconstruction. The `moderate_epistatic` truth fell within the posterior q90 interval for `10/10` parameters, and the `flatter_low_epistasis` truth did so for `8/10`, missing only `empirical_pairwise_strength` and `noise_amplitude`. This supports the claim that the remaining failure is not just an optimizer-collapse artifact.
- The main implication is that pooling across assays is part of the failure mode, but not a sufficient explanation. Even after removing the pooling confound and focusing on the strongest prior candidate assay, the model family and current summary-targeted inverse problem still do not recover the local landscape statistics convincingly enough to support `HYP-001`.

## Effect on hypothesis

- This result still weakens `HYP-001`, but indirectly rather than as a direct test of the full shared-regime claim.
- `HYP-001` is about a shared synthetic parameter regime matching multiple empirical assays. `RUN-022` fits only one assay, so it cannot by itself support that broader claim.
- It nevertheless weakens `HYP-001` because the strongest prior assay-level candidate still is not recovered convincingly enough after the main pooling confound is removed, which lowers the plausibility that a shared regime would succeed across the wider panel.

## Limitations

- The run targets only one assay and therefore does not directly test the multi-assay sharing requirement in `HYP-001`.
- The assay was chosen because it was the strongest prior deterministic candidate, which makes this a relatively optimistic diagnostic rather than a neutral sample from the panel.
- The Bayesian objective remains summary-targeted rather than a full variant-level likelihood, so residual failure could still reflect information discarded by the summaries as well as mismatch in the model family.
- The within-run branch comparison remains structurally ambiguous: different branches trade off holdout ranking, epistasis prediction, KS, and reference-to-peak behavior, so no single scalar metric fully resolves which fit is scientifically preferable.
- The synthetic-truth recovery is informative only on model-matched single-assay targets and therefore cannot establish that the model family is adequate for the real `SPTN1` assay.

## Downstream use

- Use `RES-006` as the current record of the `EXP-006` single-assay diagnostic when discussing whether pooling alone explains the remaining `HYP-001` failure.
- Use this result to constrain future claims: the evidence now supports saying that removing assay pooling helps materially, but it still does not justify saying that the current single-latent Adapt-Env family recovers even the strongest prior `SPTN1` assay convincingly.
