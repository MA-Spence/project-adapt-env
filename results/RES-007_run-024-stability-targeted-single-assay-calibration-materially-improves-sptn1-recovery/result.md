# RES-007: RUN-024 stability-targeted single-assay calibration materially improves SPTN1 recovery

## Summary

RUN-024 switched the SPTN1_CHICK single-assay Bayesian calibration to a stability-targeted synthetic readout and materially improved holdout ranking, epistasis prediction, and KS relative to RUN-022 while keeping the fitted reference far from the peak. Because the run still targets one assay rather than a shared multi-assay regime, it does not directly support HYP-001, but it does motivate the hypothesis by showing that readout-stage misspecification was a major contributor to the earlier failure.

## Generated from

- Analyses: ANA-007

## Relevant hypotheses

- Supports: None
- Weakens: None
- Refutes: None
- Motivates: HYP-001

## Evidence

- `analyses/2026-05-12_ANA-007_run-024-single-assay-stability-readout-smc-abc-calibration-review/tables/run-024_key_metrics.md`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/summary.json`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/selected_panel.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/branch_validations.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/posterior_particles.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/synthetic_truth_recovery.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/per_assay_branch_fits.csv`
- `results/RES-006_run-022-single-assay-bayesian-smc-abc-improves-on-the-shared-fit-but-still-weakens-hyp-001/metrics.json`
- `results/RES-005_run-019-bayesian-smc-abc-recovers-nonzero-epistasis-but-still-weakens-hyp-001/metrics.json`
- `experiments/2026-05-12_EXP-007_proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/config.yaml`
- `experiments/2026-05-12_EXP-007_proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/runs/RUN-024.yaml`

## Interpretation

- `RUN-024` completed the intended paired single-assay diagnostic on the same `SPTN1_CHICK_Tsuboyama_2023_1TUD` stability assay used in `RUN-022`, but it changed the calibration semantics to fit the assay through latent stability rather than collapsed synthetic fitness.
- The measurement layer remained strong. The assay-specific `mavenn` model reached test Spearman `0.908` and test NRMSE `0.542`, so the improved empirical fit is not explained by a different or degraded observation model.
- The stability-readout change mattered immediately even before the Bayesian fit. Both deterministic stability-readout controls strongly outperformed their raw-readout counterparts on single-mutant ranking, epistasis prediction, and KS. That is direct evidence that readout-stage misspecification was a major contributor to the earlier failure.
- Those deterministic gains were not sufficient by themselves because both stability-readout controls returned to a near-peak reference solution (`fraction_of_peak` about `0.992`, distance `1` to `3`) and kept `epistasis_strength = 0.0`. They improved fit metrics, but not with a scientifically clean peak interpretation.
- The Bayesian stability-targeted fit provided the strongest balanced result. The best particle reached single-mutant holdout Spearman `0.556`, double-mutant holdout Spearman `0.483`, epistasis-prediction Spearman `0.540`, and functional KS `0.233`, while keeping the fitted reference far from the peak (`fraction_of_peak 0.039`, distance `50`) and retaining nonzero `epistasis_strength = 0.0686` plus `empirical_pairwise_strength = 0.0367`.
- Relative to the earlier raw-readout Bayesian single-assay diagnostic from `RUN-022`, the best `RUN-024` Bayesian fit improved single-mutant holdout Spearman by `+0.278`, double-mutant holdout Spearman by `+0.121`, epistasis-prediction Spearman by `+0.296`, and functional KS by `-0.103`. Relative to the earlier six-assay shared Bayesian fit from `RUN-019`, the gains were larger again.
- Relative to the earlier strongest assay-specific deterministic `SPTN1` fit from `RUN-012`, the `RUN-024` Bayesian fit also established a better assay-level ceiling. It improved single-mutant holdout Spearman from `0.406` to `0.556`, double-mutant holdout Spearman from `0.164` to `0.483`, epistasis-prediction Spearman from `0.529` to `0.540`, and functional KS from `0.414` down to `0.233`, while keeping a similarly non-pathological reference location.
- The inverse problem remained well behaved under the revised semantics. The empirical SMC round diagnostics were much tighter than in `RUN-022`, and both preregistered synthetic truths fell within the posterior q90 interval for `10/10` parameters.
- The scientific implication is therefore more positive than `RES-006`, but still narrower than direct support for `HYP-001`. `RUN-024` shows that the current single-latent family can recover the strongest prior `SPTN1` assay much more convincingly once the assay is fitted through a stability-targeted readout with pairwise structure applied at the latent stability stage. It does not yet show that one shared regime will recover multiple assays under those revised semantics.

## Effect on hypothesis

- This result motivates `HYP-001` rather than supporting it directly.
- `HYP-001` is a shared-regime, multi-assay claim. `RUN-024` still fits only one assay from one assay class, so it cannot establish that broader prediction on its own.
- It nevertheless increases the plausibility of `HYP-001` by showing that a large part of the earlier empirical mismatch was caused by where the assay readout entered the causal model, not only by generic optimizer failure or by an inability of the current family to recover any empirical target at all.

## Limitations

- The run targets only one assay, `SPTN1_CHICK_Tsuboyama_2023_1TUD`, which was already the strongest prior single-assay candidate and is therefore an optimistic diagnostic rather than a neutral panel sample.
- The assay remains from the same `Tsuboyama 2023` `cDNA display proteolysis` stability class, so the result does not yet test whether the revised semantics generalize across assay families or taxa.
- The Bayesian branch labels in the exported artifacts still retain the historical `smc_abc_*_raw` names, so correct interpretation depends on the experiment config rather than the suffixes in `summary.json`.
- The durable `RUN-024` outputs are present under `data/processed`, but they are not yet DVC-tracked in the current repository state.

## Downstream use

- Use `RES-007` as the current result record for `EXP-007` when discussing whether readout-stage misspecification explains the earlier `SPTN1` single-assay failure.
- Use this result together with `RES-005` and `RES-006` to separate three failure modes: multi-assay pooling, readout-stage misspecification, and remaining shared-regime uncertainty.
- Do not cite `RES-007` as direct support that `HYP-001` is already satisfied. The revised stability-targeted calibration still needs to be tested on a broader empirical panel before that claim is warranted.
