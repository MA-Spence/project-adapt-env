# ANA-006: RUN-022 single-assay Bayesian SMC-ABC calibration review

## Purpose

Assess whether restricting the Bayesian SMC-ABC calibration to the single SPTN1_CHICK_Tsuboyama_2023_1TUD assay rescues empirical recovery enough to explain the earlier HYP-001 failure mainly as cross-assay pooling.

## Linked experiments/runs

- Experiments: EXP-006
- Runs: RUN-022

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Inputs

- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/summary.json`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/selected_panel.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/branch_validations.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/posterior_particles.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022/synthetic_truth_recovery.csv`
- `experiments/2026-05-12_EXP-006_proteingym-single-assay-smc-abc-calibration-sptn1-chick/config.yaml`
- `experiments/2026-05-12_EXP-006_proteingym-single-assay-smc-abc-calibration-sptn1-chick/runs/RUN-022.yaml`
- `results/RES-005_run-019-bayesian-smc-abc-recovers-nonzero-epistasis-but-still-weakens-hyp-001/metrics.json`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/per_assay_branch_fits.csv`

## Analysis performed

- Reviewed the reconciled run record and confirmed that `RUN-022` completed after restoring the cached `RUN-021` assay-preparation and checkpoint outputs.
- Rechecked the fixed `SPTN1_CHICK_Tsuboyama_2023_1TUD` panel selection and the assay-specific `mavenn` diagnostics to verify that the measurement layer remained well fit on the single-assay target.
- Compared the two deterministic control branches and the two Bayesian summaries on single-mutant holdout, double-mutant holdout, epistasis-prediction, KS, and reference-to-peak metrics.
- Compared the `RUN-022` Bayesian best particle with the earlier six-assay shared Bayesian result from `RUN-019` to isolate the effect of removing cross-assay pooling.
- Compared the `RUN-022` Bayesian best particle with the earlier strongest assay-specific deterministic `SPTN1` fit from `RUN-012` to determine whether the single-assay Bayesian path established a stronger assay-level ceiling.
- Inspected the posterior parameter summary and the preregistered single-assay synthetic-truth recovery outputs to distinguish inverse-problem failure from residual model-family mismatch.

## Outputs

- Figures: none
- Tables: `tables/run-022_key_metrics.md`

## Produced artifacts

- analyses/2026-05-12_ANA-006_run-022-single-assay-bayesian-smc-abc-calibration-review/tables/run-022_key_metrics.md

## Main observations

- `RUN-022` completed successfully and produced the declared durable outputs under `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022`. The run resumed from the failed `RUN-021` cache rather than starting from a cold assay preparation.
- The empirical target is scientifically narrow and explicit: one `Tsuboyama 2023` `cDNA display proteolysis` stability assay, `SPTN1_CHICK_Tsuboyama_2023_1TUD`, with `3,201` measured variants total, `1,051` single mutants, `2,150` multiple mutants, and sequence length `60`.
- The measurement layer remained strong. The assay-specific `mavenn` model reached test Spearman `0.907` and test NRMSE `0.483`, so the downstream calibration outcome is not explained by poor raw-scale observation modeling.
- The Bayesian single-assay fit no longer collapsed to zero structured effects. The best empirical particle fit `epistasis_strength = 0.0497`, `empirical_pairwise_strength = 0.0036`, and `noise_amplitude = 0.0023`, while the posterior mean retained `epistasis_strength = 0.0344` and `empirical_pairwise_strength = 0.0400`.
- Restricting the target to one assay improved several metrics relative to the earlier six-assay Bayesian shared fit from `RUN-019`. The best Bayesian single-assay fit raised single-mutant holdout Spearman from `0.195` to `0.278`, raised double-mutant holdout Spearman from `0.119` to `0.362`, and lowered functional KS from `0.483` to `0.337`, while keeping the fitted reference far from the peak (`53` mutations away, fraction of peak `0.048`).
- Those gains were not a clean rescue. Inside `RUN-022`, the `baseline_shared_raw` control achieved higher single- and double-mutant holdout Spearman (`0.340` and `0.504`) but did so with the pathological reference-at-peak artifact (`fraction_of_peak = 1.0`, distance `0`) and zero fitted epistasis. The better-behaved `predictive_richpair_shared_raw` control still beat the Bayesian best fit on single-mutant holdout Spearman (`0.310` versus `0.278`), epistasis-prediction Spearman (`0.380` versus `0.245`), and functional KS (`0.301` versus `0.337`), while the Bayesian best fit was better only on double-mutant holdout Spearman (`0.362` versus `0.283`) and on retaining nonzero epistatic structure.
- Relative to the earlier strongest assay-specific deterministic `SPTN1` fit from `RUN-012`, the Bayesian single-assay fit also remained mixed rather than dominant. It improved double-mutant holdout Spearman (`0.362` versus `0.164`) and functional KS (`0.337` versus `0.414`), but it worsened single-mutant holdout Spearman (`0.278` versus `0.406`) and epistasis-prediction Spearman (`0.245` versus `0.529`).
- The preregistered synthetic-truth recovery on the matched single-assay scaffold remained much stronger than the empirical reconstruction. The `moderate_epistatic` truth lay within the posterior q90 interval for `10/10` parameters, and the `flatter_low_epistasis` truth did so for `8/10`, missing `empirical_pairwise_strength` and `noise_amplitude`, with best-particle distances `0.230` and `0.509`.
- The scientific implication is therefore narrower than support for `HYP-001`. Removing cross-assay pooling materially improves several empirical recovery metrics, so pooling does contribute to the earlier negative result. But even on the strongest prior assay-level candidate, the Bayesian fit still does not produce a clean, dominating reconstruction across single-mutant ranking, double-mutant ranking, epistasis prediction, and KS simultaneously. The remaining mismatch is therefore not explained by pooling alone.
- On the current evidence, `RUN-022` still weakens `HYP-001`. The weakening is indirect because the run fits only one assay rather than a shared multi-assay regime, but it still reduces the plausibility that the hypothesis would succeed on the broader panel by showing that the strongest prior candidate assay is not recovered convincingly enough after the main pooling confound is removed.

## Result records created

- `RES-006`

## Hypothesis updates

- `HYP-001` is still weakened, not refuted. `RUN-022` narrows the failure mode by showing that assay pooling is part of the problem, but not a sufficient explanation for the empirical mismatch.
