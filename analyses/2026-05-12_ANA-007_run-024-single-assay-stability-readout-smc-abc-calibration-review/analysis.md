# ANA-007: RUN-024 single-assay stability-readout SMC-ABC calibration review

## Purpose

Assess whether the stability-targeted synthetic readout in RUN-024 materially improves empirical recovery on SPTN1_CHICK_Tsuboyama_2023_1TUD relative to the earlier raw-readout single-assay Bayesian fit and deterministic controls.

## Linked experiments/runs

- Experiments: EXP-007
- Runs: RUN-024

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Inputs

- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/summary.json`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/selected_panel.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/branch_validations.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/posterior_particles.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024/synthetic_truth_recovery.csv`
- `experiments/2026-05-12_EXP-007_proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/config.yaml`
- `experiments/2026-05-12_EXP-007_proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/runs/RUN-024.yaml`
- `results/RES-006_run-022-single-assay-bayesian-smc-abc-improves-on-the-shared-fit-but-still-weakens-hyp-001/metrics.json`
- `results/RES-005_run-019-bayesian-smc-abc-recovers-nonzero-epistasis-but-still-weakens-hyp-001/metrics.json`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/per_assay_branch_fits.csv`

## Analysis performed

- Reviewed the reconciled run record and confirmed that `RUN-024` completed successfully on `lab-slurm` under scheduler job `81`.
- Rechecked the fixed `SPTN1_CHICK_Tsuboyama_2023_1TUD` panel selection and the assay-specific `mavenn` diagnostics to verify that the measurement layer remained strong on the same single-assay target used by `EXP-006`.
- Compared the two raw deterministic controls, the two new stability-readout deterministic controls, and the two Bayesian summaries on single-mutant holdout, double-mutant holdout, epistasis-prediction, KS, and reference-to-peak metrics.
- Compared the best `RUN-024` Bayesian fit against the earlier `RUN-022` single-assay Bayesian result, the earlier `RUN-019` six-assay shared Bayesian result, and the strongest prior deterministic `SPTN1` fit from `RUN-012`.
- Inspected the posterior parameter summary and SMC round diagnostics to determine whether the stability-targeted Bayesian path still retained nonzero structured effects and whether the inverse problem became easier numerically.
- Verified from `EXP-007` `config.yaml` that the Bayesian calibration path used `synthetic_readout_mode: stability_margin` and `empirical_pairwise_target: stability`, because the output branch labels still retain the historical `smc_abc_*_raw` suffixes.
- Checked the preregistered synthetic-truth recovery outputs to determine whether the improved empirical fit coincided with a degraded inverse-problem calibration.

## Outputs

- Figures: none
- Tables: `tables/run-024_key_metrics.md`

## Produced artifacts

- analyses/2026-05-12_ANA-007_run-024-single-assay-stability-readout-smc-abc-calibration-review/tables/run-024_key_metrics.md

## Main observations

- `RUN-024` completed successfully and produced the declared durable outputs under `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024`, including the posterior particles, posterior round summaries, posterior parameter summary, and synthetic-truth recovery table in addition to the panel and branch-comparison artifacts.
- The empirical target remained the same narrow but well-defined assay used in `EXP-006`: one `Tsuboyama 2023` `cDNA display proteolysis` stability assay, `SPTN1_CHICK_Tsuboyama_2023_1TUD`, with `3,201` measured variants total, `1,051` single mutants, `2,150` multiple mutants, and sequence length `60`.
- The measurement layer remained strong. The assay-specific `mavenn` model reached test Spearman `0.908` and test NRMSE `0.542`, so the empirical calibration outcome is not explained by a failed observation model.
- The new deterministic stability-readout controls materially outperformed the old raw-readout controls on the same assay. `baseline_shared_stability_readout` raised single-mutant holdout Spearman from `0.340` to `0.501`, double-mutant holdout Spearman from `0.504` to `0.572`, epistasis-prediction Spearman from `0.097` to `0.520`, and lowered functional KS from `0.325` to `0.216`. `predictive_richpair_shared_stability_readout` similarly raised single-mutant holdout Spearman from `0.328` to `0.534`, epistasis-prediction Spearman from `0.427` to `0.517`, and lowered functional KS from `0.331` to `0.223`.
- Those deterministic gains were not clean scientific wins because both stability-readout controls reintroduced the near-peak reference pathology. Their fitted references sat at fraction-of-peak `0.992` with distances `3` and `1`, respectively, while both retained `epistasis_strength = 0.0`. That means the readout-stage change helps materially, but deterministic grid search alone still allows the fit to trade biological interpretability for easier recovery.
- The Bayesian stability-targeted fit achieved the strongest overall balance. The best empirical SMC particle reached single-mutant holdout Spearman `0.556`, double-mutant holdout Spearman `0.483`, epistasis-prediction Spearman `0.540`, and functional KS `0.233`, while keeping the fitted reference far from the peak at fraction-of-peak `0.039` and distance `50`, with nonzero `epistasis_strength = 0.0686` and `empirical_pairwise_strength = 0.0367`.
- The posterior mean was consistent with the best particle rather than contradicting it. It reached single-mutant holdout Spearman `0.554`, double-mutant holdout Spearman `0.486`, epistasis-prediction Spearman `0.523`, and functional KS `0.227`, while retaining nonzero posterior mean `epistasis_strength = 0.0382` and `empirical_pairwise_strength = 0.0417`. The posterior q90 intervals also stayed well away from a pure-zero explanation for those two parameters.
- Relative to the earlier raw-readout single-assay Bayesian fit from `RUN-022`, the best `RUN-024` Bayesian fit improved single-mutant holdout Spearman from `0.278` to `0.556`, improved double-mutant holdout Spearman from `0.362` to `0.483`, improved epistasis-prediction Spearman from `0.245` to `0.540`, and lowered functional KS from `0.337` to `0.233`, while keeping a similarly non-pathological reference location (`fraction_of_peak 0.039` versus `0.048`).
- Relative to the earlier six-assay shared Bayesian fit from `RUN-019`, the gains were larger again: single-mutant holdout Spearman improved from `0.195` to `0.556`, double-mutant holdout Spearman improved from `0.119` to `0.483`, epistasis-prediction Spearman improved from `0.265` to `0.540`, and functional KS fell from `0.483` to `0.233`.
- Relative to the strongest prior deterministic assay-specific `SPTN1` fit from `RUN-012`, the `RUN-024` Bayesian fit established a stronger assay-level ceiling. It improved single-mutant holdout Spearman from `0.406` to `0.556`, double-mutant holdout Spearman from `0.164` to `0.483`, epistasis-prediction Spearman from `0.529` to `0.540`, and lowered functional KS from `0.414` to `0.233`, while keeping a comparable reference fraction of peak (`0.039` in both cases).
- The SMC inverse problem also became easier numerically under the revised semantics. The final-round best distance fell to `8.014` with median distance `9.060`, versus `13.507` and `17.097` in `RUN-022`.
- The preregistered synthetic-truth recovery remained strong rather than degrading under the stability-targeted readout. Both `moderate_epistatic` and `flatter_low_epistasis` truths fell within the posterior q90 interval for `10/10` parameters, with best-particle distances `0.418` and `0.440`.
- The scientific implication is narrower than direct support for `HYP-001` but materially more positive than `RES-006`. On the current evidence, readout-stage misspecification was a major contributor to the earlier empirical failure: once the assay is fit through latent stability rather than collapsed synthetic fitness, the current family can recover the `SPTN1` target much more convincingly. But `RUN-024` still tests only one assay from one assay class, so it does not yet establish that a shared regime will recover multiple empirical assays without assay-specific fitting.

## Result records created

- `RES-007`

## Hypothesis updates

- `HYP-001` is motivated rather than directly supported. The revised stability-targeted readout rescues the strongest prior single-assay target far better than the earlier raw-readout path, but the shared-regime multi-assay claim remains untested under those revised semantics.
