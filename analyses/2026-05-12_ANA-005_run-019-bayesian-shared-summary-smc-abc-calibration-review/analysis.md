# ANA-005: RUN-019 Bayesian shared-summary SMC-ABC calibration review

## Purpose

Assess whether the resumed Bayesian SMC-ABC calibration in RUN-019 recovers the fixed six-assay ProteinGym stability panel more faithfully than the deterministic controls and whether synthetic-truth recovery supports the inference procedure.

## Linked experiments/runs

- Experiments: EXP-005
- Runs: RUN-019

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Inputs

- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/summary.json`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/selected_panel.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/mavenn_assay_metrics.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/branch_validations.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/posterior_particles.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/posterior_rounds.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/posterior_parameter_summary.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/synthetic_truth_recovery.csv`
- `experiments/2026-05-09_EXP-005_proteingym-shared-summary-smc-abc-calibration-panel/config.yaml`
- `experiments/2026-05-09_EXP-005_proteingym-shared-summary-smc-abc-calibration-panel/runs/RUN-019.yaml`
- `results/RES-004_run-012-modestly-improves-shared-raw-calibration-but-still-weakens-hyp-001/result.md`
- `results/RES-004_run-012-modestly-improves-shared-raw-calibration-but-still-weakens-hyp-001/metrics.json`

## Analysis performed

- Reviewed the reconciled run record and confirmed that `RUN-019` completed after resuming the stored `RUN-017` checkpoint.
- Rechecked the paired six-assay ProteinGym stability panel and the assay-specific `mavenn` diagnostics to verify that the measurement layer remained well fit.
- Compared the two deterministic control branches and the two Bayesian shared summaries on single-mutant holdout, double-mutant holdout, epistasis-prediction, KS, and reference-to-peak metrics.
- Inspected the posterior parameter summary and round diagnostics to determine whether the Bayesian fitter still collapsed to zero epistasis, pairwise structure, or noise.
- Evaluated preregistered synthetic-truth recovery on the same panel scaffold to distinguish inverse-problem failure from model-family mismatch.

## Outputs

- Figures: none
- Tables: `tables/run-019_key_metrics.md`

## Produced artifacts

- analyses/2026-05-12_ANA-005_run-019-bayesian-shared-summary-smc-abc-calibration-review/tables/run-019_key_metrics.md

## Main observations

- `RUN-019` completed successfully and produced the declared durable outputs under `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019`, including the posterior particles, posterior round summaries, posterior parameter summary, and synthetic-truth recovery table in addition to the panel and branch-comparison artifacts.
- The run stayed paired to the same empirical target used by `EXP-003`, `EXP-004`, and the interrupted `RUN-017`: `6` short `Tsuboyama 2023` `cDNA display proteolysis` stability assays spanning `4` taxa, `23,279` measured variants total, `6,424` single mutants, `16,855` multiple mutants, and sequence lengths from `44` to `72` residues.
- The measurement layer remained strong. Across the six assays, the assay-specific `mavenn` models reached mean test Spearman `0.890` and mean test NRMSE `0.480`, so the Bayesian outcome is not explained by a degraded observation model.
- The Bayesian fit no longer collapsed to zero structured effects. The best empirical SMC particle fit `epistasis_strength = 0.0527`, `empirical_pairwise_strength = 0.0008`, and `noise_amplitude = 0.00044`, while the posterior mean retained `epistasis_strength = 0.0429`, `empirical_pairwise_strength = 0.0050`, and `noise_amplitude = 0.0052`.
- Relative to the strongest deterministic control inside `RUN-019`, `predictive_richpair_shared_raw`, the best Bayesian shared fit improved double-mutant holdout Spearman from `0.098` to `0.119` and epistasis-prediction Spearman from `0.094` to `0.265`, and it avoided the near-peak reference artifact by moving the fitted reference fraction of peak from `0.998` down to `0.049`.
- Those gains did not rescue empirical recovery. The best Bayesian fit still reached only single-mutant holdout Spearman `0.195`, double-mutant holdout NRMSE `1.226`, functional KS `0.483`, and train NRMSE `1.081`, while the deterministic control remained better on single-mutant holdout ranking (`0.243` versus `0.195`) with essentially identical functional KS.
- The preregistered synthetic-truth recovery on the matched panel scaffold succeeded much better than the empirical reconstruction. The `moderate_epistatic` truth fell within the posterior q90 interval for `10/10` parameters, and the `flatter_low_epistasis` truth did so for `9/10`, missing only `noise_amplitude`, with best-particle distances `2.503` and `4.206`.
- The scientific implication is narrower than a positive validation of `HYP-001`. `RUN-019` shows that the Bayesian inference procedure itself can recover nonzero epistatic structure on model-matched synthetic targets and no longer fails by trivial collapse, but the shared empirical fit on the real six-assay panel still does not recover the local landscape statistics convincingly enough to support the hypothesis.
- On the current evidence, `RUN-019` still weakens `HYP-001`. It is stronger negative evidence than `RUN-017` because the resumed run completed the full empirical SMC and the preregistered synthetic-truth checks, yet the main empirical reconstruction failure persisted.

## Result records created

- `RES-005`

## Hypothesis updates

- `HYP-001` is further weakened, not refuted. The result narrows the failure mode away from pure optimizer collapse and toward a mismatch between the shared model or summary target and the empirical panel.
