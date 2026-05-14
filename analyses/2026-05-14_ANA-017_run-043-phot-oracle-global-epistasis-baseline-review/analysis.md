# ANA-017: RUN-043 PHOT oracle global-epistasis baseline review

## Purpose

Assess whether the `EXP-017` oracle additive latent score plus flexible
monotone global-epistasis baseline outperforms recent mechanistic PHOT variants,
thereby distinguishing restrictive AdaptEnv latent priors from scalar-assay
measurement or identifiability limits.

## Linked experiments/runs

- Experiments: `EXP-017`
- Runs: `RUN-043`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: python3

## Inputs

- `data/processed/proteingym-phot-oracle-global-epistasis-baseline/RUN-043/summary.json`
- `data/processed/proteingym-phot-oracle-global-epistasis-baseline/RUN-043/selected_panel.csv`
- `data/processed/proteingym-phot-oracle-global-epistasis-baseline/RUN-043/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-oracle-global-epistasis-baseline/RUN-043/model_validations.csv`
- `data/processed/proteingym-phot-oracle-global-epistasis-baseline/RUN-043/progress.json`
- `experiments/2026-05-14_EXP-017_proteingym-phot-oracle-global-epistasis-baseline/config.yaml`
- `experiments/2026-05-14_EXP-017_proteingym-phot-oracle-global-epistasis-baseline/runs/RUN-043.yaml`
- `scripts/proteingym_phot_structural_mismatch_diagnostic.py`
- `results/RES-013_run-039-stability-readout-validation-objective-improves-phot-recovery/metrics.json`
- `results/RES-014_run-040-explicit-two-trait-replicate-reproduces-phot-validation-objective-failure/metrics.json`
- `results/RES-015_run-041-flexible-monotone-two-latent-phot-diagnostic-improves-ranking-not-epistasis/metrics.json`
- `results/RES-016_run-042-residual-photophysical-trait-improves-double-ranking-not-epistasis/metrics.json`

## Analysis performed

- Confirmed from the run metadata that `RUN-043` completed successfully on
  `lab-slurm` as scheduler job `176` with exit code `0`.
- Confirmed that `RUN-043` is the collected run for `EXP-017`, which asks
  whether an oracle additive latent score plus flexible monotone global
  epistasis can outperform the recent PHOT mechanistic variants.
- Verified from `config.yaml` and the runner that this is a lightweight
  score-prediction diagnostic, not a full AdaptEnv landscape posterior. The
  workflow fits:
  - an additive ridge mutation score directly to the observed assay score
  - a linear additive-score prediction baseline
  - a monotone one-dimensional global-epistasis map over the additive score
- Reviewed the MAVE-NN assay diagnostics and selected panel to confirm that the
  same `PHOT_CHLRE_Chen_2023` activity/FACS assay was used.
- Compared the oracle validation metrics to `RUN-039` stability-targeted
  validation-objective control, `RUN-040` explicit product/gate two-trait
  replicate, `RUN-041` flexible two-latent diagnostic, and `RUN-042`
  stability-first residual-trait diagnostic.
- Recorded that functional KS, reference-to-peak geometry, posterior mass, and
  synthetic-truth recovery are not available for `RUN-043` because this
  diagnostic predicts observed assay scores rather than constructing a full
  AdaptEnv landscape.

## Outputs

- Figures: none
- Tables: `tables/run-043_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-14_ANA-017_run-043-phot-oracle-global-epistasis-baseline-review/tables/run-043_key_metrics.md`

## Main observations

- `RUN-043` completed in `13` seconds and produced the declared diagnostic
  outputs under
  `data/processed/proteingym-phot-oracle-global-epistasis-baseline/RUN-043`.
- The empirical target remained `PHOT_CHLRE_Chen_2023`, an activity/FACS assay
  with sequence length `118`, `167,529` measured variants, `2,122` single
  mutants, and `165,407` multiple mutants.
- The MAVE-NN observation layer was moderate, with test Spearman `0.682` and
  test NRMSE `0.872`.
- The diagnostic trained on `15,997` variants after applying the configured
  mutation-count and sampling limits, then evaluated `404` held-out single
  mutants and `43` held-out double mutants.
- The oracle additive linear model had near-zero single-mutant holdout Spearman
  (`0.005`) despite strong double-mutant holdout Spearman (`0.805`).
  Epistasis-prediction Spearman was only `0.252`, and epistasis KS was poor
  (`0.907`).
- The monotone global-epistasis map improved NRMSE and epistasis KS slightly,
  but it did not improve the rank-based readouts. Single-mutant Spearman stayed
  `0.005`, double-mutant Spearman was `0.803`, and epistasis Spearman fell
  slightly to `0.240`.
- Relative to the `RUN-039` stability-targeted best particle, the monotone
  oracle was far worse on single-mutant ranking (`-0.502`) and epistasis
  Spearman (`-0.395`), with worse epistasis KS (`+0.134`). It was only slightly
  better on double-mutant ranking (`+0.020`).
- Relative to the explicit product/gate two-trait `RUN-040` result, the
  monotone oracle improved double-mutant ranking by `0.551`, but it was worse
  on single-mutant ranking (`-0.194`), epistasis Spearman (`-0.025`), and
  epistasis KS (`+0.418`).
- Relative to the flexible/residual diagnostics (`RUN-041` and `RUN-042`), the
  oracle improved epistasis Spearman but remained weaker on single-mutant
  ranking and did not improve double-mutant ranking.
- The result fails the preregistered criterion for blaming restrictive AdaptEnv
  latent priors alone: a flexible additive-score plus global nonlinearity did
  not beat the recent mechanistic variants across the full readout package.
  The split metric behavior instead supports caution about identifiability from
  this scalar assay and about treating double-mutant ranking as sufficient
  evidence of epistasis recovery.

## Result records created

- `RES-017`

## Hypothesis updates

- `HYP-007` is weakened as a complete empirical recovery claim and as a simple
  "more flexible latent prior fixes PHOT" explanation. The oracle baseline is
  flexible in the global-epistasis sense, but it still does not recover the full
  primary metric package.
- The result does not refute all multiple-latent molecular phenotype models.
  It says the tested oracle additive/global-epistasis scalar baseline does not
  dominate, leaving measurement noise, scalar-assay non-identifiability, and
  missing stability-specific or sparse epistatic structure as live concerns.
