# ANA-016: RUN-042 PHOT stability-first residual-trait sparse-residual diagnostic review

## Purpose

Assess whether the `EXP-016` stability-first residual-trait and sparse
pair-residual PHOT diagnostic isolates additional readout structure beyond the
stability component and improves held-out ranking or epistasis recovery.

## Linked experiments/runs

- Experiments: `EXP-016`
- Runs: `RUN-042`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: python3

## Inputs

- `data/processed/proteingym-phot-photophysical-trait-sparse-residual-diagnostic/RUN-042/summary.json`
- `data/processed/proteingym-phot-photophysical-trait-sparse-residual-diagnostic/RUN-042/selected_panel.csv`
- `data/processed/proteingym-phot-photophysical-trait-sparse-residual-diagnostic/RUN-042/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-photophysical-trait-sparse-residual-diagnostic/RUN-042/model_validations.csv`
- `data/processed/proteingym-phot-photophysical-trait-sparse-residual-diagnostic/RUN-042/progress.json`
- `experiments/2026-05-14_EXP-016_proteingym-phot-photophysical-trait-sparse-residual-diagnostic/config.yaml`
- `experiments/2026-05-14_EXP-016_proteingym-phot-photophysical-trait-sparse-residual-diagnostic/runs/RUN-042.yaml`
- `scripts/proteingym_phot_structural_mismatch_diagnostic.py`
- `results/RES-013_run-039-stability-readout-validation-objective-improves-phot-recovery/metrics.json`
- `results/RES-014_run-040-explicit-two-trait-replicate-reproduces-phot-validation-objective-failure/metrics.json`
- `results/RES-015_run-041-flexible-monotone-two-latent-phot-diagnostic-improves-ranking-not-epistasis/metrics.json`

## Analysis performed

- Confirmed from the run metadata that `RUN-042` completed successfully on
  `lab-slurm` as scheduler job `175` with exit code `0`.
- Confirmed that `RUN-042` is the collected run for `EXP-016`, which asks
  whether, after fitting a stability-first component, an additional
  photophysical residual trait plus sparse pair residuals can isolate the
  remaining `PHOT_CHLRE` mismatch.
- Verified from `config.yaml` and the runner that this is a lightweight
  structural diagnostic, not a full AdaptEnv Bayesian posterior. The workflow
  fits:
  - a stability proxy from the alignment profile with a monotone one-dimensional
    readout
  - an additive ridge mutation score to the residual after the stability
    component
  - a monotone two-dimensional surface over stability and residual-trait
    coordinates
  - sparse pair residual effects after the two-coordinate surface
- Reviewed the MAVE-NN assay diagnostics and selected panel to confirm that the
  same `PHOT_CHLRE_Chen_2023` activity/FACS assay was used.
- Compared `model_validations.csv` to the best `RUN-039` stability-targeted
  validation-objective control, the `RUN-040` explicit product/gate two-trait
  replicate, and the `RUN-041` flexible monotone two-latent diagnostic.
- Recorded that functional KS, reference-to-peak geometry, posterior mass, and
  synthetic-truth recovery are not available for `RUN-042` because this
  diagnostic predicts observed assay scores rather than constructing a full
  AdaptEnv landscape.

## Outputs

- Figures: none
- Tables: `tables/run-042_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-14_ANA-016_run-042-phot-stability-first-residual-trait-sparse-residual-diagnostic-review/tables/run-042_key_metrics.md`

## Main observations

- `RUN-042` completed in `19` seconds and produced the declared diagnostic
  outputs under
  `data/processed/proteingym-phot-photophysical-trait-sparse-residual-diagnostic/RUN-042`.
- The empirical target remained `PHOT_CHLRE_Chen_2023`, an activity/FACS assay
  with sequence length `118`, `167,529` measured variants, `2,122` single
  mutants, and `165,407` multiple mutants.
- The MAVE-NN observation layer was moderate, with test Spearman `0.682` and
  test NRMSE `0.885`.
- The diagnostic trained on `15,997` variants after applying the configured
  mutation-count and sampling limits, then evaluated `404` held-out single
  mutants and `43` held-out double mutants.
- The stability-first monotone model reached single-mutant holdout Spearman
  `0.514` and double-mutant holdout Spearman `0.555`, but epistasis-prediction
  Spearman was negative (`-0.198`) and epistasis KS was poor (`0.907`).
- Adding the photophysical residual trait improved double-mutant holdout
  Spearman to `0.823` and improved epistasis Spearman from negative to
  approximately zero (`-0.001`), but it slightly reduced single-mutant ranking
  to `0.497` and left epistasis KS poor (`0.860`).
- The sparse pair-residual layer fitted `133` residual pair effects and used
  `1,721` photophysical mutation features, but its held-out single, double, and
  epistasis metrics were exactly the same as the residual-trait model. This
  means the sparse residual layer did not produce an observable validation
  gain on the exported held-out splits.
- Relative to the explicit product/gate two-trait `RUN-040` result, the best
  `RUN-042` residual-trait model improved single-mutant ranking by `0.299` and
  double-mutant ranking by `0.570`, but epistasis Spearman was lower by
  `0.265` and epistasis KS was worse by `0.395`.
- Relative to the stability-targeted `RUN-039` best particle, `RUN-042` was
  essentially tied on single-mutant ranking (`-0.009`) and better on
  double-mutant ranking (`+0.040`), but much worse on epistasis Spearman
  (`-0.635`) and epistasis KS (`+0.110`).
- Relative to `RUN-041`, the residual-trait model gave small ranking gains
  (`+0.024` single, `+0.016` double) but no epistasis rescue.
- The strongest interpretation is that stability-first residualization can
  expose additional ranking signal in PHOT double mutants, but the additional
  latent score and sparse residual pairs do not recover the assay's epistatic
  structure. Double-mutant ranking and epistasis recovery are therefore
  separable failure modes in this diagnostic.

## Result records created

- `RES-016`

## Hypothesis updates

- `HYP-007` is weakened as a complete empirical recovery claim. The residual
  trait improves double-mutant ranking, but the primary epistasis readout is
  still not recovered and sparse residual epistasis does not validate on the
  held-out metrics.
- The result does not refute all multi-latent causal maps. It narrows the
  mismatch: residual readout structure is real enough to improve ranking, but
  the missing PHOT structure is not isolated by this sparse pair-residual
  diagnostic.
