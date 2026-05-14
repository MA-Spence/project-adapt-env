# RES-016: RUN-042 residual photophysical trait improves double ranking but not PHOT epistasis

## Summary

On PHOT_CHLRE_Chen_2023, the EXP-016 stability-first residual-trait diagnostic improved double-mutant holdout ranking after fitting an additional photophysical residual trait, but it slightly reduced single-mutant ranking and did not rescue epistasis prediction. The sparse pair-residual layer learned 133 residual effects yet produced the same held-out metrics as the residual-trait model, so the run does not support sparse residual epistasis as the missing PHOT structure. The result weakens HYP-007 as a complete empirical recovery claim while preserving residual readout structure as a narrower contributor to the mismatch.

## Generated from

- Analyses: `ANA-016`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-14_ANA-016_run-042-phot-stability-first-residual-trait-sparse-residual-diagnostic-review/tables/run-042_key_metrics.md`
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

## Interpretation

- `RUN-042` completed the intended `EXP-016` stability-first residual-trait
  diagnostic on `PHOT_CHLRE_Chen_2023`.
- This is a lightweight score-prediction diagnostic, not a full AdaptEnv
  posterior. It fits an alignment-frequency stability proxy first, then fits
  an additive residual mutation score as a putative photophysical trait, then
  adds sparse pair residuals after a monotone two-dimensional surface.
- The observation layer was moderate. The assay-specific MAVE-NN model reached
  test Spearman `0.682` and test NRMSE `0.885`.
- The stability-only monotone component produced good single-mutant ranking
  (`0.514`) but weak double-mutant ranking (`0.555`) and negative epistasis
  Spearman (`-0.198`).
- Adding the photophysical residual trait increased double-mutant holdout
  Spearman to `0.823`, exceeding both `RUN-039` (`0.783`) and `RUN-041`
  (`0.807`) on that one metric. It did not improve the whole package:
  single-mutant ranking fell slightly to `0.497`, and epistasis Spearman
  remained essentially zero (`-0.001`) with poor KS (`0.860`).
- The sparse pair-residual layer fitted `133` pair effects but had identical
  held-out single, double, and epistasis metrics to the residual-trait model.
  The exported validation evidence therefore does not support sparse residual
  epistasis as the missing PHOT component.
- Relative to the explicit product/gate two-trait `RUN-040` result, `RUN-042`
  improves ranking but worsens epistasis recovery. Relative to the stability
  readout `RUN-039` control, it is close on single-mutant ranking and slightly
  better on double-mutant ranking, but much worse on epistasis prediction.
- The scientifically useful inference is that fitting stability first and then
  an additional residual trait is a valid diagnostic for ranking signal, but it
  is not a rigorous causal decomposition by itself. The residual trait is fit
  from the same scalar assay, and without orthogonal photophysical data the
  latent variable remains predictive rather than identified.

## Effect on hypothesis

- `HYP-007` is weakened as a complete empirical recovery claim. A residual
  photophysical trait can improve double-mutant ranking, but it does not
  recover PHOT epistasis, and the sparse residual layer does not validate on
  held-out metrics.
- The result does not refute all multiple-latent molecular phenotype models.
  It narrows the failure: the PHOT mismatch is not solved by the tested
  stability-first residual trait plus sparse pair residuals, and double-mutant
  ranking should not be treated as sufficient evidence of epistasis recovery.
- `HYP-001` is not directly updated because this remains a one-assay PHOT
  structural diagnostic, not a shared-regime multi-assay realism test.

## Limitations

- The run uses one assay, `PHOT_CHLRE_Chen_2023`.
- `RUN-042` is a diagnostic model, not a full AdaptEnv generator or Bayesian
  posterior, so it cannot assess functional KS, reference-to-peak geometry,
  posterior identifiability, or synthetic-truth recovery.
- The double-mutant holdout set is small (`43` variants), so the strong
  double-mutant ranking should be interpreted cautiously.
- The stability coordinate is a conservation-derived proxy, not an independent
  measured stability trait.
- The photophysical residual trait is fit from the same scalar assay after
  subtracting the stability model, so it is not causally identified without an
  orthogonal brightness, occupancy, abundance, or stability measurement.
- The sparse residual layer's lack of held-out gain may reflect limited
  overlap between fitted pair effects and held-out double mutants, regularized
  effect sizes, or true absence of recoverable sparse pair structure in this
  split. The current output cannot distinguish these explanations.

## Downstream use

- Use `RES-016` with `RES-015` to distinguish ranking improvements from
  epistasis recovery: both flexible/residual readout diagnostics help ranking,
  but neither recovers PHOT epistasis.
- Use `RES-016` with `RES-013` as evidence that the stability-targeted mapping
  remains the strongest current PHOT explanation when epistasis prediction is
  part of the readout package.
- Treat stability-first residualization as a useful diagnostic, not as a
  rigorous causal decomposition unless paired with orthogonal latent-trait
  measurements or stronger identifiability controls.
