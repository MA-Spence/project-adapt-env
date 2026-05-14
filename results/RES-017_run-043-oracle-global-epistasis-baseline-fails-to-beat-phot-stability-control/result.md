# RES-017: RUN-043 oracle global-epistasis baseline fails to beat the PHOT stability control

## Summary

On PHOT_CHLRE_Chen_2023, the EXP-017 oracle additive latent score plus flexible monotone global-epistasis baseline did not beat the recent mechanistic PHOT controls on the full prespecified metric package. The oracle models had near-zero single-mutant holdout Spearman, strong double-mutant ranking near 0.80, and only modest epistasis Spearman near 0.24-0.25 with poor KS. This fails the preregistered condition for blaming restrictive AdaptEnv latent priors alone and instead points back toward scalar-assay measurement limits, insufficient identifiability, or missing stability-specific/epistatic structure.

## Generated from

- Analyses: `ANA-017`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-14_ANA-017_run-043-phot-oracle-global-epistasis-baseline-review/tables/run-043_key_metrics.md`
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

## Interpretation

- `RUN-043` completed the intended `EXP-017` oracle global-epistasis diagnostic
  on `PHOT_CHLRE_Chen_2023`.
- This is a lightweight score-prediction diagnostic, not a full AdaptEnv
  posterior. It fits an additive ridge mutation score directly to the assay
  score, then compares that linear additive score against a flexible monotone
  one-dimensional global-epistasis map.
- The observation layer was moderate. The assay-specific MAVE-NN model reached
  test Spearman `0.682` and test NRMSE `0.872`.
- The oracle did not satisfy the pre-experiment prediction. The monotone
  global-epistasis model had near-zero single-mutant holdout Spearman
  (`0.005`), double-mutant Spearman `0.803`, epistasis Spearman `0.240`, and
  poor epistasis KS (`0.884`).
- The linear additive oracle was similar: single-mutant Spearman `0.005`,
  double-mutant Spearman `0.805`, epistasis Spearman `0.252`, and epistasis KS
  `0.907`.
- Compared with the `RUN-039` stability-targeted best particle, the oracle is
  much weaker on single-mutant ranking and epistasis prediction. The monotone
  oracle was lower by `0.502` on single-mutant Spearman and lower by `0.395`
  on epistasis Spearman, while epistasis KS was worse by `0.134`. It was only
  slightly better on double-mutant Spearman (`+0.020`).
- Compared with the `RUN-040` explicit product/gate two-trait result, the
  oracle improves double-mutant ranking but is worse on single-mutant ranking,
  epistasis Spearman, and epistasis KS.
- Compared with `RUN-041` and `RUN-042`, the oracle has better epistasis
  Spearman, but the single-mutant failure is severe and double-mutant ranking
  does not exceed those diagnostics.
- The strongest defensible inference is negative and diagnostic: the PHOT
  mismatch is not simply solved by allowing an additive latent score plus a
  flexible monotone global-epistasis nonlinearity. That pushes interpretation
  away from "AdaptEnv latent priors are too restrictive" as a sufficient
  explanation and back toward scalar-assay identifiability, measurement noise,
  or missing stability-specific and sparse epistatic structure.

## Effect on hypothesis

- `HYP-007` is weakened as a complete empirical recovery claim. The oracle
  baseline gives the scalar assay a flexible additive-to-observed map, yet it
  fails to beat the current PHOT stability-targeted control across the full
  primary readout package.
- The result does not refute all multiple-latent molecular phenotype models.
  It specifically weakens the claim that a flexible scalar global-epistasis
  layer is enough to expose a better latent mapping for this PHOT assay.
- `HYP-001` is not directly updated because this remains a one-assay PHOT
  structural diagnostic, not a shared-regime multi-assay realism test.

## Limitations

- The run uses one assay, `PHOT_CHLRE_Chen_2023`.
- `RUN-043` is a diagnostic score model, not a full AdaptEnv generator or
  Bayesian posterior, so it cannot assess functional KS, reference-to-peak
  geometry, posterior identifiability, or synthetic-truth recovery.
- The double-mutant holdout set is small (`43` variants), so double-mutant
  ranking differences should be interpreted cautiously.
- The additive latent score is fit from the same scalar assay and is not an
  independently measured molecular phenotype.
- The poor single-mutant ranking may reflect assay noise, split-specific
  limitations, dominance of multi-mutant training signal, or genuine
  non-identifiability of a scalar additive latent from this assay. The current
  diagnostic does not distinguish those explanations.

## Downstream use

- Use `RES-017` as the oracle check showing that a flexible additive-score plus
  monotone global-epistasis baseline does not dominate the mechanistic PHOT
  variants.
- Use `RES-017` with `RES-013` to keep the stability-targeted mapping as the
  strongest current PHOT explanation when single-mutant ranking and epistasis
  prediction are both considered.
- Use `RES-017` with `RES-015` and `RES-016` to avoid overinterpreting
  double-mutant ranking alone: several diagnostics can rank doubles well
  without recovering single-mutant ranking and epistasis together.
