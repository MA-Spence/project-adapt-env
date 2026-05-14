# RES-018: RUN-051 OXDA paired state-function diagnostic does not improve activity recovery

## Summary

On matched OXDA_RHOTO expression and activity single-mutant assays, the paired state/function branches did not improve activity holdout recovery over the function-only branch and left the fitted reference effectively at the peak. The result weakens the current HYP-007 implementation for this OXDA diagnostic without directly updating HYP-001.

## Generated from

- Analyses: `ANA-018`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None

## Evidence

- `analyses/2026-05-14_ANA-018_run-051-oxda-expression-activity-paired-latent-trait-diagnostic-review/tables/run-051_key_metrics.md`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/summary.json`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/selected_panel.csv`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/mavenn_assay_metrics.csv`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/paired_readout_correlations.csv`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/branch_validations.csv`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051.dvc`
- `experiments/2026-05-14_EXP-018_proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/config.yaml`
- `experiments/2026-05-14_EXP-018_proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/runs/RUN-051.yaml`

## Interpretation

- `RUN-051` completed the corrected `EXP-018` OXDA paired-readout diagnostic
  on `OXDA_RHOTO_Vanella_2023_expression` and
  `OXDA_RHOTO_Vanella_2023_activity`.
- The empirical pair is appropriate for the narrow state/function question:
  expression and activity are measured on the same OXDA_RHOTO sequence, both
  are FACS single-mutant assays, and their paired-readout correlation is weak
  but nonzero (Spearman `0.246` over `6387` matched variants).
- The observation layer is a limitation. The assay-specific MAVE-NN models
  reached only moderate test Spearman: `0.567` for expression and `0.542` for
  activity.
- The function-only branch remained the best branch on activity holdout
  Spearman (`0.183`) and activity functional KS (`0.520`).
- The paired state/function single-trait branch made activity recovery worse:
  holdout Spearman fell to `0.061`, holdout NRMSE rose to `0.999`, and
  functional KS was `0.523`.
- The paired state/function explicit two-trait branch nearly tied the
  function-only branch on activity holdout Spearman (`0.180` versus `0.183`),
  but it still did not improve holdout NRMSE (`0.985` versus `0.981`) or
  functional KS (`0.526` versus `0.520`).
- All branches retained `epistasis_strength = 0.0`,
  `empirical_pairwise_strength = 0.0`, and `noise_amplitude = 0.0`.
- The reference-to-peak pathology remained. Function reference fraction of
  peak was `0.9995` to `0.9999` across branches, with distance to peak only
  `1` to `6`.
- The evidence therefore does not support the preregistered prediction that
  expression-constrained state information would improve activity recovery in
  this OXDA paired-readout diagnostic.

## Effect on hypothesis

- `HYP-007` is weakened for the current implementation. The branch designed to
  instantiate a state/function split did not improve the primary single-mutant
  activity recovery metrics and did not resolve the reference geometry problem.
- `HYP-001` is not directly updated. This experiment is a focused paired
  non-KRAS readout diagnostic, not a shared-regime multi-assay realism test.

## Limitations

- Both assays contain only single mutants, so this run does not test
  double-mutant ranking, epistasis recovery, or local epistasis structure.
- The MAVE-NN observation layers are moderate rather than strong, so a failed
  observation model remains a plausible contributor to the weak branch
  performance.
- The run compares deterministic branches, not a Bayesian posterior or
  synthetic-truth recovery scaffold.
- This is one enzyme system and one paired readout pair, so it weakens but does
  not refute the broader multiple-latent molecular-phenotype direction.

## Downstream use

- Use `RES-018` as the record for `EXP-018` / `RUN-051`.
- Do not cite this OXDA paired state/function implementation as support for
  `HYP-007`; it is negative for the branch comparison that was actually run.
- Treat the evidence as limited to single-mutant activity recovery and
  state/function readout structure.
