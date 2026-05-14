# RES-019: RUN-053 KCNJ2 paired state-function diagnostic does not improve function recovery

## Summary

On matched KCNJ2_MOUSE surface trafficking and ion-conduction single-mutant assays, paired state/function branches did not improve function holdout recovery over the function-only branch and left the reference near the peak. The result weakens the current HYP-007 implementation for this KCNJ2 diagnostic, with the caveat that the function observation layer itself was weak, and does not directly update HYP-001.

## Generated from

- Analyses: `ANA-019`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None

## Evidence

- `analyses/2026-05-14_ANA-019_run-053-kcnj2-surface-function-paired-latent-trait-diagnostic-review/tables/run-053_key_metrics.md`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/summary.json`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/selected_panel.csv`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/mavenn_assay_metrics.csv`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/paired_readout_correlations.csv`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/branch_validations.csv`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053.dvc`
- `experiments/2026-05-14_EXP-019_proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/config.yaml`
- `experiments/2026-05-14_EXP-019_proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/runs/RUN-053.yaml`

## Interpretation

- `RUN-053` completed the corrected `EXP-019` KCNJ2 paired-readout diagnostic
  on `KCNJ2_MOUSE_Coyote-Maestas_2022_surface` and
  `KCNJ2_MOUSE_Coyote-Maestas_2022_function`.
- The empirical pair is appropriate for the narrow state/function question:
  surface trafficking and ion-conduction function are measured on the same
  KCNJ2_MOUSE sequence, both are FACS single-mutant assays, and their
  paired-readout correlation is weak but nonzero (Spearman `0.290` over `6789`
  matched variants).
- The observation layer is a serious limitation for the function branch. The
  assay-specific MAVE-NN models reached test Spearman `0.538` for surface
  trafficking but only `0.389` for function.
- The function-only branch remained the best branch on function holdout
  Spearman (`0.066`) and function KS (`0.460`).
- The paired state/function single-trait branch was not better: holdout
  Spearman was `0.065`, holdout NRMSE rose to `1.002`, and function KS worsened
  to `0.502`.
- The paired state/function explicit two-trait branch slightly improved
  holdout NRMSE (`0.998` versus `1.000`) but worsened the more interpretable
  ranking and distributional metrics: holdout Spearman fell to `0.021`, and
  function KS worsened to `0.511`.
- All branches retained `epistasis_strength = 0.0`,
  `empirical_pairwise_strength = 0.0`, and `noise_amplitude = 0.0`.
- The reference-to-peak pathology remained. Function reference fraction of
  peak was `0.9932` to `0.9998` across branches, with distance to peak only
  `2` to `4`.
- The evidence therefore does not support the preregistered prediction that
  surface-constrained state information would improve function recovery in
  this KCNJ2 paired-readout diagnostic.

## Effect on hypothesis

- `HYP-007` is weakened for the current implementation. The paired
  state/function branches did not improve the primary single-mutant function
  recovery metrics as a package and did not resolve the reference geometry
  problem.
- `HYP-001` is not directly updated. This experiment is a focused paired
  non-KRAS readout diagnostic, not a shared-regime multi-assay realism test.

## Limitations

- Both assays contain only single mutants, so this run does not test
  double-mutant ranking, epistasis recovery, or local epistasis structure.
- The function MAVE-NN observation layer is weak, so a failed observation model
  remains a plausible contributor to the weak branch performance.
- The run compares deterministic branches, not a Bayesian posterior or
  synthetic-truth recovery scaffold.
- This is one channel system and one paired readout pair, so it weakens but
  does not refute the broader multiple-latent molecular-phenotype direction.

## Downstream use

- Use `RES-019` as the record for `EXP-019` / `RUN-053`.
- Do not cite this KCNJ2 paired state/function implementation as support for
  `HYP-007`; it is negative for the branch comparison that was actually run.
- Treat the evidence as limited to single-mutant function recovery and
  state/function readout structure.
