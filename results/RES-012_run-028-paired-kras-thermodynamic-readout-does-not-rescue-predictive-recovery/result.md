# RES-012: RUN-028 paired KRAS thermodynamic readout improves some diagnostics but does not rescue predictive recovery

## Summary

On matched KRAS abundance and DARPin K55 binding assays, the paired thermodynamic binding readout improves functional KS and epistasis-prediction Spearman relative to raw controls, but it degrades single- and double-mutant holdout ranking, worsens abundance agreement, and leaves the fitted reference essentially at the peak. The result therefore weakens the current HYP-007 implementation rather than supporting a clean multi-latent rescue.

## Generated from

- Analyses: `ANA-012`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None

## Evidence

- `analyses/2026-05-14_ANA-012_run-028-paired-kras-abundance-binding-thermodynamic-diagnostic-review/tables/run-028_key_metrics.md`
- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/summary.json`
- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/selected_panel.csv`
- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/mavenn_assay_metrics.csv`
- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/branch_validations.csv`
- `experiments/2026-05-13_EXP-009_proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/config.yaml`
- `experiments/2026-05-13_EXP-009_proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/runs/RUN-028.yaml`

## Interpretation

- `RUN-028` completed the intended `EXP-009` matched KRAS diagnostic on
  `RASK_HUMAN_Weng_2022_abundance` and
  `RASK_HUMAN_Weng_2022_binding-DARPin_K55`.
- The observation layer was strong, especially for the binding assay. The
  assay-specific MAVE-NN models reached test Spearman `0.818` for abundance
  and `0.933` for DARPin K55 binding, so the mixed result is not primarily
  explained by a failed observation model.
- The binding-only thermodynamic branch improved over the binding-only raw
  branch on binding KS (`0.384` versus `0.506`), double-mutant holdout
  Spearman (`0.298` versus `0.195`), epistasis-prediction Spearman (`0.172`
  versus `0.106`), and epistasis KS (`0.481` versus `0.673`). This is real
  evidence that the thermodynamic readout captures some relevant structure.
- The full paired abundance-plus-binding thermodynamic branch produced the
  strongest epistasis-prediction signal in the run: epistasis-prediction
  Spearman `0.287` versus `0.104` for the paired raw branch, and functional KS
  `0.258` versus `0.381`.
- Those gains did not satisfy the preregistered prediction. The full paired
  thermodynamic branch had worse abundance KS than the paired raw branch
  (`0.832` versus `0.582`), much worse single-mutant holdout Spearman (`0.105`
  versus `0.265`), and worse double-mutant holdout Spearman (`0.144` versus
  `0.252`).
- The reference-to-peak pathology remained unresolved. All four branches placed
  the fitted reference at or near the peak, with `fraction_of_peak` between
  `0.9986` and `1.000` and distance to peak between `0` and `6`.
- The full paired thermodynamic branch pushed `empirical_pairwise_strength` to
  `0.05` while keeping `epistasis_strength = 0.0`. That is a warning sign that
  the branch is absorbing residual structure through the available pairwise
  lever rather than cleanly identifying a transferable abundance-binding latent
  decomposition.

## Effect on hypothesis

- `HYP-007` is weakened for the current implementation. The experiment moved
  from the underidentified single-readout `PHOT_CHLRE` assay to a matched KRAS
  abundance-plus-binding pair, which was the sharper test requested by the
  hypothesis. The paired thermodynamic branch improved some mechanistic
  diagnostics, but it did not improve the full primary readout package of
  single-mutant ranking, double-mutant ranking, epistasis prediction, and
  reference geometry.
- The result does not directly update `HYP-001`, because it is a focused
  matched-readout diagnostic rather than a broad shared-regime multi-assay
  realism test.

## Limitations

- `RUN-028` reports deterministic branch comparisons, not a Bayesian posterior
  with synthetic-truth recovery diagnostics. The result is therefore a clean
  branch-level diagnostic, but not a full posterior identifiability analysis.
- The paired abundance objective currently worsened abundance KS in the full
  thermodynamic branch, so the experiment does not yet show that the abundance
  readout is being used as a stable transferable latent constraint.
- Because all fitted branches retained near-peak reference artifacts, the run
  cannot distinguish whether the remaining failure is mostly readout
  misspecification, inadequate latent priors, or insufficient constraints in
  the deterministic fitter.
- This is one matched KRAS readout pair, so the result weakens but does not
  refute the broader multi-latent molecular-phenotype direction.

## Downstream use

- Use `RES-012` as the project record for `EXP-009` / `RUN-028`.
- Do not cite the current paired KRAS thermodynamic implementation as support
  for `HYP-007`. It gives partial evidence that thermodynamic structure can
  improve epistasis diagnostics, but the full empirical recovery package is
  negative.
- Treat the still-unresolved near-peak reference artifact and poor
  single/double holdout ranking as central constraints on any next KRAS
  model-design iteration.
