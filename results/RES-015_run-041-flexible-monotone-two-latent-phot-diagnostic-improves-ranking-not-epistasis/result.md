# RES-015: RUN-041 flexible monotone two-latent PHOT diagnostic improves ranking but not epistasis

## Summary

On `PHOT_CHLRE_Chen_2023`, the `EXP-015` flexible monotone
stability/activity readout diagnostic improved held-out single- and
double-mutant ranking relative to the explicit product/gate two-trait failure,
and double-mutant ranking slightly exceeded the stability-targeted
validation-objective control. It did not rescue epistasis prediction:
epistasis Spearman was near zero and KS was poor. The result suggests the
hard-coded product/gate readout contributed to ranking failure, but a flexible
two-latent monotone surface is still not a sufficient `HYP-007` rescue.

## Generated from

- Analyses: `ANA-015`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-14_ANA-015_run-041-phot-flexible-monotone-two-latent-readout-diagnostic-review/tables/run-041_key_metrics.md`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/summary.json`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/selected_panel.csv`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/model_validations.csv`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/progress.json`
- `experiments/2026-05-14_EXP-015_proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/config.yaml`
- `experiments/2026-05-14_EXP-015_proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/runs/RUN-041.yaml`
- `scripts/proteingym_phot_structural_mismatch_diagnostic.py`
- `results/RES-013_run-039-stability-readout-validation-objective-improves-phot-recovery/metrics.json`
- `results/RES-014_run-040-explicit-two-trait-replicate-reproduces-phot-validation-objective-failure/metrics.json`

## Interpretation

- `RUN-041` completed the intended `EXP-015` flexible monotone readout
  diagnostic on the same `PHOT_CHLRE_Chen_2023` scalar activity assay.
- This is not a full AdaptEnv posterior fit. It is a lightweight structural
  diagnostic that fits an alignment-frequency stability proxy, an additive
  ridge activity score, and a constrained monotone two-dimensional surface over
  those coordinates.
- The diagnostic therefore tests whether the hard-coded product/gate collapse
  was a ranking bottleneck, but it cannot assess reference-to-peak geometry or
  full functional KS.
- The observation layer was moderate. The assay-specific MAVE-NN model reached
  test Spearman `0.657` and test NRMSE `0.947`, slightly weaker than the
  immediate PHOT SMC runs but still usable for a coarse diagnostic.
- The flexible surface reached single-mutant holdout Spearman `0.473` and
  double-mutant holdout Spearman `0.807`. Relative to `RUN-040`, this improves
  single-mutant ranking by `0.275` and double-mutant ranking by `0.555`.
- Relative to the stability-targeted `RUN-039` best particle, `RUN-041` is
  slightly worse on single-mutant ranking (`-0.033`) and slightly better on
  double-mutant ranking (`+0.024`).
- The epistasis result is negative. Epistasis-prediction Spearman was only
  `0.035`, with KS `0.814`. Relative to `RUN-040`, epistasis Spearman worsened
  by `0.230`; relative to `RUN-039`, it worsened by `0.600`.
- This separates two failure modes. The hard-coded product/gate readout appears
  to be part of the held-out ranking mismatch, but relaxing that readout into a
  flexible monotone two-latent surface does not recover the epistatic structure
  of the PHOT assay.
- The strongest defensible interpretation is mixed: `RUN-041` supports readout
  misspecification as a real component of the structural mismatch, but it
  weakens the stronger claim that the current two-latent stability/activity
  direction is close to a complete empirical recovery model.

## Effect on hypothesis

- `HYP-007` is weakened as a full empirical recovery claim. The flexible
  two-latent diagnostic improved ranking but failed one of the primary readouts,
  epistasis prediction.
- The result does not refute every multiple-latent molecular phenotype model.
  It narrows the problem: readout flexibility alone is insufficient, and PHOT
  epistasis still appears mismatched or underidentified from this scalar assay.
- `HYP-001` is not directly updated because this is still a one-assay PHOT
  diagnostic, not a multi-assay realism test.

## Limitations

- The run uses one assay, `PHOT_CHLRE_Chen_2023`.
- `RUN-041` is a diagnostic model, not a full AdaptEnv generator or Bayesian
  posterior, so it cannot be compared on reference-to-peak geometry,
  posterior identifiability, or synthetic-truth recovery.
- The double-mutant holdout set is small (`43` variants), so the strong
  double-mutant ranking should be interpreted as suggestive rather than
  definitive.
- The stability coordinate is a conservation-derived proxy, not an independent
  measured stability trait.
- The activity coordinate is fit from the same scalar assay, so the latent
  decomposition is predictive rather than causally identified.

## Downstream use

- Use `RES-015` to separate readout-collapse mismatch from epistasis mismatch:
  flexible monotone two-latent readout helps held-out ranking but does not
  recover epistasis.
- Use `RES-015` with `RES-013` and `RES-014` when interpreting later PHOT
  structural-mismatch diagnostics, especially `EXP-016` and `EXP-017`.
