# RES-013: RUN-039 stability-readout validation objective strongly improves PHOT recovery relative to two-trait variants

## Summary

On PHOT_CHLRE_Chen_2023, the stability-targeted validation-objective Bayesian fit improves the full predictive package relative to the explicit two-trait validation-objective run: single-mutant holdout ranking remains strong, double-mutant holdout ranking and epistasis prediction improve sharply, KS stays low, and the fitted reference moves away from the peak. This indicates that the current PHOT multi-latent/readout implementations are still structurally mismatched, while a stability-targeted causal mapping fits the scalar activity assay much better.

## Generated from

- Analyses: `ANA-013`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-14_ANA-013_run-039-phot-stability-readout-validation-objective-control-review/tables/run-039_key_metrics.md`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/summary.json`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/selected_panel.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/branch_validations.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/posterior_parameter_summary.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/posterior_rounds.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/target_features.csv`
- `results/RES-010_run-036-explicit-two-latent-trait-validation-objective-does-not-rescue-phot-chlre-recovery/metrics.json`
- `results/RES-011_run-038-stability-readout-control-improves-phot-chlre-recovery-but-does-not-fully-rescue-the-assay/metrics.json`
- `experiments/2026-05-14_EXP-013_proteingym-phot-stability-readout-validation-objective-control/config.yaml`
- `experiments/2026-05-14_EXP-013_proteingym-phot-stability-readout-validation-objective-control/runs/RUN-039.yaml`

## Interpretation

- `RUN-039` completed the intended `EXP-013` control on
  `PHOT_CHLRE_Chen_2023`: the simpler stability-targeted calibration was fit
  under the same held-out `validation_objective` family that was used by the
  explicit two-latent `RUN-036` experiment.
- The measurement layer was comparable to the prior PHOT runs. The MAVE-NN
  observation model reached test Spearman `0.689` and test NRMSE `0.966`.
- The deterministic stability-readout branches again improved single-mutant
  ranking and KS relative to raw-readout branches, but they did not solve the
  full problem. The stronger deterministic stability branch had single-mutant
  holdout Spearman `0.505`, double-mutant holdout Spearman `0.703`,
  epistasis-prediction Spearman `0.032`, KS `0.190`, and reference fraction of
  peak `0.987`.
- The Bayesian validation-objective stability-targeted fit was much stronger.
  The best particle reached single-mutant holdout Spearman `0.506`,
  double-mutant holdout Spearman `0.783`, epistasis-prediction Spearman
  `0.634`, functional KS `0.176`, and a non-pathological reference geometry
  (`fraction_of_peak 0.323`, distance `45`).
- Relative to the explicit two-latent validation-objective fit from `RUN-036`,
  the `RUN-039` Bayesian best particle improved single-mutant ranking
  (`0.198 -> 0.506`), double-mutant ranking (`0.251 -> 0.783`), epistasis
  prediction (`0.259 -> 0.634`), KS (`0.496 -> 0.176`), and reference geometry
  (`fraction_of_peak 0.986 -> 0.323`, distance `6 -> 45`).
- Relative to the earlier summary-vector stability control from `RUN-038`,
  `RUN-039` kept the single-mutant advantage, improved double-mutant holdout
  ranking (`0.721 -> 0.783`), greatly improved epistasis prediction
  (`0.066 -> 0.634`), and moved the reference much farther from the peak, with
  only a small KS regression (`0.166 -> 0.176`).
- The posterior mean did not retain the epistasis-prediction gain
  (`-0.147`), so the posterior distribution still contains tradeoffs. The best
  particle, however, demonstrates that the configured stability-targeted model
  family can jointly satisfy ranking, KS, epistasis prediction, and reference
  geometry much better than the current PHOT two-trait implementation.
- Matched synthetic-truth recovery is absent because `validation_objective`
  mode currently disables that scaffold. The result is therefore strongest as
  an empirical model-comparison result, not as a full identifiability proof.

## Effect on hypothesis

- `HYP-007` is weakened at the current PHOT implementation level. The explicit
  two-latent model no longer has the simple objective-mismatch defense: a
  stability-targeted mapping fit with the same validation objective performs
  substantially better on the same scalar activity assay.
- The result does not falsify the broader biological idea that multiple latent
  molecular traits can matter. It says that the currently implemented PHOT
  multi-latent/readout formulation is not yet the right causal map for this
  assay, while a stability-targeted latent mapping is empirically much closer.
- `HYP-001` is not directly updated because this is still a one-assay
  diagnostic rather than a shared-regime, multi-assay realism test.

## Limitations

- The result is based on one assay, `PHOT_CHLRE_Chen_2023`.
- The observation model is only moderate by recent project standards, so the
  result should be interpreted through large metric differences rather than
  fine distinctions.
- The SMC posterior mean does not preserve the best particle's epistasis gain,
  indicating that posterior mass still spans metric tradeoffs.
- `validation_objective` mode lacks matched synthetic-truth recovery, so the
  run does not isolate inverse-problem identifiability from empirical
  structural fit.
- A strong stability-targeted fit to a scalar activity assay does not prove
  that the physical assay is purely a stability measurement.

## Downstream use

- Use `RES-013` as the control showing that the poor `RUN-036` explicit
  two-latent result is not explained by the held-out validation objective
  itself.
- Use `RES-013` together with `RES-010` and `RES-011` to narrow the PHOT
  structural-mismatch claim: the current activity/two-trait readout model is
  weaker than a stability-targeted causal mapping on this scalar assay.
