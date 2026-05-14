# RES-014: RUN-040 explicit two-trait replicate reproduces PHOT validation-objective failure

## Summary

On `PHOT_CHLRE_Chen_2023`, the `EXP-014` explicit two-trait
validation-objective replicate reproduced the `RUN-036` failure pattern almost
exactly: single- and double-mutant holdout ranking remained weak, functional KS
worsened, the fitted reference stayed near the peak, and empirical pairwise
strength again sat at the top of its prior. This strengthens the evidence that
the current PHOT explicit two-trait implementation is structurally mismatched
rather than just noisy.

## Generated from

- Analyses: `ANA-014`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-007`
- Refutes: None
- Motivates: None

## Evidence

- `analyses/2026-05-14_ANA-014_run-040-phot-explicit-two-trait-validation-objective-replicate-review/tables/run-040_key_metrics.md`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/summary.json`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/selected_panel.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/branch_validations.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/posterior_parameter_summary.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/posterior_rounds.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/target_features.csv`
- `results/RES-010_run-036-explicit-two-latent-trait-validation-objective-does-not-rescue-phot-chlre-recovery/metrics.json`
- `results/RES-013_run-039-stability-readout-validation-objective-improves-phot-recovery/metrics.json`
- `experiments/2026-05-14_EXP-014_proteingym-phot-explicit-two-trait-validation-objective-replicate/config.yaml`
- `experiments/2026-05-14_EXP-014_proteingym-phot-explicit-two-trait-validation-objective-replicate/runs/RUN-040.yaml`

## Interpretation

- `RUN-040` completed the intended direct replicate of the explicit two-trait
  PHOT validation-objective experiment on `PHOT_CHLRE_Chen_2023`.
- The run again tested the intended current two-trait mapping: built-in
  stability plus an explicit `readout` latent trait, observed scalar fitness as
  `stability_gate * trait:readout:capacity`, pairwise empirical terms targeted
  to `trait:readout`, fixed one-dimensional readout, fixed reference peak
  geometry, and a held-out validation objective.
- The measurement layer was moderate but usable. MAVE-NN reached test Spearman
  `0.677` and test NRMSE `0.878`, which is not strong enough to remove
  measurement noise as a limitation but is not a failed observation model.
- The best SMC particle reached single-mutant holdout Spearman `0.198`,
  double-mutant holdout Spearman `0.252`, epistasis-prediction Spearman
  `0.264`, functional KS `0.498`, and reference geometry still near the peak
  (`fraction_of_peak 0.986`, distance `6`).
- The replicate matched `RUN-036` almost exactly. The best-particle deltas were
  `0.000` for single-mutant holdout Spearman, `+0.001` for double-mutant
  holdout Spearman, `+0.006` for epistasis-prediction Spearman, `+0.002` for
  functional KS, approximately `0.000` for reference fraction of peak, and `0`
  for reference distance.
- This makes stochastic SMC noise an unlikely explanation for the `RUN-036`
  failure. The same model family, same assay, and same validation-objective
  family reproduced the same weak predictive package.
- The paired comparison to the `RUN-039` stability-targeted validation-objective
  control is strongly negative for the explicit two-trait mapping. Relative to
  `RUN-039`, the `RUN-040` best particle was lower on single-mutant ranking
  (`-0.308`), double-mutant ranking (`-0.531`), and epistasis prediction
  (`-0.370`), while KS was worse (`+0.322`) and reference geometry was much
  more peak-proximal (`fraction_of_peak +0.663`, distance `-39`).
- The posterior again pushed `empirical_pairwise_strength` to the top of its
  prior. The posterior mean was `0.0777` with q05 `0.0768` and q95 `0.0797`
  against an upper bound of `0.08`, and the best particle used `0.0796`.
  Therefore the model is not failing because this configured pairwise freedom
  was unused.
- The SMC rounds do not show a credible late rescue. Best distance only moved
  from `23.299` in round 0 to `23.228` in round 3, and posterior-mean metrics
  remained close to the best-particle metrics.
- The most defensible conclusion is narrow: the current explicit stability plus
  readout latent mapping is a poor empirical model for this scalar PHOT assay
  as currently implemented. This does not prove that multiple latent molecular
  traits are biologically irrelevant.

## Effect on hypothesis

- `HYP-007` is weakened. A direct replicate of the current explicit two-trait
  PHOT validation-objective model reproduced the earlier failure while the
  simpler stability-targeted validation-objective control fit the same assay
  much better.
- `HYP-001` is not directly updated. This remains a one-assay PHOT diagnostic,
  not a shared-regime multi-assay realism test.

## Limitations

- The result is based on one assay, `PHOT_CHLRE_Chen_2023`.
- The assay observation model is moderate rather than high-performing, so
  small metric differences should not be overinterpreted.
- `validation_objective` mode lacks matched synthetic-truth recovery, so this
  result does not isolate inverse-problem identifiability from empirical model
  mismatch.
- The result tests the current implemented explicit two-trait mapping and
  priors; it does not falsify all possible stability/activity causal maps or
  all non-linear combinations of additive latent traits.
- The exported SMC branch labels retain the historical `smc_abc_*_raw` suffix,
  so interpretation depends on the experiment config and not the branch label
  alone.

## Downstream use

- Use `RES-014` with `RES-010` as replicate evidence that the explicit
  two-trait PHOT validation-objective failure is reproducible.
- Use `RES-014` with `RES-013` to support the narrower structural-mismatch
  claim: on this scalar activity assay, the current PHOT explicit two-trait
  readout mapping is much weaker than the stability-targeted mapping under the
  same validation-objective family.
