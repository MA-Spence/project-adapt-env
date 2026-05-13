# ANA-009: RUN-030 latent-activity observed-fitness calibration review

## Purpose

Assess whether the `EXP-010` latent-activity observed-fitness generator improves
empirical recovery on `PHOT_CHLRE_Chen_2023` relative to both its within-run
controls and the earlier `RUN-026` biophysical-function baseline.

Amendment, `2026-05-13`: later code review showed that `EXP-010` was not yet a
clean explicit two-latent-trait fit. The experiment changed the observed
readout composition, but it did not configure `latent_trait_blocks`, it left
both generic and empirical epistasis acting on the same legacy functional
layer, and the Bayesian SMC stage still optimized the older bootstrap
summary-vector distance rather than the held-out validation objective adopted
later in `EXP-012`.

## Linked experiments/runs

- Experiments: `EXP-010`
- Runs: `RUN-030`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: `python3`

## Inputs

- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/summary.json`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/selected_panel.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/branch_validations.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/posterior_particles.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/synthetic_truth_recovery.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/summary.json`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/branch_validations.csv`
- `results/RES-008_run-026-biophysical-function-readout-does-not-improve-phot-chlre-recovery/metrics.json`
- `experiments/2026-05-13_EXP-010_proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-13_EXP-010_proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/runs/RUN-030.yaml`

## Analysis performed

- Reviewed the reconciled run record and confirmed that `RUN-030` completed
  successfully on `lab-slurm` under scheduler job `135`.
- Rechecked the fixed `PHOT_CHLRE_Chen_2023` panel selection and the assay-level
  `mavenn` diagnostics to determine whether the measurement layer changed
  materially relative to `RUN-026`.
- Verified from `EXP-010` `config.yaml` that the new Bayesian path calibrated
  `synthetic_readout_mode: fitness` while defining public scalar fitness through
  `observed_fitness_combine_mode: product` with
  `observed_fitness_terms = [stability_gate, function_capacity]`,
  `generic_epistasis_target: function`, and
  `empirical_pairwise_target: function`.
- Compared the two deterministic activity-readout branches and the two Bayesian
  summaries on single-mutant holdout, double-mutant holdout, epistasis
  prediction, KS, and reference-to-peak metrics.
- Compared `RUN-030` against `RUN-026` on the same assay to determine whether
  promoting activity to the public observed-fitness generator resolved the
  earlier `RES-008` failure mode.
- Rechecked the posterior parameter summary and synthetic-truth recovery outputs
  to determine whether any negative empirical result reflected a broken inverse
  problem or only remaining model mismatch.
- Reinterpreted the run in light of the later latent-trait code review,
  including whether `RUN-030` actually isolated an explicit readout trait,
  whether its epistasis knobs were competing on the same internal quantity, and
  whether the Bayesian target aligned with the held-out predictive metrics used
  to judge success.

## Outputs

- Figures: none
- Tables: `tables/run-030_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-13_ANA-009_run-030-latent-activity-observed-fitness-calibration-review/tables/run-030_key_metrics.md`

## Main observations

- `RUN-030` completed successfully and produced the declared durable outputs
  under
  `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030`,
  including the branch comparison table, posterior particles, posterior
  parameter summary, and synthetic-truth recovery diagnostics.
- The empirical target remained the same large functional ProteinGym assay used
  in `RUN-026`, `PHOT_CHLRE_Chen_2023`, with sequence length `118`,
  `167,529` measured variants total, `2,122` single mutants, and `165,407`
  multiple mutants. This is still a strong epistasis-relevant assay for a
  single-system diagnostic.
- The measurement layer was slightly better than in `RUN-026`, not worse. The
  assay-specific `mavenn` model reached test Spearman `0.689` and test NRMSE
  `0.921`, versus `0.680` and `0.935` in `RUN-026`. The new result therefore
  cannot be dismissed as an artifact of a degraded observation model.
- The deterministic controls in `RUN-030` were effectively identical to the
  deterministic biophysical-function controls from `RUN-026`. The strongest
  deterministic branch,
  `predictive_richpair_shared_activity_readout`, matched the earlier
  `predictive_richpair_shared_biophysical_function_readout` almost exactly on
  single-mutant holdout Spearman (`0.349`), double-mutant holdout Spearman
  (`0.551`), epistasis-prediction Spearman (`0.095`), functional KS (`0.356`),
  and reference location (`fraction_of_peak 1.000`, distance `0`). This shows
  that the generalized observed-fitness composition path largely reproduces the
  earlier deterministic activity-readout semantics on this assay.
- The later technical review explains why `RUN-030` should not be treated as a
  clean test of explicit latent-trait fitting. `EXP-010` defined public scalar
  fitness as `stability_gate * function_capacity`, but it did not configure any
  explicit named readout trait block, so the fitted `functional_sigma_base` and
  `n_functional_dims` parameters still referred to the same legacy functional
  layer used in earlier runs.
- The same review also identified a structural confound in the fitted
  epistasis controls. `generic_epistasis_target: function` and
  `empirical_pairwise_target: function` both pushed signal into the same
  underlying functional quantity, so `epistasis_strength` and
  `empirical_pairwise_strength` were partially competing levers rather than a
  clean separation between generic ruggedness and assay-matched pairwise
  effects.
- The Bayesian objective was also mismatched to the success criteria later used
  to interpret the run. `RUN-030` still used the earlier bootstrap
  summary-vector distance, whereas the project judged success mainly by
  held-out single-mutant ranking, double-mutant ranking, epistasis prediction,
  KS, and reference geometry. That objective mismatch is consistent with the
  observed pattern: strong synthetic-truth recovery on the fitted target,
  moderate KS improvement, but worse double-mutant and epistasis-prediction
  performance on the metrics that mattered scientifically.
- Within `RUN-030`, no branch achieved a clean empirical win. The strongest
  deterministic branch kept the best single-mutant holdout Spearman (`0.349`)
  and best double-mutant holdout Spearman (`0.551`) among the non-Bayesian
  fits, but it collapsed the fitted reference exactly onto the peak and had only
  weak epistasis-prediction Spearman (`0.095`).
- The best Bayesian fit under the new public-activity semantics improved one
  coarse distributional metric, functional KS (`0.253`), relative to the best
  deterministic control (`0.356`), and retained nonzero structured effects
  (`epistasis_strength 0.0424`, `empirical_pairwise_strength 0.0296`). But it
  did not improve the broader predictive package: single-mutant holdout Spearman
  fell from `0.349` to `0.294`, double-mutant holdout Spearman fell from
  `0.551` to `0.363`, epistasis-prediction Spearman fell from `0.095` to
  `-0.026`, and the fitted reference remained extremely near the peak
  (`fraction_of_peak 0.9995`, distance `5`).
- The key historical comparison is also negative. Relative to the earlier
  `RUN-026` Bayesian best fit on the same assay, the `RUN-030` Bayesian best
  fit kept single-mutant holdout essentially unchanged (`0.293` to `0.294`) and
  improved functional KS (`0.367` to `0.253`), but double-mutant holdout
  Spearman fell sharply (`0.562` to `0.363`), epistasis-prediction Spearman
  turned negative (`0.109` to `-0.026`), and the fitted reference moved from a
  non-pathological location (`fraction_of_peak 0.083`, distance `57`) back to a
  near-peak state (`0.9995`, distance `5`).
- The posterior mean was weaker than the best particle and did not define a
  robust empirical solution. It reached only single-mutant holdout Spearman
  `0.094`, double-mutant holdout Spearman `-0.052`, epistasis-prediction
  Spearman `-0.200`, and functional KS `0.612`, while the reference still
  remained very near the peak (`fraction_of_peak 0.959`, distance `4`).
- The inverse problem itself remained well behaved. The row-level
  `synthetic_truth_recovery.csv` shows that both preregistered truths fell
  within the posterior q90 interval for all `10/10` fitted parameters. This is
  a negative empirical result, not a trivial failure of the SMC machinery on
  matched synthetic data.
- The scientific implication is therefore narrower than the original read would
  suggest. `RUN-030` is a negative result for one coupled implementation:
  public-activity readout layered onto the legacy functional block, with
  competing epistasis levers and a summary-vector SMC target. It lowers KS, but
  it worsens the epistasis-relevant predictive metrics and reintroduces the
  near-peak reference pathology that the earlier Bayesian `RUN-026` fit had
  largely avoided.

## Result records created

- `RES-009`

## Hypothesis updates

- `HYP-007` is weakened at the implementation level, not cleanly falsified. The
  more direct public-activity formulation in `EXP-010` did not deliver the
  predicted improvement over either the within-run controls or the earlier
  `RUN-026` activity-readout baseline on `PHOT_CHLRE_Chen_2023`, but the run
  did not yet isolate a clean explicit two-latent-trait, validation-objective
  fit.
- `HYP-001` is not directly updated by this result. `RUN-030` is still a
  single-assay activity diagnostic rather than a shared-regime multi-assay test,
  so it does not establish or refute the broader benchmark-landscape claim.
