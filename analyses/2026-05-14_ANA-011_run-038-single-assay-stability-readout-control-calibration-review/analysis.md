# ANA-011: RUN-038 single-assay stability-readout control calibration review

## Purpose

Assess whether the simpler `EXP-011` stability-targeted single-latent control
recovers `PHOT_CHLRE_Chen_2023` better than the recent PHOT
activity-oriented formulations, and determine what that control implies for
`HYP-007`.

## Linked experiments/runs

- Experiments: `EXP-011`
- Runs: `RUN-038`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: `python3`

## Inputs

- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/summary.json`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/selected_panel.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/branch_validations.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/posterior_particles.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/synthetic_truth_recovery.csv`
- `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038/target_features.csv`
- `results/RES-007_run-024-stability-targeted-single-assay-calibration-materially-improves-sptn1-recovery/metrics.json`
- `results/RES-008_run-026-biophysical-function-readout-does-not-improve-phot-chlre-recovery/metrics.json`
- `results/RES-009_run-030-latent-activity-observed-fitness-readout-does-not-rescue-phot-chlre-recovery/metrics.json`
- `results/RES-010_run-036-explicit-two-latent-trait-validation-objective-does-not-rescue-phot-chlre-recovery/metrics.json`
- `experiments/2026-05-13_EXP-011_proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-13_EXP-011_proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/runs/RUN-038.yaml`

## Analysis performed

- Reviewed the reconciled run record and confirmed that `RUN-038` completed
  successfully on `lab-slurm` under scheduler job `163`.
- Rechecked the fixed `PHOT_CHLRE_Chen_2023` panel selection and the assay-level
  `mavenn` diagnostics to determine whether the observation layer changed
  materially relative to the earlier PHOT runs.
- Verified from `EXP-011` `config.yaml` that this run reuses the successful
  `RUN-024` stability-targeted semantics on a different assay:
  - `synthetic_readout_mode: stability_margin`
  - `empirical_pairwise_target: stability`
  - bootstrap summary-vector SMC distance rather than the later
    `validation_objective` mode
- Compared the raw-readout controls, the stability-readout deterministic
  controls, and the Bayesian summaries on single-mutant holdout, double-mutant
  holdout, epistasis prediction, KS, and reference-to-peak geometry.
- Compared `RUN-038` against the earlier PHOT Bayesian baselines from
  `RUN-026`, `RUN-030`, and `RUN-036` to determine whether the simpler control
  improves the overall activity-assay fit.
- Compared `RUN-038` against the earlier successful stability-assay result from
  `RUN-024` to determine whether the same stability-targeted formulation
  performs as cleanly on an activity readout as it did on the unfolding assay.
- Reviewed the posterior parameter summary, SMC round progression, and
  synthetic-truth recovery outputs to distinguish model mismatch from inverse-
  problem failure.

## Outputs

- Figures: none
- Tables: `tables/run-038_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-14_ANA-011_run-038-single-assay-stability-readout-control-calibration-review/tables/run-038_key_metrics.md`

## Main observations

- `RUN-038` completed successfully and produced the declared durable outputs
  under
  `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038`,
  including branch comparisons, posterior particles, posterior summaries,
  target features, and synthetic-truth recovery diagnostics.
- The empirical target remained the same large functional ProteinGym assay used
  in the earlier PHOT runs, `PHOT_CHLRE_Chen_2023`, with sequence length `118`,
  `167,529` measured variants total, `2,122` single mutants, and `165,407`
  multiple mutants. This remains a strong one-assay epistasis test.
- The measurement layer remained moderate rather than collapsing. The
  assay-specific `mavenn` model reached test Spearman `0.686` and test NRMSE
  `0.916`, which is effectively in line with the earlier PHOT runs
  (`0.680-0.689` Spearman). The control result is therefore not explained by a
  different observation-model regime.
- `EXP-011` really did test the intended simple control. Unlike the later clean
  two-trait run, this experiment deliberately reuses the same stability-targeted
  calibration semantics that worked well on the `SPTN1` stability assay:
  `synthetic_readout_mode: stability_margin` with empirical pairwise structure
  routed into latent stability.
- The target features confirm that this run used the older bootstrap
  summary-vector distance rather than the later `validation_objective` mode.
  The Bayesian fit was therefore optimized against nine assay-level summary
  features, including DFE fractions, skewness, conservation correlation, and
  epistasis variance.
- The stability-readout control mattered immediately within the run. Both
  deterministic stability-readout branches strongly outperformed the raw-readout
  controls on single-mutant holdout ranking and functional KS, and the richer
  predictive stability branch also raised double-mutant holdout Spearman to
  `0.701`. This shows that the negative PHOT outcome was not robust to reverting
  to the simpler `RUN-024`-style semantics.
- Those deterministic gains were still not a clean rescue. The deterministic
  stability-readout branches kept the fitted reference essentially at the peak
  (`fraction_of_peak 1.000` and `0.987`, distance `0` and `3`) and had weak or
  negative epistasis-prediction Spearman (`-0.250` and `0.008`).
- The best Bayesian control fit produced the strongest PHOT ranking-and-KS
  package seen so far. It reached single-mutant holdout Spearman `0.505`,
  double-mutant holdout Spearman `0.721`, and functional KS `0.166`, which are
  all better than the earlier PHOT Bayesian results from `RUN-026`, `RUN-030`,
  and `RUN-036`.
- That best Bayesian PHOT fit still did not fully solve the assay. Its
  epistasis-prediction Spearman was only `0.066`, below both `RUN-026`
  (`0.109`) and `RUN-036` (`0.259`), and the fitted reference remained too near
  the peak (`fraction_of_peak 0.867`, distance `4`) to count as a fully clean
  mechanistic recovery.
- The posterior mean did not reinforce one clean all-metric solution. It kept
  nearly the same ranking quality as the best particle and improved reference
  geometry (`fraction_of_peak 0.337`, distance `21`), but epistasis-prediction
  Spearman turned negative (`-0.156`). This is more consistent with metric
  tradeoffs than with a single robustly correct posterior basin.
- Relative to `RUN-024`, the same stability-targeted formulation does not
  perform as cleanly on the activity readout as it did on the unfolding assay.
  `RUN-038` is slightly worse on single-mutant holdout (`0.505` versus `0.556`)
  and dramatically worse on epistasis prediction (`0.066` versus `0.540`) and
  reference geometry (`fraction_of_peak 0.867`, distance `4` versus `0.039`,
  distance `50`), even though it is better on double-mutant holdout (`0.721`
  versus `0.483`) and functional KS (`0.166` versus `0.233`).
- The inverse problem remained well behaved. Both preregistered synthetic truths
  fell within the posterior q90 interval for `10/10` fitted parameters, so the
  negative parts of the PHOT result are not explained by a trivial SMC failure
  on matched synthetic data.
- The scientific implication is therefore mixed but informative. The simple
  stability-targeted control clearly outperforms the recent PHOT
  activity-oriented implementations on ranking and KS, which weakens the current
  argument that the richer PHOT formulations are adding useful mechanistic
  structure on this assay. But the control still leaves weak epistasis
  prediction and a near-peak geometry problem, so it does not fully rescue
  `PHOT_CHLRE` or show that a stability-only causal chain is the correct final
  explanation.

## Result records created

- `RES-011`

## Hypothesis updates

- `HYP-007` is weakened at the implementation level on `PHOT_CHLRE_Chen_2023`.
  The simpler stability-targeted control outperforms the recent PHOT
  activity-oriented Bayesian formulations on the main ranking and KS metrics, so
  the current richer activity/readout formulations are not yet providing the
  predicted empirical gain on this assay.
- `HYP-001` is not directly updated by this result. `RUN-038` remains a
  one-assay activity diagnostic rather than a shared-regime multi-assay test,
  and its evidence is still variant-level rather than the original summary-first
  shared-panel claim.
