# ANA-010: RUN-036 explicit two-latent-trait validation-objective calibration review

## Purpose

Assess whether the clean `EXP-012` explicit two-latent-trait formulation
improves empirical recovery on `PHOT_CHLRE_Chen_2023` relative to the earlier
PHOT implementations once the latent-trait confounders are removed and the
Bayesian fit optimizes the held-out validation objective directly.

## Linked experiments/runs

- Experiments: `EXP-012`
- Runs: `RUN-036`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: `python3`

## Inputs

- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/summary.json`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/selected_panel.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/mavenn_assay_metrics.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/branch_validations.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/posterior_particles.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/posterior_rounds.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/posterior_parameter_summary.csv`
- `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036/target_features.csv`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/summary.json`
- `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026/branch_validations.csv`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/summary.json`
- `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030/branch_validations.csv`
- `results/RES-008_run-026-biophysical-function-readout-does-not-improve-phot-chlre-recovery/metrics.json`
- `results/RES-009_run-030-latent-activity-observed-fitness-readout-does-not-rescue-phot-chlre-recovery/metrics.json`
- `experiments/2026-05-13_EXP-012_proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/config.yaml`
- `experiments/2026-05-13_EXP-012_proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/runs/RUN-036.yaml`

## Analysis performed

- Reviewed the reconciled run record and confirmed that `RUN-036` completed
  successfully on `lab-slurm` under scheduler job `162`.
- Rechecked the fixed `PHOT_CHLRE_Chen_2023` panel selection and the assay-level
  `mavenn` diagnostics to determine whether the observation layer changed
  materially relative to `RUN-026` and `RUN-030`.
- Verified from `EXP-012` `config.yaml` that the run now cleanly instantiates
  two latent traits in practice:
  - built-in stability
  - one explicit named latent block, `readout`
  - public scalar fitness defined as `stability_gate * trait:readout:capacity`
  - `generic_epistasis_target: score` with `epistasis_strength = 0`
  - `empirical_pairwise_target: trait:readout`
  - fixed `n_functional_dims = 1`
  - fixed `peak_distance_from_consensus = 2`
- Rechecked `target_features.csv` and `summary.json` to confirm that the SMC
  distance was the direct held-out validation objective, with
  `objective_core`, `objective_double`, and `objective_total` as zero-centered
  target features rather than the earlier bootstrap summary-vector target.
- Compared the deterministic two-trait controls and the Bayesian posterior
  summaries on single-mutant holdout, double-mutant holdout, epistasis
  prediction, KS, and reference-to-peak geometry.
- Compared `RUN-036` against the earlier PHOT Bayesian baselines from
  `RUN-026` and `RUN-030` to determine whether the cleaner latent-trait setup
  improved the full predictive package.
- Reviewed the posterior parameter summary and SMC round progression for signs
  of weak convergence, unused fit levers, or saturation of the remaining
  empirical pairwise term.
- Recorded the absence of matched synthetic-truth recovery in this run because
  `validation_objective` mode currently disables that scaffold.

## Outputs

- Figures: none
- Tables: `tables/run-036_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-13_ANA-010_run-036-explicit-two-latent-trait-validation-objective-calibration-review/tables/run-036_key_metrics.md`

## Main observations

- `RUN-036` completed successfully and produced the declared durable outputs
  under
  `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036`,
  including branch validations, posterior particles, posterior summaries, and
  the validation-objective target features.
- The empirical target remained the same large functional ProteinGym assay used
  in `RUN-026` and `RUN-030`, `PHOT_CHLRE_Chen_2023`, with sequence length
  `118`, `167,529` measured variants total, `2,122` single mutants, and
  `165,407` multiple mutants. This remains a strong one-assay epistasis test.
- The measurement layer was usable but only moderate. The assay-specific
  `mavenn` model reached test Spearman `0.686` and test NRMSE `0.865`, which is
  slightly below `RUN-030` on Spearman but somewhat better on NRMSE. The new
  negative result is therefore not explained by a collapsed observation layer.
- Unlike `RUN-026` and `RUN-030`, this run did test the intended clean
  two-latent setup. `EXP-012` explicitly instantiated one named `readout`
  trait block, kept stability as the other latent trait, routed only empirical
  pairwise epistasis into `trait:readout`, fixed generic epistasis to zero, and
  optimized the held-out validation objective directly rather than the earlier
  bootstrap summary-vector distance.
- The target features confirm that objective change. `RUN-036` fit directly to
  zero-centered held-out criteria,
  `objective_core`, `objective_double`, and `objective_total`, rather than to a
  summary-statistics vector that was only indirectly related to the later
  success criteria.
- The deterministic controls were already weak on this assay. The stronger
  deterministic branch for the main predictive metrics,
  `baseline_shared_two_trait_readout`, reached only single-mutant holdout
  Spearman `0.191`, double-mutant holdout Spearman `0.216`,
  epistasis-prediction Spearman `0.240`, and functional KS `0.305`, while the
  fitted reference remained essentially at the peak (`fraction_of_peak 0.99994`,
  distance `2`).
- The Bayesian fit improved some within-run predictive metrics but not the
  whole package. Relative to that stronger deterministic control, the best
  Bayesian particle increased single-mutant holdout Spearman from `0.191` to
  `0.198`, double-mutant holdout Spearman from `0.216` to `0.251`, and
  epistasis-prediction Spearman from `0.240` to `0.259`, but functional KS
  worsened from `0.305` to `0.496` and the fitted reference still remained near
  the peak (`fraction_of_peak 0.986`, distance `6`).
- The historical comparison is negative on the metrics that matter most for
  recovery. Relative to the earlier `RUN-030` Bayesian best fit on the same
  assay, `RUN-036` improved epistasis-prediction Spearman from `-0.026` to
  `0.259`, but single-mutant holdout Spearman fell from `0.294` to `0.198`,
  double-mutant holdout Spearman fell from `0.363` to `0.251`, functional KS
  worsened from `0.253` to `0.496`, and the reference remained near the peak
  (`fraction_of_peak 0.986`, distance `6` versus `0.9995`, distance `5`).
- The paired comparison to `RUN-026` is also negative overall. `RUN-036`
  improved epistasis-prediction Spearman from `0.109` to `0.259`, but
  single-mutant holdout Spearman fell from `0.293` to `0.198`, double-mutant
  holdout Spearman fell from `0.562` to `0.251`, functional KS worsened from
  `0.367` to `0.496`, and the fitted reference moved back toward the peak
  (`fraction_of_peak 0.986`, distance `6` versus `0.083`, distance `57`).
- The posterior behavior suggests that the cleaned-up model used its remaining
  epistatic freedom heavily without finding a strong global solution. The
  posterior mean for `empirical_pairwise_strength` was `0.0768` with `q95`
  `0.0796`, effectively at the top of the configured prior range `0.00-0.08`,
  while the best SMC distance improved only marginally across rounds
  (`23.346 -> 23.248`).
- The posterior mean was close to the best particle on the main predictive
  metrics and remained near the peak, so the result is not hiding a clearly
  better posterior basin that the single best particle missed.
- The absence of synthetic-truth recovery is an important limitation. Because
  `validation_objective` mode currently disables the matched synthetic-truth
  scaffold, this run cannot by itself separate structural model mismatch from
  any inverse-problem difficulty introduced by the new objective.
- Even with that limitation, the scientific interpretation is narrower and
  cleaner than for `RES-008` and `RES-009`. This is the first PHOT run that
  actually removed the earlier latent-block confounders and objective mismatch.
  Its negative empirical outcome therefore weakens the explicit two-latent
  `PHOT_CHLRE` formulation more directly than the earlier amended records did.

## Result records created

- `RES-010`

## Hypothesis updates

- `HYP-007` is weakened more directly on `PHOT_CHLRE_Chen_2023`. Once the
  explicit second trait block, single remaining epistasis lever, and held-out
  validation objective were all implemented as intended, the model still did
  not recover a strong overall predictive package on this activity assay.
- `HYP-001` is not directly updated by this result. `RUN-036` remains a
  one-assay activity diagnostic rather than a shared-regime multi-assay test.
