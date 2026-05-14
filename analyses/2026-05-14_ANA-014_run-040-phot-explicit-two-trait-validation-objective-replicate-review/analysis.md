# ANA-014: RUN-040 PHOT explicit two-trait validation-objective replicate review

## Purpose

Assess whether the `EXP-014` explicit two-trait PHOT validation-objective
replicate reproduces the `RUN-036` failure pattern and record the implication
for `HYP-007`.

## Linked experiments/runs

- Experiments: `EXP-014`
- Runs: `RUN-040`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: `python3`

## Inputs

- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/summary.json`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/selected_panel.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/branch_validations.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/posterior_particles.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/posterior_rounds.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/posterior_parameter_summary.csv`
- `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040/target_features.csv`
- `experiments/2026-05-14_EXP-014_proteingym-phot-explicit-two-trait-validation-objective-replicate/config.yaml`
- `experiments/2026-05-14_EXP-014_proteingym-phot-explicit-two-trait-validation-objective-replicate/runs/RUN-040.yaml`
- `results/RES-010_run-036-explicit-two-latent-trait-validation-objective-does-not-rescue-phot-chlre-recovery/metrics.json`
- `results/RES-013_run-039-stability-readout-validation-objective-improves-phot-recovery/metrics.json`

## Analysis performed

- Confirmed from the run metadata that `RUN-040` completed successfully on
  `lab-slurm` as scheduler job `173` with exit code `0`.
- Confirmed that `RUN-040` is the collected run for `EXP-014`, the direct
  explicit two-trait validation-objective replicate of the `RUN-036` PHOT
  experiment.
- Verified from `EXP-014` `config.yaml` that the run again used the intended
  explicit two-trait setup:
  - built-in stability plus one named latent block, `readout`
  - public scalar fitness defined as `stability_gate * trait:readout:capacity`
  - `synthetic_readout_mode: fitness`
  - `empirical_pairwise_target: trait:readout`
  - fixed `n_functional_dims = 1`
  - fixed `peak_distance_from_consensus = 2`
  - `epistasis_strength = 0`
  - SMC `distance_mode: validation_objective`
- Reviewed the MAVE-NN assay diagnostics to check whether the observation
  layer was usable and comparable to the previous PHOT runs.
- Compared deterministic two-trait controls, the best SMC particle, and the
  posterior-mean particle on single-mutant holdout, double-mutant holdout,
  epistasis prediction, KS, and reference-to-peak geometry.
- Compared the `RUN-040` best particle against the earlier `RUN-036` explicit
  two-trait validation-objective result and the `RUN-039` stability-targeted
  validation-objective control.
- Reviewed posterior parameter summaries and SMC round progression for signs
  that failure was due to stochastic search rather than model-family mismatch.
- Recorded the absence of matched synthetic-truth recovery because
  `validation_objective` mode currently disables that scaffold.

## Outputs

- Figures: none
- Tables: `tables/run-040_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-14_ANA-014_run-040-phot-explicit-two-trait-validation-objective-replicate-review/tables/run-040_key_metrics.md`

## Main observations

- `RUN-040` completed successfully and produced the declared durable outputs
  under
  `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040`.
- The empirical target remained `PHOT_CHLRE_Chen_2023`, an activity/FACS assay
  with sequence length `118`, `167,529` measured variants, `2,122` single
  mutants, `165,407` multiple mutants, and `2,298` calibration variants.
- The assay observation layer was moderate but usable. The MAVE-NN model
  reached test Spearman `0.677` and test NRMSE `0.878`. This is not a collapsed
  readout model and is close enough to the previous PHOT diagnostics that it
  does not explain the replicate failure by itself.
- The deterministic explicit two-trait branches were weak. The stronger
  deterministic branch for predictive metrics reached single-mutant holdout
  Spearman `0.191`, double-mutant holdout Spearman `0.216`, epistasis-prediction
  Spearman `0.240`, and functional KS `0.306`, while keeping the reference
  almost at the peak (`fraction_of_peak 0.99994`, distance `2`).
- The best SMC particle improved ranking only slightly over that deterministic
  branch: single-mutant holdout Spearman `0.198`, double-mutant holdout
  Spearman `0.252`, and epistasis-prediction Spearman `0.264`. Functional KS
  worsened to `0.498`, and the reference still remained near the peak
  (`fraction_of_peak 0.986`, distance `6`).
- The posterior mean was close to the best particle on the core predictive
  metrics: single-mutant holdout Spearman `0.198`, double-mutant holdout
  Spearman `0.251`, epistasis-prediction Spearman `0.255`, functional KS
  `0.483`, and reference distance `5`.
- The replicate comparison is the key result. Relative to `RUN-036`, `RUN-040`
  changed almost nothing: single-mutant holdout Spearman delta `0.000`,
  double-mutant holdout Spearman delta `+0.001`, epistasis-prediction Spearman
  delta `+0.006`, functional KS delta `+0.002`, reference fraction delta near
  `0.000`, and reference-distance delta `0`.
- Relative to the `RUN-039` stability-targeted validation-objective control,
  the explicit two-trait replicate is much weaker: single-mutant holdout
  Spearman delta `-0.308`, double-mutant holdout Spearman delta `-0.531`,
  epistasis-prediction Spearman delta `-0.370`, functional KS delta `+0.322`,
  reference fraction delta `+0.663`, and reference-distance delta `-39`.
- The posterior again pushed the remaining empirical pairwise term to the top
  of its prior range. The posterior mean for `empirical_pairwise_strength` was
  `0.0777`, with q05 `0.0768` and q95 `0.0797` against the configured upper
  bound `0.08`; the best particle used `0.0796`.
- SMC did not reveal a late better basin. The best distance moved only from
  `23.299` in round 0 to `23.228` in round 3, and the posterior mean remained
  close to the best particle.
- The result is therefore difficult to explain as stochastic SMC noise. It is
  more consistent with the current explicit two-trait PHOT latent mapping being
  structurally mismatched to this scalar activity assay, at least as currently
  parameterized.
- The main limitation is still causal resolution. Because matched synthetic
  recovery is disabled in `validation_objective` mode, `RUN-040` strengthens an
  empirical model-comparison claim but does not by itself prove which piece of
  the causal map is wrong.

## Result records created

- `RES-014`

## Hypothesis updates

- `HYP-007` is weakened. The direct replicate of the clean explicit
  two-trait validation-objective PHOT experiment reproduced the earlier failure
  almost exactly, while the simpler stability-targeted control from `RUN-039`
  fit the same scalar assay much better.
- `HYP-001` is not directly updated because this remains a one-assay PHOT
  diagnostic rather than a shared-regime, multi-assay realism test.
