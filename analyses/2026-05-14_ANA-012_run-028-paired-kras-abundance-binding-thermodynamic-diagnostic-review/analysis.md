# ANA-012: RUN-028 paired KRAS abundance-binding thermodynamic diagnostic review

## Purpose

Assess whether the EXP-009 paired KRAS abundance-plus-binding thermodynamic readout improves empirical recovery relative to binding-only and paired raw-readout controls, and determine the implication for HYP-007.

## Linked experiments/runs

- Experiments: EXP-009
- Runs: RUN-028

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `tables/run-028_key_metrics.md`

## Inputs

- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/summary.json`
- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/selected_panel.csv`
- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/mavenn_assay_metrics.csv`
- `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028/branch_validations.csv`
- `experiments/2026-05-13_EXP-009_proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/config.yaml`
- `experiments/2026-05-13_EXP-009_proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/runs/RUN-028.yaml`
- `results/RES-007_run-024-stability-targeted-single-assay-calibration-materially-improves-sptn1-recovery/result.md`
- `results/RES-008_run-026-biophysical-function-readout-does-not-improve-phot-chlre-recovery/result.md`

## Analysis performed

- Searched `analyses/`, `results/`, `registry/`, `PROJECT_STATE.md`, and the
  experiment records for existing `RUN-028` / `EXP-009` analysis or result
  records.
- Confirmed that `RUN-028` is linked to `EXP-009` and marked completed with
  scheduler job `96` and exit code `0`.
- Rechecked the `EXP-009` scientific question and prediction: the paired KRAS
  abundance-plus-binding thermodynamic readout should improve binding holdout
  ranking, double-mutant recovery, epistasis prediction, and reference-to-peak
  behavior relative to binding-only and paired raw-readout controls.
- Reviewed the selected KRAS panel, assay-specific MAVE-NN observation-layer
  diagnostics, and branch validation table.
- Compared four deterministic branches:
  - `binding_only_raw`
  - `binding_only_biophysical_binding`
  - `paired_abundance_binding_raw`
  - `paired_abundance_binding_biophysical_binding`
- Evaluated branch-level evidence using abundance KS, binding functional KS,
  single-mutant holdout Spearman, double-mutant holdout Spearman, epistasis
  prediction Spearman/KS, fitted reference-to-peak geometry, and fitted
  epistatic parameters.

## Outputs

- Figures: none
- Tables: `tables/run-028_key_metrics.md`

## Main observations

- No prior `ANA-*` or `RES-*` record analysed `RUN-028`. The run was hanging in
  the scientific record as a completed run linked to `EXP-009`, but without a
  corresponding analysis/result.
- `RUN-028` completed the intended matched KRAS panel diagnostic. The selected
  assays were `RASK_HUMAN_Weng_2022_abundance` and
  `RASK_HUMAN_Weng_2022_binding-DARPin_K55`, both length `188`, with `26,012`
  and `24,873` measured variants respectively.
- The observation layer was strong enough that the negative result is not
  primarily explained by failed assay score modeling. The abundance MAVE-NN
  model reached test Spearman `0.818` and NRMSE `0.614`; the binding MAVE-NN
  model reached test Spearman `0.933` and NRMSE `0.362`.
- The binding-only thermodynamic branch improved several metrics over the
  binding-only raw branch: binding KS improved from `0.506` to `0.384`,
  double-mutant holdout Spearman improved from `0.195` to `0.298`, epistasis
  prediction Spearman improved from `0.106` to `0.172`, and epistasis KS
  improved from `0.673` to `0.481`.
- That binding-only thermodynamic gain was not complete. Single-mutant holdout
  Spearman fell slightly from `0.308` to `0.293`, and the fitted reference
  remained essentially at the peak (`fraction_of_peak 0.999`, distance `1`).
- Adding the abundance readout without the thermodynamic binding readout did
  not rescue the model. `paired_abundance_binding_raw` improved binding KS
  relative to `binding_only_raw` (`0.381` versus `0.506`) but reduced
  single-mutant holdout Spearman (`0.265` versus `0.308`), left epistasis
  prediction essentially unchanged, and still placed the reference near the
  peak (`fraction_of_peak 0.999`, distance `6`).
- The full paired abundance-plus-binding thermodynamic branch produced the
  clearest mechanistic-looking gain on epistasis prediction and binding KS.
  Relative to the paired raw branch, binding KS improved from `0.381` to
  `0.258`, epistasis prediction Spearman improved from `0.104` to `0.287`, and
  epistasis KS improved from `0.488` to `0.428`.
- The same full paired thermodynamic branch failed the core predictive package.
  Abundance KS worsened from `0.582` to `0.832`, single-mutant holdout
  Spearman fell from `0.265` to `0.105`, and double-mutant holdout Spearman
  fell from `0.252` to `0.144`.
- The peak-geometry problem remained unresolved. All branches placed the fitted
  reference at or very near the fitted peak: `fraction_of_peak` ranged from
  `0.9986` to `1.000`, with distance `0` to `6`.
- The full paired thermodynamic branch also pushed the fitted empirical
  pairwise strength to the top tested value (`0.05`) while keeping
  `epistasis_strength = 0.0`. That is more consistent with remaining structural
  mismatch or insufficient identifiability than with a clean mechanistic
  recovery.

## Result records created

- `RES-012`

## Hypothesis updates

- `HYP-007` is weakened for the current KRAS paired-readout implementation. The
  matched abundance-plus-binding thermodynamic branch improved some
  mechanistic diagnostics, especially epistasis prediction and binding KS, but
  it did not materially outperform simpler controls on the full set of primary
  readouts specified by `HYP-007`.
- `HYP-001` is not directly updated. `RUN-028` is a focused matched-readout
  diagnostic for `HYP-007`, not a broad shared-regime empirical recovery test.
