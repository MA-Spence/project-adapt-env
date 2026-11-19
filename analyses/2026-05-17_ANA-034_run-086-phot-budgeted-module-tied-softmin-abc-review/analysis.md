# ANA-034: RUN-086 PHOT budgeted module-tied soft-min ABC review

## Purpose

Assess whether `RUN-086`, the PHOT budgeted module-tied soft-min ABC run,
supports `HYP-008` by improving calibration-selected held-out graphlet
topology, spectra, locality, and assay guardrails relative to frozen `RUN-077`
controls.

## Linked experiments/runs

- Experiments: `EXP-034`
- Runs: `RUN-086`
- Source data:
  `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086`
- DVC pointer:
  `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086.dvc`

## Scratch notes read

- `scratch/coupled_landscape_nonlinear_operator_plan.md`
- `scratch/coupled_landscape_regularisation_metrics.md`
- `scratch/model_signature_derivation.md`
- `scratch/model_literature_falsification_review.md`

These notes constrain the interpretation. Support for `HYP-008` requires a
calibration-selected model to improve held-out graphlet statistics while
showing empirical-like, local, non-diffuse innovation. Improvement on
calibration graphlets, assay fit, or test-oracle particles alone is not enough.

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-086_preflight_and_output_audit.md`
- `tables/run-086_smc_posterior_summary.md`
- `tables/run-086_graphlet_distance_comparison.md`
- `tables/run-086_selected_graphlet_metrics.md`
- `tables/run-086_assay_fit_comparison.md`
- `tables/run-086_locality_tests.md`
- `tables/run-086_hyp008_decision_table.md`

## Run and data audit

`RUN-086` completed and was collected with `labproj collect`. The promoted
output directory contains all `27` expected runtime output files and is
DVC-tracked under
`data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086.dvc`.

The run used `PHOT_CHLRE_Chen_2023` only, restored the `RUN-077` PHOT
system-specific current-model output, kept `a(x)` fixed to the fitted smooth
core, and fitted the budgeted state-difference operator
`g_tau(a,d) = a - tau log(1 + exp(-d/tau))`. It evaluated `167,529` assay
variants and `24` dense graphlets split evenly into `12` calibration and `12`
held-out test graphlets. The SMC-ABC fit used `32` particles, `3` rounds, and
`160` candidate attempts per round, with `26` posterior candidates evaluated on
the graphlet outputs.

## Key observations

- The ABC-selected best particle slightly improves calibration graphlet
  distance over the frozen current and smooth-core controls: `0.7240` versus
  `0.7318` for fitted current and `0.7454` for smooth core. This is weaker than
  `RUN-084` and far behind the empirical readout-residual diagnostic
  (`0.5428`), so it is a small calibration signal only.
- The central held-out criterion fails. The calibration-selected best particle
  has test graphlet distance `0.4060`, worse than fitted current (`0.4016`).
  It is slightly better than smooth core (`0.4119`), but that is not enough
  under the registered prediction because fitted current remains the relevant
  current-model control.
- Module controls do not support empirical module locality. The selected
  empirical-module particle is not better than the shuffled-module version on
  held-out graphlets (`0.4060` versus `0.4029`), and random-module and
  no-module variants are close (`0.4081` and `0.4144`).
- The fitted operator is effectively not using the local state mechanism. The
  selected particle has values-to-base correlation `1.0`, edge-gradient
  correlation `1.0`, innovation variance ratio `1.4e-25`, and edge/plaquette
  innovation ratios to the empirical residual of `6.0e-22` and `1.2e-21`.
  This is not empirical-like localized ruggedness; it is nearly the frozen
  base landscape.
- The PHOT graphlet package remains mismatched. Empirical held-out graphlets
  have local maxima `5.17`, accessibility `0.1975`, lag-1 autocorrelation
  `0.1905`, and WHT high-order fraction `0.3788`. The selected budgeted
  particle has local maxima `3.50`, accessibility `0.0280`, autocorrelation
  `0.0851`, and WHT high-order fraction `0.4984`.
- The assay guardrail is negative. The best assay objective is still the
  fitted-current posterior mean (`0.4631`). The best budgeted candidate is
  worse (`0.5350`), and the calibration-selected best particle is worse again
  (`0.5560`). Holdout ranking is essentially inherited from the base and
  epistasis-prediction Spearman remains negative for the selected particle
  (`-0.1745`).

## Interpretation

`RUN-086` is a valid execution of `EXP-034`, but it does not support the tested
budgeted module-tied implementation of `HYP-008`. The added budgets prevent the
large artificial innovation seen in `RUN-084`, but they do so by collapsing the
selected operator back onto the frozen base. That resolves the over-innovation
failure mode technically, not scientifically.

The best calibration-selected model does not beat fitted current on held-out
graphlets, does not improve the assay guardrail, and does not validate the
empirical module mechanism against no-module/random/shuffled-module controls.
The result therefore weakens this implementation of the local nonlinear
trait-operator hypothesis.

## Limitations

- PHOT is the only empirical system tested.
- The empirical graphlet target remains the processed `latent_phi` layer rather
  than direct exhaustive empirical measurement.
- The base landscape is frozen to the `RUN-077` PHOT fit, so this result tests
  the budgeted operator layer rather than a joint refit.
- The budgets may be too restrictive, but relaxing them would need to avoid
  returning to the `RUN-084` artificial-innovation failure mode.
- The run provides a constrained negative result for this implementation, not
  a formal refutation of every possible nonlinear local trait operator.
