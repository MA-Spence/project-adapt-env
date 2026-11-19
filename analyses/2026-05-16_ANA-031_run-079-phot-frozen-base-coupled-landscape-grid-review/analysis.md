# ANA-031: RUN-079 PHOT frozen-base coupled-landscape grid review

## Purpose

Assess whether the deterministic frozen-base PHOT coupled-landscape grid
supports `HYP-008` by improving held-out graphlet spectra, topology, and
localized innovation relative to matched controls.

## Linked experiments/runs

- Experiments: `EXP-031`
- Runs: `RUN-079`
- Source data:
  `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079`
- DVC pointer:
  `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079.dvc`

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-079_preflight_and_output_audit.md`
- `tables/run-079_empirical_target_graphlet_metrics.md`
- `tables/run-079_calibration_selected_distances.md`
- `tables/run-079_calibration_selected_by_module_source.md`
- `tables/run-079_static_control_distances.md`
- `tables/run-079_test_oracle_distances_capacity_only.md`
- `tables/run-079_relevant_test_graphlet_summaries.md`
- `tables/run-079_base_correlation_heldout.md`
- `tables/run-079_module_energy_locality_heldout.md`
- `tables/run-079_module_catalog.md`
- `tables/run-079_hyp008_decision_table.md`

## Run and data audit

`RUN-079` completed successfully with scheduler exit code `0`. The collected
output directory contains the planned graphlet catalog, module catalog,
parameter grid, graphlet metrics, WHT terms, GFT coefficients, model distances,
base-correlation diagnostics, module-energy locality tables, and summary JSON
files. The output directory is DVC-tracked and pushed.

The run used `PHOT_CHLRE_Chen_2023` only, froze the `RUN-077`/`RES-030` PHOT
current-model fit, did not rerun ABC, and enumerated `144` deterministic
coupled parameter sets over `24` dense graphlets split evenly into `12`
calibration and `12` test graphlets across dimensions `4`, `5`, and `6`. It
selected `12` empirical modules from calibration graphlet plaquette energy.

The empirical readout residual control fits a degree-3 polynomial from the
frozen `fitted_current_best` base to empirical `latent_phi`. Its calibration
`R^2` is only `0.4426`, so it is useful as a scale/readout diagnostic but not
as a biologically grounded procedural model.

## Selection rule

The grid contains both calibration and test distances. For the scientific
claim, the valid comparison is calibration-selected then held out on test
graphlets. Test-oracle rankings are included only as a capacity diagnostic.

This distinction matters. The globally best calibration parameter is an
uncoupled local-module setting with calibration distance `0.5219`, but its
held-out test distance is only `0.3951`, barely better than the frozen current
model (`0.4016`). Within candidate families, the calibration-selected
no-module coupled operator generalizes much better (`0.3137`), but that is not
evidence for localized empirical modules.

## Key observations

- The empirical PHOT test graphlets remain rugged and high-order-rich:
  sign epistasis fraction `0.6146`, RSE fraction `0.2583`, mean local maxima
  `5.17`, accessible path fraction `0.1975`, random-walk lag-1
  autocorrelation `0.1905`, WHT high-order fraction `0.3788`, and GFT
  low/mid/high fractions `0.6212`/`0.2236`/`0.1552`.
- The empirical global-readout residual is the closest test candidate
  (`0.2708`) but is a diagnostic subtraction from the empirical target, not a
  procedural generator.
- The calibration-selected no-module coupled operator beats the frozen current
  and smooth-core controls on held-out distance: `0.3137` versus `0.4016` for
  `fitted_current_best` and `0.4119` for `fitted_smooth_core_only`.
- The calibration-selected empirical plaquette-module coupled candidate does
  not show the expected local-module advantage. Its held-out distance is
  `0.3982`, worse than no-module (`0.3137`), worse than random-module
  (`0.3546`), and effectively comparable to shuffled-module (`0.3946`) and
  frozen current (`0.4016`).
- The test-oracle capacity table tells the same story in weaker form: the best
  no-module candidate is `0.3096`, the best uncoupled local-module candidate is
  `0.3251`, the best empirical-module coupled candidate is `0.3296`, and the
  best shuffled/random module candidates are `0.3385`/`0.3425`. Local empirical
  modules are not uniquely favored even when selection is allowed to look at
  the test split.
- Coupled candidates preserve the frozen base well. The calibration-selected
  no-module operator has held-out base phenotype correlation `0.982` and
  edge-gradient correlation `0.987`; the empirical-module coupled candidate
  has `0.986` and `0.988`.
- The locality diagnostic is negative. The calibration-selected empirical
  module candidate has plaquette effective module count `9.33` out of `12` and
  normalized entropy `0.942`; the empirical readout residual itself has
  effective count `9.72` and entropy `0.952`. That is diffuse, not concentrated
  local innovation. A random-module candidate is more concentrated
  (`4.62` effective modules) without being the best model.

## Interpretation

`RUN-079` is procedurally successful but scientifically negative for the
specific local-module version of `HYP-008`. A correlated nonlinear operator can
improve the frozen PHOT graphlet package without rerunning ABC, and it can do
so while preserving high base correlation. However, the improvement is carried
by the no-module operator and by an empirical readout residual diagnostic, not
by localized empirical plaquette modules.

This does not refute the broad possibility that real PHOT ruggedness arises
from nonlinear trait composition. It does weaken the current implementation
path: selecting local modules from calibration plaquette energy and injecting
them as coupled module residuals does not produce a clearly localized,
held-out-improving procedural mechanism. The inferred innovation energy remains
spread over many selected modules, and matched random/shuffled/no-module
controls are too competitive.

Relative to `RES-030`, this run does not rescue the PHOT spectral-envelope
problem as a local-module mechanism. The empirical test WHT high-order fraction
is `0.3788`; the fitted current model is `0.4984`, the
calibration-selected empirical-module coupled candidate is `0.4966`, and the
calibration-selected no-module candidate is `0.4713`. The graphlet distance
can be improved, but the localized ruggedness explanation is not yet supported.

## Limitations

- This is a deterministic PHOT-only grid, not an ABC refit and not a
  multi-system validation.
- Test-oracle rankings are only capacity diagnostics; they are not valid model
  selection evidence.
- The empirical target remains based on the `latent_phi` layer from the current
  PHOT empirical processing, not direct exhaustive empirical graphlets.
- The module construction may be too crude: many selected modules share site
  `111`, and the module-energy summaries show high effective module counts
  rather than sparse patch localization.
- The run tests one implementation of local coupled innovation, not every
  possible nonlinear trait-operator formulation.
