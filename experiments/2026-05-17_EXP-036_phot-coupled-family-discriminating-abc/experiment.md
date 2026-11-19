# EXP-036: PHOT coupled-family discriminating ABC

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008

## Question

For PHOT, does ABC fitting of the same RUN-079 coupled candidate families produce a calibration-selected held-out graphlet improvement that is specific to empirical local modules rather than no-module, random/shuffled-module, or uncoupled controls?

## Pre-experiment prediction

Support for HYP-008 requires the ABC-selected empirical-module family to beat fitted-current, smooth-core, no-module, random-module, shuffled-module, and uncoupled-module controls on held-out graphlet distance while satisfying component guardrails for WHT high-order fraction, local maxima, accessibility, lag-1 autocorrelation, SE/RSE, and GFT bands. If the RUN-089 diagnosis is correct, the fitted no-module or uncoupled controls will remain competitive and the empirical-module specificity claim will fail.

## Rationale

RUN-089 showed that the RUN-079 no-module PHOT graphlet-distance gain is stable against current and smooth controls but fails specificity because the uncoupled-module control is better. This experiment uses ABC, which previously improved calibration performance over deterministic grids, to fit the same family/control set directly under calibration graphlet objectives before evaluating held-out graphlets.

## Experimental design

- Freeze the PHOT current-model base from `RUN-077`.
- Rebuild the same complete PHOT graphlet panel and empirical plaquette-module
  catalog used by the coupled-family diagnostics.
- Fit `no_module`, `empirical_module`, `random_module`, `shuffled_module`, and
  `uncoupled_module` independently with SMC-ABC on calibration graphlets only.
- Add component guardrails to the ABC objective for SE/RSE, local maxima,
  accessibility, lag-1 autocorrelation, WHT high-order fraction, and GFT bands.
- Compare calibration-selected posterior candidates on held-out graphlets
  against fitted current, smooth core, residual-only, GNK, no-module,
  random/shuffled-module, and uncoupled controls.

## Inputs

- `RUN-077` PHOT system-specific current-model refit outputs.
- `scripts/phot_coupled_family_discriminating_abc.py`
- `scripts/phot_coupled_landscape_grid_search.py`
- `src/project_adapt_env/coupled_landscape.py`
- `src/project_adapt_env/smc_abc.py`

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- preflight.json
- progress.json
- parameter_specs.csv
- posterior_particles.csv
- posterior_rounds.csv
- evaluated_parameter_sets.csv
- graphlet_catalog.csv
- module_catalog.csv
- graphlet_metrics.csv
- walsh_terms.csv
- gft_coefficients.csv
- graphlet_summary.csv
- model_distances.csv
- top_model_distances.csv
- selected_family_distances.csv
- selected_family_pairwise_deltas.csv
- component_guardrail_summary.csv
- module_energy_rows.csv
- module_energy_by_graphlet.csv
- module_energy_summary.csv
- base_correlation.csv
- empirical_readout_fit.json

## Analysis plan

- Create an ANA record after execution.

## Completion criteria

- [x] Config recorded
- [x] Run script recorded
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- `RUN-091` completed and was collected under
  `data/processed/phot-coupled-family-discriminating-abc/RUN-091`.
- `ANA-036` and `RES-036` record the conclusion: ABC fitting improves the fixed
  held-out PHOT graphlet distance relative to current/smooth/no-module controls,
  but random-module is effectively tied, uncoupled-module is better, and the
  selected empirical-module energy remains diffuse.
