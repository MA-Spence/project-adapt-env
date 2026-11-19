# EXP-038: PHOT ABC local ruggedness objective

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008

## Question

For PHOT, does adding adjacent-edge Dirichlet-energy localization statistics to the SMC-ABC calibration objective yield coupled-family candidates whose held-out graphlet fit and ruggedness diffuseness match empirical PHOT better than current/smooth, no-module, random/shuffled-module, and uncoupled controls?

## Pre-experiment prediction

If the RUN-093 diagnosis is actionable, including adjacent-edge Dirichlet-energy autocorrelation and edge/site concentration terms in the ABC objective should select particles that better match empirical PHOT held-out ruggedness localization without losing the RUN-091 held-out graphlet-distance gain. Support for HYP-008 still requires empirical-module specificity over no-module, random/shuffled-module, uncoupled-module, residual, and GNK controls.

## Rationale

RUN-091 showed ABC can improve PHOT held-out graphlet distance but failed locality specificity, while RUN-093 showed PHOT has residual localized ruggedness that was not optimized during ABC fitting. This experiment keeps the same PHOT base and candidate/control families, but adds explicit local ruggedness diffuseness statistics to the calibration objective.

## Experimental design

- Freeze the PHOT fitted-current base and graphlet/source context from
  `RUN-077`.
- Fit the same `RUN-091` family/control set with SMC-ABC:
  `no_module`, `empirical_module`, `random_module`, `shuffled_module`, and
  `uncoupled_module`.
- Keep calibration/test separation across graphlets; ABC sees calibration
  graphlets only.
- Extend the ABC calibration distance vector with adjacent-edge Dirichlet-energy
  localization statistics: edge-energy Moran I, edge top-10 fraction, edge
  effective fraction, site top-1/effective fraction, and plaquette
  top-1/effective fraction.
- Retain the original component metrics for SE/RSE, local maxima,
  accessibility, lag-1 random-walk autocorrelation, WHT high-order fraction,
  and GFT low/mid/high energy.
- Evaluate calibration-selected particles on held-out graphlets against fitted
  current, smooth core, empirical readout residual, GNK, no-module,
  empirical-module, random/shuffled-module, and uncoupled-module controls.

## Inputs

- `RUN-077` PHOT system-specific current-model refit outputs.
- `scripts/phot_coupled_family_discriminating_abc.py`
- `scripts/phot_coupled_landscape_grid_search.py`
- `scripts/phot_localized_ruggedness_energy_diagnostic.py`
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

- Local quick validation passed with the project `.venv`; the compact run
  emitted the expected output files and confirmed that the new local ruggedness
  metrics are present in `preflight.json`, `graphlet_metrics.csv`, and
  `posterior_particles.csv`.
- Dry-run submission: `RUN-094`.
- `RUN-095` was aborted during source staging before a scheduler job id was
  assigned, then marked failed to keep the record accurate.
- Live submission: `RUN-096`, Slurm job `246`, running on `dutch` immediately
  after submission.
- `RUN-096` completed with inferred exit code `0`. Outputs were collected with
  `labproj collect` and DVC-tracked/backed up under
  `data/processed/phot-abc-local-ruggedness-objective/RUN-096`.
- `ANA-038` records the analysis and `RES-038` records the result. Adding local
  ruggedness terms to the ABC objective improves held-out distances over
  fitted current/smooth core, but the calibration-selected no-module and
  uncoupled controls beat empirical-module, and empirical-module energy remains
  diffuse across modules.
