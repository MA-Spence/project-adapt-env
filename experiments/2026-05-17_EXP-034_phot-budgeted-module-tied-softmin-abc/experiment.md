# EXP-034: PHOT budgeted module-tied soft-min ABC

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008
- Previous experiments: EXP-030, EXP-031, EXP-032, EXP-033
- Motivating results: RES-030, RES-031, RES-032, RES-033

## Question

For PHOT, does fitting an explicit module-tied state-difference field
`d(x)=b(x)-a(x)`, with finite-difference energy and state-boundary budgets, allow
a frozen-base soft-min operator to improve calibration-selected held-out
graphlet topology, spectra, and locality relative to fitted-current,
smooth-core, residual, random/shuffled-module, no-module, and GNK controls?

## Pre-experiment prediction

A supporting result requires calibration-selected posterior candidates to beat
fitted-current and smooth-core controls on held-out graphlet distance while
preserving base correlation and keeping innovation edge/plaquette energy within
the empirical residual budget. Assay-objective improvement alone remains
insufficient support for HYP-008.

## Rationale

RES-033 showed that the previous soft-min/Boltzmann ABC family could improve
calibration graphlet distance but generated very large, non-empirical
state-transition innovation and failed held-out graphlets. EXP-034 tests the
specific diagnosis that `b(x)` should be parameterised through a module-tied
`d(x)=b(x)-a(x)` field and constrained by explicit smoothness/energy budgets,
rather than by a globally correlated auxiliary landscape alone.

## Experimental Design

- Freeze the RUN-077 PHOT current-model output and use the fitted smooth-core
  landscape as `a(x)`.
- Fit calibration graphlets only; reserve held-out graphlets for scientific
  support.
- Fit the state-difference field:

```text
d(x) = bias + slope * z(a(x)) + module_scale * delta_M(x)
g_tau(a,d) = a - tau log(1 + exp(-d/tau))
```

- Select module coordinates from calibration empirical readout-residual
  plaquette energy, then evaluate no-module, shuffled-module, and random-module
  controls for the selected posterior summaries.
- Penalize high innovation variance, edge energy, plaquette energy,
  state-boundary roughness, and overly low effective temperature relative to the
  local edge scale of `d(x)`.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `scripts/phot_budgeted_softmin_coupled_landscape_abc.py`
- `src/project_adapt_env/coupled_landscape.py`

## Configuration

- `config.yaml`

## Execution

- Run via `run.sh` or `labproj submit EXP-034 --target lab-slurm`.

## Expected Outputs

- `summary.json`
- `preflight.json`
- `progress.json`
- `parameter_specs.csv`
- `module_catalog.csv`
- `posterior_particles.csv`
- `posterior_rounds.csv`
- `evaluated_parameter_sets.csv`
- `assay_fit_metrics.csv`
- `assay_dfe_metrics.csv`
- `assay_epistasis_metrics.csv`
- `assay_fit_rankings.csv`
- `graphlet_catalog.csv`
- `graphlet_metrics.csv`
- `walsh_terms.csv`
- `gft_coefficients.csv`
- `graphlet_summary.csv`
- `model_distances.csv`
- `top_model_distances.csv`
- `local_signal_by_graphlet.csv`
- `local_signal_energy_by_feature.csv`
- `local_signal_feature_energy_summary.csv`
- `local_signal_graphlet_energy_summary.csv`
- `locality_covariance_summary.csv`
- `empirical_readout_fit.json`

## Analysis Plan

Create `ANA-034` after execution. The analysis must report calibration-selected
held-out graphlet distance, topology/spectral/dynamics metrics, innovation
energy ratios against the empirical residual, state-boundary smoothness,
module-control comparisons, assay guardrail metrics, and whether HYP-008 is
supported or weakened.

## Completion Criteria

- [x] Config recorded
- [x] Run script recorded
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated before submission
- [x] PROJECT_STATE.md updated after result

## Post-experiment notes

- `RUN-086` completed and was collected under
  `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086`.
- `ANA-034` and `RES-034` record the negative interpretation: the budgeted
  module-tied operator slightly improves calibration graphlet distance but does
  not beat fitted current on held-out graphlets, does not validate module
  locality, and does not improve the assay guardrail.
