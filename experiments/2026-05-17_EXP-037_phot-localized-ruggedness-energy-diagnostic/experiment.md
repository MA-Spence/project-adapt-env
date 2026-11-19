# EXP-037: PHOT localized ruggedness energy diagnostic

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008

## Question

In PHOT empirical, frozen ABC-refit current a(x), smooth-core, RUN-091 ABC fitted coupled families, and the random/shuffled/uncoupled controls, is Dirichlet ruggedness locally clustered on neighboring graphlet edges, sites, plaquettes, or empirical modules rather than diffuse across the dense graphlet landscape?

## Pre-experiment prediction

Support for HYP-008 requires PHOT empirical residual ruggedness and the empirical-module ABC family to show stronger localized energy clustering than fitted current, smooth core, no-module, random-module, shuffled-module, and uncoupled-module controls: positive adjacent-edge Dirichlet-energy autocorrelation, lower effective edge/site/pair counts, higher top-edge/site/pair fractions, and empirical-module enrichment. If RUN-091's diffuse-module diagnosis is correct, the empirical-module family will not show specific localization beyond random or uncoupled controls.

## Rationale

RUN-091 improved held-out graphlet distance with ABC but failed mechanism specificity: random-module was fractionally better, uncoupled-module was better, and empirical-module energy was diffuse. This experiment does not introduce a new operator. It reuses the fitted RUN-091 parameters and directly tests whether the observed and fitted landscapes localize ruggedness or distribute Dirichlet energy broadly.

## Experimental design

- Freeze the PHOT graphlet/source context from `RUN-077` and rebuild the same
  graphlet/module panel used by `RUN-091`.
- Load the fitted RUN-091 ABC posterior-predictive parameter sets and the
  calibration-selected family/control candidates; do not refit ABC.
- Evaluate PHOT empirical `latent_phi`, empirical readout residual, fitted
  current `a(x)`, fitted-current posterior mean, smooth core, residual-only,
  GNK local null, and the RUN-091 no-module, empirical-module, random-module,
  shuffled-module, and uncoupled-module ABC candidates.
- For each graphlet, compute Dirichlet edge energy `(f(v)-f(u))^2`, adjacent
  edge-energy Pearson autocorrelation and Moran I on the hypercube line graph,
  edge/node/site/plaquette energy concentration, and module energy enrichment.
- Report total-landscape localization for all systems/models and
  innovation-from-base localization for the ABC coupled-family candidates.

## Inputs

- `RUN-077` PHOT system-specific current-model refit outputs.
- `RUN-091` PHOT coupled-family discriminating ABC outputs.
- `scripts/phot_localized_ruggedness_energy_diagnostic.py`
- `scripts/phot_coupled_family_discriminating_abc.py`
- `scripts/phot_coupled_landscape_grid_search.py`
- `src/project_adapt_env/coupled_landscape.py`

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- preflight.json
- localized_ruggedness_graphlet_metrics.csv
- edge_dirichlet_energy_rows.csv
- module_enrichment_by_graphlet.csv
- model_localization_summary.csv
- module_enrichment_summary.csv
- selected_vs_empirical_localization_delta.csv
- localization_metric_ranks.csv
- run091_selected_family_distances.csv
- run091_reused_parameter_sets.csv

## Analysis plan

- Analysis completed in `ANA-037`.

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

- RUN-093 completed on lab-slurm with scheduler job 244 and inferred Slurm exit
  code 0.
- Collected outputs are tracked at
  `data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093.dvc`
  (`md5: 8747ed319b4f1a5fed2230d543cc5d5b.dir`, 12 files,
  175474957 bytes) and were backed up to the labcluster DVC remote.
- ANA-037 and RES-037 record the interpretation: PHOT empirical and residual
  readout landscapes show residual localized ruggedness relative to fitted
  current/smooth-core controls, but the empirical-module ABC candidate is not
  uniquely supported over no-module, random-module, uncoupled-module, residual,
  or GNK controls.
