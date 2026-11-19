# EXP-032: PHOT hard-min coupled weakest-link landscape grid search

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008

## Question

With PHOT as the only empirical test system, does a frozen fitted smooth/additive PHOT landscape combined with a correlated offset auxiliary landscape by f(x)=min(a(x),b(x)) improve assay fit and held-out graphlet local ruggedness relative to the single-landscape base and fitted-current controls?

## Pre-experiment prediction

If weakest-link coupling is the missing local nonlinear operator, some correlated hard-min grid points should improve PHOT assay fit statistics and held-out graphlet WHT/GFT/topology while producing empirical-like localized edge/plaquette innovation energy around switching boundaries. If improvements are only diffuse, only graphlet-specific, or worsen assay fit, this implementation does not support HYP-008.

## Rationale

RES-031 showed that additive/gated no-module coupling could move PHOT graphlet distance but did not create validated localized ruggedness. The user-proposed hard-min operator is biologically different: two related essential functions are constrained by the weaker trait, so local topology may arise at trait-switch boundaries rather than from additive innovation modules. This experiment freezes the RES-030 PHOT fit and performs a deterministic grid search before any new ABC.

## Experimental design

- Use `PHOT_CHLRE_Chen_2023` as the only empirical system.
- Restore `RUN-077`/`RES-030` PHOT outputs from DVC and rebuild the fitted
  current, posterior-mean, and smooth-core landscapes. No ABC is run.
- Define the primary first landscape as `A(x) = fitted_smooth_core_only(x)`.
  The fitted current model remains a baseline/control, not the primary `A`.
- Construct auxiliary landscapes as correlated, scaled, and offset functions
  of `A` plus deterministic additive/low-order auxiliary components. Scan
  auxiliary correlation, scale, offset, pairwise auxiliary scale, pairwise
  density, and seed.
- Combine traits with the strict weakest-link operator
  `F(x) = min(A(x), B(x))`.
- Evaluate the full prepared PHOT assay for the same fit-statistic families
  used in recent PHOT refits: held-out Spearman/NRMSE/KS, single-mutant DFE,
  double-mutant fit, and double-mutant epistasis prediction.
- Evaluate complete dense PHOT graphlets for WHT, GFT, SE/RSE, local maxima,
  accessibility, random-walk autocorrelation, and model distance to empirical
  `latent_phi`, using the same graphlet metrics as `RUN-077`/`RUN-079`.
- Measure empirical local signal as empirical `latent_phi` minus a polynomial
  readout of the frozen smooth core on calibration graphlets, then compare
  empirical and model edge/plaquette energy concentration, graphlet energy
  concentration, and hard-min switching-boundary metrics.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `RUN-077` PHOT fit summary, prepared PHOT assay, MAVE-NN metrics, and PHOT
  alignment.
- Current root `src/` and `scripts/` code.

## Configuration

- See `config.yaml`.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- preflight.json
- parameter_grid.csv
- assay_fit_metrics.csv
- assay_dfe_metrics.csv
- assay_epistasis_metrics.csv
- assay_fit_rankings.csv
- graphlet_catalog.csv
- graphlet_metrics.csv
- walsh_terms.csv
- gft_coefficients.csv
- graphlet_summary.csv
- model_distances.csv
- top_model_distances.csv
- local_signal_by_graphlet.csv
- local_signal_energy_by_feature.csv
- local_signal_feature_energy_summary.csv
- local_signal_graphlet_energy_summary.csv
- locality_covariance_summary.csv
- empirical_readout_fit.json

## Analysis plan

- Create an ANA record after execution.
- First ask whether any hard-min grid point improves full-assay PHOT fit over
  both `fitted_smooth_core_only` and `fitted_current_best`. A graphlet-only
  improvement is not sufficient support if assay fit degrades.
- Then ask whether the same or calibration-selected hard-min candidates improve
  held-out graphlet WHT/GFT/topology distance relative to A-only, fitted
  current, empirical readout residual, and GNK controls.
- Finally ask whether the innovation is local and empirical-like: edge and
  plaquette energy should be concentrated across graphlets/features, and model
  hard-min switch boundaries should explain local plaquette energy better than
  diffuse no-boundary innovation.
- Treat calibration-selected versus test graphlet performance separately.

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- `RUN-081` completed successfully with scheduler exit code `0`.
- Durable outputs are DVC-tracked at
  `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081.dvc`.
- `ANA-032` and `RES-032` record the interpretation. The hard-min operator
  improved the PHOT full-assay objective but did not beat fitted-current
  held-out graphlet controls or validate empirical-like local non-diffuse
  ruggedness.
