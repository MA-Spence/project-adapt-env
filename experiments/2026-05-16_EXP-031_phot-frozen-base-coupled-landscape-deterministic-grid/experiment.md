# EXP-031: PHOT frozen-base coupled-landscape deterministic grid search

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008

## Question

With PHOT as the only empirical test system, can a frozen RES-030 current-model base plus a constrained correlated local auxiliary landscape improve held-out dense graphlet topology, spectra, and localized plaquette-energy structure without rerunning ABC?

## Pre-experiment prediction

If HYP-008 is plausible, local-module coupled candidates should improve PHOT held-out graphlet distance relative to the frozen current base, smooth core, diffuse residual, uncoupled auxiliary, shuffled/random modules, and GNK controls while preserving base correlation and concentrating innovation plaquette energy in few modules. If improvements are matched by shuffled/random/uncoupled controls, or require diffuse innovation or base-correlation loss, the coupled-landscape mechanism is not supported.

## Rationale

RES-030 showed that PHOT system-specific current refitting partially rescues graphlet dynamics but leaves a high-order spectral-envelope mismatch and weak smooth-core separation. This experiment avoids another expensive ABC run by freezing the RES-030 PHOT base and deterministically enumerating the coupled-landscape parameter space needed for a first HYP-008 machinery test.

## Experimental design

- Use `PHOT_CHLRE_Chen_2023` as the only empirical system.
- Restore the completed `RUN-077` PHOT current-model outputs from DVC and use
  them as fixed inputs. No ABC or posterior refitting is performed.
- Reconstruct the RES-030 PHOT frozen bases on the same complete dense Hamming
  graphlets: fitted current best, fitted smooth core, and posterior-mean
  sensitivity control.
- Fit a polynomial global readout from frozen current-best graphlet values to
  empirical `latent_phi` using calibration graphlets only. Define empirical
  innovation as empirical `latent_phi` minus that calibration readout.
- Select local modules from calibration empirical innovation plaquette energy,
  then evaluate held-out graphlets for empirical innovation localization.
- Deterministically grid-search coupled auxiliary landscapes over frozen
  current-best and smooth-core anchors. Compare local modules, no-module
  coupled controls, shuffled modules, random modules, uncoupled auxiliary
  controls, GNK, residual-only, and the existing frozen current/smooth models.
- Record standard graphlet metrics plus base correlation, edge-gradient
  correlation, edge-energy concentration, plaquette-energy concentration,
  effective module count, top-module energy fraction, and model distance to the
  empirical graphlet package.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `RUN-077` PHOT fit summary, prepared PHOT assay, MAVE-NN metrics, and copied
  PHOT alignment.
- Current root `src/` and `scripts/` code.

## Configuration

- See `config.yaml`.
- The deterministic grid is intentionally small enough to avoid the cost of a
  full ABC refit while still testing the coupled-landscape machinery.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- preflight.json
- graphlet_catalog.csv
- module_catalog.csv
- parameter_grid.csv
- graphlet_metrics.csv
- walsh_terms.csv
- gft_coefficients.csv
- model_distances.csv
- top_model_distances.csv
- module_energy_rows.csv
- module_energy_by_graphlet.csv
- module_energy_summary.csv
- base_correlation.csv
- empirical_readout_fit.json

## Analysis plan

- Create an ANA record after execution.
- The analysis should identify the best calibration and held-out test
  candidates, but the result must not treat calibration improvement as support
  unless it survives on held-out graphlets.
- The primary falsifiers are: shuffled/random modules match local modules;
  uncoupled auxiliary matches coupled auxiliary; improvement requires base
  correlation loss; innovation plaquette energy is diffuse; GNK matches the
  same statistics; or spectra improve while accessibility/local maxima/random
  walk statistics remain wrong.

## Completion criteria

- [ ] Config committed
- [ ] Run script committed
- [ ] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
