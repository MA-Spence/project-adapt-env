# EXP-039: PHOT nested soft-min ensemble capacity diagnostic

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008

## Question

Does adding a locked third latent residual landscape through nested independently parameterised Boltzmann/soft-min coupling improve PHOT held-out graphlet and local-ruggedness performance beyond the locked first-stage coupled model and matched negative controls?

## Pre-experiment prediction

If prior improvements mostly reflect latent residual capacity, a nested third landscape should improve held-out graphlet distance and local ruggedness statistics over the locked first-stage model, but matched random/shuffled/uncoupled controls may remain close; a mechanistic interpretation requires empirical-module specificity beyond these controls.

## Rationale

RUN-091 and RUN-096 show that secondary residual capacity improves PHOT held-out graphlet and local-ruggedness metrics, but does not identify residue-specific empirical modules. This experiment tests whether increasing ensemble capacity with a second independently parameterised coupling layer gives stable predictive improvement, and whether any gain is specific to empirical module structure or simply generic latent capacity.

## Experimental design

- Use the frozen PHOT fitted base and empirical graphlets from RUN-077.
- Lock the first soft-min/Boltzmann layer to the best RUN-084 posterior parameterisation:
  `c(x)=g(a(x),b(x))`, where `a(x)` is the fitted smooth thermostability core.
- Fit only the third latent residual landscape and second coupling layer:
  `f(x)=h(c(x),d(x))`.
- Fit independent ABC posterior populations for five second-stage families:
  no-module, empirical-module, random-module, shuffled-module, and uncoupled-module.
- Select parameters on calibration graphlets only. Evaluate the scientific decision point on held-out graphlets.
- Primary objective components: graphlet distance to empirical `latent_phi`, including sign/reciprocal sign epistasis, local maxima, accessibility, random-walk lag-1 autocorrelation, Walsh/GFT spectral fractions, adjacent-edge Dirichlet-energy autocorrelation, and edge/site/plaquette energy concentration.
- Controls: locked first-stage model, fitted current model, smooth core, current posterior mean, residual-only diagnostic, empirical readout residual, GNK null, no-module latent capacity, random-module, shuffled-module, and uncoupled-module families.
- Interpretation boundary: improvement over the locked first stage supports added latent capacity. Empirical-module specificity is required before making any claim about biologically localised module structure.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- RUN-084 posterior rank-1 soft-min parameterisation, copied into `config.yaml` as frozen first-stage parameters.

## Configuration

- `config.yaml` defines the frozen first stage, ABC priors for the fitted third layer, graphlet split policy, component guardrails, and quick-run overrides.

## Execution

- `run.sh` executes `scripts/phot_nested_softmin_ensemble_abc.py`.
- Use `labproj submit EXP-039 --target lab-slurm` for queued execution.

## Expected outputs

- nested_softmin_summary_json
- posterior_particles_csv
- heldout_model_distance_csv
- local_ruggedness_component_csv
- nested_control_comparison_csv

## Analysis plan

- Compare calibration-selected held-out graphlet distance for each nested family against the locked first-stage model and non-nested fitted controls.
- Inspect component distances to identify whether any global improvement is driven by the local ruggedness terms or only by smoother graphlet statistics.
- Inspect pairwise deltas against no-module, random-module, shuffled-module, and uncoupled-module families to test empirical-module specificity.
- Inspect module-energy and base-correlation outputs to diagnose whether any third-layer improvement is diffuse or localised.
- Create ANA and RES records after execution and update HYP-008/CLAIM-003 only if the held-out control comparison supports the inference.

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

- `RUN-098` completed with inferred Slurm exit code `0` and produced 30
  collected files under
  `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098`.
- `ANA-039` and `RES-039` record the analysis and result. The
  calibration-selected nested empirical-module family failed the primary
  held-out full-objective comparison (`0.543180`) against fitted current
  (`0.521007`), smooth core (`0.526819`), locked first-stage (`0.529590`),
  random-module (`0.482630`), and shuffled-module (`0.510961`).
- The result weakens `HYP-008` as currently operationalized and supports only a
  narrow generic-capacity interpretation for the nested posterior/negative
  controls.
