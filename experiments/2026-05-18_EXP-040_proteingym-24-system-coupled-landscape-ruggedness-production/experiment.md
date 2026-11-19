# EXP-040: ProteinGym 24-system coupled landscape ruggedness production

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Across a diverse 24-system ProteinGym multi-mutant panel, do base thermodynamic residual-epistasis, aligned two-way coupled, and aligned three-way coupled landscape fits reproduce empirical predictive and topological ruggedness behaviour on held-out variants and graphlets?

## Pre-experiment prediction

The base thermodynamic residual model should recover coarse single/double-mutant predictive structure in some systems but underfit high-order and local ruggedness summaries; aligned coupled capacity should improve bulk graphlet, DFE, walk, and spectral statistics while preserving correlation to the base function; three-way coupling may add capacity but should not be interpreted as residue-specific mechanism.

## Rationale

CLAIM-001 and CLAIM-003 support system-specific thermodynamic calibration and coupled auxiliary capacity, but almost all coupled-landscape evidence is PHOT-only. A production multi-system experiment is needed to certify whether the same model hierarchy improves empirical-like landscape topology across activity, binding, expression, organismal-fitness, and stability assays with real combinatorial mutation coverage.

## Experimental design

- Select 24 ProteinGym substitution assays with multiple-mutant coverage, capped
  at 10 Tsuboyama stability assays and spanning activity, binding, expression,
  organismal-fitness, and stability readouts.
- Fit each system independently with the current thermodynamic
  stability-margin plus stability-routed residual-epistasis ABC model. No
  pooled cross-assay objective is used.
- On each fitted base model, evaluate empirical and model dense Hamming
  graphlets for DFE, SE/RSE, magnitude epistasis, local maxima, accessibility,
  random-walk autocorrelation, local Dirichlet-energy concentration, WHT, and
  GFT summaries.
- Fit an alignment-routed two-way coupled residual layer and then an
  alignment-routed three-way layer against calibration graphlets while
  regularising the coupled function to remain correlated with the base
  thermodynamic function.
- Evaluate base, two-way coupled, and three-way coupled models on held-out
  graphlets plus held-out single and double mutants. Record random-walk and
  adaptive-walk step diagnostics where complete empirical graphlets permit the
  same measurements.
- Interpret coupled terms only as bulk ruggedness capacity routed through an
  alignment-legible scaffold, not as recovered residue-specific interaction
  terms.

## Inputs

- ProteinGym v1.3 substitution reference table and substitutions parquet. The
  runner downloads these from the configured ProteinGym URLs if they are absent
  in the staged run.
- MMseqs alignments under `data/interim/proteingym_mmseqs_alignments`; missing
  alignments are prepared by the existing ProteinGym panel helper during each
  base fit.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- resolved_config.json
- selected_systems.csv
- target_selection_audit.json
- system_fit_summary.csv
- system_variant_validation.csv
- system_graphlet_metrics.csv
- system_model_distances.csv
- system_walk_step_metrics.csv
- coupled_parameter_summary.csv
- coupling_function_stats.csv
- coupled_particles.csv
- system_walsh_terms.csv
- system_gft_coefficients.csv

## Analysis plan

- Create an ANA record after execution.

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
