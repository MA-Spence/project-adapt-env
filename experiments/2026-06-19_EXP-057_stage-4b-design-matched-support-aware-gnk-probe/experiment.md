# EXP-057: Stage 4B design-matched support-aware GNK probe

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009

## Question

Does the public additive thermodynamic core plus GNK texture show protein-scale monotone knob control when measured by design-matched, support-aware graphlets on the instantiated GNK supports?

## Pre-experiment prediction

Support-aware, design-matched graphlets will restore observable full-public versus additive-core separation as epi_weight increases, with monotone increases in high-order/specific epistasis and robust local maxima density and monotone decreases in peak-targeted reachability/autocorrelation.

## Rationale

This is the single Section 3 fix from scratch/RES-054_assessment_and_path_to_finishing_component1.md: rerun Stage 4B with support-aware sites, matched empirical dimension/alphabet, regime-stratified anchors, plateau-robust topography, and a hard core-vs-GNK separation gate before any overlap claim.

## Experimental design

This is the single Section 3 repair experiment from
`scratch/RES-054_assessment_and_path_to_finishing_component1.md`; it is not a
new overlap, joint-realism, held-out, or benchmark-claim experiment.

For each Tier-1 empirical design, the runner uses the observed combinatorial
template as the model measurement template:

- model graphlet dimension equals the empirical design dimension;
- graphlet alphabet cardinalities match the empirical alphabet, including
  20-letter Potts graphlets for GB1/TrpB/ParD and binary graphlets for
  CR9114/CR6261;
- varied sites are sampled from instantiated GNK supports or support unions,
  not from positions `0..dim-1`;
- anchors are stratified across unfolded, fold-transition, and plateau
  effective-stability margins;
- topography uses one shared tolerance, strict-greater local maxima density,
  and peak-targeted reachability.

The public observable is compared against additive-core-only and bare-GNK
controls under structured, shuffled, and random GNK support modes. The hard
gate is whether `full_public_observable` separates from
`additive_core_only_observable` by preregistered effects in high-order
fraction, sign epistasis, or local maxima density as `epi_weight` changes.

If monotonicity or the core-vs-GNK separation gate fails, the correct inference
is a Stage 4B measurement/scale-balance failure. The run must stop there; it
does not license 4C-4E reruns.

## Inputs

- `scratch/RES-054_assessment_and_path_to_finishing_component1.md`
- `scripts/design_matched_support_aware_probe_exp057.py`
- `scripts/ensemble_landscape_validation_stage4.py`
- `data/processed/ensemble-landscape-validation-datasets.dvc`
- `external/Adapt-Env/adaptenv/landscape.py`

## Configuration

- `config.yaml`

## Execution

- Dry-run `RUN-140` packaged the job.
- Live `RUN-141` attempt 1 failed by Slurm walltime as job `1948`.
- Live `RUN-141` resumed attempt 3 completed on `lab-slurm` as job `1955`.

## Expected outputs

- preflight.json
- empirical_design_templates.csv
- empirical_template_descriptor_metrics.csv
- generator_parameter_grid.csv
- support_catalog.csv
- support_aware_graphlet_audit.csv
- anchor_regime_audit.csv
- generator_descriptor_metrics.csv
- generator_knob_monotonicity.csv
- core_texture_ablation_controls.csv
- core_gnk_separation_gates.csv
- checkpoint_resume_summary.json
- summary.json

## Analysis plan

- Completed as `ANA-055`.

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

- `RUN-141` completed as an execution, but `RES-055` records a negative
  scientific result: only `318/1620` monotonicity checks passed and only
  `34/1080` positive-epistasis core-vs-GNK separation rows passed.
- Per the registered stop rule, downstream overlap, joint-realism, and held-out
  reruns remain blocked for the current public implementation.
