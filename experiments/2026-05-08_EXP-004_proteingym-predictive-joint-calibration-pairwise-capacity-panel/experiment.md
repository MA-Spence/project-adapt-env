# EXP-004: ProteinGym Predictive Joint Calibration And Pairwise Capacity Panel

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001
- Previous experiments: EXP-001, EXP-003
- Motivating results: RES-001, RES-002, RES-003

## Question

Does HYP-001 hold more cleanly on the fixed EXP-003 six-assay ProteinGym stability panel when calibration is changed from a staged summary-matching fit to a predictive joint fit and when empirical pairwise epistasis is represented with less aggressive compression?

## Pre-experiment prediction

Relative to the EXP-003 baseline, a predictive joint objective should improve single- and double-mutant holdout recovery on the shared raw-score branch, and relaxing empirical pairwise compression should further improve double-mutant and epistasis-prediction metrics while reducing collapse to zero epistasis strength. If calibration was the main bottleneck, the improved shared raw branch should approach the improved per-assay ceiling more closely than in EXP-003.

## Rationale

RES-002 showed that the Adapt-Env family spans broad local DFE envelopes, while RES-003 showed that coordinate-aware preprocessing alone does not rescue empirical calibration. This experiment is run to isolate two remaining calibration bottlenecks: the staged summary-first objective and the strongly compressed empirical pairwise epistasis representation. It fixes the panel to the exact six multi-mutant Tsuboyama stability assays used in EXP-003, then compares baseline and improved calibration branches on the same raw and latent assay constructions.

## Experimental design

- Reuse the exact six-assay Tsuboyama stability panel from EXP-003 so calibration changes can be compared on a fixed empirical target.
- Fit assay-specific `mavenn` models on raw ProteinGym scores as in EXP-003, then construct raw-score and latent-scale empirical landscapes from the same observed single and double mutants.
- Compare five preregistered branches:
  1. `baseline_shared_raw`: the EXP-003 shared raw calibration
  2. `predictive_shared_raw`: improved predictive-joint objective only
  3. `predictive_richpair_shared_raw`: improved objective plus less-compressive empirical pairwise representation
  4. `predictive_richpair_shared_latent`: same improved settings on the latent target
  5. `predictive_richpair_per_assay_raw`: improved raw per-assay fits as a ceiling comparison
- Use pair-based double-mutant holdout in the improved branches to evaluate whether richer pairwise representations generalize beyond the exact observed double variants used to construct them.

## Inputs

- ProteinGym substitution reference metadata
- ProteinGym substitutions parquet benchmark
- the six-assay EXP-003 ProteinGym stability subset defined explicitly in `config.yaml`
- real MMseqs alignments stored under `data/interim/proteingym_mmseqs_alignments`
- the branch-aware calibration runner `scripts/proteingym_predictive_joint_calibration.py`

## Configuration

- `config.yaml` records the fixed assay IDs, the retained EXP-003 baseline branch, and the improved predictive-joint and rich-pairwise branches.
- The improved branches enable the new Adapt-Env options for `functional_core_objective`, `functional_fit_strategy`, and `functional_joint_rounds`, and they relax empirical pairwise compression through `empirical_pairwise_observed_prior_blend`, `empirical_pairwise_coverage_shrink_power`, and `empirical_pairwise_normalization`.

## Execution

- Use `labproj submit` for RUN generation.

## Expected outputs

- 1) a fixed six-assay ProteinGym stability panel matching EXP-003
- 2) per-assay `mavenn` diagnostics
- 3) shared-branch calibration summaries in `branch_validations.csv`
- 4) per-assay improved-raw fits in `per_assay_branch_fits.csv`
- 5) a structured `summary.json`
- and 6) a Slurm-executable run scaffold

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

- This experiment is run because RES-003 reduced the coordinate-system confound but still left the central calibration failure unresolved.
