# EXP-005: ProteinGym Shared Summary SMC-ABC Calibration Panel

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001

## Question

Does a project-local shared SMC-ABC summary calibration recover the fixed six-assay ProteinGym stability panel more faithfully than the current deterministic Adapt-Env fitters?

## Pre-experiment prediction

If the dominant bottleneck is calibration rather than model-family expressivity, then an SMC-ABC fit against a structured per-assay summary target should improve shared raw holdout behavior relative to EXP-004 and should recover known parameters in synthetic-truth tests without collapsing immediately to zero epistasis.

## Rationale

RES-002 kept the model family scientifically live, while RES-003 and RES-004 showed that coordinate-aware preprocessing and deterministic objective redesign still leave the inverse problem weak. This experiment is run to replace the current coarse staged/grid fitting path with a project-local likelihood-free Bayesian calibration on the same six-assay ProteinGym stability panel, with synthetic-truth recovery built in before the empirical test.

## Experimental design

- Reuse the exact six-assay Tsuboyama stability panel from EXP-003 and EXP-004 so the Bayesian calibration can be compared directly to the previous deterministic fitters on the same empirical target.
- Fit assay-specific `mavenn` models on raw ProteinGym scores exactly as in EXP-003 and EXP-004, then construct shared raw-score empirical landscapes for the Bayesian and deterministic calibrations.
- Run two deterministic controls on the same panel:
  1. `baseline_shared_raw`, matching the original EXP-003 shared raw fitter
  2. `predictive_richpair_shared_raw`, matching the strongest deterministic branch from EXP-004
- Build a structured per-assay summary target from observed single- and double-mutant statistics and estimate its covariance by bootstrap.
- Fit a shared local parameter posterior with project-local SMC-ABC rather than the Adapt-Env staged grid fitter.
- Require a synthetic-truth recovery step on the same panel scaffold before interpreting the empirical posterior as evidence about HYP-001.

## Inputs

- ProteinGym substitution reference metadata
- ProteinGym substitutions parquet benchmark
- the fixed six-assay EXP-003 / EXP-004 ProteinGym stability subset defined explicitly in `config.yaml`
- real MMseqs alignments stored under `data/interim/proteingym_mmseqs_alignments`
- the deterministic Adapt-Env fitters retained from EXP-004
- the project-local Bayesian calibration runner `scripts/proteingym_bayesian_summary_calibration.py`

## Configuration

- `config.yaml` records the fixed assay IDs, the retained deterministic control branches, the SMC-ABC summary target, the parameter priors, and the synthetic-truth recovery settings.
- The Bayesian fitter is intentionally project-local and lives under `src/project_adapt_env/` rather than in `external/Adapt-Env`, because this is an experimental calibration workflow rather than a change to the external package API.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) the fixed six-assay ProteinGym stability panel
- 2) per-assay mavenn diagnostics
- 3) deterministic branch comparison metrics
- 4) SMC-ABC posterior particles and round diagnostics
- 5) synthetic-truth recovery metrics
- and 6) a structured summary.json

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

- This experiment is run because RES-004 still left HYP-001 mainly bottlenecked by how the inverse problem is solved, not obviously by lack of local model-family capacity.
