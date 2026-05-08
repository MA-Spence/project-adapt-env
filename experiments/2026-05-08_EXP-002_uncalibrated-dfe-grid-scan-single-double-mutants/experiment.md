# EXP-002: Uncalibrated DFE Grid Scan for Single and Double Mutants

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001

## Question

Without fitting to any empirical DMS landscape, what single-mutant and double-mutant DFEs does the raw Adapt-Env model family generate across a preregistered parameter grid, and does that grid span qualitatively protein-like local mutational regimes?

## Pre-experiment prediction

Across the uncalibrated grid, lower stability margins and higher epistasis strengths should shift one- and two-step DFEs toward more deleterious and lethal effects, double-mutant DFEs should be broader and more deleterious than single-mutant DFEs, and epistasis-sensitive settings should produce larger deviations from additive expectation for double mutants.

## Rationale

After RES-001 weakened HYP-001 under empirical calibration, the next rigorous question is whether the uncalibrated model family can generate a plausible envelope of one- and two-step DFEs at all. This does not directly test empirical recovery, but it is a necessary scope check: if the raw family cannot span realistic local DFE regimes before fitting, calibration failure is more likely a model-family limitation than only an objective-search problem.

## Experimental design

- Template: generic
- Sweep `stability_margin`, `functional_sigma_base`, `n_functional_dims`, and `epistasis_strength` over a preregistered uncalibrated grid.
- For each grid point and random seed, build an unconditioned Adapt-Env landscape around the default reference sequence.
- Exhaustively evaluate all single substitutions around the reference and summarize the one-step DFE.
- Sample a large panel of double substitutions around the same reference, summarize the two-step DFE, and compare each double-mutant effect to the additive expectation from the two constituent singles.
- Aggregate the resulting DFE and epistasis summaries by grid point and by parameter axis.

## Inputs

- `external/Adapt-Env`
- `config.yaml`

## Configuration

- `config.yaml` records the base landscape configuration, the parameter grid, and the single/double mutant sampling regime.
- This experiment is intentionally uncalibrated against empirical DMS data and uses no assay-specific fitting.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) a per-setting table of single- and double-mutant DFE summary statistics
- 2) aggregate summaries across parameter axes
- 3) double-mutant epistasis summaries relative to additive expectation
- and 4) a Slurm-executable experiment scaffold.

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
