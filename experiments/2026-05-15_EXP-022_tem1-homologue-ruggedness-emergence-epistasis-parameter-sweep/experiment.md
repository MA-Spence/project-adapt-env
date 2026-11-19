# EXP-022: TEM-1 homologue ruggedness emergence epistasis parameter sweep

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005

## Question

At what epistatic density, strength, rank, covariation bias, and higher-order term settings does the TEM-1 homologue stability-margin landscape stop behaving as a nearly smooth additive landscape and begin showing measurable ruggedness in pairwise epistasis and random-walk diagnostics?

## Pre-experiment prediction

If the RUN-057 smoothness arose because epistasis terms were too sparse and weak relative to the additive stability layer, then increasing pairwise density and strength, using higher-rank pairwise matrices, enriching placement toward covarying MSA positions, and adding sparse higher-order exact-match terms should increase nonzero and strong pairwise epistasis, increase random-walk step variance, reduce random-walk autocorrelation/correlation length, and potentially increase sign or diminishing-returns epistasis. If the landscape remains smooth, then the current stability-margin readout is damping or hiding epistatic ruggedness even under stress settings.

## Rationale

RES-021 supports the TEM-1 homologue consensus-distance adaptive-walk diagnostic but shows that RUN-057 used epistasis settings that left the random walks almost perfectly autocorrelated and sampled sign epistasis absent. This experiment repeats the EXP-021 TEM-1 design while sweeping stronger, denser, higher-rank, and covariation-enriched epistasis regimes to isolate where ruggedness emerges and to test HYP-005 knob directionality more directly.

## Experimental design

- Template: generic
- Configure the experiment in config.yaml.
- Run via run.sh or labproj submit.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- tem1_mmseqs_alignment.fasta
- alignment_profile_summary.json
- homologue_starts.csv
- background_neighborhoods.csv
- walk_trajectories.csv
- walk_summaries.csv
- walk_summaries.json
- epistasis_pairs.csv
- parameter_grid_summary.csv
- literature_baselines.json
- summary.json

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
