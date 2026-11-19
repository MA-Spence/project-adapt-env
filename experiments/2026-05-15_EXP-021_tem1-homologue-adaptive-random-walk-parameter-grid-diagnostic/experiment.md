# EXP-021: TEM-1 homologue adaptive and random walk parameter-grid diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005

## Question

Do TEM-1-conditioned stability-margin landscapes reproduce literature-supported adaptive-walk and random-walk behaviour when starts are real TEM-1 homologues ranked by consensus distance, and do stability-margin and epistasis-strength knobs move those behaviours directionally over realistic parameter ranges?

## Pre-experiment prediction

If the EXP-020 failure was caused by peak-fraction stratification and parameter-specific stability-threshold behaviour rather than a generic inability to produce biological walk dynamics, then TEM-1 homologue starts farther from the alignment consensus should show lower stability margins, more beneficial opportunities, larger adaptive gains, less random-walk autocorrelation, and stronger deleterious/negative epistatic constraints; increasing stability margin should increase neutrality, while increasing epistasis strength should increase pairwise epistasis and ruggedness diagnostics.

## Rationale

EXP-020 mixed adaptive-walk support with failed peak-proximal epistasis tests. Literature on TEM-1 and stability-threshold proteins suggests the relevant axis for this model is robustness or consensus/stability distance, not stored peak fraction. This experiment repeats the EXP-020 walk/neighbourhood design using a real TEM-1 reference, homologues from the MMSeqs alignment used for conditioning, random walks for ruggedness/autocorrelation diagnostics, and a realistic stability-margin and epistasis-strength parameter grid to diagnose HYP-004 and HYP-005 together.

## Experimental design

- Template: generic
- Use the reviewed TEM-1 beta-lactamase sequence `P62593` as the explicit
  AdaptEnv reference sequence.
- Build the alignment-conditioned public stability-margin AdaptEnv landscape
  from the MMSeqs/ColabFold MSA fetched for that TEM-1 query sequence.
- Rank unique homologues from the same MSA by an alignment-derived
  consensus-bitscore distance: the consensus sequence receives the highest
  expected stability score under the stability-only model, and lower-scoring
  homologues are treated as farther from the implied stability optimum.
- Select starts from near-, intermediate-, and far-consensus homologue strata
  rather than sampling artificial random backgrounds or using stored peak
  fraction.
- Repeat the EXP-020 local-neighbourhood and adaptive-walk design from those
  homologue starts, but evaluate trends against consensus distance rather than
  stored peak distance or fraction of peak fitness.
- Add unconditional random walks from the same starts to measure fitness
  autocorrelation, step-effect variance, beneficial/neutral/deleterious random
  step fractions, and correlation-length proxies.
- Sweep a realistic grid of `stability_margin` and `epistasis_strength` values
  to test whether the intended mechanistic knobs move neutrality, deleterious
  burden, pairwise epistasis, and random-walk ruggedness diagnostics in
  interpretable directions.

## Literature baseline

- TEM-1 and other marginally stable proteins show a strong connection between
  stability, robustness, and the distribution of mutational effects. For this
  stability-only model, consensus distance is therefore the more defensible
  biological axis than distance to the stored AdaptEnv peak. Sources: Stiffler
  et al. 2015, `Cell`, doi: `10.1016/j.cell.2015.01.035`; Bershtein et al.
  2006, `Nature`, doi: `10.1038/nature04806`; Bloom et al. 2006, `PNAS`, doi:
  `10.1073/pnas.0510098103`.
- Empirical and theoretical adaptation literature predicts diminishing returns
  and depleted beneficial opportunities as backgrounds approach high-fitness
  regions. Sources: Wiser et al. 2013, `Science`, doi:
  `10.1126/science.1243357`; Kryazhimskiy et al. 2014, `Science`, doi:
  `10.1126/science.1250939`; Martin and Lenormand 2008, `Genetics`, doi:
  `10.1534/genetics.108.087122`.
- Random walks are a standard way to diagnose ruggedness through fitness
  autocorrelation, correlation length, and step-effect variance. Stronger
  epistasis should reduce autocorrelation and increase ruggedness proxies.
  Sources: Kauffman and Levin 1987, `Journal of Theoretical Biology`, doi:
  `10.1016/S0022-5193(87)80029-2`; Stadler 1996, `Journal of Mathematical
  Chemistry`, doi: `10.1007/BF01165157`.

## Inputs

- `external/Adapt-Env`
- `scripts/tem1_homologue_walk_parameter_grid.py`
- `config.yaml`
- TEM-1 beta-lactamase sequence `P62593`
- MMSeqs/ColabFold alignment fetched from the TEM-1 query sequence during run
  setup and exported as a run artifact

## Configuration

- `config.yaml` records the TEM-1 sequence, MMSeqs alignment settings, public
  stability-margin AdaptEnv core, homologue start strata, walk settings,
  epistasis sampling, and the stability-margin / epistasis-strength parameter
  grid.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- `tem1_mmseqs_alignment.fasta`: the MSA used for conditional stability.
- `alignment_profile_summary.json`: consensus, conservation, and frequency
  summary for the alignment.
- `homologue_starts.csv`: selected TEM-1 homologue starts and their
  consensus-distance diagnostics.
- `background_neighborhoods.csv`: local DFE, stability, and epistasis summary
  at each start/background.
- `walk_trajectories.csv`: per-step adaptive and random-walk trajectories.
- `walk_summaries.csv` and `walk_summaries.json`: per-walk summaries.
- `epistasis_pairs.csv`: sampled pairwise double-mutant epistasis rows.
- `parameter_grid_summary.csv`: parameter-set-level diagnostic means.
- `literature_baselines.json`: machine-readable baseline expectations.
- `summary.json`: aggregate validation tests and metadata.

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
