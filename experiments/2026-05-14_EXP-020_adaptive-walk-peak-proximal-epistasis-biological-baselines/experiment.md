# EXP-020: Adaptive Walk And Peak-Proximal Epistasis Biological Baselines

## Status

- planned

## Scientific lineage

- Hypotheses: HYP-004

## Question

Do current AdaptEnv stability-margin landscapes reproduce literature-supported evolutionary-landscape behaviours: adaptive walks with decelerating gains, declining beneficial-mutation availability near local peaks, and stronger or more constraining epistasis as backgrounds approach a peak?

## Pre-experiment prediction

Across seeds and conditioned or unconditioned synthetic family backgrounds, adaptive walks should gain fitness early and then decelerate; starting closer to the stored local peak should reduce beneficial fraction and mean beneficial effect size; peak-proximal neighborhoods should show larger or more sign-constraining pairwise epistasis than low-fitness neighborhoods.

## Rationale

HYP-004 is a biological-realism test, not a machine-learning benchmark. The experiment compares synthetic walk and epistasis behaviour to qualitative expectations from empirical fitness-landscape, microbial-adaptation, and Fisher-geometric-model literature before any optimizer ranking is considered.

## Literature baseline

- Diminishing-returns epistasis is a reproducible empirical pattern in
  microbial adaptation: fitter backgrounds tend to gain less from the same or
  comparable beneficial mutations than less-fit backgrounds. The preregistered
  synthetic expectation is therefore a negative relationship between starting
  fraction of peak fitness and adaptive-walk total gain, beneficial fraction,
  and mean available beneficial effect. Sources: Kryazhimskiy et al. 2014,
  `Science`, doi: `10.1126/science.1250939`; Wiser et al. 2013, `Science`,
  doi: `10.1126/science.1243357`.
- Fisher-geometric-model theory predicts that the DFE and fixed beneficial
  effects depend on distance to an optimum. The preregistered expectation is
  depletion of beneficial effects near the relevant local peak rather than
  unlimited improvement. Source: Martin and Lenormand 2008, `Genetics`, doi:
  `10.1534/genetics.108.087122`.
- Empirical fitness landscapes show that sign epistasis and ruggedness can
  restrict accessible evolutionary paths. The preregistered expectation is
  stronger or more constraining pairwise epistasis in peak-proximal
  neighborhoods than in low-fitness neighborhoods. Sources: de Visser and Krug
  2014, `Nature Reviews Genetics`, doi: `10.1038/nrg3744`; Weinreich et al.
  2006, `Science`, doi: `10.1126/science.1123539`.

## Experimental design

- Template: generic
- Build current public-core AdaptEnv landscapes using the single
  `stability_margin` scalar readout, pairwise stability-targeted epistasis, and
  no legacy FGM readout or optimizer benchmark.
- Test both unconditioned landscapes and synthetic-alignment-conditioned
  family-like landscapes across independent seeds.
- For each landscape, sample starting backgrounds near low, mid, high, and
  near-peak fractions of the stored local peak fitness.
- From each background, run greedy and random-beneficial adaptive walks and
  record per-step fitness, cumulative gain, distance to the stored peak,
  available beneficial fraction, and available beneficial effect sizes.
- At the same backgrounds, exhaustively scan one-step neighborhoods for DFE
  summaries and sample double mutants for pairwise epistasis, sign epistasis,
  reciprocal sign epistasis, negative epistasis, and diminishing-returns
  epistasis.
- Validate against literature-derived qualitative expectations using
  preregistered rank trends with fraction of peak fitness rather than
  benchmark scores.

## Inputs

- `external/Adapt-Env`
- `scripts/adaptive_walk_peak_epistasis_baselines.py`
- `config.yaml`

## Configuration

- `config.yaml` records the AdaptEnv core configuration, conditioning modes,
  seeds, target fractions of the stored peak, walk settings, local epistasis
  sampling, and validation trend thresholds.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- `walk_trajectories.csv`: one row per adaptive-walk step.
- `background_neighborhoods.csv`: one row per sampled starting background with
  local DFE and epistasis summaries.
- `epistasis_pairs.csv`: sampled pairwise double-mutant epistasis rows.
- `walk_summaries.json`: one adaptive-walk summary per background and mode.
- `literature_baselines.json`: machine-readable baseline expectations.
- `summary.json`: validation tests and aggregate pass/fail status.

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
