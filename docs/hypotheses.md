# Hypotheses

## HYP-001: Local-Protein-Landscape-Statistics-Without-Exact-Reconstruction

- Status: superseded
- Superseded by: CLAIM-001
- Owner: matthew-spence
- Aims: AIM-001

A hybrid synthetic landscape built from a stability threshold, a low-dimensional
functional layer, and sparse structured epistasis can match empirical local
summary statistics from real protein landscapes without attempting to
reconstruct any one assay exactly.

Predictions:
- A shared synthetic parameter regime should reproduce DFE shape, neutrality,
  beneficial fraction, and fitness-dependent mutational-effect trends across a
  panel of empirical assays at the level of summary statistics.
- The match should hold without assay-specific variant-level fitting.

Why this matters:
- The benchmark goal is not exact replay of one DMS experiment, but a
  synthetic family that still looks protein-like to learning algorithms and
  evolutionary processes.

Test:
- Choose a panel of empirical assays spanning stability, binding, and activity.
- Compare synthetic landscapes to empirical assays using summary statistics
  rather than variant-by-variant fit.
- Fit broad parameter ranges, then evaluate on held-out assays or held-out
  proteins.

Primary readouts:
- single-mutant DFE shape
- neutral, deleterious, and beneficial fractions
- DFE skewness and tail weight
- fitness-dependent mutation effects
- local epistasis prevalence and sign structure

Falsification criteria:
- no shared parameter regime can reproduce the empirical summary-statistic
  envelope across multiple assay families
- realism only appears after assay-specific tuning

## HYP-002: Benchmark-Certified-Peak-Mode-Preserves-Local-Realism

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

A benchmark-facing landscape mode that certifies or tightly controls the top of
the landscape will make optimization metrics interpretable while preserving the
local DFE and epistasis properties needed for realism.

Predictions:
- A certified-top mode should reduce peak overshoot and reference-to-peak
  ambiguity.
- Local DFE, global epistasis, and neutrality metrics should remain within
  pre-specified tolerances of the biological mode.

Why this matters:
- Fraction-of-peak and regret-style metrics are not rigorous if the stored peak
  is only a local-search artifact.

Test:
- Compare the current public reference-neutral mode against a
  benchmark-certified mode such as strict peak locking or peak-preserving
  pairwise construction.
- Evaluate both modes on the same seeds and alignments.
- Require local realism metrics to remain within pre-specified tolerances.

Primary readouts:
- certified top-sequence recovery under exhaustive search on small landscapes
- overshoot rate in probe searches
- reference-to-peak distance
- DFE, global epistasis, and neutrality changes between modes

Falsification criteria:
- benchmark certification materially distorts local mutational statistics
- algorithm rankings change only because the certified mode becomes
  artificially easy

## HYP-003: Real-MSA-Conditioning-Improves-Realism-Over-Synthetic-Conditioning

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

Conditioning the landscape on real family MSAs will improve
conservation-sensitivity structure and pairwise epistasis placement more than
conditioning on internally generated synthetic alignments.

Predictions:
- Real MSA conditioning should outperform synthetic-alignment conditioning on
  held-out conservation-effect correlation.
- Real MSA conditioning should outperform synthetic-alignment conditioning on
  covariation-epistasis enrichment and natural-sequence plausibility metrics.

Why this matters:
- Current internal tests mostly show that the model is consistent with its own
  synthetic priors, not that those priors match real protein families.

Test:
- Assemble a panel of real MSAs linked to ProteinGym or other mutational
  datasets.
- Compare three conditions: no alignment, synthetic alignment, real MSA.
- Evaluate only on summary statistics not used to set the parameters.

Primary readouts:
- conservation versus mean mutational effect correlation
- enrichment of strong epistasis among high-covariation site pairs
- viability of natural-like sampled sequences
- calibration of mutation effects for consensus, natural, and off-manifold
  sequences

Falsification criteria:
- real MSA conditioning offers no improvement over synthetic conditioning
- alignment conditioning improves only trivial viability checks but not
  epistasis or sensitivity structure

## HYP-004: Benchmark-Regimes-Separate-Search-Strategies-Reproducibly

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

A scientifically useful synthetic benchmark should produce reproducible
performance separation between qualitatively different optimization strategies
across budgets, start states, and noise levels.

Predictions:
- Across seeds and budgets, strategy rankings should remain stable enough that
  random search, local mutagenesis, recombination, and uncertainty-aware active
  learning do not collapse into one indistinguishable performance class.
- The separation should persist across both noiseless and noisy evaluation
  regimes.

Why this matters:
- A landscape can look biologically plausible yet still be a poor benchmark if
  random search, exploitative search, and uncertainty-aware search are not
  reliably distinguished.

Test:
- Compare random, local mutagenesis, recombination, Bayesian optimization, and
  active-learning baselines across many seeds.
- Repeat across natural-evolution and directed-evolution presets.
- Evaluate both exact and noisy observations.

Primary readouts:
- area under the best-so-far curve versus evaluation budget
- probability that one strategy beats another across seeds
- rank stability under resampling
- sensitivity of rankings to start-library distance from the reference and peak

Falsification criteria:
- rankings are dominated by seed luck or observation noise
- regime labels such as high-throughput versus low-throughput do not yield
  stable difficulty differences

## HYP-005: Mechanistic-Knobs-Produce-Monotonic-Directional-Effects

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

The model knobs intended to represent biological concepts should move
validation metrics in the expected direction over broad parameter sweeps.

Predictions:
- Increasing stability margin should change neutral and lethal fractions in the
  expected direction.
- Functional dimensionality, epistasis strength, covariation bias, and noise
  settings should each induce directional changes in their target validation
  metrics with effect sizes larger than seed-to-seed variation.

Why this matters:
- If the knobs are not directional, the model is hard to interpret
  scientifically and hard to use for controlled benchmark generation.

Test:
- Run preregistered sweeps over stability margin, functional dimensionality,
  epistasis strength, noise amplitude, and covariation bias.
- Estimate effect sizes and monotonicity with confidence intervals over many
  seeds.

Primary readouts:
- stability margin versus neutral fraction and lethal fraction
- functional dimensionality versus beneficial fraction and adaptive
  accessibility
- epistasis strength versus local maxima count and accessible paths
- covariation bias versus enrichment of strong interactions at coupled sites

Falsification criteria:
- the expected directional effect is absent or reverses across nearby settings
- landscape behavior is dominated by seed-specific stochasticity rather than
  parameter meaning

## HYP-006: Observation-Noise-Is-Realistic-And-Does-Not-Destabilize-Conclusions

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

A useful active-learning benchmark needs an observation model that perturbs
ranking and uncertainty in a realistic way, but does not erase meaningful
differences between strategies.

Predictions:
- Empirically grounded noise models should change local rankings and search
  difficulty.
- Strategy ordering and top-k recovery should remain statistically stable over
  practical noise ranges.

Why this matters:
- The benchmark is meant for iterative decision-making pipelines, not only a
  noiseless oracle optimization task.

Test:
- Evaluate fixed, heteroscedastic, and assay-inspired noise models.
- Compare inferred and true top-k recovery, regret, and ranking stability.
- Calibrate noise scales against replicate variability from empirical assays
  where available.

Primary readouts:
- top-k recovery under noisy measurement
- Kendall or Spearman agreement between noisy and true ranking
- rank reversals among close variants
- strategy ordering stability across noise regimes

Falsification criteria:
- small noise changes cause large and erratic benchmark reordering
- the default noise model is not empirically defensible or is too weak to
  change search behavior

## HYP-007: Multiple-Latent-Molecular-Phenotypes-Improve-Empirical-Recovery

- Status: superseded
- Superseded by: CLAIM-002
- Owner: matthew-spence
- Aims: AIM-001

Replacing the current single functional latent layer with multiple explicit
latent molecular phenotypes, such as folding stability together with binding or
abundance, linked to assay readouts through nonlinear thermodynamic or
measurement maps, will improve empirical recovery and realism.

Predictions:
- On matched empirical assays, multi-phenotype models should improve held-out
  single-mutant ranking, double-mutant ranking, and epistasis-prediction
  metrics relative to the current single-latent model.
- The improvement should be largest on assays where the present fitter
  collapses to zero epistasis, near-peak reference artifacts, or weak
  double-mutant recovery.
- When orthogonal measurements such as abundance plus binding or stability plus
  activity are available, a shared latent decomposition should transfer across
  readouts better than a single scalar latent score.

Why this matters:
- `RES-005` weakens the idea that the remaining failure is mainly optimizer
  collapse. Otwinowski-style folding and binding decompositions, together with
  the Lehner-group `deepPCA`, `ddPCA`, and `MoCHI` results, suggest that much
  observed epistasis may arise from nonlinear readouts of a small number of
  additive molecular traits rather than from one generic phenotype plus ad hoc
  ruggedness.

Test:
- Build paired benchmark panels with orthogonal readouts where possible, such
  as abundance plus binding or stability plus activity, and compare the current
  model against explicit two-trait and multi-trait latent variants on held-out
  single and double mutants.
- Evaluate both assay-specific fits and limited-sharing hierarchical fits to
  determine whether the latent traits transfer across related assays or
  proteins.
- Require improvement over the current model on predictive metrics and on
  mechanistic diagnostics such as reduced reference-to-peak artifacts and
  non-degenerate latent parameter posteriors.

Primary readouts:
- single-mutant holdout Spearman or Kendall
- double-mutant holdout Spearman or Kendall
- epistasis-prediction Spearman and KS
- cross-readout transfer accuracy when fitting one shared latent decomposition
- fitted latent-trait interpretability, such as shared folding effects across
  binding readouts

Falsification criteria:
- multi-phenotype models do not materially outperform the current single-latent
  model on matched assays
- any gains come only from much larger unconstrained parameterization rather
  than from transferable latent structure
- latent decompositions fail to align across orthogonal readouts or collapse
  back to effectively one trait
