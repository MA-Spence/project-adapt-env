# EXP-059: Saturated-gate GNK mechanism calibration and held-out confirmation

## Status

- planned

## Scientific lineage

- Aim: AIM-001
- Hypothesis tested directly: HYP-005
- Previous experiment: EXP-058
- Motivating results: RES-049, RES-055, RES-057

HYP-004 and HYP-009 remain downstream. This experiment contains no empirical-realism
comparison and no agents, so it cannot directly test either hypothesis. It also does
not test functional dimensionality or reopen superseded HYP-007.

## Question

Can the private legacy two-channel product route, when the stability gate is
demonstrably saturated/open and the function-capacity channel remains non-saturated,
transmit a paired GNK-weight
intervention into observable epistasis and local topology while retaining nonzero
navigability on held-out protein-scale subspaces?

## Pre-experiment prediction

At a component-qualified open-gate operating point, increasing GNK weight from the
paired weight-zero causal null to the frozen working weight will increase q-ary
non-additive variance, sign epistasis, and local-maxima density, while reducing but
not collapsing basin convergence in held-out seeds.

## Rationale and scope

RES-057 rejected the tested Fix C configuration because its low stability margins
nearly closed the multiplicative gate. It did not test the same architecture in an
open-gate regime. EXP-059 isolates that remaining mechanistic possibility without
weakening the Stage 4C/4D/4E hard stop.

The tested product/function-capacity implementation is a private compatibility route
marked for deletion from the public AdaptEnv API. A positive result would establish
only an internally controllable candidate mechanism for the GNK-weight portion of
HYP-005. It would not establish empirical realism, biological routing privilege,
optimizer performance, public-generator readiness, or a need for more latent
functional dimensions.

## Experimental design

One queued execution contains two sequential stages. Calibration selects an operating
point using component checks alone, then tests the frozen point. Confirmation runs
only if calibration passes its mechanistic gates. Calibration and confirmation seeds
are disjoint.

### Paired intervention and controls

- Sequence length is 48 and alphabet size is 20 throughout.
- GNK weights are `0.0`, `0.3`, `0.6`, and `1.0`; `0.0` is the causal null and
  `1.0` is the preregistered working weight.
- The base-design key is `(stage, dataset_id, seed, graphlet_index)`. It excludes
  operating point, architecture, routing mode, and GNK weight. Genotypes, the
  designed-peak background, variable sites, allele choices, raw GNK terms, sampled
  edges/pairs, and basin starts therefore remain paired across interventions. The
  operating-point sweep changes only stability margin and threshold steepness.
- For reduced-cardinality confirmation axes, allele sets contain the anchor allele
  followed by the highest-stability alternatives under the frozen weight-zero core.
  This deterministic viability construction is shared across treatments; q=20 axes
  still contain the complete alphabet. The confirmation estimand is therefore
  conditional on stability-qualified, support-aware template-shaped subspaces, not
  on generic or randomly sampled protein subspaces.
- The primary architecture is the two-channel product. The mixed score/stability
  architecture is a paired comparator. An unsaturated product with the stability
  gate deliberately closed is the direct mechanistic negative control.
- Structured routing is primary. Shuffled and random routing are equivalence
  diagnostics, not causal nulls and not pass/fail requirements for the mechanism.
- Fitness neutrality uses one fixed absolute tolerance, `1e-6`, across all paired
  conditions.

### Stage 1: exact embedded calibration

For each seed `701`–`704`, one graphlet is anchored on the frozen designed-peak
background and four sites vary over all 20 alleles. The resulting `20^4 = 160,000`
genotypes and their complete one-step
edge set are enumerated. This preserves the L=48 function-capacity scaling while
making the calibration subspace exact.

Operating points are checked in the fixed order recorded in `config.yaml`:
`m14_s2`, `m14_s4`, `m18_s2`, `m18_s4`, `m22_s2`, `m22_s4`. The first point for which
all four seeds pass the gate-occupancy, function-capacity, and product-reconstruction
preflight is selected. No epistasis or topology result participates in selection.
If no point qualifies, the run records a normal scientific failure and does not run
confirmation.

At the selected point, the genotype set, q-ary factorial non-additive variance,
and local maxima are exact. Sign epistasis and deterministic steepest-ascent basin
convergence use frozen treatment-independent sampled plaquettes and starts. All
shortest paths are enumerated conditional on preregistered sampled starts and
reported as an audit; that conditional result does not define a confirmation
acceptance band.

### Stage 2: held-out confirmation

The selected operating point is frozen without retuning and applied to seeds
`1701`–`1708`. Four empirical templates—GB1, DHFR, CR9114 H1, and ParD2—supply only
subspace dimensions and alphabet cardinalities; their measured fitness values are
not used. The independent inferential unit is the generator seed:

1. summarize graphlets by the median within `(template, seed)`;
2. summarize templates by the median within seed;
3. evaluate replication over the eight seeds.

A missing or preflight-rejected confirmation base counts as a failure rather than
being dropped or replaced.

## Preregistered checks

### Preflight

- Every node in each primary viable graphlet must have stability gate `>= 0.90`.
- At least 90% of sampled dead-state nodes must have gate `<= 0.10`.
- At every weight, function capacity on viable nodes must have SD `>= 0.05`, with
  at least 65% of values strictly between `0.05` and `0.95`. This floor was frozen
  after a component-only feasibility audit across calibration seeds; no epistasis or
  topology result was inspected.
- The unweighted GNK score must have SD `>= 1e-3` on each primary graphlet, ensuring
  the support-aware intervention is nonconstant before effect metrics are evaluated.
- Reconstructed `gate * function_capacity * fitness_scale` must agree with the
  product observable to absolute tolerance `1e-8`.
- The direct negative control is shifted to achieved effective margin `-2.0`, and
  no more than 10% of its graphlet nodes may have gate above `0.10`.

### G1: component viability

At every GNK weight, observables must remain within `[0, fitness_scale]` up to
`1e-8`; the dead-state 95th percentile must be at most 10% of scale; and viable-node
fitness must have median at least 30% and SD at least 5% of scale. A fixed-background
stability titration targets achieved margins `[-2,-1,-0.5,0,0.5,1,2,4]` rather
than offsets around the already-saturated operating point. It checks implementation
only: the slope over achieved margins `[-0.5,0.5]` is positive,
saturated slope is nonnegative, and the saturated-to-transition slope ratio is at
most 0.20. These checks are not evidence for biological diminishing returns.

### G2: GNK transmission

For both non-additive variance fraction and sign-epistasis fraction, the working-
weight minus weight-zero paired delta must have overall median `>= 0.02`, positive
median in every template, and positive seed-level median in at least 7/8 seeds. The
exact one-sided sign-test must have `p <= 0.0352`, and the seed bootstrap 95% lower
bound must be greater than zero. Zero deltas count as nonpositive; ties remain in the
`Binomial(n=8, p=0.5)` denominator. The lower bound is the fifth percentile of
10,000 median bootstraps over eight seed summaries using seed `59059`. Both metrics
must pass.

### G3: topology control with retained navigability

Across weights, median Spearman correlation must be at least 0.80 overall and 0.50
within every template for non-additive variance, sign epistasis, and local-maxima
density; at least 7/8 seed summaries must have the expected direction. Undefined
correlations caused by constant series fail. At the working weight, the overall seed
median local-maxima-density delta must be at least 0.005 and median local-maxima-count
increase at least one; every template median must be positive, and at least 7/8 seed
summaries must meet both floors. For basin convergence, at least 7/8 seed summaries
and every template median must lie in `[0.05, 0.80]` and fall by at least 0.05 relative
to weight zero.

### G4: routing equivalence diagnostic

Paired seed-level percentile-bootstrap 90% intervals (10,000 resamples, seed
`59060`) compare structured routing with shuffled and random routing
against margins of `0.02` for both epistasis fractions, `0.005` for local-maxima
density, and `0.05` for basin convergence. Equivalence requires the entire interval
to lie within the symmetric margin. G4 is descriptive: failure keeps the
HYP-009 prerequisite blocked but does not overturn the EXP-059 mechanism decision.

## Decision rule

- If calibration preflight has no qualifying operating point, stop before
  confirmation and report `NO_COMPONENT_QUALIFIED_POINT`.
- Calibration G1 requires all component checks at every weight in all four seeds.
  Calibration G2 requires, for both endpoint epistasis deltas, a median of at least
  0.02 and positive deltas in 4/4 seeds; no sign test or bootstrap is used at this
  selection-independent calibration stage. Calibration G3 requires expected trends
  in 4/4 seeds, median rho at least 0.80, endpoint local-maxima-density delta at
  least 0.005 and count increase at least one, and the same working-weight basin
  band/delta used in confirmation.
- If the selected calibration point fails any of G1, G2, or G3, stop before
  confirmation and report `SELECTED_COMPONENT_QUALIFIED_POINT_FAIL`. Because the
  point was selected by a component-only fixed order, other component-qualified
  points remain untested; this outcome does not reject the entire saturated-gate
  mechanism family.
- Confirmation passes only if every required primary two-channel preflight and all
  of G1, G2, and G3 pass under structured routing.
- A passing comparator shows that the product readout is not uniquely required; it
  does not invalidate a primary mechanism pass.
- Any primary pass is reported as `LEGACY_TWO_CHANNEL_MECHANISM_DEMONSTRATED`, never
  as public benchmark readiness.

## Execution and reproducibility

- Configuration: `config.yaml`
- Entrypoint: `run.sh`
- Runner: `scripts/two_channel_viability_benchmark_exp059.py`
- Default target: `lab-slurm` (8 CPU, 64 GB RAM, 72 hours)
- All stochastic design choices use stable recorded seeds and treatment-independent
  keys. The run emits the full design catalog, component metrics, descriptor metrics,
  gate tables, frozen operating point, routing diagnostics, and final decision.

## Expected outputs

- `preflight.json`
- `design_catalog.csv`
- `component_metrics.csv`
- `descriptor_metrics.csv`
- `enumerable_landscape_metrics.csv`
- `support_catalog.csv`
- `calibration_decision.csv`
- `frozen_operating_point.json`
- `confirmation_seed_metrics.csv`
- `viability_checks.csv`
- `gnk_transmission_checks.csv`
- `topography_checks.csv`
- `routing_equivalence.csv`
- `mechanism_decision.csv`
- `summary.json`

## Analysis plan

After execution, create an ANA record that evaluates the preregistered decision
without retuning thresholds or treating graphlets as independent replicates.

## Completion criteria

- [x] Configuration and protocol preregistered
- [x] Run script and environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [x] PROJECT_STATE.md updated after execution

## Post-experiment notes

- `RUN-145` completed on `lab-slurm` as Slurm job `2038` with exit code `0`.
- The preregistered `NO_COMPONENT_QUALIFIED_POINT` branch fired: no operating
  point passed component preflight in all four calibration seeds, so no point
  was selected and confirmation was not run.
- The 15 declared outputs and three runtime logs are collected under
  `data/processed/saturated-gate-gnk-mechanism-calibration-confirmation/RUN-145`
  and tracked by
  `data/processed/saturated-gate-gnk-mechanism-calibration-confirmation/RUN-145.dvc`.
- A formal analysis/result record and hypothesis update remain pending.
