# EXP-025: TEM-1 aggregate graphlet residual ruggedness generalization diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005
- Previous experiment: EXP-024
- Motivating result: RES-024

## Question

Can sparse contact-conditioned residual fields calibrated only to aggregate
local graphlet distributions produce held-out TEM-1 dense-Hamming graphlets
with nontrivial SE/RSE, multiple local maxima, shorter random-walk
autocorrelation, and empirical-looking Walsh/GFT spectra across adaptive walks
from very different starting states, without per-graphlet hard certification
or global peak conditioning?

## Pre-experiment prediction

A contact-conditioned residual ensemble calibrated on aggregate graphlet
statistics should improve held-out SE/RSE prevalence, local maxima counts,
accessibility reduction, random-walk autocorrelation, and Walsh/GFT spectral
envelopes relative to smooth-core and exact-match controls. If improvements
appear only in calibration graphlets, only in one adaptive-walk region, or only
before embedding in the thermodynamic core, the residual design should be
treated as overfit or region-specific rather than as a scalable biological
ruggedness solution.

## Rationale

RES-024 showed that DCA clique exact-match haplotype terms provide precise
high-order Walsh coefficient control but only weak topology: lower Walsh orders
remain dominant, RSE is weak, and local maxima remain few. That result supports
the knob-control component of HYP-005 but weakens HYP-004 as a full rugged
landscape claim.

The proposed alternative is deliberately less contrived than per-graphlet
certification. The experiment tests a sparse residual-field architecture whose
parameters are selected by matching aggregate distributions over calibration
graphlets. It then evaluates the selected parameter regimes on held-out dense
Hamming graphlets sampled from adaptive walks that begin in very different
regions of the TEM-1 conditioned landscape. This directly tests whether the
architecture generalizes as a procedural landscape-family mechanism rather than
overfitting a chosen local cube.

The experiment avoids conditioning on distance to a stored or inferred peak.
The current record already shows that stored peak geometry is not a reliable
biological coordinate. All stratification uses observable local quantities:
start class, walk step, background fitness, local DFE shape, local additive
slope, neutrality, and alignment/consensus distance.

## Experimental Design

- Use the reviewed TEM-1 beta-lactamase sequence `P62593` and the same
  MMSeqs/ColabFold alignment source used by EXP-021 through EXP-024.
- Build a shared smooth thermodynamic AdaptEnv core with old generic epistasis,
  exact-match higher-order terms, deterministic roughness, and noise disabled.
- Choose dramatically different starting states:
  - the reference or consensus-like TEM-1 background;
  - high-fitness close homologues;
  - marginal-stability homologues near the folding cliff;
  - distant low-fitness homologues;
  - viable off-consensus backgrounds sampled from alignment-conditioned
    mutations where available.
- From each start, run several adaptive walks under the smooth core. Record
  the full trajectory and sample graphlet backgrounds at early, middle, late,
  and stalled/final steps. Evaluate all residual candidates and controls on
  this common graphlet panel so differences in topology are not confounded by
  different sampled neighborhoods.
- For every sampled background, construct dense binary Hamming graphlets of
  dimension `5` to `8` from DCA/contact neighborhoods and matched controls.
  Each graphlet must be evaluated after embedding in the full thermodynamic
  background, not only as an isolated residual tensor.
- Split graphlets into calibration and test sets by both walk identity and
  contact neighborhood. Held-out test graphlets must come from adaptive walks
  and DCA/contact neighborhoods not used to select residual parameters.
- Calibrate residual-field parameter regimes only against aggregate
  distributions over calibration graphlets. Do not reject or tune individual
  graphlets to pass topology targets.
- Evaluate the selected regimes on held-out graphlets and report calibration
  versus test distributions separately.

## Candidate Residual Architecture

The primary candidate is a sparse contact-conditioned residual field:

- select DCA/contact neighborhoods from the TEM-1 alignment;
- add sparse pairwise and occasional higher-order local residual coefficients;
- draw coefficient magnitudes from heavy-tailed distributions with explicit
  low-order dominance;
- allow mixed signs so that contact modules can create frustration rather than
  only smooth compensation;
- scale residual strength relative to the local one-step DFE so it can
  sometimes flip signs without dominating the global thermodynamic core;
- evaluate all readouts after embedding in the smooth core.

The experiment is not a claim that this architecture is correct. It is a test
of whether this less hand-certified residual mechanism generalizes better than
the exact-match coefficient-control mechanism from RES-024.

## Required Controls

- `smooth_core_only`: thermodynamic core with all residual and legacy
  epistasis layers disabled. This is the negative control for topology beyond
  diminishing returns.
- `res024_exact_match_haplotype`: DCA clique exact-match haplotype overlay
  matching the RES-024 mechanism and strength scale. This controls for
  high-order coefficient insertion without aggregate topology calibration.
- `dca_pair_override`: DCA state-specific pair terms matching the RES-023
  mechanism. This controls for strong selected-pair sign epistasis without
  high-order residual fields.
- `contact_residual_ensemble`: the primary aggregate-calibrated residual
  candidate.
- `shuffled_contact_residual`: same residual generator and coefficient
  distribution as the primary candidate, but with DCA/contact neighborhoods
  permuted across positions. This tests whether biological placement matters.
- `random_position_residual`: same coefficient distribution placed on
  randomly selected positions and states matched for graphlet dimension and
  local single-mutant DFE. This controls for density and strength alone.
- `gnk_local_null`: a local GNK-style non-procedural comparator matched for
  graphlet dimension, coefficient variance, and interaction order envelope.
  This is a statistical null, not a proposed procedural replacement.
- `hard_certified_graphlet_positive_control`: per-graphlet rejection sampling
  to satisfy topology targets. This is an overfit positive control and should
  not be treated as a scalable model. Its purpose is to show the maximum
  achievable topology when generalization is not required.
- `residual_only_diagnostic`: evaluate residual tensors without the smooth
  thermodynamic core, then compare against embedded results. This identifies
  cases where the core erases or manufactures apparent topology.

## Calibration Targets

Calibration targets are aggregate envelopes, not per-graphlet constraints:

- SE and RSE prevalence distributions;
- number of local maxima and local minima;
- basin-size distribution and entropy where graphlet size permits;
- fraction of accessible improving shortest paths between selected endpoints;
- random-walk lag-1 and lag-5 autocorrelation on each graphlet;
- GFT coefficient distribution by graph eigenvalue band;
- Walsh coefficient and variance distributions by order;
- high-order Walsh fraction present but not dominant;
- one-step DFE shape, neutral fraction, and deleterious/beneficial fractions
  retained within the smooth-core biological envelope.

## Primary Readouts

- Calibration and held-out test distributions for SE fraction and RSE fraction.
- Calibration and held-out test distributions for local maxima count, local
  minima count, basin count, and basin entropy.
- Calibration and held-out test accessibility summaries, including accessible
  improving shortest-path fraction between graphlet endpoints and greedy path
  length from each node.
- Calibration and held-out test random-walk autocorrelation summaries for each
  graphlet.
- Calibration and held-out test Walsh terms for every graphlet: coefficient,
  absolute coefficient, order, variance contribution, and order-level variance
  fraction.
- Calibration and held-out test GFT coefficients for every graphlet: eigenvalue,
  coefficient, absolute coefficient, squared coefficient, frequency band, and
  band-level energy fraction.
- Stratified summaries by start class, adaptive-walk step class, background
  fitness quantile, local additive slope, local DFE variance, neutral fraction,
  and graphlet dimension.

## Falsification Criteria

The primary candidate is not supported if any of the following occur:

- held-out graphlets do not improve over smooth-core and exact-match controls
  on RSE, local maxima, accessibility, and random-walk autocorrelation;
- calibration graphlets improve but held-out graphlets do not;
- improvements occur only in one adaptive-walk start class or graphlet
  dimension;
- residual-only diagnostics show topology that disappears after embedding in
  the smooth core;
- the model achieves topology only by destroying local DFE realism or making
  high-order Walsh/GFT energy dominate unrealistically;
- shuffled or random-position controls match the primary candidate, implying
  that DCA/contact placement is not contributing biological structure.

## Inputs

- `config.yaml`
- runner: `scripts/tem1_aggregate_graphlet_residual_ruggedness.py`
- `external/Adapt-Env`
- `external/landscapy`
- TEM-1 MMSeqs/ColabFold alignment fetched or cached through the runner.

## Expected Outputs

- `tem1_mmseqs_alignment.fasta`
- `adaptive_walk_starts.csv`
- `adaptive_walk_trajectories.csv`
- `graphlet_catalog.csv`
- `calibration_graphlet_metrics.csv`
- `test_graphlet_metrics.csv`
- `graphlet_walsh_terms.csv`
- `graphlet_gft_coefficients.csv`
- `graphlet_sign_epistasis_pairs.csv`
- `graphlet_accessibility.csv`
- `module_parameter_sweep.csv`
- `calibration_target_summary.json`
- `heldout_test_summary.json`
- `summary.json`

## Analysis Plan

After completion, create an ANA record that checks:

- whether the primary contact residual ensemble improves held-out graphlet
  topology relative to smooth-core, RES-023, and RES-024 controls;
- whether calibration and held-out distributions agree closely enough to rule
  out obvious local overfitting;
- whether improvements persist across distinct adaptive-walk start classes and
  walk-step strata;
- whether DCA/contact placement beats shuffled and random-position controls;
- whether GNK-style local nulls match spectra but fail or succeed on topology;
- whether WHT and GFT spectra remain empirical-looking: low-order energy
  dominant, high-order energy present but not dominant, and no deterministic
  exact-match coefficient shelves;
- whether local DFE shape remains within the smooth-core biological envelope.

## Completion Criteria

- [x] Config created
- [x] Runner implemented
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated with the planned experiment

## Post-experiment notes

- `scripts/tem1_aggregate_graphlet_residual_ruggedness.py` now implements the
  common smooth-core adaptive-walk graphlet panel, all preregistered controls,
  calibration/test split, Walsh and Hamming-graph GFT coefficient export, SE/RSE,
  local-extrema, basin, accessibility, and graph-random-walk metrics.
- A local quick smoke run completed successfully with a synthetic TEM-1-like
  alignment through `LABPROJ_QUICK=1`. This smoke run is an execution check
  only, not scientific evidence for the full experiment.
- Dry-run submission created `RUN-064`. The first live Slurm execution is
  `RUN-065`, submitted on 2026-05-15 as scheduler job `226`.
- `RUN-065` completed with scheduler exit code `0`, and the collected durable
  outputs are DVC-tracked under
  `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065.dvc`.
- `ANA-025` records the scientific review of `RUN-065`, and `RES-025` records
  the result. The primary contact residual produced held-out local topology in
  the sampled graphlets with low-order-dominated spectra, but all evaluated
  graphlets came from `reference_or_consensus_like` starts and matched controls
  did not support DCA/contact placement as the source of topology.
