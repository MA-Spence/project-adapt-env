# EXP-064: Observable-ruggedness coverage probe

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009

## Question

Can the additive-core plus binary-rank-8 extensive-GNK stability-gate observable, with beta and reference occupancy swept, specifically cover the rugged upper quantiles of the frozen Tier-1 empirical descriptor envelope?

## Pre-experiment prediction

At least one beta, p0, and GNK-amplitude region will place the structured binary-rank-8 observable within at least five of eight descriptor bands for at least 80% of rugged-end Tier-1 templates while retaining accessibility and exceeding additive, shuffled, and random controls.

## Rationale

RES-061 repairs observable texture prevalence on the fixed challenge panel, but RES-051/052 leave forward empirical-envelope coverage and specificity unresolved and used the prior representation or frozen deployment operating point. This probe tests that expressivity prerequisite without fitting.

## Experimental design

This is a paired forward-generation coverage probe, not a fit. For each
Tier-1 template and seed, the additive core, GNK field, support-aware binary
graphlets, and graphlet nodes are sampled once. GNK amplitude is an inner loop
and is excluded from the design key. The experiment sweeps:

- `beta = 1, 2, 4, 8`;
- reference occupancy `p0 = 0.5, 0.7, 0.9, 0.99`, with `mu = logit(p0)`;
- GNK amplitude `a = 0, 0.25, 0.5, 1, 2, 4`;
- `binary_r8` primary, `binary_r1` incumbent, and `dense_q20_r8` contrast;
- structured, shuffled, and random support routing;
- 30 generator seeds and all 10 frozen Tier-1 templates.

The observable is

`F(A,aG) = 10 sigmoid(mu + beta(A-Aref)/sigma_A + aG/sigma_G)`,

with `sigma_A = 2.7448641782814067`. The operating point is an
experiment-local evaluation override; the public preset and public AdaptEnv
model are unchanged.

Generated and empirical graphlets use the exact RUN-135
`metrics_from_cube` descriptor path. The frozen RUN-135 source graphlet table
supplies per-template q05, median, q90, and q95 bands. Their pooled summaries
and per-dataset medians must reproduce the published ANA-050 tables before any
scientific decision is allowed. RUN-135 has no finite random-walk
autocorrelation or edge-energy localisation value for the flat CR9114-FLUB
graphlets; those two descriptor/template cells are retained explicitly as
unavailable rather than imputed.

The rugged-end subset is the preregistered top quartile of templates by median
local-maxima density. G1 requires range coverage through the empirical
rugged-direction q90. G2 requires one configuration in at least five of eight
descriptor bands for at least 80% of rugged-end templates. G3 requires
structured binary-rank-8 rugged extremes to exceed additive-core-only,
shuffled-support, and random-support controls by the frozen descriptor-specific
margins. G4 requires the G2 configuration to retain accessible-path fraction
at least `0.05`.

The terminal decisions are `PROCEED_TO_FITTING`,
`OBSERVABLE_CONSTRUCTION_LIMIT`, `COVERAGE_NONSPECIFIC`,
`NAVIGABILITY_DEGENERATE`, or a technical/band-ingestion invalidity. Quick
mode always returns `QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION`.

## Inputs

- Frozen RUN-135 graphlet metrics, scope audit, key bands, and dataset medians.
- The DVC-tracked ensemble landscape dataset bundle.
- The EXP-057 support-aware graphlet machinery.
- The EXP-063 ranked categorical GNK field and incumbent replay.
- The EXP-062 additive-core scale and direct dimensionless sigmoid formula.

## Configuration

`config.yaml` freezes input hashes, the complete sweep, descriptor aliases and
directions, specificity margins, gate thresholds, decision precedence, quick
mode, and scope guards.

## Execution

- Focused tests:
  `python -m pytest -q tests/test_observable_ruggedness_coverage_probe_exp064.py tests/test_ranked_categorical_gnk.py tests/test_dense_descriptive.py`
- Quick smoke:
  `LABPROJ_QUICK=1 bash experiments/2026-07-24_EXP-064_observable-ruggedness-coverage-probe/run.sh`
- Queue submission:
  `labproj submit EXP-064 --target lab-slurm --dry-run`, then `--execute`.

The final inspected dry-run is `RUN-159`. `RUN-160` and `RUN-161` failed before
reaching the scheduler and have no scheduler IDs. The corrected live submission
is `RUN-162`, HomeLab Slurm job `2045` (8 CPUs, 64 GB, 24 hours).

## Expected outputs

- preflight.json
- template_catalog.csv
- support_catalog.csv
- support_aware_graphlet_audit.csv
- empirical_bands_resolved.csv
- generated_descriptor_metrics.csv
- coverage_matrix.csv
- rugged_end_reachability.csv
- specificity_table.csv
- selected_operating_region.csv
- technical_checks.csv
- checkpoint_resume_summary.json
- decision_table.csv
- summary.json

## Analysis plan

Audit source hashes, empirical-band reproduction, pairing, support hashes,
operator energy, rank, dense-code ties, incumbent replay, and the additive
synthetic descriptor control before reading G1-G4. Independently reproduce the
coverage matrix, ranked descriptor shortfalls, specificity contrasts,
rugged-end joint configurations, navigability guard, and terminal decision.

Coverage is necessary but not sufficient for per-landscape fitting. Even
`PROCEED_TO_FITTING` only unblocks EXP-065; it does not validate a fit, a
monotone knob, biological realism, the public generator, or optimizer
performance. Tier-3 assays are excluded and the Stage 4C/4D/4E hard stop is
retained under every outcome.

## Completion criteria

- [x] Config implemented
- [x] Run script implemented
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- One-template quick smoke completed locally with all technical and frozen-band
  checks passing and the required no-decision branch.
- `RUN-162` completed with inferred exit code `0`. ANA-062 independently
  reproduces all technical and G1-G4 states and the terminal
  `OBSERVABLE_CONSTRUCTION_LIMIT` decision. RES-062 records the formal
  expressivity/specificity failure and the plateau confound in the registered
  rugged-end target. EXP-065 remains blocked.
