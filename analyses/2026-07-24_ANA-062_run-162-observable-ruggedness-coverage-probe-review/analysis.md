# ANA-062: RUN-162 observable-ruggedness coverage probe review

## Purpose

Collect and independently audit completed RUN-162 against the preregistered
EXP-064 technical, empirical-band reproduction, range-coverage, joint
rugged-end, specificity, navigability, decision, and provenance gates without
post hoc retuning.

## Linked experiments/runs

- Experiment: EXP-064
- Run: RUN-162
- Direct hypothesis: HYP-005
- Motivating result: RES-061

## Collection and reproducibility audit

`RUN-162` completed on `lab-slurm` as Slurm job `2045` with inferred scheduler
exit code `0`. The full registered run contains all 10 Tier-1 templates, 30
seeds, three representation arms, three support modes, four `beta` values, four
reference occupancies, and six amplitudes. The resulting
`generated_descriptor_metrics.csv` has the expected `259,200` unique grid
rows. All 300 template-seed tasks completed, their checkpoint row totals
reproduce the aggregate table, all 14 declared top-level outputs exist, and
`stderr.txt` is empty.

The collected bundle is:

- `data/processed/observable-ruggedness-coverage-probe/RUN-162`;
- 1,517 files and 204,553,256 bytes, including 300 checkpoint directories and
  three logs;
- DVC pointer
  `data/processed/observable-ruggedness-coverage-probe/RUN-162.dvc`;
- DVC directory hash `8428ff92c79a7d5b768a97772dbce230.dir`.

The workspace file count and byte count exactly match the DVC pointer, DVC
reports the bundle up to date, and all 303 JSON files parse under a strict
non-finite-token rejection rule. The collector initially emitted a redundant
nested `logs.dvc` after already including those same three logs in the parent
bundle. Removing only that nested pointer restored a single non-overlapping,
up-to-date DVC record; no data or logs were deleted.

All `1,526/1,526` recorded technical checks pass. The four frozen ANA-050 input
hashes match the present files, all published pooled bands and dataset medians
are reproduced within `1e-6`, the known additive-cube control passes, graphlet
and support pairing are complete, binary rank 1 replays the incumbent exactly,
operator-energy errors stay below tolerance, rank-8 effective-rank checks pass,
and dense-q20 ties are absent.

The independent analysis reaggregates the 259,200 raw generated rows rather
than trusting the runner's gate tables. It reproduces every gate and the
terminal decision. Maximum absolute differences from the runner are
`1.94e-16` for coverage, `2.43e-16` for specificity, and `4.44e-16` for the
operating-region table; all Boolean states match.

The repository and AdaptEnv checkout were dirty at submission. The frozen
empirical inputs are hash-locked, but the run record does not contain hashes
for every implementation source, so this remains a broader provenance
limitation.

## Registered and independently audited decision

The registered and independent terminal decision is:

`OBSERVABLE_CONSTRUCTION_LIMIT`

| Gate | Observed | Required | Pass |
|---|---:|---:|---:|
| Execution and pairing | 1.000 | 1.000 | yes |
| Empirical-band reproduction | 1.000 | 1.000 | yes |
| G1 rugged-template range coverage | 0.333 | 0.800 | no |
| G2 rugged-end joint reachability | 0.667 | 0.800 | no |
| G3 specificity-descriptor fraction | 0.000 | >0.500 | no |
| G4 navigable rugged-end fraction | 0.667 | 0.800 | no |

EXP-065 therefore remains blocked. There is no valid operating-point prior to
hand forward for fitting.

## Range coverage

Only CR9114-H1 passes the G1 per-template range rule:

| Rugged-end template | Available descriptors | Covered descriptors | Fraction over all 8 | G1 template pass |
|---|---:|---:|---:|---:|
| CR9114-FLUB | 6 | 3 | 0.375 | no |
| CR9114-H1 | 8 | 5 | 0.625 | yes |
| CR9114-H3 | 8 | 3 | 0.375 | no |

The failure is not confined to a single descriptor in one template. Across all
10 Tier-1 templates, local-maxima density is covered in `0/10`, sign epistasis
in `2/10`, reciprocal-sign epistasis in `2/10`, Walsh high-order fraction in
`5/10`, roughness-to-slope in `6/10`, and random-walk lag-1 autocorrelation in
`6/10`. Accessibility is covered in `10/10`, and Dirichlet-energy localisation
in `8/10` (one additional template has no empirical value).

For the three preregistered rugged-end templates, the strongest exact
shortfalls are:

- local-maxima density: `0/3` covered. CR9114-FLUB and CR9114-H3 require
  density `1.0`, while the structured binary-rank-8 generator reaches only
  `0.0625`; CR9114-H1 requires rugged q90 `0.1875` and also reaches only
  `0.0625`;
- CR9114-H1 sign epistasis: reachable maximum `0.27083` versus rugged q90
  `0.70833`, shortfall `0.43750`;
- CR9114-H1 reciprocal-sign epistasis: maximum `0.02083` versus rugged q90
  `0.25`, shortfall `0.22917`;
- CR9114-H3 random-walk autocorrelation misses the decreasing-direction target
  by `0.01057`, and Dirichlet localisation misses by `0.00043`.

For CR9114-FLUB, Walsh high-order fraction and roughness-to-slope reach the
rugged q90 of zero but cannot contain the empirical median of exactly zero
because every nonzero-amplitude generated value is positive. Its random-walk
and Dirichlet descriptors are unavailable. Thus its maximum joint count is
only three.

## Joint reachability and navigability

CR9114-H1 and CR9114-H3 each have at least one structured binary-rank-8
configuration within the empirical q05-q95 band on seven of eight descriptors,
and each has at least one such configuration with accessibility at least
`0.05`. CR9114-FLUB never exceeds three in-band descriptors and has no joint or
navigable-joint pass. G2 and G4 are therefore both `2/3 = 0.667`, below the
required `0.8`.

The best row-level operating regions do not rescue the result. Both
`(beta=1, p0=0.5, a=2)` and `(beta=2, p0=0.5, a=1)` achieve a joint-row
fraction of `0.644`; the latter is much more navigable (`0.967` versus
`0.211`), but neither is a passed fitting prior.

## Specificity

G3 fails completely: none of the six increasing-direction descriptors passes
the across-template majority rule. Sign epistasis passes its per-template
specificity margin in only CR9114-H3 (`1/3`); reciprocal-sign epistasis passes
only in CR9114-FLUB (`1/3`). The other four descriptors pass in `0/3`.

The structured and structureless reachable maxima are often nearly identical:

- Walsh high-order structured-minus-shuffled minima reach `-0.00273`, with
  structured-minus-random minima `-0.00059`;
- roughness-to-slope minima are `-0.11168` versus shuffled and `-0.07816`
  versus random;
- Dirichlet-localisation differences are at most numerical-scale relative to
  the required `0.05`, and can be negative;
- local-maxima density gains are zero against additive and shuffled controls,
  and negative against the random control for CR9114-H1.

This independently reproduces the central RES-052 identifiability problem:
even where descriptor bands are reachable, the structured construction does
not separate from energy/degree-matched shuffled and random routing. The formal
terminal remains the earlier construction-limit branch because G1/G2 fail, but
the fitting experiment is also blocked by nonspecificity.

## Rugged-end target confound

The preregistered rugged-end definition is materially confounded by plateaus.
It selects the top quartile by strict local-maxima density, where a node with no
strictly fitter neighbour counts as a local maximum. Under that definition a
constant graphlet has all 16 nodes classified as local maxima and density
`1.0`, even though it has zero variance and no rugged signal.

The frozen RUN-135 source shows:

| Template | Flat graphlet fraction | Median local-maxima density | Finite random-walk/Dirichlet fraction |
|---|---:|---:|---:|
| CR9114-FLUB | 1.00000 | 1.000 | 0.00000 |
| CR9114-H1 | 0.00000 | 0.125 | 1.00000 |
| CR9114-H3 | 0.71875 | 1.000 | 0.28125 |

CR9114-FLUB has zero fitness variance in all `96/96` frozen graphlets; each has
16 local maxima, zero Walsh high-order fraction, zero sign and reciprocal-sign
epistasis, zero roughness-to-slope, and undefined random-walk and
Dirichlet-localisation metrics. CR9114-H3 is flat in `69/96` graphlets.

Consequently, the formal `OBSERVABLE_CONSTRUCTION_LIMIT` decision is valid
against the preregistered tables but is not clean evidence that the observable
cannot reach a biologically rugged upper envelope. Much of the binding
"rugged-end" failure is inability to reproduce exact flat-plateau degeneracy.
This is a flaw in the registered target construction, not an execution or
descriptor-pipeline error. The result still blocks fitting because both the
formal gate and the independent specificity test fail.

## Interpretation for HYP-005

RUN-162 weakens the scoped prediction that freeing sigmoid steepness/reference
occupancy and using the RES-061 binary-rank-8 field is sufficient to
specifically cover the frozen Tier-1 topographic envelope. Local-maxima,
nonzero sign/reciprocal-sign epistasis, and structured-versus-control
separation remain unresolved.

The evidence does not justify a blanket claim that sigmoid transmission or the
rank-8 observable is incapable of empirical ruggedness. The rugged-end target
is plateau-confounded, and two of three selected templates do have navigable
seven-descriptor in-band configurations under the looser G2 criterion. The
defensible conclusion is a formal expressivity/specificity failure of EXP-064
as registered, with the strongest biological construction-limit interpretation
withheld.

This run is not a direct test of HYP-004 adaptive-walk dynamics or HYP-009
optimizer effects. No claim-state change is warranted.

## Scope and limitations

- Coverage/expressivity probe only; forward generation; operator-free
  (`uses_exp050_transfer_function=false`); descriptor-space, with no empirical
  fitness values loaded.
- Coverage is necessary, not sufficient for per-landscape fitting. This run
  fails that prerequisite and does not attempt a fit.
- No monotone-knob claim; RES-051 remains the standing monotonicity result.
- No realism, public-generator readiness, optimizer, or active-learning claim.
  The public model is unchanged.
- Rugged/topographic claims are restricted to the 10 Tier-1 operator-free
  landscapes; Tier-3 assays remain excluded.
- Thresholds are preregistered noiseless engineering criteria, not empirical
  detection limits.
- The top-quartile local-maxima-density target conflates flat plateaus with
  rugged multi-peak structure and materially limits the biological
  interpretation.
- The source-staged repository and AdaptEnv checkout were dirty; only the
  frozen empirical inputs have complete recorded source hashes.
- The Stage 4C/4D/4E hard stop remains in force.

## Produced artifacts

- `tables/metrics.json`
- `tables/run-162_acceptance_audit.csv`
- `tables/run-162_coverage_by_template.csv`
- `tables/run-162_coverage_by_descriptor.csv`
- `tables/run-162_rugged_end_audit.csv`
- `tables/run-162_rugged_target_diagnostic.csv`
- `tables/run-162_specificity_audit.csv`
- `tables/run-162_specificity_by_descriptor.csv`
- `tables/run-162_operating_region_audit.csv`
- `tables/run-162_protocol_output_audit.csv`
- `tables/run-162_source_hash_audit.csv`
- `tables/run-162_collection_audit.csv`
