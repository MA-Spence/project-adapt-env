# RES-062: RUN-162 formally fails observable coverage and specificity with a plateau-confounded rugged target

## Summary

RUN-162 independently reproduces the preregistered
`OBSERVABLE_CONSTRUCTION_LIMIT` decision. G1 range coverage passes for only
`1/3` rugged-end templates, G2 joint reachability and G4 navigability pass for
`2/3`, and G3 specificity passes for `0/6` descriptors. EXP-065 therefore
remains blocked.

The strongest biological construction-limit interpretation is withheld. The
local-maxima-density rugged-end selector misclassifies flat graphlets as
maximally rugged: CR9114-FLUB is flat in all `96/96` frozen graphlets and
CR9114-H3 is flat in `69/96`.

## Generated from

- Analysis: ANA-062
- Experiment: EXP-064
- Run: RUN-162

## Relevant hypotheses

- Supports: None
- Weakens: HYP-005
- Refutes: None

## Evidence

- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/metrics.json`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_acceptance_audit.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_coverage_by_template.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_coverage_by_descriptor.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_rugged_end_audit.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_rugged_target_diagnostic.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_specificity_audit.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_specificity_by_descriptor.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_operating_region_audit.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_protocol_output_audit.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_source_hash_audit.csv`
- `analyses/2026-07-24_ANA-062_run-162-observable-ruggedness-coverage-probe-review/tables/run-162_collection_audit.csv`
- `data/processed/observable-ruggedness-coverage-probe/RUN-162.dvc`

## Interpretation

The run is technically valid. Slurm job `2045` exits `0`; all 300
template-seed tasks and 259,200 generated grid rows are present; all 14
declared outputs exist; all `1,526/1,526` technical checks pass; the frozen
ANA-050 hashes and published-band checks reproduce; and an independent
reaggregation matches every gate and the terminal decision to numerical
precision.

Only CR9114-H1 meets the G1 `5/8` range rule. CR9114-FLUB and CR9114-H3 cover
only `3/8` descriptors. G2 and G4 are `2/3 = 0.667`, below the required `0.8`:
CR9114-H1 and CR9114-H3 each have navigable seven-descriptor in-band
configurations, but CR9114-FLUB never exceeds three.

The broadest descriptor walls are local-maxima density (`0/10` Tier-1
templates covered), sign epistasis (`2/10`), and reciprocal-sign epistasis
(`2/10`). For CR9114-H1, the structured binary-rank-8 maximum is `0.27083`
against sign-epistasis rugged q90 `0.70833`, and `0.02083` against
reciprocal-sign q90 `0.25`. The generated local-maxima density reaches only
`0.0625` on each rugged-end template.

Specificity independently fails. No increasing-direction descriptor passes the
majority-of-rugged-templates rule. Structured maxima for Walsh high-order
fraction, roughness-to-slope, Dirichlet localisation, and local-maxima density
are often indistinguishable from or worse than matched shuffled/random
controls. Sign and reciprocal-sign specificity each pass for only one of three
templates. Thus the RES-052 identifiability problem remains even where
coverage is possible.

The rugged-end target is not a clean ruggedness target. Under the registered
strict rule, a node with no fitter neighbour is a local maximum, so every node
in a constant 4-cube is counted and the density is `1.0`. CR9114-FLUB has zero
fitness variance, 16 local maxima, zero Walsh/sign/reciprocal-sign/roughness
signal, and undefined random-walk and Dirichlet descriptors in all 96 frozen
graphlets. CR9114-H3 is flat in `71.875%` of graphlets. Much of the binding
failure is therefore inability to reproduce flat-plateau degeneracy, not
demonstrated inability to generate biological multi-peak ruggedness.

The defensible conclusion is that EXP-064 fails its registered
expressivity/specificity prerequisite. It does not establish a general
observable or sigmoid impossibility.

## Effect on hypothesis

RES-062 weakens HYP-005 narrowly. Freeing `beta` and `p0` and using the
experiment-local binary-rank-8 field is not sufficient to specifically cover
the registered Tier-1 topographic envelope: local-maxima and nonzero
sign-epistasis structure remain deficient, and structured routing does not
separate from controls.

This is not a direct test of HYP-004 adaptive-walk dynamics or HYP-009 optimizer
effects. The plateau-confounded rugged-end definition prevents a strong
biological ruggedness conclusion, so neither hypothesis is weakened or
supported by this run.

## Limitations

- Coverage/expressivity probe only; forward generation; operator-free
  (`uses_exp050_transfer_function=false`); descriptor-space, with no empirical
  fitness values loaded.
- Coverage is necessary, not sufficient for fitting. No per-landscape fitting
  was attempted.
- No monotone-knob claim; RES-051 remains the standing monotonicity result.
- The top-quartile strict local-maxima-density selector conflates flat plateaus
  with rugged multi-peak landscapes.
- CR9114-FLUB lacks finite empirical random-walk and Dirichlet-localisation
  descriptors; these were left unavailable, not imputed.
- Thresholds are preregistered noiseless engineering criteria, not empirical
  detection limits.
- The repository and AdaptEnv checkout were dirty at submission. Frozen
  empirical inputs are hash-locked, but not every implementation source has a
  recorded run-time hash.
- No realism, public-generator readiness, optimizer, or active-learning claim.
  The public model is unchanged and the Stage 4C/4D/4E hard stop remains.

## Downstream use

- Cite RES-062 for the formal EXP-064 failure and for continued lack of
  structured-versus-control specificity.
- Cite the exact local-maxima and sign/reciprocal-sign shortfalls, not a generic
  statement that all ruggedness descriptors are unreachable.
- Do not cite RUN-162 as proof that the sigmoid or rank-8 observable cannot
  generate biological ruggedness; the binding target set is plateau-confounded.
- Do not use `selected_operating_region.csv` as a fitting prior and do not
  treat EXP-065 as unblocked.
- Keep the public-model validation state and Stage 4C/4D/4E hard stop unchanged.
