# RES-043: RUN-115 bounded trait-space fast fit improves guardrails but fails PHOT full-ABC equivalence

## Summary

RUN-115 completed EXP-044 and found zero operationally equivalent bounded trait-space current-model candidates against the best available RUN-098 full-ABC PHOT comparator. The minimal model patch improved the previous gross component-guardrail failures, with 12/12 candidates passing component guardrails, but all rows failed absolute and relative aggregate-distance equivalence. Observable epistasis moved closer than previous nested-min calibrated fits, while raw latent magnitude epistasis remained tens to hundreds of units and held-out accessibility/local-energy tradeoffs persisted.

## Generated from

- Analyses: ANA-043

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/metrics.json
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_equivalence_decision.md
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_bounded_gate_test_distances.md
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_abc_comparator_test_distances.md
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_top_component_distances.md
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_observable_vs_raw_epistasis_and_landscape.md
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_best_candidate_parameters.md
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_output_collection_audit.md
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/source_tables/run-115_bounded_gate_equivalence.csv
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/source_tables/run-115_bounded_gate_distances.csv
- analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/source_tables/run-115_graphlet_summary.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs

## Interpretation

`RUN-115` completed the intended `EXP-044` test: a minimal current-model patch
bounded non-score epistatic perturbations in trait space, fitted PHOT by a
deterministic/random numerical fast search with local refinement, applied a
positive affine observable readout on calibration graphlets, and compared the
selected candidates with empirical PHOT and the previously computed `RUN-098`
full-ABC comparator.

The result does not validate the fast fitter as equivalent to the full-ABC
PHOT result. It produced `0/12` operationally equivalent rows. All selected
bounded candidates passed component-distance guardrails, but every row failed
both the absolute and relative aggregate-distance margins against the best
available full-ABC context row. The best bounded candidate test distance was
`0.4269`; the best full-ABC context row was `0.2115`.

The result is nevertheless materially different from `RES-042`. The previous
calibrated nested-min cheap protocol had `0/40` equivalent rows and failed
component guardrails for every row, with best test distance `1.1474`.
`RUN-115` reduced the best held-out distance and removed the severe
component-guardrail failures. The bounded trait-space route is therefore a
partial model-level improvement, not just another failed search protocol.

The remaining problem is calibration of the model geometry. Empirical PHOT
test-graphlet magnitude epistasis mean is `0.2926`. The best held-out bounded
candidate reaches `0.1708`, and the best calibration-objective candidate
reaches `0.1698`; both are closer than `RUN-112` but still low. Their raw
latent magnitude epistasis means are `103.65` and `67.84`, respectively. The
observable readout is therefore compressing a badly scaled latent surface
rather than demonstrating that the underlying current-model epistasis is
empirically calibrated.

Accessibility also remains unstable across the selected tradeoff. The best
held-out candidate has accessible-path fraction `0.0120` versus empirical
`0.1975`. The best calibration-objective candidate nearly matches
accessibility (`0.1986`), but has a worse aggregate test distance (`0.4744`)
and still undershoots observable magnitude epistasis.

## Effect on hypothesis

This weakens the PHOT adequacy form of `HYP-004`. The bounded current-model
patch is not yet a procedural landscape model whose fitted PHOT graphlets can
be treated as empirically matched for accessibility, epistasis, spectra, and
local-energy structure.

For `HYP-009`, the result is methodological motivation rather than support. It
shows that a model-level epistasis route can improve matched-summary behavior
more than readout-only calibration, but it still does not provide a validated
scale-up protocol for GB1/TrpB or a benchmark-ready matched landscape.

## Limitations

- PHOT remains the only empirical system in this calibration test.
- The full-ABC comparator is represented by selected `RUN-098` context rows,
  not by the full posterior distribution.
- The formal equivalence decision uses the best available full-ABC context
  row, including the empirical-readout residual diagnostic; the best bounded
  row is better than the best nested full-ABC family row, so the result should
  not be overread as a total failure of bounded trait-space routing.
- Raw WHT/GFT coefficient tables were not copied into the analysis; aggregate
  spectral summaries were copied and the full coefficient tables remain on the
  runtime volume.
- The positive affine observable readout can improve reported scale without
  fixing latent geometry. Both observable and raw summaries must remain in
  future records.

## Downstream use

- Do not use the current `EXP-044` bounded fast-fit protocol as a validated
  full-ABC replacement for complete GB1 or TrpB landscape fitting.
- Keep bounded trait-space epistasis as a promising model-level direction:
  it reduces the severe component-guardrail failures seen in `RES-042`.
- Treat the next model revision as a calibration problem in latent trait
  geometry, not only as a cheaper optimizer problem.
- Continue reporting accessibility, local energy, spectral, observable
  epistasis, and raw latent epistasis together; any readout-only success is
  insufficient for the intended procedural landscape use.
