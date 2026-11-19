# RES-041: RUN-110 cheaper PHOT nested-min fitting matches same-family full ABC but not diagnostic ceiling

## Summary

RUN-110 shows surrogate-assisted and exact-descriptive-gated cheaper fitting match the RUN-098 full-ABC selected nested-family PHOT fits across all five families under same-family equivalence margins, while deterministic summary loss and local ABC initialized from cheap fits match four of five. No cheaper protocol reaches the best empirical-readout-residual diagnostic comparator, so the result supports cheaper fitted-family PHOT protocol use but not a general claim that fast fitting is equivalent to full ABC.

## Generated from

- Analyses: ANA-041

## Relevant hypotheses

- Supports: None
- Weakens: None
- Refutes: None
- Motivates: HYP-009

## Evidence

- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/metrics.json
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_protocol_decision_table.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_same_family_equivalence.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_best_any_equivalence.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_test_distances.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_abc_comparator_test_distances.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_secondary_epistasis_test_summary.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_module_locality_test.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/tables/run-110_output_collection_audit.md
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/source_tables/run-110_fit_protocol_equivalence.csv
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/source_tables/run-110_fit_protocol_distances.csv
- analyses/2026-06-16_ANA-041_run-110-phot-cheap-fitting-equivalence-review/source_tables/run-110_summary.json
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-110/outputs

## Interpretation

`RUN-110` directly addresses the practical question behind `EXP-042`: can a
cheaper fitting protocol recover the selected full-ABC fitted-family result for
the current PHOT nested min-coupled sequential residual model family?

For same-family comparisons, the answer is yes for the two strongest cheap
protocols. `surrogate_assisted` and `exact_descriptive_gated` each pass
same-family equivalence for all five nested families. `deterministic_summary_loss`
and `local_abc_after_cheap_fit` each pass four of five families, failing only
the random-module family.

The best cheap fitted candidate is the `surrogate_assisted` random-module
candidate with held-out test distance `0.428849`. This is better than the
`RUN-098` selected full-ABC random-module comparator (`0.482630`) and better
than every selected fitted nested-family row in the `RUN-098` comparator table.

The result is not a general equivalence proof. No cheap protocol reaches the
best-any comparator row, which is the empirical readout residual diagnostic
with distance `0.211455`. That diagnostic is a ceiling showing remaining PHOT
residual structure, not a fitted nested-family ABC particle. The correct
claim is therefore same-family operational equivalence to selected full-ABC
fitted-family behaviour, not equivalence to the empirical residual diagnostic
or to the full ABC posterior distribution.

The secondary descriptive summaries preserve an important caveat. Selected
fitted nested-family candidates remain far from empirical PHOT in raw epistasis
magnitude: empirical test-graphlet magnitude epistasis mean is about `0.293`,
whereas selected fitted nested-family candidates are tens to hundreds in the
same units. The module-locality summaries also remain diffuse for the
empirical-module family, so this result does not rescue residue-level empirical
module specificity.

## Effect on hypothesis

`RES-041` motivates `HYP-009` methodologically. It identifies two tractable
protocols, `surrogate_assisted` and `exact_descriptive_gated`, that can match
the selected full-ABC fitted-family PHOT graphlet summary results without
rerunning brute-force ABC.

It does not directly support `HYP-009`'s optimizer-performance prediction,
because no PLM, active-learning, or evolutionary optimizer was tested. It also
does not strengthen the old mechanism-specific empirical-module claim:
random/shuffled/no-module controls remain central, and empirical-module
locality remains diffuse.

## Limitations

- PHOT is the only empirical system in this run.
- The comparison uses the selected `RUN-098` full-ABC candidates, not the full
  ABC posterior distribution.
- Operational equivalence is defined on the predeclared held-out graphlet
  distance and component margins; it is not parameter equality.
- The best-any comparator is an empirical readout residual diagnostic rather
  than a fitted nested-family particle.
- Secondary epistasis and module-locality readouts show remaining mismatch, so
  this result should not be cited as a complete PHOT landscape match.
- The large WHT/GFT coefficient tables remain on the runtime volume and were
  not copied into the result record.

## Downstream use

- Use `surrogate_assisted` and `exact_descriptive_gated` as the leading cheap
  fitting protocols for the GB1/TrpB scale-up attempt.
- Treat local ABC initialized from cheap fits as an optional refinement, not as
  the default, because it did not clearly outperform the cheaper deterministic
  protocols in `RUN-110`.
- Do not claim that fast fitting is generally equivalent to full ABC. The
  supported claim is narrower: on PHOT, two cheaper protocols recover
  same-family selected full-ABC fitted-family graphlet behaviour under the
  predeclared margins.
