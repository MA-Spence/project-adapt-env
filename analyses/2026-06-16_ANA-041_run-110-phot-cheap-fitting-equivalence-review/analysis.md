# ANA-041: RUN-110 PHOT cheap fitting equivalence benchmark review

## Purpose

Review RUN-110 compact outputs, compare cheaper PHOT nested-min fitting protocols against the RUN-098 full ABC comparator, and decide whether any protocol is operationally equivalent enough to justify larger landscape scale-up.

## Linked experiments/runs

- Experiments: EXP-042
- Runs: RUN-110

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-110_summary.json`
- `source_tables/run-110_preflight.json`
- `source_tables/run-110_fit_protocol_equivalence.csv`
- `source_tables/run-110_fit_protocol_distances.csv`
- `source_tables/run-110_selected_fit_protocol_candidates.csv`
- `source_tables/run-110_abc_comparator_selected_family_distances.csv`
- `source_tables/run-110_full_abc_comparator_selected_particles.csv`
- `source_tables/run-110_graphlet_summary.csv`
- `source_tables/run-110_module_energy_summary.csv`
- `source_tables/run-110_local_abc_rounds.csv`
- `source_tables/run-110_graphlet_catalog.csv`
- `source_tables/run-110_remote_output_manifest.tsv`
- `tables/metrics.json`
- `tables/run-110_output_collection_audit.md`
- `tables/run-110_protocol_decision_table.md`
- `tables/run-110_same_family_equivalence.md`
- `tables/run-110_best_any_equivalence.md`
- `tables/run-110_equivalence_by_method_reference.md`
- `tables/run-110_equivalence_by_method_family.md`
- `tables/run-110_test_distances.md`
- `tables/run-110_abc_comparator_test_distances.md`
- `tables/run-110_secondary_epistasis_test_summary.md`
- `tables/run-110_module_locality_test.md`
- `tables/run-110_local_abc_rounds.md`

## Collection Policy

`RUN-110` completed on `lab-slurm` with scheduler job `1906` and exit code
`0`. The full run output remains on the runtime volume at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-110/outputs`

This analysis deliberately did not collect the larger detailed output files.
The two largest files were `gft_coefficients.csv` (`28.7 MB`) and
`walsh_terms.csv` (`28.3 MB`), followed by `graphlet_metrics.csv` (`1.7 MB`)
and the complete candidate table (`1.0 MB`). Only compact JSON/CSV outputs
needed for the result were copied into `source_tables/`; the copied payload is
small enough for the scientific record and is listed in
`tables/run-110_output_collection_audit.md`.

## Output Audit

`RUN-110` evaluated `24` PHOT graphlets, evaluated `730` candidate parameter
sets, and selected `20` candidates for final comparison. The experiment
generated `40` equivalence rows: four fitting protocols by five candidate
families by two comparator references.

The comparison uses two reference notions:

- `same_family_full_abc`: compare a cheaper protocol's selected candidate
  against the `RUN-098` full-ABC selected candidate from the same nested family.
- `best_any_full_abc`: compare against the best row in the RUN-098 comparator
  table. In this run that row is the empirical readout residual diagnostic
  (`distance = 0.211455`), not a fitted nested-family particle.

The primary decision should therefore be read as same-family operational
equivalence to the full ABC fitted-family result. The best-any diagnostic
comparison is a useful ceiling, but failing it is not the same as failing to
replace the full ABC fitted-family search.

## Main Result

Two cheaper protocols passed same-family equivalence for all five nested
families:

- `surrogate_assisted`: `5/5` same-family equivalent.
- `exact_descriptive_gated`: `5/5` same-family equivalent.

Two protocols mostly passed:

- `deterministic_summary_loss`: `4/5` same-family equivalent, failing only the
  random-module family.
- `local_abc_after_cheap_fit`: `4/5` same-family equivalent, also failing only
  the random-module family.

No protocol was equivalent to the best-any comparator row. That best-any row is
the empirical readout residual diagnostic (`0.211455`), while the best cheaper
fitted candidate was `surrogate_assisted` random-module (`0.428849`). The best
RUN-098 selected fitted nested family was random-module (`0.482630`).

This means `RUN-110` gives positive evidence that a cheaper fitting protocol
can recover the operational behaviour of the full ABC selected nested-family
fits on PHOT. It does not show that the cheaper protocols recover the
empirical residual diagnostic, the full ABC posterior distribution, or a
general solution for larger landscapes.

## Protocol-Level Interpretation

`surrogate_assisted` is the strongest cheap protocol in this run. It passed
same-family equivalence for all five families and produced the best fitted
candidate overall: random-module test distance `0.428849`, better than the
`RUN-098` selected random-module comparator (`0.482630`).

`exact_descriptive_gated` also passed same-family equivalence for all five
families. Its best fitted candidate was no-module (`0.501041`), close to the
best local-ABC no-module candidate (`0.500039`) and better than the `RUN-098`
selected no-module comparator (`0.537199`). This is the cleanest evidence that
the `RES-040` exact-gate summaries can be used as a cheaper fitting route
without obvious loss on the original held-out graphlet summary vector.

`local_abc_after_cheap_fit` did not obviously add enough to justify itself as
the default cheap protocol. It matched same-family ABC for four families, but
failed random-module and did not beat the simpler exact-gated no-module result.
It may still be useful as a local refinement when a cheap deterministic or
surrogate fit is close but not inside tolerance.

`deterministic_summary_loss` is usable but weaker than the two best protocols:
it matched four families and failed random-module, with best test distance
`0.516987` for shuffled-module.

## Secondary Landscape Properties

The secondary graphlet summary table reports sign epistasis, reciprocal sign
epistasis, accessibility, magnitude epistasis, signed pairwise epistasis, and
diminishing-return fractions for the selected cheap candidates and for the
selected `RUN-098` full-ABC particles re-evaluated through the same reporting
path.

Those rows should be treated as descriptive checks, not as the primary
equivalence decision. They show that the fitted nested-family candidates remain
far from the empirical PHOT target in raw epistasis magnitude: empirical
test-graphlet magnitude epistasis mean is about `0.293`, while selected fitted
nested-family candidates are tens to hundreds in the same units. This reflects
the existing scale/latent-readout limitations of the PHOT nested-family
comparison and should not be hidden by the same-family distance equivalence.

The module-locality table also preserves the previous mechanism-specific
warning. Empirical-module selected candidates remain diffuse over modules on
test graphlets. For example, the `RUN-098` empirical-module selected particle
has edge effective module count `11.72/12` and edge normalized entropy
`0.995`, and the cheap empirical-module candidates have similarly diffuse
edge-module entropy. This continues to argue against interpreting empirical
module identities as recovered residue-level mechanism.

## Decision

`RUN-110` supports using `surrogate_assisted` and `exact_descriptive_gated` as
candidate cheaper PHOT fitting protocols for the nested min-coupled sequential
residual model family. These protocols are operationally equivalent to the
selected full ABC fitted-family results under the predeclared same-family
margins.

The result does not justify the stronger claim that fast fitting is generally
equivalent to full ABC. It is one PHOT-system validation, on one graphlet split,
against selected full-ABC candidates rather than the full posterior. Before
large GB1/TrpB production use, the prudent interpretation is that
`surrogate_assisted` and `exact_descriptive_gated` are now the leading cheap
protocols to try, with a smaller local-ABC refinement as a fallback rather
than the default.
