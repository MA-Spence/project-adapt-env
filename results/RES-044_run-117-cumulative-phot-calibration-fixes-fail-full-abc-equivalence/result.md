# RES-044: RUN-117 cumulative PHOT calibration fixes fail full-ABC equivalence

## Summary

RUN-117 completed the four-stage EXP-045 PHOT cumulative calibration diagnostic
and found `0/32` operationally equivalent rows. The full `A+B+C+D` stage was
the best RUN-117 stage, with held-out distance `0.4122` and observable
magnitude epistasis `0.2426`, but it still failed the best available RUN-098
comparator distance of `0.2115` and retained low accessibility and pathological
raw-scale behavior.

## Generated from

- Analyses: ANA-044

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/metrics.json
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_stage_equivalence_decision.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_stage_best_test_distances.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_bounded_gate_test_distances.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_abc_comparator_test_distances.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_top_component_distances.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_observable_vs_raw_epistasis_and_landscape.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_stage_failure_modes.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_best_candidate_parameters.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/run-117_output_collection_audit.md
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/source_tables/run-117_bounded_gate_equivalence.csv
- analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/source_tables/run-117_graphlet_summary.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-117/outputs

## Interpretation

`RUN-117` directly tests the cumulative A-through-D diagnosis from the
post-`RES-043` calibration review. It keeps the PHOT source and RUN-098
comparator fixed while adding, in order, a monotone-piecewise observation link,
a global sigmoid readout with constrained specific epistasis, an explicit
specific-epistasis/additive-scale ratio knob, and a distribution/scale
conditioned objective.

The experiment does not validate the cumulative calibration protocol. The
formal equivalence decision produced `0/32` operationally equivalent rows. All
selected rows failed the absolute and relative aggregate-distance margins
against the best available RUN-098 full-ABC context row, the
`empirical_readout_residual` comparator at test distance `0.2115`.

The stage pattern matters. Stage `A` partially preserves the `RES-043`
bounded-route improvement: all `8/8` selected rows pass component guardrails,
and the best held-out row reaches test distance `0.4235`. It also improves
observable magnitude epistasis to `0.2289` versus empirical `0.2926`. But
accessibility remains low (`0.0164` versus empirical `0.1975`), and raw latent
magnitude epistasis remains very large (`178.6`).

Stages `A+B` and `A+B+C` are clear failures under this search. The selected
rows collapse to zero-variation landscapes with zero single-mutant SD, zero
observable magnitude epistasis, zero accessibility, and no component-guardrail
passes. Adding the explicit specific-epistasis/additive-scale ratio knob in
`C` does not rescue the collapse.

The full `A+B+C+D` stage is the best of the four. It restores component
guardrail passes (`8/8`), reaches the best RUN-117 held-out distance
(`0.4122`), and improves observable magnitude epistasis to `0.2426` with a
magnitude-epistasis-to-single-mutant-SD ratio of `1.1204`, close to empirical
`1.0465`. This is also a small improvement over the best RUN-115 held-out
distance reported in `RES-043` (`0.4269`).

That improvement is not sufficient. The best `A+B+C+D` row is still nearly
twice the distance of the best available full-ABC comparator, retains low
accessibility (`0.0190` versus empirical `0.1975`), and has effectively absent
raw latent variation (`raw_magnitude_epistasis_mean_mean` about `1.16e-6`).
The result therefore shifts the failure mode rather than solving it: previous
bounded candidates had raw latent epistasis that was far too large, while the
best cumulative-stage candidate has raw latent variation that is essentially
zero.

## Effect on hypothesis

This weakens the current PHOT adequacy form of `HYP-004`. The current-model
PHOT generator, even after cumulative observation, global-readout, ratio, and
scale-conditioned objective fixes, still does not produce a held-out PHOT
landscape that can be treated as empirically matched across accessibility,
epistasis, spectra, and local-energy summaries.

For `HYP-009`, this is methodological motivation rather than support. The run
shows that distribution/scale-conditioned matching can improve the numerical
behavior of the fast fitter, but it does not yet supply a validated matched
landscape for benchmark use, and it does not test optimizer or PLM-agent
behavior.

## Limitations

- PHOT remains the only empirical system in this calibration test.
- The formal comparator is the selected best available RUN-098 context row,
  not the full RUN-098 posterior.
- The failed `A+B` and `A+B+C` stages diagnose the sampled implementation and
  fast-fit search, not every possible global-readout or ratio-constraint
  parameterization.
- The raw and observable metrics have different meanings once the global
  readout is active. Both must be retained because either raw-overlarge or
  raw-collapsed solutions can produce superficially improved observable
  summaries.
- Full WHT/GFT coefficient tables were not copied into ANA-044; aggregate
  spectral summaries and output hashes were copied, and the coefficient tables
  remain on the runtime volume.

## Downstream use

- Do not use EXP-045/RUN-117 as evidence that the current cumulative fast-fit
  calibration protocol is full-ABC equivalent for PHOT.
- Cite RUN-117 as a negative calibration boundary: the `D` objective improves
  the best held-out distance and observable epistasis relative to RES-043, but
  still fails formal equivalence and does not repair accessibility or latent
  geometry.
- Keep the raw-versus-observable distinction in any later PHOT calibration
  record. A readout-mediated observable match is insufficient if the latent
  landscape either explodes or collapses.
