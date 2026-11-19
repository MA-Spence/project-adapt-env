# RES-046: RUN-121 calibration isolation improves raw scale but fails PHOT equivalence

## Summary

EXP-047/RUN-121 improved over the RUN-119 raw-latent soft-min failure but still produced 0/32 operationally equivalent rows; the best held-out stage was A+F+B+C+D at distance 0.9155, while observable variation collapsed and all selected component-equivalence decisions failed.

## Generated from

- Analyses: ANA-046

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/metrics.json
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_abc_comparator_test_distances.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_best_candidate_parameters.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_bounded_gate_test_distances.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_graphlet_and_module_catalog_summary.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_observable_vs_raw_epistasis_and_landscape.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_output_collection_audit.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_prior_run_comparison.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_raw_latent_pathology_diagnostics.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_softmin_capability_probe_summary.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_stage_best_test_distances.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_stage_equivalence_decision.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_stage_failure_modes.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/tables/run-121_top_component_distances.md
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/source_tables/run-121_bounded_gate_equivalence.csv
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/source_tables/run-121_graphlet_summary.csv
- analyses/2026-06-17_ANA-046_run-121-phot-res045-calibration-isolation-review/source_tables/run-121_softmin_capability_probe.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-121/outputs

## Interpretation

`RUN-121` directly tests the `RES-045` calibration-isolation plan as a single
cumulative PHOT-only diagnostic. The stage order separates generator-scale
anchoring (`A`), scale-invariant objective conditioning (`A+F`), model-free
k-trait soft-min capability (`A+F+B`), sequence-conditioned k/texture fitting
(`A+F+B+C`), structured residual plus restricted sigmoid readout
(`A+F+B+C+D`), and expanded search (`A+F+B+C+D+E`).

The experiment does not validate this path as a PHOT calibration repair. The
formal equivalence decision produced `0/32` operationally equivalent rows. No
selected candidate passed the absolute or relative aggregate-distance margins
against the best available `RUN-098` full-ABC comparator, the
`empirical_readout_residual` row with test distance `0.2115`.

The best `RUN-121` held-out row was `A+F+B+C+D`, with test distance `0.9155`.
That improves on the best `RUN-119` raw-latent soft-min row from `RES-045`
(`1.2240`) by about `25%`, but it is still much worse than the best `RUN-117`
cumulative calibration row from `RES-044` (`0.4122`) and the best nested
full-ABC context row (`0.4826`). The expanded-search stage improved the
calibration objective to `0.4884`, but its best held-out row was slightly worse
than the unexpanded structured-residual stage (`0.9298` versus `0.9155`).

The stage-level pattern is informative. `A` and `A+F` stayed poor
(`1.7127` and `1.6716` best held-out distance), with no component guardrail
passes and failures dominated by plaquette-energy, sign-epistasis,
single-mutant-SD, and magnitude-epistasis components. The probe-only `A+F+B`
stage showed that the model-free additive k-trait soft-min family can fit
individual graphlets only partially: median fit correlations were about
`0.86` to `0.88`, but local maxima reached only `19%` to `35%` of the target
and magnitude epistasis only `14%` to `30%` of the target.

Later stages partly fix the specific raw-scale pathology seen in `RES-045`.
For `A+F+B+C`, `A+F+B+C+D`, and `A+F+B+C+D+E`, raw single-mutant SD and raw
magnitude epistasis are close to empirical scale (`0.97x` to `1.25x` and
`0.93x` to `1.02x`, respectively). But this does not produce a usable PHOT
landscape. The best stage reaches only `32%` of empirical raw accessibility,
while the expanded-search stage reaches `85%` of empirical raw local maxima
only with raw accessibility collapsed to `4%` of target.

The strongest remaining failure is observable readout collapse. The
stage-best rows have observable single-mutant SD only `0.08%` to `0.16%` of
empirical and observable magnitude epistasis only `0.04%` to `0.14%` of
empirical. Thus the tested k-trait/structured-residual/sigmoid family can bring
raw latent scale near the empirical range only while flattening the observable
fitness landscape, which is incompatible with the PHOT graphlet-calibration
goal.

## Effect on hypothesis

This further weakens the current PHOT adequacy form of `HYP-004`. The current
procedural route still does not match PHOT across held-out graphlet distance,
component guardrails, accessibility, local maxima, and observable epistasis
scale. `RUN-121` is less pathological than `RUN-119` on raw scale, but it does
not restore empirical observable variation or pass operational equivalence.

For `HYP-009`, this is methodological motivation rather than direct support or
refutation. The run does not test optimizer or PLM-agent behavior. It does
show that this tested calibration-isolation route should not be used as a
matched-landscape construction route for benchmark-agent comparisons.

## Limitations

- PHOT remains the only empirical system in this calibration-isolation test.
- The formal comparator is the selected best available `RUN-098` context row,
  not the full `RUN-098` posterior distribution.
- Stage `A+F+B` is intentionally probe-only; it diagnoses additive soft-min
  capability before sequence-conditioned fitting, but it does not itself have
  held-out equivalence rows.
- The failed result is specific to the tested cumulative implementation,
  candidate budget, and fast-fit search. It does not refute every possible
  latent-topology or structured-residual model.
- Full WHT/GFT coefficient tables were not copied into `ANA-046`; aggregate
  spectral summaries and SHA-256 hashes were copied, and the coefficient
  tables remain on the runtime volume.

## Downstream use

- Do not use EXP-047/RUN-121 as evidence that the RES-045
  calibration-isolation plan solves PHOT calibration.
- Cite RUN-121 as a negative boundary on the tested k-trait
  soft-min/structured-residual/sigmoid route: raw scale can be brought near
  empirical range, but only with failed held-out equivalence and collapsed
  observable variation.
- Preserve the raw-versus-observable distinction in later PHOT work. A raw
  latent-scale repair is not sufficient unless the observable landscape retains
  empirical-scale single-mutant and epistatic variation.
