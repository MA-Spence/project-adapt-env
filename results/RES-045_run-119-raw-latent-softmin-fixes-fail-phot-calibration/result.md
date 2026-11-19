# RES-045: RUN-119 raw-latent soft-min fixes fail PHOT calibration

## Summary

EXP-046/RUN-119 produced 0/24 operationally equivalent rows. The low-dimensional soft-min stage partially improved raw accessibility but left raw scale hundreds-fold too large and held-out distance 1.2240; raw-latent constraints and up-weighted accessibility worsened the fit and collapsed observable variation.

## Generated from

- Analyses: ANA-045

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/metrics.json
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_stage_equivalence_decision.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_stage_best_test_distances.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_bounded_gate_test_distances.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_abc_comparator_test_distances.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_top_component_distances.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_observable_vs_raw_epistasis_and_landscape.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_raw_latent_pathology_diagnostics.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_stage_failure_modes.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_best_candidate_parameters.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/run-119_output_collection_audit.md
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/source_tables/run-119_bounded_gate_equivalence.csv
- analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/source_tables/run-119_graphlet_summary.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-119/outputs

## Interpretation

`RUN-119` directly tests the three path-forward points recorded after
`RES-044`: replacing the random specific-`J` ruggedness source with a
low-dimensional two-additive-trait soft-min latent generator, adding raw-latent
constraints plus a restricted monotone-sigmoid readout, and then up-weighting
raw accessibility as a first-class target.

The experiment does not validate that path. The formal equivalence decision
produced `0/24` operationally equivalent rows. All selected rows failed the
absolute and relative aggregate-distance margins against the best available
`RUN-098` full-ABC context row, the `empirical_readout_residual` comparator at
test distance `0.2115`.

Stage `1` was the best EXP-046 stage but remained poor. Its best held-out
distance was `1.2240`, much worse than the best `RUN-117` A-through-D row
from `RES-044` (`0.4122`) and worse than the nested full-ABC context row
(`0.4826`). It partially moved raw topology in the intended direction, with
raw accessibility `0.0796` versus empirical `0.1975`, but raw local maxima
remained low (`2.33` versus `5.17`) and raw scale was pathological: raw
single-mutant SD was about `934x` empirical and raw magnitude epistasis about
`424x` empirical.

Stage `1+2` worsened the fit. Adding raw-latent constraints and the restricted
sigmoid readout produced best held-out distance `1.6872`, with no component
guardrail passes among selected rows. Raw accessibility rose only to `0.0958`
(`49%` of empirical), raw local maxima fell to `1.58` (`31%` of empirical),
and raw scale remained tens-fold too large. Observable variation collapsed:
observable single-mutant SD was only `4.9%` of empirical and observable
magnitude epistasis only `2.3%` of empirical.

Stage `1+2+3` did not rescue the pathology. The best selected held-out row was
effectively identical to stage `1+2`, with the same test distance (`1.6872`),
raw accessibility (`0.0958`), raw local maxima (`1.58`), and observable-scale
collapse. In this sampled fast-fit search, increasing the raw-accessibility
weight did not find a usable raw-latent PHOT geometry.

## Effect on hypothesis

This further weakens the current PHOT adequacy form of `HYP-004`. The current
procedural landscape route still does not produce a PHOT landscape matched
across accessibility, local optima, epistasis scale, spectra, and local-energy
summaries.

For `HYP-009`, this is methodological motivation rather than direct support or
refutation. The run does not test optimizer or PLM-agent behavior, but it does
show that the tested raw-latent soft-min repair is not a validated
matched-landscape construction route for benchmark use.

## Limitations

- PHOT remains the only empirical system in this calibration test.
- The formal comparator is the selected best available `RUN-098` context row,
  not the full `RUN-098` posterior.
- The failed `1`, `1+2`, and `1+2+3` stages diagnose this sampled
  low-dimensional soft-min implementation and fast-fit search, not every
  possible latent-topology model.
- The raw and observable metrics have different meanings once the monotone
  readout is active. Both are retained because an apparent observable-scale
  improvement can still be unusable if latent topology or raw scale is wrong.
- Full WHT/GFT coefficient tables were not copied into `ANA-045`; aggregate
  spectral summaries and output hashes were copied, and the coefficient tables
  remain on the runtime volume.

## Downstream use

- Do not use EXP-046/RUN-119 as evidence that the raw-latent soft-min repair
  solves the PHOT calibration pathology.
- Cite RUN-119 as a negative boundary on the tested RES-044 path-forward
  points: raw accessibility can be partially moved, but not together with
  held-out graphlet equivalence and correct raw/observable scale under this
  protocol.
- Keep the raw-versus-observable distinction in later PHOT calibration work.
  Neither readout-mediated observable matching nor direct raw-accessibility
  pressure is sufficient unless latent topology and raw scale are both matched.
