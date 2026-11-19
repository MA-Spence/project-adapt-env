# RES-042: RUN-112 calibrated cheap PHOT fit fails component-equivalence despite readout rescaling

## Summary

RUN-112 completed EXP-043 and found zero operationally equivalent calibrated cheap-fitting rows. Positive affine observable calibration reduced visible epistasis scale for selected candidates but did not match empirical epistasis amplitude, left raw latent magnitude epistasis tens to hundreds of units too large, and failed component guardrails dominated by signed pairwise epistasis and accessibility/local-summary mismatch. The protocol is therefore not an equivalent replacement for full ABC or a validated scale-up path for GB1/TrpB.

## Generated from

- Analyses: ANA-042

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/metrics.json
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/run-112_equivalence_by_method_reference.md
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/run-112_same_family_equivalence.md
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/run-112_calibrated_test_distances.md
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/run-112_top_component_failures.md
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/run-112_observable_vs_raw_epistasis_and_locality.md
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/run-112_best_epistasis_amplitude_rows.md
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/tables/run-112_empirical_module_locality_test.md
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/source_tables/run-112_fit_protocol_equivalence.csv
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/source_tables/run-112_fit_protocol_distances.csv
- analyses/2026-06-16_ANA-042_run-112-phot-epistasis-calibrated-cheap-fitting-review/source_tables/run-112_graphlet_summary.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-112/outputs

## Interpretation

`RUN-112` completed the intended `EXP-043` test: the latest PHOT nested-min
family was fitted with cheap protocols, a positive affine observable readout
was fitted on calibration graphlets, epistasis-amplitude terms were included
in the objective, and selected `RUN-098` full-ABC particles were re-evaluated
through the calibrated reporting path.

The experiment did not validate the calibrated cheap protocol. It produced
`0/40` operationally equivalent rows. This is stricter than the `EXP-042`
result: aggregate distance margins alone were often satisfied, but every
method failed the component guardrail package. For same-family comparisons,
`exact_descriptive_gated`, `local_abc_after_cheap_fit`, and
`surrogate_assisted` each passed the absolute and relative aggregate margins
for `5/5` families, while `deterministic_summary_loss` passed `4/5`; all four
methods had `0/5` same-family component-guardrail passes.

The dominant failing component was `signed_pairwise_epistasis_mean`. The
empirical signed mean is essentially zero on test graphlets, so small absolute
signed offsets generate very large scaled component distances. This means the
EXP-043 calibrated equivalence scale is not directly comparable to the older
EXP-042 same-family graphlet-distance scale.

The best calibrated-objective selected row was
`exact_descriptive_gated` / `uncoupled_module`, with test distance `1.1474`.
It still violated `magnitude_epistasis_mean` and
`signed_pairwise_epistasis_mean`, had calibrated magnitude epistasis mean
`0.0156` versus empirical `0.2926`, raw latent magnitude epistasis mean
`29.10`, and accessible-path fraction `0.0` versus empirical `0.1975`.

The row closest to the empirical calibrated magnitude epistasis mean was
`exact_descriptive_gated` / `empirical_module`, at `0.1774` versus empirical
`0.2926`. That row still retained raw latent magnitude epistasis `131.9`,
failed the signed-pairwise component, and had test distance `12.0043`.
Therefore the positive affine readout shrank observable epistasis scale but
did not solve the empirical amplitude mismatch, and it did not repair the raw
latent epistasis pathology.

Some local-energy summaries moved in the right direction. The best
calibrated-objective row had test edge top-10 energy fraction `0.3955`, close
to empirical `0.4060`. That local-energy improvement was not enough because
accessibility collapsed, site top-1 energy remained too low, and epistasis
component failures dominated the decision.

The empirical-module specificity warning remains. Empirical-module candidates
remain diffuse over modules on test graphlets: the selected `RUN-098`
empirical-module comparator has `11.72/12` effective edge modules and edge
normalized entropy `0.995`, and EXP-043 cheap empirical-module rows remain
similarly diffuse.

## Effect on hypothesis

This weakens the current PHOT adequacy form of `HYP-004`: after adding
observable-scale calibration and explicit epistasis-amplitude terms, the model
still does not provide a calibrated PHOT landscape suitable for interpreting
accessibility, neutral networks, and epistasis as empirical-like biological
features.

For `HYP-009`, this is methodological motivation rather than support. The run
shows that benchmark-facing matched landscape construction needs stricter and
better-conditioned matched-summary controls before optimizer or PLM-agent
comparisons. It does not test optimizer behaviour and does not support a
biologically routed benchmark effect.

## Limitations

- PHOT is still the only empirical system in this calibration test.
- The selected full-ABC comparator is represented by selected `RUN-098`
  particles, not by the full posterior distribution.
- The signed-pairwise epistasis mean is numerically delicate because the
  empirical mean is near zero; the component failure is real under the
  preregistered scaling, but it should be interpreted as a conditioning issue
  as well as a model-fit issue.
- The positive affine readout is a simple observation-scale correction. It
  cannot by itself demonstrate that the latent nested-min geometry is
  empirically calibrated.
- Large WHT/GFT coefficient tables and full graphlet-level metrics were not
  copied into the analysis; they remain on the runtime volume and are listed
  in the output audit.

## Downstream use

- Do not cite EXP-043 as evidence that the cheap calibrated protocol is
  equivalent to the full ABC fit.
- Do not use this exact EXP-043 protocol as the default GB1/TrpB scale-up
  route without revising the model or objective.
- Preserve the useful diagnostic distinction in future work: report both
  calibrated observable epistasis and raw latent epistasis, because readout
  rescaling can hide but not fix large latent curvature.
- Keep `surrogate_assisted` and `exact_descriptive_gated` from RES-041 as
  cheap search candidates only for same-family graphlet behaviour, not as
  validated protocols for empirical epistasis/amplitude matching.
