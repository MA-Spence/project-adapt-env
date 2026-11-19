# RES-056: RUN-129 non-additive PHOT operators fail certification

## Summary

RUN-129 completed EXP-050 but certified no non-additive MAVE-NN PHOT measurement operator: pairwise profiles were preflight-skipped at PHOT length 118, black-box/GE profiles failed under the installed MAVE-NN/Keras build interface, and 0/24 certification checks passed.

## Generated from

- Analyses: ANA-056

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004, HYP-006, HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/metrics.json
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_certification_metric_matrix.csv
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_certification_metric_matrix.md
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_fit_failure_modes.csv
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_fit_failure_modes.md
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_graphlet_catalog_summary.csv
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_graphlet_catalog_summary.md
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_injected_spectrum_summary.csv
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_injected_spectrum_summary.md
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_metric_pass_summary.csv
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_metric_pass_summary.md
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_operator_certification_status.csv
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_operator_certification_status.md
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_output_collection_audit.csv
- analyses/2026-06-22_ANA-056_run-129-phot-nonadditive-measurement-operator-certification-review/tables/run-129_output_collection_audit.md

## Interpretation

`RUN-129` completed the corrected `EXP-050` run and produced all expected
output tables. The run used `167,529` PHOT assay rows from `RUN-077`, selected
`24` dense complete graphlets, injected `176` known-spectrum terms, and tested
three spectrum profiles across two non-additive MAVE-NN operator profiles.

The certification result is negative. Neither operator profile was certified:
`0/2` operator profiles passed, `0/24` required operator-spectrum-metric checks
passed, and `0/24` certification-metric rows had recovered values. The absence
of recovered values is important. This is not evidence that the non-additive
operators recover the wrong epistasis spectrum; the current implementation did
not produce recovered non-additive landscapes for the certification metrics.

The two failure channels are explicit in `ANA-056`. The pairwise profile was
preflight-skipped for all `3/3` injected spectra because the PHOT sequence
length was `118`, exceeding the preregistered
`max_pairwise_sequence_length = 80`. The black-box/GE profile failed for all
`3/3` spectra with a `MultilayerPerceptronGPMap` / Keras `build()` signature
error in the installed MAVE-NN interface.

The Slurm job itself completed cleanly: `labproj-exit-code.txt` records `0`,
and the only log messages are TensorFlow CPU/GPU informational warnings. The
result is therefore a valid negative operator-certification record rather than
a failed scheduler run.

## Effect on hypothesis

This weakens `HYP-006` directly. The current non-additive MAVE-NN measurement
route is not certified as a support-aware, observation-aware operator for PHOT
dense graphlet recovery under the installed stack.

It weakens the current `HYP-004` route because the PHOT calibration/validation
series still lacks a certified operator that preserves injected
high-order/local-ruggedness spectra. It does not refute the biological claim
that real PHOT landscapes may contain such structure.

It weakens the current `HYP-009` benchmark-construction route because the
planned empirical/model comparison cannot inherit a certified non-additive
PHOT measurement operator from `EXP-050`. It does not test PLM-guided agents,
optimizer behaviour, or the additive-core + GNK forward generator.

## Limitations

- The result is a measurement-operator certification failure, not a biological
  PHOT ruggedness measurement.
- The pairwise profile did not run to fit; it was preflight-skipped by the
  registered sequence-length guard.
- The black-box/GE profile failure is a technical MAVE-NN/Keras interface
  failure, not a demonstrated statistical failure of a successfully fit
  black-box operator.
- No recovered landscapes exist for the certification metrics, so no
  non-additive transfer function should be estimated from this run.
- `RUN-126` and `RUN-128` are execution context only; the result is based on
  the completed corrected run, `RUN-129`.

## Downstream use

- Cite this result as showing that `EXP-050` did not certify a non-additive
  PHOT measurement operator.
- Do not cite it as evidence that non-additive MAVE-NN operators recover,
  smooth, or distort high-order spectra under PHOT.
- Do not use the two tested non-additive profiles from `EXP-050` as certified
  measurement operators for empirical/model benchmark comparisons.
