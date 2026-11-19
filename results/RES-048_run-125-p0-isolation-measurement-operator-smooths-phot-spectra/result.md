# RES-048: RUN-125 P0 isolation shows PHOT measurement operator smooths injected spectra

## Summary

RUN-125 supports the PHOT target/measurement-operator mismatch diagnosis: under the identical MAVE-NN plus dense-graphlet operator, cached/refit targets are smooth/additive-like, known injected spectra are recovered as unimodal and nearly epistasis-free, and sparse-Walsh supports are well contained in sampled graphlets.

## Generated from

- Analyses: ANA-048

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/source_tables
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/metrics.json
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_output_collection_audit.csv
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_p0_1_target_stability.csv
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_p0_1_mavenn_refit_fit_quality.csv
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_p0_2_roundtrip_transfer_key_metrics.csv
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_p0_2_profile_summary.csv
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_p0_3_support_alignment.csv
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_graphlet_catalog_summary.csv
- analyses/2026-06-18_ANA-048_run-125-phot-p0-target-operator-isolation-review/tables/run-125_phase_decisions.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-125/outputs

## Interpretation

`RUN-125` completed the three P0 isolation checks requested after `RES-047`.
The run used the PHOT assay `PHOT_CHLRE_Chen_2023`, `24` dense graphlets,
`10` P0.1 MAVE-NN refits, `36` P0.2 transfer rows, and `2048` P0.3 sparse
Walsh support terms. The analysis copied the compact output bundle into
`ANA-048` and derived the decision tables listed above.

The result supports the broad Layer-2 diagnosis in
`scratch/RES-047_miscalibration_root_cause_and_isolation_plan.md`: the target
summaries, MAVE-NN latent representation, and dense-graphlet measurement
operator are not on a common ruggedness scale. However, the observed direction
is important. The known-spectrum P0.2 assay did not show that injected order-3
signal is recovered as a high high-order fraction under the current operator.
Instead, the identical MAVE-NN plus `dense_descriptive` measurement operator
smoothed the injected spectra almost completely.

For P0.1, the cached and refit target summaries were smooth and
additive-like. The cached target mean high-order fraction was
`3.8811e-13`, and refit high-order means remained below `1.0594e-12`.
Cached local maxima averaged `1.0`. Refit accessibility was unstable across
regularisation profiles, ranging from `0.0` to `0.5417`, so accessibility
should not be treated as a stable target under the current latent-operator
configuration.

For P0.2, all three injected profiles were recovered as nearly unimodal and
nearly epistasis-free after MAVE-NN fitting and dense-graphlet measurement.
The maximum direct injected high-order fraction was `0.01393`, but the
maximum recovered high-order fraction was only `1.0091e-14`. The maximum
recovered magnitude-epistasis ratio was `4.0040e-07`, recovered local sign
epistasis was `0.0`, and recovered local maxima were exactly `1.0` for all
profiles. This shows that the current measurement path is strongly
non-identity for the ruggedness summaries being used.

For P0.3, low graphlet containment of the sparse-Walsh supports is not the
limiting explanation for the `RES-047` failure under the audited support rule.
All `2048` order-3 support terms were contained in at least one sampled
graphlet, `95.21%` were contained in calibration graphlets, and `92.92%` were
contained in test graphlets. The sampled dense graphlets therefore cover the
EXP-048 support construction well enough that the dominant failure should be
attributed to the target/measurement operator rather than to sparse support
visibility.

## Effect on hypothesis

This further weakens the current PHOT adequacy route under `HYP-004`. The
experiment was not a generator repair and does not test whether a different
landscape family could match PHOT. It shows that the current PHOT target
summaries are themselves poorly conditioned for the intended high-order and
local-ruggedness inference: cached/refit empirical targets are smooth, and
known injected spectra are suppressed by the same measurement path used to
judge generated candidates.

For `HYP-009`, this is methodological motivation rather than support. RUN-125
does not test optimizer or PLM-agent behavior, and it does not produce a
matched benchmark landscape. It does show that future benchmark-facing
construction cannot rely on the current scalar high-order/local-ruggedness
target unless the measurement operator is corrected, replaced, or explicitly
calibrated against known spectra.

## Limitations

- The result is PHOT-specific and measurement-path-specific. It should not be
  generalized to all proteins, all MAVE-NN settings, or all high-order
  landscape generators.
- P0.2 uses synthetic injected spectra, so it diagnoses the measurement
  operator rather than biological truth in the PHOT assay.
- P0.3 audits support-to-graphlet containment under the implemented EXP-048
  support source and ANA-047 module catalog. It does not prove that every
  possible support sampler is adequate.
- P0.1 accessibility is refit-sensitive, so accessibility should be read as an
  unstable target in this configuration rather than as a clean empirical
  constant.

## Downstream use

- Cite RUN-125 as evidence that the current PHOT MAVE-NN plus dense-graphlet
  operator smooths or suppresses injected high-order/local epistasis.
- Do not cite RUN-125 as support for Generator A, a repaired PHOT fit, or a
  matched benchmark landscape.
- Do not use low sparse-Walsh support containment as the primary explanation
  for the EXP-048/RES-047 failure without additional evidence; RUN-125 argues
  against that explanation for the audited support rule.
- Keep the result tied to `ANA-048` tables when interpreting the numerical
  thresholds, especially `run-125_phase_decisions.csv` and
  `run-125_p0_2_profile_summary.csv`.
