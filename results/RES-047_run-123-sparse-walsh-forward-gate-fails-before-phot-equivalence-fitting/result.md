# RES-047: RUN-123 sparse Walsh forward gate fails before PHOT equivalence fitting

## Summary

EXP-048/RUN-123 did not validate Generator A: the empirical spectrum calibration put all normalized sparse-Walsh weight on order 3, and the H1 forward sparse-Walsh probe produced 0/12 joint passes, failing accessibility, magnitude-epistasis ratio, and Walsh high-order fraction before any sequence-conditioned PHOT equivalence fit was allowed.

## Generated from

- Analyses: ANA-047

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/metrics.json
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/run-123_sparse_walsh_spectrum_calibration.md
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/run-123_forward_probe_metric_gate.md
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/run-123_forward_probe_candidate_gate.md
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/run-123_stage_gate_and_fit_decision.md
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/run-123_comparator_context.md
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/run-123_graphlet_and_module_catalog_summary.md
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/tables/run-123_output_collection_audit.md
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_summary.json
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_preflight.json
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_remote_output_manifest.tsv
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_sparse_walsh_spectrum_calibration.csv
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_sparse_walsh_forward_probe.csv
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_bounded_gate_equivalence.csv
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_abc_comparator_selected_family_distances.csv
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_graphlet_catalog.csv
- analyses/2026-06-17_ANA-047_run-123-phot-sparse-walsh-potts-spectrum-calibrated-generator-review/source_tables/run-123_module_catalog.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-123/outputs

## Interpretation

`RUN-123` directly tests the Generator A branch from
`scratch/RES-046_scalable_multipeak_generator_design.md`: additive stability
core plus sparse hashed high-order Potts/Walsh terms, with empirical
variance-by-order calibration and a forward H1 capability gate before any
sequence-conditioned fit.

The run completed successfully as an execution (`scheduler_exit_code = 0`) but
failed scientifically at the preregistered forward gate. The G-stage spectrum
calibration found that only order `3` passed the implemented noise-floor flag
among orders `2` to `5`. Normalized global sparse-Walsh weights were therefore
`0.0, 1.0, 0.0, 0.0` for orders `2, 3, 4, 5`. This means the empirical PHOT
calibration graphlets did not support a broad high-order sparse-Walsh mixture
under the EXP-048 spectrum-calibration rule.

The H1 forward generative probe then evaluated `12` sparse-Walsh candidates
across the `epi_weight` grid and high-order multiplier grid. All candidates
failed the joint gate. Each candidate passed `4/7` metrics, but every candidate
failed:

- `accessible_path_fraction`: empirical mean `0.2319`; generated central
  interval envelope `0.0` to `0.0715`.
- `magnitude_epistasis_to_single_mutant_sd_ratio`: empirical mean `1.0161`;
  generated interval envelope `1.0253` to `1.2692`.
- `walsh_high_order_fraction`: empirical mean `0.3433`; generated interval
  envelope `0.4479` to `0.7659`.

The metrics that passed consistently were local maxima, local sign epistasis,
local reciprocal sign epistasis, and diminishing-returns epistasis. Thus the
operator generated some ruggedness-like local structure, but did not jointly
match accessibility, scale-normalized epistasis, and the spectral envelope.

Because `sparse_walsh_forward_probe_passed = false`, the protocol correctly
skipped the `G+H1+FIT` sequence-conditioned numerical fast-fit stage. `RUN-123`
therefore produced `0` candidate evaluations, `0` selected candidates, and `0`
operational-equivalence rows. This is the intended gate behavior, not a missing
output problem.

## Effect on hypothesis

This further weakens the current PHOT adequacy form of `HYP-004`. EXP-048 was
designed to test whether the sparse Walsh/Potts operator could clear a forward
generative multi-peak/high-order capability gate before being trusted in the
sequence-conditioned PHOT equivalence workflow. It did not clear that gate.

For `HYP-009`, this is methodological motivation rather than support. The run
does not test optimizer or PLM-agent behavior, and it does not provide a
matched benchmark landscape. It does provide a negative boundary on this
specific Generator A construction: sparse high-order capacity alone is not
sufficient under the implemented spectrum calibration and graphlet-aligned
support procedure.

## Limitations

- The result is specific to the EXP-048 Generator A implementation,
  graphlet-aligned support pool, spectrum-calibration rule, H1 grid, and PHOT
  graphlet summaries.
- The H1 gate intentionally stopped the sequence-conditioned fast fit, so
  there is no sparse-Walsh fitted-candidate distance to compare against
  `RUN-098`. The copied `RUN-098` distances are context only.
- The run does not evaluate the FGM/H2 branch from the scratch design note.
- PHOT remains a single-assay diagnostic. This result should not be generalized
  to all sparse Walsh/Potts landscape generators or all empirical proteins.
- The order-5 empirical local fraction was nonzero, but it did not pass the
  preregistered noise-floor decision rule. The result should therefore be read
  through the registered pass flag, not by informal visual inspection of one
  spectrum column.

## Downstream use

- Do not cite EXP-048/RUN-123 as evidence that Generator A solves PHOT
  calibration or improves over `RUN-098`; no FIT candidate was generated.
- Cite RUN-123 as a negative forward-gate boundary: the implemented
  spectrum-calibrated sparse Walsh/Potts field fails accessibility,
  magnitude-epistasis-ratio, and Walsh high-order-fraction consistency before
  sequence-conditioned fitting.
- Preserve the distinction between execution success and scientific gate
  failure. The run completed successfully, but the model did not pass the
  prospective capability test.
- Keep FGM/H2 separate in any future record; EXP-048 deliberately did not test
  it.
