# RES-054: RUN-139 held-out descriptor generalisation is not validated

## Summary

RUN-139 does not validate held-out descriptor generalisation: median held-out fraction in band was 0.5, only 2/10 held-out datasets exceeded a strict majority and reached 0.6, and selected parameters collapsed to a narrow random setting.

## Generated from

- Analyses: ANA-054

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004, HYP-009
- Refutes: None
- Motivates: None

## Evidence

- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/source_tables
- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/tables/metrics.json
- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/tables/run-139_heldout_predictions_ranked.csv
- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/tables/run-139_heldout_predictions_ranked.md
- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/tables/run-139_output_collection_audit.csv
- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/tables/run-139_output_collection_audit.md
- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/tables/run-139_selected_parameter_counts.csv
- analyses/2026-06-19_ANA-054_run-139-stage-4e-heldout-generalisation-review/tables/run-139_selected_parameter_counts.md
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-139/outputs

## Interpretation

`RUN-139` completed the Stage 4E leave-one-landscape-out descriptor generalisation test across the `10` Tier-1 empirical landscapes.

The result does not validate held-out generalisation. The median held-out fraction of metrics inside empirical band was `0.500`. `6/10` held-out datasets reached at least `0.5`, but only `2/10` exceeded a strict majority and only `2/10` reached `0.6`.

The selected parameter region is also not convincing. `9/10` held-out fits selected `random__sm6.000__ew1.000__seed703`, and the remaining fit selected a closely related low-margin shuffled setting. That collapse suggests the procedure found a generic compromise point rather than a calibrated landscape-specific manifold.

The dataset-level pattern is uneven. CR9114-H3 and CR6261-H9 reached `0.625` fraction in band. TrpB, CR9114-H1, and CR6261-H1 reached only `0.375`. CR9114-FluB had a very large mean-z distance despite a nonzero fraction-in-band score, showing that at least one descriptor block remains severely mismatched.

This is mixed contact but not validation. It provides evidence that a few descriptors can be brought into range for some held-out landscapes, but it fails the planned majority-descriptor generalisation criterion.

## Effect on hypothesis

This weakens `HYP-004` because the current public generator does not predict unseen Tier-1 empirical descriptor packages with the planned reliability.

It also weakens `HYP-009` because benchmark-routing claims depend on held-out descriptor matching. Without that, downstream agent comparisons would be under-matched.

## Limitations

- The held-out test uses descriptor bands, not direct sequence-level fitness prediction.
- The grid may be too coarse to detect a better calibrated region, but the registered grid still failed the planned criterion.
- Fraction-in-band can hide severe scale mismatch, as seen in the CR9114-FluB mean-z distance.
- The result should not be interpreted as a general refutation of all future public-generator variants.

## Downstream use

- Do not cite the Stage 4 grid as held-out validated.
- Do not use these settings as finalized empirical benchmark landscapes for optimizer or PLM-agent comparisons.
- Use the held-out table to identify which assays and descriptor blocks drive the failed generalisation.
