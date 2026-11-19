# RES-052: RUN-137 calibrated overlap finds Tier-1 contact but not public-model-specific coverage

## Summary

RUN-137 finds descriptor-space contact for Tier-1 empirical bands, but overlap is not specific to the public additive-core plus GNK composite because bare GNK and additive-core controls have comparable coverage rates.

## Generated from

- Analyses: ANA-052

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004, HYP-009
- Refutes: None
- Motivates: None

## Evidence

- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/source_tables
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/metrics.json
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_best_overlap_by_dataset.csv
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_best_overlap_by_dataset.md
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_control_overlap_failure_summary.csv
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_control_overlap_failure_summary.md
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_coverage_by_model_family.csv
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_coverage_by_model_family.md
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_output_collection_audit.csv
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_output_collection_audit.md
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_posterior_by_model_family.csv
- analyses/2026-06-19_ANA-052_run-137-stage-4c-calibrated-overlap-review/tables/run-137_posterior_by_model_family.md
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-137/outputs

## Interpretation

`RUN-137` completed the Stage 4C empirical coverage and calibrated-overlap test. It produced `44820` coverage rows and `7560` posterior knob-region rows over the `135` generated parameter settings from Stage 4B.

There is real descriptor-space contact with the empirical Tier-1 envelope. The best per-dataset overlap table reports that `10/10` Tier-1 combinatorial datasets have a best candidate with at least `0.6` of considered metrics inside band, and the median best fraction is `0.688`.

The control result prevents this from becoming a public-model-specific success. The all-row within-band rate is `0.448` for `full_public_observable`, `0.452` for `bare_gnk_texture`, and similar for additive-core-only controls. Best-overlap rows for the Tier-1 datasets are often bare-GNK texture rows rather than the intended public composite.

The sparse ProteinGym order-1 records also remain weak: `4` sparse ProteinGym datasets have best overlap fraction `0.0` under the recorded one-metric comparison. That is not a higher-order problem, because those assays were not used for topographic claims; it is a reminder that even the order-1 empirical panel is not broadly covered by this grid.

The correct interpretation is boundary evidence. The generated descriptor cloud can touch several empirical Tier-1 bands, but the current test does not distinguish the additive-core plus GNK mechanism from simpler controls.

## Effect on hypothesis

This weakens the current empirical-realism reading of `HYP-004`: overlap exists, but it is not model-specific.

It also weakens the `HYP-009` prerequisite that benchmark landscapes are matched enough to support downstream optimizer or PLM-agent comparisons. If controls overlap as well as the public model, downstream performance differences could be artifacts of the chosen descriptor envelope rather than evidence about routing or biological realism.

## Limitations

- Best-candidate overlap is permissive and does not require a single robust knob region to explain multiple assays.
- Coverage is measured in descriptor space, not in sequence-level fitness prediction.
- The sparse ProteinGym rows contribute only order-1 comparisons.
- Comparable control coverage prevents mechanistic attribution to the public composite.

## Downstream use

- Do not cite `RUN-137` as validation that the public generator covers empirical landscapes in a mechanistically specific way.
- Use the posterior and coverage tables to identify which descriptors and controls drive the overlap.
- Treat the result as input to the harder Stage 4D and 4E gates, not as a standalone success.
