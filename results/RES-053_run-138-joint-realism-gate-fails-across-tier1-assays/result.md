# RES-053: RUN-138 joint realism gate fails across Tier-1 assays

## Summary

RUN-138 fails the Stage 4D hard joint-realism gate: 0 joint passes were observed, the maximum assays hit by any setting was 1 rather than the required 2, and all hits were confined to CR9114-H3.

## Generated from

- Analyses: ANA-053

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004, HYP-009
- Refutes: None
- Motivates: None

## Evidence

- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/source_tables
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/metrics.json
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_coverage_by_dataset.csv
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_coverage_by_dataset.md
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_coverage_by_metric.csv
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_coverage_by_metric.md
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_hit_dataset_counts.csv
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_hit_dataset_counts.md
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_joint_candidate_summary.csv
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_joint_candidate_summary.md
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_output_collection_audit.csv
- analyses/2026-06-19_ANA-053_run-138-stage-4d-joint-realism-review/tables/run-138_output_collection_audit.md
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-138/outputs

## Interpretation

`RUN-138` completed the Stage 4D hard joint-realism gate. The preregistered gate required a single knob setting to hit at least `0.6` of the considered metrics across at least `2` Tier-1 assays.

The gate failed. The run produced `103` candidate rows, but `0` joint passes. The maximum number of assays hit by any setting was `1`, below the required minimum of `2`. All candidate hits were confined to `cr9114_h3_phillips2021_graphfla`.

This is stronger negative evidence than Stage 4C. Stage 4C allowed per-dataset best overlap and showed that descriptor-space contact exists. Stage 4D asks whether the same knob setting can jointly match a multi-descriptor, multi-assay package. The current grid cannot.

Descriptor-level coverage shows why the failure matters. Accessibility and reciprocal-sign epistasis are easier to bring into band, but fitness scale, random-walk autocorrelation, and broad multi-assay matching are not jointly controlled. The current public model therefore cannot yet provide a matched empirical benchmark manifold for downstream optimizer or PLM-agent experiments.

## Effect on hypothesis

This weakens `HYP-004` because the current implementation does not reproduce the intended joint empirical package of ruggedness, accessibility, local/global structure, and DFE-related summaries across multiple Tier-1 assays.

It weakens `HYP-009` because downstream optimizer or PLM-guided benchmark comparisons require matched landscapes. A setting that hits only one assay does not satisfy that prerequisite.

## Limitations

- The result is a failure of the registered Stage 4D grid and thresholds, not a formal refutation of every possible additive-core plus GNK parameterization.
- The assay panel is Tier-1 and operator-free, but still finite.
- The descriptor gate is stringent by design and does not test sequence-level prediction.
- The run does not identify whether failure is caused mainly by support routing, parameter-grid resolution, readout scale, or descriptor choice.

## Downstream use

- Do not use the current Stage 4D grid as a matched benchmark landscape set.
- Do not proceed from this grid to optimizer or PLM-agent claims that require empirical matching.
- Treat CR9114-H3 as a single-assay contact case, not as evidence of broad joint realism.
