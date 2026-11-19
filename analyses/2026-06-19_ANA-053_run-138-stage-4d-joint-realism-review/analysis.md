# ANA-053: RUN-138 Stage 4D joint realism gate review

## Purpose

Collect and interpret `RUN-138` from `EXP-055`, the Stage 4D hard joint-realism gate.

This analysis asks whether any single knob setting can simultaneously hit the required multi-descriptor package across at least `2` Tier-1 assays, using the preregistered minimum fraction `0.6` of metrics inside empirical bands.

## Linked experiments/runs

- Experiments: `EXP-055`
- Runs: `RUN-138`

## Produced artifacts

- `source_tables/run-138_summary.json`
- `source_tables/run-138_preflight.json`
- `source_tables/run-138_descriptor_scope_audit.csv`
- `source_tables/run-138_posterior_knob_regions.csv`
- `source_tables/run-138_joint_realism_coverage.csv`
- `source_tables/run-138_joint_realism_candidates.csv`
- `tables/metrics.json`
- `tables/run-138_output_collection_audit.md`
- `tables/run-138_joint_candidate_summary.md`
- `tables/run-138_hit_dataset_counts.md`
- `tables/run-138_coverage_by_dataset.md`
- `tables/run-138_coverage_by_metric.md`

## Collection policy

`RUN-138` completed on `lab-slurm` as Slurm job `1945` with inferred scheduler exit code `0`. The compact runtime output bundle was collected and DVC-tracked at:

`data/processed/stage-4d-joint-multidescriptor-realism-gate/RUN-138.dvc`

The original runtime output remains at `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-138/outputs`.

## Joint realism decision

The hard gate failed. `RUN-138` produced `103` candidate rows, but `0` joint-realism passes. The maximum number of assays hit by any setting was `1`, below the preregistered minimum of `2`.

All candidate hits were to `cr9114_h3_phillips2021_graphfla`; no setting simultaneously hit a second assay. This means the current protein-scale descriptor overlap does not yet generalize even within the selected Tier-1 joint-realism panel.

## Descriptor-level pattern

The coverage table shows that accessibility is the easiest descriptor block in this gate, while fitness-scale and autocorrelation are among the hardest. The highest dataset-level coverage is `0.652`, again for the CR9114-H3 panel, not for a broad multi-assay region.

## Interpretation

`RUN-138` fails the main Part 1 joint-realism criterion. This is stronger evidence than the Stage 4C overlap boundary because it requires one knob setting to hit several descriptors and several assays simultaneously. That did not occur.

The result weakens the current public generator as a ready empirical benchmark manifold. It does not refute the entire additive-core plus GNK idea, but it shows that this exact protein-scale grid and descriptor implementation cannot yet support the multi-peaked, navigable, correct-spectrum, correct-DFE claim.

## Decision

Record this as a failed joint-realism gate. Do not use the current Stage 4D parameter grid as a matched benchmark landscape set for downstream optimizer or PLM-agent comparisons.
