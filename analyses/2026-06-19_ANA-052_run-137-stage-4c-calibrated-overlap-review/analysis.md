# ANA-052: RUN-137 Stage 4C calibrated overlap review

## Purpose

Collect and interpret `RUN-137` from `EXP-054`, the Stage 4C empirical coverage and calibrated-overlap test.

The analysis asks whether generated descriptor distributions overlap the empirical bands from Stage 4A, and whether that overlap is specific to the public additive-core plus GNK composite rather than reproduced equally by controls.

## Linked experiments/runs

- Experiments: `EXP-054`
- Runs: `RUN-137`

## Produced artifacts

- `source_tables/run-137_summary.json`
- `source_tables/run-137_preflight.json`
- `source_tables/run-137_empirical_bands.csv`
- `source_tables/run-137_descriptor_scope_audit.csv`
- `source_tables/run-137_coverage_metrics.csv`
- `source_tables/run-137_posterior_knob_regions.csv`
- `source_tables/run-137_control_overlap_failures.csv`
- `tables/metrics.json`
- `tables/run-137_output_collection_audit.md`
- `tables/run-137_best_overlap_by_dataset.md`
- `tables/run-137_coverage_by_model_family.md`
- `tables/run-137_posterior_by_model_family.md`
- `tables/run-137_control_overlap_failure_summary.md`

## Collection policy

`RUN-137` completed on `lab-slurm` as Slurm job `1944` with inferred scheduler exit code `0`. The compact runtime output bundle was collected and DVC-tracked at:

`data/processed/stage-4c-empirical-coverage-calibrated-overlap/RUN-137.dvc`

The original runtime output remains at `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-137/outputs`.

## Output audit

The run wrote `44820` coverage rows and `7560` posterior/knob-region rows over `135` parameter sets.

## Coverage and overlap

The best per-dataset overlap table gives a superficially positive coverage result for the Tier-1 combinatorial datasets: `10` datasets have a best candidate with at least `0.6` of considered metrics inside the empirical band, and `10` have at least `0.5`. The median best fraction is `0.688`.

However, the overlap is not specific to the public composite. The all-row within-band rate is `0.448` for `full_public_observable` and `0.452` for `bare_gnk_texture`, with additive-core-only controls in the same range. The best overlaps for the Tier-1 datasets are mostly bare-GNK texture rows. The sparse ProteinGym order-1 rows are not rescued: `4` sparse ProteinGym datasets have best overlap fraction `0.0` under the recorded one-metric comparison.

## Interpretation

`RUN-137` shows descriptor-space contact with the empirical envelope, but it does not validate the intended public mechanism. A coverage claim that is mostly matched by bare GNK or additive-core-only controls is not strong evidence that the additive thermodynamic core plus GNK texture jointly explains the empirical descriptor envelope.

This result is therefore a boundary rather than a success: the implementation can find overlap regions for several Tier-1 descriptor bands, but those regions are not yet discriminating enough to support the model-specific Part 1 claim.

## Decision

Do not use this run as evidence that the public generator covers empirical landscapes in a mechanistically specific way. Use the posterior/coverage tables to identify which descriptor blocks and controls must be tightened before any model-specific overlap claim.
