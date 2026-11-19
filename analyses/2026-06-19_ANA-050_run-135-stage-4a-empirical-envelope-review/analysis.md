# ANA-050: RUN-135 Stage 4A empirical ensemble descriptor envelope review

## Purpose

Collect and interpret the `RUN-135` outputs from `EXP-052`, the Stage 4A empirical reference-distribution step in `scratch/science_plan_part1_ensemble_landscape_validation.md`.

The analysis asks whether the assembled empirical panel is usable as a fixed descriptor envelope while preserving the scope guardrail from the plan: sparse ProteinGym-style assays can support order-1 summaries, but not operator-free higher-order or topographic claims.

## Linked experiments/runs

- Experiments: `EXP-052`
- Runs: `RUN-135`

## Produced artifacts

- `source_tables/run-135_summary.json`
- `source_tables/run-135_preflight.json`
- `source_tables/run-135_descriptor_scope_audit.csv`
- `source_tables/run-135_descriptor_metrics.csv`
- `source_tables/run-135_empirical_bands.csv`
- `source_tables/run-135_graphlet_metrics.csv`
- `tables/metrics.json`
- `tables/run-135_output_collection_audit.md`
- `tables/run-135_scope_audit_summary.md`
- `tables/run-135_empirical_key_bands.md`
- `tables/run-135_graphlet_dataset_medians.md`

## Collection policy

`RUN-135` completed on `lab-slurm` as Slurm job `1941` with inferred scheduler exit code `0`. The compact runtime output bundle was collected and DVC-tracked at:

`data/processed/stage-4a-empirical-ensemble-descriptor-envelope/RUN-135.dvc`

The original runtime output remains at `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-135/outputs`.

## Output audit

The run evaluated `17` datasets and wrote `340` order-1 descriptor rows, `960` graphlet rows, and `17` scope-audit rows. The higher-order/operator-free subset contains `10` Tier-1 complete or near-complete combinatorial landscapes. The sparse Tier-3 subset contains `7` datasets that were explicitly blocked from higher-order and topographic interpretation.

This is the correct scope split for the paper plan. The run did not compute epistasis, graphlet, accessibility, or peak descriptors from sparse ProteinGym records. Those records contribute only the order-1 block unless a separate measurement operator is certified.

## Empirical envelope

The Tier-1 graphlet envelope is broad. Across the sampled graphlets, the empirical median local-maxima count is `2`, the median dynamic-programming accessible-path fraction is `0`, and the median Walsh high-order fraction is `0.01639`. The 95th percentile local-maxima count reaches `16`, and the 95th percentile accessible-path fraction reaches about `0.583`.

The median graphlet sign-epistasis fraction is `0.333`, the median magnitude-epistasis-to-single-mutant-SD ratio is `0.600`, and the median random-walk lag-1 autocorrelation is `0.308`. These summaries establish an empirical descriptor envelope with both rugged and relatively navigable regimes.

## Interpretation

`RUN-135` succeeds as an empirical reference-distribution run. It provides the operator-free descriptor bands needed by later stages and prevents the PHOT/ProteinGym measurement-operator problem from leaking into higher-order claims.

The result should not be cited as evidence that the generator covers the envelope. It is the empirical target surface against which `RUN-136` through `RUN-139` are judged.

## Decision

Use the `RUN-135` bands and scope audit as the empirical descriptor reference for Stage 4B-4E. Higher-order and topographic comparisons should be restricted to the `10` Tier-1 landscapes listed in the source `summary.json`; sparse Tier-3 assays should remain order-1 only.
