# ANA-054: RUN-139 Stage 4E held-out descriptor generalisation review

## Purpose

Collect and interpret `RUN-139` from `EXP-056`, the Stage 4E leave-one-landscape-out descriptor generalisation test.

This analysis asks whether knob regions selected on the training Tier-1 landscapes predict held-out landscape descriptors within empirical bands.

## Linked experiments/runs

- Experiments: `EXP-056`
- Runs: `RUN-139`

## Produced artifacts

- `source_tables/run-139_summary.json`
- `source_tables/run-139_preflight.json`
- `source_tables/run-139_descriptor_scope_audit.csv`
- `source_tables/run-139_heldout_predictions.csv`
- `tables/metrics.json`
- `tables/run-139_output_collection_audit.md`
- `tables/run-139_heldout_predictions_ranked.md`
- `tables/run-139_selected_parameter_counts.md`

## Collection policy

`RUN-139` completed on `lab-slurm` as Slurm job `1946` with inferred scheduler exit code `0`. The compact runtime output bundle was collected and DVC-tracked at:

`data/processed/stage-4e-heldout-descriptor-generalisation/RUN-139.dvc`

The original runtime output remains at `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-139/outputs`.

## Held-out decision

The held-out result is mixed to negative. The median held-out fraction of metrics inside the empirical band was `0.500`. `6/10` held-out datasets reached at least `0.5`, but only `2` exceeded a strict majority and only `2` reached `0.6`.

The selected parameter region also collapsed to a narrow choice: `9/10` held-out fits selected `random__sm6.000__ew1.000__seed703`. This does not look like a rich, landscape-specific calibrated manifold.

## Dataset-level pattern

The strongest held-out rows were `cr9114_h3_phillips2021_graphfla` and `cr6261_h9_phillips2021_graphfla`, each at `0.625` fraction in band. TrpB, CR9114-H1, and CR6261-H1 were only `0.375`. CR9114-FluB had a very large held-out mean-z distance despite `0.429` fraction in band, indicating a severe scale/descriptor mismatch in at least one block.

## Interpretation

`RUN-139` does not satisfy the planned held-out generalisation criterion. It provides some partial overlap on individual descriptors, but the majority-of-descriptors criterion is not met for most held-out landscapes, and the selected knob region is not robustly landscape-specific.

This weakens the current Part 1 claim that the implemented protein-scale public generator can predict unseen empirical landscape descriptor packages from a calibrated envelope.

## Decision

Do not treat the Stage 4 grid as held-out validated. The result should be cited as a generalisation boundary for the current implementation, not as evidence that the generator already covers empirical Tier-1 landscapes out of sample.
