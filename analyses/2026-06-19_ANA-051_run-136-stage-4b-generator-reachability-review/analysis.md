# ANA-051: RUN-136 Stage 4B protein-scale GNK reachability review

## Purpose

Collect and interpret `RUN-136` from `EXP-053`, the Stage 4B reachability sweep for the public `additive_thermodynamic_gnk` model.

This analysis asks whether the small-L `RES-049` hashed-GNK directional knob result extends to the protein-scale sampled-subcube setting with structured, shuffled, and random routing controls.

## Linked experiments/runs

- Experiments: `EXP-053`
- Runs: `RUN-136`

## Produced artifacts

- `source_tables/run-136_summary.json`
- `source_tables/run-136_preflight.json`
- `source_tables/run-136_generator_parameter_grid.csv`
- `source_tables/run-136_generator_descriptor_metrics.csv`
- `source_tables/run-136_generator_knob_monotonicity.csv`
- `source_tables/run-136_core_texture_ablation_controls.csv`
- `tables/metrics.json`
- `tables/run-136_output_collection_audit.md`
- `tables/run-136_monotonicity_by_metric.md`
- `tables/run-136_monotonicity_by_support_mode.md`
- `tables/run-136_monotonicity_by_stability_margin.md`
- `tables/run-136_epi_weight_metric_means.md`
- `tables/run-136_failed_monotonicity_checks.md`

## Collection policy

`RUN-136` completed on `lab-slurm` as Slurm job `1943` with inferred scheduler exit code `0`. The compact runtime output bundle was collected and DVC-tracked at:

`data/processed/stage-4b-additive-core-gnk-generator-reachability/RUN-136.dvc`

The original runtime output remains at `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-136/outputs`.

## Output audit

The run evaluated `135` parameter sets and wrote `540` generator metric rows plus `54` monotonicity checks. The public model system recorded in the output is `additive_thermodynamic_gnk`.

## Directional knob result

The protein-scale reachability result is mixed to negative. Only `20/54` monotonicity checks passed (`0.370`). The best support mode was `structured`, with `8/18` checks passing; the best stability margin was `14`, also with `8/18` checks passing.

The failures are not confined to one route. Accessibility did not pass any monotonicity checks (`0/9` for `accessible_path_fraction_dp`), and reciprocal-sign epistasis also passed `0/9`. Sign epistasis and Walsh high-order fraction passed only `5/9` checks each. Local maxima passed `6/9`, so ruggedness count is more directional than accessibility, but still not a broad pass.

## Ablation/control interpretation

The full public observable and additive-core-only observable have very similar aggregate summaries across `epi_weight`, while the effective-stability-margin diagnostic remains mostly smooth and accessible until the highest epistasis weight. The bare GNK texture behaves differently but is not the public model. Thus any empirical overlap seen later must be checked for specificity to the public composite rather than credited automatically to the intended additive-core plus GNK mechanism.

## Interpretation

`RUN-136` weakens the broad protein-scale form of the `HYP-005` directional-knob claim. It does not overturn `RES-049`, because `RES-049` was a deliberately small, enumerable `L=8` validation. It does show that the current protein-scale sampled-subcube implementation does not yet give stable monotone control of the full ruggedness/accessibility/spectrum descriptor vector.

## Decision

Do not cite this run as a successful protein-scale reachability validation. Use it as a boundary: the public composite is executable and produces descriptor tables, but the registered monotonicity criterion is not met in this setting.
