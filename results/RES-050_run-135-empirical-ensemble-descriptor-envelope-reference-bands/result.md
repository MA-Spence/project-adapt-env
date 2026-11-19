# RES-050: RUN-135 empirical ensemble descriptor envelope establishes operator-free reference bands

## Summary

RUN-135 establishes the Stage 4A empirical descriptor reference: 17 datasets were audited, 10 Tier-1 operator-free landscapes can support higher-order/topographic bands, and 7 sparse Tier-3 assays are restricted to order-1 summaries.

## Generated from

- Analyses: ANA-050

## Relevant hypotheses

- Supports: None
- Weakens: None
- Refutes: None
- Motivates: HYP-004, HYP-009

## Evidence

- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/source_tables
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/metrics.json
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_empirical_key_bands.csv
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_empirical_key_bands.md
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_graphlet_dataset_medians.csv
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_graphlet_dataset_medians.md
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_output_collection_audit.csv
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_output_collection_audit.md
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_scope_audit_summary.csv
- analyses/2026-06-19_ANA-050_run-135-stage-4a-empirical-envelope-review/tables/run-135_scope_audit_summary.md
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-135/outputs

## Interpretation

`RUN-135` completed the Stage 4A empirical-envelope step for `EXP-052`. The run collected an empirical panel spanning `17` datasets and produced `340` descriptor rows, `960` graphlet rows, and `17` scope-audit rows.

The most important outcome is the scope split. The record identifies `10` Tier-1 complete or near-complete combinatorial landscapes as valid for operator-free higher-order and topographic descriptors. The sparse Tier-3 ProteinGym-style datasets are explicitly restricted to order-1 summaries. That is the right guardrail after `RES-048`: higher-order epistasis, accessibility, graphlet, and peak claims must not be inferred from sparse measurement-operator-limited assays.

The Tier-1 empirical graphlet envelope is broad enough to serve as a target for the generator stages. Across the sampled graphlets, the median local-maxima count is `2`, the median dynamic-programming accessible-path fraction is `0`, the median Walsh high-order fraction is `0.01639`, the median sign-epistasis fraction is `0.333`, and the median random-walk lag-1 autocorrelation is `0.308`. The high quantiles still include substantially rugged and accessible regimes: the 95th percentile local-maxima count is `16`, and the 95th percentile accessible-path fraction is about `0.583`.

This is a positive result for the data-foundation step only. It establishes the empirical descriptor bands and their permitted interpretation levels. It does not show that the public generator reaches those bands.

## Effect on hypothesis

This result motivates `HYP-004` and `HYP-009` by giving later generator and benchmark experiments a scoped empirical target. It does not directly support either hypothesis, because no generated landscape is compared or validated here.

It has no direct effect on `HYP-005`, because the run does not test generator knob monotonicity.

## Limitations

- The sparse Tier-3 datasets remain order-1 only; they cannot support epistasis, accessibility, peak-count, or graphlet conclusions without a certified observation operator.
- The empirical bands are descriptor summaries, not a mechanistic model fit.
- The graphlet-level bands depend on the sampled subgraphs and descriptor implementation used in Stage 4A.
- The result is not evidence of public-generator realism, held-out generalisation, optimizer discrimination, or biological routing.

## Downstream use

- Use the `RUN-135` empirical bands as the fixed Stage 4A reference for `RUN-136` through `RUN-139`.
- Restrict higher-order/topographic empirical comparisons to the Tier-1 operator-free landscapes listed in the `RUN-135` summary and scope audit.
- Treat sparse ProteinGym-style rows as order-1 calibration or distributional context only.
