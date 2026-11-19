# RES-051: RUN-136 protein-scale GNK reachability weakens broad directional knob claim

## Summary

RUN-136 weakens the broad protein-scale HYP-005 directional-knob claim: only 20/54 monotonicity checks passed, accessibility and reciprocal-sign epistasis passed 0/9 checks, and the public observable was not clearly separated from controls.

## Generated from

- Analyses: ANA-051

## Relevant hypotheses

- Supports: None
- Weakens: HYP-005
- Refutes: None
- Motivates: HYP-004, HYP-009

## Evidence

- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/source_tables
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/metrics.json
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_epi_weight_metric_means.csv
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_epi_weight_metric_means.md
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_failed_monotonicity_checks.csv
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_failed_monotonicity_checks.md
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_monotonicity_by_metric.csv
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_monotonicity_by_metric.md
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_monotonicity_by_stability_margin.csv
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_monotonicity_by_stability_margin.md
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_monotonicity_by_support_mode.csv
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_monotonicity_by_support_mode.md
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_output_collection_audit.csv
- analyses/2026-06-19_ANA-051_run-136-stage-4b-generator-reachability-review/tables/run-136_output_collection_audit.md
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-136/outputs

## Interpretation

`RUN-136` completed the Stage 4B reachability sweep for the public `additive_thermodynamic_gnk` model at protein-scale sampled-subcube resolution. The run evaluated `135` parameter settings and produced `540` generator descriptor rows plus `54` registered monotonicity checks.

The result is mixed to negative for the planned monotone-knob criterion. Only `20/54` checks passed. The failures were not minor or confined to a single nuisance descriptor: accessibility passed `0/9` monotonicity checks, reciprocal-sign epistasis passed `0/9`, sign epistasis passed `5/9`, Walsh high-order fraction passed `5/9`, random-walk autocorrelation passed `4/9`, and local maxima passed `6/9`.

The structured support mode did best, but still only passed `8/18` checks. The highest stability margin also passed `8/18`. This is not enough to claim a robust protein-scale directional control map.

The ablation controls matter. Aggregate summaries for the full public observable and additive-core-only observable were very similar across epistasis weights. Bare GNK texture behaved differently, but it is not the intended public composite. Therefore later empirical overlap cannot be credited to the additive thermodynamic core plus GNK mechanism unless it separates cleanly from these controls.

This result does not overturn `RES-049`. `RES-049` tested a deliberately small, enumerable `L = 8` hashed-GNK setting and passed there. `RUN-136` shows that the same conceptual knob map does not yet extend cleanly to the protein-scale sampled-subcube implementation.

## Effect on hypothesis

This weakens the broad protein-scale form of `HYP-005`. The implemented public model is executable and produces descriptor variation, but the registered directional-control criterion is not met.

The result still motivates `HYP-004` and `HYP-009` in a limited sense because it supplies the protein-scale generated descriptor grid and control outputs needed for overlap and benchmark-gating tests. It should not be cited as support for those hypotheses.

## Limitations

- The sweep is a sampled-subcube protein-scale diagnostic, not an exhaustive full-sequence-space landscape.
- Failure of monotonicity here is a failure of the current implementation/grid/descriptor setting, not a formal refutation of all additive-core plus GNK models.
- Control overlap means later empirical contact must be interpreted mechanistically, not only geometrically.
- The run does not test optimizers, PLM-guided agents, or empirical held-out prediction directly.

## Downstream use

- Do not cite `RUN-136` as a successful protein-scale reachability validation.
- Use it as the control-aware generated descriptor grid for Stage 4C-4E.
- Preserve the distinction between the positive small-L result in `RES-049` and the negative protein-scale monotonicity result here.
