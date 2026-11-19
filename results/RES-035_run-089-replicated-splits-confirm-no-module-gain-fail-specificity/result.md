# RES-035: RUN-089 replicated splits confirm no-module graphlet gain but fail specificity controls

## Summary

`RUN-089` shows that the `RUN-079` no-module PHOT graphlet-distance
improvement is stable across replicated held-out graphlet splits against
fitted-current and smooth-core controls. It fails the mechanism-specificity
controls: the uncoupled-module control is better overall, and
random/shuffled/empirical module families also improve against current. This
weakens the current local non-diffuse `HYP-008` mechanism while preserving a
narrow capacity result.

## Generated from

- Analyses: ANA-035
- Experiments: EXP-035
- Runs: RUN-089

## Relevant hypotheses

- Supports: None
- Weakens: HYP-008
- Refutes: None

## Evidence

- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/metrics.json
- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/tables/run-089_component_summary.md
- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/tables/run-089_family_summary.md
- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/tables/run-089_hyp008_decision_table.md
- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/tables/run-089_original_split_reference.md
- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/tables/run-089_pairwise_control_deltas.md
- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/tables/run-089_preflight_and_output_audit.md
- analyses/2026-05-17_ANA-035_run-089-phot-replicated-split-graphlet-diagnostic-review/tables/run-089_top_selection_counts.md
- results/RES-035_run-089-replicated-splits-confirm-no-module-gain-fail-specificity/metrics.json
- data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089.dvc
- data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089/summary.json
- data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089/split_family_summary.csv
- data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089/split_family_pairwise_deltas.csv
- data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089/component_summary.csv

## Interpretation

The narrow held-out graphlet-distance claim is real under this diagnostic.
Across 512 random graphlet calibration/test splits, the no-module family has
median held-out graphlet distance `0.4839`, compared with `0.5232` for fitted
current and `0.5358` for smooth core. It beats fitted current in `93.8%` of
splits and smooth core in `97.7%` of splits. The median deltas are `-0.0468`
and `-0.0577`.

That evidence is not enough for `HYP-008`. The mechanism-specificity controls
fail. The uncoupled-module control has lower median held-out distance
(`0.4469`), is best in `51.2%` of splits, and beats no-module in `63.5%` of
random splits. Random-module (`0.4793` median), empirical-module (`0.4943`),
and shuffled-module (`0.4999`) families all beat fitted current frequently, so
the result is not specific to empirical local modules or to the no-module
coupled operator.

The original `RUN-079` split is not representative of the typical effect size.
On the original split, no-module held-out distance is `0.3137`, compared with
`0.4016` for fitted current and `0.4119` for smooth core. The replicated
median no-module distance is `0.4839`, indicating that the original split
substantially overstated the advantage.

The secondary components remain biologically and topologically mismatched. The
no-module family improves some distances relative to fitted current, but it
still underestimates accessibility (`0.0525` model versus `0.2181` empirical),
local maxima (`4.08` versus `5.17`), and lag-1 autocorrelation (`0.0607`
versus `0.2031`), while overestimating WHT high-order fraction (`0.5131`
versus `0.3583`).

## Effect on hypothesis

`RES-035` weakens the current interpretation of `HYP-008`. It confirms that
the `RUN-079` no-module graphlet-distance improvement against current/smooth
controls is stable across split replicates, but the result should not be cited
as support for local non-diffuse empirical module ruggedness. The strongest
control is a single uncoupled-module landscape, and non-empirical module
families share much of the improvement.

The supported claim is narrower: the `RUN-079` candidate family contains
calibration-selectable variants that can lower PHOT held-out graphlet distance
relative to fitted current and smooth core. The unsupported claim is that this
improvement identifies the `HYP-008` mechanism.

## Limitations

- PHOT is the only empirical system tested.
- The graphlet target is the processed `latent_phi` layer.
- The experiment reuses `RUN-079` candidate metrics and does not refit the
  base model or test a new operator.
- The diagnostic tests graphlet-split robustness, not assay generalization.
- The result weakens the current `RUN-079` mechanism interpretation, not every
  possible nonlinear local trait operator.

## Downstream use

- Cite `RES-035` as evidence that the no-module `RUN-079` graphlet-distance
  gain is stable against fitted-current and smooth-core controls under
  replicated held-out graphlet splits.
- Do not cite `RUN-079`/`RES-031` alone as support for empirical local-module
  `HYP-008` mechanism; pair it with `RES-035`.
- Treat future PHOT coupled-operator claims as needing an uncoupled-landscape
  control, random/shuffled module controls, and secondary component checks.
