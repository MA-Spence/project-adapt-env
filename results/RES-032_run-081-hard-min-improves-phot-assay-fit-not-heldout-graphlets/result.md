# RES-032: RUN-081 hard-min improves PHOT assay fit but not held-out graphlet controls

## Summary

`RUN-081` is a valid deterministic PHOT hard-min weakest-link grid. It shows
that some hard-min candidates improve the full-assay PHOT objective relative to
frozen smooth-core and fitted-current controls. However, no hard-min candidate
beats the `RUN-077` fitted-current model on held-out graphlet distance, the
calibration-selected hard-min candidate does not generalize, the PHOT spectral
envelope remains mismatched, and local-signal diagnostics do not validate
empirical-like non-diffuse ruggedness. The result weakens the current hard-min
implementation of `HYP-008`.

## Generated from

- Analyses: `ANA-032`
- Experiments: `EXP-032`
- Runs: `RUN-081`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-008`
- Refutes: None

## Evidence

- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/metrics.json`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_preflight_and_output_audit.md`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_assay_fit_comparison.md`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_graphlet_distance_comparison.md`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_selected_graphlet_metrics.md`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_graphlet_energy_locality_test.md`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_feature_plaquette_locality_test.md`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_locality_covariance_test.md`
- `analyses/2026-05-16_ANA-032_run-081-phot-hard-min-weakest-link-grid-review/tables/run-081_hyp008_decision_table.md`
- `results/RES-032_run-081-hard-min-improves-phot-assay-fit-not-heldout-graphlets/metrics.json`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081.dvc`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081/summary.json`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081/assay_fit_metrics.csv`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081/model_distances.csv`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081/graphlet_metrics.csv`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081/local_signal_graphlet_energy_summary.csv`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081/local_signal_feature_energy_summary.csv`
- `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081/locality_covariance_summary.csv`

## Interpretation

The assay-level result is a narrow positive technical signal. The best hard-min
candidate, `hard_min__smooth_core_A__grid0100`, improves the assay objective to
`0.4223`, compared with smooth core `0.4814`, fitted current `0.4856`, and
posterior mean `0.4631`. The improvement does not amount to broad empirical
recovery: holdout Spearman is effectively unchanged, single-mutant Spearman is
slightly worse than controls, and epistasis-prediction Spearman remains weak
(`0.0106`).

The graphlet result is negative for the central `HYP-008` claim. The best
hard-min test distance is `0.4124`, worse than the fitted-current control
(`0.4016`) and slightly worse than smooth core (`0.4119`). The
calibration/assay-selected hard-min candidate has test distance `0.4631`, so
the calibration improvement does not generalize. The empirical readout residual
remains closest on test graphlets (`0.2708`), but it is a diagnostic
subtraction from empirical `latent_phi`, not a procedural generator.

The hard-min operator also fails to repair the PHOT spectral and topology
package. Empirical test graphlets have WHT high-order fraction `0.3788`, mean
local maxima `5.17`, accessibility `0.1975`, and random-walk lag-1
autocorrelation `0.1905`. The best assay hard-min has WHT high-order fraction
`0.5077`, local maxima `3.42`, accessibility `0.0245`, and lag-1
autocorrelation `0.0464`. This is not the empirical graphlet package.

Locality is not validated. Hard-min switch boundaries concentrate model
innovation by construction, but the signal is not empirical-like. On test
graphlets, best-assay hard-min plaquette energy is `32,593`, compared with
`15.4` for the empirical readout residual, and shared-feature covariance with
innovation/switch similarity is weak. Concentrated switch energy therefore
cannot be interpreted as evidence for biological local non-diffuse ruggedness
in this implementation.

## Effect on hypothesis

`RES-032` weakens `HYP-008` as currently operationalized by the hard-min
weakest-link model. The run does not support the claim that this nonlinear
trait operator supplies the missing local, non-diffuse PHOT ruggedness.

The result does not refute every possible version of `HYP-008`. It is a
PHOT-only deterministic frozen-base grid with one compact auxiliary-seed
design. The honest conclusion is narrower: hard-min coupling can improve an
assay objective, but this tested implementation does not pass the held-out
graphlet and locality controls required for `HYP-008` support.

## Limitations

- PHOT is the only empirical system tested.
- The run freezes the `RUN-077` PHOT base and does not perform a joint ABC
  refit.
- The empirical graphlet target is the processed `latent_phi` layer.
- The auxiliary grid uses one seed and a limited hard-min parameter range.
- Test-oracle rankings are capacity diagnostics only and are not valid
  held-out selection evidence.

## Downstream use

- Do not cite hard-min weakest-link coupling as support for `HYP-008`.
- It is valid to cite `RUN-081` as evidence that the operator can improve a
  PHOT assay objective while failing the required held-out graphlet and
  locality controls.
- Retain `RES-030`, `RES-031`, and `RES-032` together as the current PHOT
  graphlet record: system-specific current refitting partially helps, the
  first coupled/local-module grid fails locality, and hard-min coupling does
  not rescue the held-out graphlet package.
