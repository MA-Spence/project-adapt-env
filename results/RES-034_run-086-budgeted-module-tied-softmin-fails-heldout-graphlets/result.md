# RES-034: RUN-086 budgeted module-tied soft-min fails held-out graphlet and locality support

## Summary

`RUN-086` tests a budgeted module-tied soft-min state-difference field for PHOT.
It slightly improves calibration graphlet distance but does not beat
fitted-current held-out graphlets, does not validate empirical module locality,
and provides no assay-guardrail rescue. The budgeted operator mostly collapses
back onto the frozen base landscape rather than producing empirical-like local
non-diffuse ruggedness. The result weakens this implementation of `HYP-008`.

## Generated from

- Analyses: `ANA-034`
- Experiments: `EXP-034`
- Runs: `RUN-086`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-008`
- Refutes: None

## Evidence

- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/metrics.json`
- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/tables/run-086_preflight_and_output_audit.md`
- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/tables/run-086_smc_posterior_summary.md`
- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/tables/run-086_graphlet_distance_comparison.md`
- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/tables/run-086_selected_graphlet_metrics.md`
- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/tables/run-086_assay_fit_comparison.md`
- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/tables/run-086_locality_tests.md`
- `analyses/2026-05-17_ANA-034_run-086-phot-budgeted-module-tied-softmin-abc-review/tables/run-086_hyp008_decision_table.md`
- `results/RES-034_run-086-budgeted-module-tied-softmin-fails-heldout-graphlets/metrics.json`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086.dvc`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086/summary.json`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086/model_distances.csv`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086/graphlet_summary.csv`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086/assay_fit_rankings.csv`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086/local_signal_graphlet_energy_summary.csv`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086/local_signal_feature_energy_summary.csv`
- `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086/locality_covariance_summary.csv`

## Interpretation

The run is procedurally valid. It used `PHOT_CHLRE_Chen_2023` only, restored
the `RUN-077` PHOT system-specific current-model output, froze the primary
`a(x)` landscape to the fitted smooth-core base, selected `12` module
coordinates from calibration empirical readout-residual plaquette energy, and
fit the budgeted state-difference operator
`g_tau(a,d)=a - tau log(1 + exp(-d/tau))`. The collected output contains the
planned posterior, assay, graphlet, WHT/GFT, model-distance, module, and
locality tables, and the output directory is DVC-tracked.

The calibration signal is small. The ABC-selected best particle has
calibration graphlet distance `0.7240`, slightly better than fitted current
(`0.7318`) and smooth core (`0.7454`) but far worse than the empirical
readout-residual diagnostic (`0.5428`). That is only a weak technical
calibration pass.

The registered held-out support criterion fails. The selected best particle has
held-out test graphlet distance `0.4060`, worse than fitted current (`0.4016`).
It is slightly better than smooth core (`0.4119`), but fitted current remains
the relevant current-model control. No posterior particle beats fitted current
on held-out graphlets; the best posterior particle reaches `0.4033`.

The module-locality claim is not supported. The selected empirical-module
particle is worse than the shuffled-module control on held-out graphlets
(`0.4060` versus `0.4029`), and the random-module and no-module variants are
nearby (`0.4081` and `0.4144`). The selected particle also has essentially zero
innovation relative to the base: values-to-base correlation is `1.0`, edge
gradient correlation is `1.0`, innovation variance ratio to base is `1.4e-25`,
and edge/plaquette innovation ratios to the empirical residual are `6.0e-22`
and `1.2e-21`.

The PHOT graphlet package remains mismatched. Empirical held-out graphlets have
mean local maxima `5.17`, accessibility `0.1975`, random-walk lag-1
autocorrelation `0.1905`, and WHT high-order fraction `0.3788`. The selected
budgeted particle has local maxima `3.50`, accessibility `0.0280`,
autocorrelation `0.0851`, and WHT high-order fraction `0.4984`.

The assay guardrail is also negative. The best assay objective remains the
fitted-current posterior mean (`0.4631`). The best budgeted candidate is worse
(`0.5350`), and the calibration-selected particle is worse again (`0.5560`).
Epistasis-prediction Spearman remains negative for the selected particle
(`-0.1745`).

## Effect on hypothesis

`RES-034` weakens the tested budgeted module-tied soft-min implementation of
`HYP-008`. The result does not refute every possible nonlinear local
trait-operator model, but it shows that this constrained state-difference
implementation does not pass the required held-out graphlet, locality, module
control, or assay guardrail checks.

## Limitations

- PHOT is the only empirical system tested.
- The empirical target is the processed `latent_phi` graphlet layer, not direct
  exhaustive measurement of every graphlet.
- The base landscape is frozen to the `RUN-077` PHOT fit, so the result tests
  the budgeted operator layer rather than a joint model refit.
- The budget constraints may be too restrictive, but relaxing them would need
  to avoid returning to the `RUN-084` artificial-innovation failure mode.
- The result is implementation-specific; it weakens the current budgeted
  module-tied formulation rather than all possible local nonlinear operators.

## Downstream use

- Do not cite `RUN-086` as support for `HYP-008`.
- It is valid to cite `RUN-086` as evidence that explicit innovation and
  state-boundary budgets can prevent artificial over-innovation, but in this
  implementation they do so by collapsing the fitted operator onto the frozen
  base.
- Keep `RES-031`, `RES-032`, `RES-033`, and `RES-034` together as the current
  PHOT coupled-operator record: coupled operators can move calibration or assay
  objectives, but none has yet validated local, non-diffuse empirical
  ruggedness under held-out controls.
