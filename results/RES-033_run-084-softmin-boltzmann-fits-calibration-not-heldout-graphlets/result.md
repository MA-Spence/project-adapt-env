# RES-033: RUN-084 soft-min Boltzmann ABC fits calibration but not held-out graphlets

## Summary

`RUN-084` is a valid PHOT soft-min/Boltzmann two-state ABC test under
`HYP-008`. It shows that finite-temperature state competition can improve
calibration graphlet distance while preserving the frozen base landscape, and
the posterior mean slightly improves the aggregate PHOT assay objective.
However, the calibration-selected best particle and posterior mean both fail
the central held-out graphlet test: neither beats fitted-current or smooth-core
controls on test graphlet distance, and neither repairs the empirical
topology/spectral/accessibility package. Locality diagnostics show concentrated
but very large artificial innovation energy rather than empirical-like local,
non-diffuse ruggedness. The result weakens the tested soft-min/Boltzmann
implementation of `HYP-008`.

## Generated from

- Analyses: `ANA-033`
- Experiments: `EXP-033`
- Runs: `RUN-084`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-008`
- Refutes: None

## Evidence

- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/metrics.json`
- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/tables/run-084_preflight_and_output_audit.md`
- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/tables/run-084_smc_posterior_summary.md`
- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/tables/run-084_graphlet_distance_comparison.md`
- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/tables/run-084_selected_graphlet_metrics.md`
- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/tables/run-084_assay_fit_comparison.md`
- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/tables/run-084_locality_tests.md`
- `analyses/2026-05-17_ANA-033_run-084-phot-softmin-boltzmann-abc-review/tables/run-084_hyp008_decision_table.md`
- `results/RES-033_run-084-softmin-boltzmann-fits-calibration-not-heldout-graphlets/metrics.json`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084.dvc`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084/summary.json`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084/model_distances.csv`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084/graphlet_summary.csv`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084/assay_fit_metrics.csv`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084/local_signal_graphlet_energy_summary.csv`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084/local_signal_feature_energy_summary.csv`
- `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084/locality_covariance_summary.csv`

## Interpretation

The run is procedurally valid. It used `PHOT_CHLRE_Chen_2023` only, restored
the `RUN-077` PHOT system-specific output, froze the primary `a(x)` landscape
to the fitted smooth-core PHOT base, and fitted only the auxiliary/operator ABC
parameters. It produced all expected posterior, assay, graphlet, WHT/GFT,
model-distance, and locality tables. The collected output is DVC-tracked.

The calibration result is a technical positive. The ABC-selected best particle
has calibration graphlet distance `0.6100`, better than fitted current
(`0.7318`), smooth core (`0.7454`), fitted-current posterior mean (`0.7839`),
and GNK (`0.8411`). The selected particle also preserves the base landscape in
the ABC guardrails: values-to-base correlation is `0.9449`, edge-gradient
correlation is `0.9505`, and total regularization penalty is `0`.

The held-out result is negative for the registered `HYP-008` support criterion.
The calibration-selected best particle has test graphlet distance `0.4434`,
worse than fitted current (`0.4016`) and smooth core (`0.4119`). The posterior
mean is also worse (`0.4543`). One posterior particle,
`abc_particle_rank_023`, narrowly beats fitted current on held-out test
distance (`0.3964` versus `0.4016`), but its calibration distance is poor
(`0.7126`), so it is a test-oracle capacity signal rather than valid held-out
selection evidence.

The PHOT assay guardrail does not rescue the interpretation. The posterior mean
has the best aggregate assay objective (`0.4570` versus `0.4631` for fitted
current posterior mean, `0.4814` for smooth core, and `0.4856` for fitted
current), but holdout Spearman and single/double-mutant validation correlations
do not improve as a package, and epistasis-prediction Spearman remains negative
(`-0.0881`). Assay improvement was explicitly a guardrail/reporting output, not
the support criterion.

The graphlet package remains mismatched. Empirical held-out graphlets have WHT
high-order fraction `0.3788`, mean local maxima `5.17`, accessibility `0.1975`,
and random-walk lag-1 autocorrelation `0.1905`. The calibration-selected
soft-min particle has WHT high-order fraction `0.5268`, local maxima `4.25`,
accessibility `0.0037`, and autocorrelation `0.0405`. That is not the empirical
topology/spectral/evolutionary package.

Locality diagnostics are also negative. The calibration-selected particle has
test plaquette innovation energy `179,908`, compared with only `15.4` for the
empirical readout residual, and shared-feature covariance is weak
(`0.0934` plaquette-energy similarity Spearman). Concentrated state-transition
energy is therefore not interpretable as empirical-like local non-diffuse
ruggedness.

## Effect on hypothesis

`RES-033` weakens the tested soft-min/Boltzmann implementation of `HYP-008`.
The result does not refute every possible nonlinear trait-operator model, but
it shows that this finite-temperature two-state ABC parameterization does not
pass the required calibration-selected held-out graphlet and locality controls.

## Limitations

- PHOT is the only empirical system tested.
- The empirical target is the processed `latent_phi` graphlet layer, not direct
  exhaustive measurement of every dense graphlet.
- The base landscape is frozen to the `RUN-077` PHOT fit, so this result tests
  the auxiliary/operator layer rather than a joint refit.
- One posterior particle has a weak test-oracle distance advantage; this is not
  valid selection evidence but it means the operator family has some capacity.
- The auxiliary construction uses one deterministic coefficient seed and a
  compact ABC parameterization.

## Downstream use

- Do not cite `RUN-084` as support for `HYP-008`.
- It is valid to cite `RUN-084` as evidence that soft-min/Boltzmann state
  competition can fit calibration graphlets while failing calibration-selected
  held-out graphlet and locality tests.
- Keep `RES-031`, `RES-032`, and `RES-033` together as the current PHOT
  coupled-operator record: no-module coupling can move graphlet distance,
  hard-min can improve assay objective, and soft-min/Boltzmann ABC can fit
  calibration graphlets, but none yet validates local, non-diffuse ruggedness
  under the registered held-out controls.
