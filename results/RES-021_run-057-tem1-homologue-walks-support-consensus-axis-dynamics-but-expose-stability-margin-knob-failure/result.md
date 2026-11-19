# RES-021: RUN-057 TEM-1 homologue walks support consensus-axis adaptive dynamics but expose stability-margin knob failure

## Summary

RUN-057 passes 7/9 validation checks: TEM-1 homologue starts farther from the MMSeqs consensus have lower fitness, higher beneficial availability, larger adaptive-walk gains, decelerating greedy walks, and slightly lower random-walk autocorrelation; epistasis strength also increases absolute epistasis and random-step variance. The two failures are both stability-margin knob checks, and identical-start diagnostics show that nominal stability_margin values from 4 to 12 did not change realized background fitness or stability_margin_used, so the result supports the HYP-004 consensus-axis walk behaviour while weakening HYP-005 for the current stability-margin parameterization.

## Generated from

- Analyses: ANA-021

## Relevant hypotheses

- Supports: HYP-004
- Weakens: HYP-005
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-021_run-057-tem1-homologue-walk-parameter-grid-diagnostic-review/tables/run-057_key_metrics.md`
- `analyses/2026-05-15_ANA-021_run-057-tem1-homologue-walk-parameter-grid-diagnostic-review/tables/run-057_validation_tests.md`
- `analyses/2026-05-15_ANA-021_run-057-tem1-homologue-walk-parameter-grid-diagnostic-review/tables/run-057_background_summary_by_start_stratum.md`
- `analyses/2026-05-15_ANA-021_run-057-tem1-homologue-walk-parameter-grid-diagnostic-review/tables/run-057_walk_summary_by_mode_and_stratum.md`
- `analyses/2026-05-15_ANA-021_run-057-tem1-homologue-walk-parameter-grid-diagnostic-review/tables/run-057_parameter_grid_summary.md`
- `analyses/2026-05-15_ANA-021_run-057-tem1-homologue-walk-parameter-grid-diagnostic-review/tables/run-057_stability_margin_propagation_check.md`
- `results/RES-021_run-057-tem1-homologue-walks-support-consensus-axis-dynamics-but-expose-stability-margin-knob-failure/metrics.json`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/summary.json`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/literature_baselines.json`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/alignment_profile_summary.json`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/homologue_starts.csv`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/background_neighborhoods.csv`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/walk_summaries.csv`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/walk_trajectories.csv`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/epistasis_pairs.csv`
- `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057/parameter_grid_summary.csv`
- `experiments/2026-05-15_EXP-021_tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-021_tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/runs/RUN-057.yaml`

## Interpretation

- `RUN-057` successfully replaced the confounded stored-peak axis from
  `RUN-055` with a more biologically defensible TEM-1 homologue
  consensus-distance axis. The selected homologue starts span from
  consensus-like backgrounds to much more distant homologues, using the same
  MMSeqs alignment family that conditions the landscape.
- The consensus-distance results are internally coherent and consistent with
  the stability/robustness expectation encoded in `EXP-021`. Farther
  homologues from the family consensus have lower starting fitness, fewer
  neutral neighbours, more deleterious neighbours, more beneficial
  opportunities, and larger greedy adaptive gains.
- The adaptive-walk component of `HYP-004` is therefore supported more strongly
  than it was in `RUN-055`: greedy walks decelerate, high-consensus backgrounds
  exhaust beneficial moves quickly, and far-consensus backgrounds retain
  substantially more adaptive opportunity.
- The random-walk result is weaker but directionally consistent. Lag-1
  autocorrelation declines with consensus distance, but all lag-1 correlations
  remain near `0.999`, so this run supports a smooth stability-dominated
  landscape rather than strong rugged random-walk decorrelation.
- The epistasis-strength parameter behaves directionally for mean absolute
  epistasis and random-step variance. This supports that part of `HYP-005`.
- The stability-margin parameter does not behave directionally. The mean
  neutral and deleterious fractions are essentially unchanged from nominal
  margin `4` to `12`, and matched-start diagnostics show that the exported
  realized `stability_margin_used` also did not change. This is best treated as
  a parameterization or propagation failure in the current setup, not as a
  biological result about TEM-1.
- Sampled sign epistasis remains absent, and negative epistasis is rare. The
  result should not be cited as evidence that the current TEM-1 configuration
  reproduces rich sign-epistatic ruggedness.

## Effect on hypothesis

- `HYP-004` is supported for the EXP-021 consensus-axis adaptive/random walk
  diagnostic. This run shows biologically plausible adaptive-walk opportunity
  and deceleration when starts are real TEM-1 homologues ranked by alignment
  consensus distance.
- `HYP-005` is weakened as a complete mechanistic-knob claim. The
  `epistasis_strength` knob is directional, but the `stability_margin` knob
  failed both preregistered DFE-direction tests and appears not to alter the
  realized stability margin for matched starts in this run.

## Limitations

- The consensus-distance score is alignment-derived and used as a proxy for
  closeness to a stability optimum. It is not a direct empirical measurement of
  TEM-1 fitness, stability, or distance to a true global peak.
- The result uses sampled single-mutant neighbourhoods and sampled epistasis
  pairs, not exhaustive enumeration of the 286-residue TEM-1 landscape.
- Random-walk autocorrelation effects are statistically directional but small
  in magnitude because the walks remain very smooth at lag 1.
- The `stability_margin` grid cannot be interpreted as a clean biological
  sweep until the realized-margin invariance is resolved.
- The run used `epistasis_covariation_bias = 0.0`, so it does not test whether
  covariation-targeted pair placement improves biological epistasis structure.

## Downstream use

- Use `RES-021` as the completed result record for `EXP-021` / `RUN-057`.
- Cite `RUN-057` as support for the revised `HYP-004` adaptive/random walk
  behaviour on TEM-1 homologue starts.
- Do not cite `RUN-057` as support that all `HYP-005` knobs are directional.
  It specifically weakens the current stability-margin parameterization while
  supporting the epistasis-strength component.
