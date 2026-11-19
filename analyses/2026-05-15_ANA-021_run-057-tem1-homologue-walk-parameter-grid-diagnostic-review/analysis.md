# ANA-021: RUN-057 TEM-1 homologue walk parameter-grid diagnostic review

## Purpose

Assess whether RUN-057 supports HYP-004 and HYP-005 by reviewing TEM-1 homologue adaptive/random walk behaviour, consensus-distance trends, and directional parameter-grid effects.

## Linked experiments/runs

- Experiments: EXP-021
- Runs: RUN-057

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `tables/run-057_key_metrics.md`
- `tables/run-057_validation_tests.md`
- `tables/run-057_background_summary_by_start_stratum.md`
- `tables/run-057_walk_summary_by_mode_and_stratum.md`
- `tables/run-057_parameter_grid_summary.md`
- `tables/run-057_stability_margin_propagation_check.md`

## Inputs

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
- `results/RES-020_run-055-adaptive-walks-decelerate-but-peak-proximal-epistasis-fails/result.md`

## Analysis performed

- Confirmed from run metadata that `RUN-057` completed successfully as the live
  `EXP-021` execution on `lab-slurm` with scheduler job `203` and exit code
  `0`.
- Reviewed the experiment question under `HYP-004` and `HYP-005`: whether
  TEM-1-conditioned landscapes show biologically plausible adaptive/random walk
  behaviour from real homologue starts, and whether stability-margin and
  epistasis-strength knobs move diagnostics directionally.
- Checked that the run used a real TEM-1 beta-lactamase reference sequence
  (`UniProt P62593`) and a 768-sequence MMSeqs-derived alignment profile.
- Reviewed the selected homologue starts ranked by alignment-consensus
  distance. The nine starts span `68.25` to `467.16` consensus-distance bits,
  `60` to `150` consensus Hamming substitutions, and `0.790` to `0.476`
  consensus identity.
- Summarized the run across `9` parameter sets, `162` sampled backgrounds,
  `486` walk summaries, `16,355` walk-trajectory rows, and `15,552` sampled
  epistasis pairs.
- Recomputed table-level summaries from the exported CSV/JSON files using the
  preregistered validation checks in `summary.json`.
- Added a targeted propagation check for the failed stability-margin knob tests
  by comparing identical starts across `stability_margin` values at
  `epistasis_strength = 0`.

## Main observations

- The run passed `7/9` validation checks, a clear improvement over the `2/7`
  pass pattern in `RUN-055` / `RES-020`.
- The consensus-distance axis behaved as expected for a stability-conditioned
  TEM-1 landscape. Background fitness declined strongly with consensus
  distance (Spearman `-0.9609`; low-distance mean `6983`, high-distance mean
  `2595`), while beneficial-mutation availability increased strongly
  (Spearman `0.9531`; low-distance mean `0.000917`, high-distance mean
  `0.03960`).
- Greedy adaptive-walk gain increased with consensus distance (Spearman
  `0.8609`; near-consensus mean gain `576.9`, far-consensus mean gain
  `3137`), consistent with more adaptive opportunity farther from the
  consensus-like high-fitness region.
- Greedy adaptive walks showed deceleration, with median late/early gain ratio
  `0.6849`, below the preregistered maximum `0.75`.
- Random walks showed the expected direction for autocorrelation, but the effect
  was small. Lag-1 autocorrelation declined with consensus distance at Spearman
  `-0.2015`, while mean lag-1 values remained high (`0.9989` near consensus and
  `0.9987` far from consensus).
- Epistasis strength behaved directionally for the readouts tested here:
  increasing `epistasis_strength` increased mean absolute pairwise epistasis
  from effectively zero to `0.004305` and increased random-step standard
  deviation from `224.14` to `227.23`.
- Sign epistasis remained absent in the sampled pairs, and negative epistasis
  was nearly absent except for one low-frequency stratum. This run therefore
  supports broad consensus-axis walk behaviour more strongly than it supports
  rich sign-epistatic ruggedness.
- The only failed validation checks were the two stability-margin knob checks:
  neutral fraction did not increase enough with nominal `stability_margin`
  (Spearman `0.1927`, below the `0.2` threshold), and deleterious fraction did
  not decline enough (Spearman `-0.1376`).
- The stability-margin failures are not just weak biological effects. For
  identical starts with `epistasis_strength = 0`, changing nominal
  `stability_margin` from `4` to `12` left both `background_fitness` and
  `stability_margin_used` unchanged in the exported diagnostics. That indicates
  the current run does not provide a valid positive test of the
  stability-margin knob; it exposes a parameterization or propagation failure
  that matters for `HYP-005`.

## Result records created

- `RES-021`

## Hypothesis updates

- `HYP-004` is supported for the revised EXP-021 diagnostic: TEM-1 homologue
  starts farther from the MMSeqs consensus show lower fitness, more beneficial
  opportunity, larger adaptive gains, decelerating greedy walks, and slightly
  lower random-walk autocorrelation.
- `HYP-005` is weakened as a complete knob-directionality claim. The
  `epistasis_strength` knob moved its target diagnostics, but the
  `stability_margin` knob did not move neutral/deleterious DFE metrics in this
  run and appears not to have altered the realized margins for identical
  starts.
