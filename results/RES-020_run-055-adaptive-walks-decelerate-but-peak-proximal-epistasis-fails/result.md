# RES-020: RUN-055 adaptive walks decelerate but peak-proximal epistasis does not validate HYP-004

## Summary

RUN-055 supports two narrow HYP-004 expectations, adaptive-walk deceleration and declining beneficial fraction near higher-fitness backgrounds, but only 2/7 preregistered validation tests passed. Mean beneficial effects did not decline near the peak, sign epistasis was absent, negative epistasis decreased rather than increased, and synthetic-alignment near-peak backgrounds exceeded the stored peak fitness. The run therefore weakens HYP-004 as stated while preserving a narrower adaptive-walk result.

## Generated from

- Analyses: `ANA-020`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-004`
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-020_run-055-adaptive-walk-peak-proximal-epistasis-review/tables/run-055_key_metrics.md`
- `analyses/2026-05-15_ANA-020_run-055-adaptive-walk-peak-proximal-epistasis-review/tables/run-055_validation_tests.md`
- `analyses/2026-05-15_ANA-020_run-055-adaptive-walk-peak-proximal-epistasis-review/tables/run-055_background_summary_by_level.md`
- `analyses/2026-05-15_ANA-020_run-055-adaptive-walk-peak-proximal-epistasis-review/tables/run-055_greedy_walk_summary_by_level.md`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/summary.json`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/literature_baselines.json`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/background_neighborhoods.csv`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/epistasis_pairs.csv`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/walk_summaries.json`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/walk_trajectories.csv`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055.dvc`
- `experiments/2026-05-14_EXP-020_adaptive-walk-peak-proximal-epistasis-biological-baselines/config.yaml`
- `experiments/2026-05-14_EXP-020_adaptive-walk-peak-proximal-epistasis-biological-baselines/runs/RUN-055.yaml`

## Interpretation

- `RUN-055` completed the first live `EXP-020` test of `HYP-004` against
  literature-derived biological expectations for adaptive walks and
  peak-proximal epistasis.
- The adaptive-walk part is directionally consistent with real biological
  diminishing-returns adaptation: greedy walks decelerated, with median late
  gain divided by early gain `0.4206`, and beneficial-mutation availability
  declined with fraction of peak (Spearman `-0.2748`; low mean `0.06592`, high
  mean `0.05524`).
- This support is narrow. The overall preregistered validation package passed
  only `2/7` tests.
- The available beneficial-effect distribution did not behave like a bounded
  near-optimum DFE in the combined result. Mean beneficial effect was higher in
  high-fitness backgrounds than in low-fitness backgrounds (`69.12` versus
  `37.44`).
- The peak-proximal epistasis prediction failed. Mean absolute epistasis rose
  only weakly with fraction of peak (Spearman `0.151`), negative epistasis
  decreased strongly rather than increasing (Spearman `-0.6609`), and sampled
  sign epistasis was absent.
- The synthetic-alignment condition is a critical confound for interpreting
  peak proximity. Its near-peak bin averaged fraction of stored peak fitness
  `1.748`, so the stored peak is not an effective upper bound in that branch.
- The evidence therefore supports a narrower claim that the system can produce
  decelerating adaptive walks and declining beneficial availability, but it
  does not validate the stronger `HYP-004` claim about biologically plausible
  peak-proximal epistasis.

## Effect on hypothesis

- `HYP-004` is weakened as stated. The run does not justify citing the current
  synthetic landscape as reproducing the expected increase in constraining
  epistasis near a peak.
- The result does not refute the whole adaptive-walk direction because the
  walk deceleration and beneficial-fraction tests passed, especially in the
  unconditioned summaries.

## Limitations

- The combined validation is confounded by the synthetic-alignment branch,
  where the stored peak is not a reliable upper bound for the sampled
  near-peak backgrounds.
- The result uses sampled neighborhoods and sampled pairwise epistasis pairs,
  not an exhaustive landscape enumeration.
- The absence of sampled sign epistasis may reflect either the current
  landscape parameterization or the sampling design; this run alone does not
  separate those explanations.
- The biological comparison is qualitative and baseline-derived rather than a
  direct fit to a specific empirical adaptive-walk dataset.

## Downstream use

- Use `RES-020` as the analysis record for `EXP-020` / `RUN-055`.
- Do not cite `RUN-055` as support for the peak-proximal epistasis component
  of `HYP-004`.
- It is acceptable to cite the narrower adaptive-walk result: deceleration and
  declining beneficial availability are present in this run.
