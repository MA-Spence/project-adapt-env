# ANA-020: RUN-055 adaptive-walk peak-proximal epistasis biological-baseline review

## Purpose

Assess whether RUN-055 supports HYP-004 by comparing adaptive-walk dynamics, beneficial-mutation availability, and peak-proximal epistasis trends against the preregistered literature-derived biological baselines.

## Linked experiments/runs

- Experiments: `EXP-020`
- Runs: `RUN-055`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: python3

## Inputs

- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/summary.json`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/literature_baselines.json`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/background_neighborhoods.csv`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/epistasis_pairs.csv`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/walk_summaries.json`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055/walk_trajectories.csv`
- `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055.dvc`
- `experiments/2026-05-14_EXP-020_adaptive-walk-peak-proximal-epistasis-biological-baselines/config.yaml`
- `experiments/2026-05-14_EXP-020_adaptive-walk-peak-proximal-epistasis-biological-baselines/runs/RUN-055.yaml`

## Analysis performed

- Confirmed from run metadata that `RUN-055` completed successfully on
  `lab-slurm` as scheduler job `201` with exit code `0`.
- Tracked the completed output directory with DVC at
  `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055.dvc`.
- Reviewed the `EXP-020` question under `HYP-004`: whether adaptive walks,
  beneficial-mutation availability, and peak-proximal epistasis behave like
  real biological/evolutionary systems rather than serving as optimizer
  benchmarks.
- Checked the literature-derived baselines encoded in the run output:
  diminishing-returns adaptation, decelerating adaptive trajectories,
  near-optimum depletion of beneficial effects, and rugged/sign-epistatic
  constraints on accessible paths.
- Summarized the run outputs across `120` sampled backgrounds, `240` walk
  summaries, `3029` walk-trajectory rows, and `30720` sampled epistasis pairs.
- Recomputed the preregistered validation table from `summary.json` and
  compared condition-specific background and walk summaries.

## Produced artifacts

- `tables/run-055_key_metrics.md`
- `tables/run-055_validation_tests.md`
- `tables/run-055_background_summary_by_level.md`
- `tables/run-055_greedy_walk_summary_by_level.md`

## Main observations

- The run passed `2/7` preregistered validation tests.
- The adaptive-walk deceleration expectation passed: the median greedy
  deceleration ratio was `0.4206`, below the preregistered maximum of `0.75`.
- Beneficial-mutation availability declined with higher fraction of peak:
  Spearman was `-0.2748`, with low-background mean beneficial fraction
  `0.06592` and high-background mean `0.05524`.
- The total-gain test failed as implemented. Although the rank correlation
  with start fraction of peak was negative (`-0.2542`), the high-start mean
  total gain was larger than the low-start mean (`1416` versus `1200`), driven
  by the synthetic-alignment condition.
- Mean beneficial effect did not decline near the peak. The high-background
  mean was larger than the low-background mean (`69.12` versus `37.44`), and
  the Spearman trend was weakly positive (`0.058`).
- Mean absolute epistasis increased only weakly with fraction of peak
  (Spearman `0.151`), below the preregistered threshold of `0.2`.
- The constraining-epistasis predictions failed. Negative epistasis decreased
  near peak-proximal backgrounds rather than increasing (Spearman `-0.6609`;
  low mean `0.01416`, high mean `0.0008789`), and sign epistasis was absent in
  the sampled pairs.
- The unconditioned condition behaved more intuitively than the combined
  result: near-peak backgrounds had lower beneficial fraction (`0.001072`) and
  smaller greedy-walk gain (`3.586`) than low-fitness backgrounds
  (`0.01993` beneficial fraction and `68.22` total gain).
- The synthetic-alignment condition confounds the peak-proximity axis. Its
  near-peak bin averaged fraction of stored peak fitness `1.748`, so the
  stored peak is not functioning as an effective upper bound for that branch.
- The strongest supported inference is therefore narrow: adaptive walks
  decelerate and beneficial availability tends to fall with higher-fitness
  backgrounds, but the peak-proximal epistasis component of `HYP-004` is not
  validated by this run.

## Result records created

- `RES-020`

## Hypothesis updates

- `HYP-004` is weakened as stated. The experiment preserves a narrower
  adaptive-walk and beneficial-availability result, but it does not support the
  preregistered claim that epistasis becomes more constraining as backgrounds
  approach the stored peak.
