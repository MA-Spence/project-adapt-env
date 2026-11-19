# ANA-057: RUN-142 generator readout-coupling calibration review

## Purpose

Collect and interpret the completed RUN-142 (EXP-058) generator readout-coupling
calibration outputs against the preregistered EXP-058 core-guardrail, GNK-separation,
and knob-monotonicity gates, and apply the preregistered decision rule without adding
new post hoc success criteria.

## Linked experiments/runs

- Experiments: EXP-058
- Runs: RUN-142

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-142_summary.json`
- `source_tables/run-142_preflight.json`
- `tables/metrics.json`
- `tables/run-142_calibration_decision_table.csv`
- `tables/run-142_core_guardrail_summary_by_fix.md`
- `tables/run-142_gnk_separation_summary_by_fix.md`
- `tables/run-142_monotonicity_by_metric.md`
- `tables/run-142_high_plateau_collapse_fixA.md`

## Collection policy

`RUN-142` completed on `lab-slurm` as Slurm job `2003` (attempt 3, resumed from the
attempt-2 checkpoint directory) with inferred scheduler exit code `0`. The compact
runtime outputs were collected locally and DVC-tracked at:

- `data/processed/generator-readout-coupling-calibration/RUN-142/outputs.dvc`
- `data/processed/generator-readout-coupling-calibration/RUN-142/logs.dvc`

The per-parameter checkpoint tree was not copied into the scientific record; it remains
on the runtime host under
`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-142/outputs/checkpoints`.

## Execution audit

This is the full (`quick = false`) preregistered A/B/C sweep. The run evaluated all
`5022` parameter sets across `10` empirical templates, emitting `1,325,808` generator
descriptor rows and the ten expected output tables. Parameter-set counts match the
config grid: Fix A `972`, Fix B `3240`, Fix C `810`. `RUN-142` is therefore a completed
run, not a failed execution. The failure recorded below is scientific.

## Preregistered decision (calibration_decision_table.csv)

Every fix fails every registered gate; no candidate satisfies the joint requirement.

| fix | architecture | core-guardrail pass | GNK-separation pass | joint pass | monotonicity pass | all gates pass |
|---|---|---|---|---|---|---|
| A | single stability gate | 0/972 | 0/972 | 0 | 4421/9720 | False |
| B | mixed score/stability routing | 0/3240 | 0/3240 | 0 | 5151/38880 | False |
| C | two-channel product | 0/810 | 0/810 | 0 | 2311/9720 | False |

Because no fix (A, then A+B, then C) opens a passing region, the runner emits the
terminal preregistered branch:
`operator_free_combinatorial_realness_claim_fails_for_this_model_family`, and the hard
stop holds: do not rerun Stage 4C/4D/4E until the GNK separation gate passes.

## Core-guardrail gate (the binding constraint)

The guardrail fails identically for all `5022` parameter sets, and it fails on one
component: the high plateau never forms.

- low plateau present: `5022/5022` (100%);
- high plateau present: `0/5022` (0%);
- beneficial-fraction decline `rho < -0.8`: `175/5022` (~3.5%);
- rows with both plateaus present: `0`.

The observed high-plateau-regime fitness never approaches the gate. Against a
`fitness_output_scale` of `10.0` and a required `high_plateau_fitness_min >= 9.0`, the
best observed high-plateau median fitness across the entire grid is `0.0491`, and the
grid mean is `0.0032` (see `run-142_high_plateau_collapse_fixA`). In the re-centred
margin regime (`-2 ... +3`) the entire tested landscape sits near the low plateau: the
saturating observable is operated below its fold transition, so no high-fitness plateau
exists to support diminishing-returns behaviour. This is the direct, quantified
consequence of the RES-056 structural tension: the low-margin regime that exposes
epistasis is disjoint from the high-margin regime that forms a saturated plateau.

## GNK-separation gate

Because the gate requires separation *at matched core guardrail*, and the guardrail
never passes, `gnk_separation_pass = 0/5022` for all fixes. The unconditional
separation magnitudes (ignoring the guardrail conjunction) are the informative part and
differ sharply by architecture:

- Fix B (mixed routing) transmits large specific epistasis to the observable:
  `3240/3240` sets clear the high-order-fraction floor and `3240/3240` clear the
  sign-epistasis floor, with `|full - core|` high-order deltas up to `0.866` (median
  `0.708`) and sign-epistasis deltas up to `0.162`. This is the first clean evidence in
  the series that the readout-coupling defect is *fixable* for high-order and
  sign-epistasis content at protein scale.
- Fix A (single gate) is marginal: `25/972` high-order and `24/972` sign-epistasis
  passes, deltas peaking at `0.055`/`0.037`, i.e. right at the `0.03` floor.
- Fix C (two-channel product) transmits almost nothing: `0/810` on every separation
  metric, with a maximum high-order delta of `0.0021`. The product multiplies the
  collapsed near-zero stability gate, so `~0 x function_channel ~= 0`, and the function
  channel's GNK epistasis is annihilated at the observable.
- local-maxima-density separation fails in every architecture: `0/5022`, maximum
  observed delta `~7e-4` against the `5e-3` floor.

## Knob-monotonicity gate

No fix reaches a majority of checks with no wrong-direction medians
(`11,883/58,320 = 20.4%` overall). Per-metric medians (`run-142_monotonicity_by_metric`):

- Fix A: `sign_epistasis_fraction` and `reciprocal_sign_epistasis_fraction` are cleanly
  monotone in the expected direction (`99.7%` and `75.6%` pass, median `rho ~ 1.0` and
  `0.99`), but `high_order_fraction` moves the *wrong* way (median `rho -0.49`, `71%`
  wrong-direction) — more GNK reduces observable high-order content, the saturating
  low-pass signature. This wrong-direction median alone fails the gate for A.
- Fix B: broadly non-monotone (`12-16%` pass per metric; high-order and sign medians
  negative). Sweeping `gnk_epistasis_weight` and `lambda_score` jointly does not yield a
  monotone observable trend even though B transmits large-magnitude epistasis.
- Fix C: `high_order_fraction` (median `rho 0.30`) and `sign_epistasis_fraction` (median
  `rho 0.70`) trend correctly but pass weakly (`38%`, `42%`); `local_maxima_density`
  reverses (median `rho -0.58`).
- Across all fixes the topographic/accessibility knobs — `local_maxima_density`,
  `random_walk_lag1_autocorrelation`, `peak_targeted_reachability_fraction` — are poorly
  or wrongly controlled.

## Decision

Applying the preregistered decision rule: Fix A does not open a passing region, A+B do
not satisfy core guardrail and GNK separation together, and Fix C does not satisfy the
gates either. The terminal branch therefore fires —
`operator_free_combinatorial_realness_claim_fails_for_this_model_family` — and the hard
stop stands: Stage 4C/4D/4E reruns remain blocked for this implementation until the GNK
separation gate passes.

The correct reading is narrower than "the model has no epistasis" and more precise than
RES-055. RUN-142 localizes the failure to a single structural constraint: within the
tested grid, no single-observable architecture simultaneously (i) places the landscape
in a regime with a high-fitness plateau and peak-proximal diminishing returns and (ii)
transmits GNK-specific epistasis to the observable. Mixed routing (B) solves (ii)
cleanly but on a collapsed, plateau-less landscape; the product (C) inherits the
stability gate's collapse and solves neither. This is the RES-056 saturation/epistasis
tension confirmed and quantified at protein scale.
