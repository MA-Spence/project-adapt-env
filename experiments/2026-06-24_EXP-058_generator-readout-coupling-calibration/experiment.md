# EXP-058: Generator Readout-Coupling Calibration

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005, HYP-009
- Previous experiments: EXP-051, EXP-052, EXP-053, EXP-054, EXP-055, EXP-056, EXP-057
- Motivating results: RES-049, RES-050, RES-051, RES-052, RES-053, RES-054, RES-055, RES-056
- Scratch source: `scratch/RES-056_assessment_readout_coupling_and_next_steps.md`

## Question

Can the public additive thermodynamic core plus GNK generator be calibrated so
the observable simultaneously preserves core neutral-network/diminishing-return
guardrails and transmits GNK-specific epistasis through the readout?

## Pre-experiment prediction

A desaturated and re-centred single-gate readout may expose a feasible region.
If it cannot satisfy core guardrails and GNK separation together, mixed
post-readout score routing or a two-channel product observable will be required.

## Rationale

`RES-055` removed the measurement artifact from `RES-054` and made the current
failure diagnostic: the public observable routes GNK into the stability margin
and then applies a steep, near-step `stability_gate` readout while the tested
margin grid lies deep in the folded plateau. GNK is therefore present in the
latent generator but nearly annihilated at the observable. This is the same
observable-collapse failure mode seen in the PHOT calibration series, now
isolated forward and without a PHOT measurement operator.

`RES-056` is not the driver for this experiment. It did not recover landscapes
or estimate a transfer function; it records an EXP-050 tooling/certification
failure in the non-additive MAVE-NN route. Component 1 is the operator-free
combinatorial generator route, so EXP-058 keeps EXP-050 off the critical path.

## Experimental design

This is the single generator calibration experiment proposed in
`scratch/RES-056_assessment_readout_coupling_and_next_steps.md`. It is not a
Stage 4C/4D/4E rerun, not an empirical-overlap claim, and not an EXP-050 repair.

EXP-058 sweeps exactly the three proposed fix levers:

- Fix A: operate the single stability gate in a non-saturated regime by lowering
  `threshold_steepness` to the 1-4 range, re-centring the stability-margin grid
  near the fold transition (`-2` through `+3`), and raising GNK amplitude
  relative to the core margin spread.
- Fix B: route a fraction `lambda_score` of GNK post-readout through the existing
  score-space path, while routing `1 - lambda_score` through stability. This
  tests whether mixed score/stability routing can expose specific epistasis
  without losing threshold-coupled realism.
- Fix C: test the two-channel product observable, carrying core
  neutral-network/diminishing-return behaviour on the stability gate and
  GNK-specific epistasis on a non-saturated functional channel combined
  multiplicatively.

The run must evaluate each candidate at matched core guardrail, not by GNK
separation alone. It must compare the full observable against additive-core-only
and GNK-relevant controls, and it must report single-gate, mixed-routing, and
two-channel candidates in one decision table.

## Registered gates

Core guardrail:

- low and high plateaus are both present in the chosen regime;
- beneficial fraction declines toward peak with Spearman `rho < -0.8`.

GNK separation:

- at matched core guardrail, mean `abs(full - core)` must reach the registered
  effect floor for high-order fraction, sign-epistasis fraction, and local
  maxima density;
- the high-order fraction floor is at least `0.03`, matching the assessment
  document;
- sign-epistasis and local-maxima floors are explicitly registered in
  `config.yaml`.

Knob monotonicity:

- a majority of checks pass with `abs(rho) >= 0.5`;
- no median trend for a registered metric moves in the wrong direction.

## Decision rule

- If Fix A alone opens a region satisfying all three gates, the single-latent
  composite is salvageable and can proceed to a later empirical rerun.
- If A needs B, adopt mixed routing and record `lambda_score` as a model
  parameter.
- If A+B cannot satisfy core guardrail and GNK separation together, adopt Fix C
  as the pivot result: a single saturating readout cannot carry both global and
  specific epistasis in this implementation.
- If even Fix C cannot satisfy the gates, the operator-free combinatorial
  realness claim for this model family fails under this test.

Hard stop: do not rerun Stage 4C, Stage 4D, or Stage 4E until EXP-058 passes the
GNK separation gate. Without separation, empirical overlap remains
uninterpretable.

## Inputs

- `scratch/RES-056_assessment_readout_coupling_and_next_steps.md`
- `results/RES-055_run-141-fixed-support-aware-probe-fails-core-gnk-validation-gates/result.md`
- `results/RES-056_run-129-nonadditive-phot-operators-fail-certification/result.md`
- `scripts/generator_readout_coupling_calibration_exp058.py`
- `scripts/design_matched_support_aware_probe_exp057.py`
- `scripts/ensemble_landscape_validation_stage4.py`
- `external/Adapt-Env/adaptenv/landscape.py`
- `data/processed/ensemble-landscape-validation-datasets.dvc`

## Configuration

- `config.yaml` preregisters the A/B/C sweep, gates, and decision rule.

## Execution

- Initial run record: RUN-142.
- Attempt 1: submitted to `lab-slurm` as Slurm job `1991`; failed by walltime
  before producing final outputs.
- Attempt 2: resubmitted to `lab-slurm` as Slurm job `1999` after runner
  performance/checkpoint repair. It reached `297/5022` checkpointed parameter
  sets with clean stderr, but observed throughput showed the original
  `48:00:00` walltime remained too short, so it was cancelled and reconciled as
  failed.
- Attempt 3: resubmitted to `lab-slurm` as Slurm job `2003` with a
  `168:00:00` walltime request and checkpoint resume from attempt 2; completed
  with scheduler exit code `0` after evaluating all `5022` parameter sets.
- The 11 preregistered outputs and three run logs were collected under
  `data/processed/generator-readout-coupling-calibration/RUN-142`. Outputs and
  logs are DVC-tracked separately and backed up to the configured `labcluster`
  remote. The transient per-parameter checkpoint tree was not collected.

## Expected outputs

- preflight.json
- generator_calibration_grid.csv
- support_catalog.csv
- support_aware_graphlet_audit.csv
- anchor_regime_audit.csv
- generator_descriptor_metrics.csv
- core_guardrail_checks.csv
- gnk_separation_checks.csv
- knob_monotonicity_checks.csv
- calibration_decision_table.csv
- summary.json

## Analysis plan

After live execution, create one ANA record that reports the pass/fail status of
the three gates and applies the preregistered decision rule without adding new
post hoc success criteria.

## Completion criteria

- [x] Config recorded
- [x] Run script recorded
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- This planned experiment is intentionally operator-free and should not estimate
  or use an EXP-050 MAVE-NN transfer function.
- Outcome (ANA-057 / RES-057): all three fixes fail all three registered gates
  (`0/5022` core-guardrail, `0/5022` GNK-separation, `11,883/58,320` monotonicity
  with wrong-direction medians). The terminal preregistered decision fires —
  `operator_free_combinatorial_realness_claim_fails_for_this_model_family` — and the
  hard stop holds: Stage 4C/4D/4E reruns remain blocked. The binding constraint is
  that no high-fitness plateau forms in the re-centred margin regime (best observed
  high-plateau median fitness `0.049` on a `0-10` scale), so plateau-based diminishing
  returns and transmitted GNK-specific epistasis are not obtained together. Fix B mixed
  routing does transmit large observable high-order/sign epistasis (deltas up to
  `0.87`/`0.16`) but on a collapsed, plateau-less landscape; Fix C two-channel product
  transmits almost nothing because it was swept only in the non-saturated margin regime,
  so its stability gate collapsed the product.
- The collected `summary.json` and `calibration_decision_table.csv` report
  `0/5022` core-guardrail passes and `0/5022` GNK-separation passes. This is a
  collection-stage observation only; the preregistered decision still requires
  a formal analysis and result record.
