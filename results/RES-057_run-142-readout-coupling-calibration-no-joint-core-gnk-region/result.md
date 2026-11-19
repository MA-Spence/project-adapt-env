# RES-057: RUN-142 readout-coupling calibration finds no joint core-guardrail and GNK-separation region

## Summary

RUN-142 ran the full preregistered EXP-058 A/B/C readout-coupling sweep (5022 parameter
sets across 10 empirical templates). All three fixes fail all three registered gates:
0/5022 core-guardrail passes, 0/5022 GNK-separation passes, and no fix reaches a
monotonicity majority without wrong-direction medians. The terminal preregistered
decision fires — the operator-free combinatorial realness claim fails for this model
family under this test — and the hard stop holds: Stage 4C/4D/4E reruns stay blocked
until the GNK separation gate passes. The failure is localized to a single structural
constraint: no single-observable architecture simultaneously forms a high-fitness
plateau and transmits GNK-specific epistasis.

## Generated from

- Analyses: ANA-057

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004, HYP-005, HYP-009
- Refutes: None

## Evidence

- analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review/tables/metrics.json
- analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review/tables/run-142_calibration_decision_table.csv
- analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review/tables/run-142_core_guardrail_summary_by_fix.csv
- analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review/tables/run-142_gnk_separation_summary_by_fix.csv
- analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review/tables/run-142_monotonicity_by_metric.csv
- analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review/tables/run-142_high_plateau_collapse_fixA.csv
- analyses/2026-07-22_ANA-057_run-142-generator-readout-coupling-calibration-review/source_tables/run-142_summary.json
- data/processed/generator-readout-coupling-calibration/RUN-142/outputs.dvc
- data/processed/generator-readout-coupling-calibration/RUN-142/logs.dvc
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-142/outputs

## Interpretation

`RUN-142` completed cleanly on `lab-slurm` as job `2003` (attempt 3, resumed from the
attempt-2 checkpoint), evaluating all `5022` parameter sets and emitting `1,325,808`
generator descriptor rows and the ten expected tables. This is the full
(`quick = false`) sweep of the three RES-056 fixes: Fix A desaturated/re-centred single
gate, Fix B mixed score/stability routing, Fix C two-channel product.

The preregistered decision is unambiguous. No fix passes any gate jointly:
core-guardrail `0/5022`, GNK-separation `0/5022`, and monotonicity `11,883/58,320`
(`20.4%`) with wrong-direction medians present in every architecture. Per the decision
rule (A alone, then A+B, then C), the terminal branch fires:
`operator_free_combinatorial_realness_claim_fails_for_this_model_family`.

The binding constraint is the core guardrail, and it fails on one component uniformly:
the high plateau never forms (`0/5022`), while the low plateau always does
(`5022/5022`). Against a `fitness_output_scale` of `10.0` and a required
`high_plateau_fitness_min >= 9.0`, the best observed high-plateau median fitness across
the whole grid is `0.0491` (grid mean `0.0032`). In the re-centred margin regime
(`-2 ... +3`) chosen to expose epistasis, the saturating observable is operated below
its fold transition, so the entire tested landscape collapses toward the low plateau and
there is no high-fitness region to carry diminishing returns. This is the RES-056
saturation/epistasis tension confirmed and quantified: the low-margin regime that
exposes epistasis is disjoint from the high-margin regime that forms a plateau.

Two architecture-specific findings sharpen this beyond RES-055:

- Mixed routing (Fix B) does solve the epistasis-transmission problem. Every one of its
  `3240` sets clears the high-order-fraction and sign-epistasis floors, with
  `|full - core|` high-order deltas up to `0.866` (median `0.708`) and sign-epistasis
  deltas up to `0.162`. The readout-coupling defect diagnosed in RES-055/RES-056 is
  therefore fixable for high-order and sign-epistasis content at protein scale — but on a
  collapsed, plateau-less landscape, so it cannot satisfy the joint gate.
- The two-channel product (Fix C) transmits essentially nothing (`0/810` on every
  separation metric, maximum high-order delta `0.0021`). Because the product multiplies
  the collapsed near-zero stability gate, the function channel's epistasis is annihilated
  at the observable. RES-056 named Fix C the expected destination, but in this grid it
  never operated its stability gate in the saturated regime where the product would be
  well-behaved (see Limitations).

Local-maxima-density separation fails in every architecture (`0/5022`; maximum delta
`~7e-4` against a `5e-3` floor): none of the fixes move observable local-maxima
ruggedness. On monotonicity, Fix A's `sign_epistasis_fraction` and
`reciprocal_sign_epistasis_fraction` knobs are cleanly directional (`99.7%` and `75.6%`
pass, median `rho ~ 1.0`), but its `high_order_fraction` moves the wrong way (median
`rho -0.49`), which fails the no-wrong-direction-median clause on its own.

## Effect on hypothesis

- `HYP-005` is weakened for the current protein-scale public generator: no fix restores
  a majority of monotone, correct-direction knob effects, and topographic/accessibility
  knobs (local-maxima density, random-walk autocorrelation, peak reachability) are poorly
  or wrongly controlled. The narrow exception is Fix A's directional sign/reciprocal
  epistasis knobs, which are monotone but insufficient for the gate.
- `HYP-004` is weakened as an implemented public-generator realism claim: within the
  tested regime the observable cannot host a high-fitness plateau and peak-proximal
  diminishing returns together with transmitted specific epistasis, so biologically-like
  adaptive-walk dynamics are not supported by this observable.
- `HYP-009` is weakened as a benchmark-construction route: routing and optimizer
  comparisons require matched, controllable ruggedness at the observable, which this run
  does not deliver for any of the three architectures.

This result does not refute additive-core plus GNK in principle, and it is not evidence
that the model lacks epistasis. It is a negative result for the current single-observable
implementation and the specific EXP-058 A/B/C grid, with the failure now attributable to
the saturating-readout / plateau-formation tension rather than to the knob (validated at
`L=8` in RES-049), the probe (fixed in RUN-141), or the operator (RES-048/RES-056).

## Limitations

- This is an operator-free generator calibration, not a Stage 4C/4D/4E empirical-overlap
  or held-out test. It uses no EXP-050 MAVE-NN transfer function.
- Fix C was swept only over the same low re-centred margin grid (`-2 ... +3`) as A and B,
  so its stability gate was never operated in its saturated (`gate ~= 1`) plateau regime.
  The product observable therefore inherited the stability gate's collapse. The Fix C
  failure recorded here is a failure of Fix C *as configured in this grid*, not a general
  impossibility result for a two-channel product whose stability channel is operated in
  saturation while a non-saturated functional channel carries epistasis. This is a lead,
  not a new success criterion, and does not change the preregistered decision.
- The high-plateau guardrail threshold (`min fitness >= 9` on a `0-10` scale) was
  inherited from the deep-plateau Stage-4 regime; the re-centred margin grid is by
  construction unable to reach it. The joint infeasibility is real for this grid, but a
  future probe should co-vary the operating point and the plateau threshold rather than
  fixing one against the other.

## Downstream use

- Do not run or cite Stage 4C/4D/4E reruns for this implementation; the hard stop holds.
- Do not treat the single-observable additive-core plus GNK composite as validated for
  empirical benchmark landscape construction.
- Cite RUN-142 as the decisive readout-coupling result: mixed post-readout routing (Fix
  B) transmits large observable high-order and sign epistasis at protein scale, but no
  single-observable architecture in the tested grid combines that with a high-fitness
  plateau and peak-proximal diminishing returns. The natural next probe operates the
  two-channel stability gate in saturation (high margin, `gate ~= 1`) while a
  non-saturated functional channel carries GNK epistasis — the CLAIM-002 multi-latent
  direction — but that is a new experiment, not a reinterpretation of this one.
