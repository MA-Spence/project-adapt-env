# RES-027: RUN-069 fixed residual scale sweep validates the scale knob but not peak-stage unmasking

## Summary

RUN-069 corrects the RUN-067 calibration setting by scanning only fixed walk-initial residual/additive scale values. The fixed contact residual scale grid is directional: higher scale monotonically increases held-out SE, RSE, local maxima, WHT/GFT high-order fraction, residual high-order exposure, and lowers random-walk lag-1 autocorrelation. However, the result still does not validate a clean peak-proximal or late-walk unmasking mechanism, because stage-wise SE/RSE and high-order spectra are not consistently late-dominated even though smooth beneficial-effect proxies decline along the walk. The run supports HYP-005 for the fixed residual-scale knob and weakens the stronger HYP-004 biological ruggedness claim.

## Generated from

- Analyses: ANA-027

## Relevant hypotheses

- Supports: HYP-005
- Weakens: HYP-004
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/metrics.json`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_key_metrics.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_graphlet_sampling_audit.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_diminishing_returns_by_walk_stage.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_fixed_scale_test_summary.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_contact_calibration_test_transfer.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_fixed_scale_by_walk_stage.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_relative_scale_correlations.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_run067_comparison.md`
- `analyses/2026-05-15_ANA-027_run-069-tem1-fixed-scale-residual-ratio-sweep-graphlet-review/tables/run-069_validation_tests.md`
- `results/RES-027_run-069-fixed-residual-scale-sweep-validates-scale-knob-but-not-peak-stage-unmasking/metrics.json`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/summary.json`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/module_parameter_sweep.csv`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/graphlet_catalog.csv`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/calibration_graphlet_metrics.csv`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/test_graphlet_metrics.csv`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/graphlet_walsh_terms.csv`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/graphlet_gft_coefficients.csv`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069/graphlet_accessibility.csv`
- `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069.dvc`
- `experiments/2026-05-15_EXP-027_tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-027_tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/runs/RUN-069.yaml`

## Interpretation

- `RUN-069` completed successfully and produced the intended fixed-scale ratio
  sweep dataset: `30` dense graphlets, `600` graphlet/model metric rows,
  `7200` Walsh rows, `7200` Hamming-graph GFT rows, and no sign-epistasis pair
  rows. The collected output is about `4.9` MB and is DVC-tracked.
- The reduced graphlet cap no longer creates the sampling defect seen in
  `RUN-065`. Graphlets are balanced across five start classes (`6` each) and
  early, middle, and late walk stages (`10` each).
- The smooth core again shows diminishing-return proxies along the sampled
  adaptive walks: mean available beneficial effect declines from `24.56` early
  to `22.62` middle to `20.90` late, while best available effect declines from
  `213.70` to `174.33` to `162.71`.
- The smooth core still does not itself create local high-order topology:
  held-out SE is `0`, RSE is `0`, mean local maxima is `1.00`, and Walsh
  high-order fraction is effectively zero (`2.32e-29`).
- The fixed contact residual scale grid is now directional. Across fixed scales
  `1`, `4`, and `16`, held-out SE increases `0.2391 -> 0.2935 -> 0.4783`,
  RSE increases `0.0417 -> 0.0978 -> 0.1449`, local maxima increase
  `1.22 -> 1.48 -> 1.91`, and WHT/GFT high-order fraction increases
  `0.00074 -> 0.00439 -> 0.05055`.
- Calibration selects `contact_residual_ensemble__pset-002`, the highest fixed
  scale in the grid. This directly fixes the central interpretive problem in
  `RUN-067`, where the best contact residual used local-DFE scaling rather than
  the intended fixed walk-initial beneficial-effect scale.
- The selected fixed scale has the desired low-order-dominated but nonzero
  high-order spectrum: GFT/Walsh low, mid, and high fractions are `0.8447`,
  `0.1048`, and `0.0505`. This looks less like an over-rugged residual-only
  stress setting and more like a smooth background with a high-order tail.
- Random-walk lag-1 autocorrelation declines monotonically with fixed scale
  (`0.4560 -> 0.4203 -> 0.3399`). Accessible-path fraction is lower than the
  smooth core but not strictly monotone across scale (`0.1141`, `0.0960`,
  `0.1105`).
- The stronger late-stage unmasking prediction is not cleanly supported. At
  scale `16`, residual-to-available-beneficial-effect ratio rises from `3.83`
  early to `4.72` middle to `6.99` late, but SE is `0.4667`, `0.4833`, and
  `0.4458`, RSE is `0.2250`, `0.1167`, and `0.1083`, and WHT high-order
  fraction is `0.0572`, `0.0275`, and `0.0457`. The relative scale ratio
  increases late, but topology and spectra are not late-dominated.
- Compared with `RUN-067`, the corrected fixed scale sweep is a real rescue of
  the fixed-scale knob. The fixed contact setting in `RUN-067` had held-out SE
  `0.1409`, RSE `0.0216`, local maxima `1.18`, WHT high-order fraction
  `0.0043`, and residual high-order/base order-1 ratio `0.0046`; the selected
  fixed scale in `RUN-069` has SE `0.4783`, RSE `0.1449`, local maxima `1.91`,
  WHT high-order fraction `0.0505`, and ratio `0.0823`.
- The selected GNK null is not the stronger comparator in this run. It gives
  SE `0.4112`, RSE `0.0815`, local maxima `1.52`, and WHT high-order fraction
  `0.0033`, below the selected fixed contact residual on those metrics.
- Contact specificity remains unresolved. DCA pair override reaches almost the
  same held-out SE (`0.4819`) without high-order Walsh signal, while the
  selected random-position residual has comparable RSE (`0.1522`) and local
  maxima (`2.00`). The result supports scale control more strongly than it
  supports DCA/contact placement as the necessary biological cause.

## Effect on hypothesis

- `HYP-005` is supported for the fixed residual-scale knob. The intended knob
  now moves SE, RSE, local maxima, high-order WHT/GFT fraction, residual
  high-order exposure, and random-walk autocorrelation in the expected
  direction.
- `HYP-004` remains weakened as a broad biological ruggedness claim. The run
  supports a relative-scale mechanism but does not validate clean
  peak-proximal or late-walk unmasking, does not identify distance to a true
  biological peak, and does not establish contact placement as the causal
  biological structure.

## Limitations

- Only one landscape seed was used. The result establishes a directional knob
  in this diagnostic setting, not seed-robust population behaviour.
- The calibration set is small (`7` graphlets) because the reduced graphlet
  budget and held-out-start split leave most graphlets in test. The selected
  high fixed scale is plausible and supported by test trends, but calibration
  selection itself is noisy.
- Adaptive-walk stage and distance to the capped walk endpoint are proxies, not
  a known distance to a true local or global biological peak.
- Sign-epistasis pair rows were intentionally suppressed. SE/RSE fractions are
  available in graphlet metric rows, but pair-level causal inspection is not
  available from this run.
- Walsh and GFT term tables are capped. Full order-level spectral summaries are
  retained in the metric tables.
- The output is smaller than `RUN-067`, not identical in absolute storage. It
  satisfies the intended reduced-budget constraint and stays under the previous
  storage footprint, but the exact total storage is lower because the local-DFE
  primary settings were removed.
- Contact-placement controls still confound biological interpretation. The
  fixed-scale residual generator is now useful, but it is not yet a validated
  contact-conditioned procedural model of real protein ruggedness.

## Downstream use

- Use `RES-027` as evidence that the fixed residual/additive scale knob is
  operational and directional when scanned explicitly.
- Use `RES-027` to correct `RES-026`: the failure in `RUN-067` was partly a
  calibration/design issue, not proof that fixed residual scaling cannot create
  graphlet topology.
- Do not cite `RUN-069` as evidence that adaptive-walk stage alone captures
  distance to a fitness peak, or that diminishing returns automatically imply
  late-stage SE/RSE.
- Do not cite `RUN-069` as proof of contact-specific biological ruggedness.
  Matched controls remain important counterevidence.
