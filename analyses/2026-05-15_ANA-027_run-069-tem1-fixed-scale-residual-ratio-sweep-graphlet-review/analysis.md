# ANA-027: RUN-069 TEM-1 fixed-scale residual ratio sweep graphlet review

## Purpose

Review the EXP-027 fixed residual/additive scale-grid repeat, including output budget, graphlet sampling, scale-wise SE/RSE and Walsh/GFT spectra, walk-stage unmasking, calibration/test transfer, and controls relative to RUN-067.

## Linked experiments/runs

- Experiments: EXP-027
- Runs: RUN-069

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-069_key_metrics.md`
- `tables/run-069_graphlet_sampling_audit.md`
- `tables/run-069_diminishing_returns_by_walk_stage.md`
- `tables/run-069_fixed_scale_test_summary.md`
- `tables/run-069_contact_calibration_test_transfer.md`
- `tables/run-069_fixed_scale_by_walk_stage.md`
- `tables/run-069_primary_by_start_class.md`
- `tables/run-069_relative_scale_correlations.md`
- `tables/run-069_fixed_scale_monotonicity.md`
- `tables/run-069_run067_comparison.md`
- `tables/run-069_validation_tests.md`

## Inputs

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
- `results/RES-026_run-067-relative-scale-residual-exposure-adds-diagnostics-but-not-fixed-scale-peak-proximal-serse/result.md`

## Analysis performed

- Audited `RUN-069` completion, output volume, DVC tracking, and graphlet sampling
  against the preregistered `EXP-027` reduced-budget design.
- Recomputed calibration/test model summaries from the exported graphlet metric
  rows rather than relying only on `module_parameter_sweep.csv`.
- Compared the fixed contact residual scale grid
  (`residual_scale_multiplier = 1, 4, 16`) for SE, RSE, local maxima,
  accessibility, random-walk autocorrelation, Walsh/GFT high-order fraction, and
  residual high-order exposure.
- Compared the selected fixed contact residual against the smooth core, DCA pair
  override, RES-024 exact-match haplotype control, shuffled-contact and
  random-position controls, GNK local null, and residual-only diagnostic.
- Tested the intended mechanism in two parts: the fixed residual scale knob, and
  the stronger late-walk/peak-proximal unmasking prediction.
- Compared the fixed-scale result with `RUN-067`, where calibration had selected
  a local-DFE-scaled contact residual rather than the intended fixed-scale
  mechanism.

## Main observations

- `RUN-069` completed successfully with scheduler exit code `0`. The collected
  output is DVC-tracked at
  `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069.dvc`.
- The run produced a compact dataset: `30` dense graphlets, `600`
  graphlet/model metric rows, `7200` Walsh rows, `7200` Hamming-graph GFT rows,
  no sign-epistasis pair rows, and about `4.9` MB of collected data.
- The graphlet cap is now behaving as intended. Graphlets are balanced across
  five start classes (`6` each) and three walk stages (`10` early, `10` middle,
  `10` late).
- The smooth thermodynamic core again shows diminishing-return proxies along
  the sampled adaptive walks. Mean available beneficial effect declines from
  `24.56` early to `22.62` middle to `20.90` late, and best available effect
  declines from `213.70` to `174.33` to `162.71`.
- The smooth core itself remains non-rugged in the sampled graphlets: held-out
  SE is `0`, RSE is `0`, mean local maxima is `1.00`, and Walsh high-order
  fraction is effectively zero (`2.32e-29`).
- The fixed contact residual scale grid is strongly directional in held-out
  graphlets. As scale increases from `1` to `4` to `16`, SE rises from
  `0.2391` to `0.2935` to `0.4783`; RSE rises from `0.0417` to `0.0978` to
  `0.1449`; local maxima rise from `1.22` to `1.48` to `1.91`; and WHT/GFT
  high-order fraction rises from `0.00074` to `0.00439` to `0.05055`.
- The selected contact residual is now the intended fixed-scale setting:
  `contact_residual_ensemble__pset-002`, the highest fixed scale in the grid.
  This corrects the central interpretive problem in `RUN-067`, where the
  selected contact residual used local-DFE scaling.
- At the selected fixed scale, the spectrum is low-order dominated with a
  nonzero high-order tail rather than high-order dominated: GFT/Walsh low,
  mid, and high fractions are `0.8447`, `0.1048`, and `0.0505`.
- The fixed scale also shortens random-walk autocorrelation directionally:
  lag-1 autocorrelation declines from `0.4560` to `0.4203` to `0.3399` across
  scales. Accessible-path fraction is lower than the smooth core but not
  strictly monotone across the three fixed scales (`0.1141`, `0.0960`,
  `0.1105`).
- The stronger late-stage unmasking claim remains weak. For the highest fixed
  scale, residual-to-available-beneficial-effect ratio rises from `3.83` early
  to `4.72` middle to `6.99` late, but SE is `0.4667`, `0.4833`, and `0.4458`,
  RSE is `0.2250`, `0.1167`, and `0.1083`, and WHT high-order fraction is
  `0.0572`, `0.0275`, and `0.0457`. The relative-scale ratio rises late, but
  topology and spectra are not cleanly late-dominated.
- Compared with the fixed contact setting in `RUN-067`, `RUN-069` shows a much
  stronger fixed-scale effect: SE increases from `0.1409` to `0.4783`, RSE from
  `0.0216` to `0.1449`, local maxima from `1.18` to `1.91`, and WHT high-order
  fraction from `0.0043` to `0.0505`.
- The selected fixed contact residual is not dominated by the selected GNK null
  in this run. Selected GNK gives SE `0.4112`, RSE `0.0815`, local maxima
  `1.52`, and WHT high-order fraction `0.0033`, while selected contact gives
  SE `0.4783`, RSE `0.1449`, local maxima `1.91`, and WHT high-order fraction
  `0.0505`.
- Contact specificity remains unresolved. DCA pair override reaches nearly the
  same SE (`0.4819`) without high-order Walsh signal, and random-position
  residuals have comparable RSE (`0.1522`) and local maxima (`2.00`) despite a
  smaller WHT high-order tail. This means the result supports the fixed scale
  knob more than it supports DCA/contact placement as a necessary biological
  source.

## Interpretation

`RUN-069` materially changes the interpretation of the fixed residual/additive
scale mechanism relative to `RUN-067`. The previous result did not test the
intended calibration cleanly because calibration selected a local-DFE-scaled
contact residual. Here, the primary contact residual grid uses only fixed
walk-initial beneficial-effect scaling, and calibration selects the highest fixed
scale. That is a real positive result for `HYP-005`: the fixed scale knob is
directional and moves the intended graphlet statistics.

The result does not rescue the stronger `HYP-004` peak-proximal story. The
smooth core shows diminishing-return proxies, and the fixed residual scale can
make hidden high-order terms relevant, but the stage-wise data do not show a
simple late-walk monotone increase in SE/RSE or WHT high-order fraction. In
other words, the corrected causal chain is partly supported as a relative-scale
mechanism, but adaptive-walk stage is still too weak and noisy a proxy for
distance to a true biological peak.

This is also not yet strong evidence for contact-conditioned biological
specificity. The primary contact residual has a better low-order-dominated
high-order-tail spectrum than residual-only stress settings and is not beaten by
the selected GNK null, but matched controls remain close on several topology
metrics. The supported claim should therefore be limited to fixed-scale residual
control of graphlet ruggedness, not a full procedural biological ruggedness
validation.

## Result records created

- `RES-027`

## Hypothesis updates

- `HYP-005` is supported for the fixed residual-scale knob. The scale grid moves
  SE, RSE, local maxima, WHT/GFT high-order fraction, residual high-order
  exposure, and autocorrelation in the expected direction.
- `HYP-004` remains weakened as a broad biological ruggedness claim. The run
  does not validate clean peak-proximal or late-walk unmasking, and contact
  placement is not established as the causal biological structure.
