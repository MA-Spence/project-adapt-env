# RES-026: RUN-067 relative-scale residual exposure adds diagnostics but not a fixed-scale peak-proximal SE/RSE mechanism

## Summary

The scaled relative-scale graphlet repeat fixed RUN-065's sampling defect and showed that residual exposure correlates with high-order spectra and local topology, but fixed walk-initial residual scaling did not produce a clean peak-proximal SE/RSE increase and matched controls still weaken contact-specific biological interpretation.

## Generated from

- Analyses: ANA-026

## Relevant hypotheses

- Supports: HYP-005
- Weakens: HYP-004
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/metrics.json`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_key_metrics.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_graphlet_sampling_audit.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_diminishing_returns_by_walk_stage.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_selected_test_model_comparison.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_contact_calibration_test_transfer.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_walk_stage_relative_scale.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_primary_by_start_class.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_relative_scale_correlations.md`
- `analyses/2026-05-15_ANA-026_run-067-tem1-relative-scale-residual-exposure-graphlet-review/tables/run-067_validation_tests.md`
- `results/RES-026_run-067-relative-scale-residual-exposure-adds-diagnostics-but-not-fixed-scale-peak-proximal-serse/metrics.json`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/summary.json`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/module_parameter_sweep.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_catalog.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/calibration_graphlet_metrics.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/test_graphlet_metrics.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_walsh_terms.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_gft_coefficients.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_accessibility.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067.dvc`
- `experiments/2026-05-15_EXP-026_tem1-relative-scale-residual-exposure-graphlet-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-026_tem1-relative-scale-residual-exposure-graphlet-diagnostic/runs/RUN-067.yaml`

## Interpretation

- `RUN-067` completed successfully and produced the intended scaled dataset:
  `90` dense graphlets, `2250` model-graphlet metric rows, `27000` Walsh rows,
  `27000` Hamming-graph GFT rows, and no sign-epistasis pair rows. The durable
  output is about `15.8` MB and is tracked by DVC.
- The main technical defect in `RUN-065` is fixed. Graphlets are exactly
  balanced across the five requested start classes (`18` each) and across
  early, middle, and late adaptive-walk stages (`30` each).
- The smooth thermodynamic core shows diminishing-return proxies along the
  sampled walks. Mean available beneficial effect declines from `25.01` early
  to `20.29` middle to `17.36` late, and best available effect declines from
  `245.06` to `105.92` to `82.69`.
- The smooth core does not itself create SE/RSE or high-order topology in
  these graphlets. Held-out smooth-core SE is `0`, RSE is `0`, mean local
  maxima is `1.00`, and Walsh high-order fraction is effectively zero
  (`3.55e-28`).
- This directly supports the chat correction: diminishing returns in the
  thermodynamic/additive core should not be treated as automatically implying
  WHT high-order signal, sign epistasis, or reciprocal sign epistasis. In this
  run, diminishing-return proxies occur in the smooth core while SE/RSE remain
  absent.
- The selected contact residual does create held-out topology over the smooth
  core: SE `0.3511`, RSE `0.1481`, mean local maxima `2.56`,
  accessible-path fraction `0.0323`, and lag-1 random-walk autocorrelation
  `0.3634`.
- The selected contact residual also has the intended qualitative spectral
  shape: low-order dominated but not purely low-order. Held-out Walsh
  high-order fraction is `0.0855`, GFT low/mid/high energy fractions are
  `0.7197`, `0.2603`, and `0.0200`, and residual high-order to base order-1
  variance ratio is `0.2378`.
- However, the fixed-scale unmasking mechanism is not validated. The fixed
  walk-initial contact residual has weak held-out topology (`SE = 0.1409`,
  `RSE = 0.0216`, mean local maxima `1.18`) and does not show a clean
  monotone increase in SE/RSE across early, middle, and late graphlets.
- Calibration selected `contact_residual_ensemble__pset-003`, which uses
  `local_dfe_sd` scaling, not the fixed `constant_walk_initial_beneficial_effect`
  scale. That is a material result: the diagnostic machinery is useful, but
  the current model has not yet made the fixed residual scale into the primary
  successful mechanism.
- A weaker relative-scale statement is supported. Across held-out graphlets
  for the selected contact residual, residual high-order exposure relative to
  the base order-1 variance correlates strongly with Walsh high-order fraction
  (`Spearman rho = 0.980`) and moderately with RSE (`0.482`) and local maxima
  (`0.531`).
- Compared with `RUN-065`, the topology signal persists under the corrected
  multi-start sampling design. `RUN-065` had held-out SE `0.3213`, RSE
  `0.1030`, local maxima `3.54`, accessibility `0.0350`, lag-1 autocorrelation
  `0.4976`, and Walsh high-order fraction `0.0737` in reference-like graphlets
  only. `RUN-067` gives SE `0.3511`, RSE `0.1481`, local maxima `2.56`,
  accessibility `0.0323`, lag-1 autocorrelation `0.3634`, and Walsh high-order
  fraction `0.0855` across all start classes. This rescues the sampling
  generalization flaw, but not the fixed-scale mechanistic claim.
- Matched controls still weaken contact-specific interpretation. The primary
  contact residual is stronger than shuffled-contact and random-position
  residual controls, but those controls remain close enough that DCA/contact
  placement is not established as the causal source of the topology.
- GNK remains the stronger non-procedural comparator: held-out SE `0.6218`,
  RSE `0.2710`, mean local maxima `5.02`, lag-1 autocorrelation `0.0940`,
  Walsh high-order fraction `0.3511`, and residual high-order to base order-1
  variance ratio `3.49`.

## Effect on hypothesis

- `HYP-005` is supported narrowly. The revised graphlet sampler and residual
  exposure instrumentation behave directionally: residual exposure metrics
  track high-order spectra and local ruggedness, and the scaled design stays
  under the intended data-volume and graphlet limits.
- `HYP-004` remains weakened as a broad biological ruggedness claim. The run
  shows diminishing-return proxies and residual-induced local topology, but it
  does not show that diminishing returns alone creates SE/RSE, does not validate
  a fixed-scale peak-proximal unmasking mechanism, and does not establish
  contact-specific biological placement.

## Limitations

- Distance to a true global or local biological peak is still unknown. The run
  uses adaptive-walk stage and distance to the capped walk endpoint as
  observable proxies, not as a certified peak coordinate.
- The run uses one landscape seed and `90` dense graphlets. It is deliberately
  scaled for fast review and is not a large-seed robustness study.
- The calibration/test split covers all start classes but is not itself
  perfectly balanced within every start-class/stage cell.
- Sign-epistasis pair rows were intentionally suppressed to keep output small.
  SE/RSE fractions are present in graphlet metric rows, but pair-level
  inspection is not available from this run.
- Walsh and GFT term tables record only the top terms per graphlet/model. Full
  order-level spectral summaries are retained in the metric tables.
- The selected primary contact residual uses local-DFE scaling, which partly
  reintroduces graphlet-local scale dependence. That makes it a successful
  topology knob, not clean evidence for fixed absolute residual structure being
  unmasked near peaks.
- Contact placement remains confounded with strength, density, coefficient
  distribution, and graphlet context because matched placement controls remain
  close to the primary model.

## Downstream use

- Use `RES-026` as evidence that the corrected aggregate graphlet system can
  sample balanced multi-start graphlets under a small output budget and can
  record relative-scale diagnostics reproducibly.
- Use `RES-026` as support for the corrected causal language: diminishing
  returns changes the scale of smooth gains, but SE/RSE require residual or
  interaction structure that becomes relevant relative to that smooth
  background.
- Do not cite `RUN-067` as proof that fixed absolute high-order residuals are
  already the operative mechanism. The fixed-scale setting failed the stronger
  monotone SE/RSE test, and calibration preferred local-DFE scaling.
- Do not cite `RUN-067` as proof of contact-specific biological ruggedness.
  Matched controls and the stronger GNK null remain important counterweights.
