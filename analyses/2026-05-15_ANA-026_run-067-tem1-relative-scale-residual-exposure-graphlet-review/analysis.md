# ANA-026: RUN-067 TEM-1 relative-scale residual exposure graphlet review

## Purpose

Review the scaled EXP-026 relative-scale residual exposure graphlet experiment, including balanced graphlet sampling, diminishing-return proxies, fixed versus local residual scaling, Walsh/GFT spectra, topology, and controls.

## Linked experiments/runs

- Experiments: EXP-026
- Runs: RUN-067

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-067_key_metrics.md`
- `tables/run-067_graphlet_sampling_audit.md`
- `tables/run-067_diminishing_returns_by_walk_stage.md`
- `tables/run-067_selected_test_model_comparison.md`
- `tables/run-067_contact_calibration_test_transfer.md`
- `tables/run-067_walk_stage_relative_scale.md`
- `tables/run-067_primary_by_start_class.md`
- `tables/run-067_relative_scale_correlations.md`
- `tables/run-067_validation_tests.md`

## Inputs

- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/summary.json`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/module_parameter_sweep.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_catalog.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/calibration_graphlet_metrics.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/test_graphlet_metrics.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_walsh_terms.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_gft_coefficients.csv`
- `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067/graphlet_accessibility.csv`
- `experiments/2026-05-15_EXP-026_tem1-relative-scale-residual-exposure-graphlet-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-026_tem1-relative-scale-residual-exposure-graphlet-diagnostic/runs/RUN-067.yaml`
- `results/RES-025_run-065-aggregate-residuals-create-local-topology-but-fail-broad-procedural-generalization/result.md`
- `results/RES-024_run-063-dca-clique-haplotype-terms-create-controllable-high-order-walsh-signal-but-limited-topology/result.md`

## Analysis performed

- Reconciled `RUN-067` as the completed live `EXP-026` execution on
  `lab-slurm` with scheduler job `227` and scheduler exit code `0`.
- Audited output scale and graphlet sampling balance against the scaled
  `EXP-026` design.
- Recomputed review metrics from exported CSV/JSON files using
  `notebooks/analysis.py`.
- Compared the smooth thermodynamic core, the selected contact residual,
  fixed walk-initial contact residuals, RES-024 exact-match haplotype control,
  RES-023 DCA pair override, shuffled-contact and random-position placement
  controls, GNK local null, and residual-only diagnostic.
- Quantified diminishing-return proxies along the smooth-core adaptive walks:
  available beneficial fraction, mean available beneficial effect, best
  available effect, background fitness, and distance to the capped walk
  endpoint.
- Tested the chat-level causal chain explicitly: whether residual/high-order
  components become relevant because their relative scale rises as smooth
  adaptive gains decline, rather than because diminishing returns directly
  creates SE/RSE.

## Main observations

- `RUN-067` completed successfully and produced a compact review dataset:
  `90` dense graphlets, `2250` model-graphlet metric rows, `27000` Walsh rows,
  `27000` GFT rows, and no sign-epistasis pair rows. The collected output
  directory is about `15.8` MB.
- The `RUN-065` graphlet-cap failure was fixed. Evaluated graphlets were
  exactly balanced across five start classes (`18` each) and three adaptive-walk
  stages (`30` early, `30` middle, `30` late).
- Smooth-core adaptive walks show the intended diminishing-return proxies:
  mean available beneficial effect declines from `25.01` early to `20.29`
  middle to `17.36` late, and best available effect declines from `245.06` to
  `105.92` to `82.69`.
- The smooth core itself is not locally rugged in these graphlets. Held-out
  smooth-core SE is `0`, RSE is `0`, mean local maxima is `1.00`, and Walsh
  high-order fraction is effectively zero (`3.55e-28`).
- The selected contact residual creates held-out topology over the smooth core:
  SE `0.3511`, RSE `0.1481`, mean local maxima `2.56`, endpoint accessible-path
  fraction `0.0323`, and lag-1 random-walk autocorrelation `0.3634`.
- The selected contact residual has a low-order-dominated but nonzero
  high-order spectrum: held-out Walsh high-order fraction `0.0855`, GFT
  low/mid/high fractions `0.7197`, `0.2603`, and `0.0200`, and residual
  high-order to base order-1 variance ratio `0.2378`.
- The fixed walk-initial contact residual did not support the strongest form
  of the relative-scale prediction. Its held-out topology remained weak
  (`SE = 0.1409`, `RSE = 0.0216`, mean local maxima `1.18`), and SE/RSE did not
  increase monotonically across walk stages even though residual exposure rose
  late.
- Calibration selected `contact_residual_ensemble__pset-003`, which uses
  `local_dfe_sd` scaling, not `constant_walk_initial_beneficial_effect`
  scaling. This is material: the run supports residual exposure as an
  informative diagnostic, but it does not show that the fixed-scale unmasking
  mechanism is currently implemented strongly enough to be the selected primary
  mechanism.
- Within the selected contact residual, relative high-order exposure is strongly
  associated with spectral and topology metrics across held-out graphlets.
  Spearman correlation between residual high-order to base order-1 variance
  ratio and Walsh high-order fraction is `0.980`; correlations with RSE and
  local maxima are `0.482` and `0.531`.
- Contact placement has only weak support. The primary contact residual
  exceeds shuffled-contact and random-position controls, but those controls
  remain close enough that contact-specific biological structure is not
  established.
- The GNK non-procedural comparator remains stronger topologically and much
  more high-order dominated: held-out SE `0.6218`, RSE `0.2710`, local maxima
  `5.02`, Walsh high-order fraction `0.3511`, and residual high-order to base
  order-1 variance ratio `3.49`.

## Interpretation

`RUN-067` strengthens the procedural-engineering part of the record and
clarifies the science. It fixes the `RUN-065` sampling defect, keeps the data
volume small, and records the relative-scale quantities needed to avoid the
incorrect shorthand that diminishing returns automatically implies SE/RSE.

The scientific conclusion is mixed. The data are consistent with the revised
causal chain in one limited sense: when residual high-order variance becomes
large relative to the smooth order-1 background, WHT high-order fraction,
RSE, and local maxima increase. However, the stronger adaptive-walk-stage
claim is not supported. The fixed walk-initial residual setting does not
produce a clean stage-wise increase in SE/RSE as available beneficial effects
decline, and calibration prefers the local-DFE-scaled residual.

Compared with `RUN-065`, `RUN-067` is a better-designed test but a less
expansive positive result. The local topology signal persists across all start
classes, but the selected contact residual is weaker than the `RUN-065`
primary on local maxima and similar on SE/RSE. More importantly, matched
controls still challenge biological specificity, and GNK remains a stronger
non-procedural topology generator.

## Result records created

- `RES-026`

## Hypothesis updates

- `HYP-005` is supported narrowly for the revised instrumentation and residual
  exposure diagnostics: residual scale and residual high-order exposure move
  spectral/topology metrics directionally, and the graphlet sampler now behaves
  as intended under a cap.
- `HYP-004` remains weakened as a broad biological ruggedness claim. The run
  shows diminishing-return proxies and residual-induced local topology, but it
  does not validate a robust peak-proximal SE/RSE mechanism or a
  contact-specific procedural model of biological ruggedness.
