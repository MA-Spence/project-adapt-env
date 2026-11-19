# ANA-032: RUN-081 PHOT hard-min weakest-link grid review

## Purpose

Assess whether `RUN-081`, the deterministic PHOT hard-min weakest-link grid,
supports `HYP-008` by improving PHOT assay fit and held-out graphlet local
ruggedness relative to frozen `RUN-077` controls.

## Linked experiments/runs

- Experiments: `EXP-032`
- Runs: `RUN-081`
- Source data:
  `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081`
- DVC pointer:
  `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081.dvc`

## Scratch notes read

- `scratch/model_signature_derivation.md`
- `scratch/model_literature_falsification_review.md`
- `scratch/coupled_landscape_nonlinear_operator_plan.md`
- `scratch/coupled_landscape_regularisation_metrics.md`

These notes constrain the interpretation. The hard-min operator can support
`HYP-008` only if it improves held-out graphlet topology/spectra while showing
empirical-like localized innovation. A graphlet-only or assay-only improvement
is not sufficient.

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: `python3`

## Produced artifacts

- `metrics.json`
- `tables/run-081_preflight_and_output_audit.md`
- `tables/run-081_assay_fit_comparison.md`
- `tables/run-081_graphlet_distance_comparison.md`
- `tables/run-081_selected_graphlet_metrics.md`
- `tables/run-081_graphlet_energy_locality_test.md`
- `tables/run-081_feature_plaquette_locality_test.md`
- `tables/run-081_locality_covariance_test.md`
- `tables/run-081_hyp008_decision_table.md`

## Run and data audit

`RUN-081` completed with scheduler exit code `0`. The run used
`PHOT_CHLRE_Chen_2023` only, froze the `RUN-077`/`RES-030` PHOT fit, used the
smooth/additive PHOT base as `A(x)`, evaluated `f(x) = min(A(x), B(x))`, and did
not rerun ABC. The collected output contains the expected assay-fit, graphlet,
WHT, GFT, model-distance, hard-min switch, and local-signal energy tables.

The run evaluated `167,529` assay variants, `24` dense graphlets split evenly
into calibration and test graphlets, and `180` hard-min parameter sets. The
output directory is DVC-tracked and pushed.

## Key observations

- The best hard-min candidate by full-assay objective is
  `hard_min__smooth_core_A__grid0100`. Its assay objective distance is `0.4223`,
  better than smooth core (`0.4814`), fitted current (`0.4856`), and posterior
  mean (`0.4631`).
- That assay improvement is not a full empirical recovery. Holdout Spearman is
  essentially unchanged (`0.7569` hard-min versus `0.7574` fitted current),
  single-mutant Spearman is slightly worse than the controls (`0.4722` versus
  `0.4782`), double-mutant Spearman improves modestly (`0.7389` versus
  `0.7313`), and epistasis-prediction Spearman remains weak (`0.0106`).
- Held-out graphlet evidence fails the main `HYP-008` requirement. The best
  hard-min test distance is `0.4124`, worse than fitted current (`0.4016`) and
  slightly worse than smooth core (`0.4119`). The calibration/assay-selected
  hard-min candidate generalizes poorly, with test distance `0.4631`.
- The empirical PHOT test graphlets remain high-order-rich: observed
  `latent_phi` WHT high-order fraction is `0.3788`. Hard-min does not repair
  this envelope: the best assay candidate has WHT high-order fraction `0.5077`
  and the best test-distance hard-min has `0.4976`.
- Topology and dynamics are not rescued as a package. Empirical test graphlets
  have mean local maxima `5.17`, accessibility `0.1975`, and random-walk lag-1
  autocorrelation `0.1905`; the best assay hard-min has `3.42`, `0.0245`, and
  `0.0464`, respectively.
- The local-signal diagnostics do not validate empirical-like non-diffuse
  ruggedness. Hard-min innovation is concentrated at switch graphlets by
  construction, but raw test plaquette energy is far larger than the empirical
  readout residual (`32,593` for the best assay hard-min versus `15.4` for the
  empirical residual), and shared-feature covariance with innovation/switch
  similarity is weak or negative on test graphlets.

## Interpretation

`RUN-081` is procedurally valid and informative, but it does not support the
hard-min weakest-link implementation of `HYP-008`. The full-assay objective
improvement shows that a small amount of hard-min limiting can move PHOT
ranking and distributional metrics, especially double-mutant ranking and the
aggregate objective. That is a technical signal, not sufficient evidence for a
biological local nonlinear trait operator.

The decisive evidence is held-out graphlet performance. No hard-min candidate
beats the fitted current model on the test graphlet package, and the
calibration-selected hard-min candidate does not generalize. The spectral
mismatch from `RES-030` is not repaired, and local maxima/accessibility dynamics
remain weaker than the empirical graphlets.

The locality evidence is also not the expected support for `HYP-008`. Hard-min
switches necessarily concentrate innovation where the limiter changes, but the
observed concentration does not align with empirical local signal strongly
enough to support a local, non-diffuse ruggedness mechanism. In this run, the
operator can improve the assay objective without solving the held-out graphlet
or locality problem.

## Limitations

- PHOT is the only empirical system tested.
- The run is deterministic and freezes the `RUN-077` PHOT base; it is not a
  joint ABC refit.
- The empirical graphlet target is the processed `latent_phi` layer, not direct
  exhaustive experimental measurement.
- Test-oracle hard-min rankings are useful capacity diagnostics only; they are
  not valid model-selection evidence.
- The hard-min auxiliary grid is compact and uses one auxiliary seed, so this
  result weakens the current implementation rather than refuting every possible
  nonlinear trait-operator model.
