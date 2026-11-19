# RES-025: RUN-065 aggregate residuals create local topology but fail broad procedural generalization

## Summary

Aggregate residual fields created held-out local graphlet topology with low-order-dominated spectra in sampled reference-like TEM-1 graphlets, but the graphlet cap excluded other start classes and matched controls showed little DCA/contact-specific advantage.

## Generated from

- Analyses: ANA-025

## Relevant hypotheses

- Supports: HYP-005
- Weakens: HYP-004
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-025_run-065-tem1-aggregate-graphlet-residual-ruggedness-review/metrics.json`
- `analyses/2026-05-15_ANA-025_run-065-tem1-aggregate-graphlet-residual-ruggedness-review/tables/run-065_key_metrics.md`
- `analyses/2026-05-15_ANA-025_run-065-tem1-aggregate-graphlet-residual-ruggedness-review/tables/run-065_graphlet_sampling_audit.md`
- `analyses/2026-05-15_ANA-025_run-065-tem1-aggregate-graphlet-residual-ruggedness-review/tables/run-065_selected_test_model_comparison.md`
- `analyses/2026-05-15_ANA-025_run-065-tem1-aggregate-graphlet-residual-ruggedness-review/tables/run-065_primary_calibration_test_transfer.md`
- `analyses/2026-05-15_ANA-025_run-065-tem1-aggregate-graphlet-residual-ruggedness-review/tables/run-065_primary_test_by_dimension.md`
- `analyses/2026-05-15_ANA-025_run-065-tem1-aggregate-graphlet-residual-ruggedness-review/tables/run-065_validation_tests.md`
- `results/RES-025_run-065-aggregate-residuals-create-local-topology-but-fail-broad-procedural-generalization/metrics.json`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/summary.json`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/module_parameter_sweep.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/graphlet_catalog.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/calibration_graphlet_metrics.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/test_graphlet_metrics.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/graphlet_walsh_terms.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/graphlet_gft_coefficients.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/graphlet_sign_epistasis_pairs.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065/graphlet_accessibility.csv`
- `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065.dvc`
- `experiments/2026-05-15_EXP-025_tem1-aggregate-graphlet-residual-ruggedness-generalization/config.yaml`
- `experiments/2026-05-15_EXP-025_tem1-aggregate-graphlet-residual-ruggedness-generalization/runs/RUN-065.yaml`

## Interpretation

- `RUN-065` completed successfully and produced a large dense graphlet review
  dataset: `240` graphlets, `6000` model-graphlet metric rows, `720000` Walsh
  rows, `720000` Hamming-graph GFT rows, and `111000` sign-epistasis pair rows.
- The primary contact residual candidate improved held-out local topology in
  the sampled graphlets relative to the smooth core. On held-out test graphlets,
  SE increased from `0` to `0.3213`, RSE from `0` to `0.1030`, mean local
  maxima from `1.0` to `3.54`, endpoint accessible-path fraction decreased
  from `0.3097` to `0.0350`, and lag-1 random-walk autocorrelation decreased
  from `0.6455` to `0.4976`.
- The primary candidate generalized from calibration to held-out graphlets
  within the sampled region. Calibration versus test means were close: SE
  `0.3202` versus `0.3213`, RSE `0.0956` versus `0.1030`, local maxima `3.66`
  versus `3.54`, accessible-path fraction `0.0369` versus `0.0350`, and Walsh
  high-order fraction `0.0688` versus `0.0737`.
- The primary candidate also had a more plausible spectrum than the GNK local
  null or the RES-024 exact-match stress setting. Held-out Walsh high-order
  fraction was `0.0737`, and Hamming-graph GFT low/mid/high energy fractions
  were `0.7719`, `0.2197`, and `0.0084`, respectively.
- The result is not a validation of broad procedural generalization. Although
  the run generated starts and graphlet backgrounds across five start classes,
  `graphlet_catalog.csv` contains all `240` graphlets from
  `reference_or_consensus_like` starts. The graphlet cap was reached before
  graphlets from high-fitness close homologues, marginal-stability homologues,
  distant low-fitness homologues, or viable off-consensus starts were evaluated.
- DCA/contact placement is not supported as the source of the topology.
  Shuffled-contact and random-position controls closely matched the primary
  candidate on held-out SE, RSE, local maxima, accessibility, and spectral
  fractions.
- The strongest topological comparators were not the primary contact residual.
  The DCA pair override reached held-out SE `0.4441`, RSE `0.1335`, and mean
  local maxima `6.02` with essentially zero high-order Walsh tail. The GNK local
  null reached held-out SE `0.6271`, RSE `0.2770`, and mean local maxima
  `11.86`, but with a high-order Walsh fraction of `0.4532` and without being a
  procedural biological model.
- The residual-only diagnostic confirms that the smooth thermodynamic core is
  a major structural constraint. Without the core, one residual-only setting
  produced about `50.6` held-out local maxima, but embedded weak residual
  settings remained near a single maximum.

## Effect on hypothesis

- `HYP-005` is supported narrowly. Residual and pairwise strength/density
  settings can move local graphlet topology and spectral metrics in the
  expected direction in the sampled local regime.
- `HYP-004` is weakened. The experiment does not validate adaptive-walk
  ruggedness across different starting regions, and the topology that appears
  is not demonstrably contact-specific or uniquely biological.

## Limitations

- The graphlet sampling bug is material: all evaluated graphlets came from
  reference/consensus-like starts. This invalidates the preregistered
  multi-start generalization claim for this run.
- Held-out graphlets are held out by split within one start class, not across
  the intended dramatically different adaptive-walk regions.
- The selected residual controls indicate that strength, density, and pairwise
  frustration explain much of the topology; DCA/contact placement has not yet
  been shown to add biological structure.
- GNK remains a strong non-procedural topology comparator, so novelty as a
  procedural benchmark system is not established by `RUN-065`.
- The GFT implementation for dense binary Hamming graphlets uses the Walsh
  basis as the graph Fourier basis, with Laplacian eigenvalue equal to
  `2 * Walsh order`. This is appropriate for these graphlets but does not test
  arbitrary non-Hamming graph topology.

## Downstream use

- Use `RES-025` as evidence that aggregate residual calibration can create
  local held-out graphlet topology without making high-order Walsh/GFT energy
  dominant in a sampled TEM-1 reference-like region.
- Do not cite `RUN-065` as evidence that the system generalizes across
  arbitrary sequence-space regions, dramatically different adaptive-walk
  starts, or large sequence spaces.
- Do not cite `RUN-065` as evidence that DCA/contact placement is biologically
  necessary for the observed topology; matched controls currently challenge
  that interpretation.
- Treat the main structural lesson as follows: the current residual module is a
  useful topology knob, but not yet a validated contact-conditioned procedural
  model of biological ruggedness.
