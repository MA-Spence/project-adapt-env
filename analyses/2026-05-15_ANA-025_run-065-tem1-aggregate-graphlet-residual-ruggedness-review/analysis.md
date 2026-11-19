# ANA-025: RUN-065 TEM-1 aggregate graphlet residual ruggedness review

## Purpose

Review aggregate graphlet residual ruggedness generalization, held-out topology, Walsh/GFT spectra, controls, and sampling limitations.

## Linked experiments/runs

- Experiments: EXP-025
- Runs: RUN-065

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-065_key_metrics.md`
- `tables/run-065_graphlet_sampling_audit.md`
- `tables/run-065_selected_test_model_comparison.md`
- `tables/run-065_primary_calibration_test_transfer.md`
- `tables/run-065_primary_test_by_dimension.md`
- `tables/run-065_validation_tests.md`

## Inputs

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
- `results/RES-024_run-063-dca-clique-haplotype-terms-create-controllable-high-order-walsh-signal-but-limited-topology/result.md`

## Analysis performed

- Reconciled `RUN-065` as the completed live `EXP-025` execution on `lab-slurm`
  with scheduler job `226` and scheduler exit code `0`.
- Tracked the collected durable output directory with DVC at
  `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065.dvc`.
- Confirmed that the run used the reviewed TEM-1 beta-lactamase reference with
  a `768`-sequence MMSeqs/ColabFold alignment profile.
- Recomputed review metrics from the exported CSV/JSON files using
  `notebooks/analysis.py`.
- Compared the calibration-selected primary contact residual against the smooth
  core, RES-024 exact-match haplotype control, DCA pair override, shuffled and
  random-position residual controls, GNK local null, hard-certified positive
  control, and residual-only diagnostic.
- Audited the graphlet sampling frame against the preregistered requirement for
  graphlets from dramatically different adaptive-walk start classes.

## Main observations

- `RUN-065` produced `240` graphlets, `6000` model-graphlet metric rows,
  `720000` Walsh rows, `720000` Hamming-graph GFT rows, and `111000`
  sign-epistasis pair rows.
- The run selected `20` starts and `240` graphlet backgrounds across five start
  classes, with `48` backgrounds per class.
- A structural sampling failure occurred after background generation:
  `graphlet_catalog.csv` contains all `240` graphlets from
  `reference_or_consensus_like` starts. No graphlets were actually evaluated
  from high-fitness close homologues, marginal-stability homologues, distant
  low-fitness homologues, or viable off-consensus starts.
- Within the sampled reference/consensus-like graphlets, the primary
  calibration-selected contact residual (`contact_residual_ensemble__pset-002`)
  generalized from calibration to held-out test graphlets. Calibration versus
  test means were close: SE `0.3202` versus `0.3213`, RSE `0.0956` versus
  `0.1030`, local maxima `3.66` versus `3.54`, accessible-path fraction
  `0.0369` versus `0.0350`, lag-1 random-walk autocorrelation `0.5078`
  versus `0.4976`, and Walsh high-order fraction `0.0688` versus `0.0737`.
- Relative to the smooth core on held-out graphlets, the primary contact
  residual produced nonzero topology: SE increased from `0` to `0.3213`, RSE
  from `0` to `0.1030`, mean local maxima from `1.0` to `3.54`, accessible
  path fraction decreased from `0.3097` to `0.0350`, and lag-1 autocorrelation
  decreased from `0.6455` to `0.4976`.
- The primary residual spectrum remained low-order dominated with a nonzero but
  non-dominant high-order tail: held-out Walsh high-order fraction was
  `0.0737`, while Hamming-graph GFT low/mid/high energy fractions were
  `0.7719`, `0.2197`, and `0.0084`.
- The primary contact residual improved over the RES-024 exact-match control on
  held-out SE/RSE and local maxima: RES-024 exact-match had SE `0.0557`, RSE
  `0.0393`, and mean local maxima `1.92`. However, RES-024 exact-match had
  similarly low lag-1 autocorrelation (`0.4973`) and even lower endpoint
  accessible-path fraction (`0.00354`), so the comparison is not uniformly
  positive.
- DCA/contact placement was not supported as the source of the topology.
  Shuffled-contact and random-position residual controls closely matched the
  primary residual on held-out topology and spectra. Shuffled-contact test SE,
  RSE, and local maxima were `0.2981`, `0.0945`, and `3.40`; random-position
  values were `0.2801`, `0.0954`, and `3.50`.
- The DCA pair override was a stronger topological control than the primary
  contact residual while retaining effectively zero high-order Walsh tail:
  held-out SE `0.4441`, RSE `0.1335`, mean local maxima `6.02`, accessible-path
  fraction `0.00769`, and Walsh high-order fraction approximately `0`.
- The GNK local null remained the strongest topology comparator, with held-out
  SE `0.6271`, RSE `0.2770`, mean local maxima `11.86`, and lag-1
  autocorrelation `0.1998`, but it carried a large high-order Walsh fraction
  (`0.4532`) and is non-procedural by design.
- The residual-only diagnostic showed that topology can be very different
  without the smooth thermodynamic core: `residual_only_diagnostic__pset-000`
  had about `50.6` held-out local maxima but weak SE/RSE and essentially no
  comparable embedded core context. This confirms that embedding in the smooth
  core is a first-order constraint on residual ruggedness.

## Interpretation

`RUN-065` provides a useful but narrower result than `EXP-025` intended. In the
actually sampled reference/consensus-like graphlets, sparse residual fields can
create held-out SE, RSE, multiple local maxima, lower accessibility, shorter
random-walk autocorrelation, and a low-order-dominated Walsh/GFT spectrum with
a nonzero high-order tail. That is a real improvement over the smooth core and
over the exact-match RES-024 mechanism for several topology metrics.

The run does not validate the broader procedural-generalization claim. The
graphlet cap and sampling order meant that graphlets from the non-reference
start classes were never evaluated. Therefore the run does not test whether the
architecture works across dramatically different adaptive-walk regions, even
though those starts and backgrounds were generated.

The stronger scientific problem is structural rather than only procedural:
matched controls largely reproduce the primary contact residual result.
Shuffled-contact and random-position residual fields are close to the primary
candidate, and the DCA pair override is more topologically effective than the
contact residual while requiring no high-order Walsh tail. This means the
evidence supports residual strength/density as a topology knob, but it does not
show that DCA/contact placement contributes meaningful biological structure.

## Result records created

- `RES-025`

## Hypothesis updates

- `HYP-005` is supported only narrowly: residual and pairwise strength/density
  settings can move graphlet topology and spectral metrics directionally in
  the sampled local regime.
- `HYP-004` is weakened as a broad biological ruggedness/generalization claim:
  the run does not cover multiple start classes, and the strongest topology is
  either pairwise-control dominated or non-procedural rather than clearly
  contact-conditioned and biological.
