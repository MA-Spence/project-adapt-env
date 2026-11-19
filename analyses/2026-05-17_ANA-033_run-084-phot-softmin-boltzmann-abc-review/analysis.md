# ANA-033: RUN-084 PHOT soft-min Boltzmann ABC review

## Purpose

Assess whether `RUN-084`, the PHOT soft-min/Boltzmann two-state ABC run,
supports `HYP-008` by improving held-out graphlet topology, spectra, and
locality relative to frozen `RUN-077` fitted-current, smooth-core, empirical
readout residual, and GNK/RMF-style controls.

## Linked experiments/runs

- Experiments: `EXP-033`
- Runs: `RUN-084`
- Source data:
  `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084`
- DVC pointer:
  `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084.dvc`

## Scratch notes read

- `scratch/coupled_landscape_nonlinear_operator_plan.md`
- `scratch/coupled_landscape_regularisation_metrics.md`
- `scratch/model_signature_derivation.md`
- `scratch/model_literature_falsification_review.md`

These notes constrain the interpretation. A soft state-competition operator can
support `HYP-008` only if calibration-selected posterior candidates improve
held-out graphlet statistics and show empirical-like local, non-diffuse
innovation. Test-oracle improvement or assay-objective improvement alone is not
sufficient.

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-084_preflight_and_output_audit.md`
- `tables/run-084_smc_posterior_summary.md`
- `tables/run-084_graphlet_distance_comparison.md`
- `tables/run-084_selected_graphlet_metrics.md`
- `tables/run-084_assay_fit_comparison.md`
- `tables/run-084_locality_tests.md`
- `tables/run-084_hyp008_decision_table.md`

## Run and data audit

`RUN-084` completed successfully as Slurm job `239` with scheduler exit code
`0:0`. The collected output contains all planned assay, ABC, graphlet, WHT,
GFT, model-distance, and locality tables. The run used `PHOT_CHLRE_Chen_2023`
only, restored the `RUN-077` PHOT system-specific output, kept the primary
`a(x)` anchor fixed to the fitted smooth-core landscape, and fitted only the
small auxiliary/operator ABC parameterization.

The run evaluated `167,529` assay variants and `24` dense PHOT graphlets split
evenly into `12` calibration and `12` held-out test graphlets. The SMC-ABC run
used `32` particles and `3` rounds, with `160` candidate attempts per round.
The collected directory is DVC-tracked and pushed.

## Key observations

- The ABC objective did move in the intended direction on calibration
  graphlets. The calibration-selected best particle has calibration graphlet
  distance `0.6100`, beating fitted current (`0.7318`), smooth core (`0.7454`),
  posterior mean (`0.7839`), and GNK (`0.8411`), while preserving high base
  correlation (`0.9449` phenotype correlation and `0.9505` edge-gradient
  correlation in the ABC objective).
- The central held-out test fails for the calibration-selected particle. Its
  test graphlet distance is `0.4434`, worse than fitted current (`0.4016`) and
  smooth core (`0.4119`). The posterior mean is also worse (`0.4543`).
- One posterior particle, `abc_particle_rank_023`, narrowly beats fitted
  current on held-out test distance (`0.3964` versus `0.4016`), but this is not
  valid selection evidence because its calibration distance is poor (`0.7126`).
  It is a capacity signal only.
- The best assay objective is the ABC posterior mean (`0.4570`), slightly
  better than the fitted-current posterior mean (`0.4631`), smooth core
  (`0.4814`), and fitted current (`0.4856`). That improvement is narrow:
  holdout Spearman is slightly worse than controls, single-mutant Spearman is
  worse than the fitted controls, and epistasis-prediction Spearman remains
  negative (`-0.0881`).
- The empirical held-out PHOT graphlets have mean local maxima `5.17`,
  accessibility `0.1975`, random-walk lag-1 autocorrelation `0.1905`, and WHT
  high-order fraction `0.3788`. The calibration-selected best particle has
  local maxima `4.25`, accessibility `0.0037`, autocorrelation `0.0405`, and
  WHT high-order fraction `0.5268`, so the graphlet package is not rescued.
- Locality is not empirically validated. The calibration-selected best particle
  has very large test plaquette innovation energy (`179,908`) relative to the
  empirical readout residual (`15.4`) and weak shared-feature covariance
  (`0.0934` plaquette-energy similarity Spearman). Concentration alone is
  therefore not evidence of the desired biological local ruggedness.

## Interpretation

`RUN-084` is a valid execution of `EXP-033`, but it does not support the tested
soft-min/Boltzmann implementation of `HYP-008`. The run provides useful
technical information: finite-temperature state competition can improve
calibration graphlet distance while preserving the base landscape, and the
posterior includes at least one held-out-capable particle. However, the
registered evidentiary standard was calibration-selected held-out
generalization plus locality, not test-oracle capacity.

By that standard, the result is negative. The best calibration-selected
particle and posterior mean both underperform the fitted current model and
smooth core on held-out graphlet distance. The assay guardrail improvement is
too small and too isolated to count as support, especially because
epistasis-prediction remains weak and negative. The state-transition/locality
diagnostics show large artificial innovation energies without a convincing
match to empirical local signal.

## Limitations

- PHOT is the only empirical system tested.
- The empirical graphlet target is the processed `latent_phi` layer, not direct
  exhaustive experimental measurement.
- The run freezes the `RUN-077` PHOT base and fits only the auxiliary/operator
  layer, so it weakens this implementation rather than refuting every possible
  nonlinear trait-operator model.
- Test-oracle posterior particles are capacity diagnostics only and should not
  be cited as held-out support.
- The fitted auxiliary parameterization uses one deterministic auxiliary seed
  and a compact ABC parameter set.
