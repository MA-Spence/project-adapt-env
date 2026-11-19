# EXP-033: PHOT soft-min Boltzmann coupled landscape ABC

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-008
- Motivating results: RES-030, RES-031, RES-032

## Question

For PHOT, does a frozen RUN-077 base landscape combined with a correlated
auxiliary landscape through a temperature-controlled soft-min/Boltzmann
state-competition operator improve held-out graphlet topology, spectra, and
locality relative to the fitted-current, smooth-core, empirical-readout,
GNK/RMF-style, and previous hard-min controls?

## Pre-experiment prediction

If the hard-min failure in RUN-081 was mainly caused by a too-sharp
state-boundary operator, SMC-ABC should concentrate posterior mass on a finite
temperature `tau` and a correlated auxiliary landscape that improves
calibration graphlets without losing base correlation. The same posterior
samples should then improve held-out graphlet distance and preserve empirical
accessibility/local-maxima/random-walk structure. If posterior mass collapses to
near-hard-min limits, diffuse innovation, weak base correlation, or
assay-objective gains without held-out graphlet gains, this implementation does
not support HYP-008.

## Rationale

RES-030 motivates HYP-008 because system-specific current refitting improves
PHOT graphlet dynamics but leaves a high-order spectral mismatch and weak
separation from smooth core. RES-031 shows that coupled machinery can move
PHOT graphlet distance, but the tested empirical local modules fail the
locality control. RES-032 shows that a strict hard-min weakest-link operator
can improve the PHOT full-assay objective but does not beat fitted-current or
smooth-core controls on held-out graphlets and does not validate local
non-diffuse ruggedness.

The next test is therefore not a larger unconstrained fit. It is a guarded ABC
test of a small thermodynamic state-competition operator:

```text
g_tau(a,b) = -tau log(exp(-a/tau) + exp(-b/tau)).
```

As `tau -> 0`, this recovers the hard-min limit. Finite `tau` tests whether
soft state competition can generate the missing graphlet package without using
held-out graphlets or reopening the base smooth-core fit.

## Experimental design

- Use `PHOT_CHLRE_Chen_2023` as the only empirical system.
- Restore RUN-077/RES-030 PHOT outputs from DVC and rebuild the fitted-current,
  posterior-mean, and fitted smooth-core landscapes.
- Define the primary base `a(x)` as the fitted smooth-core PHOT landscape.
  Keep the fitted-current model as a control, not a freely refitted component.
- Construct a correlated auxiliary landscape `b(x)` from standardized `a(x)`
  plus deterministic additive and sparse pairwise auxiliary terms.
- Fit only the auxiliary/operator parameters with SMC-ABC:
  auxiliary correlation, auxiliary scale, auxiliary offset, auxiliary pairwise
  scale, discrete auxiliary pairwise density choice, and `log_tau`.
- Fit the ABC objective only on calibration graphlets, with penalties for low
  base correlation, low edge-gradient correlation, excessive innovation scale,
  and unused state competition.
- Include a small PHOT assay guardrail term in the ABC distance, but do not
  treat full-assay objective improvement as sufficient support for HYP-008.
- Evaluate posterior best, posterior mean, and posterior predictive particles
  on held-out graphlets.
- Report full PHOT assay fit, single-mutant DFE, double-mutant fit, and
  epistasis-prediction metrics for baselines and posterior-selected candidates.
- Compare to empirical observed graphlets, empirical global-readout residual,
  fitted smooth core, fitted current, posterior mean, fitted-current-minus-
  smooth innovation, and a GNK local null.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- RUN-077 PHOT fit summary, prepared PHOT assay, MAVE-NN metrics, and PHOT
  alignment.
- Current root `src/` and `scripts/` code.

## Configuration

- See `config.yaml`.

## Execution

- Use labproj submit for RUN generation.
- Use a dry run before live Slurm submission.

## Expected outputs

- summary.json
- preflight.json
- progress.json
- parameter_specs.csv
- posterior_particles.csv
- posterior_rounds.csv
- evaluated_parameter_sets.csv
- assay_fit_metrics.csv
- assay_dfe_metrics.csv
- assay_epistasis_metrics.csv
- assay_fit_rankings.csv
- graphlet_catalog.csv
- graphlet_metrics.csv
- walsh_terms.csv
- gft_coefficients.csv
- graphlet_summary.csv
- model_distances.csv
- top_model_distances.csv
- local_signal_by_graphlet.csv
- local_signal_energy_by_feature.csv
- local_signal_feature_energy_summary.csv
- local_signal_graphlet_energy_summary.csv
- locality_covariance_summary.csv
- empirical_readout_fit.json

## Analysis plan

- First verify procedural validity: source RUN-077 outputs restored, graphlet
  calibration/test split intact, SMC posterior rounds completed, and posterior
  mass is not unconstrained relative to the prior.
- Evaluate whether posterior-selected soft-min/Boltzmann candidates beat
  fitted-current and smooth-core controls on held-out graphlet distance.
- Evaluate whether any improvement is a coherent graphlet package: WHT/GFT
  spectra, SE/RSE, local maxima, accessibility, and random-walk autocorrelation
  must all be considered.
- Evaluate locality through edge and plaquette innovation energy
  concentration, graphlet energy concentration, soft state-transition
  fractions, and shared-feature covariance.
- Treat assay fit as a guardrail. Assay improvement without held-out graphlet
  and locality improvement weakens rather than supports the operator.

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- Local quick smoke completed with `LABPROJ_QUICK=1`. It produced the expected
  compact output set with four ABC particles, two ABC rounds, two quick PHOT
  graphlets, and the assay/graphlet/locality tables. This is execution evidence
  only, not scientific evidence for HYP-008.
- The dry-run submission is `RUN-082`.
- The first live submission record is `RUN-083`; it failed before Slurm job
  creation because the sandboxed submit environment could not resolve
  `raspi1.maspence.com`. No scheduler job was created.
- The successful live submission is `RUN-084`, submitted to `lab-slurm` as
  Slurm job `239` on `ripley` with `4` CPUs, `32` GB memory, and `08:00:00`
  walltime.
- `RUN-084` completed successfully with scheduler exit code `0:0`. The
  collected outputs are DVC-tracked and backed up under
  `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084`.
- `ANA-033` records the analysis, and `RES-033` records the result: soft-min
  Boltzmann ABC improves calibration graphlet distance and slightly improves an
  assay objective, but the calibration-selected posterior candidates do not
  beat fitted-current or smooth-core held-out graphlet controls and do not
  validate empirical-like locality.
