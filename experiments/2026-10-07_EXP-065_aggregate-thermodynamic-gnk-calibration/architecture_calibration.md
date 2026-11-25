# TrpB architecture-calibration amendment — 2026-11-21

The scientist authorized a new RUN under EXP-065 and local MacBook execution.
This is the first empirical calibration of candidate 0.7.0rc1 using the new
Optuna indirect-inference method. It does not replace RUN-175, RUN-178 or RUN-180.

## Question and scope

Can a profile-conditioned additive field plus synthetic full pairwise tables and
a sigmoid reproduce the genetic-architecture summaries inferred from measured
four-site TrpB observations, on disjoint genotypes and fresh generator seeds?
This addresses the scoped aggregate-resemblance aim (AIM-001, HYP-004/005), not
recovery of empirical interaction identities, independent biological replication,
comparative optimizer performance, ML performance, or long-sequence validity.

The architecture is anchored to Park, Metzger & Thornton (2024),
[The simplicity of protein sequence–function relationships](https://doi.org/10.1038/s41467-024-51895-5).
The local auxiliary inference follows the pinned authors' SSE plus L1 objective,
mean held-out-R² penalty selection, and post-fit reference-free gauge.
[Implementation mapping](../../external/Adapt-Env/docs/RFA_IMPLEMENTATION.md)
records the different numerical solver. TrpB was not fitted in that paper.
The outer minimum-distance objective and Gamma-normal sampling distribution are
explicit modelling choices; they are not supplied or validated by that paper.

## Frozen inputs and partition

Use only TrpB positions 183, 184, 227 and 228, in alphabet
`ACDEFGHIKLMNPQRSTVWY`. Retain the previous profile additive matrix unchanged.
Eligible genotypes are the 159,129 canonical q20 genotypes in the historical
whole-library mask. Use their two actual finite replicate observations, rather
than treating repeated appearances in historical cubes as independent samples.
A deterministic uniform permutation (seed 20261121), independent of fitness,
selects 8,192 calibration genotypes and the next 8,192 confirmation genotypes.
The two panels are disjoint. The preparation manifest records coverage, raw
sources, checksums, exact masks and rounding differences from historical means.
No imputation, fitness filtering, resampling of favourable genotypes, or selection
of preferred generator seeds is allowed.

The confirmation panel is held out from this fitting procedure. These data have
previously contributed to whole-library project analyses; it is not a previously
unseen biological dataset. The panel's complete-pair allele coverage is checked
without inspecting confirmation outcomes. Replicates and duplicate bootstrap
copies of each genotype remain in the same fixed cross-validation fold.

## Model and objective

At four sites, degree six includes all six site pairs. Each pair uses a full
centered 20×20 matrix. Each uncentered coefficient is sqrt(V) Z, where
V~Gamma(shape=alpha, rate=alpha), Z~Normal(0,1). Row/column centering and exact
mutation-RMS normalization retain the meaning of r. Alpha controls concentration,
not total interaction strength. Profile preferences remain a fixed prior, not
inferred thermodynamic free energies or empirical interaction coefficients.

Optuna serial multivariate TPE searches six controls: r, alpha, sigmoid steepness,
uniform-background mean-logit position, lower output bound and positive output
span. The complete search bounds, seeds, trial budget, inference settings and
bootstrap count are in [architecture_calibration.json](architecture_calibration.json).
One affine standardization is learned from the calibration observations only;
output parameters are also exported in native assay units.

Both measured and simulated panels receive identical auxiliary fitting, using
the same observed genotypes, replicate mask, grouped folds and penalty grid.
The loss averages squared discrepancies in 18 fixed summaries, divided by
calibration-only genotype-bootstrap variances plus the declared squared floor:
additive/pairwise variance and their ratio; within-order standardized absolute
coefficient quantiles; concentration of squared coefficient mass; observed fitness
quantiles; fitted output bounds; and lower/upper sigmoid saturation fractions.
The summaries, normalization and top-10%-mass definition are fixed in
[the method specification](../../docs/calibration/README.md). We do not compare
empirical and simulated coefficient identities.

The primary run uses `noise_mode=none`. Measured replicates retain assay noise;
simulated replicates have the same masks but share noiseless oracle values. This
is a declared intrinsic-model test, not an identified assay-noise model. The
historical shared-stop observation sensitivity is not silently inherited.

Strict numerical convergence remains required for the empirical fit, bootstrap
fits, selected-trial fits and confirmation fits. Failed simulated trials are
retained in the SQLite study with diagnostics and do not receive a fabricated
loss. Failure of an empirical/bootstrap/confirmation fit is a failed RUN, not
successful scientific calibration. A finite outer search is not proof of global
optimality or parameter identifiability. Near-best ranges are descriptive,
not posterior or confidence intervals.

## Confirmation and interpretation

After parameter selection, evaluate the untouched genotype panel with fresh
seeds. It has a separate empirical auxiliary fit, the same inference rule and
calibration-derived loss weights. Do not tune against confirmation.

Retain unfitted sign/reciprocal-sign epistasis, diminishing returns and walk
summaries on the observed panel's induced graph. Missing genotypes are not
filled in. These sparse-panel walks and peaks are not whole-library searchability
estimates, and are not directly comparable with RUN-180's whole-library graph.
Previous runs supply historical context only: the objective, input partition,
coefficient distribution and observation treatment differ. Their loss values are
not commensurate with the new architecture distance. No acceptance threshold,
new result/claim, public preset update or release is part of this execution.

## Execution and reproducibility

A retained preparation pilot times auxiliary inference using calibration data
only and tests strict convergence. It may set solver iteration and execution
budgets before the RUN is frozen; it does not choose model parameters from
confirmation outcomes. Its code, timings and failures are archived as preparation
evidence, with the final choices stated before submission.

The local helper verifies the 0.7.0rc1 candidate lock and a separate RUN source
manifest, snapshots the complete input bundle, and captures the actual Python,
NumPy, SciPy, Optuna, BLAS and host environment. The registered launch executes
from that frozen source/input snapshot with numerical-library thread counts one.
Labproj creates a nonexecuting local preview before the live RUN. A detached
supervisor retains the live logs, registered RUN identity and all final artifacts,
including failure evidence. On macOS it prevents idle system sleep during work.

The archive includes the exact input partition, auxiliary coefficients, fixed
loss weights, every Optuna trial and failure, confirmation predictions, output
parameters, sources, environment and RUN metadata. DVC tracking/backup status is
recorded separately from numerical completion. Prior archived RUNs are unchanged.

### Budget and preparation decisions

The production search retains the example's 100 Optuna trials, 20 startup trials,
four common calibration seeds (101–104), four fresh confirmation seeds
(1001–1004), and 20 genotype-bootstrap weight estimates. It retains 10-fold CV.
The calibration-only pilot showed that the example penalty grid ended at its
best value. The frozen grid therefore retains 0, 0.01, 0.1 and 1 and adds 10 and
100. A 5,000-iteration budget, function tolerance 1e-10 and at most 10 deterministic
curvature resets replace the insufficient example settings. The required
per-observation stationarity tolerance remains 1e-5. These are numerical fitting
controls, not additional synthetic landscape parameters.

With six penalties and ten folds, each complete auxiliary estimation performs
61 optimizer fits. If every trial completes, this study performs 25,986 fits:
24,400 for the outer search, 1,281 for the empirical target and 20 bootstraps,
and 305 for confirmation. This count does not include the separate preparation
pilots or internal solver iterations. The full empirical and synthetic timing
checks, failed attempts, input hashes and chosen configuration are retained in
the preparation archive referenced by the source manifest.

The local request records one CPU and 8 GB, with a nominal 24-hour allocation.
There is no local scheduler-enforced walltime; the finite trial and inner-iteration
budgets bound the computation. Timing pilots are execution estimates, not an
efficiency benchmark.
