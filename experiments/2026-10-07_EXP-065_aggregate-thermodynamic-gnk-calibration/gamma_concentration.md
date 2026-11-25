# TrpB Gamma concentration calibration with fixed observation noise

## Question

Under EXP-065 (AIM-001, HYP-004/HYP-005), does allowing unequal interaction
magnitudes improve TrpB aggregate recovery over the Gaussian fit? Authorized
25 November 2026 following RUN-188/ANA-073. The scientist accepts the observation
layer as an experimental modeling assumption. The run does not establish that
stop-control covariance identifies measurement error in functional variants.

## Only new model coordinate

Use the existing frozen Adapt-Env full-pairwise backend: each uncentered table
entry is sqrt(V) Z, with V ~ Gamma(shape=alpha, rate=alpha) and Z ~ Normal(0,1).
Thus E[V]=1 and the uncentered marginal excess kurtosis is 3/alpha. Small alpha
produces many weak and a few strong entries; entries remain signed. Double-center
the full 20×20 pair tables, then apply the existing exact mutation-RMS
normalization. The centered, normalized tables are not independent Gamma draws.
Profile additive preferences, six pair supports, reference selection, sigmoid
architecture and definition of r remain unchanged. Gaussian sampling is the
alpha→infinity limit, not exactly alpha=100. Compare to the exact Gaussian control.

Fit alpha on [0.05,100], logarithmically sampled, reusing the previously declared
Gamma bounds. Jointly fit r continuously on [0,8] and the four original readout
coordinates with their RUN-186 bounds. Scale/steepness coordinates remain logs.
No grid search or fitted noise strength. Gamma variances use the existing seeded
standard_gamma sampler; fixed streams ensure reproducibility but not a smooth
sample realization as alpha changes. Full tables are regenerated only for the
64 calibration cubes per trial; the 160,000-genotype library is evaluated only
for confirmation of the selected model. No exponential table over whole sequences.

## Matched fitting and historical controls

Use exactly the RUN-180 direct fitness objective imported from RUN-186's frozen
source: equally weighted normalized 101-quantile error, local-cube SD error and
near-floor-fraction error. No architecture, sign, cross-replicate or graph metrics
enter fitting. RUN-184's indirect-inference objective is not used.

One serial Optuna 4.9.0 TPE study: 11,712 trials, seed 58459495, 20 startup trials.
This matches RUN-188's per-arm budget. Nine queued candidates consume that budget:
the six historical r/readout settings reevaluated at alpha=1; the Gaussian winner
reevaluated at alpha=0.05 and 100; and the RUN-188 refitted r0 winner at alpha=1.
There is no exact Gaussian candidate in the finite-alpha search. The included r0
candidate bounds the attained loss by the attained null loss; the Gaussian may
still fit better. Six-dimensional proposals differ from the five-dimensional
Gaussian study despite the same sampler seed. This tests achieved fixed-budget
performance, not global optimality or a uniquely identified alpha.

Reuse RUN-188's Gaussian and r0 selected parameters without refitting either
control. Authenticate its output and runtime manifests and the nested RUN-186
archive. Reproduce all historical Gaussian objectives and latent fields before
fitting. During confirmation, reproduce every archived r0/full intrinsic value,
replicate, depth and covariance on its original cubes. Save all three arms for
paired analysis and distinguish reused controls in selected_parameters.json.

## Frozen sampling and observations

Same 32 empirical calibration cubes, 64 synthetic calibration cubes, genotype-
disjoint 32 empirical confirmation cubes, and 24 confirmation seeds with eight
assigned cubes each (16 reused and eight originally fresh). These already inspected
panels are diagnostic confirmation, not newly independent empirical validation.
Keep the shared-stop donor law, input-depth pools, observation uniforms, raw-scale
origin/scale/floor, activity cutoffs and three resolutions exactly unchanged.
No donors or observation amplitudes are refitted. Changing latent predictions may
change the conditional noise bin, so common uniforms do not mean identical errors.

## Outputs and interpretation

Save the complete SQLite study, all trial coordinates/loss components, selected
parameters, selected calibration and whole-library latent fields, cube-level
intrinsic and replicate values. Repeat RUN-188's absolute Walsh powers,
cross-replicate and shared-stop-adjusted powers, directed repeatable sign changes,
activity strata and whole-library static graph summaries on the same coverage
mask. Also save original RUN-186-style cube metrics (including diminishing-return
slopes), marginal quantiles, mean local SD and floor mass for all three arms.
These are descriptive diagnostics, not additional optimization criteria. No
learning or adaptive-walk task is run. Primary reporting should distinguish
intrinsic and observed layers and the eight originally fresh seeds from all 24.

A lower calibration loss alone does not demonstrate improved biological structure.
Assess whether confirmation improves consistently over the Gaussian control, how
large improvements are relative to between-seed spread, and whether alpha hits a
bound. No preset promotion, new adequacy threshold, or claim is part of this run.

## Execution and reproducibility

A source/input/environment manifest freezes the runner, configuration, tests and
complete RUN-188 archive (including RUN-186). The dedicated existing Linux Python
3.11.15 environment must match all pins in environment/readout-py311.lock.txt.
The new explicit gamma-concentration-slurm dispatcher leaves other defaults alone.
Submit via labproj preview then execute, using the recorded Dutch SSH transport,
one CPU, 8 GB and 24 hours. A 24-trial/two-seed smoke is SMOKE_ONLY and cannot count
as the production experiment. Failed/noncomplete trials abort; studies cannot be
overwritten or silently resumed. Archive and checksum the runtime with DVC before
submission. No package source changes are needed.
