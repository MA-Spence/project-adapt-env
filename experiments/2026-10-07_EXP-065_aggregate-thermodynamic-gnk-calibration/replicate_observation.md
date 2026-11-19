# EXP-065: replicate-aware empirical target and frozen observation comparison

The scientist authorized this as one new RUN after ANA-064. The question is:
**How much of the model–data discrepancy is reproducible biological variation,
and how much is compatible with the measurement process?** This is a diagnostic
reanalysis of RUN-166, not independent empirical confirmation or a new model fit.
No scientific success thresholds, empirical interaction reconstruction, new
latent traits, or changes to the public generator are introduced.

## Inputs and provenance

Retrieve the original Johnston et al. TrpB summary table, simplified table, and
relevant processing code from CaltechDATA DOI 10.22002/h5rah-5z170 using verified
ZIP-member byte ranges. Preserve exact members with CRC and SHA-256 provenance.
Do not execute the external notebooks. `prepare_trpb_replicates_exp065.py` produces
an explicit replicate table and freezes copies of RUN-166's empirical cubes,
observable cubes, selected parameters and manifest. All 768 prior sampled
measurements match the original replicate average within 1e-7 (import rounding).
The full public summary has 159,129 non-stop variants and 34,041 stop-containing
controls. Preserve negative measured fitness; do not impute missing genotypes.

The processing code establishes that both replicates use **averaged input
counts**, followed by per-replicate/timepoint stop-background subtraction and
maximum normalization. The final summary averages times 16, 20, 24, 28 and 36 h.
Consequently cross-replicate covariance may contain shared error. Stop variants
are control proxies, not certified identical biological phenotypes: readthrough,
other genotype effects and shared normalizations may also affect them.

## 1. Audit and stratify the original empirical target

Keep exactly RUN-166's 24 calibration and 24 confirmation cubes, with its original
calibration offset/range. Recover both published replicate values at every node.
Export individual measurements, input depth, stop count, published activity flag,
node identities and all cube-level values. Retain the rounded imported mean as an
additional provenance comparator; the primary replicate analysis uses the
unrounded original average. Assess source agreement, published activity cutoffs,
replicate differences, and covariance of control measurements.

Keep the full panel and natural cube strata: no active nodes, some active nodes,
and all 16 active nodes, using the published two-replicate detection definition.
Retain unavailable strata as zero counts; do not select additional favourable
cubes to fill them. Stratum means are descriptive and conditioned on noisy
activity classification. For synthetic noisy replicates use the identical assay
cutoffs on the imported raw scale. No claim that these cubes represent the
whole protein or a random set of enzyme sites is made.

## 2. Cross-replicate WHT and fixed-mutation slopes

For WHT coefficients c1 and c2, the raw cross-power at order k is the sum of
c1(S)*c2(S) over |S|=k. This estimates reproducible variation **plus correlated
error**, under compatible replicate scales. Under genotype-independent errors
with known per-node cross covariance q_i, subtract
`choose(4,k) * sum(q_i) / 256` from that order's cross-power. Here q_i is estimated
from training stop controls in the same input-depth quartile. Report this as a
**stop-control-adjusted sensitivity**, beside raw cross-power and ordinary WHT of
the average. Never force negative estimates to zero. Fractions are signed
estimates and are available only when their total denominator exceeds
max(resolution^2, 1e-18); they can lie outside [0,1]. Report total power and
availability with every fraction; also report ratios of aggregate powers.
High-order fractions use orders 3–4, as in RUN-166; epistatic fractions use orders
2–4. Preserve absolute powers for every order.

For each fixed 0-to-1 mutation across its eight backgrounds, calculate ordinary
slopes and the symmetric cross-replicate numerator
`(Cov(x1, effect2) + Cov(x2, effect1))/2`, divided by `Cov(x1,x2)`.
These centered covariances use divisor n=8. The same shared-error sensitivity
adds `(1-1/n)*mean(q_low)` to the numerator and subtracts it from the denominator.
This formula assumes errors at different genotypes are independent; global
normalization or correlated sequence errors can violate it. Export numerator,
denominator, correction, signed mean effect, and eligibility for every mutation.
Primary summaries include all mutations with a positive resolved denominator;
beneficial strata additionally require mean effect > resolution, using the
average of the two replicates. This classification itself is noisy and is not
claimed to provide an unbiased beneficial-only estimator. Undefined ratios remain
unavailable; no shrinkage, clipping or replacement with zero is used.

Primary resolution 0.001 and sensitivities 1e-12/0.01 reuse RUN-166 units. These
are descriptive resolutions, not estimates of assay noise. Bootstrap the 24
confirmation cubes for 2,000 draws, retaining within-cube pairing. Aggregate
cross-power fractions are formed after resampling sums. Parameters/control
estimates remain fixed; intervals omit uncertainty in the noise model and
shared-site/allele dependence.

## 3. Freeze landscapes; estimate the observation process from other measurements

Use all eight already selected RUN-166 readout/ratio arms and their exact saved
64 confirmation cubes each. No optimization or latent regeneration occurs.
Compare intrinsic outputs with the empirical average and with replicate-aware
summaries, and apply the following two observation processes with 128 Monte Carlo
draws. The generator's affine parameters remain frozen.

Exclude every genotype used in RUN-164, RUN-166 and prior engineering probes from
non-stop noise donors. Split remaining non-stop genotypes and stop controls 80/20
using a fixed RNG; only training data calibrate noise. Fitness strata use zero,
the average published activity cutoff, and the median training active fitness;
input-depth bins are training-stop count quartiles. Observed donor means act as
proxies for latent fitness when assigning noise strata; this approximation can
introduce conditioning bias and is not an identified measurement model. Each synthetic node receives
a depth sampled from training non-stop donors in its latent-fitness stratum.
Sparse donor cells (<30) fall back to the same fitness stratum, then all training
donors, with explicit counts/fallback labels. This is an estimation safeguard,
not a scientific acceptance criterion.

- **Independent-difference arm:** center the training differences f1-f2 within
  fitness/depth cells. Independently resample and randomly reflect these residuals
  for each replicate, dividing by sqrt(2). This reproduces the measured differential
  variance under an equal-error-variance assumption and assigns no shared error.
- **Shared-stop sensitivity:** add a common error to both synthetic replicates.
  Its variance is max(0, training-stop cross covariance) in the assigned depth
  bin. Resample the centered stop-control replicate means and scale them to that
  variance. Retain the independent component above. This makes the biological
  assumption about control covariance explicit, preserves a non-Gaussian residual
  shape, and avoids treating replicate disagreement as total measurement error.
  Report signed estimated covariance even if simulation requires a nonnegative
  common component. The model assumes centered, exchangeable replicate errors,
  independent errors across genotypes and transferable control covariance.

The same random uniforms pair observation draws across readouts/ratios, and the
shared arm adds to the same independent errors. No simulated values are clipped.
Holdout non-stop donors/control variants assess differential error variance and
stop covariance calibration descriptively, without tuning to improve agreement.
Export all per-node synthetic replicate observations (compressed arrays), and all
per-cube ordinary/cross/adjusted summaries. Compressed arrays include every
fixed-mutation numerator, denominator, correction, signed mean effect and
resolution-specific eligibility, with explicit row-index tables. Intrinsic reference and noise-free
readouts remain separate from measured replicas. Report seed-level summaries;
128 observation draws are not 128 independent landscapes. Do not refit the
RUN-166 joint distance to a new empirical target.

## Outputs and interpretation

Save source/input hashes, config, package versions, exclusions, donor split and
noise calibration audit; empirical replicate cubes and fixed-mutation summaries;
control diagnostics; cross-power and slope tables; simulation arrays and summaries;
conditional bootstrap summaries; and figures (PNG/PDF) showing:

1. replicate/control variation and shared covariance;
2. empirical ordinary/cross/stop-adjusted WHT and slopes, with availability;
3. frozen intrinsic/observed synthetic comparisons and activity strata.

A reduction in discrepancy after observation modelling supports a role for the
assay process under its stated assumptions. A persistent reproducible discrepancy
supports revisiting the latent spectrum or site-sampling structure. Neither
outcome proves a unique physical explanation. Retain disagreements among noise
assumptions and empty strata, rather than forcing a conclusion.

## Sources and execution

- Johnston et al. (2024), https://doi.org/10.1073/pnas.2400439121 ; Methods, stop
  controls and replicate-aware original processing motivate the measurement test.
- Data/code: https://data.caltech.edu/records/h5rah-5z170 .
- Reddy & Desai (2021), https://doi.org/10.7554/eLife.64740 ; aggregate diminishing
  returns does not uniquely identify a nonlinear physical readout.
- Sailer & Harms (2017), https://doi.org/10.1534/genetics.116.195214 ; phenotype
  scale and high-order epistasis require separate interpretation.

`run.sh` defaults to `replicate-observation`; `ADAPTENV_PROTOCOL=paired-readout`
and `initial` retain previous entry points. The shared cluster input bundle is
hash-verified before submission. Quick execution uses separate observation seeds
and reduced draws solely to test the pipeline; it is labelled SMOKE_ONLY. CPU
allocation remains 2 CPUs, 8 GB, 2 hours. No ANA/RES is created by submission.
