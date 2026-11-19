# Observation-aware global calibration

Authorized 2026-10-08 after ANA-065. This is one follow-up RUN within EXP-065.
The purpose is to find useful empirical-like global regimes for an online ML
landscape environment, with honest residual differences, not a molecular
reconstruction or universal perfection test. There are no scientific success gates.

## Question and contrast

Can the existing fixed latent architecture, calibrated with an explicit
observation process, reproduce withheld aggregate epistasis and diminishing
returns? RUN-168 added noise after fitting amplitude to noisy data, so its
increased spread does not answer this question.

Compare r=4 linear and sigmoid readouts of the same additive/GNK latent panels,
and an r=0 additive sigmoid control. Hold rank 8, reference rank 0.5, orders
2/3/4 and weights 0.35/0.45/0.20 fixed. A global map is used throughout a
landscape; no cube-specific offsets, empirical sequence effects or interaction
identities are fitted. In raw published units, linear is b+c(A+4G); sigmoid is
b+c sigmoid(mu+beta(A+rG)), c>0 and beta>0. Sigmoid fits export directly usable
ThermodynamicConfig dictionaries. These are phenomenological stability-like
coordinates, not identified free energies.

## Prospective sampling and observation assumptions

Use 64 fresh empirical calibration and 64 confirmation four-site binary cubes,
with no shared genotype between either panel or the 1,664 prior exposed nodes.
Draw uniform measured non-stop anchors and uniform alternative residues; accept
only complete, finite, unused cubes. No activity filtering. Select both panels
before estimating observation parameters and exclude their union from donors.
This is new genotype-level confirmation within the previously examined assay,
not independent assay replication.

Retain the ANA-065 observation model: independent symmetrized replicate
residuals, or those residuals plus a shared stop-control proxy. Train on 80% of
eligible donors/controls; reserve the other 20% for descriptive diagnostics.
Fitness/depth conditioning and control covariance are frozen before calibration.
Stop covariance is a sensitivity assumption: biology and shared assay processing
remain confounded. Do not tune either noise model to WHT or slopes.

Use 8 fresh calibration landscape seeds and 16 fresh confirmation seeds, 8 cubes
per seed. All arms share latent panels and random uniforms. Eight observation
draws per calibration cube form a fixed simulation objective; 32 independent
observation draws evaluate each fitted map on both panels. These repeated draws
are observation uncertainty, not additional independent landscapes. The quick
engineering mode only uses already-exposed RUN-166 empirical cubes and separate
synthetic seeds, with two iterations; it supplies no scientific estimates.

## Calibration and withheld predictions

Minimize the sum of (1) MSE at 101 equally spaced marginal fitness quantiles,
divided by empirical calibration variance; (2) squared relative error in average
within-cube population SD; and (3) squared discrepancy in CDF at the average
published activity cutoff. The latter summarizes near-floor mass in the
replicate-average distribution; the joint two-replicate activity labels remain
diagnostic. Targets are replicate-averaged measured fitness, never a supposedly
error-free corrected landscape. Common random draws keep objective comparisons
paired. Differential evolution uses 60 iterations and population multiplier 8,
without local polishing. Bounds, seeds and exact weights are fixed in YAML.
Record all evaluated parameter vectors, objective components, convergence status
and distance to bounds; a finite search does not establish a global optimum.

Withhold absolute WHT powers, equal-cube and pooled order fractions, all-mutation
and beneficial slopes, their numerators/denominators and availability, and the
distributions of cube variance/activity. Report ordinary, cross-replicate and
stop-adjusted estimators with three descriptive resolution settings. Never clip
negative corrected powers or interpret an unavailable ratio as zero. Primary
resolution is the old 0.001 normalized unit translated to raw units.

Evaluate intrinsic maps as well as observed maps, so measurement variation cannot
silently substitute for biological structure. The ANA must regenerate figures
from saved arrays in an executable notebook, audit numerical summaries, show
empirical cube and synthetic seed uncertainty, and explain material residuals.
Neither equal-cube and pooled WHT nor all-mutation and beneficial-only slopes are
interchangeable. Freeze all six fits before interpreting confirmation outcomes.

## Interpretation and literature

Useful resemblance across several withheld summaries supports these parameter
regimes for synthetic ML testing; exact equality is unnecessary. Remaining
systematic differences describe scope. A change of conclusion between observation
assumptions is measurement-limited evidence. No additional latent complexity is
required merely because one statistic differs. Four-site evidence does not imply
unchanged global statistics for every sequence length or sampling policy.

Otwinowski, McCandlish & Plotkin (2018), *Inferring the shape of global epistasis*,
https://doi.org/10.1073/pnas.1804015115 motivates distinguishing readout nonlinearity,
residual structure and measurement variation. Reddy & Desai (2021), *Global
epistasis emerges from a generic model of a complex trait*,
https://doi.org/10.7554/eLife.64740 supports testing aggregate diminishing returns
without identifying microscopic interactions. Neither paper establishes a unique
thermodynamic interpretation of these fits.
