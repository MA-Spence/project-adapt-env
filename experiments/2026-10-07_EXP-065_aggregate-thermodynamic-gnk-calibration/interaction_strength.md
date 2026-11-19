# EXP-065 follow-up — interaction strength, global statistics and adaptive walks

## Question and hypothesis

Can intermediate interaction strength provide useful joint resemblance to GB1 and
TrpB in aggregate epistasis and adaptive search? ANA-068 found longer but overly
concentrated additive-control walks, and shorter walks at r=4. Intermediate r may
balance those behaviours; a persistent trade-off is an informative outcome.

This is one non-production RUN under EXP-065. Public model architecture, rank 8,
reference rank 0.5 and order weights [0.35, 0.45, 0.20] remain fixed. Compare
r = 0, 0.5, 1, 2, 4, 8, pairing additive and interaction fields across ratios.
Use all 16 RUN-172 confirmation seeds and eight fresh seeds at every ratio. Fresh
seeds are designated before outcomes and are reported separately, without choosing
or discarding seeds by fit. Assay identities and interaction coefficients are not fit.

## Frozen empirical information and readout calibration

Reuse the exact RUN-172 calibration/confirmation cubes, full measured GB1/TrpB
parquets, calibration origins/scales and floor cutoffs. Genotype-disjoint empirical
cube splits remain intact. Whole-library graphs include calibration nodes and
noise donors, and these already-inspected assays are exploratory post-fit evidence.

For each ratio and observation assumption, fit only four sigmoid readout parameters
to the existing calibration objective: fitness quantiles, mean cube SD and floor
mass, with the same bounds and differential-evolution budget. The r=0 and r=4
settings retain their original RUN-172 fit parameters as anchors. New ratios use
the same finite search and include the two frozen endpoint parameter vectors as
additional candidate readouts; select the smallest calibration objective only.
All candidates, objective components, convergence and boundary proximity are saved.
Walks, WHT, DR and sign epistasis do not enter parameter fitting.

For TrpB, independently recalibrate the readout under the previously specified
independent-difference and shared-stop observation assumptions. Noise donor pools,
strata and covariance are reconstructed from frozen RUN-172 inputs and held fixed.
Report each intrinsic map alongside its frozen observed map. One observation per
genotype is paired by uniform random draws across ratios/assumptions; no noise is
redrawn during walks, and no noise amplitude is tuned to fit walk statistics.

## Readouts

Enumerate the four-site q20 substitution graph, restricting both empirical and
synthetic maps to the same measured-genotype mask. There is no imputation.
Exact dynamic programming gives uniform-uphill and greedy expected walk lengths,
terminal ranks and probabilities of reaching the top 1% or best measured plateau.
Exact graph reachability permits indirect paths and reversions. Compare epsilon
0, 0.001 and 0.01 times the frozen empirical calibration range. In addition,
evaluate noiseless exact topology on the latent coordinate directly; quantify
any floating-point sigmoid saturation instead of mistaking it for intrinsic topology.

Compute WHT spectra, fixed-direction DR slopes and sign epistasis on the frozen
confirmation cube geometries. Preserve pooled and equal-cube WHT summaries, DR
means and medians, and denominator availability. The observed cubes are subsets
of the same frozen full observed map, maintaining consistency with graph outcomes.
Each old seed uses its original eight cube geometries; fresh seeds cycle over the
same confirmation geometries. Single frozen observation draws are not independent
measurement-uncertainty estimates. Across-seed variability is reported separately
for reused and fresh seeds; it omits fitted-parameter and empirical uncertainty.

Neutral components remain a supporting tolerance/floor sensitivity; they are not
calibration targets. Best-genotype capture is diagnostic rather than an acceptance
condition. Primary figures show the joint trade-off between walk length, top-1%
success and aggregate epistasis, preserving all strengths and all seeds.

## Reproducibility

The input bundle preserves 24 SHA-256-verified RUN-172/ANA-068 files. A frozen
source manifest, worker environment, complete parameter traces, field arrays,
noise seeds, cube fitness values, exact graph summaries and output hashes support
replay. Twenty independent tiny graphs validate the reused C++ graph kernel.
Smoke runs are explicitly labelled and cannot be interpreted as evidence.
No RES, CLAIM, model release, universal realism claim or numerical scientific
success gate is part of this RUN.

## Scientific basis

Wu et al. (2016), https://doi.org/10.7554/eLife.16965, motivates indirect routes,
basins and contrasting uphill walk rules. Nowak & Krug (2015),
https://doi.org/10.1088/1742-5468/2015/06/P06014, motivates examining interaction
structure through walk outcomes. Johnston et al. (2024),
https://doi.org/10.1073/pnas.2400439121, supplies the TrpB empirical landscape.
