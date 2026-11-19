# EXP-065 follow-up — TrpB alignment-profile additive prior

## Authorization and comparison

On 19 November 2026 the scientist requested restoration of the profile-based
additive component and a matched repeat of the previous TrpB calibration under
this existing EXP. The comparator is the completed TrpB branch of RUN-175,
including its endpoint calibrations inherited from RUN-172. The fully synthetic
control is reused from the scientific record, not rerun. GB1 is outside this
follow-up. No result, claim, model release or scientific acceptance gate is created.

The sole scientific intervention is to replace the synthetic additive matrix
with a log-frequency profile derived from the original Johnston et al. TrpB
family-alignment frequency table. Interaction identities remain synthetic. This tests the consequences
of a biological additive prior for the same calibration and diagnostics; it does
not test PLM performance, coevolutionary interaction recovery or physical energies.

## Matched protocol

Retain the exact RUN-175 TrpB empirical table, 159,129-genotype mask, 32 calibration
and 32 genotype-disjoint confirmation cubes, calibration origin and scale, floor
cutoff, calibration seed list, 16 reused and eight fresh evaluation seeds, eight
synthetic cubes per seed, cube geometries, observation donor assignments, donor
exclusions and fixed uniform variates. Retain r = 0, 0.5, 1, 2, 4, 8; rank 8;
order weights [0.35, 0.45, 0.20]; and the three observation modes none,
independent_difference and shared_stop_sensitivity.

The reference sequence remains the original deterministic seed-specific numerical
anchor, not an asserted natural wild type. This preserves the reference and
normalized interaction field within every paired seed. Only the additive matrix
is replaced. The full-protein alignment is mapped to the assay's four physical
sites; amino-acid identities are not permuted to synthetic labels. A common global
mutation-RMS normalization preserves relative site conservation. No per-column
standardization, conservation-power multiplier, BLOSUM smoothing, empirical pair
couplings or additional fit parameter is introduced. The retained configuration
value reference_rank=0.5 applies only to synthetic matrix construction; a supplied
biological matrix bypasses that allele-rank permutation. The identical reference
sequence can therefore occupy a different percentile under the new prior. This is
part of the declared prior replacement, not an additional fitted adjustment.

Fit only the same four global sigmoid-readout parameters against fitness
quantiles, mean local SD and near-floor mass. Preserve 101 quantiles, original
bounds, differential-evolution seeds, population factor 8, 60 iterations,
tolerance 1e-5, no polishing and one optimizer worker. WHT spectra, diminishing
returns, sign epistasis and all graph outcomes remain unfitted diagnostics.

RUN-175 reused RUN-172's fitted r=0 and r=4 endpoints. Those values cannot be
reused as fitted endpoints for a different prior. Repeat their original RUN-172
calibration using its exact optimizer seeds and objective, now with the profile
additive values. Then repeat RUN-175's intermediate-ratio calibration, using
the newly fitted profile endpoints as the same two additional candidates. This
repeats the original calibration lineage rather than silently freezing readout
parameters trained on the control. Save all endpoint and intermediate traces.

Retain exact graph dynamic programming, both walk policies, upper-1% and best
plateau outcomes, intrinsic latent topology, and relative resolutions
[0, 0.001, 0.01]. Retain pooled and equal-cube WHT summaries, DR means and medians,
and denominator availability. Observations are frozen once per genotype and
paired by the existing uniforms; they are not redrawn during walks. Report all
seeds and strengths, without selecting favourable realizations.

## Profile and provenance

Use the original-study Tm9D8* / TARGET_b0.3 alignment-frequency table,
selected by the paper's published off-the-shelf-prediction methods rather than by
outcomes of this run. The four assay positions are 183, 184, 227 and 228
(one-based full-protein coordinates), with parent residues VFVS. Preserve and
verify the author-provided query, original table and source metadata.

The archived A2M alignment is unusable: independent extraction via two URL paths
confirms that the source member contains only NUL bytes with a matching ZIP CRC.
The intact author frequency table reports all amino acids and gaps, rounded to
0.001. The full-precision mutation table omits each site's wild-type frequency,
which cannot safely be recovered without its missing gap probability. A later
SSMuLA model uses a different alignment and is not substituted. These limitations
are recorded explicitly; no alignment, sequence weights or effective sequence
count is fabricated.

For each selected physical site, remove the gap category and normalize the 20
reported non-gap frequencies to obtain f_i(a). Then use the fixed probability
mixture p_i(a) = 0.99 f_i(a) + 0.01/20 and M_i(a) = log(20 p_i(a)). This is 1%
uniform probability smoothing, NOT a 0.01 pseudocount in sequence-count units.
It prevents zero log probabilities despite the published rounding. The smoothing
choice is fixed before calibration, not tuned to fitness or diagnostic results.
Source weighting and N_eff are unknown; the prior is described as an
alignment-derived frequency profile, not a newly computed exact reweighted MSA.

The prepared input manifest records source URLs and SHA-256 digests, the source
defect, query-to-physical-position mapping, canonical amino-acid order, source
precision, unknown weighting, gap handling, smoothing and exact output matrix.
The package independently supports both raw-alignment profiles and supplied
frequency profiles. The factor 20 cancels under reference centering. This is an
independent-site evolutionary-preference prior, not a fitted thermodynamic free
energy. Changed priors do not inherit CLAIM-004's numerical results automatically.

## Reproducibility and execution

Freeze the original input manifest, profile bundle, runner, reused numerical
kernels, current package source and exact configuration by SHA-256. Worker startup
must verify the manifests before calculations. Preserve the original input/source
records without overwriting them. Validation must check profile mapping and
normalization, unchanged synthetic behaviour, paired normalized interactions,
endpoint optimizer settings, and the independently verified graph kernel.
Any shortened smoke execution is labelled SMOKE_ONLY and is not scientific evidence.
Production startup also requires RUN-175's recorded Python 3.11.15 and numerical
package versions (NumPy 2.4.4, pandas 2.3.0, SciPy 1.17.1, PyArrow 20.0.0,
PyYAML 6.0.3). Local engineering smoke may use a different Python version and
records it; that exception is disabled for the scientific run.

Use the existing lab-slurm CPU resource request: four CPUs, 8 GB RAM and one hour.
On 19 November the scientist confirmed Dutch was down and instructed preparation
only, with queue submission to follow when the worker is restored. A dry-run RUN
record is not an accepted scheduler job. The source and data are prepared locally
for a self-contained submission; no replacement synthetic-control run is requested.

## Scientific basis

Johnston et al. (2024), https://doi.org/10.1073/pnas.2400439121, supplies the assay,
the alignment and its target context. Hopf et al. (2017),
https://doi.org/10.1038/nbt.3769, supplies the established independent-site versus
coupled evolutionary-model distinction. The inherited landscape and observation
limitations in CLAIM-004 and ANA-069 remain unchanged.
