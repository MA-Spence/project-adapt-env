# EXP-066: Production landscape aggregate realism across empirical systems

Status: RUN-172 completed (Slurm 2317, attempt 3); analysis pending. Protocol authorized 2026-10-08. This experiment
precedes the proposed online-policy pilot. It extends EXP-065 / ANA-066 under
AIM-001 and the limited aggregate-statistics portions of HYP-004 and HYP-005.

The question is whether a small set of parameters for the public synthetic model
can produce useful global resemblance across available empirical assay classes.
We seek representative regimes for online ML environments. We do not reconstruct
empirical interactions, identify physical energies, demand a perfect model, or
apply previous scientific success gates.

## Public method

Adapt-Env 0.3.0, GitHub `dev`, commit
`dd0bb7a52d7d92af48bbfba115331742b9a0c2eb` exposes
`ThermodynamicLandscape` / `ThermodynamicConfig`; historical public landscape names
are aliases to this method. A mutation-normalized additive latent trait plus
rank-8 synthetic GNK interactions passes through an affine sigmoid:

`fitness = b + c * sigmoid(mu + beta * (A + r G))`.

The rank, order weights (0.35, 0.45, 0.20), degree construction and reference rank
0.5 are fixed. Compare r=4 sigmoid, r=0 sigmoid, and an r=4 linear diagnostic.
The linear readout exists only inside this experiment. Every selected sigmoid
map is checked against public `evaluate_batch` for every evaluated seed. There
are no per-genotype effects, per-cube affine maps, or fitted interaction terms.
The private historical methods are deprecated; the queued snapshot cannot import
them. Public release and runner source hashes are frozen in `source_manifest.json`.

## Existing systems and suitability

Audit all 86 entries in the existing DVC-backed ensemble dataset inventory.
Attempt structural coverage sampling for all 79 nonduplicate, fixed-length
protein assay candidates, including the 64 stability panels. This includes GB1,
TrpB, CR9114 across three antigens, CR6261 across two, ParD2/3, Aakre, Ding, GFP,
HIS7, PHOT and Olson where the measured data permit the statistics. Related
antigen/assay measurements are not independent proteins.

Seven inventory entries have explicit exclusions: two DHFR nucleotide/codon
representations (not amino-acid substitution geometry), three duplicate GB1/TrpB
mirrors, AAV's mixed indel library, and the cross-protein Meltome collection.
All remaining panels are tested for complete measured cubes; do not exclude an
assay for a poor fit. Record inadequate coverage separately from an execution
error. Single-mutant-only data cannot identify the desired WHT interactions.

Read authoritative `pos*` columns where present, preserving leading-zero binary
states. Otherwise compress invariant full-sequence coordinates and record the
source positions, categorical alphabets and reference. Duplicate measured
genotypes are averaged and counted; nonfinite measurements are counted and
excluded. Never fill missing fitness values. All input parquet SHA-256 values
must match the frozen inventory. These are already processed higher-is-better
readouts; no assay-specific nonlinear transformation is fitted to the data.

## Sampling and calibration

For each assay attempt dimensions 4, 3, then 2. Request 32 calibration and 32
confirmation cubes, with no reused genotype. A fixed 20,000 attempts per dimension
and at least eight cubes per split define whether an estimator can be reported;
these are sampling/availability choices, not model success thresholds. Save every
sampling rejection count. Randomize the split after outcome-blind sampling.

Dense libraries use uniformly sampled measured anchors, sites and alternative
alleles. Sparse libraries choose sites at which a measured anchor differs from
the target; each alternative is a reversion with probability 1/2 or a uniformly
sampled alternative measured allele. This is a coverage-aware distribution over
local contexts, not uniform sampling of an unmeasured full protein landscape.
Accept only complete genotype-disjoint cubes, without looking at fitness or
activity. Higher-order power is unavailable, not zero, for dimension two.

Synthetic fields use the same categorical genotype geometries under independent
random seeds, without importing measured residue effects. Use eight calibration
seeds and 16 new confirmation seeds, eight cubes per seed. Geometry indices,
model fingerprints, latent A/G arrays and full genotype arrays are saved. The
model length is the number of variable library coordinates, not full protein
length. This experiment does not test invariance with arbitrary length.

Per assay/observation mode fit only b, c>0, mu, beta>0 (linear controls fit b,c).
One setting covers all cubes and seeds. Normalize numerical coordinates by the
calibration 1st–99th percentile interval (SD fallback); export parameters back in
raw assay units. Constant calibration panels are retained as unavailable for
parameter identification. Differential evolution has 60 iterations, population
multiplier eight, no polish, and the prospective bounds in `config.yaml`.
Retain all objective evaluations, convergence status and distances to bounds.
These finite searches yield candidates, not proven optima.

The loss is the sum of (1) 101-point marginal quantile MSE divided by empirical
variance, (2) squared relative error in average within-cube SD, and (3) squared
near-floor CDF error. For TrpB use its published replicate activity threshold
average as in EXP-065; elsewhere the calibration fifth percentile is a descriptive
near-floor proxy, not a biological inactivity threshold. WHT powers, fractions,
diminishing-returns slopes, DFE and local topology are never optimization targets.

Confirmation is disjoint within this RUN. The model family was developed using
some of these assay classes; this is not independent discovery on previously
unseen biology. Both split labels and all calibration results remain available.

## Observation and reporting

The common comparison fits observed assay values with no added noise. For assays
without a matched observation model, a fit can absorb assay noise and cannot be
interpreted as recovery of intrinsic biological epistasis. For TrpB additionally
repeat the independently estimated replicate-difference and shared-stop covariance
sensitivity branches from EXP-065. Exclude all selected cube genotypes and the
existing donor exclusions before the 80/20 donor/control split. Conditioning on
observed donor means and shared-stop errors remains an assumption. Save intrinsic
and observed predictions separately, plus paired-replicate corrected diagnostics.

For every cube compute exact WHT power by order and both equal-cube and pooled
fractions, fixed-direction effect-versus-background slopes (all and beneficial
mutations), DFE, neutral/beneficial edges, sign and reciprocal sign epistasis,
and strict/weak local peaks. Report denominators and unavailable estimates at
three resolutions: 1e-12, 1e-3, 1e-2 times the calibration fitness range. Local
peaks are not certified global optima. Avoid cross-dimension pooled rankings.

`figure_data.csv` retains one row per metric and empirical panel or synthetic
seed. Raw cubes and per-mutation numerators/denominators permit reanalysis.
A figure should show every evaluated assay, including mismatches, and separate
measurement sensitivities from intrinsic predictions. Resample empirical cubes
and synthetic seeds for uncertainty; repeated observation draws are not
independent landscapes. Suitable panels are measured versus synthetic fitness
CDFs, WHT order fractions, diminishing-returns distributions, and local topology.
No aggregate score is used to declare biological equivalence.

## Execution and durable data

Run with `labproj submit EXP-066 --target lab-slurm --execute` after a dry run.
CPU partition, four allocated CPUs, 16 GB, two-hour walltime. Fits are sequential
and numerical-library threads are capped at one for deterministic behavior;
the spare allocation permits input processing headroom. The source-only snapshot
uses an explicit staging filter to avoid unrelated cloud-backed historical files.
Inputs are separately verified from the existing shared DVC-backed bundles.

After successful execution `run.sh` verifies the output manifest and atomically
publishes the full bundle to the durable cluster
`data/processed/production_landscapes/<RUN-ID>` counterpart. Normal labproj outputs
are retained. The matching local directory is
`data/processed/production_landscapes/<RUN-ID>`; the collection script verifies
hashes before exposing it, after which it can be DVC tracked. Partial executions
are not silently published as complete. There are no scientific success gates.
No new ANA, RES or CLAIM is created by this submission request.

## Scientific basis

- [Otwinowski, McCandlish & Plotkin (2018), Inferring the shape of global epistasis](https://doi.org/10.1073/pnas.1804015115): a nonlinear mapping from an aggregate latent coordinate is a useful global-epistasis description, not proof of specific physical energies.
- [Reddy & Desai (2021), Global epistasis emerges from a generic model of a complex trait](https://doi.org/10.7554/eLife.64740): aggregate fitness-dependent effects can arise without recovering individual interactions.
- [Johnston et al. (2024)](https://doi.org/10.1073/pnas.2400439121): TrpB combinatorial landscape and the assay-specific observation context used in EXP-065.
- [Phillips et al. (2021)](https://doi.org/10.7554/eLife.71393): antibody landscapes across antigens illustrate why assay contexts should remain separate.
- [Lite et al. (2020)](https://doi.org/10.7554/eLife.60924): three-position toxin–antitoxin specificity landscapes motivate dimension-matched comparisons.

The inventory preserves per-dataset source IDs, publications and DOI provenance.

## Pre-submission validation

The outcome-blind full coverage preflight found 64 usable assay geometries: ten
four-site, two three-site and 52 two-site; 15 candidates lacked enough complete
disjoint cubes under this sampling budget. All 86 input parquet hashes match the
frozen inventory. Five representative worker smoke cases completed without
errors in 38 seconds (TrpB, ParD2, AMFR stability, PHOT and Olson). Five focused
production tests cover WHT/DR mathematics, plateaus, disjoint outcome-blind
sampling, raw-unit public API equivalence and checksum-verified publication.
These checks establish executable methods, not empirical model adequacy.

### Launcher correction before scientific execution

RUN-172 attempt 1 (Slurm 2315) failed in the automatic overlay installer because
the source-only snapshot deliberately omitted editable-package build metadata.
No scientific computation began. The corrected submission disables overlay
installation in its captured project configuration and uses the existing worker
environment validated by the smoke test. `worker_environment.json` pins Python
and the five numerical/data dependencies, and the runner verifies them before
fitting. The normal local project configuration is restored after submission.
This changes execution setup only, not sampling, fitting or scientific criteria.

RUN-172 attempt 2 was stopped after the TrpB compact mirror differed from its
replicate means by 3.17e-8 near zero, beyond the original absolute 1e-8 check.
The complete 64-cube panel was verified within one float32 epsilon. Attempt 3
uses that documented representation tolerance and records the observed maximum
difference. Primary empirical values, objectives and all statistics are unchanged;
partial attempt-2 outputs and logs are preserved separately.
