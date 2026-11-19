# EXP-065: Aggregate thermodynamic GNK calibration and confirmation

## Follow-up after RUN-164 / ANA-063

On 2026-10-07 the scientist authorized one additional production RUN under this
EXP, combining a paired linear/sigmoid comparison, conditional statistics and
embedded four-site libraries in longer landscapes. The prospective design is in
[paired_readout.md](paired_readout.md), with settings in
[paired_readout.yaml](paired_readout.yaml). It uses fresh TrpB cubes and seeds,
excludes the earlier RUN's genotypes, and imposes no numerical success gates.
The launcher now defaults to this follow-up; `ADAPTENV_PROTOCOL=initial`
reproduces the original protocol below. RUN-164's submitted source is also
preserved in its collected bundle.

## Original RUN-164 protocol

**Question:** Can the normalized rank-8 thermodynamic generator reproduce the
joint aggregate statistics of measured complete genotype subspaces?

This is one small calibration-and-confirmation experiment. It fits distributional
parameters, never empirical mutation coefficients or interaction identities.
It supersedes the proposed EXP-065 prerequisite in EXP-064 for this newly
authorized objective. Historical runs and their conclusions are unchanged.
There is no structured-versus-random routing requirement, universal-assay setting,
or inherited numerical pass gate.

## Design

Use the existing GB1, TrpB, and CR9114-H1 data. Sample 48 complete four-dimensional
binary cubes per assay from uniformly drawn measured backgrounds, variable sites,
and alternative alleles. Reject incomplete cubes and genotype overlap. The first
24 cubes calibrate and the remaining 24 confirm. Save all sampled nodes and raw
values. Exclude rows with nonfinite fitness, record their count, and average finite
duplicate measurements of a genotype before sampling. Missing genotype coordinates
are an error. No imputation, activity filtering, or preferred residue routing is allowed.

Apply one affine transformation per assay, learned from calibration values only:
subtract their minimum and divide by their range. Do not clip confirmation values.
A constant calibration panel uses a declared unit scale and is reported as flat.
The primary comparison resolution is 0.001 of that range, with sensitivity checks
at 1e-12 and 0.01. These are engineering resolutions, not measured noise floors.
Plateaus may represent biological neutrality or assay censoring; these data alone
cannot distinguish them.

Compare three families, each independently selected on calibration summaries:

- Additive stability-like trait through a sigmoid.
- Bare normalized rank-8 GNK field with an affine output scale and offset.
- Additive trait plus normalized rank-8 GNK field through a sigmoid.

The composite varies steepness, reference occupancy, interaction/additive ratio,
and two order mixtures. All fields use the same explicit seed panels; seeds are
averaged, never selected for resemblance to an assay. Four seeds calibrate and
eight new seeds confirm. Rank and support degree are fixed. See `config.yaml`
for the complete finite search, including independent baseline parameter searches.

The synthetic comparison represents the **four active sites**, with q20 alleles
and random binary restrictions. It does not pretend to reconstruct the remaining
protein. Random empirical and synthetic anchors define a sampled-background DFE;
they are not labelled wild type. Distinct assay classes may select distinct
parameter settings. No measured sequence identity enters the generator.

## Measurements and readout

Use the same aggregate_v2 functions for model and data: fitness SD, WHT variance
at orders >=2 and >=3, same-mutation/background diminishing-returns slope for
mutations whose mean effect across backgrounds is beneficial,
anchor-conditioned beneficial fraction, neutral-edge fraction, sign and reciprocal
sign epistasis, strict-peak density, and strictly increasing opposite-corner path
fraction. Output every WHT order and absolute variance as well. Flat/unresolved
spectral fractions stay unavailable. Count neutral components separately from
strict peaks. Also report all-mutation and deleterious-mutation slopes separately;
averaging opposite directions can cancel a concave-readout effect. The uncorrected
diminishing-returns slope is descriptive because
these inputs lack independent replicate errors; it is not evidence for a unique
mechanism. Single-mutant-only data are not used to estimate it.

Select parameters by energy distance between the **joint descriptor distributions**,
standardized using empirical calibration cubes only. The numerical scale floor is
0.05 in the declared normalized descriptor units. Undefined features use an
explicit availability indicator in the internal distance embedding; reported
measurements remain null, never invented zeros. Repeat comparison at the stated
resolutions with selected parameters fixed. Report held-out joint error,
descriptor errors, and variability across independent generator seeds. Do not
collapse the result into a pass/fail label.

**Prediction:** the composite will reduce held-out joint error relative to both
baselines, while generating measurable epistasis and diminishing returns in the
same realizations. If it does not, report which summaries and baselines explain
the failure. Equal performance supports choosing the simpler family. An accurate
aggregate match does not identify physical interaction terms.

## Scalability check

For lengths 32, 128, and 512 and three fresh seeds, record build time, batch-query
time, support count, and additive/GNK/observable single-mutation RMS. The mutation
operator is uniform q20 backgrounds and uniform nonself replacements. The latent
RMS should remain near one; the observable can saturate and that is reported.
This tests finite-length scale behaviour and query cost, not exact whole-protein
WHT, known global optima, or length-invariant random-subcube topology.

## Outputs and scope

`manifest.json` records effective configuration, source/input hashes, versions,
and generator fingerprints. `empirical_cubes.csv` preserves raw measurements and
the disjoint sampling audit. `calibration.csv`, `selected_parameters.json`,
`confirmation.csv`, `resolution_sensitivity.csv`, `scaling.csv`, and `summary.json`
provide the complete comparison. Quick and synthetic-smoke runs are labelled
`SMOKE_ONLY`; they use reserved engineering seeds and cannot support the prediction.
Production confirmation seeds were replaced before registration was finalized
because early engineering checks used a subset of the original panel. No production
confirmation result has been inspected.

Confirmation supports conditional recovery within these assays and sampled
subspaces. It is not held-out-protein generalization or a universal biological
prior. On 2026-10-07 the scientist authorized queue submission. The production
job requests the `lab-slurm` CPU partition with 2 CPUs, 8 GB RAM, and 2 hours.
The launcher accepts an explicit data root or uses the existing shared cluster
bundle; the effective path and every input hash are recorded in the manifest.
CSV and Parquet inputs use the same sampler; Parquet is used when CSV is absent.

## Basis

- RES-058/060/061: normalization, operating scale, and rank-8 component evidence.
- RES-062: remaining empirical shortfalls and the plateau-confounded old target.
- [Otwinowski et al. 2018](https://doi.org/10.1073/pnas.1804015115): nonlinear latent-trait maps.
- [Sailer & Harms 2017](https://doi.org/10.1534/genetics.116.195214): scale dependence of WHT epistasis.
- [Reddy & Desai 2021](https://doi.org/10.7554/eLife.64740): aggregate diminishing returns need not identify one mechanism.
- [Brookes et al. 2022](https://pubmed.ncbi.nlm.nih.gov/34937698/): GNK spectral modelling.
- [Tsuboyama et al. 2023](https://doi.org/10.1038/s41586-023-06328-6): an available empirical stability-effect prior; this first experiment explicitly uses the synthetic prior, not a fitted delta-delta-G distribution.


## Replicate-aware follow-up after ANA-064

The scientist authorized one new RUN for the frozen-model measurement diagnostic.
The prospective [protocol](replicate_observation.md) and
[configuration](replicate_observation.yaml) retain RUN-166 landscapes and empirical
cubes, recover original TrpB replicates, and compare differential-noise and
stop-control shared-error assumptions without refitting the generator or imposing
scientific success gates. The default run entry point now selects this protocol.
Earlier RUN source snapshots preserve the previous entry point.
