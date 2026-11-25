# EXP-065: Aggregate thermodynamic GNK calibration and confirmation

## 2026-11-24 joint Optuna submission

RUN-186 is running on Dutch as Slurm job **2355**, after preview RUN-185.
The frozen joint protocol optimizes continuous r in [0,8] with four readout
coordinates, using the restored RUN-180 loss and matched TrpB sampling design.
Gamma concentration is disabled. Startup verification confirms five-coordinate
trials with off-grid r values, matching source and inputs, and all 18 historical
objective replays. See [launch evidence](runs/RUN-186.launch_verification.json)
and [engineering validation](joint_readout_validation.json). Results are pending.

## Current protocol: joint continuous-r Optuna fitting, authorized 2026-11-24

The scientist requested a new TrpB RUN with r optimized continuously on [0,8]
alongside all four readout parameters, retaining RUN-180's direct fitness loss.
[The joint protocol](joint_readout.md) replaces the r grid with three five-parameter
Optuna studies, one per existing observation assumption. Gamma concentration is
disabled. The old model, inputs, masks, seeds, noise draws and confirmation design
are retained. All six old fits per mode are re-evaluated as starting candidates;
11,712 trials per mode preserve the overall 35,136-trial budget. A dedicated
Linux runtime and frozen source/input bundle support the authorized Slurm submission.

## Initial fixed-r rollback, superseded by the joint fit below

The scientist requested a rollback to the most recent useful calibration while
retaining Optuna. [The amended protocol](optuna_readout.md) restores RUN-180's
fitness-distribution/local-SD/floor objective, Gaussian profile full-pairwise
generator, reference-anchored sigmoid, original bounds, fixed r grid and all
three observation assumptions. TPE optimizes the four readout parameters, with
previous RUN-180 settings re-evaluated as candidates. The launcher now defaults
to `optuna-readout`; architecture fitting remains an explicit experimental option.
No new scientific RUN is registered or launched by this implementation. Earlier
RUN source snapshots and scientific conclusions remain unchanged.

## RUN-184 analysed in ANA-071

[ANA-071](../../analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/analysis.md) provides five reproducible figures and three executed notebooks. The selected setting does not adequately reproduce TrpB fitness spread, local sign epistasis or sparse-panel navigation; bootstrap architecture is unstable despite convergence. It is not ready to graduate as a TrpB-calibrated environment. This concerns the achieved fit, not a rejection of the model family, and does not change defaults or historical claims.

## Optuna architecture calibration authorized 2026-11-21

RUN-181 passed local preview; RUN-182 was stopped at the scientist's request,
with its partial artifacts preserved. RUN-183 passed Slurm preview. RUN-184
(job **2325**, 1 CPU, 8 GB, 24-hour limit) completed on Dutch on 2026-11-24 at
16:00:32 AEDT after about 6 h 17 min, with both exit codes zero. All 100 Optuna
attempts were processed: **37 COMPLETE, 63 FAIL** from auxiliary inference
nonconvergence. The selected trial (index 63) has calibration loss **24.0341**
and confirmation loss **27.3463** on disjoint genotypes and fresh seeds.
Selected fits pass the required convergence checks. The high numerical failure
rate limits usable search; these numbers alone do not establish model adequacy.

All 152 frozen payloads remain intact, all expected outputs are present, and
33/33 independent completion/reproducibility checks pass. See
[completion verification](runs/RUN-184.completion-verification.json),
[RUN record](runs/RUN-184.yaml), and
[execution amendment](architecture_slurm.md). The collected bundle preserves
outputs, logs, frozen runtime, exact remote transport and failure diagnostics.
Source, inputs, settings and seeds were unchanged from the cancelled local run;
no Optuna resume or preset promotion was performed.

The scientist authorized a new TrpB-only calibration RUN and local MacBook
execution using unpublished candidate 0.7.0rc1. The registered amendment is
[architecture_calibration.md](architecture_calibration.md), with the exact fitting
configuration in [architecture_calibration.json](architecture_calibration.json).
The same reference-free sigmoid-plus-pairwise inference procedure is applied to
measured and simulated data; Optuna matches architecture summaries. This changes
the fitting objective, coefficient distribution and sampling protocol relative to
RUN-180. Previous controls remain historical context, not a matched optimizer
comparison. The source, input partition and runtime environment are frozen before
submission. No new adequacy or biological-mechanism claim is made by this RUN.

## Full-pairwise TrpB follow-up authorized 2026-11-20

The scientist requested a matched repeat of profile-conditioned RUN-178 using
the refined 0.6.0rc1 full-pairwise model, with local MacBook execution. The
prospective protocol is [full_pairwise.md](full_pairwise.md) and its complete
configuration is [full_pairwise.yaml](full_pairwise.yaml). All non-model scientific
settings and inputs are retained. RUN-178 and the original synthetic RUN-175
remain archived controls and are not rerun. This follow-up does not overwrite
their evidence or create a new empirical adequacy claim.

At that stage the launcher defaulted to `full-pairwise`; historical protocols require their
recorded sources and explicit protocol selector. The local submission helper
registers a separate dry-run preview and live RUN under this EXP.

RUN-179 was the nonexecuting local preview. RUN-180 completed locally in
330.75 seconds with all 18 settings and 24 confirmation seeds. Independent
[engineering verification](runs/RUN-180.verification.json) passes: unchanged
empirical inputs, profile/additive values, references and r=0 results. The
[verifier amendment](runs/RUN-180.verification_amendment.json) documents gzip
timestamp metadata; both compressed files authenticate correctly and their
decompressed contents are identical. Original sources and outputs are unchanged.
The 346-file bundle is DVC-tracked with local cache verification; remote backup
was not attempted. This records execution, not an adequacy conclusion.

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


## 2026-11-25 Matched interaction ablation

The scientist authorized a new TrpB RUN contrasting independently refitted r=0
and full r in [0,8], with the RUN-186 shared-stop observation law fixed. See
[protocol](interaction_ablation.md) and [configuration](interaction_ablation.yaml).
Both arms receive 11,712 Optuna trials under the restored direct fitness loss.
Repeatable interaction powers and sign reversals are unfitted diagnostics.
This tests dependence on latent interactions; it cannot identify common assay
bias separately from biology. No package preset or historical claim changes.

RUN-188 is the production comparison (Slurm 2356, running on Dutch); RUN-187 is
the nonexecuting preview. Local/Linux smoke and 24 tests pass. Frozen runtime
and backup objects are verified; launch evidence is in
[runs/RUN-188.launch_verification.json](runs/RUN-188.launch_verification.json).
No production conclusion is drawn at launch.

ANA-073 completes the requested RUN-188 noise-versus-structure analysis with six
figures and an executed notebook. The r=0 observed control reproduces much of the
sign-epistatic geometry. At primary resolution, 90.4% of full-model observed
sign-epistatic faces were absent intrinsically; this does not mean 90.4% of
spectral power or biological epistasis is noise. Shared error can be repeatable;
latent pairwise terms do not uniformly improve the unfitted diagnostics.
See [ANA-073](../../analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/analysis.md). No package/default or historical claim change.

## Gamma concentration follow-up (25 November 2026)

Scientist-authorized joint alpha/r/readout calibration with unchanged shared-stop observations and direct loss; reuse RUN-188 Gaussian and r0 controls. See [protocol](gamma_concentration.md), [configuration](gamma_concentration.yaml) and [validation](gamma_concentration_validation.json). This is matched diagnostic confirmation on previously inspected TrpB panels, not new independent empirical validation.

RUN-189 is the nonexecuting preview. RUN-190 attempt 1 submitted as Slurm job **2357**, verified running on Dutch. Local tests and Linux smoke passed; source/input payloads and pinned environment verified at production startup. [Launch evidence](runs/RUN-190.launch_verification.json). No production result yet.
