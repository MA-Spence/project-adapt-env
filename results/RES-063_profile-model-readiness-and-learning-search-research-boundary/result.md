# RES-063: Profile model readiness and the learning–search research boundary

## Summary

The locked profile-conditioned model is a usable, reproducible research environment
with an interpretable latent interaction-strength parameter and characterized
four-site TrpB regimes. The evidence does not establish comparative computational
superiority, independent control of learning and search difficulty, or explanatory
transfer to biological systems. Those are the prospective questions in HYP-010,
HYP-011 and HYP-012.

This result synthesizes existing RUN-175/ANA-069 and RUN-178/ANA-070 evidence,
implementation invariants and the scientist's authorized interpretation. It reports
no new calibration, scientific experiment or ML benchmark. "Good enough" means a
working environment for those tests, not certified biological realism. The profile
choice adds marginal evolutionary information; its empirical superiority is unproven.

## Provenance and existing evidence

- Aim: AIM-001; experiment: EXP-065; runs: RUN-175 and RUN-178.
- [Synthetic control analysis, ANA-069](../../analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/analysis.md).
- [Matched profile analysis, ANA-070](../../analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/analysis.md).
- [Profile-versus-control figure](../../analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/figures/01_paired_inherited_regime.png),
  [interaction sweep](../../analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/figures/02_interaction_sweep.png),
  [observation sensitivity](../../analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/figures/03_observation_models.png).
- Those figures retain their [executed reproduction notebook](../../analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/notebooks/analysis.ipynb)
  and [independent validation notebook](../../analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/notebooks/validation.ipynb).
- [Release decision and lock](../../docs/model_release/0.5.0rc1/decision.md),
  [engineering verification](../../docs/model_release/0.5.0rc1/verification.json),
  and [frozen-source numerical replay](../../docs/model_release/0.5.0rc1/replay.json).

ANA-070 changes the additive prior and repeats the same calibration procedure;
readouts are refitted, not held fixed. It finds mixed descriptor trade-offs. At the
inherited r=0.5/shared-stop setting, profile median walk length is 8.749 against
7.459 for the synthetic control and 8.589 empirically. Pooled high-order fractions
are 0.087, 0.176 and 0.086 respectively, but only 2/8 and 4/8 paired fresh seeds
improve those two errors. Equal-cube epistasis and sign-epistasis medians move farther
from the empirical values. There is no superiority or equivalence conclusion.

The following are existing ANA-070 fresh-seed medians, at resolution 0.001 times
the empirical calibration range. Eight generator seeds are not biological replicates.

| r | No-noise equal-cube epistatic power fraction | No-noise sign epistasis | No-noise uphill steps | Shared-stop equal-cube epistatic power fraction |
|---|---:|---:|---:|---:|
| 0 | 0.129 | 0.000 | 8.805 | 0.660 |
| 0.5 | 0.211 | 0.289 | 11.247 | 0.590 |
| 1 | 0.310 | 0.359 | 10.418 | 0.539 |
| 2 | 0.462 | 0.471 | 8.195 | 0.560 |
| 4 | 0.523 | 0.497 | 6.551 | 0.624 |
| 8 | 0.545 | 0.505 | 6.055 | 0.668 |

[metrics.json](metrics.json) selects the underlying existing rows and hashes their
ANA-070 source. No-noise and shared-stop entries have separately fitted readouts.
This demonstrates useful behavior variation, not a controlled proof that a single
parameter monotonically changes all observed ruggedness metrics. Shared-stop remains
an observation sensitivity, not an identified assay-error mechanism. The public
preset is intrinsic/noiseless and does not silently implement that research layer.

## What the model parameters actually control

The latent coordinate is a(x) + r g(x). Both fields have unit single-substitution
RMS under the declared uniform q20 mutation operator. Thus r is an interpretable
interaction/additive RMS ratio, not an observed epistatic-variance percentage.

| Parameter family | Scope and limitation |
|---|---|
| Interaction ratio r | Changes relative interaction strength and can change ranks, mutation signs, peaks and paths; no universal monotonic difficulty guarantee. |
| Orders and energy weights | Allocate latent mutation energy to interaction orders, not final observable spectral fractions. |
| Support degrees | Control interaction incidence; discrete and capacity-limited. Four-site TrpB already contains all six pairs, four triples and one quadruple. |
| Representation rank | Supported values 1 and 8; changes interaction-pattern richness at fixed normalized energy. |
| Profile geometry | Changes relative site/residue preferences. No inferred biological interaction identities are introduced. |
| Sigmoid occupancy and steepness | Change magnitudes, saturation and spectra while preserving exact intrinsic ordering and zero-threshold uphill topology. |

Positive output scaling also preserves exact ordering. Thresholds, numerical
saturation and observation noise can alter measured ties and apparent topology.
Length and seed define scale and realization; batch size is computational.

For positive final profile probabilities, p'_i(a) proportional to p_i(a)^tau gives
M' = tau M plus per-site constants. Centering removes those constants and global
mutation-RMS normalization cancels tau. A common positive profile temperature is
therefore not an independent conservation knob. Site-specific concentration and
smoothing can change relative preferences. The generic r=4 config is not the
calibrated TrpB r=0.5 preset; changed configurations do not inherit that calibration.

## Interpretation: recovery, navigation and biological information

Learnability is the sample/computation requirement for a specified learner to
recover fitness under a specified query distribution, test distribution, information
access, noise and error criterion. Searchability is the success or cost of a declared
search policy, starting distribution, allowed moves and fitness target. Finding any
improvement, a local peak, a top quantile and the global optimum are distinct tasks.

| Recovery | Navigation | Candidate explanation, not an established current regime |
|---|---|---|
| Easy | Easy | Informative additive direction and a simple predictive function. |
| Easy | Hard | Compact estimable interactions with conflicting effects and local traps. |
| Hard | Easy | Detailed prediction is demanding while useful directions or good regions remain accessible. |
| Hard | Hard | Both inference and navigation are limiting under the declared budget. |

Four-quadrant coverage is a hypothesis, not a result. Exact recovery with unlimited
compute allows enumeration; the discordant easy-learning/hard-search cell requires
a computation or navigation restriction. Negligible residuals and saturated labels
must not manufacture dissociation through arbitrary accuracy or success criteria.
Separate fixed-data prediction, true-oracle navigation and model-guided optimization.
Charge every oracle query, including enumerated neighbours. Global regret requires
an enumerated or certified optimum; the model reference is not one.

Biological information used to construct the oracle differs from information available
to a learner. The exact generating profile reveals the additive component by
construction; its practical utility remains task- and regime-dependent. A direct
profile baseline is necessary before assigning extra value to a PLM.
For a controlled information intervention, a site-wise alphabet bijection of the
whole oracle f_pi(x)=f(pi^-1(x)) preserves the Hamming-graph landscape. Map reference,
starts, masks, training/test points and static noise consistently. It disrupts
compatibility with fixed biological priors without changing peaks or paths.
Permutation-closed learners retain intrinsic learnability; this intervention measures
prior-dependent learning. Generic navigation is invariant only with label-equivariant
proposals and consistent tie-breaking. Shuffling only the additive matrix is not an
exact geometry control. This wrapper is proposed, not implemented or tested here.

TrpB PLM inputs require the full author protein context and verified mapping of sites
183, 184, 227 and 228. The relabelled control is deliberately biologically misaligned,
not another natural family. Current interactions are synthetic: no result from this
control alone would demonstrate that a PLM learned empirical epistasis.

## Computational contribution and package status

At fixed alphabet, maximum order, degree and rank, storage and full-sequence scoring
are O(L); local mutation updates use incident interaction terms. Initialization also
constructs and sorts supports. Other bounded-neighbourhood NK/GNK implementations
can also have linear length scaling. Low-rank products restrict representational
capacity relative to complete interaction tables. Practical speed and memory advantages
require competent matched comparators, not an exhaustive 20^L lookup straw baseline.

The maintained Adapt-Env package has one current model family. The old model families,
_legacy package, compatibility implementations and old simulation/experiment packages
have been removed. The project facade delegates to the maintained package; historical
project analysis helpers are retained for frozen-run reproduction, not offered as
current model APIs. Fresh source and isolated-wheel tests pass 82/82 each; the wheel
contains only maintained modules and required preset/profile data. All 13 facade
exports and the 45-file release lock were checked. [package_audit.json](package_audit.json)
records these engineering checks. An ambient editable-install hook points to another
old checkout; isolated validation prevents that environment from being mistaken for
this distribution. These checks are not computational performance comparisons.

The separate project-wide research overlay is not validated by these package tests.
Its historical dependency lock omits the current Landscapy dependency closure, and
its configured Python 3.11 base conflicts with that checkout's Python >=3.12
requirement. These limitations are documented in environment/README.md; historical
RUN environments and the tested standalone model installation remain distinct.

## Effect on hypotheses and claims

- HYP-005 retains narrow component evidence and remains open as a broad monotonic
  control claim. This synthesis does not mark it newly supported or refuted.
- HYP-010: comparative computability and efficiency; prospective, untested.
- HYP-011: separable learning/search regimes; prospective, untested.
- HYP-012: biological information and conditional transitions; prospective, untested.
- CLAIM-004 retains its original synthetic-prior conditional TrpB scope. This result
  can supply its previously absent RES-level synthesis without expanding that claim
  to the profile model, ML performance or biological transfer.

The scientist's scope extension is recorded in AIM-001. This does not reinstate retired
routing implementations or overturn historical failures. No EXP or RUN is created by
this result. Synthetic interventions can explain the constructed model; explanations
of real biological phenomena require predictions tested on independent measured
fitness landscapes. Fresh generator seeds do not supply that biological validation.

## Literature and non-redundancy

- [Brookes, Aghazadeh & Listgarten (2022)](https://doi.org/10.1073/pnas.2109649118)
  connects GNK structure to Fourier sparsity and function-recovery sample complexity
  under explicit assumptions. Low-rank factor count and a sigmoid-transformed fitness
  map do not automatically inherit those recovery bounds.
- [Kaznatcheev (2019)](https://doi.org/10.1534/genetics.119.302000) establishes worst-case
  local-search difficulty for compact NK landscapes; this is not evidence that our
  typical random realizations are hard.
- [Stanton et al. (2024)](https://arxiv.org/html/2407.00236v1) supplies inexpensive
  Ehrlich sequence objectives with adjustable search difficulty and known attainable
  optima. Online evaluation and tunability are established benchmark features.
- [Greenman, Amini & Yang (2025)](https://doi.org/10.1371/journal.pcbi.1012639) separates
  model-improvement acquisition from optimization and shows that uncertainty metrics
  do not simply determine useful acquisition behavior.
- [Hopf et al. (2017)](https://doi.org/10.1038/nbt.3769) provides established evolutionary
  profile and pairwise modelling context. Marginal alignment information is useful
  infrastructure, not by itself a new model class or demonstrated advantage.

The prospective contribution is an inspectable online environment that can separate
inference, navigation and biological side-information effects. Its utility and
comparative advantage are to be tested, not asserted by this record.
