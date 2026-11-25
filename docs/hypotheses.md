# Hypotheses

2026-11-20 implementation boundary: unpublished 0.7.0rc1 adds an explicit Gamma
effect-concentration control and an Optuna indirect-inference calibration method.
These are tools for the existing HYP-004/005 and prospective HYP-011/012 questions,
not new empirical support for any hypothesis. Genetic-architecture summaries are
now fitting targets; they cannot also count as unfitted validation under this new
protocol. Historical EXP-065 objectives and conclusions are preserved. No new
experiment or empirical RUN is created by implementing this method.

2026-11-20 authorized test: EXP-065 will repeat the profile-conditioned TrpB
calibration with full pairwise interactions, preserving RUN-178's non-model
design. The repeat is relevant to HYP-004/005's scoped aggregate-statistics and
search questions. It does not test HYP-010 efficiency or HYP-011/012 learning
claims, and no new support is asserted before the outputs are analysed.

2026-11-19 model-version boundary: the new unpublished 0.6.0rc1 default uses
full pairwise tables and a sigmoid; the earlier ranked model remains in the
named calibrated presets. Historical evidence below does not validate this new
default. HYP-010/011/012 remain prospective, and EXP-067 still targets its
explicitly frozen 0.5.0rc1 construction pending a protocol amendment.

## HYP-001: Local-Protein-Landscape-Statistics-Without-Exact-Reconstruction

- Status: superseded
- Superseded by: CLAIM-001
- Owner: matthew-spence
- Aims: AIM-001

A hybrid synthetic landscape built from a stability threshold, a low-dimensional
functional layer, and sparse structured epistasis can match empirical local
summary statistics from real protein landscapes without attempting to
reconstruct any one assay exactly.

Predictions:
- A shared synthetic parameter regime should reproduce DFE shape, neutrality,
  beneficial fraction, and fitness-dependent mutational-effect trends across a
  panel of empirical assays at the level of summary statistics.
- The match should hold without assay-specific variant-level fitting.

Why this matters:
- The benchmark goal is not exact replay of one DMS experiment, but a
  synthetic family that still looks protein-like to learning algorithms and
  evolutionary processes.

Test:
- Choose a panel of empirical assays spanning stability, binding, and activity.
- Compare synthetic landscapes to empirical assays using summary statistics
  rather than variant-by-variant fit.
- Fit broad parameter ranges, then evaluate on held-out assays or held-out
  proteins.

Primary readouts:
- single-mutant DFE shape
- neutral, deleterious, and beneficial fractions
- DFE skewness and tail weight
- fitness-dependent mutation effects
- local epistasis prevalence and sign structure

Falsification criteria:
- no shared parameter regime can reproduce the empirical summary-statistic
  envelope across multiple assay families
- realism only appears after assay-specific tuning

2026-05-15 compatibility note:
- `EXP-029` was added under this historical hypothesis to test a narrower
  stale-cache confound exposed by `RES-028`: recent empirical graphlet
  comparisons used cached ProteinGym latent targets prepared under older
  epistasis/readout semantics. The experiment requires raw ProteinGym systems
  to be refit with the current smooth-core plus fixed-scale residual model
  using ABC before any empirical parameterized information is used.
- This does not undo the `CLAIM-001` supersession. Even a positive `EXP-029`
  result would support only a refit/current-semantics compatibility argument
  unless shared or hierarchical parameters generalize across assays.

## HYP-002: Benchmark-Certified-Peak-Mode-Preserves-Local-Realism

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

A benchmark-facing landscape mode that certifies or tightly controls the top of
the landscape will make optimization metrics interpretable while preserving the
local DFE and epistasis properties needed for realism.

Predictions:
- A certified-top mode should reduce peak overshoot and reference-to-peak
  ambiguity.
- Local DFE, global epistasis, and neutrality metrics should remain within
  pre-specified tolerances of the biological mode.

Why this matters:
- Fraction-of-peak and regret-style metrics are not rigorous if the stored peak
  is only a local-search artifact.

Test:
- Compare the current public reference-neutral mode against a
  benchmark-certified mode such as strict peak locking or peak-preserving
  pairwise construction.
- Evaluate both modes on the same seeds and alignments.
- Require local realism metrics to remain within pre-specified tolerances.

Primary readouts:
- certified top-sequence recovery under exhaustive search on small landscapes
- overshoot rate in probe searches
- reference-to-peak distance
- DFE, global epistasis, and neutrality changes between modes

Falsification criteria:
- benchmark certification materially distorts local mutational statistics
- algorithm rankings change only because the certified mode becomes
  artificially easy

## HYP-003: Real-MSA-Conditioning-Improves-Realism-Over-Synthetic-Conditioning

2026-11-19 implementation decision: ANA-070 compares the TrpB alignment-profile
additive prior with the recorded random-additive control, finding mixed aggregate
descriptor trade-offs. It does not test this hypothesis's conservation-effect,
covariation-placement or natural-sequence predictions. The scientist nevertheless
selected the profile-conditioned variant as the default for unpublished 0.5.0rc1
because it carries biological marginal preferences. This implementation choice
does not mark HYP-003 supported or validate ML/pLM transfer.

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

Conditioning the landscape on real family MSAs will improve
conservation-sensitivity structure and pairwise epistasis placement more than
conditioning on internally generated synthetic alignments.

Predictions:
- Real MSA conditioning should outperform synthetic-alignment conditioning on
  held-out conservation-effect correlation.
- Real MSA conditioning should outperform synthetic-alignment conditioning on
  covariation-epistasis enrichment and natural-sequence plausibility metrics.

Why this matters:
- Current internal tests mostly show that the model is consistent with its own
  synthetic priors, not that those priors match real protein families.

Test:
- Assemble a panel of real MSAs linked to ProteinGym or other mutational
  datasets.
- Compare three conditions: no alignment, synthetic alignment, real MSA.
- Evaluate only on summary statistics not used to set the parameters.

Primary readouts:
- conservation versus mean mutational effect correlation
- enrichment of strong epistasis among high-covariation site pairs
- viability of natural-like sampled sequences
- calibration of mutation effects for consensus, natural, and off-manifold
  sequences

Falsification criteria:
- real MSA conditioning offers no improvement over synthetic conditioning
- alignment conditioning improves only trivial viability checks but not
  epistasis or sensitivity structure

## HYP-004: Adaptive-Walk-Dynamics-And-Peak-Proximal-Epistasis-Are-Biological

2026-10-08 claim: CLAIM-004 records useful conditional joint TrpB aggregate and
adaptive-search resemblance from ANA-069/RUN-175. It explicitly includes the
observation assumption, seed variability and pooled-WHT mismatch. This supports
the scoped benchmark aim, not all statements of HYP-004 or arbitrary-length
biological realism. Adapt-Env 0.4.0 preserves the intrinsic model while removing
retired implementations; cleanup does not supply new scientific evidence.

2026-10-08 interaction-strength analysis: ANA-069 (EXP-065 / RUN-175) finds that
lower interaction strength improves several GB1 global statistics and walk lengths
jointly, while confirmation-cube local SD remains too large. In TrpB, strength
alone trades epistasis against search length without matching upper-1% success;
a fixed shared-stop observation sensitivity at r=0.5 gives useful coarse joint
resemblance. Noise already contributes substantial apparent epistasis at r=0,
and fresh seeds weaken the earlier r=4 success-rate agreement. This supports
conditional synthetic benchmark regimes with explicit seed and observation
variation, not identification of the biological interaction mechanism or a
universal realism claim. The finite-budget fits are not proven optima; no gate
is used to interpret these results.

2026-10-08 whole-library analysis: ANA-068 reuses the frozen EXP-066 / RUN-172
fits for GB1 and TrpB. High-fitness sets are broadly reachable, but uniform-uphill
walks are shorter in the r=4 noiseless model (median 5.27 versus 9.74 empirical
substitutions for GB1; 6.37 versus 8.59 for TrpB), and basin probabilities differ.
TrpB shared-stop observations bring upper-1% walk success closer (70% versus
65%) while further shortening walks. Neutral components depend strongly on
resolution, assay floors and observation assumptions; GB1's large exact-tie
component is entirely at zero fitness. These are descriptive post-fit comparisons
on measured four-site graphs, not independent whole-library validation or a
test of arbitrary-length invariance. No scientific gate is applied.

2026-10-08 analysis: ANA-067 (EXP-066 / RUN-172) finds conditional aggregate
resemblance, strongest for TrpB under an explicit observation sensitivity.
Systematic local-spread and sign-epistasis deficits remain across many stability
panels. Only 11 of 63 evaluated assays inform order >=3 WHT; this experiment
does not test adaptive-search performance or statistical invariance with length.

2026-10-08: EXP-066 extends the conditional aggregate-statistics comparison to
available assay classes with measured combinatorial coverage. It does not test
ML policy performance or recover empirical interaction identities.

2026-10-07 prospective scope: EXP-065 tests conditional aggregate statistical
recovery on measured complete subspaces, including epistasis, diminishing returns,
and accessibility. This is a limited test relevant to HYP-004, not validation of
full adaptive-walk dynamics or empirical interaction identities. Its independently
calibrated baselines, disjoint cubes, and fresh seeds replace historical gate
requirements for this newly authorized experiment only. No new supporting result
is claimed by implementing or registering it.

2026-10-07 analysis update: ANA-063 independently audits EXP-065/RUN-164 and
finds partial aggregate recovery. Composite has the lowest held-out point
distance for GB1, whereas bare GNK is lowest for TrpB and CR9114-H1; descriptor
mismatches and resolution-dependent rankings remain. GB1 WHT and beneficial
slopes have substantial missingness at primary resolution. This does not validate
full adaptive-walk dynamics or a universal diminishing-returns rate.

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

A biologically useful synthetic landscape should reproduce evolutionary search
behaviour seen in real rugged landscapes: random and adaptive walks should show
plausible accessibility, diminishing returns near high-fitness genotypes, and
systematic changes in epistasis as sequences approach a peak.

Predictions:
- Random walks and one-step adaptive walks should show fitness trajectories with
  early accessible gains, decelerating improvement near high-fitness regions,
  and nonzero probability of stalling on local optima.
- The beneficial fraction and mean available beneficial effect size should
  decline as sequences get closer to the relevant local peak.
- Epistasis should become more constraining near peaks, with increased
  negative or diminishing-returns epistasis and changed sign-epistasis
  prevalence relative to distant low-fitness regions.
- These patterns should be reproducible across seeds and family-like
  backgrounds without evaluating machine-learning models or benchmark strategy
  rankings.

Why this matters:
- `AIM-001` asks whether the synthetic system behaves like a biological and
  evolutionary landscape. Algorithm performance separation is a downstream
  benchmark question and is outside the current scope.
- Evolutionary walk behaviour, accessibility, and peak-proximal epistasis are
  more direct tests of biological realism than optimizer ranking.

Test:
- Sample neutral random walks and adaptive walks from reference, natural-like,
  low-fitness, and peak-proximal starting states.
- Stratify local neighbourhoods by current fitness, distance to the relevant
  local peak, and fraction of peak fitness.
- Compare walk trajectories and epistasis-distance trends to empirical or
  literature-supported qualitative expectations rather than to benchmark
  scores.

Primary readouts:
- beneficial, neutral, and deleterious fractions along walks and by distance to
  peak
- distribution of fitness gains per accepted step and adaptive-walk length
- probability of stalling, local maximum density, and accessible path count
- pairwise epistasis magnitude, sign epistasis, and diminishing-returns
  epistasis as functions of current fitness and distance to peak

Falsification criteria:
- random or adaptive walks improve or stall in ways inconsistent with
  biological rugged landscapes, such as unbounded gains near peaks or no
  deceleration in adaptive gains
- beneficial fraction and effect-size distributions do not change with
  closeness to peak
- epistasis-distance trends are absent, reversed, or dominated by arbitrary
  seed effects

Evidence update:
- `RUN-067` and `RUN-069` show smooth-core diminishing-return proxies along
  sampled TEM-1 adaptive walks, but neither run validates a clean late-walk or
  peak-proximal increase in SE/RSE or high-order WHT/GFT fraction.
- `RES-027` weakens the strong form of this hypothesis: fixed residual scale can
  expose high-order topology, but adaptive-walk stage remains a weak proxy for
  true peak distance, and contact placement is not established as the causal
  biological structure.
- `RES-028` further weakens the strong benchmark-adequacy form of this
  hypothesis. Empirical PHOT graphlets are not matched by the current embedded
  procedural candidates as a full package of spectra, SE/RSE, local maxima,
  accessibility, and random-walk statistics; the closest model is a
  residual-only diagnostic rather than the smooth-core plus contact-residual
  generator. The result falsifies the current implementation as an empirically
  validated general graphlet benchmark model, while not by itself refuting every
  broader adaptive-walk prediction in this hypothesis.
- `RES-030` partially conditions the `RES-028` conclusion after correcting the
  pooled/stale-cache comparison. Under system-specific current refitting, PHOT
  no longer reproduces the strongest wrong-evolutionary-dynamics falsifier: the
  fitted current model is the closest held-out graphlet-package candidate.
  However, the PHOT wrong-spectral-envelope falsifier still fires, and the
  current model is only marginally better than smooth core. The strong
  benchmark-adequacy form of this hypothesis therefore remains weakened.
- `RES-040` narrows the PHOT failure mode after exact complete-graphlet
  evaluation of `RUN-108`. The `RUN-039` full-ABC fit is not broadly wrong
  across sign-epistasis, local-optima, basin, spectral, or random-walk
  summaries, but all fixed branches fail local edge/site Dirichlet-energy
  concentration summaries. This weakens any claim that the current PHOT fit
  already captures local ruggedness localization while avoiding a broader
  model-family rejection.
- `RES-042` further weakens the PHOT adequacy claim for the current nested-min
  cheap-fitting implementation. Adding a positive affine observable readout
  and explicit epistasis-amplitude objective in `RUN-112` yielded `0/40`
  operationally equivalent rows, did not match empirical observable epistasis
  amplitude, left raw latent epistasis magnitudes tens to hundreds of units too
  large, and collapsed accessibility in the best calibrated candidate.
- `RES-043` adds a mixed but still negative current-model update. Bounding
  epistasis in trait space before the readout removed the worst component
  guardrail failures in `RUN-115` (`12/12` selected rows passed component
  guardrails), but produced `0/12` operationally equivalent rows against the
  best available `RUN-098` full-ABC PHOT comparator. The best held-out bounded
  candidate remained worse than the best full-ABC context row (`0.4269` versus
  `0.2115` test distance), still undershot observable magnitude epistasis
  (`0.1708` versus empirical `0.2926`), retained raw latent magnitude
  epistasis around `103.65`, and collapsed accessibility (`0.0120` versus
  empirical `0.1975`).
- `RES-044` further weakens the current PHOT adequacy claim after the
  cumulative A-through-D calibration diagnostic in `RUN-117`. The best full
  `A+B+C+D` candidate improved held-out distance slightly relative to
  `RES-043` (`0.4122` versus `0.4269`) and improved observable magnitude
  epistasis (`0.2426` versus empirical `0.2926`), but the run still produced
  `0/32` operationally equivalent rows against the best `RUN-098` comparator
  (`0.2115` test distance). Accessibility remained collapsed (`0.0190` versus
  empirical `0.1975`), while the best cumulative-stage raw latent epistasis
  moved from the previous too-large failure mode to near-zero raw variation.
- `RES-045` further weakens the current PHOT adequacy claim after the
  raw-latent soft-min path-forward diagnostic in `RUN-119`. The three-stage
  EXP-046 ablation produced `0/24` operationally equivalent rows. The best
  stage was `1`, with held-out distance `1.2240`, far worse than both the best
  `RUN-117` cumulative row (`0.4122`) and the best available `RUN-098`
  comparator (`0.2115`). It partially raised raw accessibility (`0.0796`
  versus empirical `0.1975`) but left raw local maxima low (`2.33` versus
  `5.17`) and raw scale pathologically large: raw single-mutant SD was about
  `934x` empirical and raw magnitude epistasis about `424x` empirical. Stages
  `1+2` and `1+2+3` were worse (`1.6872` held-out distance), had no component
  guardrail passes, and collapsed observable variation.
- `RES-046` further weakens the current PHOT adequacy claim after the
  RES-045 calibration-isolation diagnostic in `RUN-121`. The cumulative
  EXP-047 ablation produced `0/32` operationally equivalent rows. The best
  held-out stage was `A+F+B+C+D`, with test distance `0.9155`: about `25%`
  better than the best `RUN-119` raw-latent row (`1.2240`), but still worse
  than both the best `RUN-117` cumulative row (`0.4122`) and the best nested
  `RUN-098` full-ABC context row (`0.4826`). Later stages brought raw
  single-mutant and magnitude-epistasis scale near empirical values, but did
  so with failed component guardrails, weak accessibility, and severe
  observable readout collapse.
- `RES-047` further weakens the current PHOT adequacy claim at an earlier
  gate. The spectrum-calibrated sparse Walsh/Potts Generator A test in
  `RUN-123` completed successfully as an execution but failed the preregistered
  H1 forward capability gate before sequence-conditioned fitting. Spectrum
  calibration placed all normalized sparse-Walsh weight on order `3`, and the
  H1 probe produced `0/12` joint passes, failing accessibility, magnitude
  epistasis relative to single-mutant scale, and Walsh high-order fraction for
  every probe candidate. No PHOT equivalence fit was therefore attempted.
- `RES-048` further weakens the current PHOT adequacy route by isolating the
  target/measurement-operator layer. In `RUN-125`, cached and refit PHOT
  MAVE-NN latent targets were smooth/additive-like (`high_order_fraction`
  below `1.1e-12`; cached local maxima `1.0`), while known injected spectra
  were recovered through the same MAVE-NN plus dense-graphlet operator as
  nearly unimodal and epistasis-free (`recovered_high_order_fraction` at most
  `1.0e-14`, recovered magnitude-epistasis ratio at most `4.0e-7`, recovered
  local sign epistasis `0.0`). Sparse-Walsh support containment was not the
  bottleneck under the audited rule because all `2048` order-3 supports were
  contained in at least one sampled graphlet.
- `RES-056` further weakens the current PHOT measurement-operator route. The
  completed corrected `RUN-129` tested the preregistered non-additive
  MAVE-NN pairwise and black-box/GE operator profiles, but certified no
  operator: `0/2` profiles passed and `0/24` required checks passed. The
  pairwise profile was preflight-skipped for all three spectra at PHOT length
  `118`, while the black-box/GE profile failed under the installed
  MAVE-NN/Keras `MultilayerPerceptronGPMap build()` interface. This is a
  technical certification failure of the current route, not a biological
  measurement of PHOT ruggedness.
- `RES-055` weakens the current public additive-core plus GNK implementation
  after the fixed support-aware Stage 4B probe in `RUN-141`. The probe used
  GNK-support-aware, dimension-matched, alphabet-matched graphlets and
  regime-stratified anchors, but only `318/1620` monotonicity checks passed and
  only `34/1080` positive-epistasis core-vs-GNK separation rows passed. This
  is a current implementation/scale-balance failure, not a blanket refutation
  of every possible additive-core plus GNK model.
- `RES-057` further weakens the implemented public-generator form after the full
  EXP-058 readout-coupling sweep in `RUN-142`. Across all `5022` A/B/C parameter
  sets, `0/5022` passed the core guardrail because no high-fitness plateau formed
  in the re-centred margin regime (best observed high-plateau median fitness
  `0.049` on a `0-10` scale). No single-observable architecture combined a high
  plateau and peak-proximal diminishing returns with transmitted GNK-specific
  epistasis, so biologically-like adaptive-walk dynamics are not supported by the
  current observable. This localizes the failure to the saturating-readout /
  plateau-formation tension rather than to the knob (validated at `L=8` in
  `RES-049`) or the probe (fixed in `RUN-141`).

## HYP-005: Mechanistic-Knobs-Produce-Monotonic-Directional-Effects

2026-11-26 RES-064: RUN-188/ANA-073 shows that fixed shared-stop observation
errors add epistasis outside latent r control without accounting for all
synthetic interaction power. Sigmoid and observation effects must be separated
from latent architectural effects in interpreting knob responses. Current TrpB
replicates/stop controls do not uniquely identify biological versus shared-error
components. RUN-190/ANA-074's Gamma concentration fit gave no clear advantage.
This motivates the existing mechanistic-control question, not a monotonicity
claim or an additional experiment.

2026-10-08 analysis: ANA-067 shows useful observed regimes but substantial
assay-dependent residuals under the fixed r=4 production calibration. The r=0/r=4
comparison is a finite parameter comparison, not a complete capacity or knob
monotonicity test. TrpB observed resemblance does not identify intrinsic
thermodynamics; both observation assumptions and WHT weighting affect inference.

2026-10-08: EXP-066 tests useful empirical resemblance from the public normalized
thermodynamic generator with assay-specific global readouts, additive and linear
controls. It is evidence for parameterized aggregate regimes, not a complete
validation of this broad monotonicity hypothesis.

2026-10-07 prospective scope: the normalized rank-8 thermodynamic generator is now
implemented in AdaptEnv. EXP-065 directly measures joint held-out aggregate errors
and finite-length mutation scaling, without a universal-assay or routing-specific
success gate. This tests a useful consequence of controllable parameters, not the
full monotonic-knob hypothesis; RES-058/060/061 remain the component evidence and
past failed decisions are not retrospectively changed.

2026-10-07 analysis update: ANA-063 finds finite-length computational feasibility
and near-unit latent mutation RMS through length 512. The separate default
reference-rank-1 scaling panel nevertheless loses observable effects through
sigmoid saturation; at length 512 every sampled effect is zero in float64.
This limits observable scalability under those settings without invalidating
latent normalization. The four-site fits do not establish large-length empirical
recovery or the full monotonic-knob hypothesis.

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

The model knobs intended to represent biological concepts should move
validation metrics in the expected direction over broad parameter sweeps.

Predictions:
- Increasing stability margin should change neutral and lethal fractions in the
  expected direction.
- Functional dimensionality, epistasis strength, covariation bias, and noise
  settings should each induce directional changes in their target validation
  metrics with effect sizes larger than seed-to-seed variation.

Why this matters:
- If the knobs are not directional, the model is hard to interpret
  scientifically and hard to use for controlled benchmark generation.

Test:
- Run preregistered sweeps over stability margin, functional dimensionality,
  epistasis strength, noise amplitude, and covariation bias.
- Estimate effect sizes and monotonicity with confidence intervals over many
  seeds.

Primary readouts:
- stability margin versus neutral fraction and lethal fraction
- functional dimensionality versus beneficial fraction and adaptive
  accessibility
- epistasis strength versus local maxima count and accessible paths
- covariation bias versus enrichment of strong interactions at coupled sites

Falsification criteria:
- the expected directional effect is absent or reverses across nearby settings
- landscape behavior is dominated by seed-specific stochasticity rather than
  parameter meaning

Evidence update:
- `RES-027` supports this hypothesis for the fixed residual/additive scale knob:
  in `RUN-069`, increasing fixed contact residual scale from `1` to `4` to `16`
  raises held-out SE, RSE, local maxima, high-order WHT/GFT fraction, and
  residual high-order exposure, while shortening random-walk lag-1
  autocorrelation.
- This support is narrow. It establishes a directional procedural knob in the
  sampled graphlet diagnostic, not broad seed robustness or contact-specific
  biological causality.
- `RES-028` preserves this narrow support across the empirical falsification
  suite. In `RUN-071`, the same fixed scale grid remains directional across two
  synthetic landscape seeds, but the result should still be cited only as
  support for the residual-scale knob, not as validation of the full empirical
  procedural model.
- `RES-049` supports a narrow enumerable hashed-GNK knob map at `L = 8`, where
  epistasis weight moved ruggedness, spectrum, accessibility, sign epistasis,
  and autocorrelation metrics in the expected directions.
- `RES-055` weakens the broad protein-scale version of this hypothesis for the
  current public generator. In `RUN-141`, the design-matched support-aware
  probe passed only `318/1620` monotonicity checks (`0.1963`), with high-order
  fraction at `36/270` and local-maxima density at `28/270`. The fixed probe
  therefore did not restore a reliable directional `epi_weight` knob at
  protein scale.
- `RES-057` weakens the protein-scale form after the full EXP-058 sweep in
  `RUN-142`: `11,883/58,320` monotonicity checks passed (`20.4%`) and no fix
  reached a majority without wrong-direction medians. Topographic/accessibility
  knobs (local-maxima density, random-walk autocorrelation, peak reachability)
  were poorly or wrongly controlled, and Fix A's `high_order_fraction` moved the
  wrong way with GNK weight. The narrow exception is that Fix A's sign- and
  reciprocal-sign-epistasis knobs were cleanly monotone (`99.7%`/`75.6%` pass,
  median `rho ~ 1.0`), and Fix B mixed routing transmitted large observable
  high-order/sign epistasis (`|full - core|` deltas up to `0.87`/`0.16`) though
  not monotonically in the weight knob.
- `RES-058` supports the latent protein-scale form for an experiment-local,
  mutation-operator-normalized GNK field. In `RUN-149`, all `26/26` required
  latent gates passed across six fixed references (`L = 48-448`) and 30 GNK
  seeds. The candidate local-effect slope was `-0.00632` (95% paired
  seed-bootstrap interval `[-0.01240, 0.00645]`) versus `-0.51036` for the frozen
  Stage-4 384-term raw-amplitude control; field variance scaled with slope
  `0.99508`, and reference-shell GNK RMS slopes were `0.48196` at `0.1L` and
  `0.43034` at `0.5L`. The locally normalized fixed-support control also repaired
  local, field, and full-distance scaling, so operator normalization—not extensive
  support by itself—is sufficient for the scale-law repair. Extensive balanced
  support and bounded contrasts add complete site coverage and controlled tails.
  This support stops before sigmoid transmission and does not reverse the
  observable-level failures in `RES-055`/`RES-057`, validate the public method, or
  test HYP-004/HYP-009. One protein per length also means the seed-bootstrap
  inference is conditional on the fixed reference panel.
- `RES-059` weakens the observable protein-scale form after the registered
  sigmoid-transmission test in `RUN-151`. All 12 technical checks pass, pointwise
  scalar transfer has slope `-0.00205` (95% paired-seed interval
  `[-0.01280, 0.00703]`), and frozen-core GNK texture passes at every length with
  slope `-0.10100` (`[-0.11819, -0.08858]`). The sigmoid can therefore carry the
  RES-058 texture in a controlled transition window. With the real additive core,
  however, only about `1.4-1.5%` of complete reference single-mutant panels remain
  both-dynamic, actual-core texture has zero practical-resolution fraction at all
  six lengths, and margin-8 edge/texture practical fractions are at most
  `0.0010/0.0254`. The formal decision is `SIGMOID_TRANSFER_FAILED` (`13/16`
  required checks). This localizes the failure to common-threshold additive-core /
  steep-readout coupling rather than renewed latent scale collapse, but it leaves
  the tested observable without a broadly directional practical GNK knob. The
  experiment-local candidate, noiseless heuristic floor, one-protein-per-length
  panel, and lack of public-generator/service or optimizer testing remain binding
  scope limits.
- `RES-060` records a mixed but formally negative update from the dimensionless
  calibration and held-out confirmation in `RUN-153`. One global additive
  mutation scale (`sigma_A = 2.74486`) and the core-only selected
  `beta=4`, `mu=logit(0.99)` pass all `60/60` held-out core cells, and all
  `360/360` complete-h1 edge cells pass with median practical fraction `0.8828`
  and correlation `0.8717`. GNK-specific actual-core texture is aligned and
  non-negligible (median RMS `0.4175`, correlation `0.7322`, projection `0.9182`,
  permutation `p=0.002` throughout), but only `24/360` cells meet the
  preregistered practical-prevalence floor; the frozen-reference control meets it
  in `0/360` cells because it is evaluated at the saturated `0.99` reference
  occupancy. Exact per-landscape normalization also fails (`23/360` practical
  passes), so it is not the missing ingredient. The formal decision is
  `CALIBRATED_SIGMOID_TRANSMISSION_FAILED` (`6/9` checks); an independent
  `1e-12` equality tolerance corrects a floating-point `rho=0.80` directionality
  artifact and gives `7/9` without changing the decision. This weakens the full
  observable texture-prevalence form of HYP-005 while preserving narrow support
  for global core calibration, one-step transmission, dose directionality, and
  registered length stability on the fixed 12-protein noiseless panel. It does
  not validate the public generator or lift the Stage 4C/4D/4E hard stop.
- `RES-061` supports a narrow observable protein-scale form on the same fixed
  challenge panel after the paired categorical-representation intervention in
  `RUN-157`. The fresh binary-rank-1 arm reproduces the prevalence failure
  (`45/360` cells; median practical fraction `0.1979`), whereas binary rank 8
  passes all `360/360` cells with median practical fraction `0.4609` and a
  paired protein-level gain of `0.2624` (bootstrap lower bound `0.2480`) at
  almost unchanged texture RMS. The repair redistributes fixed second-moment
  energy away from near-zero/heavy-tailed plaquette effects rather than raising
  the dose. Dense q20 rank 1 is worse (`0/360` passes), and dense q20 rank 8
  remains below the per-length package gate, so removing binary ties is not
  sufficient. The registered decision is `BINARY_RANK8_SUFFICIENT`. This is
  experiment-local mechanistic support for a dose-directional GNK knob, not
  independent protein generalization, a minimum-rank result, or validation of
  the public method; the Stage 4C/4D/4E hard stop remains.
- `RES-062` weakens the broader forward-coverage form after the operator-free
  observable coverage probe in `RUN-162`. The independently reproduced
  decision is `OBSERVABLE_CONSTRUCTION_LIMIT`: G1 range coverage passes for
  only `1/3` rugged-end templates, G2 joint reachability and G4 navigability
  pass for `2/3`, and G3 specificity passes for `0/6` descriptors. Across all
  10 Tier-1 templates, local-maxima density is covered in `0/10`, while sign
  and reciprocal-sign epistasis are covered in only `2/10` each; structured
  maxima generally fail to separate from matched shuffled/random controls.
  This does not establish a general sigmoid or rank-8 impossibility because the
  preregistered rugged-end target is plateau-confounded: CR9114-FLUB is flat in
  all `96/96` frozen graphlets and CR9114-H3 in `69/96`, yet the strict
  no-fitter-neighbour rule assigns flat 4-cubes local-maxima density `1.0`.
  The formal expressivity/specificity prerequisite fails and fitting remains
  blocked, but the strongest biological ruggedness-limit interpretation is
  withheld. The public method and Stage 4C/4D/4E hard stop remain unchanged.

## HYP-006: Observation-Noise-Is-Realistic-And-Does-Not-Destabilize-Conclusions

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

A useful active-learning benchmark needs an observation model that perturbs
ranking and uncertainty in a realistic way, but does not erase meaningful
differences between strategies.

Predictions:
- Empirically grounded noise models should change local rankings and search
  difficulty.
- Strategy ordering and top-k recovery should remain statistically stable over
  practical noise ranges.

Why this matters:
- The benchmark is meant for iterative decision-making pipelines, not only a
  noiseless oracle optimization task.

Test:
- Evaluate fixed, heteroscedastic, and assay-inspired noise models.
- Compare inferred and true top-k recovery, regret, and ranking stability.
- Calibrate noise scales against replicate variability from empirical assays
  where available.

Primary readouts:
- top-k recovery under noisy measurement
- Kendall or Spearman agreement between noisy and true ranking
- rank reversals among close variants
- strategy ordering stability across noise regimes

Falsification criteria:
- small noise changes cause large and erratic benchmark reordering
- the default noise model is not empirically defensible or is too weak to
  change search behavior

Evidence update:
- `RES-056` weakens the current measurement-model route linked to this
  hypothesis. `EXP-050/RUN-129` did not certify either tested non-additive
  MAVE-NN PHOT operator: the pairwise profile was support/preflight-infeasible
  at sequence length `118`, the black-box/GE profile failed under the installed
  MAVE-NN/Keras interface, and `0/24` certification checks passed. The result
  does not evaluate observation noise or strategy stability, but it shows that
  the current non-additive observation/measurement layer is not yet a
  defensible benchmark operator.

## HYP-007: Multiple-Latent-Molecular-Phenotypes-Improve-Empirical-Recovery

- Status: superseded
- Superseded by: CLAIM-002
- Owner: matthew-spence
- Aims: AIM-001

Replacing the current single functional latent layer with multiple explicit
latent molecular phenotypes, such as folding stability together with binding or
abundance, linked to assay readouts through nonlinear thermodynamic or
measurement maps, will improve empirical recovery and realism.

Predictions:
- On matched empirical assays, multi-phenotype models should improve held-out
  single-mutant ranking, double-mutant ranking, and epistasis-prediction
  metrics relative to the current single-latent model.
- The improvement should be largest on assays where the present fitter
  collapses to zero epistasis, near-peak reference artifacts, or weak
  double-mutant recovery.
- When orthogonal measurements such as abundance plus binding or stability plus
  activity are available, a shared latent decomposition should transfer across
  readouts better than a single scalar latent score.

Why this matters:
- `RES-005` weakens the idea that the remaining failure is mainly optimizer
  collapse. Otwinowski-style folding and binding decompositions, together with
  the Lehner-group `deepPCA`, `ddPCA`, and `MoCHI` results, suggest that much
  observed epistasis may arise from nonlinear readouts of a small number of
  additive molecular traits rather than from one generic phenotype plus ad hoc
  ruggedness.

Test:
- Build paired benchmark panels with orthogonal readouts where possible, such
  as abundance plus binding or stability plus activity, and compare the current
  model against explicit two-trait and multi-trait latent variants on held-out
  single and double mutants.
- Evaluate both assay-specific fits and limited-sharing hierarchical fits to
  determine whether the latent traits transfer across related assays or
  proteins.
- Require improvement over the current model on predictive metrics and on
  mechanistic diagnostics such as reduced reference-to-peak artifacts and
  non-degenerate latent parameter posteriors.

Primary readouts:
- single-mutant holdout Spearman or Kendall
- double-mutant holdout Spearman or Kendall
- epistasis-prediction Spearman and KS
- cross-readout transfer accuracy when fitting one shared latent decomposition
- fitted latent-trait interpretability, such as shared folding effects across
  binding readouts

Falsification criteria:
- multi-phenotype models do not materially outperform the current single-latent
  model on matched assays
- any gains come only from much larger unconstrained parameterization rather
  than from transferable latent structure
- latent decompositions fail to align across orthogonal readouts or collapse
  back to effectively one trait

## HYP-008: Local-Non-Diffuse-Ruggedness-Requires-Nonlinear-Trait-Operators

- Status: superseded
- Superseded by: CLAIM-003, HYP-009
- Owner: matthew-spence
- Aims: AIM-001

Real protein fitness landscapes contain local, non-diffuse ruggedness generated
by nonlinear interactions among molecular traits, and the current smooth-core
plus diffuse residual model will not capture that behaviour as a biologically
grounded procedural benchmark.

Predictions:
- Empirical graphlets will show ruggedness concentrated in specific local
  sequence, structural, or trait regions rather than uniformly across sampled
  backgrounds.
- Nonlinear operators combining multiple latent molecular traits should match
  empirical SE/RSE, local maxima, accessibility, random-walk autocorrelation,
  and WHT/GFT spectra better than the current smooth-core-only or diffuse
  residual models.
- The improvement must persist after controlling for parameter count, global
  readout nonlinearity, residual-only overfitting, and stale empirical
  preprocessing.
- If the current model family is sufficient, matched single-trait
  stability-margin models or diffuse residual controls should equal or
  outperform local nonlinear-operator variants on held-out graphlets.

Why this matters:
- `RES-028` shows that empirical PHOT graphlets are more locally rugged and
  higher-order-heavy than the current embedded contact-residual generator, while
  residual-only diagnostics can get closer to empirical ruggedness magnitudes.
  That pattern suggests a structural misspecification: ruggedness may be
  localized around specific molecular-state interactions rather than spread as a
  broadly sampled residual field.
- This hypothesis is not the same as superseded `HYP-007`. `HYP-007` asked
  whether extra fitted molecular phenotypes improve empirical recovery as a
  general model. `HYP-008` asks whether local ruggedness itself is generated by
  nonlinear trait-composition operators that are sparse in sequence/trait space
  and therefore not captured by diffuse residual terms.

Test:
- Compare the current single stability-margin core, fixed-scale diffuse
  residuals, residual-only diagnostics, GNK/RMF-style controls, and explicit
  local nonlinear trait-operator variants on the same empirical and synthetic
  dense Hamming graphlets.
- Require calibration/test separation across graphlets and, where possible,
  across assay systems or start regions.
- Fit or tune the operator variants only through preregistered summary targets,
  then evaluate held-out graphlets for spectra, topology, and evolutionary
  dynamics.

Primary readouts:
- SE and RSE prevalence
- local maxima and basin-count distributions
- accessible path fraction and path-count reduction
- random-walk autocorrelation length
- WHT variance by order and high-order fraction
- GFT energy by graph eigenvalue band
- localization of ruggedness across sequence backgrounds, sites, structural
  contacts, or inferred trait regions

Falsification criteria:
- Local nonlinear trait-operator variants do not outperform matched
  single-trait or diffuse residual controls on held-out graphlets.
- Any improvement disappears after matching parameter count, fitting a flexible
  global readout, or controlling for residual-only overfitting.
- The inferred ruggedness is not localized, but instead behaves like a diffuse
  random residual field.
- A non-procedural GNK/RMF-style control matches empirical spectra, topology,
  and evolutionary dynamics as well as or better than the proposed procedural
  operator.

Supersession note:
- As of 2026-05-18, this hypothesis is superseded in its mechanism-specific
  form. The current evidence supports added coupled/auxiliary capacity for
  improving PHOT global/local ruggedness summaries over current/smooth
  controls, but it does not support empirical alignment/module locations as
  privileged or causal amino-acid interaction structures. `CLAIM-003` records
  the completed capacity claim. `HYP-009` carries the benchmark-facing question
  that remains: whether biologically routed ruggedness changes optimizer
  behaviour once landscape statistics are matched.

Evidence update:
- `RES-030` strengthens the motivation for this hypothesis but does not test it
  directly. System-specific current refitting improves PHOT graphlet-package
  alignment relative to stale or pooled comparisons, but the empirical PHOT
  spectral envelope remains too high-order-heavy and the fitted current model is
  only weakly separated from smooth core. This keeps local, non-diffuse
  ruggedness as a plausible structural misspecification to test next.
- `RES-031` weakens the first tested implementation of this hypothesis.
  A frozen-base PHOT coupled operator can improve held-out graphlet distance
  without destroying base correlation, but the improvement is carried by a
  no-module correlated operator and by an empirical readout-residual diagnostic,
  not by calibration-derived local empirical plaquette modules. The empirical
  module candidate does not beat no-module/random/shuffled controls under
  calibration-selected held-out evaluation, and module-energy diagnostics remain
  diffuse rather than locally concentrated.
- `RES-032` further weakens the tested hard-min weakest-link implementation of
  this hypothesis. `RUN-081` shows that hard-min coupling can improve the PHOT
  full-assay objective over frozen smooth/current controls, but the best
  hard-min held-out graphlet distance (`0.4124`) does not beat the fitted
  current model (`0.4016`) or smooth core (`0.4119`), the
  calibration-selected hard-min candidate does not generalize, the WHT
  high-order mismatch remains, and local-signal diagnostics do not validate
  empirical-like non-diffuse ruggedness.
- `RES-033` weakens the tested soft-min/Boltzmann implementation of this
  hypothesis. `RUN-084` shows that finite-temperature state competition can
  improve PHOT calibration graphlet distance (`0.6100` versus `0.7318` for
  fitted current and `0.7454` for smooth core) and can slightly improve the
  aggregate assay objective through the posterior mean (`0.4570`). However, the
  calibration-selected best particle fails held-out graphlets (`0.4434` versus
  `0.4016` fitted current and `0.4119` smooth core), the posterior mean also
  fails (`0.4543`), and local-state innovation is not empirical-like. A
  test-oracle posterior particle narrowly beats fitted current (`0.3964`) but
  was not calibration-selected, so it is capacity evidence rather than support.
- `RES-034` further weakens the budgeted module-tied soft-min implementation.
  `RUN-086` adds explicit state-difference and energy budgets, which prevents
  the large artificial innovation seen in `RUN-084`, but the selected operator
  effectively collapses back to the frozen base landscape. Calibration graphlet
  distance improves only slightly (`0.7240` versus `0.7318` fitted current),
  held-out graphlet distance still does not beat fitted current (`0.4060`
  versus `0.4016`), the shuffled-module control is better than the empirical
  module candidate (`0.4029` versus `0.4060`), innovation ratios to the
  empirical residual are near zero, and the assay guardrail worsens relative to
  fitted-current controls.
- `RES-035` clarifies the strongest positive signal from `RES-031` and mostly
  weakens the mechanism-specific claim. Replicated calibration/test graphlet
  splits show that the no-module `RUN-079` improvement is stable against
  fitted-current and smooth-core controls: median held-out graphlet distance is
  `0.4839` versus `0.5232` for fitted current and `0.5358` for smooth core,
  with win rates of `93.8%` and `97.7%`. But specificity controls fail. The
  uncoupled-module control is better overall (`0.4469` median distance, best in
  `51.2%` of splits, beating no-module in `63.5%`), and random, shuffled, and
  empirical module families also frequently beat fitted current. The original
  split (`0.3137` no-module distance) was unusually favorable. This supports a
  narrow capacity claim for PHOT graphlet-distance movement, not the local
  non-diffuse empirical-module mechanism required by `HYP-008`.
- `RES-036` weakens the current empirical-module ABC implementation while
  preserving a narrow capacity signal. `RUN-091` shows that ABC fitting can
  improve the fixed held-out PHOT graphlet distance for the empirical-module
  family (`0.2553`) relative to fitted current (`0.4016`), smooth core
  (`0.4119`), no-module (`0.3489`), and shuffled-module (`0.3535`). However,
  the random-module control is fractionally better (`0.255188` versus
  `0.255286`), the uncoupled-module control is clearly better (`0.2237`), and
  selected empirical-module energy remains diffuse (`9.11/12` effective
  plaquette modules, normalized entropy `0.937`). This indicates ABC-fittable
  auxiliary-landscape capacity, not validated local non-diffuse empirical
  module ruggedness.
- `RES-037` further weakens the current empirical-module implementation while
  preserving the motivation that PHOT has residual localized ruggedness.
  `RUN-093` reuses the fitted `RUN-091` ABC parameters and directly measures
  Dirichlet-energy localization. Empirical PHOT held-out graphlets are more
  edge/site concentrated than fitted current or smooth core (edge top-10
  energy fraction `0.3904` versus `0.3085`; site top-1 energy fraction
  `0.3863` versus `0.2645`), and the empirical readout residual has a local
  signal (adjacent-edge Moran I `0.0505`). But empirical calibration modules do
  not enrich empirical energy (edge-site enrichment `0.8624`, plaquette-pair
  enrichment `0.9495`), and the empirical-module ABC candidate is not uniquely
  closest to empirical localization. No-module, random-module, uncoupled,
  residual, and GNK controls match or exceed it on key localization summaries.
  This supports a residual-localization diagnosis, not the tested empirical
  local-module mechanism.
- `RES-038` tests the direct follow-up: add local Dirichlet-energy localization
  statistics to the PHOT ABC objective itself. The expanded objective improves
  calibration-selected fitted-family held-out distances over fitted current
  (`0.5210`) and smooth core (`0.5268`), and the selected families are closer
  than those controls on the Local7 ruggedness package. But empirical-module is
  not the held-out winner (`0.3653`): no-module (`0.2771`) and uncoupled-module
  (`0.2912`) are better, empirical-module has worse Local7 distance (`0.4623`)
  than no-module (`0.3268`) and uncoupled (`0.2900`), and its module energy
  remains diffuse (`9.34/12` effective plaquette modules, normalized entropy
  `0.944`). A test-oracle empirical-module particle exists (`0.2591`), but it
  was not calibration-selected. This preserves only an ABC capacity signal and
  further weakens the current empirical local-module mechanism.
- `RES-039` tests whether adding a nested third latent landscape to the locked
  `RUN-084` soft-min first stage rescues the PHOT mechanism. It does not. The
  calibration-selected empirical-module nested model has held-out full-objective
  distance `0.5432`, worse than fitted current (`0.5210`), smooth core
  (`0.5268`), locked first-stage (`0.5296`), no-module (`0.5372`),
  random-module (`0.4826`), and shuffled-module (`0.5110`), beating only the
  uncoupled-module control (`0.6401`). Local7 distances also fail empirical
  specificity: shuffled-module (`0.4199`) and no-module (`0.4766`) beat
  empirical-module (`0.5250`). The selected empirical-module energy remains
  diffuse (`10.56/12` effective plaquette modules, normalized entropy `0.972`).
  Test-oracle particles show unused capacity, but calibration does not select
  them. This further weakens nested empirical-module versions of `HYP-008`.

## HYP-009: Benchmark-Legible-Ruggedness-Routing-Affects-Optimizer-Performance

- Status: active
- Owner: matthew-spence
- Aims: AIM-001

For ML benchmark use, alignment- or covariation-routed coupled ruggedness does
not need to improve empirical landscape fit over random, shuffled, no-module,
or uncoupled capacity controls. Its role is to provide biologically legible
sequence structure that pretrained protein language models and evolution-aware
active-learning agents may exploit, while the global and local ruggedness
statistics are matched by construction.

Predictions:
- When coupled regimes are matched for global and local ruggedness statistics,
  PLM- or evolution-aware agents should perform differently on biologically
  routed landscapes than on random-routed, shuffled-routed, no-module,
  uncoupled, or diffuse-residual regimes.
- Non-biological optimizers should be less sensitive to routing identity once
  global ruggedness, local ruggedness, observation noise, and evaluation budget
  are matched.
- Any biologically routed advantage must survive replicated seeds, matched
  budget/noise conditions, and ablations where comparable ruggedness statistics
  are routed through random or shuffled modules.
- If routing identity changes all optimizers similarly, or if apparent gains
  are explained by easier landscape statistics rather than biological
  legibility, the benchmark interpretation is weakened.

Why this matters:
- The benchmark goal is to compare optimization, protein-engineering, and
  active-learning strategies on landscapes whose difficulty can be controlled
  relative to adaptive-walk or Wright-Fisher population evolution baselines.
- `RES-035` through `RES-039` show that coupled/auxiliary model capacity can
  improve empirical-like PHOT ruggedness summaries, but that the improvement
  does not require current empirical module locations. That is not a failure
  for benchmark construction. It means alignment conditioning should be treated
  as side information for agents, not as a fitted-landscape mechanism claim.
- Pretrained protein language models and evolution-aware agents can only
  express a biological prior if the benchmark exposes some biologically
  legible relationship between sequence context and ruggedness. Random,
  shuffled, no-module, uncoupled, and diffuse controls are therefore essential
  comparator regimes, not reasons to discard coupled landscapes.

Test:
- Freeze matched benchmark regimes: thermodynamic core, thermodynamic core plus
  diffuse residual, no-module coupled capacity, uncoupled capacity, random- or
  shuffled-routed coupled capacity, and alignment/covariation-routed coupled
  capacity.
- Match these regimes on the landscape summaries that define difficulty,
  including SE/RSE prevalence, local maxima, accessibility, random-walk
  autocorrelation, WHT/GFT spectra, Dirichlet-energy localization, and Local7
  distances where available.
- Run adaptive-walk and Wright-Fisher baselines plus active-learning and
  protein-engineering agents under the same sequence budgets, observation
  noise, starting states, and seeds.
- Stratify agents by prior class: generic black-box optimizers, local search
  and evolutionary baselines, supervised sequence models without biological
  pretraining, pretrained protein language models, and explicitly
  evolution-aware models.

Primary readouts:
- Best-so-far fitness, top-k recovery, and regret to the best evaluated or
  certified accessible genotype where that bound exists.
- Strategy-ranking stability across seeds, starts, budgets, and noise levels.
- Interaction between routing family and agent class, especially PLM or
  evolution-aware agents versus generic optimizers.
- Matched-landscape checks showing that performance differences are not caused
  by unmatched ruggedness, accessibility, autocorrelation, spectral, or local
  energy summaries.

Falsification criteria:
- PLM or evolution-aware agents do not perform differently on biologically
  routed landscapes than on random/shuffled/no-module controls once global and
  local ruggedness statistics are matched.
- Apparent biologically routed advantages disappear under replicated seeds,
  matched budgets, matched starts, or matched noise.
- Routing identity affects all optimizer classes similarly, implying a hidden
  landscape-statistic mismatch rather than use of a biological prior.
- The biologically routed regimes cannot be matched to random/shuffled controls
  on the agreed ruggedness summaries without changing benchmark difficulty.

Evidence update:
- `RES-035` to `RES-039` motivate this hypothesis but do not test it directly.
  They support coupled/auxiliary capacity for improving PHOT ruggedness
  summaries and show that random, shuffled, no-module, and uncoupled controls
  can match or beat current empirical-module routing. For benchmark use, that
  evidence separates fitted landscape capacity from biological legibility:
  empirical routing is retained as a possible agent-visible prior, not as a
  claim that the fitted landscape depends on specific local amino-acid modules.
- `RES-040` adds a matched-summary constraint for any future benchmark-facing
  PHOT ABC work. Before optimizer comparisons, the next fit should match the
  four exact local-energy concentration failures identified by `RUN-108`,
  rather than reusing the whole older validation objective or assuming that
  broad ABC failure has been shown.
- `RES-042` adds a negative matched-summary calibration result. The
  `RUN-112` positive-affine calibrated cheap protocol did not pass component
  equivalence and should not be used as the default GB1/TrpB scale-up route.
  It motivates better-conditioned matching of observable epistasis,
  accessibility, and local-energy summaries before any optimizer or PLM-agent
  benchmark comparison.
- `RES-043` narrows this methodological diagnosis. Bounded trait-space routing
  improved component guardrails relative to `RES-042`, but the `RUN-115`
  numerical fast fitter still failed full-comparator equivalence and left raw
  latent epistasis badly scaled. This keeps bounded routing as a model-level
  direction to refine, but it is not yet a validated tractable fitting protocol
  for complete GB1/TrpB landscape construction.
- `RES-044` adds a negative cumulative-fix result. The distribution/scale
  conditioned `A+B+C+D` stage improved the best held-out PHOT distance to
  `0.4122` and restored component guardrail passes, but it still failed the
  best `RUN-098` comparator and left accessibility and raw latent geometry
  mismatched. This motivates stricter latent-geometry matching before any
  benchmark-facing optimizer or PLM-agent comparison.
- `RES-045` adds a negative raw-latent path-forward result. The tested
  low-dimensional soft-min latent generator partially moved raw accessibility
  toward the PHOT target, but not together with matched raw scale, local
  maxima, observable variation, or held-out graphlet distance. Adding
  raw-latent constraints, a restricted sigmoid readout, and higher
  raw-accessibility weight worsened the selected fit. This should not be used
  as a matched-landscape construction route for optimizer or PLM-agent
  benchmarks.
- `RES-046` adds a negative calibration-isolation follow-up. The tested
  conditioning-first k-trait soft-min, structured-residual, and expanded-search
  route improved over `RUN-119` and partly repaired raw scale, but still failed
  all operational-equivalence decisions and collapsed observable variation to
  near-flat landscapes. This remains methodological motivation for better
  matched-summary construction, not evidence for optimizer or PLM-agent
  behavior under `HYP-009`.
- `RES-047` adds a negative forward-gate boundary for the sparse Walsh/Potts
  Generator A route. The tested graphlet-aligned sparse-Walsh construction had
  some ruggedness-like local structure, but it failed accessibility,
  scale-normalized magnitude epistasis, and Walsh high-order-fraction
  consistency before any sequence-conditioned benchmark landscape could be
  produced. This motivates keeping forward capability gates before optimizer
  or PLM-agent comparisons; it does not support a benchmark-facing matched
  landscape yet.
- `RES-048` adds a measurement-operator calibration warning for future
  benchmark-facing construction. The current PHOT MAVE-NN plus
  dense-graphlet target path suppresses known injected high-order/local
  epistasis, so the existing scalar high-order/local-ruggedness target should
  not be reused as a matched-landscape criterion without explicit correction
  or replacement. This motivates better-conditioned target calibration; it
  does not test optimizer or PLM-agent behavior.
- `RES-056` adds a negative non-additive operator-certification boundary for
  benchmark-facing PHOT construction. The completed `RUN-129` certified no
  pairwise or black-box/GE MAVE-NN operator, with all six fit-attempt rows
  failing and no recovered certification-metric values available. This blocks
  citing `EXP-050` as a certified non-additive measurement operator for
  matched empirical/model benchmark comparisons, but it still does not test
  optimizer or PLM-agent behaviour directly.
- `RES-055` blocks the current additive-core plus GNK public generator as a
  benchmark route before optimizer or PLM-agent comparison. The fixed
  support-aware probe in `RUN-141` produced support-aware graphlets and matched
  empirical template dimensions/alphabets, but failed the monotonicity and
  core-vs-GNK separation gates. Routing comparisons remain uninterpretable
  until the public observable has matched, controllable ruggedness summaries.
- `RES-057` keeps the current additive-core plus GNK public generator blocked as
  a benchmark route after the full EXP-058 readout-coupling sweep in `RUN-142`.
  No A/B/C architecture delivered matched, controllable observable ruggedness
  together with a high-fitness plateau, so routing/optimizer comparisons remain
  uninterpretable and the Stage 4C/4D/4E hard stop holds. Fix B does show that
  post-readout routing can transmit specific epistasis to the observable at
  protein scale, which is a lead for a future two-channel construction operated
  with the stability gate in saturation (the `CLAIM-002` multi-latent direction),
  but that is a new experiment rather than support from this run.


## HYP-010: Procedural-Landscapes-Provide-Competitive-Computation-And-Efficiency

- Status: active
- Owner: matthew-spence
- Aims: AIM-001
- Motivating result: RES-063

At a declared range of lengths, query workloads and interaction
complexities relevant to online landscape use, the current sparse low-rank
implementation provides a reproducible memory and/or time advantage over suitable
established implementations while retaining the specified landscape capabilities.
The claim concerns practical costs and a stated trade-off, not novelty of the
interaction-expansion/sigmoid family or universal dominance.

**Predictions and measurements.** For fixed alphabet, interaction order, degree
and rank, storage and full-sequence scoring scale approximately linearly with
length. Report construction/normalization cost, resident memory, single-query
latency, batch throughput and mutation-neighborhood evaluation separately. The
comparative advantage must survive initialization and sufficiently long workloads
to reveal amortization, rather than appear only in a favorable cached operation.

**Essential controls.** Include a competent bounded-neighborhood NK/GNK comparator;
it can also scale linearly with length when neighborhood size is fixed. Match
hardware, precision, batching and warm-up treatment, and document compiler/backend
differences. Distinguish exact numerical equivalence comparisons from comparisons
between different ensembles. Low rank and selected interaction support restrict
the ensemble; they are not a lossless acceleration of unrestricted GNK tables.
Report these representational differences and relevant descriptor coverage instead
of claiming equivalent expressivity from equal parameter counts. Dense exhaustive
lookup alone is not an adequate efficiency baseline.

**Falsifiers/limits.** No meaningful advantage in the declared workload; advantage
vanishes after accounting for construction, memory or baseline optimization; or
the gain depends on removing the capabilities whose retention was claimed. Failure
of a comparative speed claim does not refute O(L) scaling or software usefulness.
Any numerical criterion belongs in a future specified test, not a retrospective
threshold selected from benchmark outcomes.

**Maturity and biological requirement.** Scaling follows the fixed-complexity
representation, but comparative performance remains unmeasured. No real-data test
is needed for a narrowly computational claim; retaining biological resemblance
requires separate empirical evidence and cannot be inferred from throughput.

**Evidence status.** Motivated by RES-063. [EXP-067](../experiments/2026-11-19_EXP-067_matched-oracle-computability-and-efficiency/experiment.md) is the planned four-arm computational test; it has not been executed. No support is inferred from registration or complexity derivations alone.



## HYP-011: Learnability-And-Searchability-Form-Separable-Landscape-Regimes

- Status: active
- Owner: matthew-spence
- Aims: AIM-001
- Motivating result: RES-063

Under declared recovery and search tasks, the current landscape
family contains reproducible regimes spanning easy/easy, easy/hard, hard/easy and
hard/hard recovery/search behavior. Landscape structure can explain differences
between the tasks beyond a single ruggedness statistic or interaction-strength
parameter. This predicts coverage of the quadrants, not arbitrary independent
control of both axes or universal labels across algorithms and budgets.

**Predictions and measurements.** Recovery learning curves and fitness-target
attainment curves need not rank regimes identically. Interaction support, order,
amplitude and rank may affect representation complexity and navigation differently.
Measure recovery on a fixed external test distribution and search on its declared
target; retain seed/start variation and algorithm dependence. Distinguish sample
complexity from training/inference compute and from optimization compute.

**Essential controls.** Freeze information access, starting/training distributions,
budgets, target definitions and easy/hard decision rules before classifying regimes.
Use structure-aware and generic learners, plus direct-search baselines, so one
model's misspecification is not mistaken for inherent unlearnability. Include the
additive limit and noise-free latent/observable comparisons. A positive monotone
sigmoid preserves exact order and zero-threshold uphill topology; finite resolution,
noise and saturation can change measured difficulty. Near-constant readouts must
not count as successful recovery merely because raw MSE is small, or as rugged
landscapes merely because every tied point is a weak local maximum.

**Falsifiers/limits.** The declared family/tasks fail to populate one or more
quadrants reproducibly; apparent dissociation disappears with matched observation
and budget conditions; or it is explained entirely by degeneracy, information
leakage or a learner unable to represent the target. Exact recovery plus unlimited
computation permits exhaustive search, so easy recovery/hard search requires a
computational or navigation restriction. Good average regression is not uniform
accuracy at rare high-fitness variants. Fourier recovery results for GNK do not
automatically apply to low-rank parameter count or the sigmoid-transformed oracle.

**Maturity and biological requirement.** The current sweep measures landscape and
walk descriptors, not recovery curves or ML optimization. Synthetic coverage would
establish benchmark behavior only. Relating its quadrants to real protein problems
requires the same operational tasks on held-out empirical landscapes; fresh
synthetic seeds are not independent biological validation.

**Evidence status.** Motivated by RES-063; no experiment currently tests HYP-011. No support is inferred from its registration.



## HYP-012: Biological-Information-Explains-Conditional-Learning-And-Search-Transitions

- Status: active
- Owner: matthew-spence
- Aims: AIM-001
- Motivating result: RES-063

The alignment of available biological information with the target
fitness landscape helps explain when function recovery and search are jointly
easy, jointly hard, or discordant. In the current model, correctly associated
family-profile information will produce reproducible, task-dependent benefits for
information-using methods over information-destroyed controls beyond effects
explained by generic landscape scale and navigation descriptors. A stronger claim
that this explains biological phenomena requires corresponding predictive evidence
on held-out empirical systems.

**Predictions and measurements.** Benefits depend on whether the information
predicts broadly sampled fitness, accessible beneficial changes, or rare target
variants; improvement in one need not imply improvement in the others. Compare
the effect of the information intervention on both recovery and search, including
interaction with method class and query distribution. Claims concern marginal
preferences in the present implementation, not empirical interaction identities,
contact structure or coevolutionary couplings.

**Essential controls.** Separate biological information used to construct the oracle
from information supplied to a learner. The exact generating profile reveals the
additive component by construction; its practical utility remains task- and
regime-dependent. This is not evidence that natural profiles predict real assay
fitness. Include no-information and shuffled
information controls and matched generic learners. A consistent amino-acid label
permutation of the entire oracle can preserve its fitness graph while disrupting
biological symbol meaning; simple profile shuffles need not preserve difficulty
and require explicit descriptor checks. PLM comparisons require full sequence
context and verified site mapping, with pretraining/assay overlap and target leakage
documented. Four-residue labels alone do not reproduce a natural protein context.

**Falsifiers/limits.** The information intervention has no reproducible task-specific
effect; generic methods change equally because the landscape became easier; the
apparent biological-specific benefit persists unchanged after destroying the
relevant information; or predicted empirical
patterns fail on held-out real data. Synthetic effects may support a controllable
information benchmark even when the stronger biological explanation fails. Avoid
claiming that absence of one marginal profile signal implies absence of all
biological information.

**Maturity and biological requirement.** ANA-070 establishes that the restored prior
changes aggregate behavior, with mixed trade-offs. No learner, optimizer or PLM
benefit has been measured. Any statement about why real systems are learnable or
searchable remains conditional until the predicted relationship is tested using
real fitness measurements not used to select the explanation. Statistical
prediction across systems alone would not establish a causal biological mechanism.

**Additional normalization and invariance constraints.** A common positive
profile temperature applied to final positive probabilities gives
M' = tau M + per-site constants. Reference centering and mutation-RMS
normalization cancel it; it is not an effective conservation-strength knob.
Site-specific concentration or smoothing can change relative preference geometry.
Whole-oracle relabelling preserves intrinsic learnability for a permutation-closed
hypothesis class; it tests compatibility with fixed biological priors. Transform
reference, starts, test sequences, masks, static noise assignments and proposal
rules consistently. Keep biological side information fixed in the deliberately
misaligned arm; a mapped-information rescue is a separate control. Generic
comparators must be label-equivariant. Include the direct generating-profile
score before assigning additional benefit to a PLM. The relabelled arm is an
information control, not another biologically plausible protein family.

**Evidence status.** Motivated by RES-063; no experiment currently tests HYP-012. No support is inferred from its registration.



### Literature for HYP-010 through HYP-012

- Brookes, Aghazadeh & Listgarten (2022), [On the sparsity of fitness functions and
  implications for learning](https://doi.org/10.1073/pnas.2109649118): GNK
  neighborhoods, Fourier sparsity and sparse-recovery sample complexity under
  explicit sampling/noise assumptions; not a guarantee for arbitrary learners or
  the present transformed low-rank oracle.
- Kaznatcheev (2019), [Computational Complexity as an Ultimate Constraint on
  Evolution](https://doi.org/10.1534/genetics.119.302000): compact NK landscapes can
  have difficult local optimization; worst-case results do not show that typical
  current generator draws are hard.
- Schmiegelt & Krug (2014), [Evolutionary accessibility of modular fitness
  landscapes](https://doi.org/10.1007/s10955-013-0868-8): interaction architecture
  affects accessibility and path distributions; modularity is not universally an
  easy-search condition.
- Stanton et al. (2024), [Closed-Form Test Functions for Biophysical Sequence
  Optimization Algorithms](https://arxiv.org/html/2407.00236v1): Ehrlich functions
  provide computationally inexpensive, tunable search challenges with attainable
  known optima. Difficulty is assessed relative to an optimizer, not a universal
  scalar knob.
- Greenman, Amini & Yang (2025), [Benchmarking uncertainty quantification for protein
  engineering](https://doi.org/10.1371/journal.pcbi.1012639): directly distinguishes
  model-improvement acquisition from top-sequence discovery; uncertainty metrics
  and acquisition outcomes are not interchangeable, and uncertainty-aware search
  did not outperform the greedy baseline in their evaluated settings.
