# CLAIM-005: AdaptEnv adapts genetic-architecture inference into a scalable generative fitness environment

- Status: completed
- Created: 2026-11-26
- Owner: matthew-spence
- Aim: AIM-001
- Relevant hypotheses: HYP-004, HYP-005, HYP-010, HYP-011, HYP-012
- Evidence: RES-064; EXP-065 / RUN-186, RUN-188, RUN-190 / ANA-072, ANA-073, ANA-074; pinned implementation and its verification.
- Context and boundaries: RES-063. Its ranked-backend empirical results do not validate the current full-pairwise backend.
- Supersedes: none. CLAIM-004 remains a historical claim about the ranked-GNK model and RUN-175.

## Claim

**AdaptEnv adapts the additive-plus-pairwise, nonlinear sequence–function
architecture used by Park, Metzger and Thornton into a profile-conditioned,
generative fitness environment. It provides explicit controls over latent
interaction strength and structure, separates intrinsic fitness from an optional
observation layer, and evaluates arbitrary genotypes without enumerating the
combinatorial sequence space. At fixed alphabet and interaction degree, model
storage and full-sequence evaluation scale linearly with sequence length, while
the represented space grows exponentially. This permits computational studies
across widely separated synthetic sequences; biological realism across those
evolutionary distances remains unvalidated.**

## Published foundation and generative contribution

Park, Metzger and Thornton (2024), *The simplicity of protein sequence-function
relationships*, infer reference-free genetic effects with a nonlinear sigmoid
readout and regularized regression. AdaptEnv adopts their additive-plus-pairwise
functional architecture, restricted to a selected graph of site pairs.
[Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11390738/),
[DOI](https://doi.org/10.1038/s41467-024-51895-5).

The generative choices are explicit: additive preferences come from a supplied
alignment/profile, interactions are sampled as full centered 20×20 tables, and
their total strength is normalized to an interpretable mutation scale. Sparse
supports, coefficient distributions and normalization specify our ensemble.
They are not empirically inferred interaction identities. Profile frequencies
supply marginal evolutionary information, not measured folding energies or
inferred residue couplings.

“Turning inference into generation” means using the published functional form
to construct controllable synthetic landscapes. It does not mean sampling a
posterior from their inference procedure. Our selected Optuna calibration matches
fitness quantiles, local SD and floor mass; it is not their regression objective.
The separate auxiliary-inference fit in RUN-184 did not achieve an adequate TrpB
surrogate (ANA-071).

This is a methodological adaptation and implementation claim, not a claim to
invent interaction expansions, nonlinear epistasis or procedural sequence
benchmarks. For example, Stanton et al. already propose closed-form Ehrlich
functions for biophysical sequence optimization. No priority or performance
superiority over that or other generators is asserted.
[Stanton et al. (2024)](https://arxiv.org/abs/2407.00236).

## What the architecture controls

For a sequence x, write the intrinsic landscape as

\[
F(x)=b+s\,\operatorname{sigmoid}\{\mu+\beta[a(x)+r g(x)]\},\qquad s,\beta>0.
\]

Here a and g are reference-centered additive and pairwise fields, each normalized
to unit single-substitution RMS under the uniform q20 background and uniform
nonself mutation operator. Therefore r is exactly the latent interaction/additive
mutation-RMS ratio. It is not the fraction of observed epistatic variance.

The default is pairwise-only, with full pair tables: no low-rank restriction
within a selected pair. The principal controls are:

| Knob and public API | Generic default | What it changes, and how to interpret it |
|---|---|---|
| **Interaction strength r** — `interaction_ratio` | 4.0 | Exact latent interaction/additive mutation-RMS ratio. Changes background dependence and can change mutation signs, peaks and paths. r=0 removes latent pair effects, not nonlinear or observation-induced epistasis. |
| **Interaction connectivity** — `target_degrees=(d,)` | d=6 | Number of interacting partners per site, capped by available sites. Redistributes interactions across a sparse graph at fixed total mutation RMS; it is not simply another amplitude control. Four-site TrpB already includes all six pairs. |
| **Effect concentration alpha** — `interaction_concentration` | `None` (Gaussian) | Finite positive alpha samples raw coefficients as sqrt(V)Z, with V~Gamma(alpha, rate=alpha). Smaller alpha concentrates variance in fewer large raw coefficients before centering and RMS normalization. It does not add higher interaction orders or change the meaning of r. |
| **Sigmoid steepness beta** — `steepness` | 0.5 | Controls how sharply latent variation crosses the transition and how much of the map saturates. Changes magnitudes, spectral epistasis and resolution-dependent neutrality; a positive monotone readout preserves exact intrinsic ranks. |
| **Sigmoid position** — `reference_occupancy`, or optional `global_intercept` | Occupancy 0.5; intercept `None` | Moves the latent distribution relative to the sigmoid transition, changing floor/ceiling occupancy and which branch shows diminishing or increasing returns. Occupancy anchors the reference genotype; `global_intercept` instead sets the uniform-background mean logit and overrides occupancy. These are alternative coordinate conventions, not two independent fitted controls. |
| **Biological additive preferences** — supplied `alignment_profile` / `profile_prior` | Required profile input | Changes relative site/residue preferences, conservation heterogeneity and the additive direction of improvement. Profile choice and smoothing are experimental conditions, not identified biochemical energies. |
| **Output bounds** — `fitness_min`, `fitness_max` | 0, 1 | Set assay offset and positive span. They change units and effects relative to fixed thresholds/noise, but not the exact ordering of intrinsic fitness. |
| **Observation law and draw reuse** — `observe_replicates`, `observation_bank`, supplied donor pools/uniforms | Intrinsic oracle has no noise | Independent versus shared errors and fixed versus newly drawn observations change apparent epistasis and repeatability. The TrpB shared-stop law is a separate explicit layer; it is not switched on by the generic model constructor. |

**Which knobs were fitted?** RUN-186/188 jointly varied continuous r in [0,8]
and four readout coordinates: output offset, log span, sigmoid position and log
steepness (r was fixed to zero in the ablation arm). RUN-190 additionally fitted
alpha in [0.05,100] on a log scale. The profile, pair graph design and observation
law were held fixed; noise amplitude was not optimized. The achieved full TrpB
fit selected r=0.5, which differs from the generic default r=4. Gamma selected
alpha=100 without a clear improvement. Generic defaults are not the fitted TrpB
parameter set.

There is **no independent global conservation-strength knob** in the current
normalized model: multiplying every log-profile preference by the same positive
constant cancels under mutation-RMS normalization. Changing relative site
conservation, smoothing or profile composition can change the landscape.
`length` (default 100) changes the dimension; `seed` (default 0) selects a
realization. They are design/replication settings, not calibrated ruggedness
parameters. `reference_rank` selects a computational anchor and can affect
positioning under the reference-occupancy convention. Batch size changes
computation, not the mathematical map. Latent interaction order is fixed at two
in this claimed default; no higher-order or low-rank control is implied.

These controls define the latent architecture; they do not provide independent
or universally monotonic control of every observable ruggedness statistic.
The sigmoid itself can introduce magnitude and higher-order epistasis. A strictly
monotone sigmoid alone cannot create sign epistasis from an additive latent map.
In exact arithmetic it preserves the ordering of a fixed latent map.

## Observation also contributes apparent epistasis

The tested observation model has

\[
Y_j(x)=F(x)+C(x)+e_j(x),
\]

with shared C and replicate-specific e fixed per genotype in the recorded maps.
These terms add epistasis outside the latent r-controlled architecture. Requerying
the same fixed error does not average it away; shared error survives replicate
averaging. The public API makes the draws explicit rather than silently identifying
them with latent interactions.

RES-064 provides the direct evidence. In RUN-188, the r=0 arm has no intrinsic
sign epistasis but an observed median sign fraction of 62.5%. For the full arm,
the corresponding fractions are 9.9% intrinsic and 64.1% observed. Yet intrinsic
after-sigmoid fitness supplies 91.1% of pooled order-2 power in that full arm:
noise does not explain all synthetic epistasis. These are distinct metrics and
must not be treated as interchangeable noise fractions.

Observation-induced **apparent** epistasis is also a plausible contributor in
real experimental systems. Park et al. demonstrate measurement-error inflation
of inferred interactions and discuss biological or technical phenotype bounds.
This supports the concern, not the identification of our specific shared-stop
error law in nature. [Park et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11390738/).

For current TrpB data, shared genotype-level error and biological fitness are not
uniquely separable without added assumptions: replacing F by F+h and C by C−h
leaves both replicate observations unchanged. This does not imply that all noise
is unidentifiable, or that biological epistasis is mostly error. Statistical
methods can distinguish interactions under specified noise assumptions; for
example, [Carlson, Andrews and Simons (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12501145/)
use rank constraints to test specific interactions in noisy combinatorial data.

## Tractability and evolutionary-scale sequence space

Let L be the number of variable sites, q=20, P the number of selected pairs,
and **n=q^L the number of possible sequences**. The implemented intrinsic oracle
stores O(Lq+Pq²) values and evaluates one full sequence in O(L+P) operations.
With bounded degree, P=O(L): both are O(L), or **O(log_q n)** for fixed q and degree.
Evaluating m requested sequences costs O(mL); it does not require visiting the
other n−m genotypes. Arbitrary Hamming distances and starting locations are allowed.

For example, L=1000 and degree six gives 3000 pair tables: 1,200,000 float64
entries, or 9.6 MB of pair-table payload, representing a map on 20^1000 sequences.
This is a derived storage count, not a throughput or peak-memory measurement.
It excludes additive arrays, temporary allocations and runtime overhead. The
default support construction also performs sorting, giving O(L log L) graph
setup for this bounded-degree case; table generation is O(Pq²). Other graph
settings have separate construction/retry costs. All-pairs coupling makes P
quadratic in L, but still does not require sequence-space enumeration.

The costs above concern the intrinsic oracle. Offline calibration, profile
construction, graph analyses and storing sampled data have additional costs.
The observation API requires caller-managed genotype/draw identity: caching fixed
errors only for m queried genotypes need not enumerate n, but materializing a
complete observed map or all-genotype noise bank would reintroduce O(n) storage.
The HTTP service does not automatically install the research shared-stop layer.

“Evolutionary scales” here means computational access to widely diverged,
fixed-length aligned sequences in arbitrarily large **finite** combinatorial
spaces, subject to ordinary hardware limits. It does not mean unlimited sequence
length on fixed hardware, validated phylogenetic dynamics, indel evolution, or
proof that distant synthetic variants remain members of a viable protein family.
Cheap evaluation also does not guarantee easy learning or optimization. No
comparative speed/memory advantage is established; HYP-010 remains prospective.

## Evidence, implementation and scope

- [RES-064](../results/RES-064_fixed-observation-noise-epistasis-and-trpb-identifiability/result.md) synthesizes the recent full-pairwise TrpB fits, observation ablation and Gamma comparison, with links to reproducible ANA notebooks and figures. The achieved full-model loss is 0.03099517 versus 0.04207450 for r=0. Gamma alpha=100 gives 0.03103213, with no clear advantage over Gaussian. Fixed-budget searches do not establish global optima.
- [RES-063](../results/RES-063_profile-model-readiness-and-learning-search-research-boundary/result.md) supplies the research boundary: comparative efficiency, independently controlled learning/search regimes and biological transfer remain unproven.
- Implementation: Adapt-Env `0.7.0rc1`, `thermodynamic_pairwise_v3`, commit `a6b5bfc5a86efa4f1a948d1b4ca968757b3fb343`. [Full pair backend](../external/Adapt-Env/adaptenv/full_pairwise.py), [model definition](../external/Adapt-Env/docs/THERMODYNAMIC_V2.md), [public fitting and observation methods](../external/Adapt-Env/docs/CALIBRATION.md).
- [Package verification](../results/RES-064_fixed-observation-noise-epistasis-and-trpb-identifiability/package_verification.json) records 234 passing tests in Python 3.13 and installed-package replay of the archived TrpB loss and observation calculations. Historical fingerprint portability limits are retained there. These checks support implementation correctness, not biological generalization.

Key figures: [sign epistasis with and without observation](../analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/figures/01_calibration_and_sign_epistasis.pdf),
[exact power accounting](../analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/figures/03_exact_power_accounting.pdf), and
[Gamma versus Gaussian latent concentration](../analyses/2026-11-25_ANA-074_run-190-gamma-concentration-versus-gaussian-trpb-calibration/figures/06_latent_concentration.pdf).
Their reproducible notebooks and frozen-data provenance are linked through RES-064.

The current empirical calibration supports conditional four-site TrpB aggregate
resemblance, with documented residual mismatches and observation assumptions.
This claim does not promote a calibrated preset, release the package, establish
learning/search performance, or extend that empirical validation to long proteins.

## Concise manuscript wording

“We adapt the additive–pairwise nonlinear architecture used for genetic inference
by Park, Metzger and Thornton into a profile-conditioned generative fitness
environment. Explicit latent interaction controls and a separate observation
layer permit controlled synthetic experiments, while bounded-degree construction
supports on-demand evaluation across exponentially large sequence spaces without
enumeration. TrpB calibration characterizes a conditional empirical regime;
biological transfer across evolutionary distances remains to be established.”
