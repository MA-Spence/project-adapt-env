# EXP-061: Operator-normalized GNK sigmoid-transmission paired-residual validation

## Status

- completed

## Scientific lineage and scope

- Aim: AIM-001
- Hypothesis tested directly: HYP-005
- Previous experiment: EXP-060
- Motivating results: RES-055, RES-057, RES-058

EXP-060 established the mutation-operator scale of the latent GNK field before
any nonlinear readout. EXP-061 asks the next, narrower question: does that field
remain observable after the public thermodynamic sigmoid, and does any loss arise
from the sigmoid operating point rather than a renewed latent length defect?

This experiment does not test empirical fitness realism, accessibility, adaptive
walks, or optimizer performance. It cannot by itself lift the Stage 4C/4D/4E hard
stop or replace the public AdaptEnv implementation.

## Readout and estimands

The frozen readout is

`F_theta(z) = 10 sigmoid(20 (z - theta))`.

For additive stability `A(x)`, reference-centred perturbation `G(x)`, and weight
`w`, the primitive paired residual is

`R_w(x) = F_theta(A(x) + w G(x)) - F_theta(A(x))`.

The mutation and two-site texture estimands are then

- `Delta R = R(y) - R(x) = Delta F_full - Delta F_core`;
- `epsilon_R = R_ij - R_i - R_j + R_0 = epsilon_F_full - epsilon_F_core`.

The normalized pointwise secant gain is

`gamma = R / [(10 * 20 / 4) wG]`,

with its continuous limit used at `wG = 0`. The runner also records the first-
order comparator `R_lin = F'(A)wG` and the corresponding edge and plaquette
contrasts. Observable epistasis alone is not attributed to GNK: the projection on
latent `w epsilon_G`, a within-stratum permutation null, and an order-1 additive
surrogate are required because a sigmoid can create global epistasis from an
additive perturbation.

Sign preservation by a strictly increasing sigmoid and exact rank preservation in
real arithmetic are mathematical consequences, not scientific validation gates.
They are retained only as implementation and floating-point audits.

Likewise, the genotype-wise occupancy assay is a deterministic calibration of this
known transfer function. Its bounds are technical validation; passing them is not
the scientific transmission claim. Scientific decision branches are determined by
the common-threshold reference-anchored and deployment landscapes.

## Three operating-point scopes

1. `pointwise_occupancy` sets `theta` separately for each genotype so the core
   occupancy is exactly `0.10`, `0.50`, `0.90`, or `0.99`. It isolates scalar
   transfer, secant gain, and finite-shift linearity. Because it is not one common
   landscape, it supports no topology claim.
2. `reference_anchored` sets one threshold per `(reference, seed, target
   occupancy)` and freezes it across genotypes, weights, and arms. This is the
   mechanistic landscape test for `Delta R` and `epsilon_R`. Both endpoint core
   occupancies are recorded and edges are stratified as both-dynamic, one-crossing,
   same-low, or same-high; anchoring the reference at 0.5 does not force its mutants
   to remain in the transition region.
3. `deployment` uses the constructed additive-core threshold. The documented
   `.config.yaml` base defaults and `natural_single_domain` template use margin
   `8`, whereas bare `LandscapeConfig()` currently defaults to margin `10` and
   named API presets span other margins. There is therefore no universal service
   or API margin. The margin-8 base-template formula is the preregistered primary,
   and the bare-constructor margin-10 formula is a separately labelled comparator.
   Both use steepness `20`, output scale `10`, and deployment weight `0.6`.

Thresholds remain fixed across `w` within every legitimate landscape scope.
Reference centring of `G` is also frozen: it is not a harmless gauge change after a
nonlinear readout.

## Fields and controls

The primary field is the exact EXP-060 `extensive_balanced_bounded` construction.
Paired on identical sequences are:

- the same public 256-term fixed-support field on its raw global scale;
- those identical public supports and contrasts locally mutation-RMS normalized
  and reference centred;
- a deterministic order-1, reference-centred additive surrogate matched to the
  candidate's mutation-operator RMS;
- a per-candidate-plaquette additivized counterfactual that keeps `G_0`, `G_a`,
  and `G_b` fixed but sets `G_ab^0 = G_a + G_b - G_0`;
- identity-readout and `w = 0` algebraic controls.

The additive core is built with GNK, legacy epistasis, noise, and peak lock disabled.
The candidate and controls are evaluated outside that construction. One candidate
field is generated per `(reference, seed)` and scalar multiples of that same field
are used at every weight.

The primary texture contrast is
`tau_texture = epsilon_R(actual G) - epsilon_R(additivized G)`. It changes only
the double-mutant candidate value and therefore isolates transmission of latent
`epsilon_G` from the global epistasis that the same sigmoid would create with an
additive perturbation. The independent order-1 field remains a distributional
control; its observable epistasis is reported rather than dismissed merely because
its latent epistasis is zero.

The frozen routing surface is `model_system = additive_thermodynamic_gnk`,
`observed_fitness_combine_mode = stability_gate`, and
`generic_epistasis_target = stability`. Observation noise is disabled to identify
the readout mechanism.

## Intervention and probe panel

Weights are parameterized by the dimensionless sigmoid displacement
`a = 20 w sqrt(E[(Delta G_candidate)^2])`. The grid is
`a = 0, 0.25, 0.5, 1, 2, 4`; each field's `w` is derived from its analytic mutation-
operator energy before observable values are inspected. Fixed public stress weights
`0.6` and `1.0` are added, and the same numeric `w` is applied to every paired arm.

The six hash-frozen EXP-060 references span lengths 48, 93, 118, 220, 238, and
448. Seeds `2601`-`2630` are independent generator replicates. The local probe
enumerates all `19L` substitutions from each reference rather than sampling them.
Each field also includes 512 single-step pairs on uniform-q20 backgrounds so the
local conclusion is not tied only to the selected natural parent.
Reference-shell and uniform-background pairs at `0.1L` and `0.5L` retain the
EXP-060 mutation rules with 512 pairs per distance; they characterize saturation
only, because a bounded sigmoid is not required to preserve latent `sqrt(L)` growth
at extensive distance.

Two-site plaquettes are sampled on both reference and uniform-q20 backgrounds,
treatment-independently from candidate co-supported site pairs at orders 2, 3, and
4 and from uniform site-pair controls. Each plaquette
is evaluated with the real additive core and with an explicitly labelled frozen-
core mechanistic control. The frozen-core arm asks whether the sigmoid can carry
GNK curvature at all; only the real-core arm tests survival in the actual landscape.

The output cross-product is deliberately pruned without changing any primary gate:
pointwise calibration is emitted only for one-step panels; extensive-distance
reference-anchored panels retain occupancy 0.5; fractional-distance weights retain
`a = 0, 0.5, 1` and fixed `w = 0.6, 1`; and deployment plaquettes retain
`a = 0` plus fixed `w = 0.6, 1`. Reference-anchored plaquettes retain the complete
operator-scaled grid at occupancy 0.5, with `a = 0.5` retained at the other
diagnostic occupancies.

## Preregistered gates and decisions

Technical gates require deterministic field replay, treatment-independent pair and
plaquette hashes, finite inputs, exact `R_0 = 0`, occupancy calibration, public-
formula parity where direct subtraction is resolvable, identity-readout algebra to
`1e-12`, and zero latent plaquette epistasis in the additive surrogate.

Controlled scalar transfer at occupancy 0.5 and `a = 0.5` requires per-length
median normalized gain at least `0.75`, operational detection in at least 90% of
eligible pointwise residuals above a `0.10` output-unit practical resolution floor
used as a preregistered heuristic in this noiseless assay—not as a service
detection limit—
and no renewed local log-length slope larger than `0.15` in magnitude (paired-seed
95% interval inside `[-0.25, 0.25]`). At occupancies 0.1 and 0.9, median gain must be
at least `0.20`; occupancy 0.99 is diagnostic. The finite-shift linear comparator
requires normalized RMSE at most `0.25` and correlation at least `0.90`.

Reference-anchored texture at occupancy 0.5 uses `a = 0.5` for projection and
finite-shift linearity, and `a = 1.0` for the preregistered practical-resolution
fraction. It requires candidate `tau_texture` to have positive through-origin latent projection
`sum(epsilon_G tau_texture) / sum(epsilon_G^2)` in at least 80% of seeds, a paired-
seed bootstrap lower bound above zero, and a one-sided within-stratum permutation
`p <= 0.05` in at least 80% of required within-stratum tests. It also requires
correlation at least `0.50` with its finite-secant/linear comparator and a
practical-resolution fraction of at least `0.25`. Every one of those bounds is
applied separately at each expected length; the frozen-core and real-core results
are also gated separately. Their `tau_texture` RMS log-length slopes must be at
most `0.15` in absolute value, with paired-seed 95% intervals inside
`[-0.25, 0.25]`.

A missing eligible stratum cannot silently pass: every expected length-by-seed
cell is reindexed explicitly, each length requires at least 80% seed coverage,
and at least 32 eligible observations per retained row. The deployment test at
margin 8 and `w = 0.6` separately requires, at every length, at least 25%
practically resolvable local residuals, positive texture projection in at least
80% of seeds, a positive per-length bootstrap lower bound, and at least 80% of
required permutation tests passing. Direct float64 and stable analytic residuals
are both recorded so mathematical attenuation is not confused with exact numerical
ties. Only direct base-template-formula values can satisfy deployment gates; stable
residuals are explanatory and cannot rescue a numerically tied public output.

The controlled and deployment plaquette gates use the reference-background
supported strata. Uniform-q20 square backgrounds and `0.1L`/`0.5L` edge panels
report direct/stable tie, gain, projection, and practical-resolution behavior as
saturation/generalization diagnostics; they are not forced to preserve latent
square-root growth through a bounded output.

The decision branches are:

- `EXECUTION_OR_PAIRING_INVALID` for a technical failure;
- `SIGMOID_TRANSFER_FAILED` if common-threshold reference-anchored transfer fails
  after the pointwise implementation calibration passes;
- `GNK_TEXTURE_ONLY_IN_FROZEN_CORE` if texture passes only without additive-core
  variation;
- `TRANSITION_REGION_GNK_TRANSMISSION_VALIDATED_PUBLIC_POINT_SATURATED` if the
  controlled real-core mechanism passes but the margin-8 deployment point fails;
- `CANDIDATE_GNK_MARGIN8_NOISELESS_PUBLIC_READOUT_TRANSMISSION_VALIDATED` only if
  controlled transfer, real-core texture, and the margin-8 deployment gates all
  pass;
- `SIGMOID_TRANSMISSION_NOT_ESTABLISHED` for other controlled texture failures.

The margin-10 constructor comparator is mandatory for execution and interpretation
but does not substitute for the preregistered margin-8 decision.

The candidate replaces the current fixed-support generator here and observation
noise is disabled. Even the strongest label therefore validates only candidate
transmission through the current readout/routing formula at the base-template margin,
not the current generator or service end to end. Also, EXP-060's raw legacy arm used
384 terms; EXP-061 deliberately uses the current base's 256 terms and does not treat
those raw controls as identical.

## Reproducibility and analysis

- Configuration: `config.yaml`
- Entrypoint: `run.sh`
- Runner: `scripts/extensive_gnk_sigmoid_transmission_exp061.py`
- Latent field: `src/project_adapt_env/extensive_gnk.py`
- Frozen references/configuration source: EXP-060 `config.yaml`
- Default target: `lab-slurm` (8 CPU, 32 GB RAM, 24 hours)

Seed is the inferential unit. Length is paired within seed, but there is only one
protein per length, so length remains confounded with protein identity. Quick mode
is an execution smoke test and cannot issue any scientific decision. Full outputs
include design and occupancy audits, paired residual metrics, edge and epistasis
transmission summaries, topology/numerical-tie audits, scale fits, preregistered
checks, and task checkpoints.

The implementation streams the readout cross-product rather than materializing it.
Seed-level output retains the sufficient sums `n, sum(x), sum(y), sum(x^2),
sum(y^2), sum(xy)` for projections/correlations, plus threshold, endpoint
occupancies, direct/stable residual summaries, pairing hashes, and observed/null
permutation statistics. Source and implementation hashes cover both runners, the
candidate field, the public landscape, presets, and `.config.yaml`.

Interpretation should be consistent with global-epistasis work showing that a
nonlinear phenotype-to-fitness map can generate observable epistasis without latent
interactions (Sailer & Harms, *Genetics* 2017, doi:10.1534/genetics.116.195214;
Otwinowski et al., *PNAS* 2018, doi:10.1073/pnas.1804015115) and with thresholded
protein fluorescence landscapes (Sarkisyan et al., *Nature* 2016,
doi:10.1038/nature17995).

## Completion criteria

- [x] Configuration and protocol preregistered
- [x] Run script and environment recorded
- [x] Local quick smoke and focused tests pass
- [x] Dry-run record generated and inspected
- [x] Live run submitted to `lab-slurm`
- [x] Outputs collected and DVC-tracked
- [x] ANA and RES records created after completion
- [x] HYP-005 and PROJECT_STATE.md updated after analysis

## Post-experiment notes

- `RUN-150` is the inspected `lab-slurm` dry-run record. Its staged script requests
  8 CPUs, 32 GB RAM, and 24 hours and contains the declared implementation inputs.
- Live `RUN-151` was accepted by `lab-slurm` as scheduler job `2041` on
  2026-07-23 and completed with exit code `0`.
- The ten outputs, three logs, and 180 checkpoints were collected and DVC-tracked
  as a 193-file, 1,513,486,733-byte bundle. All local files match the immutable
  remote runtime SHA-256 manifest. Because the configured client mount was
  unavailable, the DVC content-addressed objects were backed up to the configured
  remote store directly over SSH and verified with a checksum-aware dry run.
- `ANA-059` independently reproduces all 16 gate states. All 12 technical checks
  pass, but reference-anchored edge coverage, actual-core texture, and margin-8
  deployment fail; frozen-core texture passes. The registered decision is
  `SIGMOID_TRANSFER_FAILED`.
- `RES-059` weakens HYP-005 at the observable protein-scale boundary while
  preserving the RES-058 latent result and the positive frozen-core mechanistic
  control. The public generator/service remains unvalidated and the Stage
  4C/4D/4E hard stop remains.
