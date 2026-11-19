# EXP-063: Paired categorical-code and representation-rank GNK texture rescue

## Status

- completed; `RUN-157` finished as Slurm job `2044` with exit code `0`,
  analysed in `ANA-061`, and recorded in `RES-061`

## Scientific lineage and scope

- Aim: AIM-001
- Hypothesis tested directly: HYP-005
- Previous experiment: EXP-062
- Motivating result: RES-060

RES-060 established that the experiment-local extensive GNK field has stable
mutation-unit semantics through the calibrated additive core and sigmoid:
all `60/60` held-out core cells and `360/360` complete-h1 edge cells passed,
GNK-specific two-site texture was aligned and dose-directional, and its RMS was
length-stable. The formal failure was practical texture prevalence: only
`24/360` actual-core cells had at least 25% of supported plaquettes above the
preregistered `0.10` output threshold.

The current field assigns one balanced binary rank-one categorical product to
each order-2/3/4 support. A nonself q20 substitution changes one such binary
contrast with probability `10/19`; a singly supported pair changes both factors
with probability `(10/19)^2 = 100/361`. This motivates, but does not establish,
the post-RES-060 diagnosis that binary allele ties and/or rank-one factorization
limit practical prevalence.

EXP-063 is a single paired representation intervention. It changes only the
categorical allele code and the number of separable components on each fixed
support. It does not change the additive core, support topology, support degree,
per-order mutation energy, global calibration, sigmoid, amplitude, probe
genotypes, or acceptance thresholds.

The candidate remains experiment-local. No empirical assay outcome is loaded,
and this experiment cannot validate DCA inference, structural contacts,
biological realism, the public AdaptEnv generator, accessibility, adaptive
walks, optimizers, or active learning. The Stage 4C/4D/4E hard stop remains
under every outcome.

## Question

At fixed support topology, mutation-operator energy, additive core, and sigmoid
readout, do dense q20 allele coding and/or rank-8 categorical tensor components
rescue the practical prevalence of transmitted GNK texture?

## Pre-experiment prediction

- If binary allele ties are sufficient to explain RES-060, `dense_q20_r1` will
  pass the frozen EXP-062 observable gates while `binary_r1` reproduces the
  prevalence failure.
- If repeated categorical directions are sufficient, `binary_r8` will pass.
- If both changes are needed, only `dense_q20_r8` will pass.
- If a richer arm repairs latent categorical density but not observable
  prevalence, representation rank is insufficient at the fixed sigmoid.
- If `binary_r1` passes on the fresh field/core seeds, the RES-060 prevalence
  failure is procedurally unstable and no rescue claim is justified.

## Frozen panel and replication

The exact 12 EXP-062 confirmation proteins are reused as a fixed failure-
challenge panel: two proteins in each of six length strata, spanning
`L=49-448`. They are crossed with five fresh additive-core seeds
`6311-6315` and six fresh GNK seeds `6321-6326`, giving 360 paired
protein/core/GNK cells per arm.

The panel was already examined in RES-060. It provides high-power paired causal
diagnosis but is not a new held-out protein sample. Any selected representation
would still require later confirmation on a new frozen panel before public
integration.

No empirical sequence-associated outcome is used. Candidate arm is excluded
from support, background, mutation, and square sampling keys so every
comparison is genotype-paired.

## Frozen additive core and readout

The experiment uses the EXP-062 confirmation construction with no calibration
or arm-specific tuning:

`F(A,G) = 10 sigmoid(mu + beta (A-A_ref)/sigma_A + a G/sigma_G)`,

with:

- `sigma_A = 2.7448641782814067`;
- `beta = 4`;
- `mu = logit(0.99) = 4.595119850134589`;
- `a = 0, 0.25, 0.5, 1, 2, 4`;
- `a=4` retained as saturation stress, outside the monotonicity-failure rule.

The primary paired residual remains

`R_a(x) = F(A(x), aG(x)) - F(A(x), 0)`.

The two-site candidate-additivized texture contrast remains

`tau = R(G_ab) - R(G_a + G_b - G_0)`,

with latent predictor

`a (G_ab-G_a-G_b+G_0) / sigma_G`.

## Paired categorical representation arms

All arms use identical balanced simple support hypergraphs, orders `(2,3,4)`,
target per-site degrees `(6,11,7)`, order mutation-energy shares
`(0.35,0.45,0.20)`, support IDs, coefficient namespaces, and reference
centering.

1. `binary_r1`: the incumbent balanced `+/-1`, one-product representation.
2. `dense_q20_r1`: one bounded, centered q20 contrast with distinct values,
   removing binary within-class mutation ties while retaining rank one.
3. `binary_r8`: eight independently hashed balanced binary product components
   per support, changing rank while retaining binary component coding.
4. `dense_q20_r8`: eight independently hashed dense q20 product components per
   support, applying both interventions.

For a support `S_t`, a rank-R term has the form

`G_t(x) = sum_r c_tr product_{i in S_t} phi_tir(x_i)`.

The implementation computes each support's exact component mutation-operator
Gram matrix. One order-level scale enforces the same registered per-order and
total expected q20 nonself single-mutation energy in every arm. Rank therefore
changes categorical geometry, not nominal dose. Support degree counts unique
supports rather than product components, and field variance includes
same-support cross-component terms.

`binary_r1` must match the existing `ExtensiveGNKField` implementation under
the same inputs. Rank is fixed globally; it cannot vary with order, protein, or
length.

## Probe and controls

Every arm is evaluated on:

- the complete `19L` reference single-mutant shell;
- candidate-supported order-2, order-3, and order-4 site pairs;
- uniformly sampled pairs excluded from the union of candidate support pairs;
- reference and uniform-q20 square backgrounds;
- 128 treatment-independent squares per stratum/background;
- the same dose grid and paired predictor/response arrays.

Required controls are:

- incumbent `binary_r1`;
- exact `a=0`;
- a mutation-RMS-matched order-1 additive surrogate;
- candidate-additivized double mutants;
- dimensionless/raw formula parity;
- frozen additive logits at the selected 0.99 reference occupancy, retained as
  an upper-shoulder saturation diagnostic;
- frozen additive logits at occupancy 0.5, added as the valid transition-window
  texture positive control.

The 0.99 frozen arm is not required to pass a texture-prevalence gate because
RES-060 established that it is saturated by construction.

## Manipulation and technical checks

Before scientific interpretation:

- all support and square hashes must match across arms;
- every arm must be centered and finite;
- per-order mutation-operator energy must match target within `1e-10`
  relative error;
- `binary_r1` must replay the incumbent implementation;
- dense-q20 nonself categorical tie fraction must be zero;
- rank-8 component audits must show effective rank at least 4;
- zero-dose, additive-surrogate, and formula-parity controls must pass.

Complete-h1 mutation RMS, field variance, latent exact/near-zero fractions,
quantiles, kurtosis, component Gram eigenvalues, and participation ratios are
reported as diagnostics rather than separately tuned targets. Complete-h1 RMS
is deliberately not forced to match: exact uniform-q20 nonself
mutation-operator energy is the controlled quantity, while reference-shell
geometry is one consequence of changing the categorical representation.

## Preregistered observable gates

The EXP-062 actual-core thresholds are retained without outcome-dependent
changes. At `a=1`, each arm requires:

- complete-h1 practical residual fraction at least 0.25;
- supported texture practical fraction at least 0.25 under the fixed `0.10`
  output floor;
- latent/observable correlation at least 0.50;
- positive projection;
- within-order one-sided permutation `p <= 0.05`;
- at least 80% of cells passing separately in every length stratum;
- hierarchical protein-bootstrap lower bound for projection above zero.

The transition-window frozen occupancy-0.5 control must independently transmit
texture. The saturated occupancy-0.99 arm and an actual/saturated RMS ratio are
diagnostic only. The rank-8 or dense-code arm must increase the protein-median
practical fraction by at least 0.05 over `binary_r1` and have a positive paired
protein-bootstrap lower bound for that improvement.

Before these arm comparisons, the shared frozen EXP-062 core must reproduce
the original high-reference, dynamic-h1, h1-effect, and extensive-distance
plateau guards in at least 80% of core cells separately in every length
stratum.

For `a=0.25,0.5,1,2`, texture RMS and practical-fraction curves retain:

- Spearman correlation at least 0.80 using a `1e-12` comparison tolerance;
- RMS growth at least 0.02 and prevalence growth at least 0.05;
- at least 80% cell/protein passes overall and within every length stratum;
- positive protein-bootstrap lower bounds for `a=2` minus `a=0.25`.

At `a=1`, the protein-median texture RMS log-length slope must have absolute
value at most 0.15 and its bootstrap interval must lie inside
`[-0.25,0.25]`.

The `0.10` and `0.25` thresholds remain noiseless engineering criteria, not
empirical detection limits.

## Decision precedence

1. `QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION`;
2. `EXECUTION_PAIRING_OR_MANIPULATION_INVALID`;
3. `FROZEN_EXP062_CORE_OPERATING_POINT_NOT_REPRODUCED`;
4. `INCUMBENT_PREVALENCE_FAILURE_NOT_REPLICATED`;
5. `DENSE_Q20_RANK1_SUFFICIENT`;
6. `BINARY_RANK8_SUFFICIENT`;
7. `ALTERNATIVE_SINGLE_FACTOR_REPAIRS_BOTH_SUFFICIENT` when both single-factor
   arms pass;
8. `COMBINED_DENSE_Q20_RANK8_REQUIRED`;
9. `APPARENT_RESCUE_NOT_PAIRED_SUPERIOR_TO_BINARY_R1`;
10. `NO_INTERVENTION_PASSES_MIDPOINT_CONTROL`;
11. `NO_INTERVENTION_PASSES_EDGE_TRANSFER`;
12. `NO_INTERVENTION_RESCUES_TEXTURE_PREVALENCE`;
13. `REPRESENTATION_RESCUE_NOT_DIRECTIONAL`;
14. `REPRESENTATION_RESCUE_NOT_LENGTH_STABLE`;
15. `REPRESENTATION_REPAIR_INSUFFICIENT_AT_FIXED_SIGMOID`.

If multiple richer arms pass, report all passing arms and prefer the
lower-component representation; no parameter is selected by optimizing a
continuous outcome.

## Reproducibility

- Configuration: `config.yaml`
- Frozen panel: EXP-062 `reference_catalog.yaml`
- Entrypoint: `run.sh`
- Runner: `scripts/categorical_representation_rank_rescue_exp063.py`
- Ranked field: `src/project_adapt_env/ranked_categorical_gnk.py`
- Default target: `lab-slurm` (8 CPU, 64 GB RAM, 24 hours)

All implementation hashes, arm-independent design hashes, exact operator-
energy audits, and task checkpoints are written before acceptance decisions.
Quick mode is an execution smoke test only.

The categorical comparison is motivated by generalized NK lookup geometry
(Brookes, Aghazadeh & Listgarten, *PNAS* 2022,
doi:10.1073/pnas.2109649118) and Potts/DCA amino-acid-specific couplings
(Ekeberg et al., *Physical Review E* 2013,
doi:10.1103/PhysRevE.87.012707). The tested random categorical arms are not
Boltzmann DCA models and do not infer couplings from an alignment.

## Expected outputs

- `preflight.json`
- `reference_catalog.csv`
- `core_guardrail_metrics.csv`
- `representation_catalog.csv`
- `representation_audit.csv`
- `latent_geometry_metrics.csv`
- `validation_edge_metrics.csv`
- `validation_texture_metrics.csv`
- `paired_arm_effects.csv`
- `factorial_effects.csv`
- `dose_response_metrics.csv`
- `scale_law_fits.csv`
- `acceptance_checks.csv`
- `decision_table.csv`
- `summary.json`

## Completion criteria

- [x] Scientific question, intervention, controls, gates, and terminal decisions preregistered
- [x] Ranked categorical field and runner implemented
- [x] Focused tests pass
- [x] Local quick smoke completes without a scientific decision
- [x] Dry-run submission generated and inspected
- [x] Live run submitted to `lab-slurm`
- [ ] Outputs collected and DVC-tracked after completion
- [ ] ANA and RES records created after completion

## Post-experiment notes

- The experiment reuses the RES-060 protein panel for a paired causal challenge
  but uses fresh core and GNK seeds. It is not new held-out protein validation.
- `RUN-154` is the initial 16-CPU dry run. `RUN-155` reached the submit host but
  Slurm rejected it before job creation because no node could satisfy 16 CPUs.
- The operational request was reduced, without changing the scientific design,
  to the previously proven 8-CPU, 64-GB, 24-hour envelope. Revised dry run
  `RUN-156` was inspected, and live `RUN-157` was accepted as Slurm job `2044`
  on `ripley`.
