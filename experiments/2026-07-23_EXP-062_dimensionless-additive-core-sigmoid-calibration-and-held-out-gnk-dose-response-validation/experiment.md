# EXP-062: Dimensionless additive-core sigmoid calibration and held-out GNK dose-response validation

## Status

- completed and analysed (`RUN-153`, `ANA-060`, `RES-060`)

## Scientific lineage and scope

- Aim: AIM-001
- Hypothesis tested directly: HYP-005
- Previous experiment: EXP-061
- Motivating results: RES-058 and RES-059

RES-058 established that the extensive candidate has a length-stable latent GNK
mutation scale. RES-059 then showed that GNK texture crosses the sigmoid when the
additive core is frozen, but not with the actual additive core: an ordinary
single mutation moves the additive score by about 55 public logits while the
sigmoid's 5%--95% dynamic band is only about 5.9 logits wide. EXP-062 tests the
narrow repair implied by that diagnosis. It calibrates the additive coordinate in
mutation units, fixes the reference operating point, and then asks whether GNK
texture is transmitted monotonically on held-out real-protein sequence contexts.

The candidate and calibrated readout are experiment-local. This experiment does
not change or validate the public model, load empirical assay outcomes, establish
biological realism, test accessibility or optimizers, or lift the Stage 4C/4D/4E
hard stop.

## Model and exact reparameterization

The primary readout is

`F(A,G) = C sigmoid(mu + beta (A-A_ref)/sigma_A + a G/sigma_G)`,

with `C=10`. Here `sigma_A` is one globally frozen calibration-panel estimate of
the complete reference-shell single-mutation RMS, and `sigma_G` is the candidate
field's analytic uniform-q20 nonself single-mutation RMS. The raw form is exactly

`k = beta/sigma_A`, `theta = A_ref - mu/k`, and
`w = a/(k sigma_G)`,

so that `k(A+wG-theta)` equals the dimensionless logit. The runner gates direct
float64 values and records the stable analytic secant only as a cancellation
diagnostic. It explicitly calls the inherited stable-residual function with
steepness one because its inputs are already logits.

The primary calibration is global: `sigma_A`, `beta`, and `mu` are frozen before
any confirmation GNK score is evaluated. A diagnostic upper-bound arm instead
sets `sigma_A` to each held-out core's own complete-h1 RMS. A rescue confined to
that arm means per-landscape normalization is required; it is not counted as a
global calibration success.

## Phase 1: additive-core-only calibration

The fixed panel has six length strata. In each stratum, the pre-existing EXP-060
reference plus one outcome-blind fixed real-protein alternative are assigned to
calibration. Two other real-protein sequences are held out for confirmation. The
short stratum spans lengths 47--52 because the source panel contains no valid
independent exact-length-48 Meltome entry. Its two confirmation sequences are the
first two under a recorded salted-hash rank over the 25 eligible canonical-q20
length-47--52 Meltome sequences, after retaining the lexicographically smallest
protein identifier for each duplicated sequence, and are weakly related to the
calibration pair.
The longer alternatives are an outcome-blind curated fixed panel, not a probability
sample. Exact sequences, source identifiers, hashes, split assignments, and these
selection limitations are frozen in `reference_catalog.yaml`.

Five additive-core seeds (`6201`--`6205`) are crossed with the 12 calibration
proteins. Their complete-h1 additive effect RMS values are pooled by one global
median. On those cells only, the experiment evaluates

- `beta = 0.5, 1, 1.5, 2, 2.25, 2.5, 3, 4, 8`;
- `mu = logit(0.95), logit(0.99)`.

A cell must have reference fitness at least 9.5, at least 50% of complete-h1
targets in occupancy `[0.05,0.95]`, h1 core-effect RMS at least 1.0 output unit,
at least 50% of `0.1L` reference-shell targets on the low plateau, and at least
90% of `0.5L` targets on the low plateau. A parameter cell is eligible only if at
least 80% of cells pass separately in every length stratum. The largest eligible
`beta`, then largest `mu`, is selected. No GNK quantity exists in the selection
table. If nothing qualifies, GNK confirmation is skipped and the decision is
`NO_CORE_OPERATING_POINT_QUALIFIES`.

These are mechanistic engineering guards, not empirical protein-fitness
calibration. The high reference, dynamic h1 shell, and extensive-distance low
plateaus also rule out solving transmission by replacing the sigmoid with a
nearly affine map.

## Phase 2: frozen held-out confirmation

The 12 held-out proteins are crossed with five new additive-core seeds
(`6211`--`6215`) and six independent GNK seeds (`6221`--`6226`), for 360
confirmation fields. Whole sequences and every seed namespace are disjoint from
calibration. Protein is the primary generalization unit; additive-core and GNK
seeds measure nested/crossed procedural variation. This remains fixed-panel,
not population-level, protein generalization.

The frozen global operating point must first reproduce the phase-1 core gates in
at least 80% of held-out cells in each length stratum. Only then is the selected
candidate tested at `a = 0, 0.25, 0.5, 1, 2, 4`. Directionality is evaluated over
`0.25`--`2`; `a=4` is a saturation stress and cannot fail monotonicity merely by
flattening.

The paired residual is

`R_a(x) = F(A(x), aG(x)) - F(A(x), 0)`.

All `19L` reference-shell substitutions are scored. Treatment-independent
reference-shell panels at `0.1L` and `0.5L` characterize extensive-distance
saturation but are not required to retain latent square-root growth through a
bounded output.

For two-site plaquettes, the primary texture contrast is evaluated as the direct
secant

`tau = R(G_ab) - R(G_a + G_b - G_0)`

at the double-mutant endpoint. Its predictor is
`a(G_ab-G_a-G_b+G_0)/sigma_G`. Candidate-supported order-2/3/4 site pairs and
uniform site-pair controls are evaluated on reference and uniform-q20 backgrounds.
Primary inference pools the reference-background supported strata within each
protein/core/GNK cell.

Required controls and diagnostic reports are a frozen additive core, the candidate-additivized double
mutant, a mutation-RMS-matched order-1 additive surrogate, exact `a=0`, direct
versus stable residuals, dimensionless/raw formula parity, steepness-only with the
margin-8 threshold, threshold-only with `k=20`, and candidate fields under raw
margin-8 and constructor-margin-10 readouts with `w=0.6`.
The steepness-only, threshold-only, and raw-margin ablations report edge transfer;
they do not by themselves establish which component is sufficient to rescue GNK
texture. The scientific decision concerns the combined calibrated operator.

## Preregistered confirmation gates

At `a=1`, each actual-core cell requires complete-h1 practical residual fraction
at least 0.25 and latent/residual correlation at least 0.50. Pooled supported
texture requires practical fraction at least 0.25, positive projection,
correlation at least 0.50, one-sided permutation `p<=0.05` with predictors
shuffled within support order inside each cell, and
actual/frozen texture RMS ratio at least 0.50. At least 80% of cells must pass
separately in each length stratum, and the hierarchical bootstrap lower bound for
texture projection must exceed zero.

Texture RMS and practical-fraction curves over `a=0.25,0.5,1,2` are computed
both per protein/core/GNK cell and after taking protein-level medians. Each curve
must have Spearman correlation at least 0.80, with preregistered growth floors of
0.02 output unit and 0.05 fraction. At least 80% of cells and proteins must pass
separately in every length stratum, at least 80% of proteins must pass overall,
and protein-bootstrap lower bounds for the `a=2` minus `a=0.25` changes must be
positive.

At `a=1`, the protein-level median texture RMS is regressed on log length. The
absolute log-log slope must be at most 0.15 and its stratified protein-bootstrap
interval must lie inside `[-0.25,0.25]`. The 0.10 practical floor and all related
bounds are preregistered noiseless mechanistic heuristics, not empirical or
service detection limits.

Decision precedence is:

1. `QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION`;
2. `EXECUTION_OR_PAIRING_INVALID`;
3. `NO_CORE_OPERATING_POINT_QUALIFIES`;
4. `CORE_CALIBRATION_DOES_NOT_GENERALIZE`;
5. `CALIBRATED_GNK_TRANSMISSION_NOT_DIRECTIONAL`;
6. `CALIBRATED_GNK_TRANSMISSION_NOT_LENGTH_STABLE`;
7. `GLOBAL_OPERATOR_CALIBRATED_GNK_TRANSMISSION_VALIDATED_NOISELESS`;
8. `PER_LANDSCAPE_CORE_NORMALIZATION_REQUIRED` if only exact-beta rescues all gates;
9. `GNK_TEXTURE_REMAINS_FROZEN_CORE_ONLY` if exact-beta also fails but the frozen
   positive control passes;
10. `CALIBRATED_SIGMOID_TRANSMISSION_FAILED` otherwise.

## Reproducibility

- Configuration: `config.yaml`
- Frozen panel: `reference_catalog.yaml`
- Entrypoint: `run.sh`
- Runner: `scripts/dimensionless_sigmoid_calibration_exp062.py`
- Candidate field: `src/project_adapt_env/extensive_gnk.py`
- Default target: `lab-slurm` (8 CPU, 64 GB RAM, 24 hours)

Every pair and plaquette design is treatment-independent and hash-replayed. Source
hashes cover the runner, inherited EXP-060/061 helpers, candidate field, and public
landscape configuration. Task checkpoints are keyed by protein, core seed, GNK
seed, effective configuration hash, and implementation hash. Quick mode is an
execution smoke test only.

The interpretation follows the distinction between latent epistasis and global
epistasis induced by a nonlinear phenotype-to-fitness map (Sailer & Harms,
*Genetics* 2017, doi:10.1534/genetics.116.195214; Otwinowski et al., *PNAS* 2018,
doi:10.1073/pnas.1804015115).

## Expected outputs

- `preflight.json`
- `reference_catalog.csv`
- `core_calibration_metrics.csv`
- `calibration_selection.csv`
- `validation_edge_metrics.csv`
- `validation_texture_metrics.csv`
- `dose_response_metrics.csv`
- `scale_law_fits.csv`
- `acceptance_checks.csv`
- `summary.json`

## Completion criteria

- [x] Configuration, panel, and protocol preregistered
- [x] Runner and focused unit tests implemented
- [x] Local quick smoke completed without a scientific decision
- [x] Dry-run `RUN-152` generated and inspected
- [x] Live `RUN-153` submitted to `lab-slurm` as Slurm job `2042`
- [x] `RUN-153` completed with scheduler exit code `0`
- [x] Outputs collected and DVC-tracked
- [x] ANA and RES records created after completion
- [x] HYP-005 and PROJECT_STATE.md updated after analysis

## Post-experiment notes

- `RUN-152` is the inspected source-staged dry-run record.
- `RUN-153` was accepted by `lab-slurm` as Slurm job `2042` and completed with
  exit code `0` on 2026-07-23.
- The collected evidence is DVC-tracked under
  `data/processed/dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-validation/RUN-153.dvc`
  (`373` files, `332,322,278` bytes,
  `02986be42c5b575ff8d7b75b3abad31c.dir`). All local files match the immutable
  shared-storage runtime by SHA-256, and the pointer has been pushed to the
  configured `labcluster` DVC remote.
- `ANA-060` reproduces the registered
  `CALIBRATED_SIGMOID_TRANSMISSION_FAILED` decision (`6/9` required checks).
  The global core calibration passes `60/60` held-out cells and complete-h1 edge
  transmission passes `360/360`, but texture practical prevalence passes only
  `24/360` actual-core and `0/360` frozen-core cells. Exact per-landscape
  normalization does not rescue texture (`23/360`).
- The formal global directionality failure is a floating-point equality
  artifact: treating `0.7999999999999999` as the registered mathematical
  boundary `0.80` with `1e-12` tolerance changes the independent count to `7/9`
  but leaves the overall decision unchanged.
- `RES-060` weakens HYP-005 at the full observable texture-prevalence boundary
  while preserving narrow held-out evidence for global core calibration,
  complete-h1 edge transfer, dose directionality, and registered length
  stability. The experiment-local/public-model scope limits and Stage 4C/4D/4E
  hard stop remain.
