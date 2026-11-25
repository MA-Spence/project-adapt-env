# TrpB matched interaction ablation under frozen shared-stop observation noise

## Question and scope

Does the r-controlled pairwise component reproduce empirical interaction structure
that a refitted profile-additive sigmoid model with the **same observation model**
cannot reproduce? This is a follow-up to ANA-072/RUN-186 under EXP-065, AIM-001 and
HYP-004/HYP-005. It addresses observation dependence, not computational scaling,
empirical coefficient recovery, or actual learning/search performance. The scientist
authorized registration and Slurm execution on 25 November 2026.

## Paired arms

1. **r0:** profile additive field, sigmoid readout, r fixed exactly to zero.
2. **full:** identical construction with all six full 20×20 Gaussian pair tables;
   jointly fit continuous r in [0,8] and the four sigmoid/readout coordinates.

Both independently refit offset, scale, position and steepness with identical
RUN-180/RUN-186 bounds and objective: the original 101-quantile discrepancy,
local-cube SD discrepancy and near-floor-fraction discrepancy, using the original
weights and numerical implementation. The frozen RUN-186 source is imported;
all six historical shared-stop objective values must replay before fitting.
The model's sigmoid can itself induce statistical epistasis at r=0. Consequently,
this control is **not** an additive model on the final fitness scale.

Each arm receives 11,712 sequential Optuna 4.9.0 TPE trials, including queued
candidates: 23,424 trials total. Seed 58459495 and 20 startup trials are fixed.
The r0 arm starts from the six historical readout vectors, all reevaluated at r=0.
The full arm starts from the six historical (r,readout) vectors plus the newly
refitted r0 winner. This ensures the full model includes the achieved null fit;
the two arms have equal total budgets but six versus seven queued candidates.
Four- and five-dimensional TPE proposals are not identical despite the same
sampler seed. This is a fixed-budget comparison, not a convergence claim or a
comparison of optimisers. A tie in training loss is not an equivalence test.
Gamma concentration and noise strength are not fitted. The package/defaults and
historical RUN-186 archive remain unchanged.

## Frozen inputs and observation assumptions

Reuse RUN-186's exact calibration fields, genotype panels, noise pools, observation
uniforms, confirmation fields and seed audit. The 32 empirical calibration cubes
and 32 confirmation cubes have disjoint genotypes. Calibration uses the same 64
synthetic cubes. Confirmation uses the original 24 generator seeds (16 reused and
eight originally fresh), eight assigned cubes per seed. Primary synthetic summaries
should use the eight originally fresh seeds; all 24 remain available. These panels
and seeds have already been inspected: this is a diagnostic ablation, not new
independent empirical validation. Full-library graph descriptors include calibration
nodes and use the same 159,129-node empirical coverage mask.

The shared-stop observation law is fixed before either fit: independent errors
come from centered replicate-difference donors conditional on fitness and depth;
common errors come from centered stop-control means rescaled to the nonnegative
stop-control replicate covariance within depth bins. Common error is shared by the
two observations of one genotype. It is independent across genotypes under this
sampling assumption and is frozen by genotype within the observed map.

Generator fields, noise uniforms and observation law are paired between arms.
Actual numerical errors need not be identical: changing predicted fitness can
change the conditional donor/depth bin. No new noise-amplitude degree of freedom
is introduced. Original replicate differences do not identify shared assay error;
Johnston's input-count normalization is shared across replicates. Stop controls
are a sensitivity assumption, not proof of the error law for functional variants.

## Unfitted diagnostics

These descriptors do not enter calibration or model selection. Save individual
cube/face values and seed identifiers, not only aggregate summaries.

- **Absolute cross-replicate interaction power** at orders 2, 3 and 4. For each
  four-site binary cube compute Walsh coefficients c1 and c2 with the same 1/16
  normalization, then sum c1[S]c2[S] by order. Independent zero-mean errors cancel
  in expectation; common errors do not. Retain signed estimates. Report order-1
  power and ordinary mean-fitness powers as context, and pooled higher-order
  fractions secondarily. Near-zero denominators must not drive the conclusion.
- **Replicate-consistent sign epistasis:** the same mutation must reverse sign
  across the same backgrounds, with both directed edge signs agreeing in both
  replicates. Reciprocal sign epistasis requires agreement of all four edges.
  Primary effect resolution is 0.001 × the fixed calibration scale; retain 0
  and 0.01 sensitivity resolutions. Do not count two opposite-orientation
  reversals as reproducible.
- **Activity strata and layers:** report all, no-active, some-active and all-active
  cubes/faces. A vertex is active only when both replicates exceed their original
  replicate-specific thresholds. Cube and face strata have different membership;
  export both denominators and retain empty strata as missing, never zero effects.
  Activity selection is observation-conditioned, not a noise-free filter. Report
  intrinsic and observed layers separately for each arm.
- **Stop-covariance sensitivity:** subtract |order| × sum_v q_v / 256 from each
  cross-power sum, where |order| denotes the number of Walsh coefficients of that
  order. This assumes independent shared errors across vertices. Use raw stop
  covariance estimates here, including negative estimates, exactly as a diagnostic
  sensitivity; sampling uses max(0,q). Keep negative adjusted powers. Do not treat
  this correction as identified biological interaction power, especially after
  activity stratification.
- **Secondary static descriptors:** exact no-beneficial-neighbor fraction, edge
  RMS and neighbor correlation on the common full-library graph. These describe
  local geometry; no learner, adaptive walk or optimisation benchmark is run.

Compare arms at matched seeds and show the empirical values against their seed
variation. Uncertainty for an eventual analysis must preserve cube dependence of
faces and distinguish model-seed variation from empirical sampling uncertainty.
Do not pool faces as independent replicates or equate a nonsignificant difference
with equality. No new arbitrary numerical graduation threshold is introduced.

## Interpretation fixed before execution

- If the full arm materially improves agreement with repeatable interaction
  structure, particularly away from the assay floor, r contributes beyond the
  additive-sigmoid-plus-noise control. This still does not establish that the
  shared-stop law is correct or all repeatable structure is biological.
- If r0 reproduces the apparent ruggedness comparably, the current observations
  do not establish that r controls the empirically matched ruggedness. The sigmoid
  and/or observation model can supply it without latent pairwise interactions.
- If both miss repeatable structure, neither arm is adequate on that diagnostic.
  If power, activity strata and signs disagree, report the disagreement rather
  than reducing it to one favourable score.

This experiment can expose noise compensating for missing structure; it cannot
uniquely distinguish shared measurement bias from stable biological effects with
these two empirical replicates. Resampling synthetic noise alone cannot resolve
that identifiability problem.

## Literature and existing record

Johnston et al. (2024), [PNAS 121, e2400439121](https://doi.org/10.1073/pnas.2400439121),
tested an additive-plus-measurement-error explanation of their TrpB landscape.
Their conclusion does not automatically transfer to this control: our sigmoid
and shared-stop law are different. ANA-065 already assessed replicate-aware
Walsh powers and the instability of stop correction near the floor; ANA-072
motivates this newly refitted matched ablation. No novel inference method is
claimed for cross-replicate powers.

## Execution and reproducibility

`interaction_ablation.yaml` fixes settings. `scripts/interaction_ablation_exp065.py`
fits both arms and exports diagnostics; `interaction_ablation_metrics_exp065.py`
contains known-answer-tested statistics. A SHA-256 manifest authenticates the new
runner/configuration and the entire immutable RUN-186 archive in a dedicated DVC
runtime bundle. The worker verifies every payload and all 36 pinned Linux package
versions before execution. Six historical objectives and every reused confirmation
seed's RUN-186 values are replayed as integrity checks.

Use `submit_interaction_ablation_slurm.py --dry-run` before `--execute`.
The helper temporarily selects the explicit ablation handoff and Dutch SSH
transport, then restores the project dispatcher/transport settings. The scientific
worker uses the existing dedicated Python 3.11.15 Linux environment, one CPU,
8 GB RAM and a 24-hour limit. Only smoke validation is executed locally. The
production two-arm fit runs through labproj on the CPU Slurm queue.
