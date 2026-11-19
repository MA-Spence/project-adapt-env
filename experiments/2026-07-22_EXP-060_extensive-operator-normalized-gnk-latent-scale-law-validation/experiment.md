# EXP-060: Extensive operator-normalized GNK latent scale-law validation

## Status

- completed

## Scientific lineage

- Aim: AIM-001
- Hypothesis tested directly: HYP-005
- Previous experiment: EXP-059
- Motivating results: RES-049, RES-055, RES-057

HYP-004 and HYP-009 remain downstream. This experiment tests neither empirical
landscape realism nor optimizer behavior.

## Question

Does an operator-normalized, extensive, balanced-support, bounded-contrast GNK
latent field preserve a controllable ruggedness scale from single substitutions
through full-protein sequence differences before any nonlinear readout is applied?

## Pre-experiment prediction

At working weight one, the candidate will maintain
`RMS(delta GNK) / RMS(delta additive) ~= 1` for the declared background mutation
operator without a length trend; its field variance will scale as `L`; and its GNK
difference RMS will scale as `sqrt(L)` when Hamming distance is a fixed fraction of
protein length. The current public fixed-term GNK control is expected to reproduce
an approximately `L^-1/2` local-effect defect.

## Scope and stopping boundary

The experiment ends at the additive stability score plus the centered GNK latent
field. The additive core is constructed with `gnk_terms = 0`; the current GNK
control is compiled only afterward for direct latent evaluation. No GNK score is
sent through AdaptEnv's sigmoid, and no observable-fitness, peak, accessibility,
walk, routing, or optimizer metric is calculated.

The candidate implementation remains local to this project. Passing EXP-060 would
qualify the latent field for a later sigmoid-transmission experiment; it would not
silently replace the public AdaptEnv implementation or establish biological realism.

## Candidate construction

For additive score `A(x) = sum_i a_i(x_i)`, the candidate field is

`G(x) = sum_r sum_{t in T_r} c_t product_{i in S_t} phi_ti(x_i)`.

- Active orders are 2, 3, and 4 with target mutation-energy shares
  `(0.35, 0.45, 0.20)`.
- The term count at order `r` is `ceil(L d_r / r)` with target site degrees
  `d = (6, 11, 7)`. Supports are unique, contain distinct sites, cover every site,
  and have per-order degree range at most one.
- Each q=20 categorical contrast is a deterministically permuted vector with ten
  `-1` and ten `+1` states. It therefore has exact uniform mean zero, variance one,
  bounded magnitude one, and uniform-nonself substitution energy `40/19`.
- Realized coefficients are rescaled separately by order so that, under a uniform
  site, uniform q=20 source, and uniform nonself replacement, the expected
  single-substitution energy is exactly
  `E[(delta G_r)^2] = v_r E[(delta A)^2]` at working weight one.
- The GNK field is centered at the recorded protein reference. All intervention
  weights are scalar multiples of the same field, so a weight `w` denotes target
  mutation RMS ratio `w` for the normalized arms.

## Arms and paired controls

1. `legacy_public_global_fixed`: the current 384-term public latent GNK with its
   raw global-field normalization and raw multiplier semantics.
2. `legacy_support_local_normalized`: the identical public terms, supports, and
   Gaussian contrasts, but globally rescaled to unit mutation RMS ratio and centered
   at the reference. This isolates local normalization from support and contrast
   changes.
3. `extensive_balanced_bounded`: the primary candidate above.

Within `(reference, seed, regime, distance)`, sequence pairs are generated once and
reused across every arm and weight. Arm and weight are excluded from the sampling
key. Weight zero is the exact causal null. Deterministic fingerprints and
configuration/code-qualified checkpoints audit replay.

## Reference panel and sampling

Six recorded ProteinGym reference sequences span lengths 48, 93, 118, 220, 238,
and 448. Their exact sequences, source dataset/field provenance, and SHA-256 digests
are frozen in `config.yaml`, so execution does not depend on a cluster-local copy of
the full Stage-4 data bundle. Thirty generator seeds (`2601`-`2630`) are independent
replicates.
For each field, 512 paired differences are sampled at Hamming distances 1, 2, 4,
8, `0.1L`, and `0.5L` in two regimes:

- `uniform_background` is the identified mutation operator used to calibrate and
  test the GNK/additive RMS ratio.
- `reference_shell` starts at the recorded natural protein and is the primary test
  of GNK texture over full-protein sequence differences.

The reference-shell GNK difference itself, not its ratio to the additive difference,
is the full-sequence scale-law estimand. Mutations away from the reference can give
the additive component a directional nonzero mean that grows as `L`; requiring a
flat GNK/additive ratio there would conflate that additive drift with GNK scaling.
Reference-shell ratios are therefore reported diagnostically. Uniform-background
ratios retain the primary operator-normalization gate.

## Preregistered gates

All required gates must pass to validate the latent candidate:

- deterministic replay, treatment-independent pair hashes, finite values, an exact
  weight-zero null, and reference centering to `1e-12`;
- complete, unique, nearly regular support at every order and exact bounded q=20
  contrast moments;
- per-length median single-substitution GNK/additive RMS ratio in `[0.85, 1.15]`,
  at least 90% of seed rows in `[0.70, 1.30]`, absolute log-length slope at most
  `0.10`, and its 95% seed-bootstrap interval inside `[-0.20, 0.20]`;
- candidate local slope at least `0.25` closer to zero than the public control, with
  the public defect expected in `[-0.70, -0.30]`;
- analytic `Var(G)` log-length slope in `[0.80, 1.20]` and coefficient of variation
  of `Var(G)/L` at most `0.20`;
- reference-shell fractional-distance GNK RMS slopes in `[0.35, 0.65]`,
  reference-shell fixed-distance absolute slopes at most `0.15`, and uniform-
  background fractional-distance ratio absolute slopes at most `0.15`;
- exact weight linearity, median empirical order-share error at most `0.10`, and
  positive order-stratified co-supported mixed differences with absolute length
  slope at most `0.15` at every order.

The q=20 local-effect kurtosis threshold is advisory and cannot alone reject the
candidate. All log-length regressions use per-length medians over the 30 independent
seeds; intervals resample seeds 10,000 times. Quick mode is explicitly a smoke test
and cannot issue a scientific pass or failure.

## Decision rule

- All required gates pass: `LATENT_GNK_EXTENSIVE_SCALE_LAW_VALIDATED`; a separate
  sigmoid-transmission experiment may then be designed.
- Technical/pairing failure: `EXECUTION_OR_PAIRING_INVALID`.
- Local normalization failure: `OPERATOR_NORMALIZATION_FAILED`.
- Local-normalized control passes but the primary support/contrast gates fail:
  `SUPPORT_OR_CONTRAST_REPAIR_INCOMPLETE`.
- Other preregistered labels localize support, contrast, full-distance, or remaining
  gate failures as recorded in `config.yaml`.

Any failure keeps sigmoid transmission blocked. A pass is evidence only for the
latent construction and weight semantics.

## Execution and reproducibility

- Configuration: `config.yaml`
- Entrypoint: `run.sh`
- Runner: `scripts/extensive_gnk_latent_scaling_exp060.py`
- Candidate: `src/project_adapt_env/extensive_gnk.py`
- Default target: `lab-slurm` (8 CPU, 64 GB RAM, 24 hours)
- Expected outputs: `preflight.json`, `reference_catalog.csv`,
  `pairing_audit.csv`, `support_audit.csv`, `contrast_audit.csv`,
  `latent_scale_metrics.csv`, `order_component_metrics.csv`,
  `scale_law_fits.csv`, `acceptance_checks.csv`, and `summary.json`.

## Analysis plan

Create an ANA record after completion and apply these gates without retuning. Treat
generator seed—not sequence pair—as the independent inferential unit. Preserve the
distinction between uniform-background operator calibration and reference-shell
full-protein behavior.

## Completion criteria

- [x] Configuration and protocol preregistered
- [x] Run script and environment recorded
- [x] Outputs generated remotely by `RUN-149`
- [x] Important outputs DVC-tracked
- [x] Analysis record created (`ANA-058`)
- [x] Result record created (`RES-058`)
- [x] Hypothesis updated (`HYP-005`)
- [x] PROJECT_STATE.md updated after execution

## Post-experiment notes

- `RUN-146` is the initial dry-run record.
- `RUN-147` reached Slurm as job `2039` but failed before scientific computation:
  the shared Stage-4 directory existed while the first declared reference CSV did
  not. No EXP-060 output was generated. The scientific design was unchanged; the
  six already-selected reference sequences and their verified hashes were frozen
  inline to remove that execution-only cluster dependency before retry.
- `RUN-148` is the revised dry-run record after the hash-verified inline-reference
  transport fix.
- `RUN-149` is the corrected live submission, accepted by `lab-slurm` as Slurm job
  `2040` and reconciled `completed` with exit code `0`. Its stdout summary reports
  all `180` tasks complete and all `26/26` required latent checks passing, with no
  stderr.
- `ANA-058` collected and independently audited all `193` output/checkpoint/log
  files; every local SHA-256 digest matches the remote runtime copy. `RES-058`
  records the scoped positive decision: mutation-operator normalization repairs the
  latent scale law, and the extensive/balanced/bounded candidate adds complete
  support and controlled tails before any sigmoid.
- The fixed-support locally normalized control also repairs local, field, and
  full-distance scaling, so the experiment does not establish that extensive
  support or bounded contrasts are individually necessary for scaling. The runner
  omitted the prose-specified `h = 1` fixed-distance gate; independent analysis of
  the existing rows gives slope `-0.001481`, which passes the same threshold and
  leaves the registered decision unchanged.
