# RES-058: RUN-149 validates mutation-normalized GNK latent scaling before sigmoid transmission

## Summary

RUN-149 passes all `26/26` preregistered required latent gates and both advisory
checks. Under the declared q=20 mutation operator, mutation-operator normalization
removes the protein-length-dependent GNK amplitude defect, while the experiment-local
extensive, balanced, bounded candidate preserves local effects, variance proportional
to length, and approximately square-root GNK differences across protein-scale
sequence distances. This supports HYP-005 only at the latent, pre-sigmoid boundary.

## Generated from

- Analysis: ANA-058
- Experiment: EXP-060
- Run: RUN-149

## Relevant hypotheses

- Supports: HYP-005
- Weakens: None
- Refutes: None

## Evidence

- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/metrics.json`
- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/run-149_collection_audit.csv`
- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/run-149_acceptance_audit.csv`
- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/run-149_local_ratio_by_length.csv`
- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/run-149_scale_law_summary.csv`
- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/run-149_candidate_support_by_length_order.csv`
- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/run-149_order_allocation_and_curvature.csv`
- `analyses/2026-07-23_ANA-058_run-149-extensive-operator-normalized-gnk-latent-scale-law-review/tables/run-149_structural_audit.csv`
- `data/processed/extensive-operator-normalized-gnk-latent-scale-law-validation/RUN-149.dvc`
- `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-149/outputs`

## Interpretation

`RUN-149` completed the full EXP-060 design: six recorded protein references at
lengths `48-448`, 30 independent GNK seeds, 180 reference/seed tasks, 2,160 paired
sequence designs, and 32,400 metric rows. All 193 collected output, checkpoint, and
log files match their remote SHA-256 digests; stderr is empty and the scheduler exit
code is `0`.

The primary local-scale result is positive. The candidate single-substitution GNK
RMS has log-length slope `-0.006316` (95% paired seed-bootstrap interval
`[-0.012404, 0.006453]`), and every per-length median GNK/additive RMS ratio lies
between `0.9912` and `1.0120`. All `180/180` seed rows fall in the registered broad
replication band. The frozen Stage-4 384-term public-form control reproduces the
predicted `L^-1/2` defect: slope `-0.510364`, interval
`[-0.535205, -0.482253]`.

The candidate also passes the global/full-distance laws. Analytic field variance has
slope `0.995076` with `CV[Var(G)/L] = 0.006416`. Reference-shell GNK difference RMS
has slopes `0.481963` at distance `0.1L` and `0.430338` at `0.5L`, within the
registered square-root band. At fixed distances `h = 2, 4, 8`, slopes are near zero
(`0.001964`, `0.004063`, `0.045973`). The runner omitted the preregistered `h = 1`
fixed-distance gate; independent analysis gives `-0.001481`, which passes the same
limit and leaves the decision unchanged. Uniform-background fractional-distance
GNK/additive ratio slopes are also near zero (`-0.005255`, `-0.006916`). Post hoc
reference-shell GNK SD slopes are `0.491618/0.495854`, confirming that texture
fluctuations retain square-root scaling despite additive directional drift in raw
reference-shell RMS ratios. The preregistered gate is an across-seed ensemble
criterion: a post hoc check finds `30/30` individual seed-wise slopes inside the
nominal `[0.35, 0.65]` band at `0.1L`, but `24/30` at `0.5L`, so the latter result
is not uniform over all GNK realizations.

The structural implementation behaves as registered: supports have no uncovered
sites or duplicates, maximum degree range one and maximum degree CV `0.01840`;
bounded q=20 contrasts have exact registered moments to numerical precision; maximum
order mutation-energy-share error is `0.00906`; and co-supported order-specific mixed differences
remain positive without a material length trend. Scalar weight response is exact to
numerical precision.

The order-share output columns are misnamed as variance fractions: they quantify
single-mutation-energy shares. Targets `0.35/0.45/0.20` imply actual uniform-
background field-variance shares `0.4667/0.4000/0.1333` at orders 2/3/4. Curvature
is likewise conditional on a co-supported pair rather than averaged over arbitrary
pairs; representative union coverage decreases from `66.9%` at `L = 48` to `10.4%`
at `L = 448`.

The fixed-support, locally normalized legacy control is the key causal boundary. It
removes the local trend (slope `-0.013145`) and also recovers field/full-distance
scaling: empirical field-variance slope `0.99397`, uniform fractional-distance GNK
slopes `0.49809/0.49518`, and reference-shell slopes `0.47074/0.47841`. Therefore
operator normalization is sufficient for the observed scale-law repair. The
extensive/balanced/bounded candidate's added value is uniform coverage and controlled
tails: local kurtosis median/95th/max are `3.086/3.600/4.142`, versus
`9.402/29.586/80.113` for the normalized frozen control; at `L = 448` the latter
leaves median uncovered-site counts `244.5/138.5/223.5` at orders 2/3/4. The bundle
passes, but individual necessity of its support and contrast changes is not shown.

## Effect on hypothesis

This supports HYP-005 narrowly: the experiment-local GNK weight has stable,
interpretable mutation-RMS semantics over the tested q=20 operator, six references,
and length range. The frozen raw-amplitude semantics are shown to be defective, while
the normalized latent candidate is directional and reproducible across generator
seeds.

RES-055 and RES-057 remain valid at the observable level. They tested GNK after the
current public nonlinear readout; EXP-060 deliberately stops before that operation.
This result therefore does not validate the current public generator or supply direct
evidence for HYP-004 or HYP-009.

## Limitations

- One reference protein is used per length, so protein identity and length are
  confounded. Seed bootstraps quantify GNK-field replication conditional on the six
  references, not generalization across a population of proteins.
- Normalization is identified for the declared uniform-site, uniform-background q=20
  nonself mutation operator; other proposal distributions are not tested.
- Extensive term counts, regularized support, and bounded contrasts were not
  individually ablated after local normalization.
- The candidate remains project-local. The public AdaptEnv implementation was not
  changed. The implementation hash covers the runner and candidate field but omits
  the declared dirty external AdaptEnv source used for the additive core/control;
  current-source metrics reproduce within numerical tolerance, so this is a
  provenance limitation rather than evidence of invalid data.
- Deterministic fingerprints replay bitwise within the run environment but are not
  cross-platform bitwise stable; reconstructed scientific metrics agree to about
  `1e-14`.
- The fixed-distance `h = 1` check was omitted from the runner's registered check
  table, although the independent slope passes.
- No sigmoid transmission, observable-fitness topology, accessibility, adaptive walk,
  empirical realism, biological routing, or optimizer performance was tested.

## Downstream use

- Cite RUN-149 as validation of the latent mutation-normalized scale law and weight
  semantics, not as validation of observable fitness or empirical realism.
- Treat the latent prerequisite for a distinct sigmoid-transmission test as satisfied.
- Keep the Stage 4C/4D/4E hard stop in place until a separately registered observable-
  transmission experiment passes its own gates.
