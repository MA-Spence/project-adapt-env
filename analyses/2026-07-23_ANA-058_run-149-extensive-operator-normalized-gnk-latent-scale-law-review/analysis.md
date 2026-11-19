# ANA-058: RUN-149 extensive operator-normalized GNK latent scale-law review

## Purpose

Collect and independently audit the completed RUN-149 latent GNK outputs against
every preregistered EXP-060 scale, support, contrast, pairing, and curvature gate,
without retuning thresholds or adding post hoc success criteria.

## Linked experiments/runs

- Experiment: EXP-060
- Run: RUN-149
- Direct hypothesis: HYP-005

## Collection and reproducibility audit

`RUN-149` completed on `lab-slurm` as Slurm job `2040` with scheduler exit code
`0`. The configured labproj client mount (`/Volumes/DutchZFS/labproj`) was not
available, so the same run-scoped promotion was performed directly over SSH from
the recorded runtime paths. The ten declared outputs, three logs, and 180
checkpoint files were copied to:

- `data/processed/extensive-operator-normalized-gnk-latent-scale-law-validation/RUN-149`

All `193/193` local files match their remote SHA-256 digests. The promoted
directory is DVC-tracked at:

- `data/processed/extensive-operator-normalized-gnk-latent-scale-law-validation/RUN-149.dvc`

`stderr.txt` is empty and `labproj-exit-code.txt` contains `0`. The run is the
full preregistered configuration (`quick = false`): six references at lengths
`48, 93, 118, 220, 238, 448`, 30 GNK generator seeds, 180 reference/seed tasks,
2,160 paired sequence designs, and 32,400 latent-metric rows. Pair hashes are
treatment-independent, all three arms use each paired design, all underlying
additive/GNK differences are finite, and the weight-zero rows are exactly zero.
Nonapplicable CSV fields contain expected NaNs. The 180 checkpointed field
fingerprints reproduce within the run environment; this is same-platform replay,
not a cross-platform bitwise guarantee.

## Preregistered decision

All `26/26` required checks pass. Both advisory checks also pass. The registered
decision is therefore:

`LATENT_GNK_EXTENSIVE_SCALE_LAW_VALIDATED`

This is a latent-stage decision. The run confirms its own scope controls: it did
not transmit GNK through a sigmoid, calculate observable fitness, accessibility,
topology, or optimizer metrics, or change the public AdaptEnv implementation.

There is one preregistration-to-runner coverage omission. The prose required
reference-shell fixed-distance checks at `h = 1, 2, 4, 8`, but the runner gated only
`h = 2, 4, 8`. Independent analysis of the existing `h = 1` rows gives slope
`-0.001481`, well inside the same `+/-0.15` limit, so the omission does not change
the decision. It is nevertheless a protocol deviation and is not counted among
the runner's `26` checks.

## Local mutation scale

Under the declared uniform-background, uniform-site, uniform-nonself q=20
mutation operator at weight one, the candidate has no material length trend:

- candidate log-length slope: `-0.006316`;
- 95% paired seed-bootstrap interval: `[-0.012404, 0.006453]`;
- per-length median `RMS(delta GNK) / RMS(delta additive)`: `0.9912-1.0120`;
- seed rows in the broad `[0.70, 1.30]` replication band: `180/180` (`100%`).

The frozen Stage-4 384-term public-form control reproduces the preregistered scale
defect:

- public-control slope: `-0.510364`;
- 95% paired seed-bootstrap interval: `[-0.535205, -0.482253]`.

The slope improvement in the candidate is `0.5040`, exceeding the registered
`0.25` requirement. Weight response is scalar-linear to numerical precision
(runner maximum relative error `5.40e-16`, minimum `R2 = 1`; independent
weight-one reconstruction maximum relative error `3.38e-15`).

The fixed-support locally normalized legacy control is an important causal
qualification. It removes the local trend (slope `-0.013145`) and also shows the
same diagnostic scale laws: empirical field-variance slope `0.99397`, uniform-
background fractional-distance GNK slopes `0.49809/0.49518`, and reference-shell
GNK slopes `0.47074/0.47841`. Thus mutation-operator normalization is sufficient
for the observed scale-law repair in this experiment. The current AdaptEnv default
uses 256 terms rather than the frozen control's 384, but fixed term count changes
the intercept rather than the expected `L^-1/2` defect. The data do not establish
that extensive support, balanced q=20 contrasts, or bounded contrasts are
individually necessary for scaling.

## Field and full-sequence scaling

The extensive candidate satisfies the registered global and full-distance laws:

- analytic `Var(G)` log-length slope: `0.995076`, 95% seed-bootstrap interval
  `[0.988914, 1.007004]`;
- `CV[Var(G)/L]`: `0.006416`;
- empirical source-field variance slope (supportive, not the analytic gate):
  `0.990399`;
- reference-shell GNK RMS slope at `0.1L`: `0.481963`, interval
  `[0.452315, 0.498355]`;
- reference-shell GNK RMS slope at `0.5L`: `0.430338`, interval
  `[0.383979, 0.494390]`;
- reference-shell GNK SD slopes at `0.1L` and `0.5L` (post hoc texture
  diagnostic): `0.491618` and `0.495854`;
- reference-shell fixed-distance slopes at `h = 2, 4, 8`: `0.001964`,
  `0.004063`, and `0.045973`;
- uniform-background fractional-distance GNK/additive ratio slopes at `0.1L`
  and `0.5L`: `-0.005255` and `-0.006916`.

The registered fractional-distance gates concern the across-seed median slope and
its paired bootstrap interval. In a post hoc seed-wise variability check, `30/30`
individual slopes at `0.1L` but `24/30` at `0.5L` lie in the nominal
`[0.35, 0.65]` square-root band. Thus the `0.5L` ensemble result passes, but is not
uniform over every GNK realization.

The fractional-distance reference-shell GNK/additive ratios decline with length
(slopes `-0.5009` and `-0.5626`). This was expected and was not the registered
full-sequence estimand: mutations away from a natural reference introduce an
additive directional component whose mean grows with distance. The registered
reference-shell estimand was GNK RMS itself, which shows the predicted
approximately square-root scaling. The post hoc SD calculation sharpens the
interpretation: GNK fluctuations retain square-root scaling even though the raw
GNK/additive RMS ratio falls as additive directional drift accumulates.

## Support, contrast, and order audit

The candidate construction passes all registered structural checks:

- no uncovered sites and no duplicate supports;
- maximum per-order degree range `1` and maximum degree CV `0.01840`;
- q=20 contrasts have exact zero mean and unit second moment to the registered
  tolerance, maximum absolute value `1`, and substitution-energy error
  `8.88e-16`;
- maximum median order mutation-energy-share error `0.00906` against the `0.10`
  limit;
- co-supported mixed-difference RMS is positive at every order (minimum
  `0.34012`), with log-length slopes `-0.00070`, `-0.04734`, and `-0.03970` for
  orders 2, 3, and 4;
- advisory median q=20 local-effect Pearson kurtosis is `3.0857`, below `6.0`.

The output columns named `target_variance_fraction` and `sample_variance_fraction`
are misnamed: the implementation computes single-mutation-energy shares. The
registered targets `0.35/0.45/0.20` therefore correspond to actual uniform-
background field-variance shares `0.4667/0.4000/0.1333` for orders 2/3/4, not
`0.35/0.45/0.20`. This does not invalidate the mutation-energy allocation gate,
but it matters for downstream interpretation of the order spectrum.

The candidate's additional value over normalization is coverage and tail control.
Its local-effect kurtosis has median `3.086`, 95th percentile `3.600`, and maximum
`4.142`, versus `9.402`, `29.586`, and `80.113` for the locally normalized frozen
control. At `L = 448`, that control leaves median uncovered-site counts
`244.5/138.5/223.5` at orders 2/3/4, while the candidate leaves none.

The curvature check is conditional on drawing a co-supported pair, not an average
over arbitrary site pairs. Representative union coverage falls from `66.9%` at
`L = 48` to `10.4%` at `L = 448`; the positive, length-stable curvature result must
be cited with that conditioning.

These observations validate the extensive, balanced, bounded candidate as a
working bundle. Because the experiment does not separately ablate extensive
term counts, regularized support, and bounded contrasts while retaining all
other candidate features, it cannot assign necessity to any one of them.

## Interpretation for HYP-005

RUN-149 supports HYP-005 narrowly. For the declared q=20 mutation operator and
the six fixed reference sequences, the experiment-local candidate gives the GNK
weight a stable mutation-RMS meaning from single substitutions through
protein-scale sequence differences. It repairs the `L^-1/2` amplitude loss of
the raw fixed-term public control and preserves order-specific nonlinear texture.

This does not reverse RES-055 or RES-057 at the observable level. Those results
tested the current public generator after nonlinear readout; EXP-060 deliberately
stops before that boundary. The result clears the latent prerequisite for a
separate sigmoid-transmission test, but it does not validate the public method or
lift the Stage 4C/4D/4E hard stop.

## Limitations

- There is one reference protein per length, so protein identity and length are
  confounded. The paired seed bootstrap quantifies GNK-field replication
  conditional on these six references, not population-wide generalization across
  proteins.
- The mutation-RMS semantics are identified for the declared uniform-background
  q=20 nonself operator. Different proposal distributions require their own
  normalization or validation.
- The candidate is experiment-local. No public AdaptEnv code or configuration was
  changed. The implementation digest covers the runner and candidate field but not
  the declared, dirty external AdaptEnv source used to construct the additive core
  and frozen control. Current source reconstructs scientific metrics within numerical
  tolerance, making this a provenance limitation rather than evidence of data error.
- The deterministic replay hash is bitwise within the execution environment, but
  locally reconstructed fingerprints differ across platform-level floating-point
  ordering even though scientific RMS values agree to about `1e-14`.
- No sigmoid transmission, observable-fitness topology, empirical assay overlap,
  biological interaction routing, accessibility, adaptive walks, or optimizer
  performance was tested. HYP-004 and HYP-009 therefore receive no direct support.

## Produced artifacts

- `source_tables/run-149_summary.json`
- `source_tables/run-149_preflight.json`
- `tables/metrics.json`
- `tables/run-149_collection_audit.csv`
- `tables/run-149_acceptance_audit.csv`
- `tables/run-149_local_ratio_by_length.csv`
- `tables/run-149_scale_law_summary.csv`
- `tables/run-149_candidate_support_by_length_order.csv`
- `tables/run-149_order_allocation_and_curvature.csv`
- `tables/run-149_structural_audit.csv`

The paired notebook and `notebooks/analysis.py` invoke
`notebooks/derive_tables.py`; rerunning that script from the project environment
regenerates all derived tables and asserts the audited decision from the DVC-tracked
RUN-149 inputs.
