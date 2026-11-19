# ANA-059: RUN-151 operator-normalized GNK sigmoid-transmission paired-residual review

## Purpose

Collect and independently audit completed RUN-151 against every preregistered
EXP-061 technical, controlled-edge, texture-transmission, and deployment gate,
preserving the frozen-core versus actual-core mechanistic boundary without post
hoc retuning.

## Linked experiments/runs

- Experiment: EXP-061
- Run: RUN-151
- Direct hypothesis: HYP-005

## Collection and reproducibility audit

`RUN-151` completed on `lab-slurm` as Slurm job `2041` with scheduler exit code
`0`. It is the full preregistered configuration (`quick = false`): six fixed
references at lengths `48, 93, 118, 220, 238, 448`, 30 generator seeds, and
`180/180` completed reference/seed tasks. The aggregate files contain the declared
`3,240` design rows, `46,080` paired-residual rows, `171,520` mutation-transmission
rows, `146,880` epistasis-transmission rows, and `112,320` topology-audit rows.

The configured client mount at `/Volumes/DutchZFS/labproj` was unavailable, so
the immutable run-scoped outputs were collected directly over SSH to:

- `data/processed/operator-normalized-gnk-sigmoid-transmission-paired-residual-validation/RUN-151`

The bundle contains ten declared outputs, three logs, and 180 task checkpoints.
All `193/193` local files match both the remote runtime SHA-256 manifest and their
DVC object MD5 digests. `stderr.txt` is empty and `labproj-exit-code.txt` contains
`0`. The DVC pointer is:

- `data/processed/operator-normalized-gnk-sigmoid-transmission-paired-residual-validation/RUN-151.dvc`
- directory hash: `5aace1d2b1f33df736417bf314cfe5bb.dir`
- data size: `1,513,486,733` bytes across `193` files

Because the mounted-path `dvc push` could not run, the exact 193 data objects and
directory manifest were copied by content address to the configured remote store
over SSH. A checksum-aware dry run then found zero remaining objects. This is a
manual content-addressed backup fallback, not a successful invocation of
`dvc push`.

The run-recorded hashes for the EXP-061 configuration, EXP-060 reference source,
both runners, candidate implementation, and the three relevant public AdaptEnv
files match the present files. The repository and AdaptEnv checkout were dirty at
submission, but the scientific source components used here are individually
hash-frozen in `preflight.json`.

## Preregistered decision

All 12 technical checks pass. Three scientific checks fail, giving `13/16`
required checks overall:

- reference-anchored controlled-edge transmission: fail;
- frozen-additive-core texture transmission: pass;
- actual-additive-core texture transmission: fail;
- direct margin-8 deployment: fail.

The independently reproduced flags are therefore `technical = true`,
`controlled_edges = false`, `frozen_texture = true`, `actual_texture = false`,
and `public = false`. The registered decision is:

`SIGMOID_TRANSFER_FAILED`

This label follows the preregistered branch ordering: controlled-edge failure is
tested before the frozen-core-only branch. It should not be read as evidence that
a sigmoid is mathematically unable to transmit GNK texture. The pointwise and
frozen-core controls show that it can. The failure is specifically broad
common-threshold transmission with the real additive core and the declared
deployment operating points.

## Pointwise scalar transfer

The per-genotype occupancy calibration passes, as expected for the known scalar
transfer function. At occupancy `0.5` and amplitude `a = 0.5`:

- the minimum per-length median normalized secant gain is `0.990650`;
- the minimum per-length practical-resolution fraction is `0.934653`;
- median residual RMS is `1.1706-1.1926` across lengths;
- median linear correlation is `0.998946` and NRMSE is `0.076443`;
- the log-length slope is `-0.002053`, with paired-seed 95% interval
  `[-0.012800, 0.007030]`.

At the `0.1/0.9` shoulders, the registered minimum median gain is `0.358714`,
also above its `0.20` bound. At occupancy `0.99`, gain is only about
`0.0395-0.0400` and the practical fraction is `0.078-0.085`, correctly exposing
plateau attenuation.

This assay changes `theta` for every genotype. It establishes transfer-function
implementation and absence of a renewed pointwise length trend, but it is not one
fitness landscape and cannot support topology or public-method claims.

## Common-threshold reference-anchored edges

Anchoring the reference core at occupancy `0.5` does not keep its single-mutant
neighbors in the sigmoid transition region. Across lengths, the median mutant
core occupancy is only `1.58e-21-2.36e-21`; `98.47-98.59%` of complete one-step
panels are one-crossing and only `1.41-1.53%` are both-dynamic.

The retained both-dynamic subset behaves well where it is large enough:

- eligible-seed median correlation is `0.912-0.929`;
- eligible-seed median practical fraction is `0.852-0.889`;
- the residual-RMS log-length slope is `0.031388`, with interval
  `[-0.075226, 0.073136]`.

However, the required `n >= 32` strata occur in only `0/30`, `4/30`, `17/30`,
`30/30`, `30/30`, and `30/30` seeds at lengths `48, 93, 118, 220, 238, 448`.
The minimum required per-length coverage is `80%`, so the gate fails. This is
partly a mechanical interaction between the absolute `n >= 32` cutoff and a
complete panel of size `19L`: the both-dynamic fraction itself is approximately
constant, while the absolute count grows with length. It is not evidence for a
biological attenuation specific to the shorter proteins.

The result is nevertheless not only a cutoff artifact. Across all one-step edges,
median correlation is just `0.113-0.125` and the median practical-resolution
fraction is `0.0148-0.0167`. Strong transmission is confined to a rare dynamic
subset rather than broadly present in the local landscape.

## GNK texture through the sigmoid

The candidate-additivized contrast
`tau_texture = epsilon_R(actual G) - epsilon_R(additivized G)` is decisive.

With the additive core frozen at the reference operating point, every registered
texture condition passes:

- all lengths have complete seed coverage and positive seed projections;
- the minimum per-length projection bootstrap lower bound is `43.9273`;
- the `a = 1` practical fraction is `0.2617-0.3926`, above the `0.25` gate;
- all `90/90` required permutation tests pass at every length;
- median linear correlation is at least `0.99746`;
- transmitted texture RMS is `0.344-0.438`;
- its log-length slope is `-0.101001`, with interval
  `[-0.118187, -0.088581]`.

Thus the sigmoid can carry the EXP-060 normalized GNK texture, including every
tested `L >= 100`, when additive-core motion is removed and the square remains in
the transition window.

With the actual additive core, the same test fails at every length. Projection is
formally positive and correlation remains almost one, but the transmitted signal
has collapsed to texture RMS `2.32e-12-1.21e-10`. The practical fraction is exactly
zero at all six lengths, only `0.289-0.467` of the required permutation tests pass,
and the slope is `-0.331937` with interval `[-1.311418, 0.577656]`.

Raw observable epistasis is not evidence of GNK transmission. Candidate
`epsilon_R` RMS remains `0.103-0.151` despite the tiny actual-core
`tau_texture`. The order-1 additive surrogate likewise has latent epsilon at
numerical zero (maximum `2.84e-14`) but observable `epsilon_R` RMS
`0.115-0.134`, with nonzero observable epistasis in `12,960` rows. These controls
show that raw `epsilon_R` is dominated by sigmoid/global-epistasis curvature and
support the preregistered additivized contrast as the relevant GNK-texture
estimand.

## Deployment operating points

The direct margin-8 base-template formula fails practical deployment. One-step
candidate residuals are exact float64 ties in about `98.5-98.7%` of edges, and the
per-length practical fractions are `0`, `0.000566`, `0.000892`, `0.000957`,
`0.000995`, and `0.000822`. Stable analytic residuals have the same practical
fractions, so avoiding direct-subtraction ties does not rescue a practically
resolved effect.

Supported plaquettes retain rare directional signal at margin 8: every length has
positive seed projection, projection bootstrap lower bounds are `0.466-0.533`,
and permutation-pass fractions are `0.844-0.944`. But `80.9-83.6%` of direct
texture rows are exact ties, and only `1.56-2.54%` exceed the `0.10` heuristic
floor, far below the required `25%`. This is evidence for rare threshold-crossing
texture, not broadly usable transmission.

The mandatory margin-10 constructor comparator is more saturated. Its direct
one-step practical fraction is zero at every length, direct ties are
`99.95-100%`, and the texture practical fraction is `0.003906` throughout. It is
diagnostic only and cannot substitute for the failed margin-8 primary.

## Protocol and output audit

Independent reaggregation reproduces all 16 pass states and the registered
decision. Five qualifications do not change it:

- The absolute controlled-edge `n >= 32` cutoff is length-confounded through the
  `19L` panel size, as described above.
- The epistasis table does not retain all square endpoint core occupancies or the
  additivized double endpoint promised by the prose. Saturation of the square's
  double endpoint is therefore inferred from the edge occupancy and
  actual-versus-frozen contrast rather than directly tabulated.
- The registered texture gate takes each seed's median of the three support-order
  projections rather than one pooled `sum_xy/sum_x2` projection. An independent
  pooled calculation leaves every positive-projection pass state unchanged.
- The runner applies the NRMSE `<= 0.25` linearization gate at occupancy `0.5`.
  Shoulder NRMSE is about `0.312-0.323`, so the prose is ambiguous if read as
  applying the same NRMSE bound at `0.1/0.9`; the separately stated shoulder gate
  is gain `>= 0.20`, which passes.
- The `finite_outputs` label is broader than its implementation: the task check
  covers raw arrays and selected RMS values, while nonapplicable CSV fields retain
  expected NaNs. No consequential nonfinite scientific value was found.

## Interpretation for HYP-005

RUN-151 preserves RES-058: there is no renewed latent or pointwise length-scaling
defect in the normalized GNK field. It also establishes a positive mechanistic
boundary: the frozen-core experiment transmits length-stable GNK texture through
the sigmoid.

It weakens HYP-005 at the intended observable protein-scale boundary. With a
common threshold and the real additive core, nearly all tested mutations leave
the dynamic window, candidate-additivized texture is not practically resolved at
any length, and the margin-8/10 deployment surfaces are saturated. The tested
observable therefore does not provide a broadly directional, practically
resolvable GNK texture knob.

This does not refute every additive-core plus GNK model, invalidate the EXP-060
latent normalization, or show that a sigmoid cannot transmit GNK texture. It
localizes the bottleneck to the interaction between additive-core motion and the
steep saturating readout under these thresholds.

## Limitations and scope

- There is one protein per length, so protein identity and length are confounded.
  Seed inference is conditional on these six fixed references.
- The additive core and GNK realization vary under the same seed-level task; the
  experiment isolates frozen versus actual core mode but does not independently
  factorialize every core realization against every GNK realization.
- The `0.10` practical floor is a preregistered heuristic for this noiseless assay,
  not an empirical service detection limit. Observation noise was disabled.
- The candidate is experiment-local and replaces the public fixed-support
  generator in the primary arm. Margin 8 is the documented base-template formula,
  not a universal API/service margin. The public AdaptEnv model was not changed.
- No empirical fitness realism, biological interaction routing, accessibility,
  adaptive walk, active-learning, or optimizer metric was tested. HYP-004 and
  HYP-009 receive no direct evidence, and the Stage 4C/4D/4E hard stop remains.

## Produced artifacts

- `source_tables/run-151_summary.json`
- `source_tables/run-151_preflight.json`
- `source_tables/run-151_remote_sha256.json`
- `tables/metrics.json`
- `tables/run-151_collection_audit.csv`
- `tables/run-151_source_hash_audit.csv`
- `tables/run-151_acceptance_audit.csv`
- `tables/run-151_pointwise_transfer_by_length.csv`
- `tables/run-151_reference_anchored_edge_by_length.csv`
- `tables/run-151_texture_transmission_by_length.csv`
- `tables/run-151_deployment_margin_comparison_by_length.csv`
- `tables/run-151_scale_law_summary.csv`
- `tables/run-151_protocol_output_audit.csv`

The paired notebook and `notebooks/analysis.py` invoke
`notebooks/derive_tables.py`. Rerunning that script regenerates every derived table
from the DVC-tracked RUN-151 bundle, verifies the persisted remote SHA-256 manifest
and DVC object digests, and asserts the registered decision.
