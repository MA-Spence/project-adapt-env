# ANA-060: RUN-153 dimensionless additive-core sigmoid calibration and held-out GNK dose-response review

## Purpose

Collect and independently audit completed RUN-153 against the preregistered
EXP-062 core-only calibration, held-out core, paired-residual transmission,
dose-directionality, length-stability, control, and provenance gates without
post hoc retuning.

## Linked experiments/runs

- Experiment: EXP-062
- Run: RUN-153
- Direct hypothesis: HYP-005

## Collection and reproducibility audit

`RUN-153` completed on `lab-slurm` as Slurm job `2042` with scheduler exit code
`0`. It is the full registered configuration (`quick = false`): 12 calibration
and 12 confirmation proteins, six length bins, five calibration core seeds,
five disjoint confirmation core seeds, and six crossed confirmation GNK seeds.
All `360/360` confirmation protein/core/GNK tasks completed. `stderr.txt` is
empty, and the aggregate tables contain the declared `17,280` edge rows and
`78,120` texture rows.

The collected run bundle is:

- `data/processed/dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-validation/RUN-153`
- `373` files and `332,322,278` bytes: ten registered outputs, 360 task
  checkpoints, and three logs
- DVC pointer
  `data/processed/dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-validation/RUN-153.dvc`
- DVC directory hash `02986be42c5b575ff8d7b75b3abad31c.dir`
- backed up successfully to the configured `labcluster` DVC remote

Every local file matches both its DVC object MD5 and the immutable shared-storage
runtime file by SHA-256 and size. A canonical sorted SHA-256 manifest over all
373 paths has digest
`c06f6fbd8c90f84462456427f537879a27db25b6573e0043f5f81327fc86ddbd`.
The run-recorded experiment config, reference catalog, seven implementation
components, effective configuration, and implementation aggregate hashes match
the present registered files. The repository and AdaptEnv checkout were dirty at
submission, but the scientific source components used here are individually
hash-frozen in `preflight.json`.

The 360 checkpoints form the exact registered `12 x 5 x 6` crossed design, with
no missing or extra task. Their edge and texture row totals reproduce the
aggregate files. One portability defect is retained in the immutable evidence:
the checkpoint files contain `401,773` Python-style literal `NaN` tokens. Python's
JSON loader accepts them, but they are not strict RFC 8259 JSON and cannot be
assumed portable to strict JSON parsers. This does not alter any collection hash
or scientific pass state.

## Registered and independently audited decision

The runner reports `6/9` required checks passing and issues:

`CALIBRATED_SIGMOID_TRANSMISSION_FAILED`

The three formal failures are:

- global actual-core GNK texture transmission;
- global actual-core dose directionality;
- the frozen-core texture positive control.

Independent reaggregation reproduces every literal runner pass state. It also
finds that the directionality failure is solely a floating-point boundary
artifact. Fifty-two practical-fraction Spearman values mathematically equal to
the registered `0.80` threshold are stored as `0.7999999999999999` and fail a
literal `>= 0.80` comparison. Applying a numerical tolerance of `1e-12` changes
the independent count to `7/9` and makes directionality pass. It does not change
the overall decision because the actual and frozen texture-prevalence failures
remain large.

## Core-only calibration and held-out confirmation

The complete-h1 calibration estimate is
`sigma_A = 2.7448641783`. Across the 60 calibration protein/core cells its
coefficient of variation is only `0.02274`, with range `2.59743-2.89401`. The
core-only selection rule chooses:

- `beta = 4`;
- `mu = logit(0.99) = 4.5951198501`;
- raw global steepness `beta/sigma_A = 1.457267`.

No GNK outcome is present in the selection table. All ten selected calibration
cells in every length bin pass the registered high-reference, dynamic-h1,
core-effect, and distant low-plateau gates.

The globally frozen calibration then passes all `60/60` held-out protein/core
cells, with a pass fraction of one in every length bin. Median held-out
reference fitness is `9.9`, median h1 dynamic fraction is `0.821-0.828`, median
h1 core-effect RMS is `4.15-4.35`, and the `0.5L` low-plateau median is one in
every bin. Realized beta under the global scale ranges from `3.8487` to `4.2042`,
already close to the exact-per-landscape value of four. This is strong
fixed-panel evidence that the additive-core operating-scale calibration
generalizes to the held-out proteins.

## Complete-h1 edge transmission

At `a=1`, every one of the `360/360` actual-core complete-h1 cells passes the
edge gate. Across cells:

- median practical-resolution fraction is `0.8828`;
- median latent/residual correlation is `0.8717`;
- median paired-residual RMS is `1.5342`.

The minimum per-bin pass fraction is one. This is a substantive repair of the
approximately `1.5%` practical one-step transmission in RES-059: mutation-unit
calibration makes the GNK perturbation broadly visible on complete reference
single-mutant panels.

## Actual-core GNK texture

The actual-core candidate-additivized texture at `a=1` is present, directional,
and non-negligible:

- median texture RMS is `0.4175`;
- median practical-resolution fraction is `0.1927`;
- median latent/residual correlation is `0.7322`;
- median projection is `0.9182`;
- every cell has positive projection, correlation at least `0.50`, and
  one-sided permutation `p=0.002`;
- the crossed-bootstrap projection lower bound is `0.8731`;
- the actual/frozen RMS ratio has median `3.286` and minimum `1.263`.

The registered transmission failure is entirely the practical-prevalence
criterion. Only `24/360` cells reach the required practical fraction of `0.25`:
23 are in the short bin and one is in `L118`; none pass in `L093`, `L220`,
`L238`, or `L448`. The per-bin texture pass fractions are therefore
`0-0.383`, far below the required `0.80`.

The correct interpretation is not that the calibrated sigmoid erases GNK
texture. It transmits an aligned texture with appreciable RMS, but that texture
does not exceed the preregistered `0.10` per-contrast floor in a sufficiently
large fraction of supported plaquettes.

## Frozen-core control and exact normalization

The frozen-core control formally fails:

- median texture RMS is `0.1317`;
- median practical-resolution fraction is `0.0807`;
- median correlation is `0.5954`;
- median projection is `0.2224`;
- projection and permutation pass in all `360/360` cells;
- correlation passes in `335/360`;
- practical prevalence passes in `0/360`.

This control is frozen at the newly selected reference occupancy of `0.99`,
deep in the sigmoid's upper shoulder. It therefore behaves as a saturation
control rather than the intended positive transmission control. RES-059's
frozen-core assay at occupancy `0.5` remains the positive evidence that the
sigmoid can transmit the normalized texture in a controlled transition window.
The preregistered EXP-062 control failure is retained, but it cannot be read as
evidence that aligned GNK texture is absent.

The exact-per-landscape beta diagnostic does not rescue texture. It gives median
RMS `0.4153`, practical fraction `0.1901`, correlation `0.7297`, and projection
`0.9095`; only `23/360` cells reach the practical-prevalence floor. Thus
per-landscape core normalization is not the missing ingredient.

The order-1 additive surrogate behaves correctly. Its latent texture is at
numerical zero (maximum absolute latent epsilon `1.42e-14`), transmitted texture
RMS is about `9.20e-16`, and practical fraction is zero, even though raw
observable epistasis RMS is about `1.88`. This again shows why raw sigmoid
epistasis is not evidence of GNK-specific texture.

## Dose response and numerical sensitivity

For the global actual-core arm, median texture RMS across the registered doses
`a = 0.25, 0.5, 1, 2` is:

`0.1085, 0.2163, 0.4175, 0.7578`.

All 360 cell-level RMS Spearman correlations are one. All `12/12` protein-median
curves pass both RMS and practical-fraction directionality, and the registered
protein-bootstrap lower growth bounds are positive (`0.6099` RMS and `0.08138`
practical fraction).

Under the runner's literal comparisons, `305/360` cells pass and the minimum
bin fraction is `0.7333`. With `1e-12` tolerance at the exact `rho=0.80`
boundary, `354/360` cells pass; per-bin fractions become `1.0, 1.0, 1.0,
0.9833, 0.9167, 1.0`, all above `0.80`. The remaining six cells genuinely miss
the practical-fraction growth floor. The scientific record should therefore
retain the machine-reported formal failure but not describe the response as
nonmonotonic.

At the excluded `a=4` saturation stress, median texture RMS continues upward to
`1.1846` while median practical fraction falls from `0.2005` at `a=2` to
`0.1771`. This is consistent with the registered expectation that a bounded
sigmoid can increase total contrast magnitude while flattening the fraction of
contrasts above a fixed practical floor.

## Length stability and raw operating-point controls

The global protein-level `a=1` texture-RMS slope is `-0.10221`, with stratified
bootstrap interval `[-0.11388, -0.09012]`. The exact arm slope is `-0.11358`,
with interval `[-0.12531, -0.10255]`. Both pass the registered slope and
interval bounds. This is tolerance-defined length stability with a modest
negative trend, not strict invariance.

The edge-only operating-point ablations remain poor:

| Arm | Median practical fraction | Median correlation |
|---|---:|---:|
| Steepness-only, margin 8 | 0.03115 | 0.2542 |
| Threshold-only, `k=20` | 0.03714 | 0.1723 |
| Candidate raw, margin 8 | 0.000566 | 0.04788 |
| Candidate raw, margin 10 | 0 | 0.03136 |

These controls support the combined dimensionless reparameterization as the
reason complete-h1 edge transfer improves. Because they report edges rather than
GNK-specific texture, they do not establish which calibration component alone
would be sufficient for texture transmission.

## Interpretation for HYP-005

RUN-153 is a mixed but formally negative result for the full observable
protein-scale form of HYP-005. It establishes four narrow positive findings on
the fixed noiseless panel:

- a single global additive mutation scale and operating point generalize to all
  held-out core cells;
- complete-h1 GNK edge transmission is broadly practical;
- the intended dose response is monotone after numerically correct treatment of
  the registered equality boundary;
- observable texture RMS meets the registered length-stability bounds.

The primary texture requirement nevertheless fails. Aligned GNK-specific
plaquette texture is measurable, but only `24/360` actual-core cells meet the
registered prevalence floor, and the preregistered frozen positive control is
itself saturated at occupancy `0.99`. Exact normalization does not rescue the
result. RES-060 should therefore weaken HYP-005 at the full observable
texture-prevalence boundary while preserving the narrower calibration, edge,
dose, and length evidence.

This result does not validate the public AdaptEnv method, establish biological
fitness realism, test empirical outcomes or noise, or support accessibility,
adaptive-walk, active-learning, or optimizer claims. The experiment-local
candidate and Stage 4C/4D/4E hard stop remain unchanged.

## Limitations

- Confirmation uses a fixed, outcome-blind panel of 12 proteins, two per length
  bin. It is held out from calibration but is not a probability sample of
  proteins, so population-level generalization is not established.
- The `0.10` contrast floor and `0.25` prevalence threshold are preregistered
  noiseless engineering heuristics, not empirical detection limits.
- The frozen-core positive control is evaluated at reference occupancy `0.99`,
  making its interpretation intrinsically confounded by upper-shoulder
  saturation.
- The formal directionality failure is an implementation-level floating-point
  comparison artifact; the tolerance-aware sensitivity is decision-preserving.
- Checkpoint files contain non-standard JSON `NaN` literals.
- Steepness-only, threshold-only, and raw-margin controls are edge reports, not
  factorial texture controls.

## Produced artifacts

- `source_tables/run-153_summary.json`
- `source_tables/run-153_preflight.json`
- `source_tables/run-153_remote_sha256.json`
- `tables/metrics.json`
- `tables/run-153_collection_audit.csv`
- `tables/run-153_source_hash_audit.csv`
- `tables/run-153_acceptance_audit.csv`
- `tables/run-153_protocol_output_audit.csv`
- `tables/run-153_core_calibration_by_length.csv`
- `tables/run-153_heldout_core_by_length.csv`
- `tables/run-153_transmission_by_length.csv`
- `tables/run-153_directionality_by_length.csv`
- `tables/run-153_dose_curve.csv`
- `tables/run-153_scale_law_summary.csv`
- `tables/run-153_ablation_by_length.csv`

The paired notebook and `notebooks/analysis.py` invoke
`notebooks/derive_tables.py`. Rerunning that script regenerates every derived
table, verifies shared-storage SHA-256 parity and DVC object digests, reconstructs
the strict registered decision, and applies the recorded `1e-12`
directionality-boundary sensitivity analysis.
