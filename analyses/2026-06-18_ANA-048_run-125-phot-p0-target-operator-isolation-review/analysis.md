# ANA-048: RUN-125 PHOT P0 target-operator isolation review

## Purpose

Collect and interpret the `RUN-125` outputs from `EXP-049`, which tested the
Phase-0 isolation plan from
`scratch/RES-047_miscalibration_root_cause_and_isolation_plan.md`.

The analysis asks three preregistered questions:

- P0.1: does the empirical PHOT target survive graphlet bootstrap and MAVE-NN
  refit perturbation?
- P0.2: does a known injected spectrum survive the identical empirical
  measurement pipeline?
- P0.3: do the sparse-Walsh generator supports and sampled
  `dense_descriptive` graphlets share the same position pool?

## Linked experiments/runs

- Experiments: `EXP-049`
- Runs: `RUN-125`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Derivation script: `notebooks/derive_tables.py`
- Kernel: `python3`

## Produced artifacts

- `source_tables/run-125_summary.json`
- `source_tables/run-125_preflight.json`
- `source_tables/run-125_graphlet_catalog.csv`
- `source_tables/run-125_p0_1_target_graphlet_metrics.csv`
- `source_tables/run-125_p0_1_target_bootstrap.csv`
- `source_tables/run-125_p0_1_mavenn_refit_metrics.csv`
- `source_tables/run-125_p0_1_replicate_audit.csv`
- `source_tables/run-125_p0_2_injected_terms.csv`
- `source_tables/run-125_p0_2_roundtrip_graphlet_metrics.csv`
- `source_tables/run-125_p0_2_roundtrip_transfer.csv`
- `source_tables/run-125_p0_2_mavenn_metrics.csv`
- `source_tables/run-125_p0_3_position_pool_audit.csv`
- `source_tables/run-125_p0_3_support_terms.csv`
- `source_tables/run-125_p0_3_support_containment_summary.csv`
- `tables/metrics.json`
- `tables/run-125_output_collection_audit.md`
- `tables/run-125_graphlet_catalog_summary.md`
- `tables/run-125_p0_1_target_stability.md`
- `tables/run-125_p0_1_mavenn_refit_fit_quality.md`
- `tables/run-125_p0_2_profile_summary.md`
- `tables/run-125_p0_2_roundtrip_transfer_key_metrics.md`
- `tables/run-125_p0_3_support_alignment.md`
- `tables/run-125_phase_decisions.md`

## Collection policy

`RUN-125` completed on `lab-slurm` as Slurm job `1933` with inferred scheduler
exit code `0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-125/outputs`

The full RUN-125 compact output bundle was small enough to copy into this
analysis: `1,209,603` bytes across `14` source files. The source payload is
audited with SHA-256 checksums in
`tables/run-125_output_collection_audit.md`. No heavy coefficient table or
alignment artifact was produced by this diagnostic.

## Output audit

The run used `167,529` PHOT assay rows from the frozen `RUN-077` source,
selected `24` complete graphlets across dimensions `4`, `5`, and `6`, and was
not a quick run. It produced:

- `10` P0.1 MAVE-NN refits: five seeds across two regularisation profiles.
- `36` P0.2 transfer rows: three known-spectrum profiles across the selected
  summary metrics.
- `2,048` P0.3 support-term rows for the audited order-3 sparse-Walsh support
  distribution.

The MAVE-NN paths did not use the quick surrogate. P0.1 refit test Spearman
values ranged from roughly `0.63` to `0.77` for the low-regularisation profile
and `0.64` to `0.70` for the high-regularisation profile, so the refits are not
just failed fits. The P0.2 synthetic fits had test Spearman values from `0.56`
to `0.73`, so the round-trip readout is interpretable but imperfect.

## P0.1 Target Stability

P0.1 does not recover a robust high-order or multi-peak empirical target under
this MAVE-NN latent operator. The cached `RUN-077` `latent_phi` target has mean
Walsh high-order fraction `3.88e-13` and mean local-maxima count `1.0` across
the 24 selected graphlets. Across all 10 MAVE-NN refits, the largest mean
Walsh high-order fraction was only `1.06e-12`, and local maxima remained
exactly `1.0`.

Accessibility is the unstable component. The cached target mean was `0.5417`
with bootstrap interval `0.3750` to `0.7083`, but refit means ranged from
`0.0` to `0.5417`. Thus the target as measured here is not a stable
high-order/multi-peak object; it is a smooth/additive-like MAVE-NN latent
surface with refit-sensitive accessibility.

The prepared PHOT assay contained no replicate-valued columns, so replicate
propagation could not be performed.

## P0.2 Round-Trip Measurement Operator

P0.2 is the most diagnostic part of the run. Known injected landscapes do not
survive the empirical measurement pipeline as spectra or local epistasis
objects.

For the `pure_order3` profile, direct graphlet measurement of the injected
landscape had mean Walsh high-order fraction only `0.0139`, mean local maxima
`2.0417`, mean accessibility `0.2033`, and magnitude-epistasis-to-single-mutant
SD ratio `0.6498`. After MAVE-NN recovery, the same profile had mean
high-order fraction `3.26e-15`, local maxima `1.0`, accessibility `1.0`, local
sign epistasis `0.0`, and magnitude ratio `4.00e-7`.

The mixed order-2/order-3 profile showed the same smoothing: recovered
high-order fraction `4.15e-15`, local maxima `1.0`, local sign epistasis `0.0`,
and magnitude ratio `2.73e-7`. The low-order control also recovered as
unimodal and nearly epistasis-free.

This is strong evidence that the empirical measurement operator is
non-identity and strongly smoothing. It does not support the scratch note's
specific predicted direction that injected order-3 should return as high
high-order fraction near `0.5+`. Instead, the measured transfer function goes
the other way: global injected order-3 terms are already weakly visible in the
local complete-graphlet summaries, and the MAVE-NN latent recovery erases the
remaining local epistasis almost completely.

## P0.3 Support-Sampler Alignment

P0.3 rejects low support containment as the limiting explanation for the
RUN-123/RES-047 local-spectrum mismatch under this audited support rule.

The dense graphlet position pool and generator cluster position pool both
contained `13` positions, and all `13` positions were shared. All `2,048`
generated order-3 supports were fully contained in the graphlet position pool
and were contained in at least one sampled graphlet. Calibration containment
was `0.9521`, test containment was `0.9292`, and each support appeared in
`8.76` sampled graphlets on average.

Thus the support and sampler pools are aligned for this run. A low containment
fraction does not explain the PHOT mis-read here.

## Interpretation

RUN-125 supports the broad Layer-2 diagnosis from the scratch note: the PHOT
calibration target and model outputs are not being measured on a common
footing. It does not support every mechanism proposed in that note.

The strongest result is P0.2. Known injected landscapes with nonzero local
epistasis and multiple local maxima recover through MAVE-NN as effectively
additive, unimodal latent landscapes. That makes the previous PHOT calibration
objective suspect: a model can appear to fail or pass depending on whether the
summary is measured directly, through local complete graphlets, or after
MAVE-NN latent recovery.

P0.1 reinforces the problem from the empirical side. Under the operator used
here, cached and refit MAVE-NN `latent_phi` summaries do not look like the
high-order/multi-peak target implied by the preceding PHOT calibration series.
High-order fraction is effectively zero and local maxima remain fixed at one.
Accessibility is refit-sensitive, so it should not be treated as a sharp scalar
target without uncertainty bands.

P0.3 narrows the diagnosis by removing one proposed failure channel. The
generated sparse-Walsh supports are not invisible to the sampled graphlets;
they are well contained. The immediate failure mode is therefore not a
support-pool mismatch in this audited implementation.

## Effect on hypotheses

For `HYP-004`, RUN-125 weakens the current PHOT adequacy route, but for a
different reason than another failed generator fit would. The issue is not
simply that the latest generator lacks capacity. The target/measurement layer
itself is not preserving the high-order, local-optima, and local-epistasis
objects that the calibration workflow is trying to tune.

For `HYP-009`, RUN-125 is methodological motivation. It does not test PLM or
optimizer agents and does not produce a matched benchmark landscape. It does,
however, provide a reason to recalibrate future benchmark-facing targets in
the recovered coordinate or to replace the measurement operator before using
PHOT graphlet summaries as optimization objectives.

## Limitations

- This is still PHOT-only. It should not be generalized to all empirical
  proteins without a second assay.
- P0.2 used the implemented synthetic spectrum profiles and graphlet operator;
  the direct local graphlet high-order fraction of the injected landscapes was
  itself low, so the round-trip result diagnoses the full local-operator plus
  MAVE-NN pipeline rather than MAVE-NN alone.
- MAVE-NN recovery quality was moderate but not perfect, especially for the
  low-order control. The qualitative conclusion is not that every recovered
  value is exact, but that high-order fraction, local maxima, local sign
  epistasis, and magnitude epistasis were collapsed by orders of magnitude.
- The run did not test a repaired generator or a corrected target objective.

## Decision

Use RUN-125 as a measurement-operator isolation result, not as a generator
equivalence result. The next PHOT calibration target should not use the old
scalar high-order/local-ruggedness objective without either correcting the
P0.2 transfer function or replacing the measurement operator. Low support
containment should not be cited as the current bottleneck, because P0.3 found
full containment under the audited support rule.
