# ANA-056: RUN-129 PHOT non-additive measurement-operator certification review

## Purpose

Collect and interpret `RUN-129` outputs from `EXP-050`, the PHOT non-additive
MAVE-NN measurement-operator certification diagnostic.

The analysis asks whether the two preregistered non-additive operator profiles
can recover injected high-order and local-ruggedness spectra under the same
dense complete-graphlet measurement path used for PHOT.

## Linked experiments/runs

- Experiments: `EXP-050`
- Runs: `RUN-129`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Derivation script: `notebooks/derive_tables.py`
- Kernel: `python3`

## Produced artifacts

- `source_tables/run-129_summary.json`
- `source_tables/run-129_preflight.json`
- `source_tables/run-129_graphlet_catalog.csv`
- `source_tables/run-129_operator_profile_catalog.csv`
- `source_tables/run-129_injected_terms.csv`
- `source_tables/run-129_roundtrip_graphlet_metrics.csv`
- `source_tables/run-129_roundtrip_transfer.csv`
- `source_tables/run-129_mavenn_metrics.csv`
- `source_tables/run-129_operator_certification_summary.csv`
- `tables/metrics.json`
- `tables/run-129_output_collection_audit.md`
- `tables/run-129_operator_certification_status.md`
- `tables/run-129_fit_failure_modes.md`
- `tables/run-129_certification_metric_matrix.md`
- `tables/run-129_metric_pass_summary.md`
- `tables/run-129_graphlet_catalog_summary.md`
- `tables/run-129_injected_spectrum_summary.md`

## Collection policy

`RUN-129` completed on `lab-slurm` as Slurm job `1937` with inferred scheduler
exit code `0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-129/outputs`

The compact output and log payload was copied into durable processed data:

- `data/processed/phot-nonadditive-measurement-operator-certification/RUN-129/outputs`
- `data/processed/phot-nonadditive-measurement-operator-certification/RUN-129/logs`

These directories are DVC-tracked as `outputs.dvc` and `logs.dvc`. The
analysis copied the source output tables into `source_tables/` and records
SHA-256 checksums in `tables/run-129_output_collection_audit.md`. The collected
payload is small: `237,462` bytes of output tables and `719` bytes of logs.

## Output audit

The run used the frozen `RUN-077` PHOT source output and recorded
`167,529` assay rows. It selected `24` dense complete graphlets and evaluated
three injected spectrum profiles:

- `pure_order3`
- `order2_order3_mixture`
- `low_order_control`

The injected landscape payload contains `176` terms. The certification gate
required four metrics to recover within the preregistered tolerance on all
three spectra for each operator:

- `walsh_high_order_fraction`
- `sign_epistasis_fraction_local`
- `n_local_maxima`
- `accessible_path_fraction`

The two tested operator profiles were `pairwise_low_regularization` and
`blackbox_ge_neural`.

## Certification result

No non-additive measurement operator was certified.

Across the two operator profiles, `RUN-129` reported `24` required
operator-spectrum-metric checks and `0` passing checks. There were `0`
certification-metric rows with recovered values because all `6/6` fit-attempt
rows failed before a recovered landscape was available. Therefore the result is
not a quantitative statement that non-additive MAVE-NN operators smooth
epistasis. It is a certification failure of the current non-additive
measurement implementation.

The failure modes were distinct:

- `pairwise_low_regularization`: all `3/3` spectra were preflight-skipped with
  `pairwise_gpmap_preflight_skip: sequence_length=118 exceeds max_pairwise_sequence_length=80`.
- `blackbox_ge_neural`: all `3/3` spectra failed with a
  `MultilayerPerceptronGPMap` / Keras `build()` signature error in the installed
  MAVE-NN stack.

The Slurm job itself did not fail. The log contains TensorFlow CPU/GPU
informational messages and `labproj-exit-code.txt` records `0`.

## Controls and interpretation

The injected-spectrum controls were present, and the direct graphlet metrics
were recorded. Those controls establish what the recovery task should have
measured, but they do not rescue the certification because neither
non-additive operator produced recovered values.

This matters scientifically because `EXP-050` was intended to repair the
measurement-operator problem exposed by `RES-048`. It does not do that. There
is no certified non-additive PHOT measurement operator from this run, and no
usable transfer function from these two non-additive profiles can be applied
symmetrically to downstream empirical/model comparisons.

The result should be interpreted narrowly. It weakens the current PHOT
measurement-operator route and the current benchmark-construction route that
depends on that route, but it does not refute biological PHOT ruggedness, the
additive-core + GNK forward generator, or optimizer/PLM-agent discrimination.

`RUN-126` and `RUN-128` are execution context only: `RUN-126` failed by OOM
before the corrected preflight behavior, and `RUN-128` was canceled before
execution. The analysis decision is based on the completed corrected run,
`RUN-129`.

## Decision

Do not cite `EXP-050` as evidence that non-additive MAVE-NN operators recover
known high-order/local-ruggedness spectra under the PHOT dense graphlet
pipeline. Cite it as a negative operator-certification result for the current
installed non-additive MAVE-NN route.
