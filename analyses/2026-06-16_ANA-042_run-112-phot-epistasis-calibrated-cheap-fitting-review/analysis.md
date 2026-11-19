# ANA-042: RUN-112 PHOT epistasis-calibrated cheap fitting review

## Purpose

Review compact RUN-112 outputs to determine whether positive affine observable calibration plus epistasis-amplitude objective makes the cheap PHOT nested-min fitting protocol behave similarly to the prior full ABC fitted-family comparator while preserving graphlet, accessibility, local-energy, raw latent epistasis, and module-locality diagnostics.

## Linked experiments/runs

- Experiments: EXP-043
- Runs: RUN-112

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-112_summary.json`
- `source_tables/run-112_preflight.json`
- `source_tables/run-112_fit_protocol_equivalence.csv`
- `source_tables/run-112_fit_protocol_distances.csv`
- `source_tables/run-112_selected_fit_protocol_candidates.csv`
- `source_tables/run-112_abc_comparator_selected_family_distances.csv`
- `source_tables/run-112_full_abc_comparator_selected_particles.csv`
- `source_tables/run-112_graphlet_summary.csv`
- `source_tables/run-112_module_energy_summary.csv`
- `source_tables/run-112_local_abc_rounds.csv`
- `source_tables/run-112_graphlet_catalog.csv`
- `source_tables/run-112_remote_output_manifest.tsv`
- `tables/metrics.json`
- `tables/run-112_output_collection_audit.md`
- `tables/run-112_equivalence_by_method_reference.md`
- `tables/run-112_same_family_equivalence.md`
- `tables/run-112_calibrated_test_distances.md`
- `tables/run-112_original_abc_comparator_test_distances.md`
- `tables/run-112_top_component_failures.md`
- `tables/run-112_observable_vs_raw_epistasis_and_locality.md`
- `tables/run-112_best_rows_landscape_deltas.md`
- `tables/run-112_best_epistasis_amplitude_rows.md`
- `tables/run-112_empirical_module_locality_test.md`
- `tables/run-112_local_abc_rounds.md`

## Collection policy

`RUN-112` completed on `lab-slurm` as Slurm job `1907` with inferred exit code
`0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-112/outputs`

This analysis copied only compact JSON/CSV outputs needed to interpret the
experiment. The largest detailed files remain on the runtime volume:
`gft_coefficients.csv` (`30.4 MB`), `walsh_terms.csv` (`29.9 MB`), and
`graphlet_metrics.csv` (`1.9 MB`). The copied payload is `227,247` bytes and
is audited in `tables/run-112_output_collection_audit.md`.

## Output audit

`RUN-112` evaluated `24` PHOT graphlets, `730` candidate parameter sets, and
selected `20` candidates. The run produced `40` equivalence rows across four
cheap fitting methods, five nested families, and two comparator references.

The two comparator notions must be kept separate:

- `same_family_full_abc` compares each cheap selected family to the selected
  `RUN-098` full-ABC particle from that same nested family after re-evaluation
  under the calibrated EXP-043 reporting path.
- `best_any_full_abc` compares each cheap selected family to the best available
  recomputed full-ABC comparator row under that same calibrated path.

The original `RUN-098` comparator distance table is copied only as contextual
provenance. It is not on the same calibrated-equivalence scale as
`fit_protocol_equivalence.csv`.

## Main result

No cheap protocol passed the preregistered operational equivalence decision.
`RUN-112` reports `0/40` operationally equivalent rows.

This is not because all aggregate distances were worse than the comparator.
Under the absolute and relative distance margins alone,
`exact_descriptive_gated`, `local_abc_after_cheap_fit`, and
`surrogate_assisted` each passed `5/5` same-family comparisons, while
`deterministic_summary_loss` passed `4/5`. The failure came from component
guardrails: every method had `0/5` same-family component-guardrail passes.

The dominant failure component was `signed_pairwise_epistasis_mean`. The
empirical test-graphlet signed mean is essentially zero, so small absolute
signed offsets become very large standardized component distances. This means
the EXP-043 objective is numerically stricter than the EXP-042 objective and
should not be read as a simple repeat of the same equivalence test.

The best calibrated-objective selected row was
`exact_descriptive_gated` / `uncoupled_module`, with test distance `1.1474`.
It still violated `magnitude_epistasis_mean` and
`signed_pairwise_epistasis_mean`, had calibrated magnitude epistasis mean
`0.0156` versus empirical `0.2926`, raw latent magnitude epistasis mean
`29.10`, and accessible-path fraction `0.0` versus empirical `0.1975`.

The row closest on calibrated magnitude epistasis was
`exact_descriptive_gated` / `empirical_module`, but it still only reached
`0.1774` versus empirical `0.2926`, retained raw latent magnitude epistasis
`131.9`, and had test distance `12.0043`. Thus the readout calibration reduced
the visible scale of epistasis relative to the raw latent values, but it did
not solve either the observable epistasis-amplitude match or the raw latent
amplitude pathology.

## Landscape diagnostics

Local edge-energy concentration improved for some best-distance rows. The best
calibrated row had test edge top-10 energy fraction `0.3955`, close to the
empirical `0.4060`. However, it missed other biologically important summaries:
accessibility collapsed to zero, site top-1 energy was too low (`0.2840`
versus empirical `0.4162`), and calibrated magnitude epistasis was much too
small.

The empirical-module locality warning also remains. On test graphlets,
empirical-module rows have high edge-module effective counts and entropy. The
selected `RUN-098` empirical-module comparator has `11.72/12` effective edge
modules and edge normalized entropy `0.995`, while the EXP-043 cheap
empirical-module rows remain similarly diffuse. This is not residue-level
module recovery.

## Decision

`RUN-112` does not support using the EXP-043 calibrated cheap-fitting protocol
as an equivalent replacement for the previous full ABC fit. It weakens the
specific idea that a simple positive affine observable readout plus epistasis
amplitude objective is enough to make the current nested-min PHOT model
usable for observable-scale accessibility, neutral-network, and epistasis
analyses.

The result preserves one narrower lesson: adding explicit calibration changes
the numerical scale and can improve selected local-energy summaries, but the
underlying latent epistasis remains badly miscalibrated and the selected
observable landscapes still fail accessibility and component guardrails.
