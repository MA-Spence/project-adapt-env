# RES-055: RUN-141 fixed support-aware probe fails core-GNK validation gates

## Summary

RUN-141 completed the design-matched, support-aware Stage 4B repair probe, but only 318/1620 monotonicity checks passed and only 34/1080 positive-epistasis core-vs-GNK separation gates passed; downstream overlap, joint-realism, and held-out reruns remain blocked for the current implementation.

## Generated from

- Analyses: ANA-055

## Relevant hypotheses

- Supports: None
- Weakens: HYP-005, HYP-004, HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/source_tables
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/metrics.json
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_output_collection_audit.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_support_template_audit.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_anchor_regime_audit.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_monotonicity_by_metric.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_monotonicity_by_support_anchor.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_monotonicity_by_dataset.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_core_gnk_gate_by_metric.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_core_gnk_gate_by_epi_weight.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_core_gnk_gate_by_dataset.csv
- analyses/2026-06-21_ANA-055_run-141-design-matched-support-aware-gnk-probe-review/tables/run-141_full_core_deltas_by_epi_weight.csv
- data/processed/stage-4b-design-matched-support-aware-gnk-probe/RUN-141/outputs.dvc
- data/processed/stage-4b-design-matched-support-aware-gnk-probe/RUN-141/logs.dvc
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-141/outputs

## Interpretation

`RUN-141` completed the single EXP-057 fix proposed in
`scratch/RES-054_assessment_and_path_to_finishing_component1.md#3`. The final
attempt computed all `135/135` parameter sets, emitted `35,640` generator
descriptor rows, wrote `8,910` support-aware graphlets, had empty stderr, and
exited with code `0`.

The probe itself is no longer the old measurement artifact. All generator
graphlets were support-aware (`support_aware_fraction = 1.0`), dimensions and
alphabet cardinalities were matched to empirical templates, and anchors covered
unfolded, fold-transition, and plateau regimes with small mean absolute target
margin errors (`0.0115`, `0.0151`, and `0.0137`).

The preregistered monotonicity gate nevertheless failed. Only `318/1620`
monotonicity checks passed (`0.1963`). Sign epistasis was the best individual
metric at `119/270` (`0.4407`), still below a majority. High-order fraction
passed only `36/270` checks (`0.1333`) with median Spearman `-0.1`; local-maxima
density passed `28/270` (`0.1037`) with median Spearman `-0.1817`.

The hard public-observable gate also failed. Across positive epistasis weights,
only `34/1080` core-vs-GNK separation rows passed (`0.0315`). All passing rows
were high-order-fraction rows. Sign-epistasis separation and local-maxima-density
separation had `0/1080` passes. Even at `epi_weight = 1.0`, only `16/270` rows
passed (`0.0593`), and the mean absolute full-minus-core high-order-fraction
delta was `0.00825`, below the preregistered `0.03` threshold.

The fixed support-aware probe therefore does not rescue the current public
additive-core plus GNK implementation. The result points to a core/GNK
scale-balance or readout-coupling problem in the implemented public observable,
not merely to the first-positions probe defect diagnosed in `RES-054`.

## Effect on hypothesis

- `HYP-005` is directly weakened for the current protein-scale public generator:
  increasing `epi_weight` did not yield broad monotone directional changes under
  the fixed support-aware measurement.
- `HYP-004` is weakened as an implemented public-generator realism claim because
  the local ruggedness, high-order, reachability, and autocorrelation controls
  do not move coherently enough to support biological-like landscape dynamics.
- `HYP-009` is weakened as a benchmark-construction route because routing and
  optimizer comparisons require matched, controllable ruggedness summaries
  before agent-visible routing can be interpreted.

This result does not refute every possible additive-core plus GNK model. It is a
negative result for the current EXP-057 grid and public observable.

## Limitations

- This is a Stage 4B generator probe, not a Stage 4C/4D/4E empirical overlap or
  held-out generalisation test.
- Empirical template descriptor rows were skipped because the empirical
  templates were near-complete but not mathematically complete. That does not
  affect the generator monotonicity gate, but it prevents interpreting this run
  as direct empirical descriptor matching.
- The result diagnoses the current amplitude/readout/grid implementation; it is
  not a formal impossibility result for all future scaling choices.

## Downstream use

- Do not run or cite downstream Stage 4C/4D/4E reruns for this implementation.
- Do not treat the public additive-core plus GNK composite as validated for
  empirical benchmark landscape construction.
- Cite `RUN-141` as the fixed-probe boundary result: after support-aware,
  design-matched measurement, the current implementation still fails the
  monotonicity and core-vs-GNK separation gates.
