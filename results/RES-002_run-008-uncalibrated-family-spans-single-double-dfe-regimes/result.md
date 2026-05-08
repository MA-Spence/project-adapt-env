# RES-002: RUN-008 shows the uncalibrated model family spans broad single- and double-mutant DFE regimes

## Summary

RUN-008 does not test empirical recovery directly, but it shows that the uncalibrated Adapt-Env family spans a broad envelope of one-step and two-step mutational regimes. Across the preregistered grid, singles were usually highly neutral, doubles were always more deleterious than singles and more lethal in most landscapes, stability margin dominated severity, and increasing epistasis strength increased double-mutant epistasis magnitude without changing single-mutant DFEs. This keeps HYP-001 scientifically live by arguing that the failure in RUN-006 is unlikely to be explained only by a complete lack of model-family capacity.

## Generated from

- Analyses: ANA-002

## Relevant hypotheses

- Supports: None
- Weakens: None
- Refutes: None
- Motivates: HYP-001

## Evidence

- `analyses/2026-05-08_ANA-002_run-008-uncalibrated-single-double-dfe-grid-review/tables/run-008_key_metrics.md`
- `data/processed/uncalibrated-dfe-grid-scan-single-double-mutants/RUN-008/summary.json`
- `data/processed/uncalibrated-dfe-grid-scan-single-double-mutants/RUN-008/single_double_dfe_per_setting.csv`
- `data/processed/uncalibrated-dfe-grid-scan-single-double-mutants/RUN-008/single_double_dfe_per_landscape.csv`
- `experiments/2026-05-08_EXP-002_uncalibrated-dfe-grid-scan-single-double-mutants/config.yaml`
- `experiments/2026-05-08_EXP-002_uncalibrated-dfe-grid-scan-single-double-mutants/runs/RUN-008.yaml`

## Interpretation

- `RUN-008` is an indirect scope check, not a realism test. It asks whether the uncalibrated Adapt-Env family can express a wide range of one-step and two-step mutational regimes before any empirical fitting.
- On that narrower question, the answer is yes. Across the full scan, singles remained mostly neutral but varied from `0.902` to `1.000` neutral fraction, while doubles ranged from almost never lethal to `32.8%` lethal and from `0.0%` to `4.0%` beneficial.
- `stability_margin` is the dominant harshness knob in this grid. Lower margins sharply increase double lethality and epistasis magnitude, while higher margins make the local neighborhood close to neutral.
- `epistasis_strength` behaves more specifically: it increases the magnitude and frequency of double-mutant epistasis without materially changing the single-mutant DFE summary.
- The scientific implication is limited but useful. `RES-001` showed that the current calibration setup does not recover external ProteinGym summary statistics well; `RUN-008` now suggests that this failure is unlikely to be explained only by a complete lack of local expressivity in the underlying model family.

## Effect on hypothesis

- This result motivates `HYP-001` by showing that the uncalibrated family can span a broad envelope of local DFE behavior worth testing against data.
- It does not directly support `HYP-001`, because the hypothesis is about recovering empirical assay statistics, and no empirical assay comparison is part of this run.

## Limitations

- The scan is centered on the model reference sequence and does not establish global landscape realism.
- Double-mutant summaries are based on `10,000` sampled doubles per landscape, not exhaustive enumeration.
- The run is unconditioned and uncalibrated, so it does not test whether real alignments or empirical fitting improve realism.
- Some controls are not cleanly monotonic, especially `n_functional_dims`, which means this scan is better at showing envelope width than at validating knob interpretability.

## Downstream use

- Use this result when interpreting `RES-001`: the current empirical failure is more likely a calibration or realism mismatch than total absence of local DFE capacity.
- Use `RES-002` as the record that `EXP-002` addressed a model-family capacity question rather than an external validation question.
