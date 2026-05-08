# ANA-002: RUN-008 Uncalibrated Single and Double DFE Grid Review

## Purpose

Review RUN-008 as an indirect scope check under HYP-001 by summarizing whether the uncalibrated Adapt-Env model family spans plausible one-step and two-step mutational effect regimes before any empirical DMS fitting.

## Linked experiments/runs

- Experiments: EXP-002
- Runs: RUN-008

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Inputs

- `data/processed/uncalibrated-dfe-grid-scan-single-double-mutants/RUN-008/summary.json`
- `data/processed/uncalibrated-dfe-grid-scan-single-double-mutants/RUN-008/single_double_dfe_per_setting.csv`
- `data/processed/uncalibrated-dfe-grid-scan-single-double-mutants/RUN-008/single_double_dfe_per_landscape.csv`
- `experiments/2026-05-08_EXP-002_uncalibrated-dfe-grid-scan-single-double-mutants/config.yaml`
- `experiments/2026-05-08_EXP-002_uncalibrated-dfe-grid-scan-single-double-mutants/runs/RUN-008.yaml`
- `results/RES-001_run-006-weakens-hyp-001-proteingym-distributional-statistics/result.md`

## Analysis performed

- Reviewed the completed Slurm run record to confirm execution status, promoted outputs, and declared parameter grid.
- Inspected the top-level run summary to quantify how often doubles were more deleterious or more lethal than singles across the full scan.
- Compared per-setting summaries across `stability_margin`, `functional_sigma_base`, `n_functional_dims`, and `epistasis_strength` to identify which controls most strongly shape local DFE behavior.
- Assessed whether the observed parameter envelope is broad enough to matter scientifically for `HYP-001`, while separating this model-capacity question from direct empirical realism.
- Interpreted the result relative to `RES-001`, which previously weakened `HYP-001` using external ProteinGym comparisons.

## Outputs

- Figures: none
- Tables: `tables/run-008_key_metrics.md`

## Produced artifacts

- analyses/2026-05-08_ANA-002_run-008-uncalibrated-single-double-dfe-grid-review/tables/run-008_key_metrics.md

## Main observations

- `RUN-008` completed successfully and produced the three declared durable outputs: `summary.json`, `single_double_dfe_per_setting.csv`, and `single_double_dfe_per_landscape.csv`.
- The experiment evaluated `81` grid settings across `4` seeds each, for `324` landscapes total, with sequence length `80` and `10,000` sampled double mutants per landscape.
- Across the full scan, doubles were more deleterious than singles in `100%` of settings and more lethal than singles in `83.3%` of settings.
- The model family spans a wide local regime envelope. Across settings, the single-mutant neutral fraction ranged from `0.902` to `1.000`, the single-mutant lethal fraction from `0.000` to `0.009`, the double-mutant lethal fraction from `0.000075` to `0.328`, the double-mutant beneficial fraction from `0.000` to `0.040`, and the mean absolute epistasis magnitude from `0.0095` to `3.126`.
- `stability_margin` was the dominant harshness control. Averaged across the other grid axes, decreasing `stability_margin` from `12` to `5` increased double lethality from `0.000075` to `0.328`, reduced double neutrality from `0.954` to `0.611`, and increased mean absolute epistasis from `0.050` to `3.081`.
- `functional_sigma_base` mainly changed the single-mutant regime rather than the double-mutant lethality baseline. At `12.0`, the mean single-mutant neutral fraction fell to `0.959` and the beneficial fraction rose to `0.0065`; at `20.0` and `28.0`, singles were almost entirely neutral and beneficial singles disappeared at this summary level.
- `epistasis_strength` had almost no effect on the single-mutant DFE, but it did increase double-mutant epistasis. The mean absolute epistasis rose from `1.115` at `0.0` to `1.143` at `0.06`, while the positive-epistasis fraction rose from `0.00043` to `0.0122`.
- `n_functional_dims` did not behave as a clean monotonic realism knob in this scan. The `4`-dimensional setting was the least neutral and most double-beneficial on average, while `2` and `8` dimensions were both more neutral.
- The strongest mean epistasis setting was `stability_margin-5__functional_sigma_base-28__n_functional_dims-2__epistasis_strength-0.06`, and the weakest was `stability_margin-12__functional_sigma_base-28__n_functional_dims-8__epistasis_strength-0`.
- This run does not test the central empirical claim in `HYP-001`, because no DMS data are involved. It does, however, show that the uncalibrated Adapt-Env family can generate a broad range of one-step and two-step mutational regimes, so the negative result in `RES-001` is unlikely to be explained only by complete lack of local DFE expressivity.

## Result records created

- `RES-002`

## Hypothesis updates

- `RES-002` motivates `HYP-001` but does not directly support, weaken, or refute it. The run is informative about model-family capacity, not external realism.
