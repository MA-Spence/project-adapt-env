# ANA-018: RUN-051 OXDA expression-activity paired latent trait diagnostic review

## Purpose

Assess whether the EXP-018 paired OXDA expression/activity diagnostic supports a two-trait state/function calibration over a function-only single-trait fit.

## Linked experiments/runs

- Experiments: `EXP-018`
- Runs: `RUN-051`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: python3

## Inputs

- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/summary.json`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/selected_panel.csv`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/mavenn_assay_metrics.csv`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/paired_readout_correlations.csv`
- `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051/branch_validations.csv`
- `experiments/2026-05-14_EXP-018_proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/config.yaml`
- `experiments/2026-05-14_EXP-018_proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/runs/RUN-051.yaml`

## Analysis performed

- Confirmed from the run metadata that `RUN-051` completed successfully on
  `lab-slurm` as scheduler job `194` with exit code `0`.
- Confirmed that `RUN-051` is the corrected live run for `EXP-018`, after the
  earlier `RUN-046` single-mutant MAVE-NN configuration failure.
- Reviewed the `EXP-018` question: whether paired OXDA expression plus enzyme
  activity supports a state/function calibration better than the function-only
  single-trait fit.
- Verified that both OXDA assays are single-mutant-only ProteinGym assays, so
  the run cannot test double-mutant ranking or epistasis recovery.
- Compared the three deterministic branches:
  `function_only_single_trait`,
  `paired_state_function_single_trait`, and
  `paired_state_function_explicit_two_trait`.
- Checked the paired-readout correlation, MAVE-NN observation-layer metrics,
  branch validation metrics, fitted parameter choices, and reference-to-peak
  geometry.

## Produced artifacts

- `tables/run-051_key_metrics.md`

## Main observations

- `RUN-051` completed and produced the declared outputs under
  `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051`.
- The durable output directory was DVC-tracked as
  `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051.dvc`.
- The empirical panel contained matched OXDA_RHOTO expression and activity
  assays, both FACS single-mutant assays on a sequence of length `364`.
- The state/function readout correlation was weak but nonzero:
  Spearman `0.246` over `6387` matched variants.
- The assay-specific MAVE-NN observation layers were only moderate, with test
  Spearman `0.567` for expression and `0.542` for activity. This limits the
  strength of any negative inference.
- The function-only branch had activity holdout Spearman `0.183`, holdout NRMSE
  `0.981`, and functional KS `0.520`.
- The paired single-trait branch was worse on activity holdout ranking
  (`0.061`) and did not improve KS (`0.523`).
- The paired explicit two-trait branch was nearly tied with the function-only
  branch on activity holdout ranking (`0.180` versus `0.183`) but did not
  improve NRMSE or KS.
- All branches fitted `epistasis_strength = 0.0`,
  `empirical_pairwise_strength = 0.0`, and `noise_amplitude = 0.0`.
- The reference-to-peak pathology remained: all branches placed the reference
  effectively at the peak, with function reference fraction of peak between
  `0.9995` and `0.9999`.
- The strongest supported inference is negative for this implementation:
  paired expression information did not improve the activity recovery package
  relative to the function-only branch.

## Result records created

- `RES-018`

## Hypothesis updates

- `HYP-007` is weakened for the current OXDA paired state/function
  implementation. The branch that actually instantiates the explicit
  state/function split does not improve the prespecified single-mutant activity
  recovery metrics.
- `HYP-001` is not directly updated because this is a focused paired-readout
  diagnostic, not a broad shared-regime empirical realism test.
