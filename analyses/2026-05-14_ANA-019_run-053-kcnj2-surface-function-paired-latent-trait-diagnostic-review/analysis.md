# ANA-019: RUN-053 KCNJ2 surface-function paired latent trait diagnostic review

## Purpose

Assess whether the EXP-019 paired KCNJ2 surface/function diagnostic supports a two-trait state/function calibration over a function-only single-trait fit.

## Linked experiments/runs

- Experiments: `EXP-019`
- Runs: `RUN-053`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: python3

## Inputs

- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/summary.json`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/selected_panel.csv`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/mavenn_assay_metrics.csv`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/paired_readout_correlations.csv`
- `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053/branch_validations.csv`
- `experiments/2026-05-14_EXP-019_proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/config.yaml`
- `experiments/2026-05-14_EXP-019_proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/runs/RUN-053.yaml`

## Analysis performed

- Confirmed from the run metadata that `RUN-053` completed successfully on
  `lab-slurm` as scheduler job `196` with exit code `0`.
- Reviewed the `EXP-019` question: whether paired KCNJ2 surface trafficking
  plus ion-conduction function supports a state/function calibration better
  than a function-only single-trait fit.
- Verified that both KCNJ2 assays are single-mutant-only ProteinGym assays, so
  the run cannot test double-mutant ranking or epistasis recovery.
- Compared the three deterministic branches:
  `function_only_single_trait`,
  `paired_state_function_single_trait`, and
  `paired_state_function_explicit_two_trait`.
- Checked the paired-readout correlation, MAVE-NN observation-layer metrics,
  branch validation metrics, fitted parameter choices, and reference-to-peak
  geometry.

## Produced artifacts

- `tables/run-053_key_metrics.md`

## Main observations

- `RUN-053` completed and produced the declared outputs under
  `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053`.
- The durable output directory was DVC-tracked as
  `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053.dvc`.
- The empirical panel contained matched KCNJ2_MOUSE surface trafficking and
  ion-conduction function assays, both FACS single-mutant assays on a sequence
  of length `428`.
- The state/function readout correlation was weak but nonzero: Spearman
  `0.290` over `6789` matched variants.
- The assay-specific MAVE-NN observation layers were asymmetric, with test
  Spearman `0.538` for surface trafficking and only `0.389` for function. This
  makes the negative branch comparison less decisive than it would be with a
  strong function observation layer.
- The function-only branch had function holdout Spearman `0.066`, holdout NRMSE
  `1.000`, and functional KS `0.460`.
- The paired single-trait branch did not improve function recovery: holdout
  Spearman was `0.065`, holdout NRMSE was `1.002`, and functional KS worsened
  to `0.502`.
- The paired explicit two-trait branch had the lowest holdout NRMSE by a very
  small margin (`0.998`) but worse ranking and distributional recovery:
  holdout Spearman fell to `0.021`, and functional KS worsened to `0.511`.
- All branches fitted `epistasis_strength = 0.0`,
  `empirical_pairwise_strength = 0.0`, and `noise_amplitude = 0.0`.
- The reference-to-peak pathology remained: all branches placed the reference
  near the function peak, with function reference fraction of peak between
  `0.9932` and `0.9998`.
- The strongest supported inference is negative for this implementation:
  paired surface information did not improve the function recovery package
  relative to the function-only branch.

## Result records created

- `RES-019`

## Hypothesis updates

- `HYP-007` is weakened for the current KCNJ2 paired state/function
  implementation. The paired branches do not improve the prespecified
  single-mutant function recovery metrics as a package, and the explicit
  two-trait branch worsens ranking and KS.
- `HYP-001` is not directly updated because this is a focused paired-readout
  diagnostic, not a broad shared-regime empirical realism test.
