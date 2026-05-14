# ANA-015: RUN-041 PHOT flexible monotone two-latent readout diagnostic review

## Purpose

Assess whether the `EXP-015` flexible monotone two-latent PHOT diagnostic
improves held-out ranking and epistasis recovery relative to the explicit
product/gate two-trait failure and the stability-targeted control.

## Linked experiments/runs

- Experiments: `EXP-015`
- Runs: `RUN-041`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Kernel: `python3`

## Inputs

- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/summary.json`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/selected_panel.csv`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/model_validations.csv`
- `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041/progress.json`
- `experiments/2026-05-14_EXP-015_proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/config.yaml`
- `experiments/2026-05-14_EXP-015_proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/runs/RUN-041.yaml`
- `scripts/proteingym_phot_structural_mismatch_diagnostic.py`
- `results/RES-013_run-039-stability-readout-validation-objective-improves-phot-recovery/metrics.json`
- `results/RES-014_run-040-explicit-two-trait-replicate-reproduces-phot-validation-objective-failure/metrics.json`

## Analysis performed

- Confirmed from the run metadata that `RUN-041` completed successfully on
  `lab-slurm` as scheduler job `174` with exit code `0`.
- Confirmed that `RUN-041` is the collected run for `EXP-015`, which asks
  whether replacing the hard-coded `stability_gate * readout_capacity` collapse
  with a flexible monotone stability/activity readout surface improves
  `PHOT_CHLRE` recovery.
- Verified from `config.yaml` and the runner that this is a lightweight
  diagnostic rather than an AdaptEnv SMC posterior. The fitted model uses:
  - an alignment-frequency stability proxy
  - an additive ridge activity score fit from sequence mutations
  - a monotone two-dimensional binned surface over the stability and activity
    coordinates
  - held-out single- and double-mutant splits by deterministic hash modulo
- Reviewed the MAVE-NN assay diagnostics and selected panel to confirm that the
  same `PHOT_CHLRE_Chen_2023` activity/FACS assay was used.
- Compared `model_validations.csv` to the best particles from `RES-014`
  (`RUN-040`, explicit product/gate two-trait replicate) and `RES-013`
  (`RUN-039`, stability-targeted validation-objective control).
- Recorded that reference-to-peak geometry and functional KS are not available
  for `RUN-041` because this diagnostic predicts observed assay scores rather
  than constructing a full AdaptEnv landscape with a searched peak.

## Outputs

- Figures: none
- Tables: `tables/run-041_key_metrics.md`

## Produced artifacts

- `analyses/2026-05-14_ANA-015_run-041-phot-flexible-monotone-two-latent-readout-diagnostic-review/tables/run-041_key_metrics.md`

## Main observations

- `RUN-041` completed in `18` seconds and produced the declared lightweight
  diagnostic outputs under
  `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041`.
- The empirical target remained `PHOT_CHLRE_Chen_2023`, an activity/FACS assay
  with sequence length `118`, `167,529` measured variants, `2,122` single
  mutants, and `165,407` multiple mutants.
- The MAVE-NN observation layer was moderate and slightly weaker than the
  immediately preceding PHOT runs: test Spearman `0.657` and test NRMSE
  `0.947`.
- The diagnostic trained on `15,997` variants after applying the configured
  mutation-count and sampling limits, and evaluated `404` held-out single
  mutants plus `43` held-out double mutants.
- The flexible monotone stability/activity surface reached single-mutant
  holdout Spearman `0.473` and double-mutant holdout Spearman `0.807`. This is
  a large ranking improvement over the explicit product/gate two-trait `RUN-040`
  result (`+0.275` single, `+0.555` double).
- The same model did not recover epistasis. Epistasis-prediction Spearman was
  `0.035`, and the epistasis KS statistic was `0.814`. Relative to `RUN-040`,
  epistasis Spearman fell by `0.230`, and epistasis KS worsened by `0.348`.
- Relative to the stability-targeted `RUN-039` control, the flexible two-latent
  diagnostic was close on single-mutant ranking (`-0.033`), slightly better on
  double-mutant ranking (`+0.024`), but much worse on epistasis-prediction
  Spearman (`-0.600`) and epistasis KS (`+0.064`).
- The result therefore supports a narrower diagnosis: the forced product/gate
  readout likely contributed to the ranking failure in `RUN-036/RUN-040`, but
  simply allowing a flexible monotone two-latent surface is not enough to
  recover PHOT epistasis.
- Because `RUN-041` has no full synthetic landscape, it cannot answer whether
  the flexible surface fixes the near-peak reference artifact. That part of the
  `EXP-015` pre-experiment prediction remains untested by this lightweight
  diagnostic output.

## Result records created

- `RES-015`

## Hypothesis updates

- `HYP-007` is weakened at the level of a full empirical recovery claim. The
  flexible two-latent readout improved held-out ranking but did not improve the
  complete primary readout package, especially epistasis prediction.
- The result does not refute all multiple-latent causal maps. It narrows the
  mismatch: readout flexibility helps ranking, but PHOT epistasis still requires
  additional structure or a different causal decomposition.
