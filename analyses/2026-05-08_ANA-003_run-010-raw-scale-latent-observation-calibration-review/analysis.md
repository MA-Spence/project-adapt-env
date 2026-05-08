# ANA-003: RUN-010 Raw-Scale Latent Observation Calibration Review

## Purpose

Review RUN-010 against HYP-001 by testing whether coordinate-aware calibration on a homogeneous multi-mutant ProteinGym stability panel improves recovery of local landscape statistics.

## Linked experiments/runs

- Experiments: EXP-003
- Runs: RUN-010

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Inputs

- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/summary.json`
- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/selected_panel.csv`
- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/mavenn_assay_metrics.csv`
- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/per_assay_latent_fits.csv`
- `experiments/2026-05-08_EXP-003_proteingym-raw-scale-latent-observation-calibration-panel/config.yaml`
- `experiments/2026-05-08_EXP-003_proteingym-raw-scale-latent-observation-calibration-panel/runs/RUN-010.yaml`
- `results/RES-001_run-006-weakens-hyp-001-proteingym-distributional-statistics/result.md`
- `results/RES-002_run-008-uncalibrated-family-spans-single-double-dfe-regimes/result.md`

## Analysis performed

- Reviewed the completed Slurm run record and confirmed that the four declared durable outputs were recollected under `data/processed`.
- Inspected the selected ProteinGym subpanel for assay homogeneity, taxa, sequence-length range, and single- versus multiple-mutant coverage.
- Evaluated the assay-specific `mavenn` fits on raw assay scores to test whether the observation-model and coordinate-system confounds from `RUN-006` were materially reduced.
- Compared the shared raw-score branch, the shared latent-phenotype branch, and the per-assay latent-branch fits using single-mutant holdout, double-mutant holdout, epistasis-prediction, and KS metrics.
- Assessed whether the preregistered `EXP-003` prediction was met and what the result implies for the current wording of `HYP-001`.

## Outputs

- Figures: none
- Tables: `tables/run-010_key_metrics.md`

## Produced artifacts

- analyses/2026-05-08_ANA-003_run-010-raw-scale-latent-observation-calibration-review/tables/run-010_key_metrics.md

## Main observations

- `RUN-010` completed successfully and produced the four declared durable outputs: `summary.json`, `selected_panel.csv`, `mavenn_assay_metrics.csv`, and `per_assay_latent_fits.csv`.
- The selected panel is scientifically cleaner than `RUN-006`: `6` short `cDNA display proteolysis` stability assays from `Tsuboyama 2023`, spanning `4` taxa (`Human`, `Prokaryote`, `Eukaryote`, `Virus`), `23,279` measured variants total, `6,424` single mutants, `16,855` multiple mutants, and sequence lengths from `44` to `72` residues.
- The assay-specific `mavenn` models fit the raw assay scores well. Mean test Spearman was `0.894` across assays (range `0.852` to `0.946`), and mean test NRMSE was `0.491`. This indicates that the raw-score coordinate and nonlinear measurement-process layer were not the dominant bottleneck in this run.
- The shared raw-score Adapt-Env branch was better than the shared latent branch on the main ranking-based validation metrics, but still weak overall. Single-mutant holdout Spearman improved to `0.232` versus `0.084`, and epistasis-prediction Spearman improved to `0.500` versus `-0.016`, yet single-mutant holdout NRMSE remained `1.075`, double-mutant holdout Spearman remained near zero at `-0.017`, and double-mutant holdout NRMSE remained `1.114`.
- The shared latent branch reverted toward a near-peak reference solution and lost most pairwise structure. It fit `n_functional_dims = 6`, `peak_distance_from_consensus = 0`, `functional_reference_fraction_of_peak = 0.897`, `empirical_pairwise_strength = 0.0`, and `epistasis_strength = 0.0`, with weak holdout metrics.
- Per-assay latent fits did not provide a convincing ceiling. Mean single-mutant holdout Spearman was only `0.136` across the `5` assays with non-null values, mean single-mutant holdout NRMSE was `1.219`, mean double-mutant holdout Spearman was `0.116`, and all `6` per-assay fits set `epistasis_strength = 0.0`.
- Only one per-assay latent fit showed moderate double-mutant rank recovery: `SPTN1_CHICK_Tsuboyama_2023_1TUD` reached single-mutant holdout Spearman `0.333` and double-mutant holdout Spearman `0.442`. The rest remained weak or near null.
- The preregistered `EXP-003` prediction was therefore not met. The measurement-layer confound was reduced, and the shared raw branch was modestly better than the shared latent branch, but the downstream Adapt-Env calibration still failed to recover local single- and double-mutant statistics cleanly.
- On the current evidence, `RUN-010` further weakens `HYP-001`. It removes the strongest coordinate-system confound left by `RUN-006`, yet the main calibration failure persists on a homogeneous stability-only panel with real multi-mutant signal.

## Result records created

- `RES-003`

## Hypothesis updates

- `HYP-001` is further weakened, not refuted. The result narrows the likely explanation for earlier failure away from assay-wise z-scoring and synthetic wild-type zeroing alone.
