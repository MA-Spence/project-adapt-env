# ANA-004: RUN-012 predictive-joint calibration and pairwise-capacity review

## Purpose

Assess whether predictive-joint calibration and reduced pairwise compression improve empirical reconstruction on the ProteinGym stability panel tested in EXP-004.

## Linked experiments/runs

- Experiments: EXP-004
- Runs: RUN-012

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Inputs

- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/summary.json`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/selected_panel.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/mavenn_assay_metrics.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/branch_validations.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/per_assay_branch_fits.csv`
- `experiments/2026-05-08_EXP-004_proteingym-predictive-joint-calibration-pairwise-capacity-panel/config.yaml`
- `experiments/2026-05-08_EXP-004_proteingym-predictive-joint-calibration-pairwise-capacity-panel/runs/RUN-012.yaml`
- `results/RES-003_run-010-weakens-hyp-001-despite-coordinate-aware-calibration/result.md`
- `results/RES-003_run-010-weakens-hyp-001-despite-coordinate-aware-calibration/metrics.json`

## Analysis performed

- Reviewed the completed Slurm run record and confirmed that the five declared durable outputs were recollected under `data/processed` and tracked with DVC.
- Rechecked the paired six-assay ProteinGym stability panel to ensure that `EXP-004` remained directly comparable to `RUN-010`.
- Evaluated the assay-specific `mavenn` fits on the raw assay scale to verify that the measurement layer remained well fit after the `EXP-004` code changes.
- Compared the four shared branches and the per-assay raw branch on single-mutant holdout, double-mutant holdout, epistasis-prediction, KS, and fitted-parameter summaries.
- Quantified the paired change relative to the `RUN-010` shared raw baseline to determine whether the predictive-joint objective and reduced pairwise compression rescued the empirical reconstruction failure.
- Assessed whether the improved fitting regime activated nonzero `epistasis_strength` or materially stronger pairwise structure.

## Outputs

- Figures: none
- Tables: `tables/run-012_key_metrics.md`

## Produced artifacts

- analyses/2026-05-09_ANA-004_run-012-predictive-joint-calibration-pairwise-capacity-review/tables/run-012_key_metrics.md

## Main observations

- `RUN-012` completed successfully and produced the five declared durable outputs: `summary.json`, `selected_panel.csv`, `mavenn_assay_metrics.csv`, `branch_validations.csv`, and `per_assay_branch_fits.csv`.
- The run remained paired to `RUN-010` on the same empirical target: `6` short `Tsuboyama 2023` `cDNA display proteolysis` stability assays spanning `4` taxa, `23,279` measured variants total, `6,424` single mutants, `16,855` multiple mutants, and sequence lengths from `44` to `72` residues.
- The measurement layer remained well fit. Across the six assays, the assay-specific `mavenn` models reached mean test Spearman `0.890` and mean test NRMSE `0.455`, so the `EXP-004` outcome is not explained by loss of raw-scale observation quality.
- No shared branch dominated the old baseline across all validation targets. `predictive_shared_raw` slightly improved single-mutant holdout Spearman to `0.202` and holdout NRMSE to `1.042`, but its double-mutant holdout Spearman worsened to `-0.197`. `predictive_richpair_shared_raw` gave the best shared-branch single-mutant holdout Spearman at `0.243`, the best shared-branch double-mutant holdout Spearman at `0.098`, and the best shared-branch functional KS at `0.463`, but its epistasis-prediction Spearman fell to `0.094`.
- Relative to the paired `RUN-010` shared raw baseline, the best `RUN-012` shared branch improved single-mutant holdout Spearman by `+0.011`, double-mutant holdout Spearman by `+0.115`, and functional KS by `-0.109`, but lost epistasis-prediction Spearman by `-0.407`.
- The best shared branch achieved those gains by moving the fitted reference close to the landscape peak rather than by recovering explicit epistatic structure. `predictive_richpair_shared_raw` fit `functional_reference_fraction_of_peak = 0.998`, `reference_distance_to_peak = 3`, `epistasis_strength = 0.0`, and `empirical_pairwise_strength = 0.0`.
- The shared latent branch remained weak overall. It fit `n_functional_dims = 6`, `peak_distance_from_consensus = 3`, `functional_reference_fraction_of_peak = 0.900`, `epistasis_strength = 0.0`, and `empirical_pairwise_strength = 0.0`, with single-mutant holdout Spearman only `0.053`.
- The per-assay raw fits used more pairwise capacity than the shared branches, but they still did not recover robust local geometry. Mean single-mutant holdout Spearman was `0.201`, mean double-mutant holdout Spearman was `0.089`, mean epistasis-prediction Spearman was `0.294`, all six fits kept `epistasis_strength = 0.0`, and `5/6` fits used nonzero `empirical_pairwise_strength`.
- Only one per-assay raw fit crossed `0.4` on single-mutant holdout rank recovery: `SPTN1_CHICK_Tsuboyama_2023_1TUD` reached single-mutant holdout Spearman `0.406`. The best double-mutant holdout rank recovery was `0.266` for `HECD1_HUMAN_Tsuboyama_2023_3DKM`.
- The preregistered `EXP-004` rescue prediction was therefore not met. The new objective and pairwise-capacity changes improved some coarse validation metrics on the shared raw branch, but they did not produce convincing recovery of variant-level single-mutant or double-mutant structure and did not activate nonzero `epistasis_strength`.
- On the current evidence, `RUN-012` further weakens `HYP-001`. It directly tests two calibration bottlenecks left by `RES-003`, yet the main reconstruction failure persists on the same well-controlled six-assay stability panel.

## Result records created

- `RES-004`

## Hypothesis updates

- `HYP-001` is further weakened, not refuted. The result reduces the plausibility that the remaining failure is explained mainly by the older staged objective or by overly aggressive empirical pairwise compression alone.
