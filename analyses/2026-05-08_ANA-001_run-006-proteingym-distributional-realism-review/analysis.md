# ANA-001: RUN-006 ProteinGym Distributional Realism Review

## Purpose

Review RUN-006 against HYP-001 by comparing the fitted shared Adapt-Env regime to ProteinGym single-mutant assay summary statistics and recording whether the observed evidence supports, weakens, or refutes the hypothesis.

## Linked experiments/runs

- Experiments: EXP-001
- Runs: RUN-006

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Inputs

- `data/processed/proteingym-dms-distributional-realism-panel/RUN-006/summary.json`
- `data/processed/proteingym-dms-distributional-realism-panel/RUN-006/selected_panel.csv`
- `experiments/2026-05-07_EXP-001_proteingym-dms-distributional-realism-panel/config.yaml`
- `experiments/2026-05-07_EXP-001_proteingym-dms-distributional-realism-panel/runs/RUN-006.yaml`
- `docs/hypotheses.md` entry for `HYP-001`

## Analysis performed

- Reviewed the completed Slurm run record to confirm execution status and declared outputs.
- Inspected the selected ProteinGym panel for assay breadth, taxonomic balance, and mutant counts.
- Extracted the fitted Adapt-Env parameters and the top-level calibration validation metrics from `summary.json`.
- Compared empirical and synthetic per-assay summary statistics, with emphasis on distribution tails, skewness, and conservation-sensitivity correlations.
- Assessed whether the resulting evidence supports the specific wording of `HYP-001`, which requires one shared regime to recover local summary statistics across assays without assay-specific fitting.

## Outputs

- Figures: none
- Tables: `tables/run-006_key_metrics.md`

## Produced artifacts

- analyses/2026-05-08_ANA-001_run-006-proteingym-distributional-realism-review/tables/run-006_key_metrics.md

## Main observations

- `RUN-006` completed successfully and produced the two declared durable outputs: `summary.json` and `selected_panel.csv`.
- The selected panel spans 8 single-mutant assays across 4 taxa (`Human`, `Prokaryote`, `Eukaryote`, `Virus`), 5 assay types, 27,249 measured single mutants, and sequence lengths from 60 to 364 residues.
- The fitted shared regime retained `n_functional_dims = 4` and `functional_sigma_base = 12.0`, but collapsed `epistasis_strength`, `empirical_pairwise_strength`, `noise_amplitude`, and `peak_distance_from_consensus` to `0.0` or `0`.
- Calibration quality was weak at the validation level: train NRMSE `0.992`, holdout NRMSE `0.998`, holdout Spearman `0.103`, and functional KS distance `0.537`.
- Across assays, the synthetic landscapes only partially recovered local score statistics. Mean absolute differences were `0.098` for the upper tail fraction above `+1 SD`, `0.045` for the lower tail fraction below `-1 SD`, `0.414` for the 5th percentile, `0.551` for the 95th percentile, `2.427` for score skewness, and `0.229` for conservation-sensitivity correlation.
- The beneficial tail was systematically compressed. The synthetic fraction above `+1 SD` was exactly `0.0` in 5 of 8 assays and remained below the empirical value in all 8 assays.
- Conservation-sensitivity direction was mostly preserved but not uniformly. The sign matched the empirical correlation in 7 of 8 assays, but `OXDA_RHOTO_Vanella_2023_activity` flipped from an empirical negative correlation (`-0.119`) to a synthetic positive correlation (`0.432`).
- Assay-level fit quality was heterogeneous. `PKN1_HUMAN_Tsuboyama_2023_1URF` showed the smallest aggregate mismatch, whereas `OXDA_RHOTO_Vanella_2023_activity` showed the largest mismatch, driven by large skewness, upper-tail, and conservation-sensitivity errors.
- On the current evidence, `RUN-006` weakens `HYP-001`. The run shows that some summary statistics can be approximated, but not that one shared parameter regime robustly recovers panel-level local landscape statistics without collapsing the richer mechanistic terms that the hypothesis invokes.

## Result records created

- `RES-001`

## Hypothesis updates

- `HYP-001` is weakened, not refuted. The current experiment tests the intended external realism question for single-mutant score distributions, but the recovered regime is too compressed and too weakly predictive to count as support.
