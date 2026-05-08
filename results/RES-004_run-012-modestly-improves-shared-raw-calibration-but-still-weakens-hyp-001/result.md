# RES-004: RUN-012 modestly improves shared raw calibration but still weakens HYP-001

## Summary

Predictive-joint and richer pairwise calibration modestly improved the best shared raw single-mutant holdout ranking and double-mutant holdout ranking on the six-assay ProteinGym stability panel, but the fitter still failed to recover strong variant-level structure and still drove epistasis_strength to zero.

## Generated from

- Analyses: ANA-004

## Relevant hypotheses

- Supports: None
- Weakens: HYP-001
- Refutes: None

## Evidence

- `analyses/2026-05-09_ANA-004_run-012-predictive-joint-calibration-pairwise-capacity-review/tables/run-012_key_metrics.md`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/summary.json`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/selected_panel.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/mavenn_assay_metrics.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/branch_validations.csv`
- `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012/per_assay_branch_fits.csv`
- `experiments/2026-05-08_EXP-004_proteingym-predictive-joint-calibration-pairwise-capacity-panel/config.yaml`
- `experiments/2026-05-08_EXP-004_proteingym-predictive-joint-calibration-pairwise-capacity-panel/runs/RUN-012.yaml`
- `results/RES-003_run-010-weakens-hyp-001-despite-coordinate-aware-calibration/metrics.json`

## Interpretation

- `RUN-012` is the first paired test of whether the two remaining calibration hypotheses left by `RES-003` were correct: that a more predictive joint objective and less aggressive empirical pairwise compression would materially rescue the empirical reconstruction on the same six-assay ProteinGym stability panel.
- The measurement layer remained strong. The assay-specific `mavenn` fits on raw assay scores again performed well, with mean test Spearman `0.890` and mean test NRMSE `0.455`, so this run does not fail because the raw observation scale became poorly modeled.
- The modified fitter improved some shared raw metrics, but not in a clean or sufficient way. The best shared branch, `predictive_richpair_shared_raw`, reached single-mutant holdout Spearman `0.243`, double-mutant holdout Spearman `0.098`, and functional KS `0.463`, which is modestly better than the paired `RUN-010` shared raw baseline on those three metrics.
- That branch did not recover stronger explicit pairwise or epistatic structure. It fit `epistasis_strength = 0.0` and `empirical_pairwise_strength = 0.0`, drove the reference close to the peak with `functional_reference_fraction_of_peak = 0.998`, and its epistasis-prediction Spearman fell to `0.094`. This means the apparent improvement is not a clean rescue of local landscape reconstruction.
- The predictive-only shared raw branch and the shared latent branch also failed to provide a convincing rescue. The predictive-only branch improved some single-mutant metrics but worsened double-mutant holdout ranking, and the latent branch remained weak on both ranking and epistasis-prediction targets.
- The per-assay raw fits partially used nonzero `empirical_pairwise_strength` in `5/6` assays, but all six still kept `epistasis_strength = 0.0`, mean single-mutant holdout Spearman was only `0.201`, and mean double-mutant holdout Spearman was only `0.089`. This does not support the idea that the remaining failure was mainly due to pooling across assays.

## Effect on hypothesis

- This result further weakens `HYP-001`.
- `RES-003` left open the possibility that the calibration failure mainly reflected the old staged objective or over-compressed pairwise fitting. `RUN-012` directly targeted those two explanations on the same empirical panel, and the main failure persisted.
- The result does not refute `HYP-001` outright because the test still uses a six-assay stability-only panel and one implemented fitter family, but it materially lowers confidence that modest objective redesign and richer pairwise capacity are enough to recover the intended empirical statistics.

## Limitations

- The panel is better controlled than `RUN-006`, but it is still limited to `6` short `Tsuboyama 2023` stability assays from one experimental platform.
- The branch comparison is informative but not exhaustive. `EXP-004` tests one concrete predictive-joint and rich-pairwise implementation, not all possible objective designs or pairwise parameterizations.
- Some improvements in the best shared branch may be partly attributable to moving the fitted reference near the peak rather than to genuinely better recovery of mutational effects.
- One assay still produced null holdout Spearman values in the stored per-assay summary, which limits the precision of per-assay aggregate comparisons.

## Downstream use

- Use `RES-004` as the current paired evidence on whether the `EXP-004` calibration changes rescue `HYP-001`; they do not.
- Use this result together with `RES-002` and `RES-003` when interpreting `HYP-001`: the model family remains locally expressive, but the present empirical inverse problem remains inadequately solved even after targeted fitter changes.
