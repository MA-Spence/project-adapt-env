# RES-005: RUN-019 Bayesian SMC-ABC recovers nonzero epistasis but still weakens HYP-001

## Summary

RUN-019 completed the resumed Bayesian SMC-ABC calibration and no longer collapsed the shared fit to zero epistasis, while preregistered synthetic-truth recovery succeeded on the same panel scaffold. Empirical holdout recovery on the six-assay ProteinGym stability panel nevertheless remained weak relative to the scientific claim, so the result still weakens HYP-001.

## Generated from

- Analyses: ANA-005

## Relevant hypotheses

- Supports: None
- Weakens: HYP-001
- Refutes: None

## Evidence

- `analyses/2026-05-12_ANA-005_run-019-bayesian-shared-summary-smc-abc-calibration-review/tables/run-019_key_metrics.md`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/summary.json`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/selected_panel.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/mavenn_assay_metrics.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/branch_validations.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/posterior_particles.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/posterior_rounds.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/posterior_parameter_summary.csv`
- `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019/synthetic_truth_recovery.csv`
- `experiments/2026-05-09_EXP-005_proteingym-shared-summary-smc-abc-calibration-panel/config.yaml`
- `experiments/2026-05-09_EXP-005_proteingym-shared-summary-smc-abc-calibration-panel/runs/RUN-019.yaml`

## Interpretation

- `RUN-019` is the first completed end-to-end test of the Bayesian shared-summary calibration in `EXP-005`. It resumed the interrupted `RUN-017` empirical SMC state, finished the final empirical round, and completed both preregistered synthetic-truth recovery branches on the same six-assay ProteinGym stability panel.
- The measurement layer remained strong. Across the six assays, the assay-specific `mavenn` fits on raw assay scores reached mean test Spearman `0.890` and mean test NRMSE `0.480`, so the empirical outcome is not explained by poor raw-scale observation modeling.
- The Bayesian fit no longer failed by trivial collapse. The best empirical SMC particle used nonzero `epistasis_strength = 0.0527`, `empirical_pairwise_strength = 0.0008`, and `noise_amplitude = 0.00044`, while the posterior mean retained `epistasis_strength = 0.0429` and `empirical_pairwise_strength = 0.0050`. This is qualitatively different from the earlier shared deterministic fits, which repeatedly drove `epistasis_strength` to zero.
- The preregistered synthetic-truth checks were substantially more successful than the empirical reconstruction. On the matched panel scaffold, the `moderate_epistatic` truth lay within the posterior q90 interval for `10/10` parameters and the `flatter_low_epistasis` truth did so for `9/10`, missing only `noise_amplitude`. That means the inference procedure itself can recover structured regimes on model-matched targets.
- The empirical fit still did not rescue the scientific claim behind `HYP-001`. The best Bayesian shared fit reached only single-mutant holdout Spearman `0.195`, double-mutant holdout Spearman `0.119`, and functional KS `0.483`. Compared with the best deterministic control inside the same run, it improved epistasis-prediction Spearman (`0.265` versus `0.094`) and avoided the near-peak reference artifact (`functional_reference_fraction_of_peak 0.049` versus `0.998`), but it was worse on single-mutant holdout ranking (`0.195` versus `0.243`) and essentially tied on functional KS.
- The main scientific implication is therefore narrower than support for the hypothesis. `RUN-019` reduces the plausibility that the failure is caused mainly by optimizer collapse or inability of the fitter to represent nonzero epistasis. The remaining mismatch is more consistent with a limitation in the shared model family, the chosen summary target, or the assumption that one shared regime can describe this empirical panel well enough.

## Effect on hypothesis

- This result further weakens `HYP-001`.
- It does not refute `HYP-001` outright because the test is still limited to `6` short `Tsuboyama 2023` stability assays from one assay class, and the Bayesian objective is only one implemented shared-summary inference strategy.
- It is stronger negative evidence than the interrupted `RUN-017` record because the resumed run completed the full preregistered workflow, including the synthetic-truth checks, and still failed to provide convincing empirical recovery of the local landscape statistics.

## Limitations

- The empirical panel is scientifically cleaner than the earlier mixed-family runs, but it is still narrow: `6` short stability assays from one study and one readout platform.
- The synthetic-truth recovery is informative about the inference procedure only because the truths are drawn from the same model family and panel scaffold. Success there does not imply that the model family is adequate for the real ProteinGym assays.
- The analysis evaluates a shared summary-targeted posterior, not a full variant-level or assay-specific reconstruction pipeline, so weak holdout recovery could still reflect what the current summaries discard as well as what the shared model cannot express.
- The durable `RUN-019` outputs are present under `data/processed`, but they are not yet DVC-tracked in the current repository state.

## Downstream use

- Use `RES-005` as the current result record for `EXP-005` when discussing the Bayesian calibration path for `HYP-001`.
- Use this result to constrain claims about external realism: the project can no longer say only that the fitter collapsed to zero epistasis, but it still cannot claim that one shared synthetic regime convincingly recovers the six-assay ProteinGym stability panel.
- Use this result together with `RES-002`, `RES-003`, and `RES-004` when interpreting `HYP-001`: the model family appears locally expressive and the Bayesian fitter works on matched synthetic targets, yet the shared empirical reconstruction remains inadequate.
