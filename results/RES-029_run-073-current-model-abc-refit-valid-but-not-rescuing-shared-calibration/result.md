# RES-029: RUN-073 current-model ABC refit is valid but does not rescue shared empirical calibration

## Summary

`RUN-073` confirms that the PHOT and SPTN1 ABC refit used current
stability-margin and stability-routed epistasis definitions rather than legacy
latent-trait/roughness settings. It also shows that the matched current-model
synthetic inverse problem is recoverable. However, the empirical two-assay
posterior remains broad and plateau-like, and the joint fit does not recover
the stronger assay-specific epistasis-prediction performance seen in prior
PHOT/SPTN1 runs. No current-model graphlet, WHT, or GFT outputs were produced,
so this result resolves only the stale-cache concern for the ABC refit itself;
it does not overturn the `RUN-071` empirical graphlet falsification result.

## Generated from

- Analyses: `ANA-029`
- Experiments: `EXP-029`
- Runs: `RUN-073`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-001`
- Refutes: None

## Evidence

- `analyses/2026-05-16_ANA-029_run-073-current-model-two-assay-abc-refit-review/metrics.json`
- `analyses/2026-05-16_ANA-029_run-073-current-model-two-assay-abc-refit-review/tables/run-073_preflight.md`
- `analyses/2026-05-16_ANA-029_run-073-current-model-two-assay-abc-refit-review/tables/run-073_branch_validation.md`
- `analyses/2026-05-16_ANA-029_run-073-current-model-two-assay-abc-refit-review/tables/run-073_posterior_parameters.md`
- `analyses/2026-05-16_ANA-029_run-073-current-model-two-assay-abc-refit-review/tables/run-073_synthetic_truth_recovery.md`
- `analyses/2026-05-16_ANA-029_run-073-current-model-two-assay-abc-refit-review/tables/run-073_prior_fit_comparison.md`
- `results/RES-029_run-073-current-model-abc-refit-valid-but-not-rescuing-shared-calibration/metrics.json`
- `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073.dvc`
- `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073/summary.json`
- `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073/current_model_preflight.json`
- `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073/branch_validations.csv`
- `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073/posterior_parameter_summary.csv`
- `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073/synthetic_truth_recovery.csv`

## Interpretation

- `RUN-073` completed successfully with scheduler exit code `0`, and its
  durable outputs are DVC-tracked under
  `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073`.
- The current-model preflight passed. The observed readout was
  `stability_margin`, generic epistasis was routed to latent stability, and the
  fitted parameter set excluded the older FGM, designed-peak, deterministic
  roughness, noise, and higher-order exact-match settings.
- Matched synthetic truth recovery passed: both preregistered current-model
  synthetic truths placed all `20/20` fitted truth parameters inside posterior
  q90 intervals.
- The empirical SMC posterior did not concentrate strongly. The final best
  objective distance was `10.8228`, close to the round-0 best distance
  `10.7975`, with the overall round-best minimum only slightly lower at
  `10.7028`.
- The best empirical branch achieved functional KS `0.1862`, single-mutant
  holdout Spearman `0.4767`, double-mutant holdout Spearman `0.7691`, and
  epistasis-prediction Spearman `0.1681`.
- Compared with prior assay-specific fits, the joint current-model two-assay
  fit remains weak on epistasis prediction: `RUN-024` SPTN1 reached `0.5402`,
  and `RUN-039` PHOT reached `0.6344` on the same reported metric.
- The PHOT and SPTN1 target summaries differ substantially: PHOT has observed
  epistasis variance `0.0361` and beneficial fraction `0.0094`, while SPTN1
  has epistasis variance `0.6427` and beneficial fraction `0.2303`. This is a
  plausible explanation for the broad shared posterior, but it is not evidence
  that the shared model is validated.
- No dense current-model graphlet, Walsh, GFT, local maxima, accessible-path,
  or random-walk autocorrelation outputs are present. Therefore `RUN-073` does
  not directly retest the central `RES-028` graphlet falsifiers.

## Effect on hypothesis

`RES-029` weakens `HYP-001` in its historical shared-regime form. The result is
not a total failure of the pipeline: current semantics are enforced, and
synthetic truth recovery works. The failure is the empirical inference. A
single shared PHOT+SPTN1 current-model posterior is not a validated empirical
calibration and does not restore the stronger benchmark-realism claim that
`HYP-001` originally aimed at.

This result also keeps the `RES-028` graphlet conclusion active. The stale-cache
concern remains real for future parameterized graphlet comparisons, but this
run does not provide the missing current-model graphlet evidence needed to
reinterpret `RES-028`.

## Limitations

- The run used only `PHOT_CHLRE_Chen_2023` and
  `SPTN1_CHICK_Tsuboyama_2023_1TUD`, not the broader current-refit panel
  originally described in `EXP-029`.
- The result is based on ABC summary-vector calibration, not held-out graphlet
  calibration and testing.
- The PHOT MAVE-NN layer is only moderate (`test Spearman = 0.6823`), so PHOT
  target features may remain noisier than SPTN1 target features.
- The result is not evidence against all assay-specific or hierarchical current
  model variants. It is evidence against treating this restricted shared
  two-assay posterior as a rescued empirical calibration.

## Downstream use

- Use `RUN-073` as the current-semantics PHOT+SPTN1 ABC refit record, not as a
  graphlet/WHT/GFT validation result.
- Do not use old cached PHOT/SPTN1 latent targets as current-model empirical
  parameter sources unless explicitly labelled as stale-cache controls.
- Keep `RES-028` as the operative graphlet falsification result until a
  current-model refit emits directly comparable graphlet and spectral outputs.
