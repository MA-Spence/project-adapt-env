# ANA-029: RUN-073 current-model two-assay ABC refit review

## Purpose

Assess whether `RUN-073` resolves the stale empirical-cache concern for the
current smooth-core plus fixed-scale residual model by refitting the relevant
ProteinGym assays from raw data under current ABC semantics.

## Linked experiments/runs

- Experiments: `EXP-029`
- Runs: `RUN-073`
- Source data: `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073`

## Produced artifacts

- `metrics.json`
- `tables/run-073_preflight.md`
- `tables/run-073_panel.md`
- `tables/run-073_branch_validation.md`
- `tables/run-073_smc_rounds.md`
- `tables/run-073_posterior_parameters.md`
- `tables/run-073_synthetic_truth_recovery.md`
- `tables/run-073_target_features.md`
- `tables/run-073_prior_fit_comparison.md`

## Analysis scope

`RUN-073` is a restricted PHOT plus SPTN1 current-model ABC refit. It does not
execute the full broader `EXP-029` plan of refitting KRAS and six Tsuboyama
controls, nor does it emit posterior predictive dense graphlet, Walsh, or GFT
outputs. The analysis therefore addresses whether the ABC refit itself used
current semantics and how the resulting empirical fit compares with prior
assay-specific fits. It cannot replace the `RUN-071` graphlet falsification
suite.

## Checks performed

- Verified the run completed with scheduler exit code `0`.
- Verified the current-model preflight passed and the fitted parameters were
  limited to the current stability-margin/stability-routed epistasis set:
  `stability_scale`, `stability_margin`, `threshold_steepness`,
  `blosum_blend`, `stability_conservation_power`, `n_epistatic_terms`,
  `epistasis_strength`, `empirical_pairwise_strength`,
  `pairwise_epistasis_rank`, and `epistasis_covariation_bias`.
- Verified the observed readout was `stability_margin`, with generic epistasis
  targeted to latent stability.
- Summarised the empirical MAVE-NN assay layer, target features, SMC rounds,
  posterior parameter spread, synthetic truth recovery, and branch-validation
  metrics.
- Compared the best `RUN-073` current two-assay branch against prior
  assay-specific PHOT and SPTN1 fits.

## Key observations

- The procedural guard did what it was supposed to do. `RUN-073` did not fit
  the legacy FGM/peak-distance/noise/roughness/higher-order exact-match
  parameter set that motivated the stale-cache concern.
- Synthetic truth recovery is not the limiting failure mode here. Both matched
  current-model synthetic truths recovered all `20/20` truth parameters inside
  posterior q90 intervals.
- The empirical posterior is broad and weakly concentrated. The best SMC
  distance was already `10.7975` in round 0, the overall minimum round best was
  `10.7028`, and the final best branch was `10.8228`. The posterior q05-q95
  intervals span much of the allowed prior range for the residual parameters:
  `epistasis_strength` `0.0039-0.0912`,
  `empirical_pairwise_strength` `0.0038-0.0869`, and
  `n_epistatic_terms` `31.95-392.7`.
- The joint empirical fit is mixed. The best branch has functional KS `0.1862`,
  single-mutant holdout Spearman `0.4767`, double-mutant holdout Spearman
  `0.7691`, and epistasis-prediction Spearman `0.1681`.
- The PHOT and SPTN1 targets are structurally different under the current
  summary vector. PHOT has low observed epistasis variance (`0.0361`) and very
  low beneficial fraction (`0.0094`), while SPTN1 has much larger epistasis
  variance (`0.6427`) and beneficial fraction (`0.2303`). This makes a single
  shared current-model posterior hard to interpret as a validated empirical
  regime.
- Against prior assay-specific runs, `RUN-073` does not rescue the empirical
  calibration claim. It is competitive on double-mutant ranking but weaker than
  `RUN-024` on SPTN1 epistasis prediction (`0.1681` versus `0.5402`) and weaker
  than `RUN-039` on PHOT epistasis prediction (`0.1681` versus `0.6344`).
- No dense graphlet/WHT/GFT output is present in the collected run directory.
  Therefore the run resolves only the stale-cache concern for this PHOT+SPTN1
  ABC refit, not the full empirical graphlet target-envelope question raised by
  `RES-028`.

## Interpretation

`RUN-073` supports a narrow procedural conclusion: the current ABC pathway can
be run on raw PHOT and SPTN1 ProteinGym assays without silently using legacy
epistasis definitions, and the inverse problem is recoverable on matched
current-model synthetic data. That is useful because it removes one procedural
objection to the empirical calibration pathway.

The empirical scientific conclusion is more negative. The current shared
two-assay posterior is broad, plateau-like, and does not match the best prior
assay-specific empirical fits on the critical epistasis-prediction metric. It
also provides no new graphlet spectral or topological evidence. This means the
stale-cache explanation is not enough to treat the current model as validated:
`RES-028` remains the relevant evidence for graphlet ruggedness until a
current-model refit actually emits comparable held-out graphlet, WHT, GFT,
accessibility, autocorrelation, SE, and RSE statistics.

## Limitations

- The run used only `PHOT_CHLRE_Chen_2023` and
  `SPTN1_CHICK_Tsuboyama_2023_1TUD`, not the broader panel described in the
  original `EXP-029` design.
- The PHOT MAVE-NN layer is only moderately predictive
  (`test Spearman = 0.6823`, `test NRMSE = 0.9249`), so PHOT-derived target
  summaries remain noisy compared with the SPTN1 layer
  (`test Spearman = 0.9066`, `test NRMSE = 0.5504`).
- The summary-vector distance combines two quite different assays; the result
  should not be interpreted as an assay-specific optimum for either PHOT or
  SPTN1.
- No current-model graphlet, Walsh, GFT, local maxima, basin, accessible-path,
  or random-walk autocorrelation tables were produced.
- The large double-holdout NRMSE values in `branch_validations.csv` are not
  interpreted here because the Spearman and KS metrics are the stable
  comparison quantities used in the relevant prior results.
