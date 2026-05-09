# Project state

- Status: active aim, hypotheses recorded
- Active aim: `AIM-001 Create-Model-Landscape`

## Current scientific state

- The parent project still has one active aim and had no project-local hypotheses before this review.
- The practical scientific state currently lives in `external/Adapt-Env`, which already contains a synthetic landscape implementation, a population-genetics layer, and a validation suite aimed at qualitative protein-like behavior.
- The current model claim that is best supported is narrow: AdaptEnv is a biologically informed synthetic environment for benchmarking iterative search methods, not an exact reconstruction of a real protein fitness landscape.

## What the current quick validations support

- `experiment_dfe_characterisation --quick` supports protein-like local DFE structure.
  - unconditioned consensus neutral fraction: `0.677`
  - conditioned consensus neutral fraction: `0.906`
  - beneficial fraction increases with distance from peak: supported
  - DFE left skew: supported
- `experiment_alignment_conditioning --quick` supports a real effect of conditioning, though not yet external biological validity.
  - conservation-sensitivity Pearson: conditioned `0.385`, unconditioned `-0.067`
  - strong-covariation pairs show higher epistasis than random pairs: `0.593` versus `0.241`
  - sampled natural-like sequences remain viable: `1.0`
- `experiment_global_epistasis --quick` supports negative global epistasis across tested FGM dimensionalities.
  - mean Pearson correlations by dimensions `2/4/8`: `-0.557`, `-0.579`, `-0.536`
  - the intended monotonic strengthening with dimensionality is not supported
- `experiment_wf_dynamics --quick` supports qualitative adaptation and diversity turnover but not a strong long-term dynamics match.
  - early diversity increase: supported
  - baseline deceleration target: failed
  - truncation selection outperforms proportional selection in the tested setup
- `experiment_synthetic_config_realism --quick` supports most qualitative preset checks but not all regime claims.
  - alignment sequences viable across presets: supported
  - negative global epistasis majority: supported
  - smooth-bottom / rugged-top majority: supported
  - high-throughput beating low-throughput: failed in the quick run
  - viral earlier improvement than bacterial: failed in the quick run
- `RUN-006` for `EXP-001` now provides the first external ProteinGym panel test under `HYP-001`.
  - 8 assays selected across 4 taxa with 27,249 single mutants total
  - real MMseqs alignments were generated and saved under `data/interim/proteingym_mmseqs_alignments`
  - the best shared fitted regime set `epistasis_strength`, `empirical_pairwise_strength`, and `noise_amplitude` to zero
  - holdout performance was weak (`Spearman 0.103`, `NRMSE 0.998`)
  - the run therefore weakens `HYP-001` rather than supporting it

## Main scientific weaknesses at the current state

- Peak interpretation is still the most important weakness for benchmark use.
  - `experiment_peak_sensitivity_noise_regimes --quick`: noise-free peak retention failed for `0/4` presets, implying the issue is mostly topology rather than noise.
  - `experiment_peak_sensitivity_pairwise_modes --quick`: an internal peak-preserving pairwise mode improved no-noise peak distance in `4/4` presets, but that mode is not the public default.
- Many validation experiments depend on the stored `peak_sequence`, which is local-search derived rather than globally certified.
- Most experiment helpers validate alignment conditioning using internally generated synthetic alignments in `external/Adapt-Env/experiments/common.py`, so those tests establish internal consistency more strongly than external biological realism.
- Several key knobs are not yet behaving as fully interpretable scientific controls.
  - FGM dimensionality does not cleanly strengthen diminishing returns
  - Wright-Fisher deceleration is weaker than intended
  - some preset-level ranking claims are configuration-sensitive
- The current validation suite is largely qualitative and internal.
  - it lacks preregistered external targets drawn from empirical benchmark panels
  - it does not yet cleanly separate "biological realism" from "benchmark usefulness"

## Record update from this review

- Project-local hypotheses have now been recorded in `docs/hypotheses.md` and the `labproj` registry.
- `EXP-001` has been created under `HYP-001` to assemble a real ProteinGym DMS panel, fetch MMseqs alignments from assay wild-type sequences, and fit a shared synthetic regime against external distributional statistics.
- `ANA-001` records the scientific review of `RUN-006`, and `RES-001` records that the present evidence weakens `HYP-001`.
- `EXP-002` has been created under `HYP-001` as an indirect scope check on the uncalibrated model family.
  - it scans `stability_margin`, `functional_sigma_base`, `n_functional_dims`, and `epistasis_strength`
  - it measures one-step and two-step DFEs around the reference sequence without DMS calibration
  - it records double-mutant deviations from additive expectation
  - `RUN-008` completed successfully on `lab-slurm` under scheduler job `31`
  - the uncalibrated family spans a broad one-step and two-step DFE envelope, with doubles always more deleterious than singles across the scanned settings
  - `stability_margin` dominates local harshness, while `epistasis_strength` mainly changes double-mutant epistasis magnitude
  - this keeps `HYP-001` scientifically live as a realism target, but does not support it directly because no empirical comparison is involved
- `ANA-002` records the scientific review of `RUN-008`, and `RES-002` records that the result motivates `HYP-001` without resolving it.
- `EXP-003` has been created under `HYP-001` to resolve the main interpretability flaw left by `RES-001`.
  - it restricts the panel to homogeneous multi-mutant ProteinGym stability assays instead of mixing assay families
  - it adds `mavenn` to the run environment and fits assay-specific latent observation maps from raw assay scores
  - it replaces assay-wise z-scoring and synthetic wild-type zeroing with assay-specific raw-scale and latent wild-type anchors inferred from `mavenn`
  - it compares a shared Adapt-Env fit on raw scores with non-affine observation models against a shared latent-scale fit and per-assay latent fits
  - `RUN-010` completed successfully on `lab-slurm` under scheduler job `32`
  - the promoted durable outputs are now stored under `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010`
  - the selected panel contained `6` homogeneous Tsuboyama `cDNA display proteolysis` stability assays spanning `23,279` measured variants, including `16,855` multiple mutants
  - assay-specific `mavenn` models fit the raw assay scale well, with mean test Spearman `0.894`
  - the shared raw-score Adapt-Env branch was modestly better than the shared latent branch, but both remained weak on holdout recovery and both shared fits retained `epistasis_strength = 0.0`
  - per-assay latent fits also remained poor, so removing z-scoring and synthetic wild-type zeroing did not rescue the empirical calibration claim
  - `ANA-003` records the scientific review of `RUN-010`, and `RES-003` records that the result further weakens `HYP-001`
- `EXP-004` has now been created under `HYP-001` to test the two remaining calibration bottlenecks left after `RES-003`.
  - it fixes the empirical target to the exact six-assay `EXP-003` Tsuboyama stability panel so branch comparisons are paired on the same observed landscapes
  - it retains the `EXP-003` shared raw branch as a baseline control
  - it adds new predictive-joint calibration branches that use a DFE-and-predictive core objective instead of the earlier summary-first staged fit
  - it adds new rich-pairwise branches that relax empirical pairwise compression through less aggressive prior blending, coverage shrinkage, and normalization
  - it compares shared raw, shared latent, and per-assay raw fits under the improved settings
  - the required external code changes were committed and pushed to `external/Adapt-Env` branch `feat/predictive-joint-calibration` at commit `8ce716211160f83c0d4c3eaa13a9e0ffefde1814`
  - `RUN-012` completed successfully on `lab-slurm` under scheduler job `33`
  - the promoted durable outputs are now stored under `data/processed/proteingym-predictive-joint-calibration-pairwise-capacity-panel/RUN-012` and tracked with DVC
  - the run remained paired to the `EXP-003` six-assay Tsuboyama stability panel with `23,279` measured variants, including `16,855` multiple mutants
  - assay-specific `mavenn` models again fit the raw assay scale well, with mean test Spearman `0.890`
  - the modified shared fitter improved some coarse shared raw metrics relative to the `RUN-010` shared raw baseline, with best shared-branch single-mutant holdout Spearman `0.243`, double-mutant holdout Spearman `0.098`, and functional KS `0.463`
  - those gains did not rescue variant-level reconstruction: the best shared branch still fit `epistasis_strength = 0.0` and `empirical_pairwise_strength = 0.0`, moved the reference very near the peak, and had weak epistasis-prediction Spearman `0.094`
  - per-assay raw fits used some nonzero empirical pairwise strength in `5/6` assays, but all six still fit `epistasis_strength = 0.0`, and mean single- and double-mutant holdout Spearman remained only `0.201` and `0.089`
  - `ANA-004` records the scientific review of `RUN-012`, and `RES-004` records that the result still weakens `HYP-001`
- The project-local `src/` package has been renamed from the placeholder `scientific_project` to `project_adapt_env`.
  - the new package now contains project-side experiment code instead of leaving all workflow logic in ad hoc scripts
  - `src/project_adapt_env/smc_abc.py` implements a reusable SMC-ABC engine
  - `src/project_adapt_env/adapt_env_bayes.py` implements Adapt-Env-specific summary targeting, bootstrap covariance estimation, and synthetic-truth recovery helpers
  - `src/project_adapt_env/proteingym_panel.py` factors out ProteinGym panel assembly, MMseqs alignment fetching, and `mavenn`-anchored empirical landscape construction from the earlier scripts
- `EXP-005` has now been created under `HYP-001` to test whether a likelihood-free Bayesian calibration can recover the six-assay ProteinGym stability panel more faithfully than the current deterministic fitters.
  - it keeps the exact fixed `EXP-003` / `EXP-004` six-assay Tsuboyama stability panel so the Bayesian calibration is paired against the existing evidence
  - it retains two deterministic controls: `baseline_shared_raw` and `predictive_richpair_shared_raw`
  - it builds a structured per-assay summary target from observed single-mutant and double-mutant statistics and estimates the target covariance by bootstrap
  - it fits a shared posterior over `stability_scale`, `stability_margin`, `blosum_blend`, `stability_conservation_power`, `functional_sigma_base`, `n_functional_dims`, `peak_distance_from_consensus`, `epistasis_strength`, `empirical_pairwise_strength`, and `noise_amplitude`
  - it includes preregistered synthetic-truth recovery on the same panel scaffold before the empirical posterior is interpreted as evidence about `HYP-001`
  - the runner is `scripts/proteingym_bayesian_summary_calibration.py`, and the experiment record lives under `experiments/2026-05-09_EXP-005_proteingym-shared-summary-smc-abc-calibration-panel`
  - local code validation passed through `python -m compileall` and `pytest tests/test_smc_abc.py`, while a local `EXP-005 --quick` workflow smoke was blocked only by the absence of `mavenn` in the current interpreter
  - `RUN-014` is now running on `lab-slurm` under scheduler job `55`
