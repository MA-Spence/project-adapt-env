# EXP-029: ProteinGym empirical refit with current residual-ruggedness model

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001
- Directly linked results: RES-005, RES-006, RES-007, RES-025, RES-026,
  RES-027, RES-028
- Note: HYP-001 is superseded by CLAIM-001. This experiment does not reopen the
  original strong shared-regime claim as supported; it tests whether the recent
  empirical falsification suite used stale empirical targets that are no longer
  valid for the current model semantics.

## Question

Do current-model ABC refits of the relevant ProteinGym assays remove the stale-cache confound exposed by RUN-071 and yield empirical graphlet and spectral targets that can be fairly compared to the new residual-ruggedness model?

## Pre-experiment prediction

If RUN-071's mismatch is partly a stale empirical-target artifact, then refitting PHOT, SPTN1, KRAS and matched stability/activity controls under the current smooth-core plus fixed-scale residual model should change empirical latent targets and graphlet spectra enough that model-vs-empirical distances improve relative to RES-028 while retaining held-out predictive checks. If mismatch persists after refit, it is evidence against the current model rather than against old cached preprocessing.

## Rationale

RES-028 deliberately used cached prepared ProteinGym outputs from older calibration and latent-trait definitions, which makes it a valid falsification of comparisons to those cached targets but not a definitive test of the current model after the residual-scale and epistasis semantics were revised. This experiment rebuilds the empirical targets and ABC posteriors under the current model before using any parameterized empirical information.

## Experimental design

This is a current-model empirical refit experiment. It is not allowed to use
old cached `latent_phi`, `mavenn_yhat`, or fitted empirical parameter values as
targets for the current model, except in the explicit stale-cache replay
control. The primary inference is made only after raw ProteinGym assay data are
refit under the current smooth-core plus fixed-scale residual model with ABC.

The experiment has four phases:

1. Synthetic truth recovery under current semantics. Simulate assay panels from
   the current model, including the smooth thermodynamic core, nonlinear
   readout, and fixed-scale residual terms. Run the full ABC pipeline on these
   synthetic panels and require recovery of the fitted summaries and posterior
   predictive graphlet statistics before interpreting empirical fits.
2. Empirical current-model ABC refit. Refit the relevant ProteinGym systems
   from raw substitutions data. Primary systems are `PHOT_CHLRE_Chen_2023`,
   `SPTN1_CHICK_Tsuboyama_2023_1TUD`, and
   `RASK_HUMAN_Weng_2022_binding-DARPin_K55`, because these are the empirical
   graphlet systems used in RUN-071. The six-assay Tsuboyama stability panel
   from EXP-005 is retained as a context/control panel for HYP-001.
3. Posterior predictive graphlet extraction. For each fitted empirical system,
   sample posterior predictive landscapes and extract dense local Hamming
   graphlets. Record the same statistics in the ABC calibration set and in
   held-out test graphlets: SE, RSE, local maxima, accessible paths, random-walk
   autocorrelation, Walsh variance by order, GFT variance by frequency band,
   single/double DFE summaries, and global-readout-corrected residual spectra.
4. Direct comparison against stale cached targets and structural controls. Replay
   the exact cached graphlet targets used by RES-028 as a negative control, then
   compare current-model posterior predictive graphlets against smooth-core-only,
   residual-only, shuffled-contact residual, random-position residual, and GNK
   local-null controls.

The critical comparison is not whether the current model can exactly reproduce
one empirical assay. The test is whether refitting raw empirical systems under
the current model changes the empirical target envelope enough to remove the
stale-cache confound while preserving held-out predictive checks. If the current
model still cannot match spectra, SE/RSE, topology, and accessibility after this
refit, the RES-028 structural falsification remains valid.

## Inputs

- Raw ProteinGym v1.3 substitutions metadata and substitutions parquet.
- ProteinGym systems requiring current-model ABC refit:
  `PHOT_CHLRE_Chen_2023`, `SPTN1_CHICK_Tsuboyama_2023_1TUD`,
  `RASK_HUMAN_Weng_2022_binding-DARPin_K55`,
  `HECD1_HUMAN_Tsuboyama_2023_3DKM`,
  `SDA_BACSU_Tsuboyama_2023_1PV0`,
  `MYO3_YEAST_Tsuboyama_2023_2BTT`,
  `POLG_PESV_Tsuboyama_2023_2MXD`, and
  `CSN4_MOUSE_Tsuboyama_2023_1UFM`.
- Prior cached empirical outputs from RUN-024, RUN-026, RUN-028, and RUN-071
  only for stale-cache replay and comparison. They are not valid current-model
  parameter sources.
- Current Adapt-Env model implementation and project-local ABC calibration
  utilities.

## Configuration

- `config.yaml` records the assay panel, current-model parameter priors,
  ABC summaries, synthetic-truth checks, graphlet metrics, stale-cache replay
  controls, and falsification criteria.
- The planned runner is intentionally separate from the older
  `proteingym_bayesian_summary_calibration.py` semantics unless that script is
  updated to emit current-model residual and graphlet summaries.

## Execution

- Run `scripts/proteingym_current_model_abc_refit.py`, which first checks that
  the EXP-029 config uses only the PHOT and SPTN1 assays and does not fit
  legacy FGM, designed peak, deterministic roughness, or higher-order
  exact-match parameters.
- After the guard passes, the wrapper delegates to the same
  `scripts/proteingym_bayesian_summary_calibration.py` SMC-ABC pathway used in
  prior ProteinGym calibration experiments.
- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- current_model_preflight.json
- synthetic_truth_recovery.csv
- empirical_refit_summary.json
- assay_refit_metrics.csv
- selected_panel.csv
- mavenn_assay_metrics.csv
- target_features.csv
- posterior_particles.csv
- posterior_rounds.csv
- posterior_parameter_summary.csv
- branch_validations.csv
- progress.json
- progress.log

## Analysis plan

- Create an ANA record after execution.
- First verify synthetic truth recovery. If the ABC pipeline cannot recover
  current-model synthetic data, do not interpret empirical posterior mismatch as
  biological evidence.
- Compare current-model empirical refits to stale cached graphlet targets from
  RES-028. A substantial shift in empirical spectral and SE/RSE targets would
  support the stale-cache-confound explanation; no shift would preserve the
  RES-028 falsification.
- Evaluate held-out graphlets separately from calibration graphlets. A fit that
  only improves calibration graphlets is overfit and does not support HYP-001 or
  CLAIM-001.
- Compare current-model posterior predictive graphlets against smooth-core-only,
  residual-only, shuffled-contact, random-position, and GNK controls using the
  same joint distance family as RES-028.
- Reassess HYP-001 only in its historical context. Even if this experiment
  improves empirical alignment, it can at most support the reduced
  assay-specific interpretation captured in CLAIM-001 unless a shared or
  hierarchical parameter structure generalizes across assays.

Falsification criteria:

- Synthetic truth recovery fails under current semantics.
- Refit empirical graphlets remain closest to smooth-core-only, residual-only,
  GNK, random-position, or shuffled-contact controls rather than the embedded
  current model.
- Current-model refits do not materially change the empirical target envelope
  relative to the stale cached targets used by RES-028.
- Refits improve calibration graphlets but fail held-out graphlets.
- High-order WHT/GFT signal disappears after fitting an appropriate nonlinear
  global readout, implying the residual ruggedness signal was scale
  confounded.
- Assay-specific tuning is still required with no stable shared or hierarchical
  structure, preserving the original HYP-001 falsification even if individual
  systems can be fit.

## Completion criteria

- [ ] Config committed
- [ ] Run script committed
- [ ] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
