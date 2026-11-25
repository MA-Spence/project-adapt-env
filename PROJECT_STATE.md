# Project state

## 2026-11-26 RES-064 — fixed-noise epistasis and public methods

- Created RES-064 linking EXP-065 RUN-186/188/190 and ANA-072/073/074. Fixed
  genotype-level observation errors add epistasis outside latent r control;
  they do not account for all synthetic epistasis. Current TrpB replicates and
  stop controls do not uniquely identify biological versus shared-error
  contributions. This is conditional nonidentifiability, not impossibility in
  every experimental design. The synthetic decomposition is exactly known.
- Gamma concentration gave no demonstrated advantage over Gaussian. Gaussian
  remains the default; no release, calibrated-preset promotion or new fit.
- AdaptEnv source commit `a6b5bfc5a86efa4f1a948d1b4ca968757b3fb343` is pushed
  to `origin/dev`; the project submodule and external-package lock pin it.
- Public AdaptEnv APIs now expose the direct aggregate objective, serial joint
  Optuna r/readout fitting, optional alpha fitting, fixed-r ablation, explicit
  replicate errors and frozen banks. The separate experimental auxiliary-
  inference method remains clearly distinguished. Immutable RUN snapshots
  were not changed.
- Package validation: 234 tests pass in Python 3.13; six new API tests and
  archived-method parity pass in the production Python 3.11 environment. The
  installed wheel passes the same archived parity check. Historical ranked
  fingerprint fixture differences across environments are documented in
  RES-064/package_verification.json; no historical fixtures were changed.
- Validation replays fixed independent/shared observations and banks exactly,
  reconstructs Gaussian/Gamma calibration fields, and reproduces all four loss
  components for six historical starts and selected RUN-186/188/190 settings
  (maximum difference 2.8e-17). See RES-064/public_api_parity.json.


## ANA-074 — Gamma concentration adds no demonstrated TrpB advantage (25 November 2026)

- Completed RUN-190 analysis with seven visually checked PNG/PDF figures and a ten-code-cell executed notebook. Gamma selects alpha=100, r=0.5 and the Gaussian readout: loss 0.03103213 versus 0.03099517 (0.12% higher). No later Optuna proposal beats the queued candidate.
- Selected latent pair-table correlation with Gaussian has fresh-seed median 0.99875; top 1% squared coefficient mass is 8.53% versus 8.40%. Confirmation changes are small and mixed: pooled high-order fractions remain about 0.20 versus empirical 0.086; sign prevalence is close, while activity, absolute powers and seed variation remain limitations.
- Static difficulty proxies remain similar. No learner/search task, new independent empirical validation, model promotion or new result/claim was performed. This does not prove equivalence or a global optimum.
- Verified 65 output and 536 frozen files, 11,712 completed trials, exact Gaussian/r0 replay, independent powers/quantiles/slopes/sign checks, and reconstruction of 48 full latent pair fields. RUN-190 is locally collected, DVC tracked and backed up; all 543 unique cache objects verify locally and remotely.
- [Report](analyses/2026-11-25_ANA-074_run-190-gamma-concentration-versus-gaussian-trpb-calibration/analysis.md) and [executed notebook](analyses/2026-11-25_ANA-074_run-190-gamma-concentration-versus-gaussian-trpb-calibration/notebooks/analysis.ipynb). ANA-074 archive `93cbf1674a6c316f03db00a735123a9b.dir` is backed up with all cache objects verified locally and remotely; its 42 canonical payload files match the archive. Verification is recorded with the analysis.

## RUN-190 TrpB Gamma concentration calibration completed (25 November 2026)

RUN-190 (EXP-065; Slurm 2357 on Dutch) completed successfully in 8,809.31 seconds (~2 h 27 min). All 11,712 trials completed with no failures and all 24 confirmation seeds finished; both exit codes are zero, stderr is empty, and all 65 manifest-listed output files passed SHA-256 verification. Reconciled through labproj. Completion evidence: `experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-190.completion-verification.json`.

The selected Gamma fit retains the queued alpha=100, r=0.5 candidate (trial 7): loss 0.03103213 versus the reused exact Gaussian control 0.03099517 (~0.12% higher). No new Optuna candidate improved this warm start. This achieved fixed-budget calibration provides no gain from stronger concentration of interaction magnitudes; it is not a global optimality or equivalence proof. Confirmation diagnostics analysed in ANA-074 above; no model promotion.

Design: joint alpha [0.05,100] (log), r [0,8] and four readout coordinates; unchanged direct loss, fixed shared-stop noise, panels and seeds; RUN-188 Gaussian/r0 controls reused. RUN-189 was the nonexecuting preview. Local 27-test suite and Linux smoke passed; frozen runtime DVC directory `705126f8fe5da52a44a94e6459a918d3.dir` has 470 unique cache objects verified locally and in backup. Production outputs have now been collected and DVC archived with verified project backup as part of ANA-074.

## 2026-11-25 ANA-073 — RUN-188 observation noise versus structured interactions

- Created the requested figure-led ANA with six PNG/PDF figures and an executed
  notebook (11 code cells), paired source, exact noise-component reconstruction
  and independent cube-power/face-classification checks. No refitting or new draws.
- At primary resolution, r0 intrinsic sign epistasis is zero but observed median
  is 62.5%; full r=0.5 moves from 9.9% intrinsic to 64.1% observed, versus TrpB
  65.1%. Median 90.4% of full-model observed sign-epistatic faces were absent
  intrinsically (72.3% at zero resolution; 98.3% at coarser resolution).
- Noise dominance depends on the metric: intrinsic after-sigmoid fitness supplies
  91.1% of pooled full-model order2 power, 61.3% of order3, and 16.3% of order4.
  Shared error supplies 72.3% of order4 power; signed cross terms are retained.
  These are synthetic power accounts, not identified biological noise fractions.
- Full-model calibration loss is 26.3% lower, but unfitted diagnostic advantages
  are mixed. Shared noise produces repeatable reversals; sparse all-active strata
  and negative near-floor stop-adjusted powers leave biological interpretation
  unresolved. The observed map is not established as principally r-controlled.
- All 41 output payloads and 460 frozen payloads verify; all 48 observation replays
  match exactly. RUN-188 is now locally collected, DVC tracked and backed up;
  all 460 unique local/remote storage objects verify. No model promotion, new
  result/claim, or ML/search experiment was performed.
- Completed ANA-073 archive is DVC tracked and backed up; all 41 unique
  local/remote objects verify and canonical files match the archive byte-for-byte.
- [ANA-073 report](analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/analysis.md)
  and [notebook](analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/notebooks/analysis.ipynb).

## 2026-11-25 RUN-188 completed successfully

- Both calibration studies completed all 11,712 trials (23,424 total) with no
  failures. All 24 confirmation seeds completed. Scientific runtime was
  12,758.31 seconds (3 h 32 min 38 sec); scientific and runner exit codes are zero,
  stderr is empty, and all 41 manifest-listed output files passed SHA-256 checks.
- The r=0 control retains loss 0.04207450018; the full model retains r=0.5 and
  loss 0.03099516773 (26.3% lower). Both winners are historical queued candidates;
  no new Optuna proposal improved them. This is a calibration result, not proof
  that the observation model is justified or that latent interactions account for
  the empirical ruggedness.
- Replicate-aware diagnostics are complete (800 cube rows, 57,600 face rows,
  1,164 summary rows and 291 graph rows), but have not received a new scientific
  analysis. No ANA, result, claim or package change is made in this status check.
- Scheduler retention no longer includes job 2356; completion is supported by
  both saved zero exit codes, finished studies and authenticated final outputs.
  Evidence: [completion snapshot](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-188.completion-status.json).
  Production outputs remain on Dutch; full collection/archival is pending.

## 2026-11-25 RUN-188 progress at 13:25 AEDT

- Slurm 2356 remains RUNNING on Dutch after 1 h 19 min. The r=0 control
  has completed 9,512/11,712 trials (81.2% of that arm; 40.6% of the combined
  fitting budget), with no failed trials and empty stderr.
- Best control loss remains 0.04207450018, the initial queued candidate.
  Full-r fitting and post-fit diagnostics are still pending. No scientific
  comparison can yet be concluded from this fitting progress.
- Evidence: [progress snapshot](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-188.progress-20261125T022546Z.json).

## 2026-11-25 RUN-188 progress at 12:42 AEDT

- Slurm 2356 remains RUNNING on Dutch after 35 min 53 sec. The r=0 control
  has completed 4,580/11,712 trials (39.1% of that arm; 19.6% of the combined
  fitting budget), with no failed trials and empty stderr.
- Best r=0 calibration loss remains 0.04207450018, the initial queued candidate.
  The full-r arm and post-fit replicate-aware diagnostics have not started.
  Trial progress is not total pipeline completion or evidence of adequacy.
- Read-only scheduler/log/SQLite evidence is retained in
  [progress snapshot](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-188.progress-20261125T014222Z.json).

## 2026-11-25 RUN-188 — matched interaction ablation running

- Registered the authorized TrpB comparison under EXP-065 as RUN-188 and submitted
  Slurm job 2356. Verified RUNNING on Dutch with one CPU, 8 GB RAM and a 24-hour
  limit. RUN-187 is the nonexecuting submission preview.
- Independently refit the profile-additive sigmoid control at r=0 and the full
  pairwise model with continuous r in [0,8], each with 11,712 Optuna trials. The
  RUN-186 shared-stop law, calibration loss, panels and generator/noise seeds are
  frozen. Noise strength and Gamma concentration are not fitted.
- Unfitted diagnostics measure absolute cross-replicate interaction powers,
  matching sign reversals, activity strata and intrinsic/observed graph geometry.
  Previously inspected confirmation panels are reused; shared assay error remains
  unidentified. No learning/search tasks or package promotion is performed.
- All 24 targeted tests, local and Linux end-to-end smoke runs passed. Worker
  verified 460 frozen payloads and 36 exact package pins; all six historical
  objectives replayed, the r0 fit has started and initial stderr is empty.
  All 404 local and project backup DVC objects were checksum-verified.
- See [protocol](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/interaction_ablation.md),
  [RUN record](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-188.yaml)
  and [launch verification](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-188.launch_verification.json).
  No production scientific result is available at launch.

## 2026-11-25 ANA-072 — RUN-186 TrpB behaviour, architecture and difficulty proxies

- Created the requested figure-led analysis under EXP-065/RUN-186: four notebooks,
  seven figures with PDF copies, reviewable Python sources and numerical tables.
  Fresh-kernel execution and independent mathematical/scientific review pass.
- RUN-186 retains all three RUN-180 warm starts after 35,136 joint Optuna trials.
  The shared-stop observed setting partially resembles TrpB: median mean 0.02274
  versus 0.02115, SD 0.06132 versus 0.05998; local SD and upper-tail mismatches remain.
  Shared-stop remains an observation sensitivity assumption.
- Equal-cube and sign-epistasis summaries look relatively close, but pooled
  higher-order variance is 0.202 versus 0.086 empirically. The latent model remains
  profile additive plus full Gaussian pair tables; no empirical coefficient
  inference or Gamma fitting was performed in RUN-186.
- At the primary resolution, shared-stop observed no-beneficial-neighbour and
  inaccessible-top-1% fractions are 0.595% and 0.958%, versus TrpB 0.514% and 1.190%.
  Intrinsic median landscapes are much less locally trapping; neutral connectivity
  still differs. These are static proxies, not measured learning or optimisation
  performance. No such tasks were run for this ANA.
- Primary summaries use eight originally fresh seeds; all 24 are retained. These
  panels/seeds have already been inspected, full-library comparisons include
  calibration nodes, and observed maps use error frozen by genotype. This is
  conditional evidence for a four-site analogue, not broad graduation or a new
  independent biological validation. No preset, result or claim was changed.
- RUN-186 has now been collected with its frozen runtime, DVC tracked and backed
  up; all 393 unique local and remote objects verify. Analysis reproduction and
  archive evidence are linked from
  [ANA-072](analyses/2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies/analysis.md).
  The completed ANA archive is DVC tracked and backed up; all 67 unique local
  and remote objects verify. Earlier status notes below are historical.

## 2026-11-25 RUN-186 completed successfully

- Slurm job 2355 finished on Dutch at 03:42 AEDT on 25 November. Scientific
  runtime was 7 h 44 min. Both scientific and runner exit codes are zero, stderr
  is empty, and all 86 manifest-listed output files match their SHA-256 checksums.
  labproj reconciles completion as INFERRED_EXIT_CODE because the scheduler no
  longer retains the job and Slurm accounting is disabled.
- All three studies completed 11,712 trials each (35,136 total), with no failed
  trials. None improved its historical RUN-180 warm start: no-noise r=4,
  loss 0.13755328; independent-difference r=2, loss 0.06074944; shared-stop
  sensitivity r=0.5, loss 0.03099517. This reports fixed-budget execution,
  not convergence or global optimality.
- All 24 confirmation seeds completed; outputs include 363 graph rows, 363
  cube-summary rows and 72 latent-graph rows. These diagnostics have not received
  a new scientific adequacy analysis. Full output archival/collection is not part
  of this status check; durable outputs remain on Dutch. No preset is changed.
- Completion evidence and checksum verification are retained in
  EXP-065/runs/RUN-186.completion-status.json and RUN-186.completion-integrity.json.
  Earlier progress snapshots below remain historical.

## 2026-11-24 RUN-186 progress at 20:08 AEDT

- Slurm job 2355 remains RUNNING on Dutch after 10 min 21 sec. The first
  observation study (no noise) has completed 867 of 11,712 trials, with no failed
  trials; 867/35,136 is 2.47% of the total optimization budget. Stderr is empty.
- Best calibration loss remains the historical warm start: 0.13755328 at r=4.
  No improvement has been observed yet. The other two observation studies and
  the 24-seed confirmation/graph diagnostics remain pending. Trial percentage
  is not overall job completion or an adequacy assessment.
- Timestamped queue, logs and read-only SQLite evidence is retained in
  EXP-065/runs/RUN-186.progress-20261124T090830Z.json.
  Slurm accounting storage is disabled; current RUNNING state is verified by
  squeue/scontrol, with no completion exit code or final manifest yet.

## 2026-11-24 Joint continuous-r Optuna fitting running as RUN-186

- Scientist authorized a new TrpB-only RUN under EXP-065, jointly optimizing
  continuous r in [0,8] and the four readout parameters with the restored RUN-180
  direct fitness loss. There is no r fitting grid and Gamma concentration stays
  disabled. Three studies preserve the observation assumptions and total 35,136
  trials; six historical fits per mode are warm starts. Old data, masks, seeds,
  noise draws, bounds and post-fit diagnostics remain fixed.
- Added the joint fitter, updated the default dispatch, and prepared an isolated
  Dutch Python 3.11.15 runtime with all 36 dependency pins. 114 targeted tests pass;
  the full local eight-trial smoke completes three settings/two confirmation seeds.
- Initial Linux preflight exposed only floating-point pair-field differences:
  maximum 1.33e-15 on calibration and 2.22e-15 across 24×160,000 confirmation values;
  additive values and references match exactly. The joint replay audit now allows
  rtol 1e-13 / atol 1e-14, with actual differences recorded, while retaining old loss
  replay tolerances and exact input authentication. The old strict audit remains
  the default for the superseded fixed-r adapter. No model values are altered.
- The complete 339-payload source/input runtime is frozen and DVC tracked under
  data/processed/trpb-joint-readout-slurm-runtime. The prelaunch amendment and
  numerical/environment evidence are retained in EXP-065/joint_readout_remote_environment.
  Dutch preflight passes all 18 historical objective replays (maximum component
  difference 6.94e-18), source/dependency verification and native graph checks.
  All 323 unique local and remote DVC objects verify. RUN-185 is the successful
  nonexecuting preview; RUN-186 is running on Dutch as Slurm job 2355, with
  4 CPUs, 8 GB memory and a 24 h limit. Worker startup reverified the frozen files,
  dependencies and all 18 historical objectives. Its first study has 89 completed
  trials including 83 off-grid r values; all five coordinates are present and
  stderr is empty. See EXP-065/runs/RUN-186.launch_verification.json.
  No completed calibration or model-graduation conclusion is available yet.

## 2026-11-24 RUN-180 calibration protocol restored with Optuna

- Scientist requested rollback to the most recent useful calibration while
  retaining Optuna. EXP-065 now defaults to `optuna-readout`, with its amendment
  in [optuna_readout.md](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/optuna_readout.md).
  No new scientific RUN was registered or launched, and no public preset changed.
- Restored the exact RUN-180 101-quantile, mean local-cube-SD and near-floor
  objective, affine units, parameter bounds, fixed r grid, all three observation
  modes, Gaussian profile full-pairwise fields and reference-anchored sigmoid.
  Optuna TPE optimizes four readout parameters separately for each of 18 settings,
  with 1,952 trials per setting including authenticated historical candidates.
  Every same-setting RUN-180 fit must replay and remains a candidate, so selected
  calibration loss cannot regress. This is not a confirmation-performance promise.
- Existing model features and the experimental architecture API remain available;
  the restored project protocol does not fit Gamma concentration or auxiliary
  coefficients. Historical RUN-180/RUN-184 archives and fitting functions remain
  intact. The shared pipeline gained explicit fitting/validation/provenance hooks.
- Engineering verification: 48 targeted tests pass; all 18 previous objective
  components replay exactly. Current model fields reproduce eight calibration
  seeds and all 24 confirmation seeds across 160,000 genotypes each. The bounded
  eight-trial smoke completes all 18 settings and two confirmation seeds; selected
  parameters and corresponding graph/cube diagnostics reproduce RUN-180.
  One additional historical test still assumes obsolete implicit ranked/synthetic
  defaults; its incompatibility predates this rollback and is recorded explicitly.
- [Validation](experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/optuna_readout_validation.json)
  records source hashes, replay details, commands, smoke scope and limitations.
  Previous confirmation panels have already been inspected; this rollback is not
  new independent biological validation, a model graduation, or a test of Optuna
  superiority. Shared-stop remains an observation sensitivity assumption.

## 2026-11-24 ANA-071 — RUN-184 does not yet support graduation

- Scientist requested a figure-led ANA with reproducible notebooks and an explicit
  assessment of TrpB landscape recovery. Created ANA-071 under EXP-065/RUN-184.
  Three fresh-kernel notebooks executed; five PNG/PDF figures show native fitness,
  independent local diagnostics, all 18 architecture targets, inference stability
  and the complete search. All seeds and failures are retained; no refitting.
- Selected model fitness SD is about 27 times too small; 87.72% of measured means
  lie outside its attainable output interval. Confirmation sign epistasis is
  15.4–19.4% versus 65.9%, reciprocal sign epistasis 0.7–1.8% versus 27.2%.
  Walks and peaks also differ on the identical sparse confirmation graph.
- The empirical coefficient model independently predicts confirmation observations
  at R²=0.720 (0.735 against genotype means). This predictive value does not carry
  over to the synthetic coefficient generator. All 20 bootstrap architectures shift
  radically from the original estimate despite numerical convergence, limiting
  interpretation of fixed weighting scales. The cause remains unresolved.
- Verdict: reproducible as a synthetic oracle, but this fitted setting is not
  ready as a TrpB-calibrated environment. No retrospective numerical gate is added,
  and no rejection of the full model family is inferred. Earlier favorable claims
  apply to different versions/protocols. No default, preset, RES or CLAIM changed.
- [ANA-071 report](analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/analysis.md)
  links figures, source notebooks, verification and reproduction commands. Archive
  and backup verification are recorded with ANA-071; original RUN-184 is immutable.

## 2026-11-24 TrpB RUN-184 completed and verified

- Job 2325 finished on Dutch at 16:00:32 AEDT after starting at 09:43:36
  (about 6 h 17 min). Both runner and scientific exit codes are zero; all 11
  expected outputs are present. Slurm no longer retains the job, so labproj
  records INFERRED_EXIT_CODE from runtime evidence, not a queried final state.
- All 100 Optuna attempts were processed: 37 COMPLETE and 63 FAIL due to
  auxiliary inference nonconvergence. Failed attempts and diagnostics are
  retained; they consumed the fixed search budget and limit usable search.
- Best trial index 63 has calibration loss 24.0341 and confirmation loss 27.3463
  on disjoint genotypes with four fresh seeds. Required empirical, successful
  trial and confirmation fits pass the strict convergence checks.
- Independent completion/reproducibility audit passes 33/33 checks: frozen
  inputs/source intact, SQLite and trial records consistent, losses recomputed,
  and native-parameter prediction replay agrees within 5.3e-18. All 35 remote
  output/log file hashes match the local collection. This establishes execution
  and reproducibility, not model adequacy or superiority; no preset is promoted.
- See RUN-184.completion-status.json and RUN-184.completion-verification.json.
  Outputs, logs, frozen runtime and transport evidence are collected under
  data/processed/aggregate-thermodynamic-gnk-calibration/RUN-184. Its 225-file
  DVC archive is backed up; all 178 unique local and remote cache objects verify
  (RUN-184.storage-verification.json). Earlier progress below is historical.

## 2026-11-24 TrpB RUN-184 — historical 15:29 progress check

- Latest check at 15:29 AEDT: job 2325 remains RUNNING on Dutch after 5 h 45 min,
  using approximately one CPU and 151 MiB resident memory. Empirical inference
  and all 20 bootstrap fits completed with convergence flags true.
- Optuna has finished 90 of its 100 attempted trials: 36 COMPLETE, 54 FAIL, with
  the 91st attempt running. All failures report auxiliary inference nonconvergence.
  Forty failed trials' full fits reached the 5,000-iteration cap; fourteen report
  function-tolerance termination, which alone does not establish convergence of
  the full fit and every selected-penalty CV fold. Failed attempts consume budget.
- The 20th successful trial finished at 13:22 AEDT, enabling adaptive TPE sampling.
  Best valid objective is now 24.03 (trial index 63), improved from 76.68 at the
  previous check. This is a calibration discrepancy, not an adequacy verdict.
  Confirmation and final outputs are not yet available. The 60% numerical failure
  rate still limits usable search; it does not itself establish biological inadequacy.
- Source, convergence requirements and budget remain unchanged. Queue, process,
  logs and database evidence is retained in RUN-184.progress-20261124T042858Z.json.
  Earlier progress snapshots and trial-failure evidence are preserved.

## 2026-11-24 TrpB RUN-184 — historical startup check

- The scientist repaired controller hostname resolution and reconfigured Slurm.
  Job 2325 / RUN-184 started on Dutch at 09:43:36 AEDT. Dutch was MIXED and
  Ripley IDLE after recovery; the later RESUME request was unnecessary.
- At 10:37 AEDT the calibration was actively computing: about 54 minutes elapsed,
  one thread using approximately 99% of a CPU and 110 MiB resident memory.
  Worker startup verified all 152 frozen payloads and 31 dependency pins;
  stderr was empty. No source, configuration, seed or fitting-budget change.
- Only initial input/configuration/provenance outputs exist. The frozen method
  writes its empirical architecture after the empirical fit and all 20 bootstrap
  fits; it has not yet created the Optuna study. Exact bootstrap progress and a
  completion percentage are not observable from current instrumentation. No
  numerical convergence or scientific adequacy conclusion is available.
- RUN metadata was reconciled with Slurm. See RUN-184.progress-20261124.json;
  earlier pending-state observations remain historical records.

## 2026-11-22 Slurm status — RUN-184 still pending

- At 14:04 AEDT, Slurm reports one queued job (2325 / RUN-184), PENDING
  with zero elapsed execution and no running jobs. Dutch and Ripley are reachable
  by SSH with active slurmd services, but both are DOWN+NOT_RESPONDING in Slurm.
- The raspi-1 controller is active. Both workers carry NO NETWORK ADDRESS FOUND;
  controller hostname lookups return no worker addresses. This indicates a
  controller-side name-resolution problem. No cluster settings, services, job
  configuration or node states were changed. See RUN-184 queue-status-20261122.

## 2026-11-21 TrpB calibration moved to Slurm — RUN-184 queued

- Scientist requested stopping local RUN-182 and queuing the same calibration.
  The verified worker received SIGTERM; worker and supervisor exited. RUN-182
  is cancelled, with its 190-file immutable partial collection DVC-tracked and
  backed up. All 149 unique local and remote cache objects verify. Its original
  automatic failure record is retained; missing final results reflect user
  interruption, not evidence of numerical or scientific failure.
- RUN-183 passed the nonexecuting Slurm preview. RUN-184 was accepted as Slurm
  job **2325** on the CPU partition (1 CPU, 8 GB, 24 hours). At 2026-11-21
  01:52:44 UTC it was PENDING because the required worker was unavailable;
  Dutch is reachable by SSH but marked DOWN in Slurm. No cluster service or node
  state was changed. The job has not begun scientific execution.
- The self-contained runtime preserves all 152 RUN-182 prelaunch payloads and
  its original manifest; all staged scientific bytes, dispatcher and environment
  evidence verify remotely. Its DVC bundle and all 131 unique local/remote cache
  objects verify. Source, inputs, fitting budget and seeds are unchanged; this
  is a fresh run, with no unsupported Optuna resumption.
- A dedicated Linux Python 3.11.15 environment matches all 31 original package
  pins; imports and structural input/configuration validation pass. Linux/BLAS
  differences are recorded; numerical bitwise identity with Mac execution is not
  assumed. Queue acceptance is not completion or an adequacy claim.
- See EXP-065 `architecture_slurm.md` and RUN-184 remote verification. The
  original dispatcher, project overlay, SSH transport and stage rules were
  restored after submission. Final outputs await queued execution.

## 2026-11-21 TrpB Optuna architecture calibration — initial RUN-182 launch

Historical launch record below; RUN-182 was subsequently cancelled and replaced
by queued RUN-184 as recorded above.

- Scientist authorized a new RUN under existing EXP-065 and local MacBook
  execution of the 0.7.0rc1 indirect-inference calibration, restricted to TrpB.
- Prepared deterministic, disjoint 8,192-genotype calibration and confirmation
  panels from the historical 159,129-genotype mask, retaining both actual assay
  replicates. The fixed profile matrix is unchanged. Every site-pair/allele cell
  is observed; the split uses genotype IDs and a fixed seed, not fitness values.
- `architecture_calibration.md` records the changed objective, six controls,
  noiseless primary synthetic observations, strict inference convergence and
  historical-data limitations. A calibration-only timing/convergence pilot is
  retained before freezing the execution budget. No confirmation result is used
  for preparation, model selection or budget selection.
- The source lock and 29 CLI/input/transport checks pass. Raw author replicate
  columns, disjointness, every pair/allele combination, grouped folds and replay
  were verified. The 19-file input bundle is DVC-tracked and backed up; all 20
  local and remote cache objects verify.
- Calibration-only pilots showed that the example 1,000-iteration limit and
  function tolerance could stop before strict stationarity. The final protocol
  retains the 1e-5 stationarity requirement, raises max iterations to 5,000,
  tightens function tolerance to 1e-10 and allows 10 curvature resets. The penalty
  grid retains the four example values and adds 10 and 100 after the original
  optimum hit its boundary. Full empirical, first bootstrap and fixed synthetic
  fits now converge; all earlier attempts are preserved, including the pilot
  adapter's initial uint8-index failure (production already converts to int64).
- The run retains 100 Optuna trials, 4 calibration seeds, 4 fresh confirmation
  seeds and 20 bootstrap weight estimates. Its 10-fold, six-penalty inner search
  entails up to 25,986 optimizer fits. Concurrent pilot timings are approximately
  104–138 seconds per complete auxiliary estimation; an overnight-scale workload
  is anticipated, not a guaranteed finish time or performance benchmark.
- Preview RUN-181 passed. Live RUN-182 started locally at approximately 12:39
  AEDT, using a frozen 124-file scientific source/protocol/preparation manifest
  and complete input snapshot. Actual output configuration, initial input/source
  hashes, continuing supervisor heartbeat and restored project overlay were
  verified. The native Python 3.11.15 environment and dependency versions are
  captured. Idle sleep is inhibited while the detached supervisor runs.
- The 45-file preparation bundle is also DVC-tracked and backed up. The supervisor
  will preserve success or failure artifacts in the RUN-182 bundle; final DVC
  tracking remains pending until outputs are immutable. See RUN-182 launch
  verification for paths. No ANA/RES/CLAIM, completed calibration, preset update,
  publication or Git push is claimed. The locked 0.7.0rc1 source and historical
  RUNs remain unchanged. Final strict refresh/check pass with 376 green,
  22 historical/prospective yellow and zero red checks; whitespace checks pass.

## 2026-11-20 Gamma concentration and Optuna architecture fitting implemented

- Scientist authorized Gamma-normal concentration and the proposed Thornton-
  anchored indirect-inference fitting method. Unpublished candidate 0.7.0rc1
  adds `interaction_concentration` and seed-independent `global_intercept`.
  Alpha changes per-entry variance concentration; r retains exact mutation-RMS
  meaning. None retains the Gaussian/reference-occupancy construction exactly.
  No arbitrary finite alpha is promoted to a calibrated default.
- New `adaptenv.inference` follows the pinned authors' sparse sigmoid SSE+L1
  objective, grouped cross-validation and post-fit zero-sum gauge. The Python
  split-variable solver differs from their OWL-QN executable; source mapping,
  original MIT notice and explicit differences are documented. An independently
  executed upstream R objective/gradient fixture agrees. Stationarity checks
  prevent premature function-tolerance stops being reported as converged.
- New `adaptenv.calibration` fits r, alpha, steepness, mean-logit position, lower
  and span with serial multivariate Optuna TPE. Empirical and synthetic panels
  use the same inference rule, masks and genotype-grouped folds. Eighteen
  architecture/readout summaries use fixed calibration-only genotype-bootstrap
  weights. Seeds are never fitted. Noiseless and explicit pooled-replicate-
  residual sensitivities are supported; historical shared-stop is not implied.
- Artifacts retain inputs, source snapshots, versions, weights, coefficients,
  every trial/failure, boundary distances, native-unit oracle parameters and
  fresh-seed confirmation. Separate confirmation genotypes must be disjoint;
  tests prove changing confirmation labels does not alter calibration. Unfitted
  local sign epistasis, diminishing returns and walks concern the observed
  panel's induced graph, not an invented whole-library graph. Near-best trial
  ranges are descriptive, not posterior/confidence regions.
- Source plus CLI: 239 tests pass. Clean installed wheel: 228 tests pass; no
  skipped tests in either Python 3.13 verification. Documentation examples and
  upstream-equivalence, failure, leakage and native-parameter replay checks pass.
  Same-interpreter comparisons preserve all 28 checked old maps/fingerprints
  against frozen RUN-180 source. Seven old fixed fingerprints differ under
  Python 3.11 because of Python float-sum semantics; frozen source has the same
  differences and the unchanged original fixtures pass under Python 3.13.
- Engineering sampling with alpha=0.2 at 1,000 sites retains 3,000 full pair
  tables (9.6 MB) and scores 256 sequences in about 0.01 s locally. This is not
  an HYP-010 comparison or long-sequence biological validation. Offline all-pairs
  inference/CV is much more expensive: example optimization entails 16,400
  auxiliary fits before bootstrap/confirmation, so it requires a timing pilot.
- `docs/calibration` supplies the exact protocol, validated CLI and example;
  `docs/model_release/0.7.0rc1` preserves engineering evidence and source lock.
  Separate Python 3.11/3.13 validation dependency snapshots leave RUN-180's
  environment unchanged. Historical candidates and scientific outputs remain
  frozen. No new empirical RUN, fitted TrpB parameters, adequacy claim, package
  publication, deployment or Git push was performed.
- The 69-file candidate archive is DVC-tracked locally (517,118 bytes; directory
  MD5 `8d7dfd3d3c5662e38fe8ba5008c87d2a.dir`). All archive payloads and 70 cache
  objects verify. The 67-file source/verification lock passes. Labproj also
  completed its configured remote backup; all 70 remote objects verify under
  the available DutchZFS volume. The archive pointer is staged by labproj; no
  commit was made. Strict refresh/check pass: 374 green, 22 existing/prospective
  yellow and zero red; source and project whitespace checks pass.

## 2026-11-20 matched full-pairwise TrpB calibration — RUN-180 completed

- Scientist authorized a new RUN under EXP-065 and local MacBook execution of
  the refined 0.6.0rc1 model. `full_pairwise.md` / `full_pairwise.yaml` preserve
  RUN-178's profile, data/masks, eight calibration seeds, 24 confirmation seeds,
  18 settings, optimizer rules, observation models and unfitted diagnostics.
  Only the interaction construction changes; RUN-178 is the direct control and
  synthetic RUN-175 is retained without rerunning either control.
- Local preview RUN-179 passed; live RUN-180 completed with return code zero in
  330.75 seconds using the same pinned Python/numerical environment and compiler.
  All 18 settings, 24 seed audits, 2,163 graph rows, 2,163 cube-summary rows and
  144 latent graph rows are present. The full budget was used, with four workers.
- The 84-file source lock and 133 input/control files verified before execution;
  all 50 focused preflight tests passed. Independent post-run engineering checks
  pass. All 351 checked numerical invariants match exactly, including additive
  values, references, empirical inputs and r=0 fits/graph/cube outcomes. The
  intended interaction fields differ. The original scientific sources remain
  unchanged; the 51-file package candidate lock still passes.
- Initial verification flagged the noise membership checksum: gzip timestamps
  changed while decompressed donor data remained byte-identical. The original
  failure report is retained. An explicit v2 verifier authenticates each file's
  own hash and requires identical decompressed bytes and remaining model fields;
  six new and seven retained verifier tests pass. See RUN-180.verification.json
  and RUN-180.verification_amendment.json in EXP-065/runs.
- Source, complete inputs, environment, logs and outputs are collected in the
  DVC-tracked RUN-180 bundle: 346 files, 119,576,625 bytes; 326 unique cache objects
  and all collected SHA-256 hashes verified locally. Remote backup was not
  attempted. Post-run verification reports and correction are recorded alongside
  the RUN; the collected bundle retains its original prelaunch sources.
- No ANA/RES/CLAIM, adequacy conclusion, package preset recalibration, release or
  Git push was created. Prior calibration evidence and the locked candidate are
  retained unchanged. Final strict labproj refresh/check pass: 374 green,
  22 existing/prospective yellow and zero red; source and record whitespace
  checks pass.

## 2026-11-19 full-pairwise default — unpublished 0.6.0rc1

- Scientist requested full q20 pair tables and pairwise-only latent terms as the
  default, anchored to Park–Metzger–Thornton's interaction expansion plus sigmoid.
  The profile additive prior, exact mutation-RMS meaning of r, reference anchoring
  and sigmoid are retained. Pair tables are centered Gaussian matrices with all
  361 pure-pair degrees of freedom; there is no rank-eight/binary restriction.
- Python, service and HTTP defaults now identify `thermodynamic_pairwise_v3`.
  Bounded site degree six is retained for tractability; each selected pair has
  a full table. At four sites all six pairs are present. All-pairs topology can
  be selected explicitly, with quadratic rather than linear length costs.
- All six calibrated TrpB presets explicitly retain `ranked_gnk`, original
  orders/weights/degrees/rank and `thermodynamic_gnk_v2`. Regression checks across
  four seeds per preset preserve components, mutation effects, references and
  fingerprints exactly. The RUN-178 packaged profile and historical evidence
  remain unchanged. The new default is uncalibrated; no RUN was created or run.
- 148 source tests and 148 clean installed-wheel tests pass, plus 11 documentation
  examples and all 13 project-facade exports. Exhaustive q20 checks at two and
  three sites verify exact normalization, including overlapping pairs. A second
  implementation review found no material mathematical or integration defect.
- A reproducible local engineering check evaluated 256 diverse uniform-q20
  sequences at L=100/500/1,000. At L=1,000: 3,000 pair tables, 9.6 MB table payload,
  about 0.104 s construction and 0.0105 s batch scoring; whole-process peak RSS
  was about 88.5 MB. This is not a comparative HYP-010 benchmark or long-sequence
  biological validation. Fixed generic readout settings can saturate at large L.
- The new decision, reproducible check, verification and 51-file source lock are
  in `docs/model_release/0.6.0rc1`. The previous candidate lock/archive is retained
  as historical provenance. No package release, tag, deployment or Git push was
  performed. EXP-067 still targets the frozen ranked candidate; its disabled
  harness requires a protocol amendment before benchmarking this new default.
- The 64-file candidate archive is DVC-tracked and its local cache hashes verify;
  remote backup was not attempted. Strict labproj refresh/check pass with 372
  green, 22 existing/prospective yellow and zero red checks. The external-package
  lock accurately records the source checkout as modified; the candidate manifest
  identifies the tested bytes. Source and record whitespace checks pass.

## 2026-11-19 HYP-010 computational experiment registered

- Created planned EXP-067 through labproj, linked to AIM-001, HYP-010, EXP-065 and
  RES-063. Four arms compare the frozen factorized oracle with its exact local-table
  expansion, and conventional GNK tables with an identical counter-based on-demand
  map. Ehrlich is excluded. No RUN, execution, performance result or ANA is created.
- The protocol/config and complexity_notes.md define n=q^L as all possible complete
  sequences and Q as actual queries. At fixed alphabet/order/degree/rank all arms
  have O(log n) storage/full-query cost. NK enumerates local tables, not necessarily
  all n states. Current support sorting/retries and retained Gram matrices are
  explicitly included in setup/memory derivations; no unconditional linear setup
  or universal asymptotic advantage is claimed.
- The design fixes four lengths, ten paired seeds and three process repetitions,
  a 9,600-query workload, memory/time envelopes, correctness gates and comparison
  rules. A/B have exact realized normalization; C/D have ensemble-expected
  normalization. Cross-pair capability differences limit the resulting claims.
- This is a prospective protocol, with a guarded run.sh. Benchmark adapters,
  counter implementation, isolated dependency lock, host and resource monitor must
  be frozen before a RUN. The historical project overlay is disabled for this
  experiment. The unpublished model candidate and its source lock are unchanged.
- Labproj refresh updates external-package clean-state flags and the derived lock
  header; it does not change pinned commits or resolve the historical environment.
- Strict refresh/check pass (372 green, 22 historical/prospective yellow, zero red).
  Protocol/config checks validate pairing, normalization algebra, table-byte counts,
  input hashes and links; shell syntax is valid and all 45 candidate-lock files
  still verify. No timing, correctness benchmark or model recalibration was run.

## 2026-11-19 result synthesis, benchmark hypotheses and Git synchronization

- Created completed RES-063 through labproj as a synthesis of ANA-069/070 and
  the current implementation. It records usable TrpB regimes and precise latent
  controls, with mixed profile trade-offs and no new empirical experiment.
  Existing figures/notebooks remain canonical; metrics.json selects 48 existing
  ANA-070 aggregate rows with the source hash. CLAIM-004 now links this result
  without expanding its original synthetic-prior claim.
- Created active, untested HYP-010 (computability and comparative efficiency),
  HYP-011 (separable learnability/searchability regimes), and HYP-012 (biological
  information and conditional learning/search transitions). Each states task
  definitions, controls and falsifiers. AIM-001 records the scientist-authorized
  extension; real biological explanations require independent measured fitness
  evidence. No experiment or run is created or submitted by this change.
- The package audit passes 82 source and 82 isolated installed-wheel tests,
  verifies all 13 project facade exports, and confirms the 45-file release lock.
  The distribution contains only maintained modules; retired model families are
  absent. Historical project research helpers remain separate for frozen-run
  reproduction. An unrelated old editable hook in base Python can expose another
  checkout's legacy modules, so isolated wheel checks establish the distribution
  boundary. RES-063/package_audit.json and its logs preserve the fresh audit.
- Adapt-Env commit 745c00b3a369da33e7d7c7b9850ccbc212fac714 is pushed to origin/dev,
  including the earlier consolidation commit and the profile default. The project
  gitlink and external-package lock point to that clean revision. Git synchronization
  of this scientific project is authorized on its existing main branch.
- Publication preparation excludes local .labproj caches and duplicate experiment
  outputs/logs. Thirty small historical DVC pointers were restored to cited original
  paths from identical .legacy copies (29 already indexed plus RUN-086); raw historical
  payloads were not moved. The release-candidate DVC pointer is explicitly included.
  All 207 cloud-only record files and 1,842 historical Git objects were materialized;
  record file sizes were preserved and the candidate contents were scanned before
  staging. Strict record checks pass with zero red records (370 green, 23 yellow).
  Commit checks preserve generated evidence bytes and apply a 10 MiB size limit
  to versioned analysis figures/notebooks/tables; ordinary files retain the
  original size check. Raw experiment payloads remain DVC-managed.
- The full project research overlay remains a historical, incomplete dependency
  snapshot: its lock omits Landscapy's install entry/closure, and STACK-000 Python
  3.11 conflicts with the pinned Landscapy requirement >=3.12. This is now explicit
  in environment/README.md and the lock header. The tested standalone model and
  frozen RUN environments are separate; no untested full-stack upgrade was made.
- Git synchronization does not publish a package release, create a tag or deploy a
  service. DVC remote backup remains a separate pending operation while Dutch is
  unavailable. Existing historical labproj warnings remain visible; new hypotheses
  intentionally have no linked experiment yet.

## 2026-11-19 profile-conditioned default — unpublished 0.5.0rc1

- Scientist selected the alignment/profile additive variant as the model default
  because it carries biological marginal preferences. ANA-070 records mixed
  descriptor trade-offs, not consistent superiority or statistical equivalence;
  no new empirical claim or ML/pLM benefit is inferred from this choice.
- Adapt-Env is now locally versioned **0.5.0rc1**, an unpublished release candidate.
  The Python aliases, service and HTTP API require alignment/profile input by
  default. Gaussian additive effects require explicit `additive_prior="synthetic"`
  or a named synthetic control. Advanced supplied matrices/effect samples retain
  accurate source labels and are not automatically called biological profiles.
- The new `trpb_run178_r05_profile_shared_stop_intrinsic` preset bundles the exact
  four original frequency rows and matching RUN-178 readout. The full factory
  loads both; the config-only factory cannot discard the profile. All calibrated
  presets reject non-seed parameter, profile and reference substitutions. The
  synthetic-prior RUN-170/RUN-175 presets remain explicitly identified controls.
  Runtime noise remains off, and reference anchors remain computational.
- All 82 maintained package tests pass from source and from the installed clean
  wheel; 11 documentation examples and the project API delegation check pass.
  A first wheel exposed obsolete modules in an ignored historical build cache;
  the final wheel was built from clean source staging and contains no retired
  modules. The old build cache was preserved under `.labproj/tmp`.
- Pinned-Python replay covers all 24 recorded seeds for each prior, with 1,020
  latent sample sequences and 128 recorded cube sequences per seed. Current
  intrinsic values match the frozen implementations exactly on these queries;
  differences from saved outputs are at most 2.23e-15. Configuration fingerprints
  change because the prior-selection policy is now explicit. No calibration or
  empirical simulation was rerun. Historical RUN/ANA artifacts remain unchanged.
- `docs/model_release/0.5.0rc1/decision.md`, `manifest.json`, `verification.json`
  and `replay.json` record the decision, 45-file source/evidence lock and checks.
  `lock.py --check` verifies the candidate. The 50-file source/wheel/evidence
  archive is DVC-tracked at `data/processed/model-release-candidates/0.5.0rc1.dvc`
  with local payload/cache hashes verified; remote backup awaits Dutch.
  No release, tag, publication or deployment was performed.

## 2026-11-19 TrpB profile comparison — ANA-070

- Completed the requested figure-led analysis of local RUN-178 against the
  recorded RUN-175 synthetic control under EXP-065. ANA-070 contains five
  PNG/vector-PDF figures, two executed notebooks (11 code cells, five inline
  figures), paired-source Python files, full comparison tables and independent
  numerical/provenance audits. The synthetic control was not rerun.
- At the inherited r=0.5/shared-stop setting and resolution 0.001, fresh-seed
  median walk length changes 7.459 to 8.749 versus measured 8.589, but only 2/8
  paired seeds improve. Pooled high-order power fraction changes 0.176 to 0.087
  versus 0.086 measured, with 4/8 seeds closer. Equal-cube epistasis, sign
  epistasis and diminishing-returns medians move farther from TrpB. These are
  mixed trade-offs, not consistent profile superiority. All strengths, noise
  assumptions, seed cohorts and resolutions remain visible in figures/tables.
- The profile supplies biological amino-acid preferences followed by the same
  readout calibration procedure. It does not add biological interaction terms
  or replace seed-specific computational reference anchors with wild type.
  Shared-stop observations remain a sensitivity assumption; fresh seeds do not
  provide independent biological validation. Whole-library graph descriptors
  include calibration nodes. No ML/pLM or active-learning performance was tested.
- Audits verify 238 output hashes, all 36 selected objectives, 42 finite-budget
  fit traces, 4,326 independent WHT summaries and 24 paired full interaction
  fields. Every optimizer search exhausted the original 60-iteration budget;
  no convergence or optimality claim is made. macOS/Linux platform differences
  and rounded-profile source limitations remain explicit.
- The 45-file analysis archive is DVC-tracked at
  `data/processed/aggregate-thermodynamic-gnk-calibration/ANA-070.dvc`, with all
  collected payloads and local cache objects verified. Remote backup awaits
  Dutch. EXP-065 follow-up is completed and links ANA-070. No RES or CLAIM was
  requested or created; CLAIM-004 retains its synthetic-prior scope.

## 2026-11-19 authorized local TrpB profile execution

- Scientist requested execution on the local MacBook while retaining the
  scientific record. This supersedes the preparation-only hold for the local
  target; Slurm submission remains held and RUN-176 stays a nonexecuting preview.
- RUN-178 completed the full-budget local run under EXP-065 in 416.51 scientific
  seconds (return code 0), after nonexecuting local preview RUN-177, with unchanged frozen source,
  inputs and calibration settings. An isolated Python 3.11.15 environment uses
  the five exact recorded numerical package versions. macOS arm64 and Apple
  Clang replace the Linux/Zig execution platform; this is recorded explicitly
  and is not asserted to give bit-for-bit cross-platform equivalence.
- Local execution used labproj RUN records plus a separately saved prelaunch
  source/input snapshot, captured logs and manually verified collection because
  the current local target lacks automatic source staging/collection support.
  The eight local Python 3.11.15 runner/input tests and 20 freshly compiled
  graph-kernel cases passed. RUN-178 verified its source/input hashes and exact
  numerical versions on startup. All 18 settings, 24 seeds (16 reused + 8 fresh),
  2,163 graph rows, 2,163 cube rows, 144 latent rows and 99 output hashes were
  independently verified. All 24 complete normalized GNK fields match RUN-175
  within 4.45e-15, and all 38 empirical graph numerical columns match exactly
  at three resolutions. The synthetic control was not recalibrated.
- The 208-file source/input/environment/log/output bundle is DVC-tracked at
  `data/processed/aggregate-thermodynamic-gnk-calibration/RUN-178.dvc`, with local
  cache and collected-file hashes verified. Detailed execution, protocol and
  storage audits accompany the RUN record. PyArrow CPU-feature/cache probe
  permission warnings are retained in stderr; execution and numerical checks
  passed. Remote backup remains pending while Dutch is offline. Local execution
  fulfills the requested run; no duplicate Slurm job was submitted. No new
  empirical adequacy assessment, ANA, RES or CLAIM is created by completion.

## 2026-11-19 TrpB alignment-profile repeat preparation

- Scientist authorized restoring the profile-based additive component and a
  matched repeat under EXP-065, reusing the recorded synthetic control. The
  comparator is the TrpB branch of RUN-175, with its RUN-172 endpoint-calibration
  lineage. Six strengths, three observation modes, all original calibration and
  evaluation seeds, masks, cube geometries, donor pools, optimization settings
  and unfitted diagnostics are preserved; only the additive prior changes.
- The inherited r=0 and r=4 readouts require the original endpoint fitting to be
  repeated under the profile prior before anchoring intermediate-ratio fits.
  Retaining synthetic-prior endpoint parameters would not repeat calibration.
- Implemented optional raw-alignment and frequency-profile additive adapters in
  the Python package, with frequency-profile input through the public service
  and HTTP API. The original computational methods and synthetic defaults are
  unchanged. The matched runner freezes 42 source files, the original numerical
  kernels, input identities and full calibration settings by SHA-256.
- The original-study TARGET_b0.3 A2M archive member contains only NUL bytes.
  The frozen prior therefore uses the authors' complete alignment-frequency
  table, rounded to 0.001, at physical sites 183/184/227/228 (parent VFVS).
  Non-gap conditioning and a fixed 1% uniform probability mixture precede
  log-frequency conversion; weighting and effective sequence count are unknown.
  Source defects, alternative exclusions and exact mapping are recorded in
  `data/processed/trpb-profile-additive-inputs/prior_manifest.json`.
- All 49 maintained package tests, five matched-run tests and three profile-input
  tests pass. The actual-input SMOKE_ONLY execution passed all 18 settings with
  two seeds in 40.81 seconds, including 20 independent graph-kernel checks.
  It is engineering evidence only. Scientific startup requires the exact
  recorded RUN-175 Python/numerical package versions; local smoke used Python
  3.13.11. No empirical result is available and CLAIM-004 remains unchanged.
- RUN-176 is the prepared, nonexecuting labproj dry-run under EXP-065: four CPUs,
  8 GB and one hour, with both input bundles included in source staging. The
  profile bundle is DVC-tracked and its 12 local cache objects verified. Remote
  backup and worker-side environment/compiler checks await Dutch restoration.
  The preparation record is `experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-176.preparation.json`.
- HomeLab gateway is reachable, but Dutch reports no route to host and the
  gateway lacks the shared runtime paths. Scientist confirmed Dutch is down and
  instructed preparation only, with submission to follow after restoration.
  No job has been submitted for this follow-up.

## 2026-11-19 TrpB package-readiness and literature review

- Follow-up source/history audit: the synthetic-prior thermodynamic candidate
  was introduced on 7 October while the previous public environment remained
  available. Commit dd0bb7a (8 October, 11:17 AEDT, v0.3) switched public aliases
  and service to that candidate, moving alignment/query conditioning to private
  legacy code. Commit 37d23c6 (8 October, 22:25 AEDT, v0.4) removed that archive;
  the latter cleanup did not newly change the candidate's conditioning.
  This was a scope/architecture transition, not a demonstrated rejection of
  additive consensus/conservation information. RES-035--039 constrain PHOT
  module specificity; CLAIM-003 explicitly retained biological side information
  for PLM benchmarking. HYP-003's direct conditioning question remains unsettled.
  The documented authorization for removing retired methods does not by itself
  establish explicit approval of the loss of that benchmark capability.
- A minimal proposed restoration is a reweighted, pseudocount-smoothed MSA
  profile h_i(a)=log p_i(a), supplied through the existing additive_matrix
  constructor argument with a mapped natural or consensus reference. Existing
  global mutation-RMS normalization preserves relative site conservation; no
  per-column standardization is proposed. This retains current sparse GNK
  scoring and sigmoid, with synthetic interactions remaining synthetic. It is
  an evolutionary-preference prior, not measured free energy or assay fitness;
  changed priors do not inherit the existing TrpB calibration automatically.
  Alignment and consensus utilities remain installed, although a stale module
  docstring still implies automatic landscape conditioning. No adapter, API
  change, empirical run or revised scientific claim was implemented in this audit.
- Follow-up advice on moving to active learning: consolidate the established
  model-family positioning; quantify initialization, memory, full scoring and
  mutation-batch performance against bounded-neighbourhood GNK if claiming an
  efficiency advantage; freeze the generator, seed distribution, observation
  semantics and scored target. Reproduce the advertised TrpB regime through
  that interface if retaining its observed-map claim. These are bounded closure
  tasks, not a requirement for another broad empirical-fitting campaign.
  Algorithm discrimination and ranking robustness belong to the next stage.
  Existing HYP-009 concerns retired biological routing/module constructions and
  should not silently govern the consolidated generator's active-learning test.
  Random synthetic amino-acid effects do not supply the evolutionary structure
  needed to interpret PLM performance as biological competence. This is advice
  only: no new protocol, run, hypothesis or model implementation was created.
- Scientist clarified that the requested single-system review concerns Johnston
  et al.'s TrpB, not Papkou et al.'s DHFR. Reviewed AIM-001, the hypotheses,
  historical results, CLAIM-004, ANA-063--069, current package source and primary
  literature. This is an assessment, not a new experiment or scientific claim.
- Assessment: the evidence supports freezing a scoped synthetic research package
  for landscape analysis and ML calibration against its declared synthetic ground
  truth. It does not establish transfer of ML uncertainty or optimizer rankings
  to experimental TrpB. RUN-170 supplies genotype-disjoint, withheld-descriptor
  checks; RUN-175's highlighted r=0.5 regime remains a post-sweep illustration.
- The strongest empirical example includes the shared-stop observation
  sensitivity. The packaged preset contains its intrinsic map only, and public
  Gaussian observations are a different process. Static graph observations are
  frozen per genotype; independent noise redrawn on online queries changes the
  learning/search problem. CLAIM-004's pooled-spectrum mismatch, seed variation
  and four-site/length-transfer restrictions remain in force.
- The existing Park et al. interaction-expansion plus sigmoid relationship
  remains the closest whole functional-form precedent. Low-rank factors specify
  a constrained ensemble, not a lossless acceleration of independent-table GNK.
  O(L) evaluation at fixed interaction complexity does not establish a speed
  advantage over existing bounded-neighbourhood implementations.
- Read [Stanton et al. (2024), Ehrlich functions](https://arxiv.org/html/2407.00236v1):
  Markov feasibility constraints and quantized spaced motifs define a different
  synthetic optimization geometry with a known attainable optimum. It supports
  the legitimacy of useful benchmark abstraction without molecular reconstruction,
  not validation of AdaptEnv. [SLIP](https://github.com/google-research/slip)
  already supplies MSA/Potts-based synthetic landscapes, and
  [FLEXS](https://arxiv.org/abs/2010.02141) already supplies model-guided sequence
  exploration infrastructure. The online-environment category is established;
  the candidate contribution is the particular controllable generator and its
  documented computational and empirical properties.
- Fresh verification: all 33 maintained Adapt-Env package tests pass on the
  current 37d23c6 checkout in 6.18 seconds, using installed dependencies and
  disabling unrelated pytest plugin autoload. No package code was changed.
  Labproj status, strict refresh and strict check completed; refresh reported
  green=363, yellow=21, red=1 and many historical data-tracking warnings, so this
  is not a claim of clean repository-wide compliance. No new EXP/RUN/RES/CLAIM.

## 2026-10-27 scientist prioritizes implementation novelty and sparse exploration

- Scientist clarified that the intended contribution is computational efficiency
  across large, sparsely sampled sequence spaces, with intermediate/phylogenetic
  adaptation as an application, rather than a new biological landscape class.
  Recorded this in AIM-001 without replacing its existing empirical constraints
  or creating a hypothesis, experiment, result or claim.
- Current model fits that mathematical scope: one coherent deterministic map on
  all q20 sequences of a configured length, no reference-distance evaluation
  cutoff, and no requirement to enumerate or observe nearby genotypes. Sparse
  exploration and a sparse interaction representation are different concepts.
- Source inspection confirms O(L) field storage and full scoring at fixed
  alphabet, maximum order, degree and rank. Mutation-effect batches use incident
  terms after one O(L) baseline calculation; the existing API recomputes that
  baseline per call, so no persistent O(1)-per-step walk implementation is claimed.
  Initialization has additional support construction/sorting costs.
- Bounded-neighbourhood NK/GNK already allows on-demand evaluation; avoiding
  enumeration alone is not new. Rank restrictions reduce local tensor storage
  from q^k entries to O(R*k*q) factors but alter the random landscape ensemble.
  They define an exactly evaluated subclass, not a lossless implementation of
  unrestricted GNK. Computational novelty and comparative speed are not proven
  by the functional-form relationship or asymptotic O(L) scaling alone.
- The implementation can supply a stationary fixed-length fitness map for walks
  or branching simulations across distant backgrounds. This is not itself a
  phylogenetic/population model or evidence of biological validity over those
  distances. Existing saturation and length-transfer limitations remain relevant;
  mutation-RMS normalization does not guarantee stable observable dynamic range.
  No source implementation or empirical evidence changed in this discussion.

## 2026-10-27 closest published formulation and literature foundation

- In response to the scientist's literature question, comparison beyond the
  package reference list identifies [Park, Metzger and Thornton (2024),
  The simplicity of protein sequence-function relationships](https://doi.org/10.1038/s41467-024-51895-5)
  as the closest complete functional form among the papers examined. Their
  reference-free categorical interaction expansion followed by a bounded sigmoid
  (Figure 1b) contains the current intrinsic equation as a constrained instance:
  their genetic score s corresponds to logit(p0)+beta*(a+r*g).
  This is functional-form inclusion, not identity of inference methods or random
  ensembles. Their main empirical fits emphasize additive/pairwise terms; their
  general expansion allows higher orders.
- AdaptEnv specifies a procedural ensemble using sparse supports, default orders
  2-4, rank-8 categorical products and mutation-operator normalization, whereas
  Park et al. infer coefficients from measured genotype-phenotype data. The broad
  construction of specific interactions inside a sigmoid therefore has published
  precedent. Novelty of the particular generator or its empirical calibration is
  not established by this literature comparison.
- [Carlson, Andrews and Simons (2025)](https://doi.org/10.1073/pnas.2509444122)
  explicitly combine additive and pairwise latent effects with a monotone readout
  in Equations 2a-b and use two-state thermodynamic simulations with specific
  interactions. [Faure and Lehner (2024), MoCHI](https://doi.org/10.1186/s13059-024-03444-y)
  provides another inference framework combining nonlinear mappings and explicit
  interactions, including higher-order terms.
- Foundational roles remain distinct: [Brookes et al. (2022)](https://doi.org/10.1073/pnas.2109649118)
  for GNK random fields and interaction spectra;
  [Otwinowski, McCandlish and Plotkin (2018)](https://doi.org/10.1073/pnas.1804015115)
  for additive latent-trait global epistasis;
  [Otwinowski (2018)](https://doi.org/10.1093/molbev/msy141)
  for thermodynamic folding/binding readouts and the limits of a single two-state
  coordinate; [Reddy and Desai (2021)](https://doi.org/10.7554/eLife.64740)
  for interaction-driven diminishing returns under stated statistical assumptions.
  These papers do not establish physical free energies, biological support
  topology, the chosen order weights or arbitrary-length empirical realism for
  AdaptEnv. Existing CLAIM-004 scope is unchanged; no new scientific record or
  experiment was created.

## 2026-10-27 formal classification of the current model

- Scientist requested classification and analytical relationships, including RMF.
  Inspected the unchanged Adapt-Env 0.4.0 source at
  `37d23c64b15a26cafa41232940feabc218505cd5`. The intrinsic model is a sigmoid
  transform of an additive field plus a sparse, low-rank, mixed-order categorical
  interaction field. It has a generalized multibody Potts/spin-glass Hamiltonian
  representation and a GNK-type local-factor structure; its constrained contrasts,
  rank and exact energy normalization do not reproduce the usual independent
  Gaussian lookup-table GNK ensemble.
- With p=(f-f_min)/(f_max-f_min), logit(p)=logit(p0)+beta*(a+r*g) exactly
  in real arithmetic. Default latent interaction orders are 2, 3 and 4; hence
  logit fitness has no uniform-product Fourier interactions above order four.
  The sigmoid can generate higher-order magnitude epistasis on the fitness scale,
  but preserves exact rankings, mutation signs, peaks and uphill paths. Statements
  about invariant walks require rank-only rules and zero tolerance; observation
  noise, finite resolution and numerical saturation can change them.
- This is not standard RMF: the residual is correlated across genotypes and the
  additive matrix does not impose one common Hamming-distance slope. Increasing r
  does not remove the interaction field's correlations. At r=0 the model is an
  additive latent-trait global-epistasis model. The current construction is not
  standard Fisher geometric adaptation to a multidimensional phenotypic optimum,
  and the occupancy equation alone does not identify physical free energies.
- Analytical consequence under the implemented uniform-q20 nonself mutation
  operator: for genotype-mean-centered pure-order g_k, P*g_k equals
  (1-20*k/(19*L))*g_k. With order mutation-energy shares w_k summing to one,
  Var(g)=19*L/40*sum_k(w_k/k), and the correlation averaged over uniform directed
  one-mutation edges is 1-1/(2*Var(g)). Default weights at L>=4 give
  sum_k(w_k/k)=0.375 and correlation 1-2.80701754386/L. These are intrinsic
  uniform-background identities, not empirical or reference-shell claims.
- Direct arithmetic checks at lengths 4, 9 and 100, seed 7, beta=0.1 and r=0.5:
  implemented versus derived Var(g) differs by at most 2.14e-14; inverse-logit
  recovery on 64 uniform sequences per length (NumPy generator seed 19) differs
  by at most 4.45e-16. No new empirical experiment, RES or CLAIM was created.
- Literature checked: [Neidhart et al. (2014), RMF](https://pmc.ncbi.nlm.nih.gov/articles/PMC4196622/);
  [Brookes et al. (2022), GNK and spectra](https://pmc.ncbi.nlm.nih.gov/articles/PMC8740588/);
  [Otwinowski et al. (2018), global epistasis](https://pmc.ncbi.nlm.nih.gov/articles/PMC6094095/);
  [Neidhart et al. (2013), amplitude spectra](https://arxiv.org/abs/1301.1923).
  Labproj status, strict refresh and strict check each timed out at 12 seconds;
  no successful whole-record validation is asserted.

## 2026-10-08 current model consolidated — Adapt-Env 0.4.0 / CLAIM-004

- Scientist authorized removal of the retired methods, clear model/API documentation,
  and a claim about useful TrpB behaviour. Adapt-Env dev is committed locally at
  `37d23c64b15a26cafa41232940feabc218505cd5`; the external working tree is clean. No push or deployment was performed.
- Removed 65 tracked legacy files: the entire `_legacy` package, compatibility
  imports, dependent old simulation/experiment packages and their tests/docs.
  Frozen scientific RUN sources and ignored historical result data are preserved.
  The package and project model namespace now advertise only the current model.
- Wrote model equations, normalization and scaling assumptions, Python/HTTP API,
  configuration/preset documentation, migration instructions and TrpB evidence.
  Added `trpb_run175_r05_shared_stop_intrinsic`, preserving RUN-175's exact intrinsic
  readout. The public Gaussian observation option remains distinct from the separate
  replicate-derived shared-stop research sensitivity.
- CLAIM-004 supports conditional approximate joint TrpB aggregate/search resemblance,
  directly citing ANA-069/RUN-175 and preserving the negative observations from
  ANA-067/068. It includes equal-cube and pooled WHT, seed ranges, observation
  dependence and limits on cross-assay/length transfer. No new experiment or RES.
- Verification: 33 maintained package tests (including the built wheel), 40 current
  project tests, four executable Python documentation examples and 36 local doc links.
  Replay covers 437,184 values / 576 configuration-seed pairs. Fitness and fingerprints
  match the frozen source exactly when evaluated locally; stored worker values agree
  within 1.07e-14. Wheel/sdist omit archived packages and compatibility modules.
- Verification, source manifest and removal inventory are in `claims/CLAIM-004_evidence`.
  Broad labproj checks still time out on historical records. Older calibration/SMC
  integration tests require their earlier source environment, as documented.

## 2026-10-08 EXP-065 interaction-strength follow-up — RUN-175

- Scientist authorized the paired strength sweep motivated by ANA-068, queue
  execution, active monitoring through completion and a figure-led notebook ANA.
  Ratios are 0, 0.5, 1, 2, 4, 8; each GB1/TrpB setting uses 16 original and eight
  fresh seeds. Only the original global readout objective is fit; WHT, DR and
  adaptive walks remain unfitted diagnostics. Observation assumptions stay fixed.
- RUN-174 is the dry preview. RUN-175 attempt 1 (Slurm 2321) failed before science
  because the worker lacked g++. The failed attempt is preserved. A project-local
  pinned Zig C++ compiler passed all 20 independent graph checks; shared Python
  environments and model architecture are unchanged.
- RUN-175 attempt 2 completed as Slurm 2322 on dutch in 558.29 seconds of
  scientific computation, with runtime exit code 0 and empty stderr. Slurm accounting
  is unavailable; terminal evidence and all output hashes are preserved.
- The 24 frozen input files are SHA-256-verified on the worker and DVC-tracked as
  `data/processed/interaction-strength-exp065-inputs.dvc`; all 26 remote objects
  are MD5-verified. Local smoke reproduces 30 ANA-068 graph rows with zero
  discrepancy. No scientific result is inferred from these implementation checks.

- ANA-069 contains eight executed notebook cells and five PNG/PDF figures. All 390
  ANA-068 endpoint graph rows reproduce exactly; all 24 objectives and 1,732 WHT
  summary rows are independently verified. The sweep produces 864 synthetic maps,
  2,598 graph summaries and 288 exact latent-coordinate summaries.
- Fresh-seed GB1 r=0.5 improves joint resemblance: uniform walk length 8.20 versus
  9.74 empirical; WHT epistatic fraction 0.546 versus 0.555; DR median -0.637 versus
  -0.690. Local SD remains 2.64 times empirical. TrpB noiseless r=1 matches walk
  length more closely but top-1% success remains 99.5% versus 65.0% empirical.
- TrpB shared-stop r=0.5 gives 7.46 substitutions, 71.4% top-1% success and WHT
  fraction 0.630 versus empirical 8.59, 65.0%, 0.688. This is a conditional coarse
  benchmark regime; shared noise also generates much apparent epistasis at r=0.
  Fresh-seed top-1% success has a 52–85% 10–90% range. Earlier r=4 shared-stop
  agreement is seed-sensitive (70.1% reused versus 87.5% fresh median success).
- RUN data, source and inputs are DVC-tracked: RUN-175, 186 files / 151,681,986
  bytes, directory MD5 b33c8fdb785160b4df3f30d0a26c4711.dir. All 184 unique remote
  objects on dutch verified. No public model change, new RES/CLAIM or scientific
  success gate. These are four-site post-fit comparisons, not arbitrary-length
  or independently validated biological realism.

## 2026-10-08 GB1/TrpB whole-library dynamics — ANA-068

- ANA-068 reuses EXP-066 / RUN-172 parameters and 16 frozen seeds per condition
  without refitting. Nine executed notebook cells produce six PNG/PDF figures
  comparing adaptive walks, accessibility, neutral components, observation
  sensitivity and missing-genotype sensitivity. No new RUN, RES, CLAIM or
  scientific gate is introduced.
- Empirical graphs retain 149,361 GB1 and 159,129 TrpB genotypes. Synthetic maps
  use identical masks. Computation covers 130 landscapes, 390 resolution-specific
  graph summaries and 532,480 sampled walks; exact dynamic programming gives walk
  expectations. Graph calculations take 139 seconds locally. Twenty independent
  small-graph cases pass; all sampled walks terminate before the cap.
- High-fitness regions are broadly accessible, but routes and walk outcomes
  differ. Uniform-uphill expected lengths are 9.74 empirical versus 5.27 model
  median for GB1 and 8.59 versus 6.37 for TrpB (r=4, no added noise). TrpB's
  upper-1% walk success is 65% empirical, 96% noiseless and 70% shared-stop
  observed; that observation branch further shortens walks to 4.87 substitutions.
- GB1's exact neutral component of 29,477 genotypes is entirely at the zero-fitness
  floor. Nonzero-tolerance components are strongly resolution/observation dependent
  and can span many edge tolerances; they do not establish shared useful function.
  The whole-library analysis includes earlier calibration nodes and is descriptive
  post-fit evidence, not independent validation or arbitrary-length evidence.
- Derived data: `data/processed/production_landscapes/ANA-068.dvc`, directory MD5
  `4784c138999f69f1b522198f43220458.dir`, 403 files / 173,163,035 bytes. Local
  cache objects and all output hashes are verified. Remote backup is pending:
  the lab gateway reports no route to `dutch`. Broad labproj checks time out on
  historical cloud files; targeted source, mathematical, notebook and artifact
  verification is recorded in the ANA.

## 2026-10-08 RUN-172 production analysis — ANA-067

- ANA-067 provides the executed figure-led analysis of EXP-066 / RUN-172 at
  `analyses/2026-10-08_ANA-067_run-172-production-landscape-aggregate-realism/`.
  Ten notebook code cells produce six inline figures, PNG/PDF exports and exact
  estimates with 1,500-draw cube/seed bootstrap intervals. No scientific gates,
  new RES or CLAIM, or empirical interaction reconstruction are introduced.
- Coverage is 63 evaluated assays (13 function/binding, 50 stability): 52 two-site,
  two three-site and nine four-site panels. Only 11 inform order >=3 WHT.
  The derived inventory restores three omitted final exclusion rows without
  modifying the frozen RUN; all 86 inventory entries are accounted for.
- The fixed r=4 sigmoid offers conditional useful regimes, not uniformly close
  recovery across assay classes. Median synthetic/empirical mean cube SD is
  0.696 for function/binding and 0.270 for stability. Stability median empirical
  versus synthetic sign epistasis is 0.394 versus 0.078; mean DR is -0.398 versus
  -0.004. PHOT/HIS7 also have substantial WHT/sign-structure shortfalls.
- TrpB r=4 shared-stop observed predictions resemble empirical equal-cube
  epistasis (0.689 versus 0.688), high-order fraction (0.293 versus 0.316),
  mean DR (-0.896 versus -0.880), SD (0.0358 versus 0.0379) and sign epistasis
  (0.638 versus 0.651). Pooled high-order agreement is weaker (0.202 versus
  0.086). Observation assumptions materially contribute: intrinsic sign
  epistasis is 0.212 and DR is -0.451. This is not intrinsic mechanism recovery.
- DR means can be unstable: GB1's empirical mean is +4.023 while its median is
  -0.480. Resolution sensitivity, slope availability and medians are retained.
  Calibration uses finite candidate searches (76/195 report convergence);
  testing r=0 and r=4 does not exhaust the public model's parameter space.
- The complete local data and exact execution source are now DVC-tracked at
  `data/processed/production_landscapes/RUN-172.dvc`: directory MD5
  `521d0e1b381ead592b0f47504e0f13ba.dir`, 1,842 files / 227,667,682 bytes.
  All 1,840 unique remote DVC objects were MD5-verified after SSH backup.
- Independent audit verifies 1,819 output hashes, 12 source hashes, 179,310
  statistic rows, 23,016 empirical nodes, 230,400 cached intrinsic values,
  public API replay on every assay and all 195 selected objective minima.
  Sixty-five tiny-power WHT fraction entries show bounded roundoff sensitivity
  only at the smallest resolutions; primary-resolution conclusions are unchanged.
  Broad labproj checks timed out on historical cloud-backed records; targeted
  data, source, numerical, notebook and artifact checks are recorded in ANA-067.

## 2026-10-08 RUN-172 production execution completed

- RUN-172 attempt 3 (Slurm 2317) completed in 4,416.3 seconds (73 min 36 sec),
  exit 0, empty stderr. It evaluated 63 assays and produced 195 fitted settings;
  there were no assay execution errors. One geometrically eligible CR9114 FluB
  panel had constant calibration values, and 15 candidates lacked sufficient
  sampled disjoint cubes.
- The 1,819-file output manifest and completed bundle were automatically published
  to `/tank/data/labproj/projects/project-adapt-env/data/processed/production_landscapes/RUN-172`.
  Published and runtime manifests match. Local collection/DVC tracking and
  scientific analysis have not yet been performed for this production RUN.
- Inventory bookkeeping caveat: `suitability.csv` has 83 of the 86 inventory rows.
  The runner did not flush the final three explicitly excluded FLIP entries
  (duplicate GB1, mixed-indel AAV and cross-protein Meltome). All 79 candidate
  decisions and all 63 fitted assays are present. Preserve the frozen bundle and
  include these exclusions in the later analysis. No scientific conclusion is
  inferred from successful execution alone.

## 2026-10-08 public release and production EXP-066

- Adapt-Env 0.3.0 is pushed to GitHub `dev` at
  `dd0bb7a52d7d92af48bbfba115331742b9a0c2eb`. The public landscape aliases,
  presets, service and API now use the normalized rank-8 thermodynamic model.
  Historical methods live under private `_legacy` with deprecation warnings.
  Default evaluation is noiseless; calibrated TrpB observation assumptions are
  explicit experimental branches, not implicit public-model noise.
- Fixed feasible short-system balanced-support construction; established successful
  seeds retain their supports. External package suite: 87 passed, one skipped;
  public/observation checks: 14 passed; root mathematical tests: 41 passed.
  Built wheel installed and imported separately. The external checkout is clean.
- EXP-066 is the authorized systematic production comparison before an online ML
  pilot. It audits all 86 existing inventory entries, structurally samples all
  79 nonduplicate fixed-length protein candidates (including 64 stability panels),
  and records coverage exclusions and all unfavorable comparisons. No scientific
  success gates or actual interaction reconstruction are attempted.
- Calibration uses only fitness quantiles, near-floor mass and mean local SD;
  confirmation withholds WHT, DR, DFE and topology on disjoint empirical cubes
  and fresh generator seeds. Dimension-matched geometry and empirical readout
  context remain explicit. This does not establish all-length statistical
  invariance or uniquely identified thermodynamics.
- The run publishes complete, hash-verified figure-ready bundles to the durable
  cluster counterpart of `data/processed/production_landscapes/<RUN-ID>`;
  collection and DVC tracking use the same local namespace. A source-only snapshot
  and explicit SHA-256 manifests avoid unrelated cloud-backed history.
- RUN-172 attempt 3 is running as Slurm 2317 on `dutch` (four CPUs, 16 GB,
  two-hour allocation). RUN-171 and RUN-173 were dry previews. Attempt 1 failed
  before science in editable overlay setup; attempt 2 was stopped after a TrpB
  compact-mirror float32 precision check. Both corrections and attempt history
  are recorded; the scientific design and data are unchanged.
- All 86 input parquet hashes were verified. Geometry preflight found 64 usable
  assays (10 four-site, two three-site, 52 two-site), with 15 additional candidates
  lacking enough sampled disjoint cubes. Five production tests and the five-assay
  worker smoke passed. The full TrpB panel passed the corrected rounding check.
- Submission verification confirms all 12 frozen source hashes and the exact
  worker versions, with scientific fitting in progress and empty stderr. Outputs
  will publish automatically after completion. No production conclusion is yet
  asserted. Broad labproj checks timed out on historical cloud-backed records;
  targeted source, numerical, package and input checks are recorded.

## 2026-10-08 RUN-170 analysis — ANA-066

- RUN-170 completed as Slurm 2314 in 60.86 seconds, exit 0, empty stderr.
  ANA-066 supplies eight executed notebook code cells and five inline figures with
  PNG/PDF exports. Exact analysis lives at
  `analyses/2026-10-08_ANA-066_run-170-observation-aware-global-calibration/`.
- Six observation-aware global fits used only measured fitness quantiles,
  near-floor mass and average cube SD. The r=4 sigmoid calibration losses are
  0.0475 independent-error / 0.0080 shared-stop, versus 1.0157 / 0.4871 linear.
  Sigmoid searches exhausted their finite budget; some scales approach the upper
  bound. These are useful candidates, not proven optima or identified energies.
- Fresh confirmation gives empirical / r=4 independent / r=4 shared-stop
  equal-cube high-order fractions 0.267 / 0.200 / 0.283 and pooled fractions
  0.121 / 0.118 / 0.139. All-epistasis fractions are 0.660 / 0.591 / 0.683
  equal-cube and 0.496 / 0.483 / 0.525 pooled. Absolute power varies by seed.
- Empirical all-mutation slope is -0.599 (interval [-0.949,-0.084]); beneficial
  slope -0.335 ([-1.057,+0.695]). The prior approximately -1.03 beneficial value
  is not a precise universal target. Observed r=4 slopes are -0.667 / -0.890
  all-mutation and -0.474 / -0.800 beneficial across the two noise assumptions.
- Intrinsic r=4 all/beneficial slopes are -0.476/-0.124 under independent-error
  calibration and -0.062/+1.018 under shared-stop calibration. Observed
  diminishing returns must not be conflated with intrinsic biological structure.
  Additive shared-stop observed high-order fraction is 0.274 but intrinsic is
  0.047, compared with 0.218 for r=4 shared-stop.
- Activity mismatch persists: empirical active vertices 5.76%, no-active cubes
  59.4%; r=4 shared-stop gives 14.4%, 41.7%. Additive shared-stop gives 6.52%,
  60.3% but much less intrinsic high-order structure. No setting dominates all
  descriptors. No scientific success gate or equivalence claim is imposed.
- Recommendation: keep these frozen regimes and test the actual online ML task
  in one bounded pilot (random, local and one ML policy; fixed budgets), scoring
  intrinsic fitness and noisy observations separately. Independent-error r=4 is
  the clearer intrinsic-fitness starting point; shared-stop is an assay
  sensitivity. Add latent complexity only if a residual distorts relevant ML
  behavior. This is advice, not another submitted RUN. Four-site evidence does
  not establish global statistical invariance with length or new-assay transfer.
- Notebook audit checks 103 output hashes, 212 remote hashes, seven source
  hashes, 2,048 fresh empirical genotypes, 3,712 donor exclusions, 114,432
  statistic rows, 1,179,648 replayed observations and all six selected losses.
  Intervals resample 64 empirical cubes or 16 synthetic seeds, not 32 repeated
  observation draws as independent landscapes. All five figures were inspected.
- The exact bundle has 215 files / 166,435,332 bytes, DVC directory MD5
  `1f928b0b95c7ec39570765c6a051f8c9.dir`; all 208 unique objects were backed up
  and MD5-verified. SSH collection replaced the unavailable client share. Broad
  strict checks timed out on historical cloud-backed records; targeted notebook,
  source, data and artifact checks passed. No RES or CLAIM was created.

## 2026-10-08 RUN-170 completed; ANA-066 in progress

- RUN-170 attempt 2 completed as Slurm job 2314 in 60.86 seconds, exit code 0,
  empty stderr. Attempt 1 failed in source packaging before scheduler submission;
  narrowing the source snapshot resolved this. RUN-169 was the dry preview.
- Six global parameter fits were frozen before confirmation. The new panels have
  64 calibration and 64 confirmation cubes, disjoint from each other and all old
  exposed nodes; 3,712 genotypes were excluded from observation donors.
- Source/data checks, 12 numerical tests and an old-panel smoke test passed.
  `labproj run reconcile RUN-170` confirmed completion. Collection uses SSH because
  the configured client NFS/SMB share is unavailable.
- ANA-066 is auditing the complete outputs and generating five figures through an
  executable notebook. No additional experiment, result or claim is authorized.

## 2026-10-08 observation-aware calibration authorized

- The scientist authorized the ANA-065 proposed follow-up as one new RUN within
  EXP-065, including submission, active monitoring and a figure-led executable ANA.
- The prospective protocol is `observation_calibration.md` and YAML under EXP-065.
  Six global fits cross r=4 linear/sigmoid and r=0 sigmoid with two frozen noise
  assumptions. Calibration uses marginal measured fitness, near-floor mass and
  mean within-cube SD; WHT allocation, slopes and conditional distributions are
  withheld. Fresh disjoint panels precede donor estimation.
- Focus is useful global resemblance for online ML testing, with explicit residual
  differences and no scientific success gates. No specific interactions are fitted.
- Broad labproj status/refresh/check each timed out after eight seconds on existing
  records. The network target doctor passed with the known staged-root visibility
  warning; targeted validation and exact source/data checks are used for this RUN.

## 2026-10-08 scientific interpretation after ANA-065

- Asked for interpretation and a next experiment, prioritize an observation-aware
  calibration of the existing model before changing the latent interaction
  spectrum. This is advice and a proposed experimental direction, not an
  authorized new EXP/RUN, implementation or queue submission.
- The equal-cube WHT fraction remains 0.2668 ordinary / 0.2671 cross-replicate,
  whereas pooled power gives 0.1025 / 0.0972. Three cubes contain 74.7% of ordinary
  total power. Alongside 12/24 no-active empirical cubes versus 8/64 intrinsic
  sigmoid cubes at r=4, this motivates a hypothesis about the distribution of
  local activity and variation. It does not establish that operating-regime
  mismatch explains the residual, or that shared error is absent.
- RUN-168 deliberately added observation noise after freezing affine scales
  fitted in RUN-166 to noisy measured mean within-cube SD. Some spread inflation
  is therefore expected. Its SD increase does not demonstrate that a model
  calibrated with the observation process included must have excessive spread.
  The old RUN remains a valid frozen-parameter diagnostic; no result is revised.
- Proposed question: can the current global readout and observation process,
  calibrated to the distribution of measured fitness and overall within-cube
  variation, predict held-out WHT power allocation and fixed-mutation slopes?
  Keep rank, supports and order weights fixed; use the r=4 linear/sigmoid pair
  and an additive-latent sigmoid control. Fit only global readout/affine settings,
  never per-cube offsets or empirical interaction identities. Hold WHT order
  fractions and diminishing-returns slopes out of calibration.
- Retain independent-difference and shared-stop observation assumptions as
  separate analyses; estimate their parameters from replicate/control training
  data rather than fitting them to erase aggregate discrepancies. Sample fresh
  genotype-disjoint calibration/confirmation cubes before noise calibration,
  exclude their genotypes from noise donors, use fresh synthetic seeds, and
  preserve the original complete-cube sampling rule without activity filtering.
  Prior assay-wide exposure still limits claims of wholly independent replication.
- Compare absolute WHT order powers, equal-cube and pooled fractions, slope
  numerators/denominators and availability; also predict the distributions of
  cube variance and activity composition. Calibration can use marginal fitness
  shape and average within-cube variation while leaving these conditional
  distributions as checks. Do not use the unstable stop-subtracted ratios as a
  uniquely identified biological target or introduce scientific success gates.
- If matching observed marginals also predicts the withheld aggregate summaries,
  that supports an operating-regime/measurement explanation within this family.
  Residual conditional spectral or slope discrepancies would support changing
  latent effect heterogeneity or order structure. Conclusions that reverse
  between noise assumptions remain measurement-limited rather than identifying
  a unique biological mechanism.
- Literature grounding: Otwinowski, McCandlish & Plotkin (2018),
  https://doi.org/10.1073/pnas.1804015115 , separates global nonlinearity and
  measurement/error components; Reddy & Desai (2021),
  https://doi.org/10.7554/eLife.64740 , shows that diminishing returns can emerge
  from pervasive microscopic interactions. These motivate the controls and
  interpretation, not a claim that their particular fitted models are adopted
  here or that this proposal identifies actual molecular free energies.

## 2026-10-07 RUN-168 analysis — ANA-065

- RUN-168 (Slurm 2313) completed with exit code 0, empty stderr and 103.16 seconds
  of production execution. `labproj run reconcile RUN-168` confirmed completion.
  ANA-065 provides its figure-led report and canonical executed notebook at
  `analyses/2026-10-07_ANA-065_run-168-replicate-aware-statistics-and-observation-models/`.
  All eight code cells execute without errors; five inline figures are also
  exported as PNG/PDF. No RES, CLAIM or further experiment was created.
- Independent audit verifies 162 remote file hashes, 73 output hashes, eight
  submission source hashes, 14 input/source hashes, all 768 empirical measurements,
  193,170 activity flags and exclusion of 1,664 exposed genotypes from donors.
  A separate butterfly WHT/covariance implementation verifies 394,752 synthetic
  statistic rows and replays 4,210,688 observation values from 131,584 saved
  replicate-pair cubes. It also checks 144 empirical statistic rows, 1,728 mutation
  rows, 1,440 synthetic summary rows and all 108 original bootstrap rows.
- Only 32/384 confirmation vertices exceed both published activity criteria.
  Twelve cubes contain no active vertices, twelve contain some, and none contain
  all active vertices. Training/holdout control covariance agrees descriptively;
  donor holdout/training differential-variance ratios span 0.784–1.060. This tests
  repeatability of control summaries, not identification of pure assay error.
- Mean per-cube high-order WHT fraction is 0.2668 ordinary and 0.2671 cross-replicate:
  independent replicate disagreement alone does not explain this discrepancy.
  The ratio of pooled powers is instead 0.1025 ordinary and 0.0972 cross; three
  cubes contain 74.7% of ordinary total power. These are different estimands;
  pooled agreement cannot retrospectively establish adequate equal-cube recovery.
- Stop-control adjustment gives pooled high-order ratio 0.0632, but every no-active
  cube has negative adjusted total power. Only 12/24 individual spectra remain
  resolved; four fractions lie outside [0,1]. Their mean 0.6781 is unstable and
  cannot be interpreted as an increase in biological epistasis. Shared inputs,
  normalization, control biology and covariance transfer remain confounded.
- Ordinary/cross/stop-adjusted beneficial slopes are -1.029/-1.882/-0.958, with
  47/47/16 eligible mutations. The conditional stop-adjusted interval is
  [-2.415,+0.538]. The corresponding all-mutation slopes are -0.923/-1.225/-0.617.
  Availability and ratios change with measurement assumptions; no unique
  error-free diminishing-returns rate is identified.
- At frozen r=4, shared-stop observation changes sigmoid mean per-cube WHT fraction
  from 0.1442 to 0.1885 and beneficial slope from -0.500 to -0.692, while SD rises
  from 0.05990 to 0.06867 versus measured 0.05644. Linear changes similarly.
  Observation noise explains part of some descriptor differences while increasing
  spread. Activity composition remains different; 128 observation draws reuse the
  same 64 cubes per arm and are not independent landscapes. No refitting or
  scientific acceptance gates were introduced.
- Collected exact outputs, logs and submitted source over SSH because the client
  share was unavailable. DVC tracks 165 files / 224,508,438 bytes with directory
  MD5 `b71af13a0354c4a7a8a65cb77e8f019a.dir`; all 159 unique cache objects were
  backed up and MD5-verified in the existing project store. The original TrpB
  replicate bundle remains the second explicit data dependency.
- `labproj analysis sync ANA-065 --write` generated the text companion and indexed
  20 artifacts. Notebook execution, all five inline images, ten figure exports,
  report links and source hashes were checked. Broad status/refresh/check and the
  data wrapper timed out on historical cloud-backed records; direct DVC tracking
  succeeded. No repository-wide strict-check pass is claimed.
- Completion commands saved ANA-065 and EXP-065 as completed before follow-up
  checks timed out. Targeted metadata and registry reads verify that RUN-168,
  ANA-065 and EXP-065 are all completed.

## 2026-10-07 replicate-aware follow-up — RUN-168 submitted

- The scientist authorized a new RUN under EXP-065 for the direction proposed
  after ANA-064. `replicate_observation.md` and its YAML freeze the original
  48 empirical cubes and eight selected RUN-166 arms (64 confirmation cubes each).
  The task is diagnostic reanalysis and observation modelling, not a new generator
  fit or independent empirical confirmation; no scientific success gates apply.
- Recovered exact original TrpB summary tables and processing code from
  CaltechDATA DOI 10.22002/h5rah-5z170 using ZIP byte ranges, member CRC and SHA-256.
  All 768 RUN-166 measurements match the published replicate mean within 1e-7.
  The derived activity definition matches all 193,170 published flags. Both
  replicates share averaged input counts; stop controls also covary strongly.
  Cross-replicate estimators therefore cannot be labelled error-free.
- New runner exports ordinary/cross/stop-adjusted WHT powers and fixed-mutation
  slope components, activity-stratified summaries, conditional bootstraps, frozen
  synthetic observation arrays and three PNG/PDF figures. Noise estimation excludes
  1,664 exposed genotypes and reserves donor/control holdouts. Signed corrected
  powers and unresolved ratios are retained. Stop covariance is an explicit
  sensitivity assumption that can include control biology and normalization.
- Eight numerical controls pass, including heteroscedastic shared-noise correction
  and donor holdout isolation. The separate-seed, two-draw SMOKE_ONLY pipeline
  completed in 19 seconds. Frozen r=4 intrinsic summaries reproduce RUN-166.
- New 15-file input bundle is DVC-tracked as
  `data/processed/trpb-replicate-observation-exp065.dvc`
  (`ad4e31a4a5704520d0a4be2cab9fa281.dir`, 90,025,877 bytes).
  Broad labproj status/refresh/check and data-track commands timed out while
  accessing historical cloud-backed records; direct DVC tracking succeeded.
  This does not establish a clean global strict check. Remote verification passed
  for all 14 declared prepared/source hashes and all 16 DVC objects; production
  submission is confirmed below. No ANA/RES/CLAIM created.
- RUN-167 is the nonexecuting preview. RUN-168 attempt 1 was accepted as Slurm
  job 2313 and verified RUNNING on dutch, requesting 2 CPUs, 8 GB and 2 hours.
  Startup reports the production 128-draw protocol with empty stderr; eight staged
  source/config/protocol/test/pointer hashes match the local files. The submission
  audit is `runs/RUN-168.submission-verification.json` within EXP-065. The temporary
  source-staging filter is archived with the RUN and the original project filter
  was restored. No scientific conclusion is drawn from the smoke outputs.

## 2026-10-07 scientific interpretation after ANA-064

- In response to the scientist's request for interpretation and a next step,
  prioritize a replicate-aware check of the empirical target before increasing
  model complexity. This is a proposed direction, not a new EXP/RUN submission.
  The small primary-resolution sigmoid advantage does not identify a physical
  thermodynamic mechanism; local observable activity through length 512 and
  recovery of empirical aggregate structure remain separate findings.
- Rechecking Johnston et al. (2024), Methods and Fig. 1, establishes that the
  original TrpB study used two replicate flasks, averaged their fitness estimates,
  and explicitly controlled noisy low-fitness measurements in its epistasis
  analyses. Its replicate-specific activity cutoffs were approximately 0.054.
  The prepared table used by RUN-166 lacks replicate columns; this is not evidence
  that the original measurements lack replicates. Public data/code are deposited
  at https://data.caltech.edu/records/h5rah-5z170 (archive contents not yet audited).
  Primary source: https://www.pnas.org/doi/10.1073/pnas.2400439121 .
- An exploratory read-only calculation on saved RUN-166 confirmation values finds
  that 348/384 (90.625%) are below 0.05 on the imported raw fitness scale, with
  median 0.011294. This is not the paper's two-replicate activity classification;
  it identifies a reason to examine assay-floor sensitivity. The original paper
  finds epistasis exceeding its noise null, so no conclusion that TrpB epistasis
  is merely measurement noise follows. Shared initial sampling/normalization
  can correlate replicate errors and must be assessed before cross-replicate
  estimators are described as error-free.
- A second exploratory calculation uses the 64 saved confirmation GNK cubes
  alone: the linear pure-G (r approaching infinity) mean high-order WHT fraction
  is 0.134164, compared with 0.112693 at r=4 and 0.266817 in the empirical panel.
  This is a fixed-panel limiting diagnostic, not an independent confirmation or
  a bound on all nonlinear/readout/order-weight configurations. It suggests
  limited benefit from increasing only the interaction/additive ratio; persistent
  reproducible discrepancy would instead motivate changes in aggregate order
  weights and the sampling/geometry of interacting sites.
- The proposed minimal follow-up would retain frozen model parameters and compare
  raw versus replicate-aware empirical WHT power and fixed-mutation slopes, then
  compare the same synthetic landscapes through an observation process estimated
  from replicate/control measurements. Retain the full sampled panel alongside
  explicit activity strata. Do not infer specific empirical interaction identities
  or fit measurement noise simply to erase a descriptor discrepancy.
- Reddy & Desai (2021), https://doi.org/10.7554/eLife.64740 , provides the relevant
  theoretical caution: diminishing returns can emerge from widespread latent
  interactions without a saturating readout. The new TrpB comparison therefore
  supports useful aggregate modelling more directly than a unique physical
  explanation. No previously completed ANA or RES is retrospectively revised.

## 2026-10-07 RUN-166 analysis — ANA-064

- Created and completed `ANA-064` for the three-part EXP-065 follow-up. Its
  canonical `notebooks/analysis.ipynb` executes all eight code cells without
  errors and displays five figures, also exported as PNG/PDF. The report is
  `analyses/2026-10-07_ANA-064_run-166-paired-readouts-and-embedded-neighbourhoods/analysis.md`;
  its README documents restoration and reproduction. No RES or CLAIM was created.
- Independently checked all 20 output/log hashes, nine submitted-source hashes,
  both input hashes, 768 empirical measurements, 40 calibration losses, and 24
  confirmation distances. Replayed 3,312 observed cubes from saved latent fields,
  recalculated 9,936 statistic rows and verified 105,984 oriented mutation edges.
  Genotype exclusion, seed separation, exact pairing and frozen transfer
  parameters are audited. This is saved-latent replay, not cross-platform
  regeneration of GNK coefficients.
- At the calibration-selected ratio r=4, primary held-out distance is 0.2784
  linear versus 0.2617 sigmoid. The paired difference is -0.0166, with a
  conditional 95% bootstrap interval [-0.0647, +0.0734]. Interaction strength
  dominates the readout comparison within this grid; the selected ratio is its
  upper boundary. Mean high-order WHT fractions remain 0.1127 linear / 0.1442
  sigmoid versus 0.2668 measured. Mean beneficial slopes are -0.577 / -0.500
  versus -1.029 measured. Observable SD is much closer after amplitude calibration,
  but absolute distances cannot be compared directly with RUN-164's embedding.
- Operating occupancy changes slopes, WHT fractions, absolute power and neutrality
  together. These are fixed-parameter interventions, not extra fitted models.
  The primary readout advantage is small and uncertain; resolution sensitivity
  and unavailable spectra/slopes are shown explicitly. Empirical measurement
  variance and shared-error slope bias remain uncorrected because no independent
  replicate errors are available.
- Embedded libraries preserve observable reference-neighbourhood variation
  through length 512. At reference rank 0.5 all 24/24 reference-anchored cubes
  per site policy are resolved, versus only 6/24 or 7/24 for uniform backgrounds.
  Uniformly sampled sites lose latent high-order structure at large length;
  a shared-support control retains some, without reproducing the four-site
  TrpB statistics. Rank-1 uniform backgrounds are unresolved at every length
  and numerically flat at 128/512. These results delimit sampling and reference
  dependence; they do not establish arbitrary-length statistical invariance.
- Collected outputs, logs and exact submitted source via SSH because the client
  share was unavailable. DVC tracks the 83-file bundle (54,423,275 bytes), directory
  MD5 `c128be3fda98517848f7d5533b5d2f52.dir`; all 82 unique cache objects including
  the directory manifest were backed up and hash-verified in the existing
  project cluster store. Temporary DVC bypassed the cloud-backed project
  environment without changing its configuration.
- `labproj analysis sync` generated the text companion and indexed 20 artifacts.
  ANA-064 and EXP-065 metadata/registry statuses are completed. Completion commands
  saved those statuses before their global follow-up checks timed out; initial
  status/refresh/check attempts had the same historical-file problem. Notebook
  execution, all five inline images, ten figure exports, links and targeted record
  statuses are verified; no repository-wide strict-check pass is claimed.

## 2026-10-07 EXP-065 paired-readout follow-up — RUN-166

- Completion check: RUN-166 has exited successfully with code 0. All 17 declared
  outputs, including the three PNG/PDF figure pairs, exist on the cluster. They
  were subsequently collected and analysed in ANA-064 (see above).
- Submitted production `RUN-166` on attempt 1 as Slurm job `2312`, CPU partition,
  2 CPUs, 8 GB RAM and a two-hour limit. Scheduler startup is confirmed RUNNING
  on `dutch`; the job has reached the embedded length-128 panels with empty
  stderr. Staged runner/config SHA-256 hashes match the local submitted files.
  The source-only staging filter is retained with the RUN, and the original
  project staging filter was restored. No scientific conclusion is inferred
  from progress logs; the original RUN-164/ANA-063 remain completed.
- The scientist authorized one further production RUN under EXP-065 containing
  all three proposed comparisons: paired linear/sigmoid readouts of the same
  latent field, conditional WHT/mutation statistics and operating-position
  diagnostics, and four-site libraries embedded at lengths 32, 128 and 512.
  `paired_readout.md` and `paired_readout.yaml` in the EXP directory define the
  prospective design. The primary assay is TrpB; no new EXP, ANA or RES is created.
- Both readouts receive calibration-only affine amplitude adjustment. Four
  latent ratios yield four paired contrasts, with fixed rank-8 order weights.
  New empirical cubes exclude all 768 RUN-164 TrpB genotypes and 128 genotypes
  exposed by the engineering smoke run. Confirmation and embedding seeds are
  separate from calibration and engineering seeds. No numerical scientific
  acceptance gates are added.
- The embedded comparison freezes the calibration-selected sigmoid and its
  linear partner at the same ratio. It crosses length, reference rank (0.5/1),
  anchor Hamming distance (0/4/16/uniform), and site selection (uniform/shared
  order-four support), with library/seed identifiers retained. The shared-support
  arm is a sampling control, not an inferred TrpB interaction map.
- Verification: 27 focused tests pass, and an end-to-end measured-data smoke
  generates all three parts and three PNG/PDF figures using reserved engineering
  seeds. The dry-run submission is RUN-165. Slurm target checks and numerical/
  plotting imports on the x86 worker succeed; the ARM submit host cannot execute
  that worker environment directly. Repository-wide status/refresh/check still
  time out on historical cloud-backed files; no global strict-check pass is claimed.

## 2026-10-07 RUN-164 analysis — ANA-063

- Created and completed `ANA-063` for `EXP-065/RUN-164`. Its canonical
  `notebooks/analysis.ipynb` executed all seven code cells without errors and
  produces five scientific figures in PNG/PDF plus audited comparison tables.
  `analyses/2026-10-07_ANA-063_run-164-aggregate-thermodynamic-gnk-confirmation/analysis.md`
  contains the interpretation; its README documents reproduction.
- Collected all eight outputs, three logs and the submitted source snapshot over
  SSH because the shared filesystem was not mounted. Verified output/log/source
  hashes and the three input-table hashes. Independently reproduced measured
  WHT, beneficial-mutation slopes, neutral edges and strict peaks for all 144
  cubes, all nine joint distances, and calibration/confirmation separation.
  The 60-file bundle is DVC-tracked; all 60 unique cache objects are backed up
  and hash-verified on the project's existing cluster store. A temporary DVC
  environment and SSH transfer avoided the unavailable project environment and
  local network mount, without changing permanent environment/remote settings.
- Composite held-out distances are 0.9565 (GB1), 0.5750 (TrpB), and 0.3109
  (CR9114-H1). It has the lowest point distance for GB1; bare GNK is lowest for
  the other two assays. The composite's small GB1 advantage over additive and
  CR9114-H1 difference from GNK are uncertain in a paired conditional bootstrap.
  The general predicted advantage over both baselines is not supported.
- Recovery remains partial: composite fitness spread is too large and mean
  high-order WHT fractions are too low in all three assays. GB1 provides only
  14/24 resolved spectra and 6/24 beneficial-slope summaries. Model rankings
  change with resolution. No empirical interaction identities, universal rates
  of diminishing returns, or physical free energies are inferred.
- Construction and paired-query costs are feasible through length 512 and latent
  mutation RMS stays near one. However, the separate default-parameter scaling
  panel saturates: observable RMS is already around 1e-49–1e-45 at length 32.
  Stable RMS reveals nonzero effects at length 128 despite recorded zeros from
  squaring underflow; all sampled effects are zero in float64 at length 512.
  This is not a large-length transfer test of the fitted four-site models.
- Frozen-source scaling replay matches recorded RMS numerically but not the
  original coefficient fingerprints across platforms; both are retained. The
  analysis does not claim bitwise replay or arbitrary-length observable adequacy.
  It adds no historical success gates and creates no RES/CLAIM record.
- `labproj analysis sync` updated the notebook companion and 16-artifact index;
  `analysis complete` wrote completed status before its repository-wide check
  stalled on historical cloud-backed files. Global status/refresh/check attempts
  likewise timed out. Targeted analysis audits passed; no full strict-check pass
  is claimed. `EXP-065` is also marked completed after run and analysis completion;
  its status command saved the update before hitting the same global I/O issue.

## 2026-10-07 EXP-065 queue submission

- At the scientist's request, submitted the full EXP-065 configuration to
  `lab-slurm`: live `RUN-164`, Slurm job `2311`, CPU partition, 2 CPUs, 8 GB RAM,
  and a 2-hour limit. The scheduler confirmed `RUNNING` on `dutch`; the job has
  emitted comparison progress for all three assays and stderr was empty at the
  startup check. No scientific conclusion is inferred from those progress lines.
- The job subsequently completed while submission verification was finishing.
  Targeted reconciliation records `RUN-164` as completed with exit code `0`.
  All eight declared output files exist remotely; `summary.json` reports
  `DESCRIPTIVE_CONFIRMATION_COMPLETE` with nine assay/family comparisons. Outputs
  have not yet been collected or scientifically analysed, and no ANA/RES record
  was created in this submission task.
- `RUN-163` is the successful dry run. The first two `RUN-164` attempts failed
  while assembling the local source archive with cloud-filesystem timeout errno 60,
  before any scheduler submission. Attempt 3 staged only the experiment's required
  source and declared inputs and was accepted. The exact staging filter is retained
  at `experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-164.staging-filter.txt`;
  the original project `.labproj-stageignore` was restored after submission.
- Added an explicit data-root argument and CSV/Parquet input support for the
  cluster's existing data bundle. All three cluster Parquet SHA-256 hashes match
  the local copies, and sampled nodes/measurements match their CSV versions.
  The 20 focused generator/measurement/runner tests pass after these launcher
  changes. Staged runner and generator hashes match the local submitted source.
- Record-wide status/refresh/check operations still time out on historical files;
  targeted `labproj run reconcile RUN-164` succeeds. Submission required no change
  to scientific parameters, confirmation seeds, or the full experiment budget.

## 2026-10-07 implementation and EXP-065

- Implemented `adaptenv.ThermodynamicLandscape` and `ThermodynamicConfig` in
  `external/Adapt-Env/adaptenv/thermodynamic.py`. The candidate combines a
  stability-like additive trait, the extensive operator-normalized binary-rank-8
  field from EXP-060/063, and a dimensionless sigmoid readout. Shared GNK code
  now lives in the installable package; project imports remain compatible.
  Batch and incremental mutation scoring avoid peak searches and sequence-space
  enumeration. Existing `FitnessLandscapeEnv` and service endpoints are unchanged.
- The additive prior is explicitly synthetic, with an optional independently
  sampled empirical stability-effect distribution. This is biological motivation,
  not a claim that inferred latent values are physical folding energies. See
  `external/Adapt-Env/docs/THERMODYNAMIC_V2.md` for the equation, costs, literature,
  and limitations. RES-058/060/061 motivate the components; empirical adequacy
  of this combined implementation remains unestablished.
- Measurement version `aggregate_v2` separates strict peaks from neutral components,
  reports unresolved WHT fractions as unavailable, distinguishes total epistatic
  variance (orders >=2) from high-order variance (>=3), and requires explicit DFE
  direction. Same-mutation slopes separate beneficial/deleterious strata and offer
  independent-replicate covariance correction. EXP-065's inputs lack the independent
  replicate errors needed for that correction, so its slopes are descriptive.
- Created `EXP-065` through labproj as a planned, local experiment. It compares
  additive-plus-sigmoid, bare GNK, and composite families on GB1, TrpB, and CR9114-H1,
  using 24 calibration and 24 genotype-disjoint confirmation cubes per assay,
  four calibration seeds and eight fresh confirmation seeds. Selection uses joint
  descriptor distance; outputs retain descriptor errors, resolution sensitivity,
  and finite-length scaling. No empirical interaction identities are fitted.
- The scientist explicitly authorized replacing arbitrary prior gate requirements
  for this new protocol. EXP-065 therefore does not require routing specificity,
  one universal setting, or a numerical pass label. Historical results and decisions
  are preserved; the review below records the earlier state before implementation.
- Verification: 46 focused tests pass (new generator/statistics/runner plus existing
  dense-descriptor and GNK tests). Synthetic and measured-data quick runs complete;
  the final measured quick run covers all three assays, writes all eight declared
  outputs, and is explicitly `SMOKE_ONLY`. CR9114-H1 has 442 nonfinite fitness rows,
  excluded under a declared policy with counts in the sampling audit. Complete
  measured cubes only are used. This is engineering verification, not confirmation
  of the scientific prediction. No production RUN, analysis, or result was created.
  A separate length-512 check built 4,310 supports in 5.86 seconds and evaluated
  one substitution at each of 512 sites in 0.232 seconds; all effects were finite.
  These single-machine timings demonstrate an executable path, not a scaling-law
  estimate or empirical validation.
- Record maintenance: `labproj status` printed EXP-065 among 65 experiments but
  did not terminate within 45 seconds. `refresh --write --strict` and
  `check --strict` also timed out at 45 seconds, so no repository-wide strict
  check pass is claimed. EXP-065 YAML, all 18 declared input paths, Python syntax,
  and the executable shell entry point were checked directly; its wrapper completes
  a measured-data quick run. Engineering seeds are now separate from the production
  panel; fresh production confirmation seeds replace those seen in early quick checks.

## Prior scientific results

- Status: active aim, hypotheses recorded; `EXP-041` completed with `ANA-040`
  and `RES-040`; `EXP-042` completed with `ANA-041` and `RES-041`;
  `EXP-043` completed with `ANA-042` and `RES-042`; `EXP-044` completed with
  `ANA-043` and `RES-043`; `EXP-045` completed with `ANA-044` and `RES-044`;
  `EXP-046` completed with `ANA-045` and `RES-045`; `EXP-047` completed with
  `ANA-046` and `RES-046`; `EXP-048` completed with `ANA-047` and `RES-047`;
  `EXP-049` completed with `ANA-048` and `RES-048`; `EXP-051` completed with
  `ANA-049` and `RES-049`; `EXP-052` through `EXP-056` completed Step 4A-4E
  ensemble landscape validation with `ANA-050` through `ANA-054` and `RES-050`
  through `RES-054`; `EXP-057` completed the design-matched support-aware Stage
  4B repair probe with `ANA-055` and `RES-055`; `EXP-050` completed the
  non-additive PHOT measurement-operator certification run with `ANA-056` and
  `RES-056`; `EXP-058/RUN-142` completed the registered readout-coupling
  generator sweep on attempt 3 as Slurm job `2003` after two failed attempts,
  analysed in `ANA-057` with `RES-057`. All three fixes fail all three
  preregistered gates (`0/5022` core-guardrail, `0/5022` GNK-separation), the
  terminal preregistered decision fires, and the Stage 4C/4D/4E hard stop holds.
  `EXP-059` is registered as the HYP-005-only saturated-gate two-channel
  viability benchmark; `RUN-143` is its successful dry-run submission record,
  while live `RUN-144` failed before reaching Slurm because the sandbox could
  not resolve the SSH host. The network-enabled retry `RUN-145` completed on
  `lab-slurm` as Slurm job `2038` with exit code `0` and took the preregistered
  `NO_COMPONENT_QUALIFIED_POINT` branch. Its outputs and logs are collected,
  DVC-tracked, and backed up. `EXP-060` is registered as the HYP-005-only
  extensive, mutation-operator-normalized GNK latent scale-law test. After
  initial and revised dry runs (`RUN-146`, `RUN-148`) and a pre-computation
  staging failure (`RUN-147`, Slurm job `2039`), corrected live `RUN-149` was
  submitted to `lab-slurm` as Slurm job `2040` and completed with exit code
  `0`. `ANA-058` collected and independently audited the complete output bundle,
  and `RES-058` records all `26/26` required latent gates passing. The evidence
  supports mutation-operator-normalized GNK scaling only before sigmoid
  transmission. `EXP-061/RUN-151` has now completed the paired-residual
  sigmoid-transmission test on `lab-slurm` as Slurm job `2041` with exit code `0`.
  `ANA-059` verifies the 193-file DVC bundle and reproduces the registered
  `SIGMOID_TRANSFER_FAILED` decision (`13/16` required checks). Pointwise scalar
  transfer and frozen-core texture pass, but common-threshold edge coverage,
  actual-core texture, and margin-8 deployment fail. `RES-059` therefore weakens
  HYP-005 at the observable protein-scale boundary without invalidating the
  RES-058 latent result or showing that sigmoid transmission is impossible in a
  controlled transition window. The public generator/service remains unvalidated
  and the Stage 4C/4D/4E hard stop remains. `EXP-062` is now registered as the
  dimensionless additive-core calibration and held-out GNK dose-response test;
  its runner and fixed real-protein calibration/confirmation panel pass focused
  tests and a local quick smoke. Dry-run `RUN-152` was inspected, and live
  `RUN-153` completed on `lab-slurm` as Slurm job `2042` with exit code `0`;
  its 373-file output/log/checkpoint bundle is collected, SHA-256 verified,
  DVC-tracked, backed up, and analysed in `ANA-060`. `RES-060` records the formal
  `CALIBRATED_SIGMOID_TRANSMISSION_FAILED` decision (`6/9` required checks):
  global core calibration passes `60/60` held-out cells and complete-h1 edge
  transfer passes `360/360`, but only `24/360` actual-core and `0/360`
  frozen-core cells meet the registered texture-prevalence floor. Exact
  per-landscape normalization does not rescue texture. The recorded
  directionality failure is a floating-point equality artifact and passes under
  a `1e-12` audit tolerance (`7/9`), without changing the overall result.
  HYP-005 is weakened at the full observable texture-prevalence boundary; the
  public generator/service remains unvalidated and the Stage 4C/4D/4E hard stop
  remains. `EXP-063/RUN-157` has now completed the paired categorical-code and
  representation-rank repair test on `lab-slurm` as Slurm job `2044` with exit
  code `0`. `ANA-061` independently reproduces all 30 recorded pass states and
  the preregistered `BINARY_RANK8_SUFFICIENT` decision. `RES-061` records that
  binary rank 8 raises median practical texture prevalence from `0.198` to
  `0.461` and passes `360/360` cells at nearly unchanged texture RMS, whereas
  dense q20 rank 1 worsens prevalence and dense q20 rank 8 misses the
  per-length gate. This supports HYP-005 narrowly for an experiment-local GNK
  representation on the reused challenge panel; the public method remains
  unchanged and the Stage 4C/4D/4E hard stop remains.
- Active aim: `AIM-001 Create-Model-Landscape`

## 2026-10-07 review before implementation

- The scientist clarified the immediate target: a synthetic, procedurally scalable,
  biologically informed generator calibrated to aggregate empirical statistics,
  including WHT variance partitions and diminishing returns. Recovering actual
  residue-specific interactions is explicitly outside this target. Aggregate
  statistical matching does not require one universal setting across assay classes.
- Recommended candidate, not an implemented or validated change: retain one
  stability-like latent coordinate with a nonlinear readout and incorporate the
  mutation-normalized extensive, binary-rank-8 GNK construction from EXP-060/063.
  RES-060/061 support component calibration and transmission, not joint empirical
  adequacy or public-generator readiness. RES-062 remains a recorded failure.
- Code inspection and small deterministic diagnostic checks confirm that
  `dense_descriptive.local_extrema_metrics` counts all 16 vertices of a constant
  4-cube as maxima. Separate unresolved/neutral plateaus from strict peaks before
  rebuilding the empirical target. The Stage-4 `directed_edge_effects` function
  includes both signs of every edge, giving a symmetric edge-effect distribution;
  this must not stand in for a reference- or population-conditioned DFE.
- The next protocol should test joint descriptor recovery rather than require
  structured routing to exceed shuffled/random routing. Retain these as comparison
  controls, consistent with CLAIM-003/HYP-009. Any revised protocol must explicitly
  supersede the relevant gates; this review does not retroactively pass old runs
  or lift their recorded hard stops.
- Scientific qualifications: observed WHT energy includes nonlinear-readout effects;
  an order spectrum alone does not determine landscape topology. A strictly
  monotone scalar readout preserves latent edge signs and strict topology before
  numerical or experimental resolution limits. Continuous concavity alone does not
  imply a single maximum on a discrete Hamming graph. Single-mutant-only data do
  not identify diminishing returns across genetic backgrounds.
- Proposed shortest work sequence: repair and freeze measurement/sampling rules;
  consolidate the tested candidate behind one generator implementation with bounded
  evaluation costs; calibrate a small parameter family against joint empirical
  summaries and confirm with independent simulation seeds and appropriately held-out
  empirical evidence. These are review recommendations, not new experiment records,
  completed model changes, or claims of empirical validation.
- Review verification used source inspection and direct deterministic checks. The
  project pytest invocation was interrupted while importing dependencies before
  tests ran; `labproj status` did not complete, and bounded refresh/check attempts
  timed out. No test-suite or strict-record-check pass is claimed.

## Current scientific state

- The parent project still has one active aim and had no project-local hypotheses before this review.
- The practical scientific state currently lives in `external/Adapt-Env`, which already contains a synthetic landscape implementation, a population-genetics layer, and a validation suite aimed at qualitative protein-like behavior.
- The current model claim that is best supported is narrow: AdaptEnv is a biologically informed synthetic landscape analogue for studying local statistics and evolutionary dynamics, not an exact reconstruction of a real protein fitness landscape and not yet a machine-learning benchmark claim.
- `RES-040` sharpens the PHOT model-misspecification picture: the existing
  `RUN-039` full-ABC branches are not broadly wrong across the exact
  complete-graphlet package, but they consistently miss local edge/site
  Dirichlet-energy concentration. Any new PHOT ABC target should start from
  those exact descriptive failures, not from the full older validation vector.
- `RES-047` is the prior PHOT calibration update: the EXP-048 sparse
  Walsh/Potts Generator A test completed as an execution but failed the
  preregistered H1 forward capability gate. The empirical spectrum calibration
  put all normalized sparse-Walsh weight on order `3`; the H1 probe produced
  `0/12` joint passes, failing accessibility, magnitude-epistasis ratio, and
  Walsh high-order fraction for every candidate. The sequence-conditioned fit
  was therefore correctly skipped, so there is no sparse-Walsh equivalence row
  to compare against `RUN-098`.
- `RES-048` is the current PHOT measurement-operator update. EXP-049/RUN-125 completed
  the P0 target-operator isolation motivated by
  `scratch/RES-047_miscalibration_root_cause_and_isolation_plan.md`: P0.1
  showed cached/refit PHOT MAVE-NN latent targets are smooth/additive-like,
  P0.2 showed the identical MAVE-NN plus `dense_descriptive` operator
  smooths known injected spectra into nearly unimodal and epistasis-free
  landscapes, and P0.3 showed sparse-Walsh support containment is not the
  limiting explanation under the audited support rule. This weakens the
  current PHOT adequacy route and motivates explicit measurement-operator
  calibration before using these scalar targets for benchmark construction.
- `RES-049` is the current hashed GNK forward-generator update. EXP-051/RUN-127
  supports the WS-B knob-map claim in the narrow model-native setting: all
  `18/18` preregistered monotonicity checks passed on enumerable binary `L=8`
  coordinates, and `epi_weight = 0.6` gave a joint multi-peaked but navigable
  regime. It supports `HYP-005` for directional knob control, while motivating
  rather than supporting `HYP-004`/`HYP-009` because it does not yet test
  empirical assay overlap, a certified measurement operator, or optimizer
  behaviour.
- `RES-056` completes the EXP-050 non-additive measurement-operator
  certification follow-up to `RES-048`. The corrected `RUN-129` produced all
  expected tables and exited cleanly, but certified no pairwise or black-box/GE
  MAVE-NN PHOT measurement operator: `0/2` operators and `0/24` certification
  checks passed. Pairwise fits were preflight-skipped at PHOT length `118`, and
  black-box/GE fits failed under the installed MAVE-NN/Keras
  `MultilayerPerceptronGPMap build()` interface. This weakens the current
  measurement-operator route under `HYP-004`, `HYP-006`, and `HYP-009`, without
  refuting biological PHOT ruggedness or optimizer/agent hypotheses.
- `RES-050` through `RES-054` are the current ensemble-landscape validation
  update for Step 4A-4E. `RES-050` establishes the empirical reference envelope
  and the required scope split: `10` Tier-1 operator-free landscapes can support
  higher-order/topographic descriptors, while `7` sparse Tier-3 assays remain
  order-1 only. The later stages are boundary or negative for the present
  public generator: `RES-051` weakens the broad protein-scale `HYP-005`
  monotone-knob claim (`20/54` checks passed), `RES-052` finds Tier-1
  descriptor contact but not public-model-specific coverage, `RES-053` fails
  the hard joint-realism gate (`0` joint passes), and `RES-054` does not
  validate held-out descriptor generalisation (median held-out fraction in band
  `0.5`, only `2/10` strict-majority/`0.6` held-out successes).
- The `external/Adapt-Env` public model surface has now been updated to make
  the additive thermodynamic stability core plus deterministic hashed GNK
  texture the default and base system. Public presets, service state, package
  exports, tests, and docs now point at `model_system =
  additive_thermodynamic_gnk`, `observed_fitness_combine_mode =
  stability_gate`, and explicit `gnk_*` controls. Older FGM, pairwise-J,
  score-space, exact-match, deterministic-noise, sparse-Walsh alias, and
  peak-design paths remain only as private legacy compatibility routes that
  warn they will be deleted soon when used.
- `RES-055` is the current fixed-probe Component 1 update. `EXP-057/RUN-141`
  completed the design-matched support-aware Stage 4B repair probe after the
  first live attempt failed by Slurm walltime and the final resumed attempt
  completed as Slurm job `1955`. The probe successfully used support-aware,
  dimension-matched, alphabet-matched graphlets and regime-stratified anchors,
  but the fixed gates failed: only `318/1620` monotonicity checks passed and
  only `34/1080` positive-epistasis core-vs-GNK separation rows passed. This
  weakens the current protein-scale public additive-core plus GNK implementation
  under `HYP-005`, `HYP-004`, and `HYP-009`, and blocks downstream Stage
  4C/4D/4E reruns for this implementation.
- `EXP-058` is now the registered Component 1 calibration response to the
  `RES-055` readout-coupling diagnosis. It is operator-free and sweeps the
  three fixes specified in
  `scratch/RES-056_assessment_readout_coupling_and_next_steps.md`: Fix A
  desaturates and re-centres the single stability gate, Fix B tests mixed
  score/stability GNK routing, and Fix C tests the two-channel product
  observable. The registered decision requires core guardrail, GNK separation,
  and knob-monotonicity gates simultaneously, and it keeps Stage 4C/4D/4E
  reruns blocked until GNK separation passes. `RUN-142` attempt 1 failed by
  Slurm walltime as job `1991`; the runner now reuses per-parameter
  environments, writes per-parameter checkpoints, and uses a finite broad Fix C
  function channel. Attempt 2 job `1999` reached `297/5022` checkpointed
  parameter sets with clean stderr, but observed throughput showed another
  `48:00:00` timeout was likely, so it was cancelled and reconciled failed.
  Attempt 3 completed on `lab-slurm` as job `2003` with scheduler exit code
  `0`, evaluating all `5022` parameter sets. The 11 expected outputs and three
  logs are collected under
  `data/processed/generator-readout-coupling-calibration/RUN-142`, DVC-tracked,
  and backed up. `ANA-057` and `RES-057` now record the outcome: all three
  fixes fail all three preregistered gates (`0/5022` core-guardrail, `0/5022`
  GNK-separation, `11,883/58,320` monotonicity with wrong-direction medians), so
  the terminal preregistered decision fires
  (`operator_free_combinatorial_realness_claim_fails_for_this_model_family`) and
  the Stage 4C/4D/4E hard stop holds. The binding constraint is that no
  high-fitness plateau forms in the re-centred margin regime (best observed
  high-plateau median fitness `0.049` on a `0-10` scale): plateau-based
  diminishing returns and transmitted GNK-specific epistasis are not obtained
  together. Fix B mixed routing does transmit large observable high-order/sign
  epistasis (`|full - core|` deltas up to `0.87`/`0.16`) but on a collapsed,
  plateau-less landscape; Fix C two-channel product transmits almost nothing
  because it was swept only in the non-saturated margin regime, so its stability
  gate collapsed the product. `RES-057` weakens `HYP-004`, `HYP-005`, and
  `HYP-009` and flags a saturated-stability-gate two-channel construction as the
  next probe.
- `EXP-059` implements that scoped next probe using the private legacy
  two-channel product construction. A single sequential run first selects the
  first component-qualified operating point from the frozen `m14/m18/m22` by
  steepness `2/4` grid using only stability-gate occupancy, function-channel
  capacity, raw-GNK variance, reconstruction, and dead/closed controls. Only if
  that calibration passes does it test held-out seeds `1701-1708` on frozen
  GB1, DHFR, CR9114, and ParD2-shaped graphlets. Confirmation is explicitly
  conditional on stability-ranked, viability-enriched support-aware subspaces;
  it is not evidence for random genotype neighbourhoods, empirical fitness
  realism, the public model, or optimizer performance. The `w=0` null,
  structured primary route, mixed-score comparator, shuffled/random routing
  diagnostic, seed-level G1-G4 gates, exact-sign test, and Stage 4C/4D/4E hard
  stop are preregistered. Focused validation passes (`11` tests plus quick
  calibration, resume, and forced-confirmation executions). `RUN-143` staged
  the complete payload in dry-run mode; `RUN-144` is retained as a failed
  pre-scheduler transport attempt caused by sandbox DNS resolution. `RUN-145`
  completed on `lab-slurm` as Slurm job `2038` with exit code `0`. No operating
  point passed component preflight in all four calibration seeds: pass counts
  were `0/4`, `0/4`, `1/4`, `1/4`, `2/4`, and `2/4` in the fixed order. At
  both `m22` points, gate occupancy was `1.0` in every seed, but the
  function-capacity interior fraction fell below the `0.65` floor for seed
  `701` at `w=1` (`0.64685`) and seed `703` at `w=0.6` and `w=1` (worst case
  `0.535994`). Capacity SD, raw-GNK variance, reconstruction, dead-state, and
  closed-gate controls passed. The registered early stop therefore fired
  before descriptor calculation or held-out confirmation; G1-G4 and the
  GNK-to-epistasis/topology prediction are unevaluated, not failed. The 15
  outputs and three logs are collected under
  `data/processed/saturated-gate-gnk-mechanism-calibration-confirmation/RUN-145`,
  DVC-tracked, and backed up.
- `EXP-060` isolates the first proposed repair layer before any nonlinear
  readout: extensive order-specific term counts, unique nearly regular
  supports, bounded balanced q=20 contrasts, exact normalization under the
  uniform non-self mutation operator, and reference centering. It compares the
  candidate against both the public legacy scaling and a support-localized
  diagnostic control across six real-protein reference lengths (`48`-`448`),
  30 generator seeds, and fixed and fractional Hamming distances in uniform-
  background and reference-shell regimes. No sigmoid, observable fitness,
  accessibility, or optimizer is included, and the public AdaptEnv model is
  unchanged by this experiment. A separate sigmoid-transmission experiment
  remains blocked unless all registered latent gates pass. `RUN-147` failed
  before computation because a cluster-side Stage-4 reference file was absent;
  the same six hash-verified sequences were then frozen inline without changing
  the scientific design. Revised live `RUN-149` completed on `lab-slurm` as
  Slurm job `2040` with exit code `0`; its runner summary reports `180/180`
  tasks and `26/26` required latent gates complete and passing. `ANA-058` verifies
  all 193 collected files against their remote SHA-256 digests and independently
  reproduces the registered decision; the bundle is DVC-tracked. `RES-058`
  supports `HYP-005` narrowly: the candidate local slope is `-0.00632` versus
  `-0.51036` for the frozen Stage-4 raw-amplitude control, analytic field-variance
  slope is `0.99508`, and reference-shell fractional-distance GNK RMS slopes are
  `0.48196/0.43034`. The locally normalized fixed-support control also repairs
  local, field, and full-distance scaling, identifying operator normalization as
  sufficient for the scale repair; extensive balanced support and bounded
  contrasts add coverage and tail control but were not necessity-ablated. The
  runner omitted the prose-specified fixed-distance `h = 1` gate, but independent
  analysis gives slope `-0.001481`, inside the same limit. The result remains
  latent-only, conditional on one fixed protein per length, and does not overturn
  the observable-level `RES-055`/`RES-057` failures or lift the Stage 4C/4D/4E
  hard stop.
- `EXP-061` is the completed readout-layer follow-up to `RES-058`. Its primitive
  estimand is the paired residual `R_w(x) = F(A(x)+wG(x)) - F(A(x))`; mutation
  and two-site contrasts use `Delta R` and `epsilon_R`, respectively. It reuses
  the six frozen EXP-060 proteins and 30 generator seeds, enumerates every
  reference single substitution, and pairs the candidate with fixed-support raw
  and locally normalized controls plus a mutation-RMS-matched order-1 surrogate.
  Three scopes are kept distinct: genotype-wise occupancy calibration tests the
  scalar transfer function only, reference-anchored thresholds test edge and
  plaquette transmission in one landscape, and deployment thresholds test the
  documented margin-8 base preset plus the bare-constructor margin-10 comparator.
  Stable analytic and direct float64 residuals distinguish attenuation from exact
  plateau ties. A controlled pass with deployment failure is preregistered as a
  transition-window-only result, not public validation; all empirical and
  optimizer hard stops remain in force under every EXP-061 branch. `RUN-151`
  completed all `180/180` reference/seed tasks with exit code `0`; the 193-file
  output/log/checkpoint bundle is DVC-tracked and remotely backed up by verified
  content-addressed transfer. `ANA-059` reproduces all 16 checks: all 12 technical
  gates pass, but the controlled-edge, actual-core texture, and margin-8 deployment
  gates fail. The fixed `n >= 32` edge rule is partly length-confounded through
  the `19L` panel size, yet all-edge practical transmission remains only
  `1.48-1.67%`. Frozen-core texture passes at every length (slope `-0.10100`, 95%
  interval `[-0.11819, -0.08858]`), whereas actual-core texture has zero practical
  fraction throughout and margin-8 edge/texture practical fractions never exceed
  `0.0010/0.0254`. `RES-059` weakens HYP-005 for this observable construction:
  the sigmoid can carry normalized texture in a controlled dynamic window, but
  real-core motion and deployment saturation prevent a broadly resolved knob.
- `EXP-062` tests the specific operating-scale repair implied by `RES-059` without
  searching on GNK outcomes. Phase 1 takes one median complete-h1 additive-effect
  RMS over 12 calibration proteins and five core seeds, then selects the largest
  eligible `beta` and `mu` from an additive-core-only grid while requiring a high
  reference, dynamic h1 targets, and low plateaus at `0.1L` and `0.5L` in every
  length stratum. Phase 2 freezes that choice before evaluating 12 disjoint
  held-out proteins under five new core seeds crossed with six new GNK seeds. The
  primary estimands are paired edge residuals and candidate-additivized plaquette
  texture over `a=0.25-2`; `a=4` is saturation stress. A per-landscape exact-beta
  arm is diagnostic only. Frozen-core, order-1 additive-surrogate, zero-dose,
  steepness-only, threshold-only, raw margin-8/margin-10, and algebraic-parity
  controls are preregistered. The candidate remains experiment-local, no empirical
  outcome is loaded, and every public-model, realism, accessibility, optimizer,
  and Stage 4C/4D/4E hard stop remains in force. Focused tests pass (`47` across
  EXP-060/061/062 helpers), the full core-only calibration audit deterministically
  selects `beta=4`, `mu=logit(0.99)` without GNK data, and quick mode completes
  end to end while correctly issuing no scientific decision. Inspected dry-run
  `RUN-152` and live `RUN-153` preserve the registered 8-CPU, 64-GB, 24-hour
  resource request; `RUN-153` was accepted as Slurm job `2042` and completed
  with exit code `0`. Its ten registered outputs, 360 task checkpoints, and
  three logs are collected, DVC-tracked, and backed up as a 373-file,
  332,322,278-byte bundle with directory hash
  `02986be42c5b575ff8d7b75b3abad31c.dir`; every local file matches
  shared-storage SHA-256 and DVC MD5. `ANA-060` independently
  reproduces the formal `6/9` decision. The core-only choice
  (`sigma_A=2.74486`, `beta=4`, `mu=logit(0.99)`) passes all `60/60` held-out
  core cells, and all `360/360` complete-h1 edge cells pass. Actual-core texture
  remains aligned (median RMS `0.4175`, correlation `0.7322`, projection
  `0.9182`) but only `24/360` cells reach the `0.25` practical-prevalence floor;
  the frozen control reaches it in `0/360` because its `0.99` reference
  occupancy is saturated. Exact beta is not a rescue (`23/360`). A `1e-12`
  sensitivity changes the directionality gate from fail to pass by treating
  stored `0.7999999999999999` Spearman values as the registered mathematical
  boundary `0.80`, giving `7/9` but leaving the overall decision unchanged.
  `RES-060` therefore weakens the full observable protein-scale HYP-005 claim
  while retaining narrow fixed-panel evidence for global calibration, edge
  transfer, dose directionality, and bounded length stability.
- `EXP-063` is the registered one-step representation intervention motivated by
  `RES-060`. It freezes the EXP-062 additive core, global sigmoid operating
  point, support topology, per-order uniform-q20 nonself mutation energy, dose
  grid, probe genotypes, and gates, while crossing dense q20 versus balanced
  binary allele coding with rank 1 versus rank 8 on the same supports and base
  coefficient draws. The exact 12 EXP-062 confirmation proteins are reused as a
  paired failure-challenge panel with fresh core/GNK seeds; this is mechanistic
  development evidence, not independent protein confirmation. The runner adds
  exact Gram normalization including same-support cross-component terms,
  incumbent binary-rank-1 replay, genuinely unsupported-pair controls, a
  midpoint-sigmoid positive control, paired protein-level rescue estimates, and
  factorial coding/rank contrasts. The candidate is experiment-local and does
  not modify the public AdaptEnv implementation. All 48 focused tests pass and
  quick mode completes all declared outputs with the required
  `QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION` branch. Repository-wide project
  tests have 127 passes plus three unrelated existing failures (two legacy
  calibration-mode assertions and a sandbox-blocked Ray process probe);
  optional `external/landscapy` collection also lacks `cogent3`. `RUN-154` and
  revised `RUN-156` are inspected dry runs. `RUN-155` was rejected before Slurm
  job creation because a 16-CPU node was unavailable; resource-only retry
  `RUN-157` completed as Slurm job `2044` with 8 CPUs, 64 GB, 24 hours, and exit
  code `0`. Its 15 outputs, 360 strict-JSON checkpoints, and three logs are
  collected and locally DVC-tracked as a 378-file, 456,516,694-byte bundle with
  directory hash `7f2e06c1c337fa49b0e9c4006be443ff.dir`; local/cache MD5 and
  shared-runtime SHA-256 audits pass. The configured redundant DVC push failed
  with a partial remote (`155/378` content objects and no directory object), so
  that backup cannot reconstruct the bundle, although the local cache and
  shared-runtime copy remain intact. `ANA-061` reproduces all 30 pass states and
  the registered `BINARY_RANK8_SUFFICIENT` decision. The fresh `binary_r1`
  arm reproduces the failure (`45/360` texture passes; median practical fraction
  `0.1979`), while `binary_r8` passes `360/360` with median practical fraction
  `0.4609` and paired gain `0.2624` (bootstrap lower `0.2480`) at nearly
  unchanged texture RMS. Rank 8 redistributes fixed RMS energy away from
  near-zero, heavy-tailed plaquette effects. `dense_q20_r1` fails `0/360` and
  its midpoint control; `dense_q20_r8` improves but fails the long-bin absolute
  gate. `RES-061` therefore supports HYP-005 narrowly for the experiment-local
  binary-rank-8 construction on this reused 12-protein panel, not as minimum-rank
  evidence, independent protein generalization, or public-model validation. The
  Stage 4C/4D/4E hard stop remains.
- `EXP-064` is now registered as the immediate forward-generation coverage gate
  after `RES-061`. It asks whether the experiment-local additive-core plus
  binary-rank-8 extensive-GNK observable can cover the frozen RUN-135 Tier-1
  ruggedness envelope when `beta` and reference occupancy `p0` are swept rather
  than fixed at the saturated deployment preset. It uses all 10 operator-free
  Tier-1 templates and excludes the seven Tier-3 assays from topographic claims.
  The generated and empirical graphlets share the exact Stage-4A
  `metrics_from_cube` descriptor path. Per-template q05/median/q90/q95 bands are
  resolved from the frozen RUN-135 graphlet table, with a hard reproduction
  check against the published ANA-050 pooled bands and dataset medians before
  any scientific reading. The paired design excludes GNK amplitude from the
  design key and compares structured binary-rank-8 coverage against `a=0`
  additive, shuffled-support, and random-support controls. G1-G4 separately test
  range coverage, single-configuration rugged-end reach, specificity, and
  nonzero accessibility. A pass only returns `PROCEED_TO_FITTING` and unblocks
  the separate proposed EXP-065 fitting experiment; it is not a fit, monotone-
  knob, realism, public-generator, or optimizer result. The public AdaptEnv
  implementation is unchanged and the Stage 4C/4D/4E hard stop remains under
  every outcome. The one-template quick smoke emitted all declared outputs,
  reproduced the frozen bands, passed pairing/energy/rank/incumbent-replay
  checks, and correctly returned `QUICK_SMOKE_ONLY_NO_SCIENTIFIC_DECISION`.
  The final dry-run record is `RUN-159`. `RUN-160` and `RUN-161` are retained
  as failed pre-scheduler transport attempts with no scheduler job ID; diagnosis
  showed that source staging was traversing and archiving prior multi-gigabyte
  processed outputs. `.labproj-stageignore` now excludes generated processed
  directories while retaining the required DVC pointer and formal inputs.
  The corrected live record `RUN-162` was submitted to HomeLab Slurm as job
  `2045` with 8 CPUs, 64 GB RAM, and a 24-hour limit. It completed with inferred
  exit code `0`; all 1,517 files (204,553,256 bytes), 300 checkpoints, and
  259,200 generated grid rows are collected and DVC-tracked under
  `data/processed/observable-ruggedness-coverage-probe/RUN-162.dvc` with
  directory hash `8428ff92c79a7d5b768a97772dbce230.dir`. All 1,526
  technical checks pass, and ANA-062 independently reproduces every gate and
  the `OBSERVABLE_CONSTRUCTION_LIMIT` decision. G1 passes for `1/3` rugged-end
  templates, G2 and G4 for `2/3`, and G3 for `0/6` descriptors, so EXP-065
  remains blocked. Local-maxima density is covered in `0/10` Tier-1 templates;
  sign and reciprocal-sign epistasis are each covered in only `2/10`; and
  structured reachable maxima do not separate from matched shuffled/random
  controls. RES-062 weakens the broader forward-coverage form of HYP-005 but
  withholds the strongest biological construction-limit reading because the
  registered rugged-end selector is plateau-confounded: CR9114-FLUB is flat in
  `96/96` frozen graphlets and CR9114-H3 in `69/96`, while the strict
  no-fitter-neighbour rule labels a flat 4-cube as local-maxima density `1.0`.
  The public implementation is unchanged and the Stage 4C/4D/4E hard stop
  remains.

## 2026-06-18 EXP-050/EXP-051 forward benchmark generator plan

- Following `scratch/RES-048_forward_plan_tunable_benchmark_generator.md`,
  `EXP-050` and `EXP-051` were created as the two immediate workstreams:
  `EXP-050` for non-additive measurement-operator certification and `EXP-051`
  for direct knob-map validation of a hashed GNK benchmark generator.
- `EXP-051` implements the intended replacement for the failing residual
  epistasis module as a deterministic hashed GNK texture layer, while retaining
  the additive stability/Fuji backbone and saturating readout. `RUN-127`
  completed on `lab-slurm` as Slurm job `1935` and produced the preregistered
  enumerable small-L outputs, including `graphlet_metrics.csv`,
  `monotonicity_checks.csv`, `support_catalog.csv`,
  `hamming_layer_metrics.csv`, `threshold_epistasis_metrics.csv`, and
  `core_readout_guardrail_checks.csv`.
- `ANA-049` copied the full compact `RUN-127` output bundle and derived the
  EXP-051 acceptance summaries. `RES-049` records the scoped result: all
  direction checks passed and the retained core/readout still gives declining
  beneficial fraction and readout plateaus at `epi_weight = 0`, but the
  threshold sign-epistasis enrichment claim is mixed rather than a clean pass.
- `EXP-050` reuses the RES-048 round-trip harness with non-additive operator
  profiles. The first live attempt, `RUN-126`, failed by OOM during full
  pairwise MAVE-NN construction. The runner was then updated to preflight-skip
  infeasible full-length pairwise profiles and record them as
  `fit_skipped_preflight` instead of allowing an OOM kill. The 64 GB pending
  resubmission `RUN-128` was cancelled before execution because it would wait
  for `ripley`; it is marked failed in the run registry.
- The corrected EXP-050 run, `RUN-129`, completed on `lab-slurm` as Slurm job
  `1937`. It produced all expected WS-A output tables. The pairwise profile is
  recorded as preflight-infeasible at PHOT length 118, and the black-box/GE
  MAVE-NN profile failed under the installed Keras/MAVE-NN interface with a
  `MultilayerPerceptronGPMap.build` signature error. This is a technical
  operator-certification failure, not support for an empirical benchmark
  target. `ANA-056` collected the compact run outputs and logs under
  `data/processed/phot-nonadditive-measurement-operator-certification/RUN-129`,
  and `RES-056` records the scoped negative result: no tested non-additive
  PHOT measurement operator was certified.

## 2026-06-18 Step 4A-4E ensemble landscape validation submissions

- Following `scratch/science_plan_part1_ensemble_landscape_validation.md` and
  the public additive thermodynamic core plus hashed GNK texture protocol,
  `EXP-052` through `EXP-056` were created for Step 4A empirical descriptor
  envelopes, Step 4B generator reachability, Step 4C calibrated overlap, Step
  4D joint realism, and Step 4E held-out descriptor generalisation.
- The runner is `scripts/ensemble_landscape_validation_stage4.py`. It uses the
  active processed bundle `data/processed/ensemble-landscape-validation-datasets`
  and blocks higher-order/topographic descriptors for sparse ProteinGym-style
  assays that do not have a complete mutation graph, preventing epistasis claims
  from single-mutation-only or sparse records.
- The planned controls are embedded in the experiment configs: additive core
  only (`epi_weight = 0`), public observable additive-core plus GNK texture,
  bare GNK texture controls, structured/shuffled/random GNK routing controls,
  and measurement-operator scope checks for sparse Tier-3 assays.
- Dry-run submission records are `RUN-130` through `RUN-134`. Live submissions
  completed on `lab-slurm`: `RUN-135` (`EXP-052`, Slurm job `1941`), `RUN-136`
  (`EXP-053`, job `1943`), `RUN-137` (`EXP-054`, job `1944`), `RUN-138`
  (`EXP-055`, job `1945`), and `RUN-139` (`EXP-056`, job `1946`). The first
  `RUN-135` staging attempt was interrupted before Slurm handoff because it was
  sending local data directories; attempt 2 used the slim source staging rules
  and submitted successfully. A duplicate manual fallback job `1942` was
  cancelled.
- Remote jobs use the shared Slurm-side dataset bundle at
  `/tank/data/labproj/runtime/project-adapt-env/shared/ensemble-landscape-validation-datasets`.
  `.labproj-stageignore` now excludes bulky local data directories from source
  staging while leaving source code, configs, and DVC pointer files in the
  submitted project bundle.
- `ANA-050` through `ANA-054` collected the compact outputs for `RUN-135`
  through `RUN-139` into analysis records and DVC-tracked the corresponding
  processed run bundles under `data/processed/stage-4a-empirical-ensemble-descriptor-envelope`,
  `data/processed/stage-4b-additive-core-gnk-generator-reachability`,
  `data/processed/stage-4c-empirical-coverage-calibrated-overlap`,
  `data/processed/stage-4d-joint-multidescriptor-realism-gate`, and
  `data/processed/stage-4e-heldout-descriptor-generalisation`.
- `RES-050` records Stage 4A as a successful empirical reference-envelope and
  scope-audit step only. `RES-051` records Stage 4B as a protein-scale
  reachability boundary, not a successful monotone-knob validation. `RES-052`
  records Stage 4C as descriptor-space contact without public-model specificity.
  `RES-053` records Stage 4D as a failed joint-realism gate, and `RES-054`
  records Stage 4E as failed held-out descriptor generalisation. These results
  should not be used to justify downstream optimizer or PLM-agent benchmark
  claims from the current Stage 4 grid.
- `EXP-057` has now completed the Section 3 repair from
  `scratch/RES-054_assessment_and_path_to_finishing_component1.md`. `ANA-055`
  collected and analysed the compact `RUN-141` output bundle under
  `data/processed/stage-4b-design-matched-support-aware-gnk-probe/RUN-141`.
  `RES-055` records the negative fixed-probe result: support-aware measurement
  worked, but the current public composite did not restore monotone protein-scale
  `epi_weight` control or core-vs-GNK separation. Do not use this implementation
  for downstream overlap, joint-realism, held-out, optimizer, or PLM-agent
  benchmark claims.
- `EXP-058` has been created from
  `scratch/RES-056_assessment_readout_coupling_and_next_steps.md` as the single
  preregistered generator calibration to settle the readout-coupling defect.
  It sweeps non-saturated single-gate calibration, mixed post-readout score
  routing, and the two-channel product observable in one design. Its registered
  gates require the core guardrail and GNK-specific separation to hold together,
  with a hard stop on downstream Stage 4C/4D/4E reruns until separation passes.
  `RUN-142` attempt 1 failed by Slurm walltime as job `1991`. The runner was
  repaired to avoid repeated graphlet-level environment construction, to write
  per-parameter checkpoints, and to keep the Fix C product channel finite rather
  than zero-width. `RUN-142` attempt 2 ran as Slurm job `1999` and was
  cancelled/reconciled failed after `297/5022` checkpointed parameter sets,
  because live throughput showed the original `48:00:00` walltime remained too
  short. With the walltime raised to `168:00:00`, `RUN-142` attempt 3 completed
  as Slurm job `2003` with exit code `0`. The registered output bundle and logs
  have been collected and backed up with DVC; the transient checkpoint tree was
  intentionally left in runtime storage.

## 2026-05-18 coupled-landscape benchmark interpretation

- `CLAIM-003` has been reframed around the defensible result: fitted coupled or auxiliary landscape capacity improves PHOT global/local ruggedness summaries over current fitted and smooth-core controls, but the improvement does not rely on empirical alignment/contact/module locations as specific amino-acid interaction structures.
- The shuffled, random, no-module, and uncoupled controls should be read as interpretation controls, not as a general negative result for coupled landscapes. They weaken empirical module specificity and residue-level mechanism claims; they do not weaken the benchmark use case where the point is to construct matched landscape regimes with controllable ruggedness.
- `HYP-008` is now superseded in its mechanism-specific form. `HYP-009` records the active benchmark-facing hypothesis: biologically routed ruggedness may affect PLM or evolution-aware optimizer performance once global and local landscape statistics are matched against random/shuffled/no-module/diffuse controls.
- Alignment or covariation conditioning is retained only as biologically legible side information for agents. It is not currently claimed to improve landscape-statistic fit over generic capacity controls.

## 2026-05-18 EXP-040 ProteinGym production panel

- `EXP-040` has been created under `HYP-004` and `HYP-009` to test, across a 24-system multi-mutant ProteinGym panel, whether the base thermodynamic residual model, an alignment-routed two-way coupled model, and an alignment-routed three-way coupled model reproduce empirical predictive and topological landscape summaries.
- The selected panel is deliberately not a Tsuboyama-only stability panel. It includes `14` non-Tsuboyama systems spanning fluorescent/activity, binding, expression, and organismal-fitness readouts, plus `10` capped Tsuboyama stability controls.
- The primary non-Tsuboyama systems are `PHOT_CHLRE_Chen_2023`, `GFP_AEQVI_Sarkisyan_2016`, three Somermeyer GFP-family assays, two `SPG1_STRSG` binding assays, two KRAS readouts, `GRB2_HUMAN_Faure_2021`, `HIS7_YEAST_Pokusaeva_2019`, two `F7YBW8_MESOW` assays, and `GCN4_YEAST_Staller_2018`.
- The production runner is `scripts/proteingym_multi_system_coupled_landscape_production.py`. It writes per-system resumable checkpoints and skips completed `system_summary.json` outputs on restart, so timeout recovery should resume rather than restart the full 24-system panel.
- The runner records held-out single/double/triple/high-order predictive metrics, graphlet distances, DFE summaries, local ruggedness/autocorrelation, random/adaptive walk epistasis summaries, WHT/GFT coefficients, coupled-particle summaries, and base-versus-coupled function correlation statistics.
- A local import check and target-selection preflight passed. The preflight selected exactly `24` systems with `10` Tsuboyama assays and confirmed `PHOT_CHLRE_Chen_2023` as the first anchor system.
- The first live submission, `RUN-100`, failed before execution because Slurm rejected the initial `96G` memory request as unschedulable on the homelab `all` partition. The experiment memory request was reduced to `86G`, the largest conservative request fitting the available `90,000 MB` node while preserving the `168:00:00` walltime.
- `RUN-102`, submitted to `lab-slurm` as Slurm job `249` on `2026-05-18`, failed after completing the PHOT base and coupled fits. The fatal error was a `KeyError: 'split'` in the reused GFT helper during walk-step summarisation because EXP-040 passed the original graphlet record rather than the graphlet payload containing the helper's required `split` field.
- The production runner has been fixed to pass the graphlet payload into walk summarisation and to seed a fresh run from `resume.source_run_id` so failed-run outputs are copied into the next run directory before processing.
- `RUN-104`, submitted to `lab-slurm` as Slurm job `250` on `2026-05-18`, successfully seeded from `RUN-102`, completed `PHOT_CHLRE_Chen_2023` with graphlet/GFT/walk outputs, and had progressed to `GFP_AEQVI_Sarkisyan_2016` before being cancelled because the empirical memory use did not justify the `86G` request.
- The current resumed production run is `RUN-106`, submitted to `lab-slurm` as Slurm job `251` on `2026-05-18` with `12` CPUs, `24G` memory, and a `7-00:00:00` limit. The config now seeds from `RUN-104`; logs confirm `RUN-106` copied `systems` and `system_configs`, skipped the completed PHOT system, and resumed at `GFP_AEQVI_Sarkisyan_2016`.

## 2026-06-15 EXP-041 exact PHOT descriptive failure gate

- `EXP-041` was completed under `HYP-004` and `HYP-009` to enforce the rule that ABC should only follow an exact descriptive analysis showing which model summaries fail.
- The new runner `scripts/phot_exact_descriptive_failure_analysis.py` evaluates complete PHOT graphlets from the existing `RUN-039` prepared assay in binary hypercube node order, without stochastic random-walk sampling or ABC resampling.
- The benchmark compares exact empirical `latent_phi` graphlet summaries against fixed `RUN-039` branches: SMC-ABC best particle, SMC-ABC posterior mean, baseline shared stability-readout, and predictive rich-pair stability-readout.
- Outputs include `exact_model_failures.csv`, `exact_model_distances.csv`, and `abc_descriptive_gate.csv`; only components exceeding the preregistered failure threshold are eligible targets for any later ABC refit.
- Regression coverage for exact graphlet ordering, completeness, additive no-epistasis behavior, deterministic walk summaries, and failure ranking now lives in `tests/test_dense_descriptive.py`.
- `RUN-108` completed on `lab-slurm`. `ANA-040` copied only compact summary
  tables into the analysis record; the larger exact graphlet tables remain on
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-108/outputs`.
- `RES-040` records the result: all four fixed `RUN-039` branches fail the
  exact gate only for `edge_energy_top10_fraction`,
  `edge_energy_effective_fraction`, `site_energy_effective_fraction`, and
  `site_energy_top1_fraction`. Sign-epistasis, local-optima, basin, spectral,
  random-walk, and Moran-I summaries do not cross the severe failure threshold.

## 2026-06-16 EXP-042 PHOT nested-min cheap fitting equivalence benchmark

- `EXP-042` was created under `HYP-004` and `HYP-009` to test whether cheaper
  fitting protocols can recover the same operational PHOT nested min-coupled
  residual landscape behavior as the full `RUN-098` ABC comparator.
- The experiment keeps the latest sequential residual nested-min family fixed
  and compares deterministic summary-loss search, surrogate-assisted fitting,
  exact-descriptive-gated fitting, and local ABC initialized around cheap fits.
- Primary equivalence uses the original `RUN-098` held-out graphlet summary
  distance vector so the comparison is apples-to-apples. Additional PHOT
  properties, including sign/reciprocal/magnitude epistasis and diminishing
  returns, are reported as secondary descriptive readouts with the selected
  `RUN-098` full-ABC particles re-evaluated through the same reporting path.
- `RUN-110` completed on `lab-slurm` as Slurm job `1906` with exit code `0`.
  `ANA-041` copied only the compact source tables needed for interpretation;
  the larger WHT/GFT coefficient tables remain on
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-110/outputs`.
- `RES-041` records the result: `surrogate_assisted` and
  `exact_descriptive_gated` pass same-family equivalence against the selected
  `RUN-098` full-ABC fitted-family candidates for all five nested families,
  while deterministic summary loss and local ABC initialized from cheap fits
  pass four of five. No cheap protocol reaches the best empirical-readout
  residual diagnostic comparator.
- The supported conclusion is narrow but useful for scale-up: cheaper protocols
  can recover selected full-ABC fitted-family PHOT graphlet behaviour under
  the preregistered same-family margins, but `RUN-110` does not prove general
  equivalence to the full ABC posterior or solve the empirical residual
  diagnostic mismatch.

## 2026-06-16 EXP-043 PHOT epistasis-calibrated cheap fitting protocol

- `EXP-043` was created under `HYP-004` and `HYP-009` to test whether the cheap
  PHOT nested-min protocol can resolve the `RES-041` epistasis-amplitude failure
  by fitting an explicit positive affine observable readout and promoting
  pairwise epistasis-amplitude summaries into the fitting objective.
- The runner now fits the readout on calibration graphlets, applies it before
  objective scoring and final reporting, and keeps raw latent epistasis in
  `raw_*` columns so the analysis can distinguish a readout-scale fix from a
  true latent-geometry fix.
- The exact descriptive gate still includes the four `RES-040` local edge/site
  Dirichlet-energy failures, now augmented with magnitude/signed pairwise
  epistasis and accessibility summaries.
- A local quick smoke test completed with the calibrated observable epistasis
  metrics on the empirical scale while preserving raw latent reporting columns.
  Focused regression tests for the benchmark runner passed.
- `RUN-111` records the dry-run submission. `RUN-112` completed on `lab-slurm`
  as Slurm job `1907` with inferred exit code `0`. Runtime outputs remain at
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-112/outputs`; `ANA-042`
  copied only compact summary/source tables and left large WHT/GFT tables on
  the runtime volume.
- `RES-042` records the result: the calibrated cheap protocol failed the
  preregistered equivalence decision (`0/40` operationally equivalent rows).
  Aggregate distance margins often passed, but every method failed component
  guardrails, dominated by signed pairwise epistasis. The best calibrated
  candidate (`exact_descriptive_gated` / `uncoupled_module`) had test distance
  `1.1474` but magnitude epistasis mean `0.0156` versus empirical `0.2926`,
  raw latent magnitude epistasis `29.10`, and accessibility `0.0` versus
  empirical `0.1975`. The closest row by calibrated epistasis amplitude still
  retained raw latent magnitude epistasis `131.9`. This weakens the current
  PHOT adequacy form of `HYP-004` and motivates better-conditioned matched
  summaries for `HYP-009`, but it does not refute the broader benchmark
  construction aim.

## 2026-06-16 EXP-044 PHOT bounded trait-space current-model fast fit

- `EXP-044` was created under `HYP-004` and `HYP-009` to test a minimal
  model-level correction to the PHOT epistasis-amplitude failure: route
  epistatic signals through bounded trait space before the readout rather than
  relying on unbounded score/residual-min curvature.
- The core AdaptEnv patch adds optional `trait_space_epistasis_bound` and
  `trait_space_epistasis_bound_mode` fields. Defaults preserve historical
  behavior; when enabled, non-score epistatic perturbations are transformed by
  a symmetric `tanh` or clip bound before they perturb stability/function.
- The new runner `scripts/phot_trait_space_bounded_gate_fast_fit.py` rebuilds
  actual current-model PHOT landscapes from the frozen `RUN-077` source,
  performs deterministic/random numerical fast fitting over stability,
  epistasis, BLOSUM, and bounded trait-space parameters, and compares selected
  candidates against empirical PHOT graphlets plus the `RUN-098` full-ABC
  comparator distance table.
- A local quick smoke completed successfully under `scratch/exp044_quick`,
  writing the expected compact outputs and readable empty-table headers.
  Focused regression tests passed:
  `.venv/bin/pytest tests/test_calibration_modes.py tests/test_phot_trait_space_bounded_gate_fast_fit.py -q`.
- `RUN-113` is the dry-run record. The first live staging attempt, `RUN-114`,
  was interrupted before Slurm submission to avoid sending unnecessary local
  cache/data payloads and was marked failed. `.gitignore` now excludes local
  `.dvc/cache`, `.dvc/tmp`, `.pytest_cache`, and `scratch/` directories.
- `RUN-115` is the live submission on `lab-slurm` as Slurm job `1909` with
  `8` CPUs, `32G` memory, and a `12:00:00` walltime. Runtime outputs are under
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs`.
- `RUN-115` completed with inferred Slurm exit code `0`. `ANA-043` copied only
  compact source/summary tables into the analysis record and left the larger
  WHT/GFT coefficient tables on the runtime volume.
- `RES-043` records the result: `0/12` selected bounded trait-space candidates
  were operationally equivalent to the best available `RUN-098` full-ABC PHOT
  comparator. All selected bounded candidates passed component-distance
  guardrails, which is an improvement over `RES-042`, but every row failed the
  absolute and relative aggregate-distance margins. The best held-out bounded
  candidate had test distance `0.4269` versus the best full-ABC context row
  `0.2115`, observable magnitude epistasis `0.1708` versus empirical `0.2926`,
  raw latent magnitude epistasis `103.65`, and accessibility `0.0120` versus
  empirical `0.1975`.
- The current scale-up conclusion is therefore negative: do not treat the
  `EXP-044` bounded fast-fit protocol as a validated full-ABC replacement for
  GB1/TrpB. Keep bounded trait-space routing as a promising model-level
  direction, but the remaining problem is latent geometry calibration rather
  than only optimizer cost.

## 2026-06-17 EXP-045 PHOT cumulative epistasis-calibration fix ablation

- `EXP-045` was created under `HYP-004` and `HYP-009` to follow
  `scratch/RES-043_epistasis_calibration_review.md` as a controlled PHOT
  diagnostic. It uses the same frozen PHOT source (`RUN-077`), the same
  `RUN-098` full-ABC comparator table, and the numerically efficient
  deterministic/random fast-fit runner from recent PHOT experiments.
- The experiment produces four result sets in one run: `A`,
  `A+B`, `A+B+C`, and `A+B+C+D`. These cumulatively add a monotone-piecewise
  observation link, a global sigmoid stability-margin readout with constrained
  specific epistasis, a fitted specific-epistasis/additive-scale ratio knob,
  and objective conditioning based on magnitude-to-single-mutant-SD ratio plus
  signed-epistasis distribution summaries instead of the near-zero signed mean.
- A quick smoke run completed under `scratch/exp045_quick/outputs` with all four
  stages represented in `summary.json` and selected-candidate tables. Focused
  regression tests passed:
  `.venv/bin/pytest tests/test_calibration_modes.py tests/test_phot_trait_space_bounded_gate_fast_fit.py tests/test_phot_nested_fit_protocol_benchmark.py -q`.
- `RUN-116` is the dry-run submission record. `RUN-117` is the live submission
  on `lab-slurm` as Slurm job `1911`, requesting `8` CPUs, `32G` memory, and
  `12:00:00` on the `all` partition. `RUN-117` completed with inferred exit
  code `0`; runtime outputs remain under
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-117/outputs`.
- `ANA-044` copied compact RUN-117 source tables and derived stage-level
  equivalence, component, raw/observable epistasis, and collection-audit
  summaries. Large WHT/GFT coefficient tables remain on the runtime volume.
- `RES-044` records the result: EXP-045 failed full-ABC equivalence with
  `0/32` operationally equivalent rows. Stage `A` and full `A+B+C+D` passed
  component guardrails but failed aggregate margins; stages `A+B` and
  `A+B+C` collapsed to zero-variation landscapes. The best full-stage row had
  test distance `0.4122`, observable magnitude epistasis `0.2426`, and
  accessibility `0.0190`.

## 2026-06-17 EXP-046 PHOT raw-latent soft-min pathology fix ablation

- `EXP-046` was created under `HYP-004` and `HYP-009` to follow
  `docs/RES-044_epistasis_calibration_followup.md` as a controlled PHOT
  diagnostic. It uses the same frozen PHOT source (`RUN-077`), the same
  `RUN-098` full-ABC comparator table, and the same deterministic/random
  numerical fast-fit runner used by recent PHOT experiments.
- The experiment produces three result sets in one run: `1`, `1+2`, and
  `1+2+3`. These cumulatively add a low-dimensional two-additive-trait
  `soft_min_values` latent generator with random specific epistasis disabled,
  raw-latent objective constraints on raw single-mutant SD, raw local-maxima
  count, and raw accessibility plus a restricted monotone-sigmoid observation
  readout, then up-weight raw accessibility as a first-class calibration
  target.
- The runner now reports raw accessibility and raw local-maxima columns in the
  selected-candidate graphlet tables and records `raw_latent_objective_*`
  components in `bounded_gate_candidates.csv`, so the analysis can separate
  observable readout fit from latent-geometry repair.
- A quick smoke run completed under `scratch/exp046_quick/outputs` with all
  three stages represented in `summary.json` and selected-candidate tables.
  Focused regression tests passed:
  `.venv/bin/pytest tests/test_calibration_modes.py tests/test_phot_trait_space_bounded_gate_fast_fit.py tests/test_phot_nested_fit_protocol_benchmark.py -q`.
- `RUN-118` is the dry-run submission record. `RUN-119` is the live submission
  on `lab-slurm` as Slurm job `1917`, requesting `8` CPUs, `32G` memory, and
  `12:00:00` on the `all` partition. `labproj refresh --write --strict`
  updated `RUN-119` to `completed`.
- `ANA-045` copied compact RUN-119 source tables and derived stage-level
  equivalence, component, raw/observable epistasis, raw-latent pathology, and
  collection-audit summaries. Large WHT/GFT coefficient tables remain on the
  runtime volume.
- `RES-045` records the result: EXP-046 failed PHOT calibration with `0/24`
  operationally equivalent rows. Stage `1` was the least bad row (`1.2240`
  held-out distance) but retained raw scale hundreds-fold too large; stages
  `1+2` and `1+2+3` had `1.6872` held-out distance, no component guardrail
  passes, and collapsed observable variation.

## 2026-06-17 EXP-047 PHOT RES-045 calibration isolation cumulative ablation

- `EXP-047` was created under `HYP-004` and `HYP-009` to follow
  `scratch/RES-045_epistasis_calibration_isolation_plan.md` as a controlled
  PHOT-only diagnostic. It keeps the same frozen PHOT source (`RUN-077`), the
  same `RUN-098` full-ABC comparator table, and the recent deterministic/random
  fast-fit reporting path.
- The experiment runs cumulative stages in the plan's conditioning-first order:
  `A` anchors the soft-min generator scale by z-scoring output before a single
  output-scale knob; `A+F` decouples scale-invariant geometry from the raw
  single-mutant scale term; `A+F+B` records a model-free k-trait additive
  soft-min capability probe and skips duplicate sequence-conditioned fitting;
  `A+F+B+C` sweeps k-trait dimensionality and auxiliary texture;
  `A+F+B+C+D` adds a small sparse structured residual and a restricted
  monotone-sigmoid readout; `A+F+B+C+D+E` increases random and local-refinement
  search.
- The runner `scripts/phot_trait_space_bounded_gate_fast_fit.py` now supports
  generic k-trait soft-min latent generation, optional generator-scale
  normalization, sparse structured residual texture, seeded initial candidates,
  probe-only cumulative stages, and `softmin_capability_probe.csv` output.
- Local validation before submission passed:
  `.venv/bin/pytest tests/test_phot_trait_space_bounded_gate_fast_fit.py -q`
  and `LABPROJ_QUICK=1 ... run.sh`, which wrote all expected EXP-047 quick
  outputs under `scratch/exp047_quick/outputs`.
- `RUN-120` is the dry-run submission record. `RUN-121` is the live submission
  on `lab-slurm` as Slurm job `1918`, requesting `8` CPUs, `32G` memory, and
  `12:00:00` on the `all` partition. `RUN-121` completed with scheduler exit
  code `0`; runtime outputs remain under
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-121/outputs`.
- `ANA-046` copied compact RUN-121 source tables and derived stage-level
  equivalence, component, raw/observable epistasis, capability-probe,
  prior-run-comparison, and collection-audit summaries. Large WHT/GFT
  coefficient tables remain on the runtime volume.
- `RES-046` records the result: EXP-047 improved over the `RES-045`
  raw-latent soft-min failure but did not solve PHOT calibration. It produced
  `0/32` operationally equivalent rows. The best held-out stage
  `A+F+B+C+D` reached distance `0.9155`; raw scale was near empirical, but
  accessibility remained weak and observable single-mutant/epistatic variation
  collapsed to near-zero.

## 2026-06-17 EXP-048 PHOT sparse Walsh Potts spectrum-calibrated generator probe

- `EXP-048` was created under `HYP-004` and `HYP-009` to follow the Generator
  A plan in `scratch/RES-046_scalable_multipeak_generator_design.md`: an
  additive stability core plus sparse hashed high-order Potts/Walsh terms,
  with structured graphlet-aligned support clusters and empirical spectrum
  calibration before any sequence-conditioned fit.
- The `external/Adapt-Env` landscape implementation now has an opt-in sparse
  Walsh/Potts epistatic field. Defaults preserve historical behavior; when
  enabled, deterministic hashed supports and per-site mean-zero/unit-variance
  contrasts generate a generic epistatic signal that can be routed through the
  existing trait-space/readout path.
- The PHOT fast-fit runner now implements the section-6 Generator A protocol:
  `G` estimates the empirical graphlet Walsh variance-by-order and bootstrap
  noise floor, maps local order weights to global sparse-field order weights
  by Monte Carlo over the support pool, `H1` runs the forward sparse-Walsh
  capability probe, and the sequence-conditioned numerical fast fit proceeds
  only after the capability gate in the full run.
- The experiment explicitly excludes the FGM/H2 branch from this run. That is
  a scope decision, not a scientific conclusion about FGM.
- A local quick execution under `scratch/exp048_quick/outputs` wrote the new
  `sparse_walsh_spectrum_calibration.csv` and
  `sparse_walsh_forward_probe.csv` tables plus the downstream fit/equivalence
  tables. The quick run is only an execution smoke because it uses two
  graphlets and failed the tiny forward probe.
- Focused regression tests passed:
  `.venv/bin/pytest tests/test_calibration_modes.py tests/test_phot_trait_space_bounded_gate_fast_fit.py -q`.
- `RUN-122` is the dry-run submission record. `RUN-123` is the live submission
  on `lab-slurm` as Slurm job `1930`, requesting `8` CPUs, `32G` memory, and
  `12:00:00` on the `all` partition. `RUN-123` completed with inferred
  scheduler exit code `0`; runtime outputs remain under
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-123/outputs`.
- `ANA-047` copied compact RUN-123 source tables and derived the
  spectrum-calibration, H1 forward-probe, stage-gate, comparator-context, and
  output-audit summaries. The PHOT alignment FASTA remains on the runtime
  volume because it is not needed for the gate decision.
- `RES-047` records the result: EXP-048 did not validate Generator A. The
  H1 forward gate produced `0/12` joint passes; accessibility, magnitude
  epistasis relative to single-mutant scale, and Walsh high-order fraction
  failed for every probe candidate, so no sequence-conditioned PHOT
  equivalence fit was attempted.

## 2026-06-18 EXP-049 PHOT miscalibration P0 target-operator isolation

- `EXP-049` was created under `HYP-004` and `HYP-009` as the direct Phase-0
  follow-up to `RES-047` and
  `scratch/RES-047_miscalibration_root_cause_and_isolation_plan.md`.
- The runner `scripts/phot_miscalibration_p0_isolation.py` implements all
  three P0 checks in one experiment: P0.1 target bootstrap and MAVE-NN refit
  stability, P0.2 known-spectrum synthetic assay round-trip through MAVE-NN
  and dense complete graphlets, and P0.3 sparse-Walsh support containment
  against the sampled `dense_descriptive` graphlets.
- P0.1 is configured for the scratch-note requirement of at least five
  MAVE-NN seeds and at least two regularisation profiles. P0.2 uses pure
  order-3, mixed order-2/order-3, and low-order control injected spectra. P0.3
  uses the EXP-048 sparse-Walsh support source and the ANA-047 module catalog.
- Local validation passed for syntax, focused tests
  (`tests/test_phot_miscalibration_p0_isolation.py`), and a quick-mode smoke
  run to `/tmp/exp049-quick.aXIQqM` that wrote P0.1/P0.2/P0.3 outputs.
- `RUN-124` is the dry-run submission record. `RUN-125` is the live submission
  on `lab-slurm` as Slurm job `1933`, requesting `8` CPUs, `32G` memory, and
  `12:00:00` on the `all` partition. It completed with inferred scheduler
  exit code `0` and wrote the expected output set under
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-125/outputs`.
  The full run used `24` graphlets, `10` P0.1 MAVE-NN refits, `36` P0.2
  transfer rows, and `2048` P0.3 support terms. Stderr showed TensorFlow CPU
  startup and memory-allocation warnings but no fatal error.
- `ANA-048` copied the compact RUN-125 output bundle into the scientific
  record and derived decision tables for P0.1 target stability, P0.2
  known-spectrum round-trip transfer, and P0.3 support containment.
- `RES-048` records the result: the PHOT measurement path smooths/suppresses
  injected high-order and local epistasis, while sparse-Walsh support
  containment is adequate under the audited support rule. This weakens the
  current PHOT adequacy route under `HYP-004` and motivates better-conditioned
  target calibration for `HYP-009`; it does not support a generator repair or
  benchmark-ready matched landscape.

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
  - the first execution, `RUN-014`, was terminated and marked `failed` after it remained effectively single-process and made no durable progress through the SMC stage after roughly `20` hours on `lab-slurm`
  - the project-local Bayesian fitter was then refactored for HPC execution
  - `src/project_adapt_env/smc_abc.py` now supports batch-parallel Ray execution, resumable JSON checkpoints, and batch-level progress callbacks
  - `src/project_adapt_env/proteingym_panel.py` now caches per-assay `mavenn` preparation artifacts under the run output so interrupted runs do not retrain every assay from scratch
  - `scripts/proteingym_bayesian_summary_calibration.py` now writes atomic outputs, emits live progress logs and ETA estimates, restores cached empirical targets, and resumes synthetic-truth recovery incrementally
  - the `EXP-005` overlay environment now includes `ray==2.55.1`, and the run entrypoint exports unbuffered Python with single-threaded BLAS defaults for multi-worker CPU use
  - `RUN-016` then failed at the start of the Ray SMC stage because worker actors could not import `project_adapt_env` from the staged source tree
  - the Ray executor now propagates the staged `src/` path into worker `PYTHONPATH`, and the Slurm entrypoint now exports that path explicitly before launch
  - local code validation now includes `python -m compileall`, `pytest tests/test_smc_abc.py` with checkpoint-resume and Ray-worker import coverage (`4 passed`), and a Ray-backed SMC smoke test on a toy simulator
  - a full local `EXP-005 --quick` workflow smoke remains blocked in the current workstation interpreter because `mavenn` and `tensorflow` are not installed locally, but those dependencies are recorded in the Slurm overlay
  - `RUN-017` later consumed its full `24:00:00` walltime and was terminated by Slurm at `2026-05-11 09:11:50 AEST` while still in empirical SMC round `3`
  - the final durable `RUN-017` checkpoint captured completed empirical rounds `0`, `1`, and `2`, plus `344/512` proposal-pool attempts in round `3`
  - `scripts/proteingym_bayesian_summary_calibration.py` now supports seeding a new run from a prior run's cached `checkpoints/` directory via a `resume.source_run_id` config entry
  - `EXP-005` is now configured for a resumed `12`-worker Ray backend on the `gpu` Slurm partition so the job runs on `ripley` rather than `dutch`
  - the dry-run record for that resumed submission is `RUN-018`
  - the resumed execution was `RUN-019`, which was submitted to `lab-slurm` as scheduler job `67`, restored the `RUN-017` checkpoint/cache into its own output directory, and completed successfully on `2026-05-12`
  - the promoted durable outputs are now stored under `data/processed/proteingym-shared-summary-smc-abc-calibration-panel/RUN-019`
  - `RUN-019` completed the final empirical SMC round and both preregistered synthetic-truth recovery branches on the six-assay paired panel; the measurement layer remained strong with mean `mavenn` test Spearman `0.890`
  - the Bayesian fit no longer collapsed to zero epistasis: the best particle fit `epistasis_strength = 0.053` and `empirical_pairwise_strength = 0.0008`, while the posterior mean retained `epistasis_strength = 0.043`
  - synthetic-truth recovery was strong on the matched panel scaffold: the `moderate_epistatic` truth lay within the posterior q90 interval for `10/10` parameters, and the `flatter_low_epistasis` truth did so for `9/10`, missing only `noise_amplitude`
  - empirical recovery nevertheless remained weak. The best Bayesian shared fit reached only single-mutant holdout Spearman `0.195`, double-mutant holdout Spearman `0.119`, and functional KS `0.483`; it improved epistasis-prediction Spearman relative to the best deterministic control (`0.265` versus `0.094`) and avoided the near-peak reference artifact, but it did not rescue the shared empirical reconstruction claim
  - `ANA-005` records the scientific review of `RUN-019`, and `RES-005` records that the result still weakens `HYP-001` while shifting the main failure away from pure optimizer collapse and toward shared-model mismatch with the empirical panel
  - `EXP-006` has now been created under `HYP-001` to isolate whether the remaining `RES-005` failure is driven mainly by pooling across assays rather than by mismatch between the model family and an individual empirical target
  - it reuses the project-local Bayesian SMC-ABC workflow from `EXP-005` but restricts the empirical target to one assay, `SPTN1_CHICK_Tsuboyama_2023_1TUD`, chosen because it was the strongest prior per-assay candidate in the earlier deterministic records
  - the experiment keeps the same `12`-worker Ray backend and `48:00:00` walltime request on the `gpu` Slurm partition, but reduces requested memory to `24G` to avoid over-reserving queue capacity for the single-assay diagnostic
  - the dry-run record for this submission is `RUN-020`
  - the first live execution was `RUN-021`, submitted to `lab-slurm` as scheduler job `78` on `2026-05-12`
  - `RUN-021` failed during empirical SMC round `3` after `468/512` proposal-pool evaluations when the Ray backend raised `ActorUnavailableError` with a keepalive watchdog timeout; the partial checkpoint and assay-preparation cache remain available under `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-021`
  - the project-local SMC executor now retries recoverable Ray actor failures by rebuilding the actor pool and rerunning the same parameter batch instead of aborting the whole run
  - the resumed execution was `RUN-022`, submitted to `lab-slurm` as scheduler job `80` on `2026-05-12`; it restored the `RUN-021` output cache via `resume.source_run_id` and completed successfully
  - the promoted durable outputs are now stored under `data/processed/proteingym-single-assay-smc-abc-calibration-sptn1-chick/RUN-022`
  - `RUN-022` kept the measurement layer strong on `SPTN1_CHICK_Tsuboyama_2023_1TUD` with `mavenn` test Spearman `0.907`
  - the best Bayesian single-assay fit improved over the earlier six-assay shared Bayesian fit on single-mutant holdout Spearman (`0.278` versus `0.195`), double-mutant holdout Spearman (`0.362` versus `0.119`), and functional KS (`0.337` versus `0.483`), while retaining nonzero `epistasis_strength = 0.0497`
  - that single-assay improvement did not rescue the broader realism claim: the fit still did not dominate the strongest prior assay-specific deterministic `SPTN1` fit from `RUN-012`, and no branch achieved a clean reconstruction across ranking, epistasis prediction, and reference-to-peak behavior simultaneously
  - `ANA-006` records the scientific review of `RUN-022`, and `RES-006` records that the result still weakens `HYP-001` while narrowing the failure mode: cross-assay pooling matters, but it is not a sufficient explanation for the remaining empirical mismatch
  - `HYP-007` has now been recorded to capture a mechanistic alternative motivated by `RES-005`: the model may need multiple explicit latent molecular phenotypes, such as folding stability together with binding or abundance, rather than one generic functional layer plus sparse epistasis
  - this hypothesis is motivated by Otwinowski-style folding and binding decomposition and by the Lehner-group `deepPCA`, `ddPCA`, and `MoCHI` line of work, but it is not yet linked to a project experiment
  - the project-local `external/Adapt-Env` package now supports assay-class-specific calibration readouts through `synthetic_readout_mode = fitness | stability | stability_margin`, and empirical pairwise terms can now be targeted either to final score space or directly to the latent stability trait before the folding gate
  - this change is a minimal mechanistic test of the `RES-006` diagnosis: the current family may be expressive enough, but it likely places the assay nonlinearity and sparse empirical epistasis at the wrong stage of the causal model for Tsuboyama-style proteolysis data
  - regression coverage for the new behavior now lives in `tests/test_calibration_modes.py`, and local validation passed with `pytest tests/test_calibration_modes.py tests/test_smc_abc.py -q` (`7 passed`)
  - `EXP-007` has now been created under `HYP-001` to test that readout-stage misspecification directly on the same `SPTN1_CHICK_Tsuboyama_2023_1TUD` single-assay target used by `EXP-006`
  - it keeps the earlier `baseline_shared_raw` and `predictive_richpair_shared_raw` branches as paired controls, and adds `baseline_shared_stability_readout` plus `predictive_richpair_shared_stability_readout`, which calibrate against `stability_margin` and route empirical pairwise corrections into latent stability
  - because the calibration semantics changed, the experiment is configured as a fresh run rather than a resume from `RUN-021` or `RUN-022`
  - the dry-run record for this submission is `RUN-023`
  - the live execution was `RUN-024`, submitted to `lab-slurm` as scheduler job `81` on `2026-05-12`, and it completed successfully
  - the promoted durable outputs are now stored under `data/processed/proteingym-single-assay-stability-readout-smc-abc-sptn1-chick/RUN-024`
  - `RUN-024` kept the measurement layer strong on `SPTN1_CHICK_Tsuboyama_2023_1TUD` with `mavenn` test Spearman `0.908`
  - the new deterministic stability-readout controls materially improved recovery over the earlier raw-readout controls, but both returned to a near-peak reference solution (`fraction_of_peak` about `0.992`), which kept the peak-interpretation problem alive in the deterministic path
  - the best Bayesian stability-targeted fit provided the strongest balance so far on this assay: single-mutant holdout Spearman improved from `0.278` in `RUN-022` to `0.556`, double-mutant holdout Spearman improved from `0.362` to `0.483`, epistasis-prediction Spearman improved from `0.245` to `0.540`, and functional KS fell from `0.337` to `0.233`, while the fitted reference remained far from the peak (`fraction_of_peak 0.039`, distance `50`) and `epistasis_strength` stayed nonzero at `0.0686`
  - relative to the earlier strongest assay-specific deterministic `SPTN1` fit from `RUN-012`, the `RUN-024` Bayesian fit also established a stronger assay-level ceiling, improving single-mutant ranking, double-mutant ranking, epistasis prediction, and KS simultaneously
  - preregistered synthetic-truth recovery remained strong under the revised stability-targeted readout, with both registered truths falling within the posterior q90 interval for `10/10` parameters
  - `ANA-007` records the scientific review of `RUN-024`, and `RES-007` records that the result motivates `HYP-001`: readout-stage misspecification was a major contributor to the earlier failure, but the broader shared-regime multi-assay claim remains untested under the revised semantics
  - the project-local `external/Adapt-Env` package now also supports an explicit stability-plus-function assay readout through `synthetic_readout_mode = stability_function`, additional readout parameters `readout_stability_midpoint`, `readout_stability_slope`, and `readout_function_exponent`, and a new `empirical_pairwise_target = function` route that pushes empirical pairwise signal into the functional layer before the assay readout is formed
  - local validation for that change passed via `python -m py_compile` on the modified modules and direct component-level import checks in the workstation environment; a full local ProteinGym end-to-end smoke remains unavailable here because the current interpreter does not include the full `pandas` / `mavenn` runtime used on Slurm
  - `EXP-008` has now been created under both `HYP-001` and `HYP-007` to test whether that explicit biophysical-function readout improves empirical recovery on one multi-mutant functional ProteinGym assay, `PHOT_CHLRE_Chen_2023`
  - the assay was chosen as a single-system activity target with substantial epistasis signal (`2,122` single mutants and `165,407` multiple mutants in the ProteinGym reference), so the new latent readout can be tested without cross-assay pooling
  - `EXP-008` keeps the same single-assay SMC-ABC workflow and raw-readout control branches used in the recent `SPTN1` round, but adds paired `baseline_shared_biophysical_function_readout` and `predictive_richpair_shared_biophysical_function_readout` branches
  - the dry-run submission for `EXP-008` is `RUN-025`
  - the live execution was `RUN-026`, submitted to `lab-slurm` as scheduler job `95` on `2026-05-12`, and it completed successfully
  - the promoted durable outputs are now stored under `data/processed/proteingym-single-assay-biophysical-function-readout-smc-abc-phot-chlre/RUN-026`
  - `RUN-026` used the intended `PHOT_CHLRE_Chen_2023` functional activity assay, with `2,122` single mutants and `165,407` multiple mutants, so the single-system test had substantial epistasis coverage
  - the assay-specific `mavenn` model was usable but only moderate by recent project standards, with test Spearman `0.680` and test NRMSE `0.935`
  - the explicit stability-plus-function readout did not produce the predicted across-the-board gain over the paired raw-readout controls: the best Bayesian fit improved double-mutant holdout Spearman to `0.562`, but single-mutant holdout Spearman fell to `0.293`, epistasis-prediction Spearman fell to `0.109`, and functional KS worsened to `0.367` relative to the strongest raw-readout control within the same run
  - the deterministic biophysical-function branches also collapsed the fitted reference exactly onto the peak (`fraction_of_peak 1.000`, distance `0`), so the richer readout did not resolve the peak-interpretation problem in the deterministic path
  - synthetic-truth recovery remained strong at the row level, with both preregistered truths inside the posterior q90 interval for `13/13` fitted parameters, although the aggregate `summary.json` truth-count block underreports those counts and should not be treated as authoritative
  - `ANA-008` records the scientific review of `RUN-026`, and `RES-008` records that this first functional activity-assay test weakens `HYP-007`; the richer readout improved some metrics but did not outperform the raw-readout controls on the full empirical package
  - a literature-grounded review now suggests that `RES-008` is more consistent with causal misspecification of the chosen assay than with a clean falsification of the multi-latent direction itself
  - `PHOT_CHLRE` fluorescence likely mixes folded fraction, FMN occupancy, and photophysical brightness, so the `EXP-008` readout was both underidentified and structurally mismatched to the implemented two-latent stability-plus-function map
  - the literature points toward matched orthogonal abundance-plus-binding assays and explicit thermodynamic folded/unfolded/bound readouts as a better causal test, consistent with the Otwinowski GB1 decomposition, the Lehner-group `MoCHI` framework, and Weng et al.'s KRAS `ddPCA` analysis
  - `external/Adapt-Env` now also supports `synthetic_readout_mode = stability_binding`, which exposes a thermodynamic binding readout with additional parameters `readout_function_midpoint` and `readout_function_slope`
  - `external/Adapt-Env` now also has a backward-compatible latent-trait block system for benchmark generation itself: one scalar observed fitness can be composed from multiple named internal FGM-style trait blocks through configured `observed_fitness_terms`, and both empirical and generic epistatic signal can now be routed into score space, latent stability, or a named latent trait before collapse
  - `EXP-009` has now been created under `HYP-007` as that sharper test, using the matched KRAS ProteinGym assays `RASK_HUMAN_Weng_2022_abundance` (`3,066` single and `22,946` multiple mutants) and `RASK_HUMAN_Weng_2022_binding-DARPin_K55` (`3,084` single and `21,789` multiple mutants)
  - the new runner, `scripts/proteingym_paired_biophysical_calibration.py`, compares binding-only versus paired abundance-plus-binding fits, and raw collapsed readouts versus the new thermodynamic binding readout
  - local static validation passed via `python3 -m py_compile`; deeper local runtime smoke remains unavailable in the workstation interpreter because it lacks the `numpy` / `pandas` / `mavenn` stack used on Slurm
  - the dry-run submission for `EXP-009` is `RUN-027`
  - the live execution was `RUN-028`, submitted to `lab-slurm` as scheduler job `96` on `2026-05-13`, and it completed successfully
  - the promoted durable outputs are now stored under `data/processed/proteingym-paired-kras-abundance-binding-thermodynamic-diagnostic/RUN-028`
  - `RUN-028` confirmed that both KRAS readout layers were learnable by the assay model: abundance MAVE-NN test Spearman was `0.818`, and binding MAVE-NN test Spearman was `0.933`
  - the paired thermodynamic binding branch improved binding KS from `0.381` to `0.258` and epistasis-prediction Spearman from `0.104` to `0.287` relative to the paired raw-readout branch, but it also worsened abundance KS from `0.582` to `0.832`, single-mutant holdout Spearman from `0.265` to `0.105`, and double-mutant holdout Spearman from `0.252` to `0.144`
  - all `RUN-028` branches still placed the fitted reference essentially at the peak, so the thermodynamic KRAS paired-readout diagnostic does not resolve the peak-interpretation problem
  - `ANA-012` records the scientific review of `RUN-028`, and `RES-012` records that the mixed improvement pattern weakens the current `HYP-007` implementation rather than supporting a clean multi-latent rescue
  - `EXP-010` has now been created under `HYP-001` and `HYP-007` as a direct follow-up to `EXP-008`, returning to the same multi-mutant `PHOT_CHLRE_Chen_2023` activity assay but replacing the old special readout with the new generalized observed-fitness composition path in `Adapt-Env`
  - its base generator now makes public scalar fitness the product of `stability_gate` and latent `function_capacity`, and routes both generic and empirical epistatic signal through the latent function layer before collapse (`observed_fitness_combine_mode = product`, `observed_fitness_terms = [stability_gate, function_capacity]`, `generic_epistasis_target = function`)
  - because the current single-assay Bayesian runner supports one base generator per experiment, `EXP-010` treats `EXP-008` / `RES-008` as the historical baseline rather than pretending to compare incompatible generators inside one posterior fit
  - the new deterministic branches are `baseline_shared_activity_readout` and `predictive_richpair_shared_activity_readout`, and the Bayesian fit now calibrates against `synthetic_readout_mode = fitness` so the target is the scalar activity readout itself rather than a secondary latent proxy
  - the dry-run submission for `EXP-010` is `RUN-029`
  - the live execution was `RUN-030`, submitted to `lab-slurm` as scheduler job `135` on `2026-05-13`, and it completed successfully
  - the promoted durable outputs are now stored under `data/processed/proteingym-single-assay-latent-activity-observed-fitness-smc-abc-phot-chlre/RUN-030`
  - the assay-specific `mavenn` layer was slightly stronger than in `RUN-026`, with test Spearman `0.689` and test NRMSE `0.921`
  - the deterministic `RUN-030` activity-readout branches were effectively identical to the deterministic biophysical-function branches from `RUN-026`, implying that the generalized observed-fitness composition path largely reproduced the earlier deterministic activity-readout behavior on `PHOT_CHLRE`
  - the best Bayesian `RUN-030` fit improved functional KS relative to `RUN-026` (`0.253` versus `0.367`) but did not rescue the broader empirical package: single-mutant holdout Spearman remained essentially unchanged (`0.294` versus `0.293`), double-mutant holdout Spearman fell sharply (`0.363` versus `0.562`), epistasis-prediction Spearman turned negative (`-0.026` versus `0.109`), and the fitted reference moved back to a near-peak state (`fraction_of_peak 0.9995`, distance `5`) rather than the earlier non-pathological location (`0.083`, distance `57`)
  - row-level synthetic-truth recovery remained strong, with both preregistered truths inside the posterior q90 interval for `10/10` fitted parameters
  - `ANA-009` records the scientific review of `RUN-030`, and `RES-009` records that this public-activity observed-fitness reformulation still does not rescue `PHOT_CHLRE` recovery, further weakening `HYP-007` without directly updating `HYP-001`
  - `EXP-011` has now been created under `HYP-001` and `HYP-007` as a direct control for the negative `PHOT_CHLRE` activity results
  - it keeps the same `PHOT_CHLRE_Chen_2023` assay and single-assay SMC-ABC workflow used in `EXP-008` and `EXP-010`, but reverts the generator and calibration semantics to the successful `RES-007` stability-targeted single-latent formulation (`synthetic_readout_mode = stability_margin`, `empirical_pairwise_target = stability`)
  - the control question is whether the negative `PHOT_CHLRE` results mainly reflect over-parameterization or misrouting in the newer activity-oriented formulations, versus a more fundamental mismatch between the assay and a stability-only latent causal chain
  - the dry-run submission for `EXP-011` is `RUN-031`
  - the first live submission was `RUN-032`, submitted to `lab-slurm` as scheduler job `158` on `2026-05-13`, but it failed immediately because `EXP-011` still had the generic scaffold `run.sh` placeholder rather than the project-specific Bayesian calibration launcher
  - `EXP-011/run.sh` has now been replaced with the same `proteingym_bayesian_summary_calibration.py` launcher used in the recent single-assay rounds
  - the corrected live execution was `RUN-033`, resubmitted to `lab-slurm` as scheduler job `159` on `2026-05-13`
  - `RUN-033` no longer counts as a missing-launch failure: it completed panel preparation and all four deterministic control branches, then reached `300/384` particles in empirical SMC round `0` before the process exited without a model-level traceback
  - because the paired `EXP-012` PHOT run later showed the same Ray-backed failure pattern under the `12`-worker, `24` GB execution profile and completed after reducing to `6` workers with `48` GB, `EXP-011` has now been resubmitted with that same execution-only change while leaving the scientific model and control question unchanged
  - the refreshed dry-run/live resubmission is `RUN-037` / `RUN-038`
  - `RUN-038` was submitted to `lab-slurm` as scheduler job `163` on `2026-05-13`, and it completed successfully
  - the promoted durable outputs are now stored under `data/processed/proteingym-single-assay-stability-readout-control-smc-abc-phot-chlre/RUN-038`
  - the simpler stability-readout control materially improved the recent PHOT Bayesian activity formulations on single-mutant holdout Spearman (`0.505`), double-mutant holdout Spearman (`0.721`), and functional KS (`0.166`), making it the strongest Bayesian PHOT fit so far on those metrics
  - the run still did not provide a clean rescue of the PHOT activity assay: epistasis-prediction Spearman remained weak (`0.066`), and the best fit still sat too near the peak (`fraction_of_peak 0.867`, distance `4`)
  - compared with the earlier successful stability-targeted `SPTN1` unfolding run from `RES-007`, the same stability-targeted formulation does not transfer cleanly to this activity readout even though some coarse predictive metrics improved
  - row-level synthetic-truth recovery remained strong, with both preregistered truths inside the posterior q90 interval for `10/10` fitted parameters
  - `ANA-011` records the scientific review of `RUN-038`, and `RES-011` records that this direct control weakens `HYP-007` further by showing that the PHOT negatives are not just an over-parameterization artifact, while still stopping short of a clean activity-assay rescue
  - subsequent code review of the negative PHOT activity runs found that `EXP-008` and `EXP-010` did not yet fit an explicit second latent trait block; both changed readout semantics, but neither configured `latent_trait_blocks`, so the generator still fell back to the legacy single functional block
  - `ANA-009`, `RES-008`, and `RES-009` have therefore been amended to narrow their interpretation: they weaken the specific PHOT implementations that were actually run, not the explicit two-latent-trait direction in its intended form, and they now record the additional confounders that the old runs still used competing epistasis levers on the same legacy functional block and optimized the older summary-vector SMC objective
  - the calibration code has now been tightened so legacy functional parameter updates (`functional_sigma_base`, `n_functional_dims`, and related fields) also update the primary explicit latent trait block when `latent_trait_blocks` are configured, preventing the next latent-trait experiment from fitting parameters that the generator ignores
  - the project-local Bayesian runner now also supports `abc.distance_mode = validation_objective`, which uses the same held-out predictive objective terms already used in deterministic calibration instead of only a bootstrap summary-vector distance
  - `EXP-012` has now been created under `HYP-001` and `HYP-007` as the explicit replacement test on the same `PHOT_CHLRE_Chen_2023` assay
  - it uses exactly two latent trait layers in practice: the built-in stability layer plus one explicit named `readout` trait block
  - the public scalar fitness is formed as `stability_gate * trait:readout:capacity`
  - the experiment removes the earlier competing generic epistasis lever by fixing `epistasis_strength = 0` and fitting only `empirical_pairwise_strength`, routed to `trait:readout`
  - it also drops the earlier free FGM-complexity search by fixing the explicit readout block to `1` latent dimension and fitting only its width through `functional_sigma_base`, which now aliases to the primary explicit trait block
  - the Bayesian SMC fit for `EXP-012` is configured to optimize the held-out validation objective directly, and preregistered synthetic-truth recovery is disabled for this run because that scaffold remains implemented only for the older summary-vector distance mode
  - the dry-run submission for `EXP-012` is `RUN-034`
  - the first live submission, `RUN-035`, failed during the Ray-backed SMC phase after deterministic controls completed
  - node-local Ray logs on `ripley` recorded repeated `RAY_WORKER_FAILURE` events followed by a dead raylet / GCS-heartbeat chain, which is consistent with a resource-level Ray runtime failure under the original `12`-worker, `24` GB execution profile rather than a scientific-model exception in the Python fit code
  - `EXP-012` has therefore been resubmitted with execution-only changes: the Ray backend was reduced to `6` workers with batch size `6`, and the job memory request was increased to `48` GB while leaving the scientific model, latent-trait structure, and validation objective unchanged
  - the new dry-run/live resubmission uses `RUN-036`, submitted to `lab-slurm` as scheduler job `162` on `2026-05-13`, and it completed successfully
  - the promoted durable outputs are now stored under `data/processed/proteingym-single-assay-explicit-two-latent-trait-validation-smc-abc-phot-chlre/RUN-036`
  - `RUN-036` was the first clean PHOT activity test that actually instantiated the intended explicit two-trait setup: built-in stability plus one named `readout` latent block, public scalar fitness formed as `stability_gate * trait:readout:capacity`, only `empirical_pairwise_strength` routed into `trait:readout`, and the Bayesian fit optimized directly on the held-out validation objective rather than the earlier summary-vector distance
  - the assay-specific `mavenn` layer remained usable but only moderate, with test Spearman `0.686` and test NRMSE `0.865`
  - the clean two-trait fit improved epistasis-prediction Spearman relative to the earlier PHOT Bayesian runs (`0.259` versus `0.109` in `RUN-026` and `-0.026` in `RUN-030`), but it still did not rescue the broader predictive package: single-mutant holdout Spearman fell to `0.198`, double-mutant holdout Spearman to `0.251`, functional KS worsened to `0.496`, and the fitted reference remained near the peak (`fraction_of_peak 0.986`, distance `6`)
  - the posterior pushed the only fitted epistatic lever, `empirical_pairwise_strength`, to the top of its prior range without finding a strong solution, which is more consistent with remaining structural mismatch than with underuse of epistatic capacity
  - `ANA-010` records the scientific review of `RUN-036`, and `RES-010` records that this cleaner explicit two-latent-trait validation-objective formulation still does not rescue `PHOT_CHLRE` recovery, weakening `HYP-007` more directly on this assay while not directly updating `HYP-001`
  - `EXP-013` through `EXP-017` have now been created under `HYP-007` to isolate the remaining PHOT structural mismatch after `RES-010` and `RES-011`
  - `EXP-013` reruns the PHOT stability-readout control under the same held-out validation objective used by `EXP-012`, so the stability-control advantage is not confounded with objective mismatch
  - `EXP-014` records a paired explicit two-trait validation-objective replicate/control of `EXP-012`
  - `EXP-015` tests whether a flexible monotone two-latent stability/activity readout surface beats the forced `stability_gate * readout_capacity` collapse
  - `EXP-016` tests a stability-first residual-trait diagnostic with sparse pair residuals, explicitly separating a stability component from additional PHOT readout structure
  - `EXP-017` is the oracle global-epistasis baseline: an additive latent score plus a flexible monotone MAVE-NN/MoCHI-style nonlinearity, intended to distinguish restrictive AdaptEnv latent priors from measurement noise or scalar-assay non-identifiability
  - the lightweight PHOT diagnostic runner for `EXP-015` to `EXP-017` is `scripts/proteingym_phot_structural_mismatch_diagnostic.py`; static validation passed with `.venv/bin/python -m py_compile`, while full local runtime smoke remains unavailable because the workstation environment lacks the Slurm overlay's full `pandas` / `mavenn` runtime
  - live Slurm submissions are now present in the scheduler: `RUN-039` for `EXP-013` ran as scheduler job `172` and completed successfully, and `RUN-040` for `EXP-014` ran as scheduler job `173` and completed successfully; later structural-mismatch runs are tracked in their own run metadata and remain outside the `RUN-040` analysis
  - the promoted durable outputs for `RUN-039` are stored under `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039`
  - `RUN-039` answered the `EXP-013` control question: the PHOT stability-targeted mapping retained and strengthened its advantage when fitted with the same held-out `validation_objective` family used by the explicit two-latent `RUN-036` experiment
  - the `RUN-039` Bayesian best particle reached single-mutant holdout Spearman `0.506`, double-mutant holdout Spearman `0.783`, epistasis-prediction Spearman `0.634`, functional KS `0.176`, and a non-pathological reference geometry (`fraction_of_peak 0.323`, distance `45`)
  - relative to the explicit two-latent validation-objective `RUN-036` best fit, `RUN-039` improved single-mutant ranking by `0.308`, double-mutant ranking by `0.531`, epistasis prediction by `0.376`, functional KS by `0.320` in the favorable direction, and lowered reference fraction of peak by `0.663`
  - relative to the earlier summary-vector stability control from `RUN-038`, `RUN-039` preserved the single-mutant advantage, improved double-mutant holdout ranking by `0.062`, improved epistasis-prediction Spearman by `0.568`, and greatly improved reference geometry, with only a small KS regression (`0.166` to `0.176`)
  - `ANA-013` records the scientific review of `RUN-039`, and `RES-013` records that this weakens the current PHOT `HYP-007` implementation by showing that the poor explicit two-latent result is not explained by the validation objective itself
  - the promoted durable outputs for `RUN-040` are stored and DVC-tracked under `data/processed/proteingym-phot-explicit-two-trait-validation-objective-replicate/RUN-040`
  - `RUN-040` answered the `EXP-014` replicate question: the explicit two-trait validation-objective PHOT model reproduced the `RUN-036` failure pattern almost exactly, making stochastic SMC noise a weak explanation for the earlier negative result
  - the `RUN-040` MAVE-NN observation layer reached test Spearman `0.677` and test NRMSE `0.878`
  - the `RUN-040` Bayesian best particle reached single-mutant holdout Spearman `0.198`, double-mutant holdout Spearman `0.252`, epistasis-prediction Spearman `0.264`, functional KS `0.498`, and near-peak reference geometry (`fraction_of_peak 0.986`, distance `6`)
  - relative to the earlier explicit two-trait validation-objective `RUN-036` best fit, `RUN-040` was essentially unchanged: single-mutant ranking delta `0.000`, double-mutant ranking delta `0.001`, epistasis-prediction delta `0.006`, functional KS delta `0.002`, and reference-distance delta `0`
  - relative to the `RUN-039` stability-targeted validation-objective control, `RUN-040` was substantially weaker: single-mutant ranking delta `-0.308`, double-mutant ranking delta `-0.531`, epistasis-prediction delta `-0.370`, functional KS delta `+0.322`, reference fraction-of-peak delta `+0.663`, and reference-distance delta `-39`
  - the `RUN-040` posterior again pushed `empirical_pairwise_strength` to the top of its configured prior (`mean 0.0777`, q05 `0.0768`, q95 `0.0797`, upper bound `0.08`) without finding a strong solution
  - `ANA-014` records the scientific review of `RUN-040`, and `RES-014` records that this replicate weakens the current PHOT `HYP-007` implementation by showing reproducible explicit two-trait structural mismatch under the same validation-objective family
  - the promoted durable outputs for `RUN-041` are stored and DVC-tracked under `data/processed/proteingym-phot-flexible-monotone-two-latent-readout-diagnostic/RUN-041`
  - `RUN-041` answered the `EXP-015` diagnostic question by replacing the forced product/gate collapse with a lightweight flexible monotone stability/activity readout surface fit on the same `PHOT_CHLRE_Chen_2023` assay
  - the `RUN-041` MAVE-NN observation layer was moderate, with test Spearman `0.657` and test NRMSE `0.947`
  - the flexible monotone two-latent diagnostic reached single-mutant holdout Spearman `0.473` and double-mutant holdout Spearman `0.807`, improving those ranking metrics over the explicit product/gate two-trait `RUN-040` result by `0.275` and `0.555`
  - the same diagnostic did not rescue epistasis: epistasis-prediction Spearman was only `0.035`, with epistasis KS `0.814`, which is worse than both `RUN-040` and the `RUN-039` stability-targeted control
  - because `RUN-041` is a score-prediction diagnostic rather than a full AdaptEnv landscape posterior, it cannot evaluate functional KS, reference-to-peak geometry, posterior identifiability, or synthetic-truth recovery
  - `ANA-015` records the scientific review of `RUN-041`, and `RES-015` records a mixed result that weakens `HYP-007` as a complete empirical recovery claim while preserving readout-collapse misspecification as a narrower contributor to the PHOT structural mismatch
  - the promoted durable outputs for `RUN-042` are stored and DVC-tracked under `data/processed/proteingym-phot-photophysical-trait-sparse-residual-diagnostic/RUN-042`
  - `RUN-042` answered the `EXP-016` diagnostic question by fitting a stability-first monotone component, fitting an additive photophysical residual trait on the residuals, and then adding sparse pair residual effects after the two-coordinate surface
  - the `RUN-042` MAVE-NN observation layer was moderate, with test Spearman `0.682` and test NRMSE `0.885`
  - the stability-first monotone model reached single-mutant holdout Spearman `0.514` and double-mutant holdout Spearman `0.555`, but epistasis-prediction Spearman was negative (`-0.198`) with poor KS (`0.907`)
  - adding the photophysical residual trait improved double-mutant holdout Spearman to `0.823`, slightly above `RUN-039` and `RUN-041`, but it slightly reduced single-mutant ranking to `0.497` and left epistasis Spearman essentially zero (`-0.001`) with poor KS (`0.860`)
  - the sparse pair-residual layer fitted `133` residual pair effects and used `1,721` photophysical mutation features, but it produced identical held-out single, double, and epistasis metrics to the residual-trait model, so the exported validation evidence does not support sparse residual epistasis as the missing PHOT structure
  - because `RUN-042` is a score-prediction diagnostic rather than a full AdaptEnv landscape posterior, it cannot evaluate functional KS, reference-to-peak geometry, posterior identifiability, or synthetic-truth recovery
  - `ANA-016` records the scientific review of `RUN-042`, and `RES-016` records a mixed result that weakens `HYP-007` as a complete empirical recovery claim while preserving stability-first residualization as a useful diagnostic for ranking signal
  - the promoted durable outputs for `RUN-043` are stored and DVC-tracked under `data/processed/proteingym-phot-oracle-global-epistasis-baseline/RUN-043`
  - `RUN-043` answered the `EXP-017` diagnostic question by fitting an oracle additive latent score to the observed PHOT assay and then applying a flexible monotone global-epistasis nonlinearity
  - the `RUN-043` MAVE-NN observation layer was moderate, with test Spearman `0.682` and test NRMSE `0.872`
  - the oracle additive-linear baseline reached single-mutant holdout Spearman `0.005`, double-mutant holdout Spearman `0.805`, epistasis-prediction Spearman `0.252`, and epistasis KS `0.907`
  - the oracle monotone global-epistasis baseline reached single-mutant holdout Spearman `0.005`, double-mutant holdout Spearman `0.803`, epistasis-prediction Spearman `0.240`, and epistasis KS `0.884`
  - relative to the `RUN-039` stability-targeted validation-objective control, the monotone oracle was much worse on single-mutant ranking (`-0.502`) and epistasis prediction (`-0.395`), only slightly better on double-mutant ranking (`+0.020`), and worse on epistasis KS (`+0.134`)
  - because `RUN-043` does not beat all mechanistic variants, it fails the preregistered condition for blaming restrictive AdaptEnv latent priors alone and points back toward scalar-assay measurement limits, insufficient identifiability, or missing stability-specific/sparse structure
  - `ANA-017` records the scientific review of `RUN-043`, and `RES-017` records a negative oracle-baseline result that weakens `HYP-007` as a complete empirical recovery or more-flexible-latent-prior explanation for the PHOT mismatch
  - `EXP-017` has been resubmitted as live `RUN-048` on `lab-slurm` with scheduler job `192`; at submission time it was queued pending priority.
  - the stale live `EXP-018` attempt `RUN-046` was reconciled to `failed`, and `EXP-018` has been resubmitted as live `RUN-051` on `lab-slurm` with scheduler job `194`; at submission time it was running on `ripley`.
  - `EXP-018` and `EXP-019` have now been created under `HYP-007` as paired non-KRAS single-mutant diagnostics for the multiple-latent-trait direction.
  - `EXP-018` tests `OXDA_RHOTO_Vanella_2023_expression` as the state readout against `OXDA_RHOTO_Vanella_2023_activity` as the function readout.
  - `EXP-019` tests `KCNJ2_MOUSE_Coyote-Maestas_2022_surface` as the state readout against `KCNJ2_MOUSE_Coyote-Maestas_2022_function` as the function readout.
  - both experiments use the new lightweight runner `scripts/proteingym_paired_trait_calibration.py`, comparing `function_only_single_trait`, `paired_state_function_single_trait`, and `paired_state_function_explicit_two_trait` branches.
  - both experiments are explicitly limited to single-mutant recovery and paired-readout structure; they do not test double-mutant ranking or epistasis recovery because both ProteinGym assay pairs lack multiple mutants.
  - the first live submissions, `RUN-046` for `EXP-018` and `RUN-047` for `EXP-019`, failed immediately because MAVE-NN requires a linear GE nonlinearity when only single mutants are present.
  - the experiment configs have been corrected to use `ge_nonlinearity_type: linear` for these single-mutant-only paired assays.
  - the corrected dry-run records are `RUN-049` for `EXP-018` and `RUN-050` for `EXP-019`.
  - the corrected live OXDA execution is `RUN-051`, submitted to `lab-slurm` as scheduler job `194`; a duplicate pending OXDA submission, `RUN-052` / scheduler job `195`, was canceled to avoid running the same experiment twice.
  - the corrected live KCNJ2 execution is `RUN-053`, submitted to `lab-slurm` as scheduler job `196`; at submission time it was queued pending priority.
  - `RUN-051` completed successfully, and its promoted durable outputs are now stored and DVC-tracked under `data/processed/proteingym-oxda-expression-activity-paired-latent-trait-diagnostic/RUN-051`
  - `RUN-051` answered the `EXP-018` OXDA paired-readout diagnostic question on matched `OXDA_RHOTO_Vanella_2023_expression` and `OXDA_RHOTO_Vanella_2023_activity` single-mutant assays
  - the OXDA expression/activity readouts were only weakly correlated (`Spearman 0.246` over `6387` matched variants), while the assay-specific MAVE-NN layers were moderate (`0.567` expression and `0.542` activity test Spearman)
  - the paired state/function branches did not improve activity recovery over the function-only branch: function-only activity holdout Spearman was `0.183`, paired single-trait was `0.061`, and paired explicit two-trait was `0.180`
  - activity functional KS also did not improve (`0.520` function-only, `0.523` paired single-trait, `0.526` explicit two-trait), and all branches retained near-peak reference geometry with function reference fraction of peak `0.9995` to `0.9999`
  - because both OXDA assays contain only single mutants, `RUN-051` does not test double-mutant ranking or epistasis recovery
  - `ANA-018` records the scientific review of `RUN-051`, and `RES-018` records that this OXDA paired state/function diagnostic weakens the current `HYP-007` implementation without directly updating `HYP-001`
  - `RUN-053` completed successfully, and its promoted durable outputs are now stored and DVC-tracked under `data/processed/proteingym-kcnj2-surface-function-paired-latent-trait-diagnostic/RUN-053`
  - `RUN-053` answered the `EXP-019` KCNJ2 paired-readout diagnostic question on matched `KCNJ2_MOUSE_Coyote-Maestas_2022_surface` and `KCNJ2_MOUSE_Coyote-Maestas_2022_function` single-mutant assays
  - the KCNJ2 surface/function readouts were weakly correlated (`Spearman 0.290` over `6789` matched variants), while the assay-specific MAVE-NN layers were asymmetric (`0.538` surface and `0.389` function test Spearman)
  - the paired state/function branches did not improve function recovery over the function-only branch as a package: function-only function holdout Spearman was `0.066`, paired single-trait was `0.065`, and paired explicit two-trait was `0.021`
  - function KS also did not improve (`0.460` function-only, `0.502` paired single-trait, `0.511` explicit two-trait), and all branches retained zero `epistasis_strength`, zero `empirical_pairwise_strength`, zero `noise_amplitude`, and near-peak reference geometry with function reference fraction of peak `0.9993`, `0.9998`, and `0.9932`
  - because both KCNJ2 assays contain only single mutants, `RUN-053` does not test double-mutant ranking or epistasis recovery
  - `ANA-019` records the scientific review of `RUN-053`, and `RES-019` records that this KCNJ2 paired state/function diagnostic weakens the current `HYP-007` implementation with a function-observation-layer caveat and without directly updating `HYP-001`
  - Based on the accumulated negative evidence for FGM-style functional latent traits, `external/Adapt-Env` now treats FGM components as legacy-only implementation paths rather than public biological benchmark controls.
  - The root package no longer exports the FGM trait-block config classes, public presets and docs no longer advertise FGM knobs, API/service state omits legacy FGM config fields, and runtime code logs a warning that FGM components are not empirically biological in this system and will be removed in a future release when those paths are instantiated or explicitly supplied.
  - Local validation of the AdaptEnv package after this interface change passed with `PYTHONPATH=. pytest tests -q` (`71 passed`, `6 skipped`) and `python3 -m compileall adaptenv`.
  - The `external/Adapt-Env` public landscape core has now been tightened further around the single thermodynamic stability-margin latent trait.
  - The default scalar readout is `stability_margin`, generic and empirical pairwise effects now target stability by default, deterministic internal roughness defaults to zero, and nonzero designed-peak distance is no longer a public/default control.
  - The old FGM scalar path, observed-fitness composition path, score-space epistasis, named trait routing, higher-order exact-match epistasis, deterministic hash roughness, strict peak locking, peak-search controls, and peak-distance design are retained only as legacy/backward-compatible implementation paths tagged with `# TODO: Remove in next version`.
  - Public presets, service state, and documentation no longer advertise those legacy controls; supplying them through the service or preset overrides logs a warning that the legacy AdaptEnv core is not empirically biological in this system and will be removed in a future release.
  - Local validation after the broader legacy quarantine passed with `PYTHONPATH=. pytest tests -q` (`71 passed`, `6 skipped`) and `python3 -m compileall adaptenv`.
  - `CLAIM-001` and `CLAIM-002` have been created as the first formal claims in the scientific record, superseding `HYP-001` and `HYP-007`.
  - `CLAIM-001` records the reduced interpretation of `HYP-001`: empirical-like local statistics are recoverable for a single isolated family or assay system at a time, but the record does not support a universal shared parameter regime across heterogeneous systems.
  - `CLAIM-002` records the reduced interpretation of `HYP-007`: empirical biophysical fitting is supported as a single thermodynamic stability or stability-margin latent trait, while additional fitted latent molecular traits have not improved the full empirical recovery package.
  - `HYP-004` has been reframed away from benchmarking optimizer or machine-learning strategy separation. It now tests whether random walks, adaptive walks, beneficial-effect availability, and epistasis change with closeness to a peak in a biologically plausible way.
  - `EXP-020` has been created under `HYP-004` to test adaptive-walk dynamics and peak-proximal epistasis against literature-derived biological expectations rather than optimizer rankings.
  - The experiment preregisters expectations from empirical diminishing-returns adaptation, Fisher-geometric-model DFE theory, and empirical rugged/sign-epistatic fitness landscapes: adaptive gains should decelerate, beneficial availability should decline near a peak, and epistasis should become more constraining near peak-proximal backgrounds.
  - The runner is `scripts/adaptive_walk_peak_epistasis_baselines.py`, with outputs `walk_trajectories.csv`, `background_neighborhoods.csv`, `epistasis_pairs.csv`, `walk_summaries.json`, `literature_baselines.json`, and `summary.json`.
  - A reduced local quick run completed successfully as a smoke test; it supported adaptive-gain deceleration and declining beneficial fraction but did not force a positive peak-proximal epistasis result.
  - The dry-run submission is `RUN-054`, and the first live execution is `RUN-055`, submitted to `lab-slurm` as scheduler job `201` on `2026-05-14`; a later labproj refresh reconciled `RUN-055` to `completed`.
  - The promoted durable outputs for `RUN-055` are now stored and DVC-tracked under `data/processed/adaptive-walk-peak-proximal-epistasis-biological-baselines/RUN-055`.
  - `RUN-055` passed only `2/7` preregistered validation tests for `HYP-004`: adaptive walks decelerated (`median greedy deceleration ratio 0.4206`) and beneficial fraction declined with fraction of peak (`Spearman -0.2748`), but the peak-proximal epistasis expectations did not validate.
  - Mean beneficial effects did not decline near the peak (`69.12` high-background mean versus `37.44` low-background mean), mean absolute epistasis rose only weakly (`Spearman 0.151`), negative epistasis decreased rather than increased (`Spearman -0.6609`), and sampled sign epistasis was absent.
  - The synthetic-alignment branch also confounds the peak-proximity axis because its near-peak bin averaged fraction of stored peak fitness `1.748`, so the stored peak is not acting as a reliable upper bound for that branch.
  - `ANA-020` records the scientific review of `RUN-055`, and `RES-020` records that the run weakens `HYP-004` as stated while preserving a narrower adaptive-walk and beneficial-availability result.
  - A literature-grounded diagnosis of the `RUN-055` failure is that the test mixed several axes that are not equivalent in the current public stability-margin model.
  - In a stability-threshold protein model, the biologically meaningful axis is expected to be robustness or consensus/stability distance, not fraction of a stored AdaptEnv peak; this is especially important because the stored peak is not a certified biological optimum and is not needed for the current `HYP-004` question.
  - The original peak-proximal epistasis prediction may also be too broad for a TEM-1-like stability model: robust consensus-like sequences can sit on a neutral plateau, while stronger deleterious and negative-epistatic effects may emerge as random drift or homologue distance exhausts stability margin.
  - `EXP-021` has therefore been created under both `HYP-004` and `HYP-005` as a direct diagnostic repeat of the `EXP-020` design on real TEM-1 beta-lactamase sequence `P62593`.
  - `EXP-021` uses the TEM-1 MMSeqs/ColabFold alignment for conditional stability, starts adaptive and random walks from real TEM-1 homologues selected by alignment-consensus bitscore distance, and validates walk/neighbourhood trends against consensus distance rather than stored peak fraction.
  - The experiment adds unconditional random walks to record fitness autocorrelation, step-effect variance, correlation-length proxies, and beneficial/neutral/deleterious random-step fractions, while retaining adaptive walks, local DFE summaries, and sampled pairwise epistasis summaries from `EXP-020`.
  - The parameter grid sweeps realistic `stability_margin` values (`4.0`, `8.0`, `12.0`) and `epistasis_strength` values (`0.0`, `0.025`, `0.05`) to test whether HYP-005 knob directions explain the failed `RUN-055` peak-proximal epistasis suite.
  - A local quick smoke of `scripts/tem1_homologue_walk_parameter_grid.py` completed using a quick-mode-only synthetic TEM-1-like alignment to avoid waiting on the public MMSeqs service; the live run still fetches and records the real MMSeqs alignment.
  - The dry-run submission is `RUN-056`; the live execution is `RUN-057`, submitted to `lab-slurm` as scheduler job `203` on `2026-05-15`, and later reconciled to `completed` with scheduler exit code `0`.
  - The promoted durable outputs for `RUN-057` are now stored and DVC-tracked under `data/processed/tem1-homologue-adaptive-random-walk-parameter-grid-diagnostic/RUN-057`.
  - `RUN-057` passed `7/9` validation checks. The TEM-1 homologue consensus-distance axis behaved as expected: background fitness declined with consensus distance (`Spearman -0.9609`), beneficial fraction increased (`Spearman 0.9531`), greedy adaptive total gain increased (`Spearman 0.8609`), greedy walks decelerated (`median ratio 0.6849`), and random-walk lag-1 autocorrelation declined weakly (`Spearman -0.2015`).
  - The `epistasis_strength` knob moved target diagnostics directionally: mean absolute epistasis increased from effectively zero to `0.004305`, and random-step standard deviation increased from `224.14` to `227.23`.
  - The `stability_margin` knob failed its two preregistered DFE-direction tests. Neutral fraction and deleterious fraction were essentially unchanged across nominal margins, and matched-start checks showed that changing nominal `stability_margin` from `4.0` to `12.0` did not change realized `background_fitness` or `stability_margin_used` for the tested identical starts.
  - `ANA-021` records the scientific review of `RUN-057`, and `RES-021` records that this supports the `HYP-004` consensus-axis adaptive/random-walk diagnostic while weakening `HYP-005` as a complete directional-knob claim because the stability-margin parameterization did not behave as intended.
  - `EXP-022` has now been created under `HYP-004` and `HYP-005` as a direct ruggedness-emergence follow-up to `RES-021`.
  - It repeats the `EXP-021` TEM-1 homologue start design but replaces the weak two-value epistasis-strength grid with explicit epistasis regimes spanning an additive control, the `RUN-057` sparse baseline, stronger sparse pairwise terms, dense pairwise terms, higher-rank pairwise matrices, covariation-enriched pair placement, very dense covarying high-rank pairwise terms, and one higher-order exact-match stress setting.
  - The runner `scripts/tem1_homologue_walk_parameter_grid.py` now records all varied epistasis parameters in the output tables and adds ruggedness-specific diagnostics: `ruggedness_parameter_score`, nonzero and strong epistasis fractions, epistasis variance, random-walk lag-1 and lag-5 autocorrelation, random-walk correlation length, and random-step class fractions.
  - A local quick smoke of the updated runner completed successfully with the quick synthetic TEM-1-like alignment. The quick smoke already showed increased random-step variance and lower autocorrelation in the stress settings, but it is not scientific evidence because it used the synthetic smoke alignment and only three reduced parameter sets.
  - The dry-run submission for `EXP-022` is `RUN-058`; the live execution is `RUN-059`, submitted to `lab-slurm` as scheduler job `206` on `2026-05-15`, and `labproj run reconcile RUN-059` currently reports it remains running.
  - `RUN-059` later completed successfully with scheduler exit code `0`, and its promoted durable outputs are now stored, DVC-tracked, and pushed under `data/processed/tem1-homologue-ruggedness-emergence-epistasis-parameter-sweep/RUN-059`.
  - `RUN-059` passed `10/13` validation checks. The consensus-distance adaptive-walk axis from `RUN-057` persisted: farther-from-consensus TEM-1 homologues had lower fitness, higher beneficial availability, larger greedy adaptive gains, and greedy walks still decelerated.
  - Stronger epistasis settings produced clear sampled pairwise epistasis: ruggedness score correlated with mean absolute epistasis (`Spearman 0.9636`), epistasis variance (`0.9636`), and nonzero epistasis fraction (`0.8146`), with the higher-order stress regime reaching mean absolute epistasis `0.6497` and nonzero epistasis fraction `0.1510`.
  - The same sweep did not produce strong random-walk ruggedness. Random-step standard deviation did not increase directionally with ruggedness score (`Spearman 0.1152`), lag-1 autocorrelation did not decline (`Spearman 0.0303`), lag-1 autocorrelation remained between `0.99928` and `0.99943`, and sign epistasis remained nearly absent even in the strongest regime (`0.000289`).
  - `ANA-022` records the scientific review of `RUN-059`, and `RES-022` records the mixed result: strong epistasis knobs support the direct `HYP-005` pairwise-epistasis subcomponent, but the result weakens the richer `HYP-004` interpretation that the current TEM-1 configuration produces biologically rich random-walk ruggedness or frequent constraining/sign epistasis.
  - `EXP-023` has now been created under `HYP-004` and `HYP-005` as a targeted architecture diagnostic following `RES-022`.
  - It disables the old generic low-rank pairwise epistasis layer, covariation-biased random pair placement, higher-order exact-match terms, and noise, then constructs sparse state-specific pairwise terms from TEM-1 MMSeqs DCA/covariation signals and routes them into latent stability through `pairwise_epistasis_override`.
  - The new runner is `scripts/tem1_dca_specific_epistasis_hamming_graphs.py`; it selects a concrete local background by realized effective stability margin, chooses DCA-coupled alternative residue states with enriched alignment support and preferably mildly deleterious single effects, enumerates the resulting `2^k` binary Hamming graph, and analyses it with `external/landscapy`.
  - Expected outputs include the DCA-selected pair table, the exact pairwise override `.npz`, the dense Hamming node table, pair epistasis table, Walsh terms, graph metrics, and a portable landscapy bundle under `landscapy_hamming_landscape/`.
  - A local quick smoke completed successfully using the quick synthetic TEM-1-like alignment: it built a 16-node landscapy Hamming cube, detected nonzero DCA-specific pair epistasis, produced one sign-epistatic selected pair in the smoke graph, and wrote the portable landscapy bundle. The smoke is an execution check only, not scientific evidence.
  - The dry-run submission for `EXP-023` is `RUN-060`; the live execution is `RUN-061`, submitted to `lab-slurm` as scheduler job `209` on `2026-05-15`.
  - `RUN-061` later completed successfully with scheduler exit code `0`, and its promoted durable outputs are now stored, DVC-tracked, and pushed under `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061`.
  - `RUN-061` confirms the DCA-specific architecture works as an isolated lever: old generic epistasis is disabled, the complete `2^10` Hamming cube contains `1024` nodes and `5120` edges, the five selected DCA-enriched pairs have mean absolute pair epistasis `57.0791`, and selected-pair sign epistasis fraction is `1.0000`.
  - The result remains scientifically mixed. Random-walk lag-1 autocorrelation drops to `0.7795`, but reciprocal sign epistasis is absent, first-order Walsh terms still explain `0.9911` of cube variance, and the graph has only one local maximum with no accessible improving shortest path from the selected anchor to the all-selected-alternative genotype.
  - `ANA-023` records the scientific review of `RUN-061`, and `RES-023` records that the result supports `HYP-005` for the DCA-specific epistasis architecture subcomponent while weakening `HYP-004` as a full multi-peak/topological-ruggedness claim.
  - `EXP-024` has now been created under `HYP-004` and `HYP-005` as the next architecture diagnostic, targeting the missing higher-order Walsh signal rather than reciprocal sign epistasis.
  - It disables the old generic epistasis layers and overlays sparse DCA/contact-graph clique exact-match haplotype terms on observed TEM-1 MSA states, testing expected orders `2` through `5` with dense `2^(n+1)` local Hamming grids around each selected term.
  - The primary readout is controllability of target-order Walsh coefficients and variance fractions, especially fifth-order Walsh signal in order-5 maps at the upper haplotype-term strength. Pairwise sign epistasis and reciprocal sign epistasis are counted in each map as secondary diagnostics.
  - The runner is `scripts/tem1_dca_clique_haplotype_epistasis.py`, with outputs including selected DCA haplotype terms, clique rejection diagnostics, dense Hamming node tables, Walsh terms, sign/RSE pair counts, and order-by-strength summaries.
  - A local quick smoke completed successfully using a quick-mode synthetic TEM-1-like alignment. In the smoke run, zero-strength maps had no target-order lift and the high-order overlay produced the expected target Walsh signal.
  - The dry-run submission for `EXP-024` is `RUN-062`; the live execution is `RUN-063`, submitted to `lab-slurm` as scheduler job `223` on `2026-05-15`, and later reconciled to `completed` with scheduler exit code `0`.
  - The promoted durable outputs for `RUN-063` are now stored, DVC-tracked, and backed up under `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063`.
  - `RUN-063` confirms that DCA clique exact-match haplotype terms are a direct high-order Walsh-control lever: all `48` requested terms were selected, zero-strength controls had no meaningful target-order signal, and target coefficients matched `strength / 2^order` with maximum error `2.842e-13`.
  - At the upper strength, order-5 maps show the preregistered fifth-order signal with mean target coefficient `93.75` and fifth-order variance fraction `0.0294`.
  - The same result remains biologically limited. In upper order-5 maps, lower Walsh orders still dominate variance (`0.2344` order 1, `0.2945` order 2, `0.2945` order 3, `0.1472` order 4, `0.0294` order 5), sign epistasis is weak (`0.0611`), RSE is weaker (`0.0194`), and maps average only `1.5` local maxima.
  - `ANA-024` records the scientific review of `RUN-063`, and `RES-024` records that the result supports `HYP-005` for the DCA clique/haplotype exact-match coefficient lever while weakening `HYP-004` as a full rugged/evolutionary-landscape claim.
  - `EXP-025` has now been created under `HYP-004` and `HYP-005` as the preregistered follow-up to the `RES-024` diagnosis.
  - The experiment tests a less hand-certified procedural alternative: sparse contact-conditioned residual fields calibrated only to aggregate local graphlet distributions, then evaluated on held-out dense Hamming graphlets sampled along adaptive walks from dramatically different TEM-1 starting states.
  - The design explicitly avoids global peak-distance conditioning because the current record shows stored peak geometry is not a reliable biological coordinate. It stratifies graphlets by observable local quantities instead: start class, walk step, background fitness, local DFE shape, local additive slope, neutrality, and alignment/consensus distance.
  - Required controls in `EXP-025` include the smooth thermodynamic core alone, the `RES-024` exact-match haplotype overlay, the `RES-023` DCA pair override, shuffled-contact and random-position residual controls, a GNK-style local non-procedural null, a residual-only embedding diagnostic, and a hard-certified per-graphlet positive control used only to expose overfitting.
  - The primary planned outputs record calibration and held-out test distributions separately for SE, RSE, local maxima, basin/accessibility metrics, random-walk autocorrelation, Walsh terms by order, and GFT coefficients by graph eigenvalue band.
  - The EXP-025 runner `scripts/tem1_aggregate_graphlet_residual_ruggedness.py` has now been implemented. It uses a common smooth-core adaptive-walk panel, then evaluates all residual candidates and controls on identical dense Hamming graphlets so topology changes are attributable to the residual/control model rather than graphlet resampling.
  - The EXP-025 implementation records calibration and held-out test graphlet metrics, Walsh terms, Hamming-graph GFT coefficients, sign/RSE pair rows, accessibility summaries, local extrema, basins, random-walk autocorrelation, and calibration-selected parameter summaries.
  - A local EXP-025 quick smoke completed successfully with `LABPROJ_QUICK=1` and a synthetic TEM-1-like alignment. The smoke produced `240` graphlets and `1440` model-graphlet metric rows, confirming execution and output shape only; it is not scientific evidence for the full real-alignment run.
  - The dry-run submission for `EXP-025` is `RUN-064`. The first live execution is `RUN-065`, submitted to `lab-slurm` as scheduler job `226` on `2026-05-15`; after reconciliation it is recorded as running.
  - `RUN-065` completed successfully with scheduler exit code `0`, and its durable outputs are now stored, DVC-tracked, and backed up under `data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065`.
  - `RUN-065` generated a large aggregate graphlet review dataset: `240` dense graphlets, `6000` model-graphlet metric rows, `720000` Walsh rows, `720000` Hamming-graph GFT rows, and `111000` sign-epistasis pair rows.
  - The primary contact residual candidate improved held-out local topology in the sampled graphlets relative to the smooth core: held-out SE increased from `0` to `0.3213`, RSE from `0` to `0.1030`, mean local maxima from `1.0` to `3.54`, accessible-path fraction decreased from `0.3097` to `0.0350`, and lag-1 random-walk autocorrelation decreased from `0.6455` to `0.4976`.
  - The same primary residual preserved a low-order-dominated spectrum with a nonzero high-order tail: held-out Walsh high-order fraction was `0.0737`, and Hamming-graph GFT low/mid/high energy fractions were `0.7719`, `0.2197`, and `0.0084`.
  - The main scientific limitation is material. Although `RUN-065` generated starts and graphlet backgrounds across five start classes, all `240` evaluated graphlets came from `reference_or_consensus_like` starts, so the run does not test the intended generalization across dramatically different adaptive-walk regions.
  - Matched controls also challenge the biological/contact-specific interpretation: shuffled-contact and random-position residual controls closely matched the primary contact residual, while the DCA pair override produced stronger topology without a high-order Walsh tail, and the GNK local null produced still stronger but non-procedural high-order-heavy topology.
  - `ANA-025` records the scientific review of `RUN-065`, and `RES-025` records the result: aggregate residuals can create local held-out topology in the sampled reference-like graphlets, but `RUN-065` weakens `HYP-004` as a broad ruggedness/generalization claim while only narrowly supporting `HYP-005` for residual and pairwise topology knobs.
  - The aggregate graphlet runner now makes the relative-scale mechanism explicit. Residual parameter sets can use a fixed walk-initial beneficial-effect scale rather than always scaling to the current local DFE, and every graphlet/model row now records base WHT variance by order, residual WHT variance by order, residual high-order absolute variance, residual-to-base order-1 variance, residual-to-available-beneficial-effect ratios, distance to the capped walk endpoint, and walk-initial DFE scales.
  - The same runner now applies balanced round-robin graphlet-background ordering before the graphlet cap, so a capped run cannot silently consume only the first start class as `RUN-065` did. Walsh/GFT term output can also be capped while preserving full order-level summaries in metric rows.
  - `EXP-026` has been created under `HYP-004` and `HYP-005` as the scaled repeat of `EXP-025`. It tests the revised causal chain: fixed residual/high-order components should become more visible as smooth adaptive gains diminish, rather than diminishing returns directly creating SE/RSE.
  - `EXP-026` deliberately reduces output volume: one TEM-1 landscape seed, one start per class, graphlet dimensions `4` to `6`, two contact-source families, at most `90` graphlets, no full pair-row output, and only the top Walsh/GFT terms per graphlet/model. The full order-level WHT/GFT and relative-exposure summaries remain in the metric tables.
  - A local quick smoke for `EXP-026` completed successfully with synthetic TEM-1-like alignment. It produced balanced graphlets across all five start classes under the cap (`8` graphlets per class in quick mode), reduced term tables (`2880` Walsh rows and `2880` GFT rows), and confirmed that the new relative-scale metric columns are populated. The smoke is execution evidence only, not a scientific result.
  - The dry-run submission for `EXP-026` is `RUN-066`. The live execution is `RUN-067`, submitted to `lab-slurm` as scheduler job `227` on `2026-05-15`; it completed successfully with scheduler exit code `0`.
  - `RUN-067` produced the intended scaled dataset and corrected the `RUN-065` sampling failure: `90` graphlets, `2250` model-graphlet metric rows, `27000` Walsh rows, `27000` GFT rows, no sign-pair rows, and about `15.8` MB of DVC-tracked durable output under `data/processed/tem1-relative-scale-residual-exposure-graphlet-diagnostic/RUN-067`.
  - The evaluated graphlets are exactly balanced across the five start classes (`18` each) and adaptive-walk stages (`30` early, `30` middle, `30` late). Smooth-core walks show diminishing-return proxies: mean available beneficial effect declines from `25.01` early to `20.29` middle to `17.36` late, while best available effect declines from `245.06` to `105.92` to `82.69`.
  - The smooth core still has no SE/RSE or high-order topology in held-out graphlets (`SE = 0`, `RSE = 0`, mean local maxima `1.00`, Walsh high-order fraction `3.55e-28`), so the result supports the corrected causal language: diminishing returns alone does not create WHT high-order signal or SE/RSE.
  - The selected contact residual (`contact_residual_ensemble__pset-003`) creates held-out topology over the smooth core across the balanced graphlet panel: SE `0.3511`, RSE `0.1481`, mean local maxima `2.56`, accessible-path fraction `0.0323`, lag-1 autocorrelation `0.3634`, Walsh high-order fraction `0.0855`, and GFT low/mid/high fractions `0.7197`, `0.2603`, and `0.0200`.
  - The stronger fixed-scale unmasking prediction is not supported. The fixed walk-initial contact residual has weak held-out topology (`SE = 0.1409`, `RSE = 0.0216`, mean local maxima `1.18`), does not show a monotone stage-wise SE/RSE rise, and calibration selects local-DFE scaling rather than fixed walk-initial scaling.
  - `ANA-026` records the scientific review of `RUN-067`, and `RES-026` records the result: relative-scale exposure diagnostics are useful and residual topology generalizes beyond the `RUN-065` reference-like sampling defect, but the current model does not validate a fixed-scale peak-proximal SE/RSE mechanism or contact-specific biological ruggedness. `RES-026` narrowly supports `HYP-005` and keeps `HYP-004` weakened.
  - `EXP-027` has been created under `HYP-004` and `HYP-005` as the fixed-scale residual ratio sweep repeat of `EXP-026`, motivated directly by `RES-026`.
  - The experiment removes the local-DFE-scaled primary residual setting and scans only `constant_walk_initial_beneficial_effect` residual scaling with fixed multipliers `1.0`, `4.0`, and `16.0`. In this mode, the multiplier is applied to the walk-initial mean available beneficial effect rather than the graphlet-local DFE.
  - To keep total compute and storage comparable with `EXP-026`, the graphlet budget is reduced from `90` to `30`, equal to the previous budget divided by the three scale-grid values. The same compact-output policy is retained: no sign-pair rows, capped Walsh/GFT term rows, and full order-level/statistical summaries in the metric rows.
  - The aggregate graphlet runner now balances graphlet tasks, not only graphlet backgrounds, before applying a cap. A quick synthetic smoke of `EXP-027` produced `15` graphlets spread evenly across the five start classes (`3` each) and three walk stages (`5` early, `5` middle, `5` late), confirming that the reduced cap still samples the intended mechanism axis.
  - The dry-run submission for `EXP-027` is `RUN-068`. The live execution is `RUN-069`, submitted to `lab-slurm` as scheduler job `228` on `2026-05-15`; reconciliation records it as completed with scheduler exit code `0`.
  - `RUN-069` outputs have now been collected, DVC-tracked, and backed up under `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069`. The DVC pointer is `data/processed/tem1-fixed-scale-residual-ratio-sweep-graphlet-diagnostic/RUN-069.dvc`.
  - `RUN-069` produced the intended compact fixed-scale sweep dataset: `30` dense graphlets, `600` graphlet/model metric rows, `7200` Walsh rows, `7200` Hamming-graph GFT rows, no sign-epistasis pair rows, and about `4.9` MB of collected output. Graphlets are balanced across five start classes (`6` each) and three walk stages (`10` early, `10` middle, `10` late).
  - The smooth thermodynamic core again shows diminishing-return proxies along sampled adaptive walks: mean available beneficial effect declines from `24.56` early to `22.62` middle to `20.90` late, and best available effect declines from `213.70` to `174.33` to `162.71`. The smooth core still has no held-out SE/RSE or high-order topology (`SE = 0`, `RSE = 0`, mean local maxima `1.00`, Walsh high-order fraction `2.32e-29`).
  - The fixed contact residual scale grid is directional in held-out graphlets. Across fixed scales `1`, `4`, and `16`, SE rises `0.2391 -> 0.2935 -> 0.4783`, RSE rises `0.0417 -> 0.0978 -> 0.1449`, local maxima rise `1.22 -> 1.48 -> 1.91`, WHT/GFT high-order fraction rises `0.00074 -> 0.00439 -> 0.05055`, and lag-1 random-walk autocorrelation declines `0.4560 -> 0.4203 -> 0.3399`.
  - Calibration now selects the intended fixed-scale contact residual (`contact_residual_ensemble__pset-002`, the highest fixed scale), correcting the central `RUN-067` interpretive problem where the selected contact residual used local-DFE scaling.
  - The selected fixed scale has a low-order-dominated spectrum with a nonzero high-order tail: GFT/Walsh low/mid/high fractions are `0.8447`, `0.1048`, and `0.0505`. Compared with the fixed contact setting in `RUN-067`, `RUN-069` raises SE from `0.1409` to `0.4783`, RSE from `0.0216` to `0.1449`, local maxima from `1.18` to `1.91`, and WHT high-order fraction from `0.0043` to `0.0505`.
  - The stronger late-stage unmasking prediction is still not cleanly supported. At the highest fixed scale, residual-to-available-beneficial-effect ratio rises from `3.83` early to `4.72` middle to `6.99` late, but SE (`0.4667`, `0.4833`, `0.4458`), RSE (`0.2250`, `0.1167`, `0.1083`), and WHT high-order fraction (`0.0572`, `0.0275`, `0.0457`) are not late-dominated.
  - Matched controls continue to limit biological/contact-specific interpretation. The selected GNK null is not stronger in this run, but the DCA pair override reaches nearly the same SE without high-order Walsh signal, and random-position residuals have comparable RSE/local maxima. `ANA-027` records the analysis and `RES-027` records the result: `RUN-069` supports `HYP-005` for the fixed residual-scale knob while keeping `HYP-004` weakened as a broad biological ruggedness claim.
  - Scratch analytical notes have been created under `scratch/model_signature_derivation.md` and `scratch/model_literature_falsification_review.md` to formalize the current structural interpretation.
  - The scratch derivation separates the model into a smooth thermodynamic/geometric core plus sparse fixed-scale residuals: the core explains nonlinear global epistasis, diminishing returns, stability-threshold neutrality, and deleterious tails, while residual terms explain SE/RSE, local maxima, reduced accessibility, shortened random-walk autocorrelation, and a low-order-dominated WHT/GFT spectrum with a high-order tail.
  - The literature review treats the mechanism as a claim to falsify. It finds that the ingredients are grounded in FGM/global-epistasis, stability-threshold, GNK/Walsh-sparsity, and SE/RSE topology literature, but the exact AdaptEnv residual generator remains a procedural design choice rather than a derived or validated biological contact-ruggedness model.
  - `EXP-028` has been created under `HYP-004` and `HYP-005` as an explicit empirical falsification suite for the smooth-core plus fixed-scale residual ruggedness mechanism.
  - The experiment directly tests six breakpoints: smooth-core-only empirical fit, no residual-scale effect, wrong empirical spectral envelope, control equivalence, wrong evolutionary dynamics, and scale confounding after nonlinear global-readout correction.
  - The implemented runner is `scripts/empirical_graphlet_falsification_suite.py`. It runs a compact RUN-069-compatible synthetic graphlet suite, extracts complete empirical ProteinGym binary Hamming graphlets from cached PHOT, SPTN1, and KRAS prepared assays, computes matched SE/RSE/topology/WHT/GFT statistics, fits an additive-plus-polynomial global-readout correction within empirical graphlets, and writes direct synthetic-vs-empirical model-distance summaries.
  - A local quick smoke completed successfully with quick synthetic graphlet limits and three PHOT empirical graphlets. The smoke verified the entrypoint and output shape only; it is not scientific evidence for the registered experiment.
  - The dry-run submission for `EXP-028` is `RUN-070`. The live execution is `RUN-071`, submitted to `lab-slurm` as scheduler job `230` on `2026-05-15`; it completed with scheduler exit code `0`.
  - `RUN-071` outputs have been collected, DVC-tracked, and backed up under `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071`. The DVC pointer is `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071.dvc`.
  - `RUN-071` produced `60` synthetic graphlets and `56` empirical graphlets. The primary PHOT empirical graphlets have SE `0.5979`, RSE `0.2524`, mean local maxima `4.9583`, WHT high-order fraction `0.3944`, and GFT low/mid/high fractions `0.6056`, `0.2472`, and `0.1473`.
  - Four registered falsifiers did not break the model: smooth-core-only empirical fit is rejected, fixed residual scale is directional across scales `1`, `4`, and `16`, matched GNK/random/shuffled controls do not beat the contact residual under the registered distance criterion, and nonlinear global-readout correction does not remove the high-order signal.
  - Two falsifiers are material. The primary empirical PHOT graphlets violate the expected modest high-order-tail spectral envelope because WHT high-order fraction is `0.3944`, above the preregistered `0.35` ceiling; and no embedded procedural candidate matches empirical spectra, SE/RSE, local maxima, accessibility, and random-walk statistics as a package.
  - The closest model to the primary empirical target is the residual-only diagnostic (`distance = 0.7387`), not the embedded contact residual (`0.8432`). The contact residual still beats smooth core (`1.9990`), GNK (`1.2414`), random-position residual (`1.2853`), and shuffled-contact residual (`1.3002`), but it underestimates empirical local maxima and WHT high-order fraction.
  - `ANA-028` records the falsification-focused analysis of `RUN-071`, and `RES-028` records the result: the fixed residual-scale knob remains supported narrowly under `HYP-005`, but the strong current residual-ruggedness benchmark-adequacy claim is falsified and `HYP-004` is further weakened.
  - `EXP-029` has been created under the historical `HYP-001` lineage as a compatibility/refit experiment linked to `RES-005`, `RES-006`, `RES-007`, and `RES-025` through `RES-028`.
  - The experiment records the stale-cache concern raised after `RES-028`: prior empirical graphlet targets from PHOT, SPTN1, and KRAS use cached prepared assay outputs and older latent/readout semantics, so they are not valid parameter sources for the current smooth-core plus fixed-scale residual model.
  - `EXP-029` requires raw ProteinGym systems to be refit with current-model ABC before empirical parameters are used. Old cached `latent_phi`, `mavenn_yhat`, and prepared assay outputs are allowed only as stale-cache replay controls.
  - The experiment remains scoped by `CLAIM-001`: even if the refit improves empirical alignment, it does not by itself restore the original strong shared-regime version of `HYP-001`; shared or hierarchical generalization across assays would still have to be demonstrated.
  - `EXP-029` has now been configured for a restricted current-model ABC refit using only `PHOT_CHLRE_Chen_2023` and the previously used Tsuboyama stability assay `SPTN1_CHICK_Tsuboyama_2023_1TUD`.
  - The new guarded runner is `scripts/proteingym_current_model_abc_refit.py`. It delegates to the existing SMC-ABC calibration pathway but blocks legacy latent-trait, peak-distance, noise, FGM, deterministic roughness, and higher-order exact-match parameters from being fitted or silently activated.
  - The fitted parameter set is limited to current stability-margin and stability-routed pairwise epistasis semantics: stability scale/margin/threshold/conservation terms plus sparse empirical pairwise controls targeted to latent stability.
  - The dry-run submission is `RUN-072`. The live execution is `RUN-073`, submitted to `lab-slurm` as scheduler job `232` on `2026-05-16` with `STACK-000`, `8` CPUs, `64` GB memory, and `24:00:00` walltime.
  - `RUN-073` completed successfully with scheduler exit code `0`, and its durable outputs are now DVC-tracked under `data/processed/proteingym-current-model-abc-refit-residual-ruggedness/RUN-073`.
  - The current-model preflight passed: the run used `observed_fitness_combine_mode = stability_margin`, routed generic epistasis to latent stability, and fitted only the current stability-margin/stability-routed residual parameter set rather than the older FGM, designed-peak, deterministic roughness, noise, or higher-order exact-match parameters.
  - Matched synthetic truth recovery passed for the restricted current-model inverse problem, with both synthetic truths recovering all `20/20` fitted truth parameters inside posterior q90 intervals.
  - The empirical PHOT+SPTN1 shared ABC fit remains weak as a validation result. The final best objective distance was `10.8228`, close to the round-0 best `10.7975`, indicating little SMC concentration; the best branch had functional KS `0.1862`, single-mutant holdout Spearman `0.4767`, double-mutant holdout Spearman `0.7691`, and epistasis-prediction Spearman only `0.1681`.
  - Compared with prior assay-specific runs, `RUN-073` does not rescue empirical epistasis prediction: `RUN-024` SPTN1 reached epistasis-prediction Spearman `0.5402`, and `RUN-039` PHOT reached `0.6344`.
  - The run did not emit current-model dense graphlet, WHT, GFT, local maxima, accessibility, or random-walk autocorrelation outputs. It therefore resolves the stale-cache concern only for the restricted PHOT+SPTN1 ABC refit, not for the full empirical graphlet target-envelope question raised by `RES-028`.
  - `ANA-029` records the analysis of `RUN-073`, and `RES-029` records the result: the current-semantics ABC refit is procedurally valid, but it does not restore a shared empirical calibration under `HYP-001` and does not overturn the `RUN-071` graphlet falsification.
  - `HYP-008` has now been created under `AIM-001`: local non-diffuse ruggedness may require nonlinear operators over multiple latent molecular traits, and the current smooth-core plus diffuse residual model may not capture that behaviour as a biologically grounded procedural benchmark.
  - This is deliberately distinct from superseded `HYP-007`. The new hypothesis is not "more fitted traits will improve all empirical recovery"; it is a sharper structural claim that empirical SE/RSE, local maxima, accessibility, autocorrelation, and WHT/GFT spectra may be localized around specific sequence/trait regions rather than spread as a diffuse residual field.
  - `HYP-008` is motivated by `RES-028`, but constrained by `CLAIM-001` and `CLAIM-002`: any future support must beat matched single-trait stability-margin, diffuse residual, residual-only, global-readout, and GNK/RMF-style controls on held-out graphlets, not merely add unconstrained capacity.
  - A review of `RUN-073` after `RES-029` found a scientific design error: the run pooled `PHOT_CHLRE_Chen_2023` and `SPTN1_CHICK_Tsuboyama_2023_1TUD` into one ABC target, which conflicts with the system-specific calibration interpretation now recorded in `CLAIM-001` and the recent empirical record.
  - The same review found that `RUN-073` SMC sampling was not productively using the fourth empirical round. Round best distances were approximately `10.7975`, `10.9120`, `10.7028`, and `10.8228`, and the posterior remained broad rather than contractive, so the next run reduces the empirical budget to three rounds, 48 particles, one simulation replicate, and 5/6 pool multipliers.
  - `EXP-030` has been created under `HYP-001`, `HYP-004`, and `HYP-005` to correct the pooling error directly. It fits PHOT and SPTN1 as two independent current-model landscapes in one run, using fresh raw ProteinGym preprocessing and the same current stability-margin/stability-routed epistasis semantics.
  - The new runner, `scripts/proteingym_system_specific_current_graphlet_refit.py`, delegates each assay to the existing ABC calibration pathway, then computes per-system graphlet SE/RSE, local maxima, accessibility, random-walk autocorrelation, WHT, and GFT statistics for empirical graphlets, fitted current models, smooth-core controls, residual-only diagnostics, GNK-style controls, and nonlinear global-readout residual controls.
  - `EXP-030` explicitly retests the `RES-028` fired falsifiers separately for each system: wrong spectral envelope and wrong embedded evolutionary/topological dynamics. It also retains the scale-confounding and smooth-core-only screens needed to interpret those falsifiers.
  - A preflight graphlet search on the `RUN-073` prepared assays found that PHOT has complete 4D, 5D, and 6D graphlets under the planned setting, while SPTN1 has complete 2D graphlets only and no complete 3D or 4D graphlets up to a top feature pool of 120. Therefore the high-order spectral falsifier may be formally non-testable for SPTN1 under the complete-graphlet criterion.
  - Scratch planning for the future `HYP-008` coupled-landscape direction has been recorded under `scratch/coupled_landscape_nonlinear_operator_plan.md` and `scratch/coupled_landscape_regularisation_metrics.md`. The proposed model keeps a frozen system-specific base landscape, adds correlated auxiliary landscapes, localizes innovation energy in modules or patches, combines traits through a simple nonlinear operator, and regularizes against loss of base correlation, diffuse innovation, DFE distortion, and GNK-like spectral overfitting.
  - `RUN-075` failed after completing the full PHOT system-specific ABC fit. The failure was an implementation bug in EXP-030 graphlet post-processing: graphlet dictionaries passed to the shared GFT writer omitted the `split` field required by `gft_rows_from_walsh`, producing `KeyError: 'split'`.
  - The EXP-030 runner has been patched to include `split` in graphlet metadata and to reuse completed per-system fit outputs from a prior run when `resume.source_run_id` and `reuse_completed_system_fits` are configured. This allows the replacement run to continue from `RUN-075`'s completed PHOT checkpoint instead of rerunning the PHOT ABC fit.
  - The replacement dry-run is `RUN-076`. The live resumed execution is `RUN-077`, submitted to `lab-slurm` as scheduler job `236`. Early remote logs confirm that `RUN-077` reused the completed PHOT fit from `RUN-075`, generated PHOT graphlet WHT/GFT outputs past the previous crash point, and then started the SPTN1 system-specific ABC fit.
  - `RUN-077` completed successfully with scheduler exit code `0`, and the collected outputs are DVC-tracked under `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077`.
  - `RUN-077` produced current system-specific PHOT and SPTN1 ABC fits plus graphlet metrics, Walsh terms, GFT coefficients, model-distance tables, and falsification-test summaries. The run used independent PHOT and SPTN1 landscapes, passed the current-model legacy-parameter guard, and emitted the graphlet evidence missing from pooled `RUN-073`.
  - The PHOT held-out graphlets partially reverse the strongest `RES-028` dynamics failure: the fitted current model is the closest held-out candidate for the joint graphlet package (`distance = 0.4016`), ahead of smooth core (`0.4119`), posterior mean (`0.4408`), GNK (`0.7594`), and residual-only diagnostic (`2.6756`).
  - The PHOT spectral-envelope falsifier still fires. Held-out PHOT graphlets have WHT high-order fraction `0.3788`, above the preregistered `0.35` ceiling, and the high-order signal is not removed by nonlinear global-readout correction. The fitted current model is also only marginally separated from smooth core, so the residual mechanism is not cleanly validated.
  - The SPTN1 graphlets are smooth 2D controls: SE/RSE are `0`, local maxima average `1`, and high-order WHT/GFT falsifiers are not testable under the complete-graphlet criterion.
  - `ANA-030` records the `RUN-077` graphlet falsifier analysis, and `RES-030` records the mixed result: system-specific current refitting repairs part of the PHOT topological/dynamical mismatch but does not rescue the empirical spectral envelope or the strong `HYP-004` benchmark-adequacy claim; it further motivates `HYP-008`.
  - `EXP-031` has been created under `HYP-008` as a cheap deterministic pre-ABC machinery test for the coupled-landscape idea. It freezes the fitted PHOT base landscape from `RUN-077`, uses calibration graphlets to select localized empirical innovation modules from plaquette energy, and grid-searches coupled auxiliary/module parameters rather than jointly refitting both landscapes.
  - The new shared implementation is `src/project_adapt_env/coupled_landscape.py`, and the registered runner is `scripts/phot_coupled_landscape_grid_search.py`. The runner rebuilds PHOT fitted-current, posterior-mean, smooth-core, empirical, empirical-readout-residual, GNK, shuffled-module, random-module, uncoupled-module, no-module, and coupled candidates on identical complete graphlets, then writes graphlet WHT/GFT/topology metrics, base-correlation summaries, model distances, and module-locality energy tables.
  - `EXP-031` explicitly tests whether localized coupled innovation improves the still-firing PHOT graphlet falsifiers from `RES-030` without an expensive ABC run. It also checks the central structural risk for `HYP-008`: a coupled model that only improves global distance by diffuse innovation, base decorrelation, or GNK-like spectral overfitting should fail the module-locality and control-comparison summaries.
  - Local verification completed before submission: Python byte-compilation passed for the new module and runner, `LABPROJ_QUICK=1` produced the expected compact outputs, and a full local smoke completed with `24` PHOT graphlets, `12` selected modules, `144` coupled parameter sets, and about `131` MB of outputs. These smokes are execution checks only, not scientific evidence.
  - The dry-run submission for `EXP-031` is `RUN-078`. The live execution is `RUN-079`, submitted to `lab-slurm` as scheduler job `237` on `2026-05-16`; the registered Slurm request is `4` CPUs, `32` GB memory, and `04:00:00` walltime. Immediately after submission, `squeue` showed job `237` running on `ripley`.
  - `RUN-079` completed successfully with scheduler exit code `0`, and the collected outputs are DVC-tracked and pushed under `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079`.
  - `RUN-079` produced the planned deterministic PHOT grid outputs: `24` dense graphlets, `12` empirical modules, `144` coupled parameter sets, graphlet WHT/GFT/topology metrics, model distances, base-correlation summaries, and module-energy locality tables.
  - `ANA-031` records the `RUN-079` review and `RES-031` records a negative result for the current local-module version of `HYP-008`. The calibration-selected no-module coupled operator improves held-out PHOT graphlet distance over frozen current (`0.3137` versus `0.4016`) while preserving base correlation, but the calibration-selected empirical plaquette-module coupled candidate is not competitive (`0.3982`) against no-module, random-module (`0.3546`), shuffled-module (`0.3946`), or frozen-current controls.
  - The module-locality diagnostic is also negative. The empirical-module coupled candidate has plaquette effective module count `9.33/12` and normalized entropy `0.942`, close to the empirical readout residual's diffuse `9.72/12` and `0.952`. The result weakens `HYP-008` as currently operationalized: graphlet distance can be moved by a coupled operator, but not yet by a validated local, non-diffuse module mechanism.
  - `EXP-032` has been created under `HYP-008` as the hard-min weakest-link follow-up to `RES-031`. It freezes the `RUN-077` PHOT fit, defines the primary first landscape as the fitted smooth/additive PHOT base, constructs correlated and offset auxiliary landscapes, and combines them by `f(x) = min(a(x), b(x))`.
  - The new runner is `scripts/phot_hard_min_coupled_landscape_grid_search.py`, with shared hard-min and switch-boundary utilities added to `src/project_adapt_env/coupled_landscape.py`. The runner records full PHOT assay fit metrics, DFE summaries, double-mutant epistasis prediction metrics, dense graphlet WHT/GFT/topology metrics, model distances, hard-min limiter-switch diagnostics, and edge/plaquette local-signal concentration summaries.
  - Local verification completed before submission: Python byte-compilation passed, and `LABPROJ_QUICK=1` completed with `2500` PHOT assay variants, `2` quick graphlets, `8` hard-min grid points, and all expected output tables. The quick smoke is execution evidence only, not scientific evidence.
  - The dry-run submission for `EXP-032` is `RUN-080`. The live execution is `RUN-081`, submitted to `lab-slurm` as scheduler job `238` on `2026-05-16`; the registered Slurm request is `4` CPUs, `32` GB memory, and `04:00:00` walltime. Immediately after submission, `squeue` showed job `238` running on `ripley`.
  - `RUN-081` completed successfully with scheduler exit code `0`, and the collected outputs are DVC-tracked and backed up under `data/processed/phot-hard-min-coupled-weakest-link-grid/RUN-081`.
  - `RUN-081` produced the planned deterministic PHOT hard-min grid outputs: `167,529` assay variants, `24` dense graphlets split evenly into calibration and test sets, `180` hard-min parameter sets, assay fit tables, graphlet metrics, WHT/GFT tables, model distances, hard-min switch diagnostics, and local-signal energy summaries.
  - The best hard-min full-assay objective improved over frozen controls (`0.4223` versus `0.4814` smooth core, `0.4856` fitted current, and `0.4631` posterior mean), but the improvement was narrow: holdout Spearman was essentially unchanged, single-mutant Spearman was slightly worse than controls, and epistasis-prediction Spearman remained weak (`0.0106`).
  - The central held-out graphlet test was negative. The best hard-min test graphlet distance was `0.4124`, which did not beat fitted current (`0.4016`) or smooth core (`0.4119`), and the calibration/assay-selected hard-min candidate had worse held-out distance (`0.4631`).
  - Hard-min coupling did not repair the PHOT spectral/topology mismatch: empirical test WHT high-order fraction was `0.3788`, while the best assay hard-min was `0.5077`; empirical test graphlets had mean local maxima `5.17`, accessibility `0.1975`, and lag-1 autocorrelation `0.1905`, while the best assay hard-min had `3.42`, `0.0245`, and `0.0464`.
  - `ANA-032` records the `RUN-081` review and `RES-032` records the result: hard-min weakest-link coupling can improve the PHOT assay objective, but it does not pass the held-out graphlet or locality controls required to support `HYP-008`; it therefore weakens the current hard-min implementation of the local nonlinear trait-operator hypothesis.
  - `EXP-033` has been created under `HYP-008` as the guarded soft-min/Boltzmann follow-up to `RES-032`. It freezes the `RUN-077` PHOT current-model output, uses the fitted smooth-core PHOT landscape as the primary `a(x)` anchor, constructs a correlated auxiliary landscape `b(x)`, and fits only a small ABC parameterization of the two-state operator rather than reopening the base landscape fit.
  - The registered operator is `g_tau(a,b) = -tau log(exp(-a/tau) + exp(-b/tau))`, with `tau` fitted through `log_tau`. The ABC parameterization includes auxiliary correlation, auxiliary scale, auxiliary offset, auxiliary pairwise scale, a discrete auxiliary pairwise density choice, and `log_tau`; held-out graphlets remain reserved for scientific support.
  - The shared implementation in `src/project_adapt_env/coupled_landscape.py` now includes soft-min, two-state Boltzmann, and soft-state diagnostic utilities. The registered runner is `scripts/phot_softmin_coupled_landscape_abc.py`, which writes PHOT assay guardrail metrics, graphlet topology/WHT/GFT metrics, model-distance summaries, posterior particles/rounds, and locality/state-transition diagnostics.
  - Local verification completed before submission: Python byte-compilation passed, and `LABPROJ_QUICK=1` completed with `4` ABC particles, `2` ABC rounds, `2` quick PHOT graphlets, `2500` assay variants, and the expected compact output tables. The quick smoke is execution evidence only, not scientific evidence.
  - The dry-run submission for `EXP-033` is `RUN-082`. The first live submission record, `RUN-083`, failed before Slurm job creation because the sandboxed submit environment could not resolve `raspi1.maspence.com`; it has no scheduler job id.
  - The successful live execution is `RUN-084`, submitted to `lab-slurm` as scheduler job `239` on `2026-05-17`; the registered Slurm request is `4` CPUs, `32` GB memory, and `08:00:00` walltime. Immediately after submission, `squeue` showed job `239` running on `ripley`.
  - `RUN-084` completed successfully with scheduler exit code `0:0`, and the collected outputs are DVC-tracked and backed up under `data/processed/phot-softmin-boltzmann-coupled-landscape-abc/RUN-084`.
  - `RUN-084` produced the planned PHOT soft-min/Boltzmann ABC output set: `167,529` assay variants, `24` dense graphlets split evenly into calibration and held-out test sets, `32` ABC particles, `3` SMC rounds, `26` posterior candidates evaluated on graphlets, assay-fit tables, graphlet metrics, WHT/GFT tables, model distances, and local-signal/state-transition summaries.
  - The ABC-selected best particle improved calibration graphlet distance (`0.6100`) over fitted current (`0.7318`), smooth core (`0.7454`), fitted-current posterior mean (`0.7839`), and GNK (`0.8411`) while preserving base correlation in the ABC objective. This is a technical calibration signal, not held-out support.
  - The central held-out graphlet test was negative. The calibration-selected best soft-min particle had test distance `0.4434`, worse than fitted current (`0.4016`) and smooth core (`0.4119`); the posterior mean was also worse (`0.4543`). One posterior particle narrowly beat fitted current on held-out distance (`0.3964`), but it had poor calibration distance (`0.7126`) and is therefore a test-oracle capacity signal only.
  - The assay guardrail was only narrowly positive. The posterior mean improved aggregate assay objective (`0.4570` versus `0.4631` posterior mean control, `0.4814` smooth core, and `0.4856` fitted current), but holdout and single/double-mutant validation correlations did not improve as a package, and epistasis-prediction Spearman remained negative (`-0.0881`).
  - Locality did not validate HYP-008. The calibration-selected best particle generated very large held-out plaquette innovation energy (`179,908`) relative to the empirical readout residual (`15.4`) and weak shared-feature covariance (`0.0934` plaquette-energy similarity Spearman), so concentrated state-transition energy cannot be interpreted as empirical-like local non-diffuse ruggedness.
  - `ANA-033` records the `RUN-084` review and `RES-033` records the result: soft-min/Boltzmann ABC can fit calibration graphlets and slightly improve an assay objective, but it does not pass the calibration-selected held-out graphlet and locality controls required to support `HYP-008`; it weakens this tested implementation.
  - `EXP-034` has now been created under `HYP-008` as the direct budgeted follow-up to `RES-033`. It keeps the `RUN-077` PHOT base frozen, fits the state-difference field `d(x)=b(x)-a(x)` rather than a free globally correlated auxiliary landscape, selects module coordinates from calibration empirical readout-residual plaquette energy, and uses `g_tau(a,d)=a - tau log(1 + exp(-d/tau))`.
  - The new shared implementation in `src/project_adapt_env/coupled_landscape.py` adds deterministic module-residual coefficients, module-tied difference-field evaluation on graphlets and assay variants, soft-min-from-difference evaluation, finite-difference energy metrics, state-boundary roughness metrics, and reusable budget penalties.
  - The registered runner is `scripts/phot_budgeted_softmin_coupled_landscape_abc.py`. Its ABC objective still fits calibration graphlets only, reports the PHOT assay objective as a guardrail, and now penalizes excessive innovation variance, empirical-residual-relative edge/plaquette energy, rough state boundaries, and overly low effective temperature relative to the local edge scale of `d(x)`.
  - Local verification before submission passed by byte-compiling the new module and runner, then running `LABPROJ_QUICK=1` with temporary outputs. The quick smoke emitted the expected compact posterior, graphlet, assay, module, and locality tables; this is execution evidence only, not scientific support.
  - The dry-run submission for `EXP-034` is `RUN-085`. The live submission is `RUN-086`, submitted to `lab-slurm` as Slurm job `241` on `2026-05-17`; the registered request is `4` CPUs, `32` GB memory, and `08:00:00` walltime. After reconciliation, `RUN-086` was recorded as running, and `squeue` showed job `241` running on `ripley`.
  - `RUN-086` completed, and its outputs have now been collected with `labproj` and DVC-tracked under `data/processed/phot-budgeted-module-tied-softmin-abc/RUN-086`.
  - `RUN-086` produced the planned budgeted module-tied soft-min ABC output set: `167,529` PHOT assay variants, `24` dense graphlets split evenly into calibration and held-out test sets, `32` ABC particles, `3` SMC rounds, `26` posterior candidates evaluated on graphlets, module catalogs, assay-fit tables, graphlet WHT/GFT/topology metrics, model distances, and locality/energy summaries.
  - The budgeted operator produced only a weak calibration signal. The ABC-selected best particle improved calibration graphlet distance slightly over frozen controls (`0.7240` versus `0.7318` fitted current and `0.7454` smooth core), but it did not beat the empirical readout-residual diagnostic (`0.5428`).
  - The central held-out graphlet test remained negative. The calibration-selected best particle had held-out distance `0.4060`, worse than fitted current (`0.4016`), though slightly better than smooth core (`0.4119`). The shuffled-module control was also better than the empirical-module selected particle (`0.4029` versus `0.4060`), so module-specific locality is not supported.
  - The budget constraints prevented the large artificial innovation seen in `RUN-084` by collapsing the selected operator almost exactly onto the frozen base: values-to-base correlation was `1.0`, edge-gradient correlation was `1.0`, innovation variance ratio to base was `1.4e-25`, and edge/plaquette innovation ratios to the empirical residual were `6.0e-22` and `1.2e-21`.
  - The assay guardrail did not rescue the run. The best assay objective remained the fitted-current posterior mean (`0.4631`), while the best budgeted candidate was worse (`0.5350`) and the calibration-selected best particle worse again (`0.5560`).
  - `ANA-034` records the `RUN-086` review and `RES-034` records the result: explicit budgeted module-tied soft-min ABC does not pass the held-out graphlet, module-control, locality, or assay-guardrail checks required to support `HYP-008`; it weakens this tested implementation.
  - `EXP-035` has now been created under `HYP-008` as a replicated held-out split diagnostic rather than another operator variant. It reuses the completed `RUN-079` PHOT graphlet candidate metrics, keeps the same fitted-current, smooth-core, no-module, empirical-module, random-module, shuffled-module, uncoupled-module, residual-only, and GNK candidate families, and resamples calibration/test graphlet assignments many times while stratifying by graphlet dimension.
  - The registered runner is `scripts/phot_replicated_split_graphlet_diagnostic.py`. For each split, it selects the best candidate within each family using calibration graphlet distance only, then reports held-out graphlet distance and secondary SE/RSE, local maxima, accessibility, lag-1 autocorrelation, WHT high-order, and GFT energy components.
  - The purpose is to decide whether the `RUN-079` no-module held-out improvement (`0.3137` versus `0.4016` fitted current on the original split) is stable or a split artifact before proposing another nonlinear operator. Local quick validation passed and emitted the planned summary, split-assignment, selected-distance, family-summary, pairwise-delta, component-summary, and selection-count tables.
  - The dry-run submission for `EXP-035` is `RUN-087`. The first live submission record, `RUN-088`, failed before Slurm job creation because the sandboxed submit environment could not resolve `raspi1.maspence.com`; it has no scheduler job id. The approved network/SSH resubmission is `RUN-089`, submitted to `lab-slurm` as Slurm job `242` on `2026-05-17` with `2` CPUs, `16` GB memory, and `01:00:00` walltime. Remote Slurm showed job `242` running on `ripley` immediately after submission.
  - `RUN-089` completed, and its outputs have now been collected with `labproj` and DVC-tracked under `data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089`.
  - `RUN-089` produced the planned replicated split diagnostic outputs: `512` random calibration/test graphlet splits plus the original `RUN-079` split, `24` graphlets stratified by graphlet dimension, `144` candidate parameter rows, split assignments, selected distances, family summaries, pairwise deltas, component summaries, and selection counts.
  - The no-module `RUN-079` graphlet-distance improvement is stable against fitted-current and smooth-core controls across replicated splits: no-module median held-out distance is `0.4839`, versus `0.5232` fitted current and `0.5358` smooth core, with win rates of `93.8%` and `97.7%`.
  - The mechanism-specific interpretation fails. The uncoupled-module control is better overall (`0.4469` median held-out distance), is best in `51.2%` of splits, and beats no-module in `63.5%` of random splits; random, shuffled, and empirical module families also beat fitted current frequently.
  - The secondary component package remains mismatched: no-module still underestimates accessibility, local maxima, and lag-1 autocorrelation and overestimates WHT high-order fraction. The original `RUN-079` no-module split distance (`0.3137`) was unusually favorable compared with the replicated median (`0.4839`).
  - `ANA-035` records the `RUN-089` review and `RES-035` records the result: the no-module graphlet-distance gain is stable against current/smooth controls, but specificity controls fail, so the result weakens the current local non-diffuse mechanism interpretation of `HYP-008` while preserving only a narrow capacity claim.
  - `EXP-036` has now been created under `HYP-008` as a discriminating ABC follow-up to `RES-035`, not a new operator. It freezes the `RUN-077` PHOT base, rebuilds the same PHOT graphlet/module context, and fits the RUN-079-style no-module, empirical-module, random-module, shuffled-module, and uncoupled-module families independently with SMC-ABC.
  - The registered runner is `scripts/phot_coupled_family_discriminating_abc.py`. Its ABC objective uses calibration graphlet aggregate distance plus component guardrails for SE/RSE, local maxima, accessibility, lag-1 autocorrelation, WHT high-order fraction, and GFT energy bands; held-out graphlets remain reserved for the scientific comparison.
  - Local verification before submission passed by byte-compiling the new runner and running `LABPROJ_QUICK=1` with temporary outputs. The quick smoke emitted the expected posterior, distance, selected-family, component-guardrail, graphlet, module, and locality tables; this is execution evidence only, not scientific support.
  - The dry-run submission for `EXP-036` is `RUN-090`. The live submission is `RUN-091`, submitted to `lab-slurm` as Slurm job `243` on `2026-05-17`; the registered request is `4` CPUs, `32` GB memory, and `04:00:00` walltime. A direct `squeue` check showed job `243` running on `ripley`.
  - `RUN-091` completed successfully, and its outputs have now been collected with `labproj`, DVC-tracked, and pushed under `data/processed/phot-coupled-family-discriminating-abc/RUN-091`. The DVC pointer is `data/processed/phot-coupled-family-discriminating-abc/RUN-091.dvc` (`bd1334082ba0a381c3bec0afc8f5a1a1.dir`, 29 files, 114,905,657 bytes).
  - `RUN-091` produced the planned PHOT discriminating ABC output set: `24` complete graphlets split into `12` calibration and `12` held-out test graphlets, fitted no-module, empirical-module, random-module, shuffled-module, and uncoupled-module ABC families, plus selected-family distances, pairwise control deltas, component guardrails, raw graphlet metrics, posterior summaries, and module-energy locality tables.
  - The ABC empirical-module family improves the fixed held-out graphlet distance over fitted current (`0.2553` versus `0.4016`), smooth core (`0.4119`), no-module (`0.3489`), and shuffled-module (`0.3535`). This confirms that ABC fitting can move the PHOT graphlet metric more effectively than the previous grid-search selection on this split.
  - The mechanism-specific interpretation still fails. Random-module is effectively tied and fractionally better (`0.255188` versus `0.255286`), the uncoupled-module control is clearly better (`0.2237`), and the selected empirical-module energy is diffuse rather than locally concentrated (`9.11/12` effective plaquette modules, normalized entropy `0.937`).
  - `ANA-036` records the `RUN-091` analysis and `RES-036` records the result: ABC fitting improves the fixed held-out PHOT graphlet distance but fails random/uncoupled specificity and module-locality checks, so it weakens the current empirical local non-diffuse `HYP-008` mechanism while preserving only a narrow capacity claim.
  - `EXP-037` has now been created under `HYP-008` as a localized-ruggedness energy diagnostic rather than another operator or refit. It reuses the completed `RUN-091` ABC fitted parameter sets and the frozen `RUN-077` PHOT graphlet/source context.
  - The registered runner is `scripts/phot_localized_ruggedness_energy_diagnostic.py`. It rebuilds dense PHOT graphlets, evaluates the empirical PHOT landscape, empirical readout residual, fitted current `a(x)`, fitted-current posterior mean, smooth core, residual-only, GNK local null, and RUN-091 no-module/empirical-module/random-module/shuffled-module/uncoupled-module candidates.
  - The primary diagnostics are edge Dirichlet energy `(f(v)-f(u))^2`, adjacent-edge energy Pearson autocorrelation and Moran I on the hypercube line graph, edge/node/site/plaquette energy concentration, and empirical/random/shuffled/uncoupled module enrichment. ABC candidates are reported both as total landscapes and innovation-from-base fields.
  - Local verification before submission passed by byte-compiling the new runner and running `LABPROJ_QUICK=1` with temporary outputs. The quick smoke emitted the planned preflight, graphlet localization, edge energy, module enrichment, selected-vs-empirical delta, rank, and reused-parameter tables.
  - The dry-run submission for `EXP-037` is `RUN-092`. The live submission is `RUN-093`, submitted to `lab-slurm` as Slurm job `244` on `2026-05-17`; the registered request is `2` CPUs, `16` GB memory, and `01:00:00` walltime. A direct `squeue` check showed job `244` running on `dutch`.
  - `RUN-093` completed successfully with inferred Slurm exit code `0`. Its outputs have now been collected with `labproj`, DVC-tracked, and backed up under `data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093`. The DVC pointer is `data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093.dvc` (`8747ed319b4f1a5fed2230d543cc5d5b.dir`, 12 files, 175,474,957 bytes).
  - `RUN-093` produced the planned PHOT localized-ruggedness output set: 24 graphlets split into 12 calibration and 12 held-out test graphlets, 80 reused `RUN-091` ABC parameter sets, edge Dirichlet-energy rows, graphlet localization summaries, module enrichment summaries, selected-vs-empirical localization deltas, and metric-rank tables.
  - The empirical PHOT graphlet layer is locally more edge/site concentrated than the fitted current/smooth core: held-out edge top-10 energy fraction is `0.3904` versus `0.3085`, site top-1 energy fraction is `0.3863` versus `0.2645`, and edge effective fraction is `0.4085` versus `0.4558`. The empirical readout residual also has a local signal, with adjacent-edge Moran I `0.0505` and edge effective fraction `0.3903`.
  - The empirical-module mechanism is not supported. Empirical PHOT energy is not enriched in empirical calibration modules (`0.8624` edge-site enrichment and `0.9495` plaquette-pair enrichment), and the selected empirical-module ABC candidate is not uniquely closest to empirical localization. No-module is closest for edge Moran I and edge effective fraction; random-module is closest for edge top-10 and site effective fraction; uncoupled and shuffled innovations are competitive or stronger on module-enrichment controls.
  - `ANA-037` records the `RUN-093` analysis and `RES-037` records the result: PHOT has residual localized ruggedness that the current fitted core under-represents, but the tested empirical-module mechanism does not explain it. This weakens the current operational form of `HYP-008` while preserving a narrower residual-localization diagnosis.
  - `EXP-038` has now been created under `HYP-008` as the direct ABC-objective follow-up to `RES-036` and `RES-037`. It keeps the same PHOT base, graphlet split, and no-module/empirical-module/random-module/shuffled-module/uncoupled-module family set as `RUN-091`, but adds local Dirichlet-energy diffuseness metrics to the calibration objective.
  - The shared runner `scripts/phot_coupled_family_discriminating_abc.py` now computes adjacent-edge Dirichlet-energy Moran I/Pearson autocorrelation plus edge, node, site, and plaquette concentration summaries for every graphlet metric row. `EXP-038` includes edge-energy Moran I, edge top-10/effective fractions, site top-1/effective fractions, and plaquette top-1/effective fractions in the ABC distance vector.
  - Local verification before submission passed by byte-compiling the modified runner and running `LABPROJ_QUICK=1` with the project `.venv`. The quick smoke emitted the expected compact output set and confirmed that `preflight.json`, `graphlet_metrics.csv`, and `posterior_particles.csv` contain the local ruggedness objective components. The first live staging attempt, `RUN-095`, was aborted before scheduler-job creation and marked failed to avoid recording a non-submitted job as running.
  - The dry-run submission for `EXP-038` is `RUN-094`. The live submitted execution is `RUN-096`, submitted to `lab-slurm` as Slurm job `246` on `2026-05-17`; the registered request is `4` CPUs, `32` GB memory, and `04:00:00` walltime. A direct `squeue` check showed job `246` running on `dutch`.
  - `RUN-096` completed successfully with inferred exit code `0`; its outputs have now been collected with `labproj`, DVC-tracked, and backed up under `data/processed/phot-abc-local-ruggedness-objective/RUN-096`. The DVC pointer is `data/processed/phot-abc-local-ruggedness-objective/RUN-096.dvc` (`4cd3edc3c50e33ecde7c5c7e32e9ec2a.dir`, 29 files, 115,223,180 bytes).
  - `RUN-096` produced the planned PHOT local-ruggedness-objective ABC output set: `24` graphlets split into `12` calibration and `12` held-out test graphlets, five independently fitted ABC families, seven added local Dirichlet-energy objective metrics, selected-family distances, pairwise deltas, component guardrails, graphlet summaries, posterior summaries, and module-energy summaries.
  - The added local ruggedness objective improves selected fitted families over fitted current and smooth core on the expanded held-out objective: no-module is `0.2771`, uncoupled-module `0.2912`, empirical-module `0.3653`, random-module `0.3668`, shuffled-module `0.3795`, fitted current `0.5210`, and smooth core `0.5268`. The local-only component package also improves relative to fitted current/smooth core (`0.6746`), with uncoupled-module `0.2900` and no-module `0.3268`.
  - The mechanism-specific result is negative. The calibration-selected empirical-module candidate loses to no-module by `0.0882` and uncoupled-module by `0.0741` on held-out full objective distance, has worse Local7 distance (`0.4623`) than no-module and uncoupled controls, and remains diffuse across empirical modules (`9.34/12` effective plaquette modules, normalized entropy `0.944`). A test-oracle empirical-module particle has held-out distance `0.2591`, but it was not calibration-selected and is therefore capacity evidence only.
  - `ANA-038` records the `RUN-096` analysis and `RES-038` records the result: adding local ruggedness statistics to the ABC objective gives a real fitted-family capacity signal over fitted current/smooth core, but it does not validate the empirical local-module mechanism required by `HYP-008`.
  - `CLAIM-003` has been created and completed to record the positive scoped interpretation of `RES-035` through `RES-038`: structured high-order residual or auxiliary capacity improves PHOT held-out graphlet/local-ruggedness benchmarks over the current fitted model. The claim explicitly does not assert residue-level identification of true epistatic interaction terms; DCA/module/mutation-feature routing is treated as an interpretable scaffold for representing learnable high-order ruggedness signal in benchmarks and downstream supervised models.
  - `EXP-039` has been created under `HYP-008` as a nested soft-min/Boltzmann ensemble capacity diagnostic. It freezes the first soft-min layer to the RUN-084 rank-1 parameterisation, then fits a third latent residual landscape through `c(x)=g(a(x),b(x)); f(x)=h(c(x),d(x))` with independently parameterised second-stage no-module, empirical-module, random-module, shuffled-module, and uncoupled-module ABC families.
  - The registered runner is `scripts/phot_nested_softmin_ensemble_abc.py`. It reuses the frozen RUN-077 PHOT graphlet context, fits only on calibration graphlets, evaluates held-out graphlet/local-ruggedness distance, writes posterior/round summaries, selected-family and pairwise-control distance tables, module-energy summaries, base-correlation diagnostics, and a locked first-stage transform report.
  - Local verification before submission passed by byte-compiling the new runner and running `LABPROJ_QUICK=1` with temporary outputs. The quick smoke emitted the expected compact output set and confirmed the first-stage transform is fit once over the experiment graphlet-node pool rather than per graphlet.
  - The dry-run submission for `EXP-039` is `RUN-097`. The live submitted execution is `RUN-098`, submitted to `lab-slurm` as Slurm job `247` on `2026-05-18`; the registered request is `4` CPUs, `32` GB memory, and `06:00:00` walltime. `RUN-098` later completed with inferred Slurm exit code `0`; its outputs have now been collected with `labproj`, DVC-tracked, and backed up under `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098`. The DVC pointer is `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc` (`e1e92177585cf2d09ed03bb3495b6ae4.dir`, 30 files, 206,813,520 bytes).
  - `RUN-098` produced the planned nested soft-min ensemble output set: `24` complete PHOT graphlets split into `12` calibration and `12` held-out test graphlets, five independently fitted second-stage ABC families, selected-family distance tables, pairwise deltas, component guardrails, posterior summaries, module-energy locality summaries, base-correlation diagnostics, and the locked first-stage transform report.
  - The nested empirical-module mechanism is negative. The calibration-selected empirical-module candidate has held-out full-objective distance `0.5432`, worse than fitted current (`0.5210`), smooth core (`0.5268`), locked first-stage (`0.5296`), no-module (`0.5372`), random-module (`0.4826`), and shuffled-module (`0.5110`), beating only uncoupled-module (`0.6401`). The best calibration-selected fitted family is the random-module negative control (`0.4826`).
  - The local-ruggedness and locality diagnostics also fail specificity. Local7 distance is `0.5250` for empirical-module, worse than shuffled-module (`0.4199`) and no-module (`0.4766`), and the selected empirical-module energy is diffuse (`10.56/12` effective plaquette modules, normalized entropy `0.972`). Test-oracle posterior particles show unused capacity (`0.3953` no-module and `0.4092` empirical-module), but those are not calibration-selected support.
  - `ANA-039` records the `RUN-098` analysis and `RES-039` records the result: nested third-landscape soft-min capacity does not rescue the empirical local-module implementation of `HYP-008`; it further supports treating prior improvements as generic capacity/selection effects rather than residue-local mechanism evidence.
