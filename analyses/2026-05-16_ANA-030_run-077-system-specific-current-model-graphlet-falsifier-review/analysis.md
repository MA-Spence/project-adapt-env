# ANA-030: RUN-077 system-specific current-model graphlet falsifier review

## Purpose

Assess whether system-specific current-model PHOT and SPTN1 refits overturn the RES-028 graphlet falsifiers after correcting the RUN-073 pooling and stale empirical-comparison problem.

## Linked experiments/runs

- Experiments: `EXP-030`
- Runs: `RUN-077`
- Source data:
  `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077`
- DVC pointer:
  `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-077_preflight_and_output_audit.md`
- `tables/run-077_graphlet_catalog_summary.md`
- `tables/run-077_fit_summary.md`
- `tables/run-077_branch_validation.md`
- `tables/run-077_empirical_primary_graphlet_metrics.md`
- `tables/run-077_test_model_graphlet_metrics.md`
- `tables/run-077_model_distances.md`
- `tables/run-077_falsifier_summary.md`

## Run and data audit

`RUN-077` completed successfully with scheduler exit code `0`. The collected
output directory contains 55 files, including `summary.json`,
`falsification_tests.json`, system-specific ABC summaries, branch validation
tables, dense graphlet metrics, Walsh terms, and GFT coefficients. The output
directory is DVC-tracked.

The current-model preflight passed. The run used two independent
system-specific landscapes for `PHOT_CHLRE_Chen_2023` and
`SPTN1_CHICK_Tsuboyama_2023_1TUD`, observed fitness was combined as
`stability_margin`, generic epistasis was routed to stability, and the fitted
parameters were limited to the current stability and residual-epistasis
parameter set. This directly fixes the central procedural defect in `RUN-073`:
there is no pooled PHOT+SPTN1 posterior and no reliance on legacy cached
epistasis definitions.

## Scope relative to previous falsifiers

This analysis treats `RUN-077` as a direct retest of the relevant break points
from `RES-028`, conditioned by the procedural limitation documented in
`RES-029`.

The six falsifiable points are assessed as follows:

- Smooth-core-only empirical fit: partly unresolved. In PHOT, the current model
  is closest on test graphlets, but only slightly closer than the smooth-core
  control (`0.4016` versus `0.4119` distance). That is not strong evidence that
  the residual/current epistasis layer is doing indispensable work. In SPTN1,
  smooth-core rejection is formally passed, but the retained graphlets are only
  2D and topologically trivial.
- No residual-scale effect: not directly retested as a scale grid in this run.
  `RUN-077` compares fitted current, smooth-core-only, residual-only, and GNK
  controls, but it does not replace the dedicated residual-scale evidence from
  `RUN-069`/`RES-027` and `RUN-071`/`RES-028`.
- Wrong spectral envelope: reproduced for PHOT. Held-out PHOT graphlets have
  WHT high-order fraction `0.3788`, above the preregistered non-dominant-tail
  ceiling used by the falsifier. SPTN1 is not testable for this point because
  only 2D graphlets survived the dense-graphlet construction.
- Control equivalence: GNK and residual-only controls do not match the fitted
  current model on test distances. For PHOT, GNK distance is `0.7594` and
  residual-only distance is `2.6756`, compared with `0.4016` for the fitted
  current model. Smooth-core equivalence remains the main control concern.
- Wrong evolutionary dynamics: not reproduced for PHOT. The fitted current
  model is the best held-out PHOT candidate across the joint graphlet package
  of spectra, SE/RSE, local maxima, accessibility, and random-walk
  autocorrelation. The SPTN1 runner flags this falsifier because the posterior
  mean branch beats the best particle on the test metric, but that should be
  interpreted as branch-selection instability in a 2D smooth case, not as
  evidence for missing high-order ruggedness.
- Scale confounding: not supported for PHOT. The high-order signal is retained
  after the nonlinear global-readout correction. SPTN1 is not informative for
  this control because the high-order signal is undefined/absent at graphlet
  dimension 2.

## Key observations

- The graphlet panel contains 32 dense graphlets: 24 PHOT graphlets across
  dimensions 4, 5, and 6, and 8 SPTN1 graphlets at dimension 2.
- The primary PHOT empirical graphlets remain rugged. On held-out graphlets,
  sign epistasis fraction is `0.6146`, reciprocal sign epistasis fraction is
  `0.2583`, local maxima average `5.17`, accessible path fraction is `0.1975`,
  WHT high-order fraction is `0.3788`, and GFT energy is split
  low/mid/high as `0.6212`/`0.2236`/`0.1552`.
- The PHOT fitted current model is the closest held-out candidate:
  `0.4016` distance versus `0.4119` for smooth core, `0.4408` for posterior
  mean, `0.7594` for GNK, and `2.6756` for residual-only.
- The small PHOT separation between current and smooth-core controls is the
  main scientific weakness. This is not a clean validation of the residual
  epistasis mechanism; it is a partial rescue of the evolutionary/topological
  graphlet package after system-specific refitting.
- The SPTN1 graphlets are essentially smooth in the retained local Hamming
  graphs: SE and RSE are both `0`, local maxima average `1`, and WHT
  high-order fraction is `0` by construction. These outputs are useful as a
  negative/smooth system control but cannot test the high-order spectral
  claims.
- The fitted empirical branch validation remains mixed. PHOT has holdout
  Spearman `0.5033` and double-holdout Spearman `0.7180`, but epistasis
  prediction Spearman is only `0.0577`. SPTN1 has holdout Spearman `0.5539`,
  double-holdout Spearman `0.4888`, and epistasis prediction Spearman
  `0.5264`.

## Interpretation

`RUN-077` partially overturns the strongest negative reading of `RES-028`.
Under current semantics, with PHOT and SPTN1 fit separately, the embedded
current model is no longer beaten by GNK or residual-only controls on held-out
PHOT graphlet package distance. That matters because it shows that the stale
cache and pooled-calibration objections were not cosmetic.

The result does not validate the model as a general graphlet benchmark
generator. The PHOT empirical spectral-envelope falsifier still fires, and the
fitted current model is only marginally better than smooth core. Thus the
supported conclusion is narrower: system-specific current refitting repairs
part of the topological/dynamical mismatch, but the present smooth-core plus
diffuse residual structure still does not explain the empirical high-order
spectral envelope convincingly.

This strengthens the motivation for `HYP-008`. If local ruggedness in empirical
protein landscapes is concentrated around particular sequence, structural, or
trait regions, then a diffuse residual layer may fit some aggregate topology
without reproducing the empirical localization and spectral structure that the
PHOT graphlets expose.

## Limitations

- The SPTN1 dense graphlets are only 2D, so they cannot test the high-order
  WHT/GFT part of the falsifier.
- `RUN-077` is not a residual-scale sweep. It should not be used as the primary
  evidence for monotonic residual-scale effects.
- The PHOT empirical layer is still only a fitted MAVE-NN/latent readout over
  ProteinGym data, not direct exhaustive assay truth for every sampled
  graphlet.
- The current-vs-smooth-core PHOT separation is small enough that conclusions
  about the residual mechanism should remain cautious.
