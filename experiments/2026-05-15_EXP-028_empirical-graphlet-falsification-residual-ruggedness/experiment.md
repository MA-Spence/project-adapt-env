# EXP-028: Empirical graphlet falsification suite for residual ruggedness model

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005

## Question

Do the six explicit falsifiers for the smooth-core plus fixed-scale residual ruggedness mechanism hold when empirical ProteinGym graphlets, nonlinear readout controls, matched residual controls, and RUN-069-style synthetic graphlets are evaluated in one paired diagnostic suite?

## Pre-experiment prediction

If the model structure is scientifically defensible, empirical graphlets should show low-order-dominated spectra with a measurable high-order tail; smooth-core-only and additive global-readout controls should not match empirical SE/RSE/topology as well as residual models; increasing fixed residual scale should remain directional across synthetic seeds; contact residuals should outperform random-position, shuffled-contact, and GNK controls on joint empirical distance; and high-order signal should not disappear after global-readout correction.

## Rationale

RES-027 supports the fixed residual-scale knob procedurally but leaves six explicit falsifiers unresolved: smooth-core-only empirical fit, no residual-scale effect, wrong empirical spectral envelope, control equivalence, wrong evolutionary dynamics, and scale confounding. This experiment tests those breakpoints directly while keeping output columns comparable to RUN-067 and RUN-069.

## Experimental design

- Run one compact RUN-069-compatible synthetic graphlet suite with the same
  candidate model families: smooth core, RES-024 exact-match haplotype control,
  RES-023 DCA pair override, fixed-scale contact residual, shuffled-contact
  residual, random-position residual, GNK local null, and residual-only
  diagnostic.
- Add a second TEM-1 synthetic landscape seed relative to RUN-069 so the
  residual-scale directional test is not a one-seed artifact.
- Extract complete measured ProteinGym binary Hamming graphlets from cached
  prepared assay data. `PHOT_CHLRE_Chen_2023` supplies the primary 4-, 5-, and
  6-dimensional empirical graphlets; `SPTN1_CHICK_Tsuboyama_2023_1TUD` and
  `RASK_HUMAN_Weng_2022_binding-DARPin_K55` provide lower-dimensional
  stability/binding controls.
- Evaluate empirical graphlets on both `mavenn_yhat` and `latent_phi` score
  scales. The latent scale is the primary score for falsification because it is
  intended to remove part of the nonlinear global-readout effect.
- Fit an additive latent plus polynomial global readout inside every empirical
  graphlet and record how much WHT high-order fraction remains in the residual.
- Compare synthetic model/control distributions to empirical graphlet
  distributions using the same metrics as RUN-067/RUN-069: SE, RSE, local
  maxima, accessibility, random-walk autocorrelation, WHT high-order fraction,
  and GFT high-frequency fraction.

The six explicit falsifiers are tested as follows:

- Smooth-core-only empirical fit: smooth-core or additive-global-readout fits
  match empirical graphlet topology/spectra as well as the residual model.
- No residual-scale effect: the fixed residual scale grid fails to move SE,
  RSE, local maxima, WHT/GFT high-order fraction, or autocorrelation
  directionally.
- Wrong spectral envelope: empirical graphlets do not show low-order dominance
  with a measurable high-order tail.
- Control equivalence: random-position, shuffled-contact, or GNK controls match
  the empirical target as well as contact residuals.
- Wrong evolutionary dynamics: spectra match but accessibility, local maxima,
  and random-walk statistics fail.
- Scale confounding: empirical high-order WHT signal mostly disappears after
  additive-global-readout correction.

## Inputs

- Cached ProteinGym prepared assay panels from prior completed runs:
  `PHOT_CHLRE_Chen_2023`, `SPTN1_CHICK_Tsuboyama_2023_1TUD`, and
  `RASK_HUMAN_Weng_2022_binding-DARPin_K55`.
- TEM-1 reference sequence and MMSeqs alignment workflow used by EXP-027.
- RUN-067/RUN-069-compatible graphlet metric definitions from
  `scripts/tem1_aggregate_graphlet_residual_ruggedness.py`.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- falsification_tests.json
- synthetic_empirical_model_distances.csv
- synthetic_empirical_best_by_model.csv
- synthetic_graphlet_suite/module_parameter_sweep.csv
- synthetic_graphlet_suite/test_graphlet_metrics.csv
- synthetic_graphlet_suite/graphlet_walsh_terms.csv
- synthetic_graphlet_suite/graphlet_gft_coefficients.csv
- empirical_graphlet_suite/empirical_assay_audit.csv
- empirical_graphlet_suite/empirical_graphlet_metrics.csv
- empirical_graphlet_suite/empirical_graphlet_summary.csv
- empirical_graphlet_suite/empirical_global_readout_confounding.csv
- empirical_graphlet_suite/empirical_global_readout_summary.csv

## Analysis plan

- Create an ANA record after execution.
- Treat every falsifier as capable of weakening the model. In particular, a
  positive result for the scale knob is insufficient if empirical graphlets are
  smooth-core-explainable, if GNK/random controls match empirical statistics as
  well as contact residuals, or if high-order signal disappears after nonlinear
  readout correction.
- Compare this experiment's live RUN output directly to RES-026/RES-027 values for the synthetic
  fixed-scale contact residual grid.

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [x] Outputs generated
- [x] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- 2026-05-15: Implemented
  `scripts/empirical_graphlet_falsification_suite.py` and wired EXP-028
  `run.sh` to execute it. The experiment deliberately uses cached empirical
  prepared assays to avoid retraining MAVE-NN and to make the falsification
  suite directly reproducible from existing project data.
- 2026-05-15: Local quick smoke completed successfully with quick synthetic
  graphlet limits and a three-graphlet PHOT empirical panel. The smoke confirms
  execution and output shape only; it is not evidence for the full experiment.
- 2026-05-15: RUN-071 completed with scheduler exit code `0`; output was
  collected, DVC-tracked, analysed in `ANA-028`, and interpreted in `RES-028`.
  The result preserves the fixed residual-scale knob but falsifies the strong
  current empirical benchmark-adequacy claim.
