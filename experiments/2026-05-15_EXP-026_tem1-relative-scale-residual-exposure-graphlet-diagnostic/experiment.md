# EXP-026: TEM-1 relative-scale residual exposure graphlet diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005
- Previous experiment: EXP-025
- Motivating result: RES-025

## Question

Do fixed-scale residual and high-order components become more topologically relevant as adaptive-walk diminishing returns reduce the smooth thermodynamic/additive gain, when dense graphlets are balanced across start classes and walk stages?

## Pre-experiment prediction

If the relative-scale mechanism is correct, absolute residual high-order variance should remain controlled while residual exposure relative to available beneficial effects or base order-1 variance increases in later walk stages, accompanied by increased WHT high-order fraction, SE/RSE, local maxima, and reduced accessibility compared with the smooth core. Local-DFE-scaled residuals and placement controls test whether this is merely strength scaling rather than a biologically meaningful procedural mechanism.

## Rationale

RUN-065 showed that aggregate residuals can create local topology but failed the multi-start sampling design and did not make the relative-scale causal chain explicit. This scaled follow-up fixes graphlet-cap sampling, keeps outputs small, and directly records whether residual/high-order terms are unmasked as smooth adaptive gains diminish.

## Experimental design

- Use the same TEM-1 reference and MMSeqs/ColabFold alignment source as
  `EXP-025`, but cap the panel at one landscape seed, one start per start
  class, three graphlet dimensions (`4`, `5`, `6`), two contact-source
  families, and at most `90` graphlets.
- Fix the `RUN-065` sampling flaw by ordering graphlet backgrounds with
  balanced round-robin sampling over both start class and adaptive-walk stage
  before applying the graphlet cap.
- Run smooth-core adaptive walks from the five existing start classes:
  reference/consensus-like, high-fitness close homologue, marginal-stability
  homologue, distant low-fitness homologue, and viable off-consensus.
- Sample graphlet backgrounds from early, middle, and final/late walk stages.
  The run does not assert true global distance to peak; it records walk step,
  distance to the capped walk endpoint, available beneficial fraction, mean
  available beneficial effect, and best available effect as observable
  diminishing-return proxies.
- Compare two residual scaling modes:
  - `constant_walk_initial_beneficial_effect`, where residual coefficients are
    scaled from the initial walk's beneficial-effect scale and therefore do
    not automatically shrink with local DFE changes;
  - `local_dfe_sd`, the `RUN-065`-style control where residual strength scales
    with each local graphlet/DFE.
- Record explicit relative-scale readouts: base WHT variance by order,
  residual WHT variance by order, residual high-order absolute variance,
  residual-to-base order-1 variance ratio, residual-to-available-beneficial
  effect ratio, and total high-order-to-base order-1 variance ratio.
- Keep output volume small by recording only top Walsh/GFT terms per
  graphlet/model while retaining full order-level variance summaries in
  metric rows.

## Inputs

- `config.yaml`
- runner: `scripts/tem1_aggregate_graphlet_residual_ruggedness.py`
- `external/Adapt-Env`
- `external/landscapy`
- TEM-1 MMSeqs/ColabFold alignment fetched or cached through the runner.

## Configuration

- `config.yaml` contains the complete reduced design and quick-mode smoke
  overlay.

## Execution

- Use `labproj submit EXP-026 --target lab-slurm --dry-run` before live
  submission, then `--execute` for the registered run.

## Expected outputs

- tem1_mmseqs_alignment.fasta
- adaptive_walk_starts.csv
- adaptive_walk_trajectories.csv
- adaptive_walk_graphlet_backgrounds.csv
- start_candidate_report.csv
- summary.json
- module_parameter_sweep.csv
- graphlet_catalog.csv
- calibration_graphlet_metrics.csv
- test_graphlet_metrics.csv
- graphlet_walsh_terms.csv
- graphlet_gft_coefficients.csv
- graphlet_sign_epistasis_pairs.csv
- graphlet_accessibility.csv
- calibration_target_summary.json
- heldout_test_summary.json

## Analysis plan

- Check that graphlets span all requested start classes and walk-stage classes
  under the cap.
- Stratify the selected primary residual and controls by walk-stage class,
  available beneficial fraction, mean available beneficial effect, best
  available effect, and distance to capped walk endpoint.
- Test the causal chain explicitly: the relative-scale model is supported only
  if residual high-order absolute variance is controlled while residual
  exposure and total WHT/SE/RSE/topology readouts rise as smooth adaptive gains
  diminish.
- Compare against smooth core, exact-match RES-024 control, DCA pair override,
  shuffled-contact residual, random-position residual, GNK local null, and
  residual-only embedding diagnostic.
- Treat a result as negative if topology rises only because absolute residual
  variance increases, if local-DFE scaling explains the same pattern, or if
  shuffled/random placement matches the primary contact residual.

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
