# EXP-044: PHOT bounded trait-space current-model fast fit

## Status

- completed; reviewed in `ANA-043` and `RES-043`

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Does a minimal bounded trait-space epistasis route in the current AdaptEnv model, fitted by a numerical fast fitter, produce PHOT graphlet accessibility, epistasis, spectral, and local-energy summaries comparable to empirical PHOT and the previously computed full-ABC PHOT comparator?

## Pre-experiment prediction

Bounding epistatic perturbations before they enter stability should reduce the gross epistasis-amplitude pathology seen in residual min-coupled fits; success still requires held-out accessibility, sign/reciprocal/magnitude epistasis, diminishing returns, graphlet spectra, and local-energy concentration to fall near empirical PHOT and not materially worse than the full-ABC comparator.

## Rationale

EXP-042/EXP-043 showed that cheaper fitting of the nested residual-min family did not rescue PHOT component failures, suggesting the failure is in the model route rather than only in optimizer cost. This experiment tests the smallest model-level correction: keep current AdaptEnv structure, but bound epistatic signals in trait space before the readout, then fit by deterministic/random numerical search rather than broad ABC.

## Experimental design

- Patch the current AdaptEnv model so epistatic signals routed to causal trait
  space can be symmetrically bounded before they perturb stability/function.
  Historical behavior is preserved when `trait_space_epistasis_bound` is unset.
- Reuse RUN-077 as the frozen PHOT prepared-assay and graphlet source.
- Rebuild actual current-model PHOT landscapes under a deterministic/random
  numerical search over stability, epistasis, BLOSUM, and bounded trait-space
  parameters. This is not a nested residual-min graphlet model and is not broad
  SMC-ABC.
- Fit a positive affine observable readout on calibration graphlets so distances
  are computed on the empirical PHOT `latent_phi` scale; retain `raw_*`
  summaries to expose whether any apparent improvement is only readout scaling.
- Evaluate selected candidates on calibration/test/all graphlet splits for
  accessibility, local maxima, random-walk autocorrelation, Walsh/GFT spectral
  mass, local edge/site/plaquette energy concentration, magnitude/sign/
  reciprocal epistasis, and diminishing returns.
- Compare selected bounded current-model candidates with the RUN-098 full-ABC
  PHOT comparator distance table. Because this is a different model family,
  equivalence is judged against the best available full-ABC comparator, not
  parameter identity.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
- `external/Adapt-Env/adaptenv/landscape.py`
- `scripts/phot_trait_space_bounded_gate_fast_fit.py`
- `scripts/phot_coupled_family_discriminating_abc.py`

## Configuration

- `config.yaml` defines the frozen PHOT source, full-ABC comparator, graphlet
  summaries, guardrails, bounded trait-space parameter search, and quick-mode
  smoke-test override.

## Execution

- Run through `labproj submit EXP-044 --target lab-slurm --execute`.
- `run.sh` pulls RUN-077/RUN-098 DVC targets if they are missing and then runs
  `scripts/phot_trait_space_bounded_gate_fast_fit.py`.

## Expected outputs

- summary.json
- preflight.json
- bounded_gate_candidates.csv
- selected_bounded_gate_candidates.csv
- bounded_gate_distances.csv
- bounded_gate_equivalence.csv
- graphlet_metrics.csv
- graphlet_summary.csv
- walsh_terms.csv
- gft_coefficients.csv
- abc_comparator_selected_family_distances.csv

## Analysis plan

- Create an ANA record after execution comparing calibrated and raw epistasis
  amplitude, accessibility, neutral/path structure proxies, local spectral
  summaries, local-energy concentration, and held-out distance against empirical
  PHOT and RUN-098 full ABC.
- Treat a calibrated-only improvement without raw-scale/local-energy improvement
  as an observation/readout fix, not a model-level procedural rescue.
- Do not copy large runtime artifacts unless they are required for the ANA;
  compact summaries and exact runtime paths are sufficient for the first pass.

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [x] Outputs generated
- [ ] Important outputs DVC-tracked
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- Focused regression coverage before submission:
  `.venv/bin/pytest tests/test_calibration_modes.py tests/test_phot_trait_space_bounded_gate_fast_fit.py -q`.
- Local quick smoke completed before submission with `LABPROJ_QUICK=1`, writing
  disposable outputs under `scratch/exp044_quick/outputs`. The smoke exercised
  graphlet construction, numerical candidate evaluation, metric summaries, and
  comparator-table loading.
- `RUN-113` is the dry-run record created by
  `labproj submit EXP-044 --target lab-slurm --dry-run`.
- `RUN-114` was an interrupted live staging attempt. It had no scheduler job id
  and was marked failed after interruption; the interruption avoided continuing
  a potentially oversized source-staging transfer.
- `.gitignore` now excludes local `.dvc/cache`, `.dvc/tmp`, `.pytest_cache`,
  and `scratch/` directories to reduce future staged payloads.
- `RUN-115` is the live submission created by
  `labproj submit EXP-044 --target lab-slurm --execute`. It was submitted to
  `lab-slurm` as Slurm job `1909` with `8` CPUs, `32G` memory, and a `12:00:00`
  walltime. The runtime output directory is
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs`.
- `RUN-115` completed with inferred Slurm exit code `0`. `ANA-043` copied the
  compact source tables needed for interpretation and left `gft_coefficients.csv`,
  `walsh_terms.csv`, and the PHOT alignment on the runtime volume. The compact
  copied payload is audited in
  `analyses/2026-06-16_ANA-043_run-115-phot-bounded-trait-space-fast-fit-review/tables/run-115_output_collection_audit.md`.
- `RES-043` records the result: bounded trait-space current-model fitting
  improved the previous component-guardrail failure mode, but produced `0/12`
  operationally equivalent rows against the best available `RUN-098` full-ABC
  PHOT comparator. The best held-out bounded candidate had test distance
  `0.4269` versus the best full-ABC context row at `0.2115`; raw latent
  magnitude epistasis remained tens to hundreds of units too large.
- DVC tracking of full runtime outputs was intentionally not performed during
  this first-pass review. The relevant compact tables are retained in
  `ANA-043`; full coefficient-level outputs remain at the RUN-115 runtime path.
