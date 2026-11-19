# EXP-045: PHOT cumulative epistasis calibration fix ablation

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Do the RES-043 epistasis-calibration fixes A through D, added cumulatively, resolve the PHOT magnitude-epistasis/accessibility tradeoff under the existing numerical fast-fit protocol?

## Pre-experiment prediction

Fix A should test whether a monotone observation link alone can raise observable magnitude epistasis without collapsing accessibility; if not, the cumulative A+B, A+B+C, and A+B+C+D stages should show whether global-nonlinearity sourcing, an explicit specific-epistasis variance-ratio knob, and scale-invariant/distribution-conditioned objectives are needed.

## Rationale

RES-043 showed that bounded trait-space routing improved component guardrails but left raw latent epistasis far too large and accessibility unstable. The scratch RES-043 calibration review argues that the defect is a combination of affine readout ratio invariance, sourcing magnitude epistasis from specific J, lack of an explicit epistasis/additive variance-ratio target, and poor objective conditioning around near-zero signed means. This experiment tests those fixes cumulatively on PHOT before any broader scale-up.

## Experimental design

- Reuse the EXP-044 PHOT fast-fit runner and frozen inputs so the comparison is
  paired to `RUN-115`: PHOT graphlets and prepared assay from `RUN-077`, plus
  the `RUN-098` full-ABC comparator table.
- Run four cumulative fit stages in a single experiment:
  - `A`: replace the positive-affine observation calibration with a fitted
    monotone-piecewise link.
  - `A+B`: add a global sigmoid stability-margin readout so magnitude and
    diminishing-returns epistasis can arise from a monotone global
    nonlinearity while keeping specific pairwise `J` smaller.
  - `A+B+C`: add an explicit fitted specific-epistasis/additive-scale ratio
    knob and score the realized magnitude-epistasis/single-mutant-SD ratio.
  - `A+B+C+D`: remove the near-zero signed-pairwise-mean objective from the
    final conditioning and instead compare signed-epsilon distributional
    summaries plus scale-invariant epistasis metrics.
- Fit each stage against `PHOT_CHLRE_Chen_2023` only, using the same
  deterministic/random numerical fast search with one local-refinement round
  used by recent PHOT fast-fit experiments.
- Select the top candidates within each cumulative stage, not globally across
  stages, so the run yields four interpretable result sets.
- Report observable and raw latent metrics together, including accessibility,
  local maxima, sign/reciprocal-sign epistasis, magnitude epistasis, raw
  magnitude epistasis, magnitude-to-single-mutant-SD ratio, Walsh/GFT spectra,
  and local edge/site/plaquette energy summaries.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
- `scratch/RES-043_epistasis_calibration_review.md`
- `external/Adapt-Env/adaptenv/calibration.py`
- `external/Adapt-Env/adaptenv/landscape.py`
- `scripts/phot_trait_space_bounded_gate_fast_fit.py`

## Configuration

- `config.yaml` defines the frozen PHOT source, RUN-098 comparator, graphlet
  summaries, cumulative stage overrides, fast-fit search budgets, and quick
  smoke-test override.

## Execution

- Run through `labproj submit EXP-045 --target lab-slurm --execute`.
- `run.sh` pulls RUN-077/RUN-098 DVC targets if missing, then runs
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

- Create an ANA record after execution comparing all four cumulative stages.
- The decisive readouts are whether observable magnitude epistasis can approach
  empirical PHOT (`0.2926`) while accessibility remains near empirical PHOT
  (`0.1975`), and whether raw latent magnitude epistasis and the
  magnitude-to-single-mutant-SD ratio stop showing the RES-043 scale defect.
- Treat an improvement in observable readout without raw ratio/accessibility
  improvement as a readout fix only, not a validated procedural landscape fix.
- Keep the interpretation PHOT-local until the same calibration behavior is
  shown on a second assay.

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- The implementation follows the advisory `scratch/RES-043_epistasis_calibration_review.md`
  as a diagnostic test, not as established evidence.
- Local quick smoke completed under `scratch/exp045_quick/outputs`, with four
  staged selected rows (`A`, `A+B`, `A+B+C`, `A+B+C+D`) and the new
  scale-invariant epistasis columns present.
- Focused regression tests passed:
  `.venv/bin/pytest tests/test_calibration_modes.py tests/test_phot_trait_space_bounded_gate_fast_fit.py tests/test_phot_nested_fit_protocol_benchmark.py -q`.
- `RUN-116` is the dry-run record created by
  `labproj submit EXP-045 --target lab-slurm --dry-run`.
- `RUN-117` is the live submission created by
  `labproj submit EXP-045 --target lab-slurm --execute`. It was submitted to
  `lab-slurm` as Slurm job `1911` with `8` CPUs, `32G` memory, and a
  `12:00:00` walltime.
