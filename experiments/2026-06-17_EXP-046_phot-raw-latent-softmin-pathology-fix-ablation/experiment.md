# EXP-046: PHOT raw-latent soft-min pathology fix ablation

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Do the three RES-044 path-forward fixes, phased in cumulatively, repair PHOT raw-latent accessibility and scale pathologies while retaining the working A/C/D observable calibration protocol?

## Pre-experiment prediction

If the dominant failure is latent geometry, adding a low-dimensional soft-min latent should improve raw accessibility relative to EXP-045, raw-latent constraints should prevent collapse or explosion, and making raw accessibility first-class should isolate whether accessibility can be fitted without losing magnitude/ratio calibration.

## Rationale

RES-044 showed that A/C/D improved observable magnitude and ratio but left raw accessibility roughly an order of magnitude too low, with the monotone readout unable to repair topology. This experiment keeps the PHOT source graphlets, RUN-098 comparator, and recent deterministic numerical fast-fit protocol fixed while cumulatively testing the three RES-044 path-forward points.

## Experimental design

- Reuse the EXP-045 PHOT fast-fit inputs and output schema so the comparison is
  paired to `RUN-117`: PHOT graphlets and prepared assay from `RUN-077`, plus
  the `RUN-098` full-ABC comparator table.
- Run three cumulative RES-044 path-forward stages in a single experiment:
  - `1`: replace random specific-`J` ruggedness with a low-dimensional
    two-additive-trait `soft_min_values` latent generator while retaining the
    working EXP-045 A/C/D observable metric conditioning.
  - `1+2`: add raw-latent objective constraints on raw single-mutant SD, raw
    local-maxima count, and raw accessibility, and restrict the fitted
    observation readout to a calibration-fitted two-parameter monotone sigmoid.
  - `1+2+3`: keep `1+2` and up-weight raw accessibility so it becomes a
    first-class calibration target rather than a diagnostic afterthought.
- Fit each stage against `PHOT_CHLRE_Chen_2023` only, using the same
  deterministic/random numerical fast search with one local-refinement round
  used by recent PHOT fast-fit experiments.
- Select the top candidates within each cumulative stage, not globally across
  stages, so the run yields three interpretable result sets.
- Report observable and raw latent metrics together, including raw
  accessibility and raw local maxima.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
- `docs/RES-044_epistasis_calibration_followup.md`
- `src/project_adapt_env/coupled_landscape.py`
- `scripts/phot_trait_space_bounded_gate_fast_fit.py`

## Configuration

- `config.yaml` defines the frozen PHOT source, RUN-098 comparator, graphlet
  summaries, cumulative stage overrides, fast-fit search budgets, raw-latent
  objective terms, and quick smoke-test override.

## Execution

- Run through `labproj submit EXP-046 --target lab-slurm --execute`.
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

- Create an ANA record after execution comparing all three cumulative stages.
- The decisive readouts are whether raw accessibility moves toward empirical
  PHOT while raw single-mutant SD, raw local maxima, observable magnitude,
  magnitude-to-SD ratio, and held-out graphlet distance remain interpretable.
- Treat an improvement in observable readout without raw accessibility/local
  topology improvement as another readout fit, not a validated procedural
  landscape fix.
- Keep the interpretation PHOT-local until the same calibration behavior is
  shown on a second assay.

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

- The implementation follows the advisory
  `docs/RES-044_epistasis_calibration_followup.md` as a diagnostic test, not as
  established evidence.
- Local quick smoke completed under `scratch/exp046_quick/outputs`, with three
  staged selected rows (`1`, `1+2`, `1+2+3`) and raw-latent objective columns
  present in `bounded_gate_candidates.csv`.
- Focused regression tests passed:
  `.venv/bin/pytest tests/test_calibration_modes.py tests/test_phot_trait_space_bounded_gate_fast_fit.py tests/test_phot_nested_fit_protocol_benchmark.py -q`.
- `RUN-118` is the dry-run record created by
  `labproj submit EXP-046 --target lab-slurm --dry-run`.
- `RUN-119` is the live submission created by
  `labproj submit EXP-046 --target lab-slurm --execute`. It was submitted to
  `lab-slurm` as Slurm job `1917` with `8` CPUs, `32G` memory, and a
  `12:00:00` walltime.
- `labproj refresh --write --strict` updated `RUN-119` to `completed`; outputs
  remain on the remote runtime volume.
- `ANA-045` copied the compact RUN-119 output tables and derived stage-level
  equivalence, raw-latent pathology, component, and collection-audit summaries.
  The large WHT/GFT coefficient tables remain on the runtime volume.
- `RES-045` records the result: EXP-046 produced `0/24` operationally
  equivalent rows. Stage `1` was the least bad stage but had held-out distance
  `1.2240` and raw scale hundreds-fold too large; stages `1+2` and `1+2+3`
  worsened held-out fit and collapsed observable variation.
- No new DVC target was promoted for `RUN-119`; the compact decision-relevant
  JSON/CSV outputs are copied and hash-audited inside `ANA-045`, and full
  outputs remain under the remote runtime path.
