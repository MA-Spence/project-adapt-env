# EXP-043: PHOT epistasis-calibrated cheap fitting protocol

## Status

- completed; reviewed in `ANA-042` and `RES-042`

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Can a cheap PHOT nested-min fitting protocol with an explicit positive affine readout calibration and epistasis-amplitude objective remove the raw pairwise epistasis magnitude mismatch while preserving graphlet ruggedness, accessibility, and local-energy summaries?

## Pre-experiment prediction

Adding a calibrated observable readout plus pairwise epistasis-amplitude summaries to the cheap fitting objective will reduce observable-scale epistasis mismatch relative to RUN-098/RUN-110 selected fitted families, but raw latent epistasis may remain large; success requires calibrated held-out epistasis amplitude, accessibility, graphlet distance, and RES-040 local-energy failures to remain within guardrails.

## Rationale

RES-041 showed that cheap protocols reproduce the selected full-ABC fitted-family PHOT result, but the full-ABC result itself has pairwise epistasis amplitudes tens to hundreds of units versus about 0.293 empirically. This experiment keeps the nested-min family fixed and tests whether cheap fitting with an explicit readout calibration and epistasis-amplitude targets can produce a usable observable landscape before scaling to complete GB1 and TrpB landscapes.

## Experimental design

- Keep the latest PHOT nested min-coupled sequential-residual model family fixed.
  This experiment does not rerun the broad full SMC-ABC search.
- Reuse RUN-077 as the frozen PHOT prepared-assay and graphlet source, and RUN-098
  selected full-ABC particles as the comparator particles.
- Fit the same second-stage residual landscape parameterization through the cheap
  fitting protocol implemented in `scripts/phot_nested_fit_protocol_benchmark.py`.
- For each candidate, fit a positive affine observable readout on calibration
  graphlets, mapping latent nested-min outputs onto empirical PHOT `latent_phi`.
  Apply that readout before objective scoring and final reporting.
- Promote pairwise epistasis-amplitude summaries into the objective:
  `magnitude_epistasis_mean`, `magnitude_epistasis_median`,
  `magnitude_epistasis_sd`, `signed_pairwise_epistasis_mean`,
  sign/reciprocal sign epistasis, diminishing returns, and accessibility.
- Keep the RES-040 local-energy failure terms in the exact descriptive gate:
  edge/site Dirichlet-energy concentration plus the epistasis-amplitude terms.
- Report raw latent epistasis as `raw_*` columns so the analysis can distinguish
  an observable-readout fix from a true latent-geometry fix.
- Re-evaluate the selected RUN-098 full-ABC particles through the same calibrated
  reporting path before comparing held-out distances.
- Interpretation boundary: success means the observable calibrated landscape is
  usable for PHOT-style accessibility/neutral-network analysis. It does not prove
  the raw latent nested-min state has empirically correct pairwise curvature.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
- `scripts/phot_nested_fit_protocol_benchmark.py`
- `scripts/phot_nested_softmin_ensemble_abc.py`
- `scripts/phot_coupled_family_discriminating_abc.py`
- `src/project_adapt_env/smc_abc.py`

## Configuration

- `config.yaml` defines the positive affine readout calibration, epistasis-amplitude
  objective terms, exact descriptive gate, fitted families, cheap-fitting budgets,
  equivalence margins, and quick-mode smoke-test overrides.

## Execution

- Run through `labproj submit EXP-043 --target lab-slurm --execute`.
- `run.sh` pulls the required RUN-077 and RUN-098 DVC targets if missing, then
  runs `scripts/phot_nested_fit_protocol_benchmark.py` with the EXP-043 config.

## Expected outputs

- summary.json
- preflight.json
- fit_protocol_candidates.csv
- selected_fit_protocol_candidates.csv
- fit_protocol_distances.csv
- fit_protocol_equivalence.csv
- graphlet_metrics.csv
- graphlet_summary.csv
- local_abc_rounds.csv
- local_abc_parameter_specs.csv
- walsh_terms.csv
- gft_coefficients.csv
- module_energy_summary.csv
- abc_comparator_selected_family_distances.csv
- full_abc_comparator_selected_particles.csv

## Analysis plan

- Create an ANA record after execution comparing calibrated observable epistasis
  amplitude, raw latent epistasis amplitude, accessibility, local-energy
  concentration, graphlet distance, and recomputed full-ABC comparator distances.
- Treat a reduction only in calibrated metrics as a readout-scale fix, not as
  evidence that the latent nested-min state itself has empirically calibrated
  pairwise curvature.
- Do not move large WHT/GFT tables into the analysis unless needed; compact
  summaries and remote runtime paths are sufficient for the first result record.

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

- Local quick smoke completed before submission with `LABPROJ_QUICK=1`, writing
  disposable outputs under `scratch/exp043_quick/outputs`. The smoke confirmed
  that calibrated observable epistasis metrics and raw latent `raw_*` epistasis
  metrics are both present in `graphlet_summary.csv`.
- Focused regression coverage for the benchmark runner passed with
  `.venv/bin/pytest tests/test_phot_nested_fit_protocol_benchmark.py -q`.
- `RUN-111` is the dry-run record created by
  `labproj submit EXP-043 --target lab-slurm --dry-run`.
- `RUN-112` is the live submission created by
  `labproj submit EXP-043 --target lab-slurm --execute`. It was submitted to
  `lab-slurm` as Slurm job `1907` with `8` CPUs, `32G` memory, and a `12:00:00`
  walltime. The runtime output directory is
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-112/outputs`.
- `RUN-112` completed with inferred Slurm exit code `0`. `ANA-042` copied only
  compact source tables into the scientific record and left the large WHT/GFT
  coefficient tables on the runtime volume.
- `RES-042` records the scientific result: the calibrated cheap-fitting
  protocol failed component equivalence (`0/40` operationally equivalent rows),
  did not rescue observable epistasis amplitude or accessibility, and retained
  large raw latent epistasis amplitudes. The result is a negative/mixed
  calibration result rather than support for using this protocol as the default
  GB1/TrpB scale-up route.
