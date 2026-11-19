# EXP-041: PHOT exact descriptive failure gate for full ABC fit

## Status

- completed

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Before running any new ABC search, which exact dense PHOT graphlet summaries does the existing RUN-039 full ABC fit fail relative to the empirical latent_phi landscape?

## Pre-experiment prediction

The RUN-039 full ABC fit will retain useful low-order predictive structure but will fail exact dense-subspace topology or spectral/local-energy summaries that were not part of the original validation objective.

## Rationale

The current ABC workflow is expensive and can obscure model misspecification by optimizing a broad summary vector. This experiment first computes deterministic complete-graphlet summaries on PHOT and emits a ranked failure ledger; a later ABC objective is justified only for summaries shown to fail descriptively.

## Experimental design

- Load the existing `RUN-039` PHOT prepared assay, selected panel, MAVE-NN
  wildtype metrics, and full validation-objective ABC summary.
- Enumerate complete observed PHOT graphlets deterministically from the most
  frequent mutation features for dimensions 4, 5, and 6.
- Evaluate every node in binary hypercube order, with no random walk sampling
  and no ABC resampling.
- Compare exact empirical `latent_phi` summaries against four fixed `RUN-039`
  branches: SMC-ABC best particle, SMC-ABC posterior mean, baseline shared
  stability-readout, and predictive rich-pair stability-readout.
- Report a ranked component-failure table with scaled signed errors for
  topology, accessibility, sign epistasis, Walsh/GFT spectral allocation, basin
  structure, and local Dirichlet-energy concentration.
- Emit `abc_descriptive_gate.csv`; a later ABC experiment is eligible only for
  summary components that exceed the preregistered exact descriptive failure
  threshold.

## Inputs

- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/summary.json`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/selected_panel.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/checkpoints/panel_preparation/PHOT_CHLRE_Chen_2023/prepared_assay.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/checkpoints/panel_preparation/PHOT_CHLRE_Chen_2023/mavenn_metrics.json`
- `data/interim/proteingym_mmseqs_alignments/PHOT_CHLRE_Chen_2023.fasta`
- `experiments/2026-05-14_EXP-013_proteingym-phot-stability-readout-validation-objective-control/config.yaml`
- `scripts/phot_exact_descriptive_failure_analysis.py`
- `src/project_adapt_env/dense_descriptive.py`

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- summary.json
- graphlet_catalog.csv
- exact_graphlet_values.csv
- exact_graphlet_metrics.csv
- exact_model_failures.csv
- exact_model_distances.csv
- abc_descriptive_gate.csv

## Analysis plan

- Inspect `exact_model_distances.csv` to rank fixed full-ABC and baseline
  candidates by exact descriptive failure.
- Inspect `exact_model_failures.csv` to identify the specific summary
  components that fail and the direction of failure.
- Treat `abc_descriptive_gate.csv` as the only admissible source of target
  summaries for any later ABC refit. If no components exceed the gate, the
  current model is not justified for another ABC pass on PHOT.

## Completion criteria

- [x] Config recorded
- [x] Run script recorded
- [x] Environment recorded
- [x] Outputs generated
- [x] Important compact outputs recorded without collecting large exact tables
- [x] Analysis record created
- [x] Result record created
- [x] Hypothesis updated
- [x] PROJECT_STATE.md updated

## Post-experiment notes

- `RUN-108` completed on `lab-slurm`.
- `ANA-040` reviewed compact copied outputs from the remote runtime directory
  without running full `labproj collect`; `exact_graphlet_values.csv` and
  `exact_graphlet_metrics.csv` remain on
  `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-108/outputs`.
- `RES-040` records the main finding: the fixed `RUN-039` PHOT branches do not
  broadly fail the exact complete-graphlet package, but all fail four
  local-energy concentration summaries. Later ABC is therefore justified only
  for those summaries unless another exact diagnostic adds targets.
