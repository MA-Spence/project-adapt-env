# EXP-003: ProteinGym Raw-Scale Latent Observation Calibration Panel

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001

## Question

Does HYP-001 hold more cleanly when assay-coordinate confounds are removed by fitting homogeneous multi-mutant ProteinGym stability assays on raw scores, using assay-specific latent observation models and non-affine measurement transforms instead of assay-wise z-scoring and wild-type zeroing?

## Pre-experiment prediction

Within a homogeneous multi-mutant stability panel, a shared Adapt-Env regime fit on raw scores with assay-specific nonlinear observation models and on MAVE-NN-derived latent phenotypes should outperform the original z-scored baseline from EXP-001, while per-assay latent fits should still exceed the shared fit and thereby reveal the remaining gap attributable to true cross-protein heterogeneity.

## Rationale

RES-001 weakened HYP-001, but that run mixed assay families, z-scored each assay, fixed wild type to zero, and excluded double-mutant signal from calibration. RES-002 then showed that the raw Adapt-Env family is expressive enough to span broad local DFE regimes. This experiment is therefore run to isolate a more specific question: whether the main failure in RES-001 came from coordinate-system and observation-model mismatch rather than from lack of landscape expressivity. It does so by using a homogeneous ProteinGym stability subpanel with real multi-mutant measurements, learning assay-specific latent observation maps from raw scores, and comparing shared versus per-assay Adapt-Env calibration on those cleaner coordinates.

## Experimental design

- Template: generic
- Assemble a homogeneous ProteinGym panel restricted to `cDNA display proteolysis`
  stability assays with real multi-mutant measurements and sequence lengths short
  enough for tractable latent-phenotype fitting.
- Use each assay wild-type sequence as an MMseqs query and save the returned real
  alignment FASTA under `data/interim/proteingym_mmseqs_alignments`.
- Fit an assay-specific `mavenn` additive global-epistasis model on raw assay
  scores to infer a latent phenotype scale and an assay-specific nonlinear
  measurement process.
- Use the fitted `mavenn` model to estimate the wild-type location on the raw
  assay scale and on the latent phenotype scale instead of z-scoring the assay
  and forcing wild type to zero.
- Build two empirical branches for Adapt-Env calibration from the same observed
  single- and double-mutant sequences:
  1. a raw-score branch with non-affine observation models enabled in Adapt-Env
  2. a latent-phenotype branch using `mavenn` inferred `phi` values
- Fit one shared Adapt-Env regime across the panel on each branch, then fit
  per-assay latent-branch calibrations as a ceiling comparison.
- Record fitted parameters, validation metrics, and per-assay latent-model
  diagnostics to the run output directory.

## Inputs

- ProteinGym substitution reference metadata
- ProteinGym substitutions parquet benchmark
- wild-type query sequences from the selected stability assays
- `mavenn` assay-specific latent phenotype fits derived from raw ProteinGym
  observations

## Configuration

- `config.yaml` records the homogeneous panel-selection filter, MMseqs settings,
  `mavenn` fit controls, Adapt-Env base configuration, and separate calibration
  branches for raw-score shared fit, latent shared fit, and per-assay latent
  fits.
- The current implementation resolves the coordinate issue by replacing
  assay-wise z-scoring and synthetic wild-type zeroing with assay-specific
  `mavenn` estimates of the wild-type location and latent phenotype scale.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) a reproducible homogeneous ProteinGym stability panel and metadata table
- 2) MMseqs-derived real alignment FASTAs saved under `data/interim`
- 3) per-assay `mavenn` latent-phenotype diagnostics from raw scores
- 4) a shared raw-score Adapt-Env calibration summary with assay-specific
  nonlinear observation models
- 5) a shared latent-phenotype Adapt-Env calibration summary
- 6) per-assay latent-branch calibration summaries
- and 7) a Slurm-executable experiment scaffold

## Analysis plan

- Create an ANA record after execution.

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
