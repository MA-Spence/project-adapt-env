# EXP-001: ProteinGym-DMS-Distributional-Realism-Panel

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001

## Question

Can Adapt-Env recover local distributional statistics from real ProteinGym DMS assays when alignments are built from assay wild-type sequences using the MMseqs API?

## Pre-experiment prediction

A shared synthetic parameter regime should match empirical DFE and mutational-effect summary statistics across a panel of ProteinGym substitution assays better than naive synthetic controls, but with assay-to-assay residual mismatch.

## Rationale

HYP-001 should be tested first against external DMS assays because benchmark usefulness depends on reproducing real sample statistics before any stronger claim about exact landscape structure is attempted.

## Experimental design

- Template: generic
- Download the ProteinGym substitution reference table and substitutions parquet.
- Select a balanced panel of single-mutant assays filtered by taxon, sequence
  length, and mutant count.
- Use each assay wild-type sequence as the MMseqs query sequence and save the
  returned alignment as FASTA under `data/interim/proteingym_mmseqs_alignments`.
- Standardize each assay's DMS scores within assay, set the wild type to zero,
  and expose the assay through an empirical-landscape adapter.
- Use Adapt-Env calibration machinery to fit one shared synthetic parameter
  regime across the selected assay panel.
- Write the selected panel, fitted parameters, validation metrics, and per-assay
  comparison summaries to the run output directory.

## Inputs

- ProteinGym substitution reference metadata
- ProteinGym substitution parquet benchmark
- wild-type query sequences from the selected DMS assays

## Configuration

- `config.yaml` records the panel-selection filter, MMseqs settings, Adapt-Env
  base configuration, and calibration search grids.
- The current implementation assumes assay-wise z-scored DMS scores and uses a
  synthetic wild-type fitness of zero for calibration.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) a reproducible ProteinGym assay panel and metadata table
- 2) downloaded real DMS assay files
- 3) MMseqs-derived alignment FASTAs saved under data/interim
- 4) empirical versus synthetic summary-statistic comparisons
- and 5) a Slurm-executable experiment scaffold.

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
