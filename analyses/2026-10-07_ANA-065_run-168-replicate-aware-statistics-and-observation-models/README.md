# ANA-065: RUN-168 replicate-aware statistics and observation models

Read [analysis.md](analysis.md) for the figure-led scientific account. The canonical
executable analysis is [notebooks/analysis.ipynb](notebooks/analysis.ipynb). It audits
the raw inputs and frozen simulations, regenerates all analysis tables and five
PNG/PDF figure pairs, and retains inline figures in its executed output.

## Reproduce

From the repository root, restore both versioned data bundles:

```bash
dvc pull data/processed/aggregate-thermodynamic-gnk-calibration/RUN-168.dvc
dvc pull data/processed/trpb-replicate-observation-exp065.dvc
```

The restore commands require the configured labcluster mount to be available;
the SSH storage location is recorded below. The first is the completed
output/log/submitted-source bundle. The second retains
the original TrpB replicate and control data, selected public processing code,
exclusions, and frozen RUN-166 inputs. No additional RUN-166 download is required.
Neither external notebook is executed. The analysis refuses to skip its input
checks when a bundle is missing or has changed.

Use Python 3.13 (the analysis execution version) and install the pinned packages:

```bash
python3 -m venv /tmp/adaptenv-ana065
/tmp/adaptenv-ana065/bin/python -m pip install -r analyses/2026-10-07_ANA-065_run-168-replicate-aware-statistics-and-observation-models/notebooks/requirements.txt
OPENBLAS_NUM_THREADS=1 /tmp/adaptenv-ana065/bin/python analyses/2026-10-07_ANA-065_run-168-replicate-aware-statistics-and-observation-models/notebooks/execute.py
```

Alternatively select a kernel with those dependencies and run every notebook cell
in order. `execute.py` creates a temporary kernel specification tied to its own
Python interpreter; it does not install or modify a user-wide kernel. It executes
from the repository root, overwriting the canonical notebook only on success.
Allow approximately 2 GB of RAM and a few minutes. Relative repository discovery
supports execution from the root or an analysis subdirectory.

After edits and successful execution:

```bash
labproj analysis sync ANA-065 --write
```

This updates the percent-format `notebooks/analysis.py` text companion and indexes
figures/tables. Do not edit the text companion as a separate canonical analysis.

## Scientific scope

This is descriptive reanalysis of the same 48 empirical cubes and eight frozen
RUN-166 arms, not an independent empirical confirmation or new model fit. Eight
synthetic seed clusters provide 64 saved confirmation cubes per arm. The 128
observation draws reuse those cubes; they are not 128 independent landscapes.
Empirical cube-bootstrap intervals condition on noise estimates and omit shared
site/allele dependence. The report separates equal-cube mean WHT fractions from
ratios of pooled powers, keeps signed corrections and unavailable ratios, and
retains empty activity strata.

Observation comparisons assume centered exchangeable errors, independent errors
between genotypes and transfer from measured donor means to synthetic latent
fitness. The shared-stop branch additionally assigns stop cross covariance to
common error; control biology and shared normalization are not identified away.
Training/holdout control agreement tests repeatability of those control summaries,
not their causal interpretation or transfer to functional variants. No scientific
success gates, specific empirical interaction inference, RES or CLAIM are added.

## Data provenance

The RUN bundle is `data/processed/aggregate-thermodynamic-gnk-calibration/RUN-168`:
**165 files, 224,508,438 bytes**, DVC directory MD5
`b71af13a0354c4a7a8a65cb77e8f019a.dir`.
Its **159 unique cache objects**, including the directory manifest, are backed up
to the existing project lab store and individually MD5-verified. The input bundle
is `trpb-replicate-observation-exp065.dvc`, directory MD5
`ad4e31a4a5704520d0a4be2cab9fa281.dir` (15 files, 90,025,877 bytes).
See `tables/data_bundle_audit.json` for the collection/backup audit.

`labproj run reconcile RUN-168` confirmed completion. The collection preview found
no client-visible NFS/SMB share, so collection used SSH from the existing gateway
`matthew-spence@raspi-1`. Runtime files remain under
`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-168` and submitted source
under `/tank/data/labproj/staging/project-adapt-env/RUN-168/project`.
The existing DVC object store is
`/tank/data/labproj/projects/project-adapt-env/data/files/md5`.

The project-local DVC wrapper timed out while accessing historical cloud-backed
dependencies. Tracking used temporary DVC 3.63.0 with pathspec 0.12.1, preserving
the repository's environment and remote configuration. Broad labproj
status/refresh/check commands also timed out; no repository-wide strict-check pass
is claimed. Targeted notebook, provenance and record checks are documented.

## Audit and artifacts

`notebooks/analysis_tools.py` implements the independent butterfly WHT,
covariance calculations, donor reconstruction, observation replay and figures.
It never imports the submitted runner. Audited totals include 162 remote file
hashes, 73 output hashes, eight submission hashes, 14 input/source hashes, 768
empirical measurements, 193,170 activity flags, 131,584 synthetic replicate-pair
cubes, 394,752 synthetic statistic rows and 4,210,688 saved observation values.
The original empirical bootstrap is reproduced, and the first/last observation
draw and intrinsic aggregate summaries are independently recomputed. All
per-cube synthetic statistics are checked, not just these summary samples.

- `figures/01_…` through `05_…`: measurement/control diagnostics, WHT weighting,
  slopes/availability, frozen observation comparisons, and activity strata.
- `tables/empirical_summary.csv`: both empirical splits, all three resolutions,
  estimators and activity strata, including zero counts.
- `tables/empirical_bootstrap.csv`: original intervals plus mean-per-cube WHT and
  order-power intervals, preserving unavailable draws.
- `tables/frozen_model_summary.csv`: conditional observation means/ranges for all
  arms, resolutions, strata and estimators.
- `tables/observation_contrasts.csv`: paired changes from each frozen intrinsic arm.
- `tables/seed_summaries.csv`: original seed-level summaries. In noisy modes counts
  include 128 repeated observations of each of eight cubes per seed; these are
  repeated measurements, not additional landscapes.
- `tables/control_diagnostics.csv`, `donor_holdout_diagnostics.csv`: held-out
  checks, with no subsequent noise tuning.
- `tables/audit.json`, `metrics.json`, `data_bundle_audit.json`: numerical checks,
  scoped headline findings, versions and data verification.

Original per-vertex replicates, per-cube statistics, mutation components and masks
remain unchanged in the DVC RUN bundle. The prior ANA-063 and ANA-064 records are
not retrospectively edited.

Completion commands saved ANA-065 and EXP-065 as completed before follow-up
checks timed out. Targeted metadata and registry reads confirm completed status
for RUN-168, ANA-065 and EXP-065.
