# ANA-064 — RUN-166 paired readouts and embedded neighbourhoods

Read [analysis.md](analysis.md) for the scientific interpretation and five figures.
The canonical [analysis.ipynb](notebooks/analysis.ipynb) has been executed from
start to finish: **8 code cells, no errors**, with all five figures displayed
inline. Every figure is also exported as a PNG and vector PDF.

## Reproduce

From the project root, create or select a Python 3.13 environment and run:

```bash
python -m pip install -r analyses/2026-10-07_ANA-064_run-166-paired-readouts-and-embedded-neighbourhoods/notebooks/requirements.txt
python analyses/2026-10-07_ANA-064_run-166-paired-readouts-and-embedded-neighbourhoods/notebooks/execute_notebook.py
```

Alternatively, select that environment as the notebook kernel and run all cells.
The execution helper uses nbclient and retains cell outputs. The Jupyter kernel
needs permission to open local sockets. Temporary Jupyter/Matplotlib files live
outside the repository. The notebook does not submit a job or regenerate model
fields: it audits and analyses the frozen RUN data.

Two data dependencies must be restored. With the project's configured DutchZFS
share mounted, use:

```bash
dvc pull data/processed/aggregate-thermodynamic-gnk-calibration/RUN-166.dvc
dvc pull data/processed/ensemble-landscape-validation-datasets.dvc
```

The latter supplies
`data/processed/ensemble-landscape-validation-datasets/datasets/trpb_johnston2024_four_site_graphfla/variants.parquet`.
The notebook deliberately requires this original prepared table to verify the
saved empirical measurements and input hash. It does not silently skip that
audit when the table is unavailable. All analysis dependencies, including
PyArrow for Parquet input, are pinned in `notebooks/requirements.txt`.

## Data and verification

The RUN bundle contains the 17 original outputs, three logs, remote SHA-256
manifest, completed RUN metadata, and the exact submitted source snapshot.
Its DVC pointer is
`data/processed/aggregate-thermodynamic-gnk-calibration/RUN-166.dvc`, directory
MD5 `c128be3fda98517848f7d5533b5d2f52.dir`: **83 files, 54,423,275 bytes, 82 unique
cache objects including the directory manifest**. Local files and cache objects
were independently hash-verified. All 82 cache objects were also backed up and
hash-verified on the existing project cluster store; the verification is recorded
in `tables/data_bundle_audit.json`.

Collection used SSH because `labproj collect --dry-run` found the client mount
unavailable. The original outputs remain at
`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-166`, and submitted source
at `/tank/data/labproj/staging/project-adapt-env/RUN-166/project` on the existing
`raspi-1` cluster gateway. The project-local DVC environment timed out while
loading cloud-backed dependencies. Tracking therefore used temporary DVC 3.63.0
with pathspec 0.12.1; backup uses the existing store at
`/tank/data/labproj/projects/project-adapt-env/data/files/md5` over SSH. No
permanent Python environment or remote configuration was changed.

The independent audit reconstructs all 3,312 observed cubes from saved additive
and GNK values, and recalculates the ten distance descriptors plus WHT quantities
for 9,936 resolution-specific rows. It verifies 105,984 mutation edges, all 40
candidate calibrations/losses, all 24 held-out distances, the 768 empirical node
measurements, all declared source/input hashes, and genotype/seed separation.
This is numerical replay of saved latents, not bitwise model reconstruction.

The ten-feature distance is an energy V-statistic with empirical calibration
standardization and explicit missingness coordinates. Descriptive uncertainty
uses 2,000 resamples, seed 664001, preserving synthetic seed clusters and pairing.
Calibration is not repeated. The report explicitly discusses dependence,
measurement error, missingness, resolution, and the limited three-seed embedded
panel. No scientific acceptance thresholds were added.

## Artifacts

- `notebooks/analysis.ipynb`: canonical executed analysis, methods and figures.
- `notebooks/analysis.py`: labproj-generated percent-format text companion.
- `notebooks/analysis_tools.py`: independent numerical audits and figure code.
- `figures/01_…` through `05_…`: paired readouts, conditional statistics, aggregate
  fingerprint, embedded neighbourhoods, and saturation/availability (PNG/PDF).
- `tables/distances.csv` and `contrasts.csv`: all resolutions and latent ratios,
  with paired conditional bootstrap intervals.
- `tables/descriptors.csv`: primary-resolution means, medians, available counts
  and conditional intervals for every readout/ratio arm.
- `tables/operating_summary.csv`: fixed-parameter p₀ interventions.
- `tables/embedded_means.csv`, `embedded_seed_means.csv`, and
  `embedded_availability.csv`: both reference ranks and both site policies.
- `tables/audit.json`, `metrics.json`, and `data_bundle_audit.json`: verification,
  environment versions, source hash, findings and data provenance.
- The original per-cube/per-edge/per-latent data and calibration-defined
  conditional bins remain in the DVC bundle, without overwriting RUN outputs.

After edits, rerun the canonical notebook, then use
`labproj analysis sync ANA-064 --write` to update its text companion and artifact
index. Repository-wide status/refresh/check commands timed out on historical
cloud-backed files. Targeted analysis audits are recorded separately; no global
strict-check pass is claimed. No new RES or CLAIM accompanies this ANA.

Completion commands saved ANA-064 and EXP-065 as completed before their global
follow-up checks timed out. Targeted metadata and registry checks confirmed the
completed statuses and linked artifacts.
