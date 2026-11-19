# ANA-063 — RUN-164 aggregate thermodynamic GNK analysis

Read [analysis.md](analysis.md) for the scientific interpretation and five
figures. Open [notebooks/analysis.ipynb](notebooks/analysis.ipynb) for executable
methods, retained cell outputs and all figures. The notebook was executed from
start to finish: **7 code cells, no errors**, using Python 3.13.11.

## Reproduce

From the project root, use a Python environment with the pinned dependencies:

```bash
python -m pip install -r analyses/2026-10-07_ANA-063_run-164-aggregate-thermodynamic-gnk-confirmation/notebooks/requirements.txt
python analyses/2026-10-07_ANA-063_run-164-aggregate-thermodynamic-gnk-confirmation/notebooks/execute_notebook.py
```

Alternatively, select that environment as the notebook's Python kernel and run
all cells. The execution helper uses `nbclient` and retains outputs in the ipynb;
its local kernel needs permission to open localhost sockets. Temporary Jupyter
and Matplotlib state lives outside the repository. The analysis requires no
queue access or full experiment rerun.

The notebook reads
`data/processed/aggregate-thermodynamic-gnk-calibration/RUN-164/`. This bundle
contains the eight original run outputs, three logs, remote SHA-256 manifest,
completed RUN metadata, and the submitted source snapshot. The source snapshot
is used for nine finite-length scaling replays; the current working-tree
generator is not imported. Original source and output hashes are checked before
analysis. Original full assay tables are checked against their hashes when
available; the stored sampled measurements suffice to reproduce the analysis.

The bundle was collected over SSH from `raspi-1` because the shared filesystem
was not mounted locally. Original run files remain at
`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-164`, with submitted source
at `/tank/data/labproj/staging/project-adapt-env/RUN-164/project`.

The 60-file bundle is DVC-tracked by
`data/processed/aggregate-thermodynamic-gnk-calibration/RUN-164.dvc`
(directory MD5 `42453e6e1a305fdf4b34fb9dd2d40092.dir`). Its 60 unique cache objects,
including the directory manifest, were backed up to the existing project store
on `raspi-1` and hash-verified there. With the configured `DutchZFS` share mounted,
restore it using:

```bash
dvc pull data/processed/aggregate-thermodynamic-gnk-calibration/RUN-164.dvc
```

The project-local DVC executable failed to load cloud-backed dependencies during
collection. Tracking used a temporary DVC 3.63.0 environment with pathspec 0.12.1;
backup copied only this bundle's verified content-addressed objects over SSH to
`/tank/data/labproj/projects/project-adapt-env/data/files/md5`. No existing project
Python environment or remote configuration was changed.

## Artifacts and provenance

- `notebooks/analysis.ipynb`: canonical, executed notebook.
- `notebooks/analysis.py`: labproj-generated text companion.
- `notebooks/analysis_tools.py`: independent audits, bootstrap and plotting code.
- `figures/01_…` through `05_…`: PNG previews and vector PDF exports.
- `tables/aggregate_summary.csv`: held-out distances and conditional intervals.
- `tables/model_contrasts.csv`: paired composite-minus-baseline comparisons.
- `tables/descriptor_summary.csv`: means, medians, discrepancies and availability.
- `tables/scaling_replay.csv`: original timings, stable RMS and numerical caveats.
- `tables/audit.json` and `tables/metrics.json`: audit results, source hash, package
  versions, quantitative findings and limitations.

Bootstrap seed is 63001 with 2,000 resamples. Parameters remain fixed at the
RUN-164 selections. Resampling preserves synthetic seed clusters; empirical
cubes are genotype-disjoint but can share sites/alleles. Intervals are conditional
descriptions and exclude fitting uncertainty. Cross-platform scaling replay
reproduces checked statistics numerically, while coefficient fingerprints differ;
the report and output table explicitly retain that limitation.

After editing the canonical notebook, regenerate its text companion and artifact
index with `labproj analysis sync ANA-063 --write`. Repository-wide
status/refresh/check commands timed out on historical cloud-backed files during
this task. The completion command saved completed status before its global
check stalled; targeted status and artifact checks confirm ANA-063 and EXP-065
are completed. No full strict-check pass is claimed.
