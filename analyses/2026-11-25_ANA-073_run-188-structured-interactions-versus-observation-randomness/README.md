# Reproduce ANA-073

Read [analysis.md](analysis.md) for the figure-led interpretation. The canonical
[executed notebook](notebooks/analysis.ipynb) contains the audits, all calculations,
six figures and interpretation. [analysis.py](notebooks/analysis.py) is its paired
percent-cell source. No fitting, new observation draws, or ML/search task is run.

From the project root, restore the immutable input archive if absent:

```bash
labproj data restore data/processed/aggregate-thermodynamic-gnk-calibration/RUN-188 --pull
```

Use Python 3.11.15 and the exact analysis environment in
[environment.lock.txt](notebooks/environment.lock.txt). The narrower direct
requirements are also retained in [requirements.txt](notebooks/requirements.txt).
The validated local interpreter is `.labproj/envs/ana070-py311/bin/python`.
To reproduce in a separate environment, create a Python 3.11.15 virtual environment
and install the environment lock, then use that environment's Python below.

```bash
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .labproj/envs/ana070-py311/bin/python \
  analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/notebooks/execute.py
```

`execute.py` creates a temporary kernelspec pointing to its own interpreter and
executes the notebook from a fresh kernel at the project root. Jupyter requires
local loopback sockets. The notebook may also be opened and run top-to-bottom
with this environment. It finds the project root from the working directory.
It requires the complete RUN-188 archive, including its frozen runtime and embedded
RUN-186 source/data; mutable package code is not used.

The notebook overwrites its own derived figures/tables. It authenticates the
input manifest against an embedded SHA-256 pin and verifies every listed payload.
The only regenerated random numbers are the original uniform arrays from recorded
noise seeds, used to replay the existing observations exactly. Independent-only
and shared-only layers subtract known saved components; they are not refits.

Figures are saved as PNG and vector PDF. `tables/seed_layer_stratum_summary.csv`
contains all 24 seeds and three resolutions, including empty strata and signed
powers. The main figures use eight originally fresh seeds. Activity membership for
intermediate layers is fixed to the corresponding arm's fully observed labels;
it is not comparable to the original RUN's intrinsic activity stratification.
Original intrinsic powers and observed face classifications are verified exactly.

All code cells must have execution counts and no error outputs. The six embedded
figure images must agree with the exported files; see
[tables/artifact_validation.json](tables/artifact_validation.json).

Durable archives:

- RUN inputs/output/source: `data/processed/aggregate-thermodynamic-gnk-calibration/RUN-188`
- Completed analysis: `data/processed/aggregate-thermodynamic-gnk-calibration/ANA-073`
- Backup verification: [storage_verification.json](storage_verification.json)
- RUN backup verification: [run_storage_verification.json](run_storage_verification.json)

The analysis archive excludes its own `storage_verification.json` sidecar to avoid
a circular archive checksum. It includes the report, executed notebook, paired
source, figures, tables, metadata and environment locks. The existing project DVC
remote is `/Volumes/DutchZFS/labproj/projects/project-adapt-env/data`.
