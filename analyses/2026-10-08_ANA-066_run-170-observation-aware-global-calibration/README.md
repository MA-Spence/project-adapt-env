# ANA-066: RUN-170 observation-aware global calibration

Start with [the figure-led report](analysis.md) or
[the executed notebook](notebooks/analysis.ipynb). All five figures are generated
by the notebook, with PNG/PDF exports under `figures/`. The full uncertainty table,
audit and portable parameter settings are under `tables/`.

From the repository root, restore the two DVC inputs if required:

```bash
dvc pull data/processed/aggregate-thermodynamic-gnk-calibration/RUN-170.dvc
dvc pull data/processed/trpb-replicate-observation-exp065.dvc
```

The configured client share must be mounted for `dvc pull`. The existing cluster
store is `/tank/data/labproj/projects/project-adapt-env/data/files/md5`.
This run's files were collected over SSH because the client mount was unavailable.
The run bundle contains exact submitted source, outputs, logs, input hashes,
remote SHA256 checksums and reconciled RUN metadata.

Use Python 3.13 with `notebooks/requirements.txt` (plus SciPy 1.17.1) or a compatible
scientific environment. Create a dedicated venv; do not modify the project base
environment solely for this analysis:

```bash
python -m pip install -r analyses/2026-10-08_ANA-066_run-170-observation-aware-global-calibration/notebooks/requirements.txt
python analyses/2026-10-08_ANA-066_run-170-observation-aware-global-calibration/notebooks/execute.py
labproj analysis sync ANA-066 --write
labproj analysis sync ANA-066 --check
```

`execute.py` registers a temporary kernel using its own Python interpreter and
writes the executed notebook only after successful execution. `analysis_tools.py`
uses saved inputs and arrays; it does not invoke the submitted fitting runner or
alter selected parameters. `analysis.py` is the labproj-generated text companion.

Production used Python 3.11.15 / NumPy 2.4.4 on the CPU worker; analysis used Python
3.13 / NumPy 2.4.4. Numerical identities are audited across the two environments.
The 32 observation draws reuse 128 confirmation cubes from 16 landscape seeds.
The uncertainty table resamples seeds, not observation draws as independent units.
All scientific comparisons are descriptive; there are no acceptance gates.
