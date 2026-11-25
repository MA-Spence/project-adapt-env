# ANA-074: RUN-190 Gamma concentration comparison

Start with [the figure-led report](analysis.md) or the
[executed notebook](notebooks/analysis.ipynb). The notebook audits the frozen data,
recomputes statistics, and writes all seven PNG/PDF figures and supporting tables.
The paired [Python source](notebooks/analysis.py) is generated from the notebook
by `labproj analysis sync ANA-074 --write`.

## Reproduce from the project root

Python 3.11.15 was used. Create an environment and install the pinned requirements:

```bash
python3.11 -m venv /tmp/ana074-reproduction
/tmp/ana074-reproduction/bin/python -m pip install -r analyses/2026-11-25_ANA-074_run-190-gamma-concentration-versus-gaussian-trpb-calibration/notebooks/requirements.txt
labproj data restore data/processed/aggregate-thermodynamic-gnk-calibration/RUN-190 --pull
MPLCONFIGDIR=/tmp/ana074-mpl OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /tmp/ana074-reproduction/bin/python analyses/2026-11-25_ANA-074_run-190-gamma-concentration-versus-gaussian-trpb-calibration/notebooks/execute.py
```

The project's configured DVC remote must be mounted/reachable for restoration.
The full original analysis environment is recorded in
[environment.lock.txt](notebooks/environment.lock.txt); package versions used by
the notebook are also emitted to [analysis_environment.json](tables/analysis_environment.json).
The execution helper creates a temporary kernelspec for its own interpreter and
runs the notebook from the repository root, then saves its executed output.

Inputs are entirely within the RUN-190 bundle: outputs, logs, frozen runtime and
the nested RUN-188/RUN-186 controls. The SHA-256 pins in the notebook identify
these exact inputs. No editable Adapt-Env source, network query, fit, new seed or
new noise draw is used. Full pairwise fields are reconstructed algebraically from
saved data. All outputs are regenerated; no hidden notebook state is required.

Figures use eight originally fresh seeds as primary; tables retain all 24.
Resolution sensitivity and individual paired differences are explicit. “Fresh”
refers to the original design: these panels/seeds have since been inspected.
No inferred biological interaction identities, learner performance or optimization
policy outcomes are claimed.
