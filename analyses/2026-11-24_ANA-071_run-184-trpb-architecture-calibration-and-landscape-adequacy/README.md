# ANA-071 — RUN-184 TrpB landscape adequacy

Start with **[analysis.md](analysis.md)** for five figures and the graduation assessment.
The achieved calibration does not adequately reproduce TrpB's fitness distribution,
local epistasis or sparse-panel navigation. Computational reproducibility passes.

| Notebook | Role |
|---|---|
| [analysis.ipynb](notebooks/analysis.ipynb) | Figures 1–2: native fitness distribution, high tail, mutation effects, sign epistasis, walks, peaks and replicate/resolution sensitivity. |
| [architecture.ipynb](notebooks/architecture.ipynb) | Figures 3–5: all 18 targets, loss decomposition, empirical coefficient prediction, bootstrap instability, search and six parameters. |
| [validation.ipynb](notebooks/validation.ipynb) | Frozen hashes, panels, objectives, strict convergence evidence, native oracle replay and empirical-summary reconstruction. |

All three notebooks are executed, with figures embedded. Percent-format `.py`
sources are paired with their cells. Helpers regenerate the complete [tables](tables/)
and [PNG/PDF figures](figures/) from frozen inputs.

## Reproduce

Use Python **3.11.15** from the repository root. Restore RUN-184 if absent:

```sh
labproj data restore data/processed/aggregate-thermodynamic-gnk-calibration/RUN-184 --pull
python3.11 -m venv .labproj/envs/ana071-py311
.labproj/envs/ana071-py311/bin/python -m pip install -r analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/notebooks/environment.lock.txt
MPLBACKEND=Agg MPLCONFIGDIR="$PWD/.labproj/tmp/ana071-mpl" IPYTHONDIR="$PWD/.labproj/tmp/ana071-ipython" .labproj/envs/ana071-py311/bin/python analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/notebooks/execute.py
.labproj/envs/ana071-py311/bin/python analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/notebooks/verify_artifacts.py
```

The executor creates a temporary kernel specification for that interpreter and
executes validation, landscape and architecture notebooks in fresh kernels.
Local Jupyter loopback sockets must be allowed. Original execution used the
existing `ana070-py311` environment with the recorded package versions; that
name does not select different scientific inputs.

After editing paired `.py` sources, regenerate before executing:

```sh
.labproj/envs/ana071-py311/bin/python analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/notebooks/build_notebooks.py
```

No internet download, Optuna optimization or auxiliary refitting is performed.
Validation imports the model only from RUN-184's hashed source snapshot and
replays every confirmation prediction. `landscape_tools.py` independently
constructs observed edges/faces and exact uphill expectations;
`architecture_tools.py` reconstructs weights, losses and held-out empirical
coefficient predictions using NumPy/SciPy.

## Interpretation and provenance

Confirmation genotypes are disjoint for this run, but come from the already-studied
TrpB assay. Four seeds are generator realizations, not biological replicates or
confidence intervals. Synthetic observations are noiseless. Graph results concern
the observed panel only. Architecture summaries are fitting targets; local graph
metrics, upper-tail percentiles and auxiliary held-out prediction are additional
assessments. No retrospective numerical adequacy gate is imposed.

The [RUN record](../../experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-184.yaml)
and [protocol](../../experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/architecture_calibration.md)
retain parameters, masks, seeds, failures and execution evidence. Previous RUNs
remain historical context, not matched controls for this changed protocol.

`tables/validation_checks.json`, `tables/landscape_verification.json`,
`tables/architecture_findings.json` and `tables/artifact_validation.json` retain
computational checks. `artifact_manifest.json` identifies analysis outputs.
The analysis copy is DVC-tracked under
`data/processed/aggregate-thermodynamic-gnk-calibration/ANA-071`.
Post-archive backup verification is in `storage_verification.json`, outside that
copy to avoid circular identities. Defaults, presets and existing claims are unchanged.
