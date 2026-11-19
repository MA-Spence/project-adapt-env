# ANA-067 — RUN-172 production landscape analysis

Open [`notebooks/analysis.ipynb`](notebooks/analysis.ipynb) for the executed,
figure-led analysis, or [`analysis.md`](analysis.md) for the interpretation.
All figures export to PNG and PDF in `figures/`; exact estimates, bootstrap
intervals, corrected suitability inventory and audit tables are in `tables/`.

```bash
python -m pip install -r notebooks/requirements.txt
python notebooks/execute.py
```

Run from this directory, or provide the full path to `execute.py`. The executor
uses the invoking Python interpreter through a temporary Jupyter kernel spec.
The notebook locates the repository from its working directory. It reads the
DVC-restored `data/processed/production_landscapes/RUN-172` bundle and its frozen
source snapshot, independently audits it, then regenerates every table/figure.
It neither refits models nor submits jobs. The paired `analysis.py` and
`analysis_tools.py` make all notebook operations inspectable.

Restore data if necessary with:

```bash
labproj data restore data/processed/production_landscapes/RUN-172 --pull
```

The verified DVC object store is on the cluster at
`/tank/data/labproj/projects/project-adapt-env/data/files/md5`. The configured
client share was unavailable during this analysis, so collection/backup used
SSH; `tables/data_bundle_audit.json` records object verification. The exact
published RUN remains unchanged. Three omitted inventory exclusion rows are
restored only in `tables/suitability_complete.csv`.

Figures:

1. Coverage, calibration and local spread.
2. Withheld global statistics across all evaluated assays.
3. Every-assay residual atlas.
4. Named fitness-distribution and local-spread examples.
5. TrpB observation sensitivities, including intrinsic/observed separation.
6. Resolution and diminishing-returns denominator sensitivity.
