# ANA-069 — RUN-175 interaction-strength sweep

[Read the figure-led analysis](analysis.md) or open the [executed notebook](notebooks/analysis.ipynb).

- Parent: EXP-065 (non-production); RUN-175 attempt 2, Slurm 2322, completed.
- Five figures in PNG and vector PDF; eight executed notebook cells.
- All ratios and both seed cohorts retained, with observation and resolution sensitivities.
- Numerical tables and independent verification results are in `tables/`.

## Reproduce

From the project root, restore `data/processed/aggregate-thermodynamic-gnk-calibration/RUN-175.dvc` if needed, install the versions in `notebooks/requirements.txt`, then run:

```bash
python analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/notebooks/execute.py
```

Alternatively, execute the notebook top to bottom. The analysis imports source frozen in the RUN bundle; endpoint replay also reads ANA-068's graph-summary CSV. `analysis.py` is the paired text representation. No cluster access or rerun of model fitting is required to regenerate the figures. Bands describe seed variation, not confidence intervals or acceptance gates.
