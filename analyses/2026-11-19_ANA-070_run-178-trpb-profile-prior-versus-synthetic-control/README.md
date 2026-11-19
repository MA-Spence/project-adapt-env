# ANA-070 — RUN-178 TrpB profile prior versus synthetic control

Start with **[analysis.md](analysis.md)** for the figure-led scientific interpretation.

- [analysis.ipynb](notebooks/analysis.ipynb): executed notebook rebuilding all five figures and comparison tables.
- [validation.ipynb](notebooks/validation.ipynb): executed provenance, calibration, interaction-field and independent WHT checks.
- [Figures](figures/): five PNGs with matching vector PDFs.
- [Tables](tables/): every paired seed result, aggregate comparisons, fit parameters and audit outputs.
- [Metadata](metadata.yaml): EXP-065, RUN-178 and archived RUN-175 control links.

The result is mixed: some median agreements improve, several local epistasis measures worsen, and individual-seed improvements are inconsistent. The analysis does not rerun either calibration or simulation.

## Reproduce

Run from the project root with Python **3.11.15**. Restore the frozen RUN-175 and RUN-178 DVC bundles if absent. Both are under `data/processed/aggregate-thermodynamic-gnk-calibration/`; their pointers preserve the original inputs, outputs and numerical source. Remote retrieval requires the configured backup host to be available.

```sh
python3.11 -m venv .labproj/envs/ana070-py311
.labproj/envs/ana070-py311/bin/python -m pip install -r analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/notebooks/environment.lock.txt
MPLBACKEND=Agg MPLCONFIGDIR="$PWD/.labproj/tmp/ana070-mpl" IPYTHONDIR="$PWD/.labproj/tmp/ana070-ipython" .labproj/envs/ana070-py311/bin/python analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/notebooks/execute.py
```

The executor creates a temporary kernel specification using that exact interpreter, then executes validation and analysis in fresh kernels. Standard local Jupyter loopback sockets must be permitted. It saves outputs only after a notebook executes successfully. Alternatively, select the same environment in Jupyter and run each notebook from the repository or notebook directory.

The paired `.py` files are readable percent-format notebook sources. After editing a source, regenerate the notebooks before executing:

```sh
.labproj/envs/ana070-py311/bin/python analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/notebooks/build_notebooks.py
```

`data_tools.py` constructs tables directly from the frozen output CSVs. `plot_tools.py` draws the figures. `audit_tools.py` verifies the archived sources and inputs, independently recomputes fitted objectives and WHT summaries, and reconstructs paired latent fields. It reuses the verified archived observation kernel; it does not rerun the optimizer or whole-graph dynamic programming. `verify_artifacts.py` checks executed notebooks, their paired sources, embedded figure identities and report links, then writes the artifact manifest. `requirements.txt` lists the principal pinned packages and `environment.lock.txt` records the complete analysis environment.

No empirical adequacy threshold or aggregate ranking score is added. “Closer” counts smaller absolute error per paired seed with a 1e-10 numerical tie tolerance. All settings, missing denominators and both seed cohorts remain in the tables. Seed bands are variation across eight computational seeds, not biological confidence intervals.

## Files and storage

| Figure | Contents | Vector export |
|---|---|---|
| 1 | Individual paired seeds at the inherited setting | [PDF](figures/01_paired_inherited_regime.pdf) |
| 2 | All six interaction strengths | [PDF](figures/02_interaction_sweep.pdf) |
| 3 | All three observation models | [PDF](figures/03_observation_models.pdf) |
| 4 | Seed cohort and resolution sensitivity | [PDF](figures/04_sensitivity.pdf) |
| 5 | Calibration, confirmation spread and availability | [PDF](figures/05_calibration_and_availability.pdf) |

The analysis archive is `data/processed/aggregate-thermodynamic-gnk-calibration/ANA-070.dvc`. It contains the report, notebooks, helper source, environment lock, figures and tables. `artifact_manifest.json` records SHA-256 hashes; `storage_verification.json` verifies the DVC pointer, collected payloads and local cache. Metadata and storage sidecars are kept outside the archive to avoid circular hashes. Remote backup remains pending while Dutch is unavailable.
