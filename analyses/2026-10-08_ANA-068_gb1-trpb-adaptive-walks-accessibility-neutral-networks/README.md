# ANA-068 — Whole-library GB1 and TrpB dynamics

Open [notebooks/analysis.ipynb](notebooks/analysis.ipynb) for the executed analysis
and six figures; [analysis.md](analysis.md) provides the interpretation.

```bash
python -m pip install -r notebooks/requirements.txt
python notebooks/execute.py
```

Run from this ANA directory or use the full path to the executor. It selects the
invoking interpreter as a temporary Jupyter kernel. A C++17 compiler is required:
macOS Command Line Tools (`clang++` and `xcrun`) or Linux `g++`. The native graph
kernel is compiled into the operating-system temporary directory, never the data
bundle. It is tested against independent Python/SciPy small-graph calculations.
The notebook recomputes the full analysis; typical calculation time is a few minutes.

Restore `data/processed/production_landscapes/ANA-068.dvc` and `RUN-172.dvc` first
if absent. The ANA bundle includes the two original empirical parquet files;
RUN-172 supplies the frozen source, parameters, seeds and TrpB noise-calibration
provenance. Existing processed input files can bootstrap the ANA bundle on its
first execution. No optimization or job submission occurs in the notebook.

`protocol.md` defines the estimands and scope. `analysis_tools.py` implements
calculation and verification; `graph_kernel.cpp` provides implicit-adjacency graph
algorithms; `plot_tools.py` renders figures; `analysis.py` is the notebook's paired
text source. `tables/` contains exact per-seed results and audits. Error bars and
bands represent variation across 16 fixed model seeds, not biological confidence
intervals. Every scientific comparison remains descriptive, without success gates.
