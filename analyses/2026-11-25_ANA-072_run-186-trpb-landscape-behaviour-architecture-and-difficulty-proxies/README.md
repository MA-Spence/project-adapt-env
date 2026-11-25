# ANA-072 — RUN-186 TrpB landscape adequacy

Read [analysis.md](analysis.md) for the figure-led interpretation. This analysis
assesses average fitness behaviour, observable genetic architecture and structural
proxies for learning and optimisation difficulty. It performs no refitting,
learner training, optimisation task or new adaptive-walk benchmark.

| Executed notebook | Reproducible outputs |
|---|---|
| [validation.ipynb](notebooks/validation.ipynb) | Figure 7: the joint Optuna search; archive/source hashes, all trials, objective and model replay, historical-control agreement. |
| [analysis.ipynb](notebooks/analysis.ipynb) | Figures 1–2: native fitness distributions, local variation, near-floor mass, mutation magnitudes and diminishing returns. |
| [architecture.ipynb](notebooks/architecture.ipynb) | Figures 3–4: equal-cube and pooled interaction-order spectra; latent, intrinsic and observed epistasis. |
| [difficulty.ipynb](notebooks/difficulty.ipynb) | Figures 5–6: neighbouring-fitness correlation, local traps, reachability, resolution and neutral connectivity. |

Each notebook has a paired percent-format `.py` source. Helpers regenerate the
numerical [tables](tables/) and all seven [PNG/PDF figure pairs](figures/).
The notebooks contain their executed figures and checks.

## Reproduce from frozen inputs

Run from the project root with **Python 3.11.15**. The following restores the
complete DVC-tracked RUN archive, creates an isolated environment, builds all
four notebooks from source, executes them in fresh kernels and checks artifacts:

```sh
labproj data restore data/processed/aggregate-thermodynamic-gnk-calibration/RUN-186 --pull
python3.11 -m venv .labproj/envs/ana072-py311
.labproj/envs/ana072-py311/bin/python -c 'import sys; assert sys.version_info[:3] == (3, 11, 15), sys.version'
.labproj/envs/ana072-py311/bin/python -m pip install -r analyses/2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies/notebooks/environment.lock.txt
.labproj/envs/ana072-py311/bin/python analyses/2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies/notebooks/build_notebooks.py validation analysis architecture difficulty
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR="$PWD/.labproj/tmp/ana072-mpl" IPYTHONDIR="$PWD/.labproj/tmp/ana072-ipython" PYTHONDONTWRITEBYTECODE=1 .labproj/envs/ana072-py311/bin/python analyses/2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies/notebooks/execute.py validation analysis architecture difficulty
.labproj/envs/ana072-py311/bin/python analyses/2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies/notebooks/check_artifacts.py
```

The configured DVC remote must be reachable for `--pull`; no pull is needed when
the complete archive is already present. Installing dependencies requires package
access. The scientific analysis itself downloads nothing and does not need Slurm.
The executor binds temporary Jupyter kernel specifications to the invoked Python;
local Jupyter loopback sockets must be allowed.

[environment.lock.txt](notebooks/environment.lock.txt) records the complete
macOS execution environment. [requirements.txt](notebooks/requirements.txt) gives
the direct scientific/notebook package pins for another compatible environment;
that shorter list does not lock all transitive dependencies. Original execution
used the existing `ana070-py311` interpreter, with Python 3.11.15 and the versions
in [analysis_environment.json](tables/analysis_environment.json). The environment
name does not select a different RUN or scientific input. Cross-platform replay
allows the documented floating-point tolerances; it does not require identical
rendered PNG or notebook bytes.

To rebuild a single notebook, pass its stem to `build_notebooks.py` and then to
`execute.py`. Rebuilding clears its saved execution, so rerun it before checking.
After scientific-source changes, rerun all affected notebooks and the checker.
For a final read-only check of an archived analysis, append `--check` to the
checker command. That validates current artifacts against the recorded hashes
without overwriting `tables/artifact_validation.json`.

## Inputs and boundaries

The [RUN-186 DVC pointer](../../data/processed/aggregate-thermodynamic-gnk-calibration/RUN-186.dvc)
restores outputs, execution logs, collection checks and the complete frozen runtime.
Its `frozen_runtime/source` includes the model/observation code, empirical fitness
and genotype mask, profile, calibration and confirmation panels, donor assignments,
original configuration, historical RUN-178/RUN-180 control data and source manifests.
Nothing is taken from a currently edited package or reconstructed by refitting.
The [RUN record](../../experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-186.yaml)
and [joint-fitting protocol](../../experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/joint_readout.md)
retain the original execution and calibration design.

The validation notebook checks collection hashes, all 86 manifest-listed outputs
and all 339 frozen runtime payloads. It replays objectives and all 24 model fields.
The average/difficulty reconstruction restores fixed observation maps from the
recorded donors and seeds and checks all 120 maps against archived cube values.
Observation error is static by genotype in these archived maps; it is not newly
sampled on each oracle query. Intrinsic model values remain separate. Shared-stop
error is a sensitivity assumption, not an identified biological mechanism.

The main figures show the eight originally fresh seeds, which were already used
in RUN-180; tables retain all 24. Seed ranges are not confidence intervals, and
local-cube ranges also reflect archived panel assignments. Whole-library results
reuse known measured genotypes, including calibration nodes. Confirmation cubes
were excluded from fitting but previously inspected. No independent biological
replication, specific interaction-identity recovery or ML-policy transfer is
established. Static graph proxies do not measure actual learner/search performance.

[validation_checks.json](tables/validation_checks.json),
[average_reconstruction.json](tables/average_reconstruction.json),
[architecture_validation.json](tables/architecture_validation.json) and
[difficulty_verification.json](tables/difficulty_verification.json) retain numerical
checks. `tables/artifact_validation.json` records notebook/source parity, execution,
figure/table coverage and report-link checks. These checks concern reproducibility,
not a numerical biological-adequacy acceptance gate.
