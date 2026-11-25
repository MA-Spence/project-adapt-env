Project-specific overlay environment guidance.

## RUN-180 readout rollback with Optuna (2026-11-24)

`readout-py311.lock.txt` records the Python 3.11.15 engineering environment for
the restored EXP-065 readout workflow. It retains RUN-180's numerical package
versions and adds Optuna 4.9.0 and its dependencies. Local checks use the existing
`indirect-py311` interpreter with `trpb-py31115/lib/python3.11/site-packages` on
`PYTHONPATH`; the two historical environments are unmodified. A clean environment
can install this combined lock for subsequent registered execution. The runner
checks Python, NumPy, SciPy, pandas, PyArrow, PyYAML and Optuna versions and records
the optimizer storage dependencies. This does not refresh the historical overlay.

## Architecture-calibration implementation (2026-11-20)

The new optional package extra `adaptenv[calibration]` uses Optuna 4.x. Engineering
validation uses separate local environments; the completed RUN-180 environment
was not modified. `calibration-py311.lock.txt` records the Python 3.11.15 fitting
test environment. `calibration-py313.lock.txt` records the Python 3.13.11 dependency
closure used for full source and installed-wheel checks, including the original
historical fingerprint fixtures. These are validation environments, not a new
registered biological RUN or a replacement for the historical overlay lock.

Exact ranked fingerprints depend on Python's float-summation semantics. All 28
checked old maps still match frozen RUN-180 source in the same Python 3.11 process;
the original fixtures pass unchanged under Python 3.13. Use recorded interpreters
for exact historical replay. Candidate verification excludes ambient editable
import hooks explicitly. See `docs/model_release/0.7.0rc1`.

## Current validity boundary (2026-11-19)

`requirements.in` records accumulated research dependencies. The existing
`requirements.lock.txt` is a historical overlay snapshot, **not a complete lock
of the current input**: it predates the Landscapy install entry and its dependency
closure. Its labproj header records current external source identities only.
The configured STACK-000 Python 3.11 also does not satisfy the pinned Landscapy
checkout's Python >=3.12 requirement. Refreshing metadata does not resolve either
issue. Do not use this overlay as a validated environment for the current model.

The maintained model is independently installable from `external/Adapt-Env`.
RES-063 records historical source and isolated-wheel tests under Python 3.13.
The current full-pairwise candidate's checks and source lock are in
`docs/model_release/0.6.0rc1`; `0.5.0rc1` retains the prior ranked model's evidence.
Reproduce historical experiments with their recorded RUN environments and source
snapshots. Rebuilding the full research overlay requires an explicit compatible
base stack and dependency resolution; no broad dependency upgrade was performed
as part of this source and scientific-record synchronization.

Rules:
1. Use the approved base stack by default.
2. Use requirements.in only for small project-specific additions.
3. Regenerate requirements.lock.txt after changing requirements.in.
4. Do not add GPU frameworks such as torch/JAX/CUDA here.
5. Do not commit virtual environments.
6. If a dependency is needed by multiple projects, request a registered stack update.
