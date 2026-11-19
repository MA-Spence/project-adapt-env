Project-specific overlay environment guidance.

## Current validity boundary (2026-11-19)

`requirements.in` records accumulated research dependencies. The existing
`requirements.lock.txt` is a historical overlay snapshot, **not a complete lock
of the current input**: it predates the Landscapy install entry and its dependency
closure. Its labproj header records current external source identities only.
The configured STACK-000 Python 3.11 also does not satisfy the pinned Landscapy
checkout's Python >=3.12 requirement. Refreshing metadata does not resolve either
issue. Do not use this overlay as a validated environment for the current model.

The maintained model is independently installable from `external/Adapt-Env`.
RES-063 records successful source and isolated-wheel tests under Python 3.13;
the release candidate's source lock is in `docs/model_release/0.5.0rc1`.
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
