# EXP-065 local execution of the frozen TrpB profile repeat

On 19 November 2026 the scientist requested execution on the local MacBook,
retaining the scientific record, while Dutch remains unavailable. This authorizes
the local target; it does not resume the held Slurm submission. RUN-176 remains
the earlier nonexecuting Slurm preview. New local preview and execution RUN IDs
are created through labproj under EXP-065.

The scientific protocol remains `alignment_prior.md` and `alignment_prior.yaml`.
No source in `alignment_prior_source.json`, input, seed, optimizer setting,
observation assumption, fitted parameter or diagnostic is changed. This is the
full calibration, not a shortened smoke run. The existing RUN-175 synthetic
control remains the comparator and is not recalibrated.

## Execution environment

An isolated project environment at `.labproj/envs/trpb-py31115` uses Python
3.11.15, NumPy 2.4.4, pandas 2.3.0, SciPy 1.17.1, PyArrow 20.0.0 and PyYAML
6.0.3, exactly matching the recorded numerical version guard. The actual local
platform is macOS on arm64. The graph source is compiled with the local Apple
Clang toolchain rather than the Linux worker's project-local Zig toolchain.
Compiler, platform, installed distributions and numerical-library configuration
are saved with the run. Fresh kernel validation checks the compiled graph code.

Matching versions and algorithmic settings is not a claim of bit-for-bit
cross-platform equality. Floating-point libraries and compiler differences may
change final digits or optimization trajectories. The original protocol's
scientific intervention is unchanged; the execution platform is an explicitly
recorded reproducibility qualification.

The local labproj template does not activate STACK-000. `PYTHON` therefore
explicitly selects the isolated environment, while the actual environment
record takes precedence over the nominal stack label. The runner uses four
workers and one BLAS/OMP thread per worker. Local CPU, memory and walltime
metadata are requests, not operating-system enforced resource limits.

## Record and collection

`submit_alignment_prior_local.py` verifies all frozen source and input hashes,
saves a pre-execution source/input snapshot, and uses `labproj submit EXP-065
--target local --dry-run` followed by a separate `--execute` invocation.
The scientific entry point and configuration remain unchanged. Temporary
submission settings are restored on exit.

The local labproj executor runs synchronously and stores captured stdout/stderr
inside its RUN YAML. It does not stage a source snapshot or populate the log
files itself. The helper preserves the captured streams in the declared log
paths and manually collects verified outputs, logs, the pre-execution snapshot,
actual environment record and RUN sidecars into
`data/processed/aggregate-thermodynamic-gnk-calibration/RUN-ID/`.
This explicit collection is required because the current `labproj collect`
path mapper expects remote-storage settings absent from the local target.

Completion requires the scientific return code to be zero, the full-budget
manifest to identify `MATCHED_PRIOR_REPEAT`, all expected settings and seeds,
all declared outputs, and matching hashes for every manifest-listed file.
The collected bundle is then tracked with `labproj data track --no-backup`;
remote backup remains pending while Dutch is offline. Technical completion
does not constitute an empirical adequacy conclusion or a new CLAIM/RES.
