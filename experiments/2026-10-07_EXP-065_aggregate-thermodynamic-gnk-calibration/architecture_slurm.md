# TrpB architecture calibration: Slurm execution amendment

On 2026-11-21 the scientist requested stopping local RUN-182 and submitting the
same calibration to the remote Slurm queue. RUN-182 was interrupted with SIGTERM
(runner exit 143); its worker and supervisor exited. Its partial collection is
preserved unchanged, including the original automatic failure classification.
The current RUN record clarifies user cancellation. No completed calibration or
numerical failure is inferred from missing final artifacts.

The remote execution restarts from the beginning in a new RUN under EXP-065.
Optuna resumption is not used because the implementation does not restore the
sampler RNG state. The scientific method remains
`thornton_auxiliary_optuna_indirect_inference_v1`, described in the original
[protocol](architecture_calibration.md).

## Frozen scientific inputs and settings

`data/processed/trpb-architecture-slurm-runtime` contains the original RUN-182
prelaunch manifest and its 152 payloads, copied byte-for-byte: source, configuration,
inputs and original Mac environment evidence. The source manifest SHA-256 is
`a93780a166d61ed95919376ce7a7ce7e70bff015fbd0fb728fc411bbd005bbe7`;
the input manifest SHA-256 is
`3762eec819d6777e37d995957ec2e34182f12c1ae69fa3bb7953ba76f358a1b4`.

The model, profile, both 8,192-genotype panels and actual assay replicates,
100 Optuna trials, 20 startup trials, 20 bootstrap fits, four calibration and four
confirmation seeds, search bounds, inference penalty grid, grouped folds,
convergence rules and all seeds are unchanged. No confirmation result is used
to prepare this execution. Historical controls are not rerun.

## Execution changes

The remote interpreter is a dedicated Python 3.11.15 environment at
`/tank/data/labproj/runtime/project-adapt-env/envs/architecture-0.7.0rc1-py311`.
It uses the exact original dependency-version pins. Linux x86_64 and its BLAS
implementation replace macOS ARM and its numerical libraries. Identical scientific
inputs do not guarantee bitwise-identical optimization trajectories across these
platforms. Actual remote versions, platform and numerical libraries are recorded.

The job requests the existing CPU partition, one CPU, 8 GB and 24 hours, with
numerical libraries restricted to one thread. The Slurm walltime is enforced;
the earlier local 24-hour request was only descriptive. No model simplification
or budget reduction is introduced. A timeout must be reported as incomplete.

The submission helper routes through Dutch because the gateway lacks the shared
runtime mount. It stages only required transport files and the self-contained
frozen runtime. Its temporary dispatcher explicitly selects `architecture-slurm`
so the selector survives SSH submission. Temporary project configuration, stage
rules and dispatcher changes are restored after submission. The queued worker
verifies the original manifest and dependency pins before invoking the original
frozen scientific CLI. No dependency installation occurs inside the job.

Slurm currently reports Dutch DOWN. Submission is authorized while the job waits
in the queue; this action does not authorize restarting cluster services or
changing node state. Scheduler acceptance, pending reason and remote staging
checks belong in the new RUN record. Queue acceptance is not execution completion.
