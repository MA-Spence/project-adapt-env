# Joint continuous-r TrpB calibration with Optuna

Authorized 2026-11-24 under EXP-065. The scientist requested removal of the r
grid, joint floating-point optimization on [0,8], and a new Slurm RUN. The
current project selector is `joint-readout`; the Slurm worker uses the same
frozen source through `joint-readout-slurm`. Gamma concentration remains disabled
(`interaction_concentration=None`); alpha is not fitted.

## Scientific design

Retain RUN-180's profile, full Gaussian 20 × 20 pair tables, six site pairs,
reference-anchored sigmoid, native fitness observations, calibration and
confirmation cube partitions, generator seeds, observation donors and fixed
noise draws, affine units, readout bounds, genotype mask and unfitted diagnostics.
The exact loss remains the [restored direct loss](optuna_readout.md):
101 quantiles including endpoints, variance-normalized quantile discrepancy,
relative mean within-cube SD discrepancy, and near-floor mass discrepancy.
Each term has unit weight. No auxiliary coefficient inference or bootstrap
weighting is introduced.

Run three independent serial Optuna TPE studies: `none`,
`independent_difference`, and `shared_stop_sensitivity`. These remain separate
observation assumptions; their fitness outputs are not pooled or selected using
confirmation. Every study jointly optimizes five floating-point coordinates:

| Coordinate | Search range |
|---|---|
| Interaction RMS ratio r | [0,8], continuous linear sampling |
| Standardized output offset | [-2,2] |
| Log standardized span | log([0.001,10]) |
| Reference sigmoid position | [-14,8] |
| Log sigmoid steepness | log([0.05,8]) |

Use 11,712 trials per mode, including 20 startup trials and six historical warm
starts. This retains the prior planned total of 35,136 trials (6 × 1,952 × 3).
The six archived (r, readout) fits per mode are authenticated and re-evaluated
before being queued as candidates. They are not a fitting grid. Subsequent r
proposals use `suggest_float('r', 0, 8)` with no step or categorical choice.
The selected calibration loss cannot exceed the best re-evaluated historical
candidate for that observation mode. This is not a confirmation guarantee or
an independent comparison of optimizer efficiency.

Each study returns one fitted r/readout combination, giving **three selected
models**, replacing the old 18 fixed-r settings. Confirm these against the same
32 withheld empirical cubes, 16 reused and eight originally fresh generator
seeds. Evaluate whole-library graph and cube summaries at the selected values;
latent graph summaries use the unique selected r values. The full library mask
contains 159,129 measured genotypes. Whole-library graphs are post-fit diagnostics,
not fully held-out data. Previously inspected confirmation panels are not new
independent biological validation. Shared-stop remains a sensitivity assumption.

## Execution and provenance

Configuration: [joint_readout.yaml](joint_readout.yaml). Source hashes:
[joint_readout_source.json](joint_readout_source.json). Source and input payloads
are frozen in `data/processed/trpb-joint-readout-slurm-runtime` before submission.
The prelaunch manifest records every payload; the worker authenticates all bytes
and exact dependency pins before fitting. Historical RUN-180 and earlier controls
are reused and remain immutable. Trial databases, per-trial CSV histories,
all five coordinates, separate loss components and warm-start provenance are saved.
No study resumes from an old database. Budget completion does not imply convergence.

Submit through labproj's `lab-slurm` target via Dutch, using a dedicated Python
3.11.15 environment with the [readout dependency lock](../../environment/readout-py311.lock.txt).
Request four CPUs, 8 GB RAM and 24 hours on the CPU partition. Each Optuna study
runs serially with one numerical-library thread; four threads evaluate post-fit
confirmation seeds. Use the existing Zig C++ graph compiler and record its version
and checksum. Linux x86 replaces RUN-180's Mac ARM execution; numerical replay
and actual numerical-library configuration are checked and recorded before launch.
The initial remote preflight measured a maximum calibration interaction-field
difference of 1.33e-15. The joint protocol accepts field replay at rtol=1e-13
and atol=1e-14, records actual differences, and retains the historical objective
checks at rtol=1e-12 and atol=1e-14. Generated values are neither rounded nor
replaced. Exact integer references, input bytes and dependency pins remain required.

`--preflight` replays all 18 old calibration objectives without fitting.
`--quick` is an eight-trial/two-confirmation-seed engineering smoke mode, labelled
`SMOKE_ONLY`. Neither replaces the registered scientific protocol. The submission
helper defaults to a nonexecuting labproj preview; `--execute` submits a new RUN.
