# ANA-055: RUN-141 design-matched support-aware GNK probe review

## Purpose

Collect and interpret RUN-141 fixed Stage 4B support-aware probe outputs against the scratch RES-054 core-vs-GNK separation and monotonicity gates.

## Linked experiments/runs

- Experiments: EXP-057
- Runs: RUN-141

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-141_summary.json`
- `source_tables/run-141_preflight.json`
- `source_tables/run-141_checkpoint_resume_summary.json`
- `tables/metrics.json`
- `tables/run-141_output_collection_audit.md`
- `tables/run-141_support_template_audit.md`
- `tables/run-141_anchor_regime_audit.md`
- `tables/run-141_monotonicity_by_metric.md`
- `tables/run-141_monotonicity_by_support_anchor.md`
- `tables/run-141_monotonicity_by_dataset.md`
- `tables/run-141_core_gnk_gate_by_metric.md`
- `tables/run-141_core_gnk_gate_by_epi_weight.md`
- `tables/run-141_core_gnk_gate_by_dataset.md`
- `tables/run-141_full_core_deltas_by_epi_weight.md`

## Collection policy

`RUN-141` completed on `lab-slurm` as Slurm job `1955` with inferred scheduler
exit code `0`. The compact runtime outputs were collected locally and
DVC-tracked at:

- `data/processed/stage-4b-design-matched-support-aware-gnk-probe/RUN-141/outputs.dvc`
- `data/processed/stage-4b-design-matched-support-aware-gnk-probe/RUN-141/logs.dvc`

The checkpoint tree was not copied into the scientific record. It remains on
the runtime host at
`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-141/outputs/checkpoints/generator_parameters`.

## Execution audit

The final attempt computed all `135/135` parameter sets, emitted `35,640`
generator descriptor rows, and wrote `8,910` support-aware graphlets. The local
log audit found `135` stdout lines, empty stderr, and exit code `0`.

`RUN-141` is therefore a completed run, not a failed execution. The failure is
scientific: the fixed probe did not pass the preregistered gates required before
re-running overlap, joint-realism, or held-out tests.

## Probe audit

The experiment did implement the main measurement fixes requested in
`scratch/RES-054_assessment_and_path_to_finishing_component1.md#3`:

- support-aware graphlets were used for all graphlets (`support_aware_fraction = 1.0`);
- dataset dimensions and alphabets were matched to the empirical templates;
- anchors covered unfolded, fold-transition, and plateau regimes with small
  average margin error (`0.0115`, `0.0151`, and `0.0137` respectively);
- topography was reported through strict local-maxima density and peak-targeted
  reachability.

This means the negative gate result should no longer be dismissed as the old
first-positions probe artifact from `RES-051`.

The empirical template descriptor rows were all marked
`skipped_incomplete_template`, so this run should not be read as a direct
empirical-overlap experiment. It is a fixed Stage 4B generator probe.

## Monotonicity gate

Only `318/1620` monotonicity checks passed (`0.1963`). The best individual
metric was sign-epistasis fraction, with `119/270` checks passing (`0.4407`).
The metrics most central to the Stage 4B public-generator gate did not pass:

- high-order fraction: `36/270` (`0.1333`), median Spearman `-0.1`;
- local-maxima density: `28/270` (`0.1037`), median Spearman `-0.1817`;
- random-walk lag-1 autocorrelation: `37/270` (`0.1370`);
- peak-targeted reachability: `58/270` (`0.2148`).

Pass rates were not rescued by routing mode or anchor regime. Structured routing
was only slightly higher than random or shuffled routing (`0.2111` versus
`0.1889` for both controls), and plateau anchors were still below a majority
pass rate (`0.2389`).

## Core-vs-GNK separation gate

The hard preregistered gate failed. Across positive epistasis weights, only
`34/1080` rows passed core-vs-GNK separation (`0.0315`).

All passing rows came from the high-order-fraction gate. The other two
preregistered separation components had no positive rows:

- high-order fraction: `34/1080` (`0.0315`);
- sign-epistasis fraction: `0/1080`;
- local-maxima density: `0/1080`.

The separation also did not become convincing at high epistasis weight. Pass
rates by `epi_weight` were `0/270` at `0.15`, `2/270` at `0.30`, `16/270` at
`0.60`, and `16/270` at `1.00`.

The effect sizes were correspondingly small. At `epi_weight = 1.0`, the mean
absolute full-minus-core high-order-fraction delta was `0.00825`, below the
`0.03` gate threshold. The sign-epistasis mean absolute delta was `0.000743`,
and the local-maxima-density mean absolute delta was `0.0000169`.

## Decision

`RUN-141` does not restore protein-scale monotone knob control or public
full-vs-core separation for the current additive thermodynamic core plus GNK
implementation.

The correct interpretation is narrower than "GNK is biologically wrong" but
stronger than "the old probe was flawed." The fixed, support-aware probe now
sees GNK supports and still finds that the public observable rarely separates
from the additive core by the preregistered effects. That points to a
core/GNK scale-balance or readout-coupling failure in the current public
implementation.

Per the EXP-057 stop rule, downstream Stage 4C/4D/4E overlap, joint-realism,
and held-out reruns remain blocked for this implementation.
