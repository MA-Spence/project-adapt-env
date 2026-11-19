# ANA-046: RUN-121 PHOT RES-045 calibration isolation review

## Purpose

Review RUN-121 compact outputs to determine whether the cumulative RES-045 isolation stages repair PHOT raw scale, accessibility, local maxima, spectral texture, sequence-map, and search failures.

## Linked experiments/runs

- Experiments: EXP-047
- Runs: RUN-121

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-121_summary.json`
- `source_tables/run-121_preflight.json`
- `source_tables/run-121_remote_output_manifest.tsv`
- `source_tables/run-121_bounded_gate_candidates.csv`
- `source_tables/run-121_selected_bounded_gate_candidates.csv`
- `source_tables/run-121_bounded_gate_distances.csv`
- `source_tables/run-121_bounded_gate_equivalence.csv`
- `source_tables/run-121_graphlet_summary.csv`
- `source_tables/run-121_graphlet_metrics.csv`
- `source_tables/run-121_abc_comparator_selected_family_distances.csv`
- `source_tables/run-121_softmin_capability_probe.csv`
- `source_tables/run-121_graphlet_catalog.csv`
- `source_tables/run-121_module_catalog.csv`
- `tables/metrics.json`
- `tables/run-121_stage_equivalence_decision.md`
- `tables/run-121_stage_best_test_distances.md`
- `tables/run-121_bounded_gate_test_distances.md`
- `tables/run-121_abc_comparator_test_distances.md`
- `tables/run-121_top_component_distances.md`
- `tables/run-121_observable_vs_raw_epistasis_and_landscape.md`
- `tables/run-121_raw_latent_pathology_diagnostics.md`
- `tables/run-121_softmin_capability_probe_summary.md`
- `tables/run-121_best_candidate_parameters.md`
- `tables/run-121_stage_failure_modes.md`
- `tables/run-121_output_collection_audit.md`
- `tables/run-121_graphlet_and_module_catalog_summary.md`
- `tables/run-121_prior_run_comparison.md`

## Collection policy

`RUN-121` completed on `lab-slurm` as Slurm job `1918` with scheduler exit
code `0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-121/outputs`

This analysis copied the compact JSON/CSV outputs needed to interpret the
experiment. The copied source payload is `3,820,656` bytes and is audited in
`tables/run-121_output_collection_audit.md`. The detailed coefficient tables
were not copied because their aggregate spectral summaries are already present
in `graphlet_summary.csv` and `graphlet_metrics.csv`: `gft_coefficients.csv`
(`29.2 MB`) and `walsh_terms.csv` (`28.5 MB`) remain on the runtime volume.
The PHOT alignment FASTA also remains on the runtime volume because it is not
needed for the stage-level decision.

## Output audit

`RUN-121` evaluated `24` PHOT graphlets, `415` sequence-conditioned candidate
parameter sets across five fitted cumulative stages, selected `32` candidates,
and produced `32` equivalence rows against the best available `RUN-098`
full-ABC comparator. Stage `A+F+B` was intentionally probe-only and therefore
has no selected held-out candidate rows; it produced `48` model-free k-trait
soft-min capability probe rows.

The comparator used for the formal decision is the same best available
`RUN-098` row used in recent PHOT bounded-gate analyses:
`empirical_readout_residual` /
`latent_phi__minus_poly_fitted_current_best`, with test distance `0.2115`.
For context, the best nested full-ABC family row in the copied comparator table
is `random_module`, with test distance `0.4826`.

## Main result

No EXP-047 stage passed operational equivalence. `RUN-121` reports `0/32`
operationally equivalent rows, and no selected candidate passed the absolute
or relative aggregate-distance margin against the best available `RUN-098`
full-ABC comparator.

The best `RUN-121` held-out row came from stage `A+F+B+C+D`,
`A+F+B+C+D__local_refinement_round_1__0090109041`, with test distance
`0.9155`. This is a `25.2%` improvement over the best `RUN-119` raw-latent
path-forward row (`1.2240`), but it is still worse than the best `RUN-117`
cumulative row from `RES-044` (`0.4122`), worse than the nested full-ABC
context row (`0.4826`), and far from the formal best available full-ABC
comparator (`0.2115`).

## Stage interpretation

Stages `A` and `A+F` were negative. Generator-scale anchoring alone gave a
best held-out distance of `1.7127`; adding the scale-invariant geometry plus
raw single-mutant-scale objective reduced that only slightly to `1.6716`. Both
stages had `0` component guardrail passes, and their failures were dominated
by plaquette-energy concentration, sign epistasis, single-mutant SD, and
magnitude-epistasis components.

The model-free `A+F+B` capability probe shows that per-graphlet additive
k-trait soft-min fits can produce moderate direct fit correlations but not the
target local ruggedness package. Across k values `2`, `3`, `5`, and `8`, the
median fit correlation stayed around `0.86` to `0.88`, but local maxima reached
only `19%` to `35%` of the target and magnitude epistasis reached only `14%`
to `30%` of the target. This argues against treating the sequence-conditioned
search as the only problem: the tested additive soft-min family is already
limited in the model-free probe.

Stage `A+F+B+C` was the first substantial sequence-conditioned improvement.
Re-enabling fitting with k-trait and auxiliary-texture sweeps reduced the best
held-out distance to `1.0124`. Stage `A+F+B+C+D`, which adds sparse structured
residual texture and a restricted sigmoid readout, gave the best held-out row
at `0.9155`. The expanded `A+F+B+C+D+E` search lowered the calibration
objective to `0.4884`, but its best held-out row was slightly worse (`0.9298`)
than the unexpanded structured-residual stage.

## Failure diagnosis

EXP-047 partially fixes the specific raw-scale pathology from `RES-045`, but
not in a scientifically usable way. In the later fitted stages, raw
single-mutant SD and raw magnitude epistasis are close to the empirical PHOT
scale (`0.97x` to `1.25x` and `0.93x` to `1.02x`, respectively). However,
raw accessibility remains poor: the best stage reaches only `32%` of the
empirical accessible-path fraction. The expanded-search stage reaches `85%`
of empirical local maxima but collapses raw accessibility to `4%` of the
target.

The more serious failure is observable readout collapse. The stage-best rows
have observable single-mutant SD only `0.08%` to `0.16%` of empirical and
observable magnitude epistasis only `0.04%` to `0.14%` of empirical. Thus the
later stages can put raw latent scale near the empirical range only while
compressing the observable landscape almost flat. That is incompatible with
the PHOT calibration goal, even before considering the failed held-out
equivalence margins.

## Decision

`RUN-121` does not validate the RES-045 calibration-isolation path as a PHOT
repair. It improves over the failed `RUN-119` raw-latent soft-min diagnostic
and identifies a useful boundary: raw scale can be brought near the empirical
range under the tested k-trait/structured-residual/sigmoid family. But the
same family still fails held-out graphlet equivalence, fails every selected
component-equivalence decision, undershoots accessibility, and collapses
observable variation.

The current PHOT adequacy claim remains weakened. For benchmark construction,
EXP-047 is methodological evidence only: it rules out this tested
calibration-isolation route as a validated matched-landscape construction
route, but it does not test optimizer or PLM-agent behavior.
