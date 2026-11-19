# ANA-045: RUN-119 PHOT raw-latent soft-min pathology fix review

## Purpose

Review RUN-119 compact outputs to determine whether the three cumulative EXP-046 RES-044 path-forward fixes repair PHOT raw-latent accessibility and scale pathologies while retaining observable calibration.

## Linked experiments/runs

- Experiments: EXP-046
- Runs: RUN-119

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-119_summary.json`
- `source_tables/run-119_preflight.json`
- `source_tables/run-119_remote_output_manifest.tsv`
- `source_tables/run-119_bounded_gate_candidates.csv`
- `source_tables/run-119_selected_bounded_gate_candidates.csv`
- `source_tables/run-119_bounded_gate_distances.csv`
- `source_tables/run-119_bounded_gate_equivalence.csv`
- `source_tables/run-119_graphlet_summary.csv`
- `source_tables/run-119_graphlet_metrics.csv`
- `source_tables/run-119_abc_comparator_selected_family_distances.csv`
- `source_tables/run-119_graphlet_catalog.csv`
- `source_tables/run-119_module_catalog.csv`
- `tables/metrics.json`
- `tables/run-119_stage_equivalence_decision.md`
- `tables/run-119_stage_best_test_distances.md`
- `tables/run-119_bounded_gate_test_distances.md`
- `tables/run-119_abc_comparator_test_distances.md`
- `tables/run-119_top_component_distances.md`
- `tables/run-119_observable_vs_raw_epistasis_and_landscape.md`
- `tables/run-119_raw_latent_pathology_diagnostics.md`
- `tables/run-119_best_candidate_parameters.md`
- `tables/run-119_stage_failure_modes.md`
- `tables/run-119_output_collection_audit.md`
- `tables/run-119_graphlet_and_module_catalog_summary.md`

## Collection policy

`RUN-119` completed on `lab-slurm` as Slurm job `1917` with inferred exit code
`0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-119/outputs`

This analysis copied the compact JSON/CSV outputs needed to interpret the
experiment. The copied source payload is `2,516,588` bytes and is audited in
`tables/run-119_output_collection_audit.md`. The detailed coefficient tables
were not copied because their aggregate spectral summaries are already present
in `graphlet_summary.csv` and `graphlet_metrics.csv`: `gft_coefficients.csv`
(`19.3 MB`) and `walsh_terms.csv` (`18.8 MB`) remain on the runtime volume.
The PHOT alignment FASTA also remains on the runtime volume because it is not
needed for the stage-level decision.

## Output audit

`RUN-119` evaluated `24` PHOT graphlets, `219` candidate parameter sets across
the three cumulative stages, selected `24` candidates, and produced `24`
equivalence rows against the best available `RUN-098` full-ABC comparator.

The comparator used for the formal decision is the best available `RUN-098`
row in the copied comparator distance table:
`empirical_readout_residual` /
`latent_phi__minus_poly_fitted_current_best`, with test distance `0.2115`.
For context, the best nested full-ABC family row in the same table is
`random_module`, with test distance `0.4826`.

## Main result

No EXP-046 stage passed the operational-equivalence decision. `RUN-119`
reports `0/24` operationally equivalent rows. No selected candidate passed the
absolute or relative aggregate-distance margin against the best available
`RUN-098` full-ABC comparator.

The best `RUN-119` held-out row came from stage `1`,
`1__local_refinement_round_1__4052355884`, with test distance `1.2240`. This
is worse than the best `RUN-117` cumulative A-through-D row from `RES-044`
(`0.4122`), worse than the nested full-ABC context row (`0.4826`), and far
from the formal best available full-ABC comparator (`0.2115`).

## Stage interpretation

Stage `1`, the low-dimensional two-additive-trait soft-min latent generator
with the working EXP-045 observable conditioning, was the best EXP-046 stage
but still failed. Its best held-out row reached raw accessibility `0.0796`,
about `40%` of empirical PHOT (`0.1975`), and raw local maxima `2.33`, about
`45%` of empirical (`5.17`). That partial topology movement did not solve
scale: raw single-mutant SD was `252.5`, about `934x` empirical, and raw
magnitude epistasis was `124.1`, about `424x` empirical. Observable
magnitude epistasis was only `0.0647`, about `22%` of empirical.

Stage `1+2`, which added raw-latent objective constraints and restricted the
observation readout to a fitted monotone sigmoid, was worse. Its best held-out
distance was `1.6872`; raw accessibility rose only to `0.0958` (`49%` of
empirical), while raw local maxima fell to `1.58` (`31%` of empirical).
The raw scale remained too large: raw single-mutant SD was `13.44` (`50x`
empirical) and raw magnitude epistasis was `7.01` (`24x` empirical).
Observable variation nearly collapsed, with observable single-mutant SD only
`4.9%` of empirical and observable magnitude epistasis only `2.3%` of
empirical.

Stage `1+2+3`, which up-weighted raw accessibility as a first-class target,
did not rescue the fit. The selected best held-out row and aggregate
diagnostics were effectively identical to stage `1+2`: test distance `1.6872`,
raw accessibility `0.0958`, raw local maxima `1.58`, and the same raw-scale
and observable-collapse pattern. In this sampled fast-fit search, increasing
the raw-accessibility weight did not find a better PHOT-matched latent
geometry.

## Decision

`RUN-119` does not validate the three RES-044 path-forward fixes as a PHOT
calibration repair. The low-dimensional soft-min latent generator partially
moves raw accessibility toward the empirical target, but it leaves raw scale
pathologically large and the held-out graphlet distance poor. Adding
raw-latent constraints, a restricted sigmoid readout, and higher
raw-accessibility weight worsens held-out fit and collapses observable
variation without matching raw topology.

The result sharpens the failure diagnosis. EXP-045 showed that an observable
readout can improve magnitude and ratio summaries while leaving latent
geometry unusable. EXP-046 shows that directly targeting raw latent topology
with the tested low-dimensional soft-min generator and fast-fit constraints
does not produce a usable PHOT landscape either. The failure is now both
topological and scale-calibration: raw accessibility remains below empirical,
raw local maxima are too few, and raw epistatic/single-mutant scales remain
orders of magnitude too large.
