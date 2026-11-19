# ANA-044: RUN-117 PHOT cumulative epistasis calibration fix review

## Purpose

Review RUN-117 compact outputs to determine whether the four cumulative
EXP-045 PHOT calibration stages, `A`, `A+B`, `A+B+C`, and `A+B+C+D`, repair
the epistasis calibration and magnitude failures diagnosed after `RES-043`.
The formal decision compares each selected current-model candidate with the
best available `RUN-098` full-ABC PHOT comparator, using the same held-out
graphlet-equivalence framing as recent PHOT records.

## Linked experiments/runs

- Experiments: EXP-045
- Runs: RUN-117

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-117_summary.json`
- `source_tables/run-117_preflight.json`
- `source_tables/run-117_remote_output_manifest.tsv`
- `source_tables/run-117_bounded_gate_candidates.csv`
- `source_tables/run-117_selected_bounded_gate_candidates.csv`
- `source_tables/run-117_bounded_gate_distances.csv`
- `source_tables/run-117_bounded_gate_equivalence.csv`
- `source_tables/run-117_graphlet_summary.csv`
- `source_tables/run-117_graphlet_metrics.csv`
- `source_tables/run-117_abc_comparator_selected_family_distances.csv`
- `source_tables/run-117_graphlet_catalog.csv`
- `source_tables/run-117_module_catalog.csv`
- `tables/metrics.json`
- `tables/run-117_stage_equivalence_decision.md`
- `tables/run-117_stage_best_test_distances.md`
- `tables/run-117_bounded_gate_test_distances.md`
- `tables/run-117_abc_comparator_test_distances.md`
- `tables/run-117_top_component_distances.md`
- `tables/run-117_observable_vs_raw_epistasis_and_landscape.md`
- `tables/run-117_best_candidate_parameters.md`
- `tables/run-117_stage_failure_modes.md`
- `tables/run-117_output_collection_audit.md`
- `tables/run-117_graphlet_and_module_catalog_summary.md`

## Collection policy

`RUN-117` completed on `lab-slurm` as Slurm job `1911` with inferred exit code
`0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-117/outputs`

This analysis copied the compact JSON/CSV outputs needed to interpret the
experiment. The copied source payload is `2,633,098` bytes and is audited in
`tables/run-117_output_collection_audit.md`. The detailed coefficient tables
were not copied because their aggregate spectral summaries are already present
in `graphlet_summary.csv` and `graphlet_metrics.csv`: `gft_coefficients.csv`
(`26.2 MB`) and `walsh_terms.csv` (`25.8 MB`) remain on the runtime volume.
The PHOT alignment FASTA also remains on the runtime volume because it is not
needed for the stage-level decision.

## Output audit

`RUN-117` evaluated `24` PHOT graphlets, `312` candidate parameter sets across
the four cumulative stages, selected `32` candidates, and produced `32`
equivalence rows against the best available `RUN-098` full-ABC comparator.

The comparator used for the formal decision is the best available `RUN-098`
row in the copied comparator distance table:
`empirical_readout_residual` /
`latent_phi__minus_poly_fitted_current_best`, with test distance `0.2115`.
For context, the best nested full-ABC family row in the same table is
`random_module`, with test distance `0.4826`.

## Main result

No cumulative fix stage passed the preregistered operational-equivalence
decision. `RUN-117` reports `0/32` operationally equivalent rows. No selected
candidate passed the absolute or relative aggregate-distance margin against
the best available full-ABC comparator.

The best `RUN-117` held-out row was the full `A+B+C+D` stage,
`A+B+C+D__local_refinement_round_1__0671141027`, with test distance `0.4122`.
That is a small improvement over the best `RUN-115` bounded candidate reported
in `RES-043` (`0.4269`) and is better than the best nested full-ABC family
context row (`0.4826`), but it remains far from the formal best available
full-ABC comparator (`0.2115`). The stronger equivalence claim therefore
fails.

## Stage interpretation

Stage `A`, the monotone-piecewise observation link alone, preserved the
improved component-guardrail behavior from `RES-043`: all `8/8` selected rows
passed component guardrails. The best held-out row had test distance `0.4235`,
observable magnitude epistasis `0.2289` versus empirical `0.2926`, and
single-mutant SD `0.2030` versus empirical `0.2704`. This is closer on
observable scale than `RUN-115`, but it still has low accessibility
(`0.0164` versus `0.1975`) and raw latent magnitude epistasis remains very
large (`178.6`).

Stages `A+B` and `A+B+C` failed more severely. Their selected best held-out
rows had test distance `2.9206`, zero accessible paths, zero single-mutant SD,
zero observable magnitude epistasis, and zero raw magnitude epistasis. No
selected row passed component guardrails in either stage. The global-readout
and constrained-specific-epistasis settings sampled in these stages therefore
collapsed the selected landscape variation rather than repairing calibration.
Adding the explicit specific-epistasis/additive-scale ratio knob in `C` did
not rescue this failure under the selected fast-fit search.

Stage `A+B+C+D` was the best stage, and the distribution/scale-conditioned
objective fixed the worst numerical conditioning failure seen in `A+B` and
`A+B+C`. All `8/8` selected rows passed component guardrails, and the best
held-out row reached observable magnitude epistasis `0.2426`, single-mutant SD
`0.2172`, and magnitude-epistasis-to-single-mutant-SD ratio `1.1204`, close
to empirical `1.0465`. However, accessibility remained low (`0.0190` versus
`0.1975`), local-energy concentration remained mismatched, and raw latent
variation was effectively absent (`raw_magnitude_epistasis_mean_mean` about
`1.16e-6`). This is not a calibrated latent current-model geometry.

## Decision

`RUN-117` does not validate the cumulative A-through-D calibration fixes as a
full-ABC-equivalent PHOT fitting protocol. The `D` objective materially
improves the numerical behavior of the cumulative model relative to `A+B` and
`A+B+C`, and it slightly improves the best held-out distance over `RES-043`,
but the final stage still fails the aggregate-equivalence margins and retains
a major accessibility/local-energy mismatch.

The experiment sharpens the model diagnosis. The problem is not only that raw
latent epistasis is too large, as in `RUN-115`; one attempted correction can
instead drive raw latent variation to near zero while preserving some
observable epistasis through the readout. For the intended procedural
landscape use, neither extreme is sufficient. The model still needs a latent
geometry that simultaneously supports observable epistasis amplitude,
single-mutant scale, accessibility, local energy, and spectra without relying
on readout-induced artifacts.
