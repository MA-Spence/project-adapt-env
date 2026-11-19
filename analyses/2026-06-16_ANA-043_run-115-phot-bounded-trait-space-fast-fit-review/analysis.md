# ANA-043: RUN-115 PHOT bounded trait-space fast-fit review

## Purpose

Review RUN-115 compact outputs to determine whether bounded trait-space epistasis in the current AdaptEnv model, fitted by a numerical fast fitter, produces PHOT graphlet, accessibility, epistasis, diminishing-return, spectral, and local-energy summaries comparable to empirical PHOT and the RUN-098 full-ABC comparator.

## Linked experiments/runs

- Experiments: EXP-044
- Runs: RUN-115

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-115_summary.json`
- `source_tables/run-115_preflight.json`
- `source_tables/run-115_remote_output_manifest.tsv`
- `source_tables/run-115_bounded_gate_candidates.csv`
- `source_tables/run-115_selected_bounded_gate_candidates.csv`
- `source_tables/run-115_bounded_gate_distances.csv`
- `source_tables/run-115_bounded_gate_equivalence.csv`
- `source_tables/run-115_graphlet_summary.csv`
- `source_tables/run-115_graphlet_metrics.csv`
- `source_tables/run-115_abc_comparator_selected_family_distances.csv`
- `source_tables/run-115_graphlet_catalog.csv`
- `source_tables/run-115_module_catalog.csv`
- `tables/metrics.json`
- `tables/run-115_output_collection_audit.md`
- `tables/run-115_equivalence_decision.md`
- `tables/run-115_bounded_gate_test_distances.md`
- `tables/run-115_abc_comparator_test_distances.md`
- `tables/run-115_top_component_distances.md`
- `tables/run-115_observable_vs_raw_epistasis_and_landscape.md`
- `tables/run-115_best_candidate_parameters.md`
- `tables/run-115_graphlet_and_module_catalog_summary.md`

## Collection policy

`RUN-115` completed on `lab-slurm` as Slurm job `1909` with inferred exit code
`0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs`

This analysis copied the compact JSON/CSV outputs needed to interpret the
experiment. The copied payload is `940,677` bytes and is audited in
`tables/run-115_output_collection_audit.md`. The detailed coefficient tables
were not copied because their aggregate spectral summaries are already present
in `graphlet_summary.csv` and `graphlet_metrics.csv`: `gft_coefficients.csv`
(`6.4 MB`) and `walsh_terms.csv` (`6.2 MB`) remain on the runtime volume.

## Output audit

`RUN-115` evaluated `24` PHOT graphlets, `102` bounded current-model candidate
parameter sets, selected `12` candidates, and produced `12` equivalence rows
against the best available `RUN-098` full-ABC comparator.

The comparator used for the formal decision is the best available RUN-098 row
in the copied comparator distance table:
`empirical_readout_residual` /
`latent_phi__minus_poly_fitted_current_best`, with test distance `0.2115`.
For context, the best nested full-ABC family row in the same table is
`random_module` with test distance `0.4826`.

## Main result

No selected bounded trait-space candidate passed the preregistered operational
equivalence decision. `RUN-115` reports `0/12` operationally equivalent rows.

The failure is not the same failure mode as `RUN-112`. All selected bounded
candidates passed the component-distance guardrail (`12/12` rows), whereas
`RUN-112` failed component guardrails for every row. The remaining failure is
aggregate held-out distance versus the best available full-ABC comparator: all
rows failed both the absolute and relative margins.

The best held-out bounded candidate was
`bounded_trait_space__local_refinement_round_1__4177169950`, with test
distance `0.4269`. That is much closer than the `RUN-112` best calibrated
cheap row (`1.1474`), but it is still worse than the best available full-ABC
context row (`0.2115`). It is slightly better than the best nested full-ABC
family context row (`0.4826`), so the result should not be read as a complete
failure of the bounded route. It is a failure of the stronger equivalence
claim against the best available comparator.

## Landscape diagnostics

The bounded route improved the previous observable epistasis-amplitude failure
but did not fully calibrate the landscape. Empirical PHOT test-graphlet
magnitude epistasis mean is `0.2926`. The best held-out bounded candidate has
observable magnitude epistasis mean `0.1708`; the best calibration-objective
candidate has `0.1698`. This is closer than the `RUN-112` best row
(`0.0156`), but still low by about `0.12` on the empirical scale.

Raw latent epistasis remains badly miscalibrated. The best held-out bounded
candidate has raw magnitude epistasis mean `103.65`, and the best
calibration-objective candidate has raw magnitude epistasis mean `67.84`,
versus empirical raw/observable magnitude `0.2926`. The positive affine
readout therefore hides a large latent-scale mismatch rather than proving that
the underlying current-model trait geometry is empirically calibrated.

Accessibility remains tradeoff-dependent. The best held-out candidate has
accessible-path fraction `0.0120` versus empirical `0.1975`. The best
calibration-objective candidate matches accessibility well (`0.1986`), but
has worse aggregate test distance (`0.4744`) and still undershoots observable
epistasis magnitude.

Diminishing-returns epistasis is not the main failure in this run. The best
held-out candidate has diminishing-returns fraction `0.3997` versus empirical
`0.3872`. Local-energy and spectral summaries are mixed: the best held-out
candidate has edge top-10 energy fraction `0.3476` versus empirical `0.4060`,
site top-1 energy `0.3149` versus `0.4162`, and high-order Walsh fraction
`0.4785` versus `0.3788`.

## Decision

`RUN-115` does not validate the bounded trait-space numerical fast fitter as an
equivalent replacement for the previous full-ABC PHOT result. It does,
however, narrow the model diagnosis: bounding epistasis before the readout can
remove the worst component-guardrail failures and produce a substantially
closer PHOT held-out landscape than the calibrated nested-min cheap protocol.

The remaining gap is model calibration, not merely search cost. The observable
readout rescales a current-model landscape whose raw latent epistasis remains
tens to hundreds of units too large, and the best distance/accessibility
solutions are not the same parameterization.
