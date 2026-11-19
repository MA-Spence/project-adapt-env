# ANA-047: RUN-123 PHOT sparse Walsh Potts spectrum-calibrated generator review

## Purpose

Review RUN-123 compact outputs to determine whether the spectrum-calibrated sparse Walsh/Potts Generator A clears the H1 forward generative capability gate and should be promoted to sequence-conditioned PHOT equivalence fitting.

## Linked experiments/runs

- Experiments: EXP-048
- Runs: RUN-123

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-123_summary.json`
- `source_tables/run-123_preflight.json`
- `source_tables/run-123_remote_output_manifest.tsv`
- `source_tables/run-123_sparse_walsh_spectrum_calibration.csv`
- `source_tables/run-123_sparse_walsh_forward_probe.csv`
- `source_tables/run-123_bounded_gate_equivalence.csv`
- `source_tables/run-123_abc_comparator_selected_family_distances.csv`
- `source_tables/run-123_graphlet_catalog.csv`
- `source_tables/run-123_module_catalog.csv`
- `tables/metrics.json`
- `tables/run-123_sparse_walsh_spectrum_calibration.md`
- `tables/run-123_forward_probe_metric_gate.md`
- `tables/run-123_forward_probe_candidate_gate.md`
- `tables/run-123_stage_gate_and_fit_decision.md`
- `tables/run-123_comparator_context.md`
- `tables/run-123_graphlet_and_module_catalog_summary.md`
- `tables/run-123_output_collection_audit.md`

## Collection policy

`RUN-123` completed on `lab-slurm` as Slurm job `1930` with inferred scheduler
exit code `0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-123/outputs`

This analysis copied the compact JSON/CSV outputs needed to interpret the
stage gate. The copied source payload is `54,267` bytes across `17` files and
is audited in `tables/run-123_output_collection_audit.md`. The PHOT alignment
FASTA was not copied because it is not needed for the sparse-Walsh gate
decision.

## Output audit

`RUN-123` evaluated `24` PHOT graphlets, split evenly into `12` calibration and
`12` test graphlets across dimensions `4`, `5`, and `6`. It recorded `48`
candidate support modules across empirical-calibration plaquette, random,
shuffled, and uncoupled empirical-module sources.

The run produced `4` sparse-Walsh spectrum-calibration rows and `84` H1
forward-probe rows. It produced no sequence-conditioned candidate evaluations,
no selected candidates, and no operational-equivalence rows. This is not a
collection failure: the `G+H1+FIT` stage was blocked because the preregistered
H1 forward gate did not pass.

## Main result

The spectrum-calibrated sparse Walsh/Potts Generator A did not clear the H1
forward generative capability gate. The summary reports
`sparse_walsh_forward_probe_passed = false`, `0` joint-pass rows across `12`
probe candidates, and therefore `0` sequence-conditioned fast-fit evaluations.

The G-stage empirical spectrum calibration was sharply concentrated. Of orders
`2` through `5`, only order `3` passed the preregistered noise-floor flag. The
normalized global sparse-Walsh order weights were therefore:

- order `2`: `0.0`
- order `3`: `1.0`
- order `4`: `0.0`
- order `5`: `0.0`

This is already an important constraint on the intended Generator A
interpretation. The empirical calibration graphlets did not support a broad
high-order mixture under the implemented noise-floor criterion. The H1 probe
therefore tested a mostly third-order sparse field, not a broad order-2-to-5
Potts/Walsh texture.

## H1 failure diagnosis

Every H1 probe candidate passed exactly `4/7` target metrics and failed the
joint gate. The stable failure metrics were:

- `accessible_path_fraction`: empirical mean `0.2319`, model central-interval
  envelope `0.0` to `0.0715`; the generated landscapes were much less
  accessible than PHOT.
- `magnitude_epistasis_to_single_mutant_sd_ratio`: empirical mean `1.0161`,
  model interval envelope `1.0253` to `1.2692`; the generated landscapes
  overshot the empirical ratio for every probe candidate.
- `walsh_high_order_fraction`: empirical mean `0.3433`, model interval
  envelope `0.4479` to `0.7659`; despite the order-3-only calibrated weights,
  the generated fields retained too much high-order Walsh fraction in the
  reporting metric.

The metrics that did pass for all candidates were local maxima, local sign
epistasis, local reciprocal sign epistasis, and diminishing-returns epistasis.
That means Generator A has some qualitative ruggedness capacity under this
implementation, but it does not jointly match the accessibility, scale-normalized
epistasis, and spectral envelope required by the prospective protocol.

## Comparator context

The copied `RUN-098` comparator table is context only because no sparse-Walsh
FIT candidates were generated. The best available `RUN-098` context row remains
`empirical_readout_residual` with test distance `0.2115`; the best nested
full-ABC family row remains `random_module` with test distance `0.4826`.

No `RUN-123` row can be compared to those values. The correct interpretation is
not that the sparse-Walsh sequence-conditioned fit lost to `RUN-098`; the
sequence-conditioned fit was never reached because the forward gate failed.

## Decision

`RUN-123` does not validate Generator A as the next PHOT repair. It falsifies
the pre-fit condition in the EXP-048 protocol: the spectrum-calibrated sparse
Walsh/Potts field did not pass the forward generative capability probe. Under
the preregistered integration plan, promotion to sequence-conditioned fitting
was appropriately skipped.

The result is narrower than a general rejection of sparse Walsh/Potts
generators. It tests the specific implementation and calibration in EXP-048:
graphlet-aligned support clusters, empirical variance-by-order spectrum
calibration, the implemented noise-floor rule, the H1 grid over `epi_weight`
and high-order multipliers, and the PHOT graphlet summaries used in this
project. It also does not evaluate the FGM/H2 branch, which EXP-048 explicitly
left out of scope.

For `HYP-004`, this further weakens the current PHOT adequacy path: the new
sparse high-order operator does not even pass the model-free forward gate on
PHOT accessibility, epistasis-ratio, and spectral summaries. For `HYP-009`,
the run is a negative boundary on this particular Generator A construction
rather than direct evidence about optimizer or PLM-agent behavior.
