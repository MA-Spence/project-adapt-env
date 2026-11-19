# ANA-049: RUN-127 hashed GNK benchmark generator knob-map review

## Purpose

Collect and interpret the `RUN-127` outputs from `EXP-051`, the WS-B
forward-generator test proposed in
`scratch/RES-048_forward_plan_tunable_benchmark_generator.md`.

The analysis asks whether the hashed GNK texture layer, inserted inside the
retained additive stability/Fuji core and saturating readout, produces the
expected direct knob maps on enumerable small coordinates:

- core/readout guardrails at `epi_weight = 0`;
- monotone effects of increasing `epi_weight` on local maxima, accessibility,
  spectrum, sign epistasis, reciprocal sign epistasis, and autocorrelation;
- a joint regime that is both multi-peaked and still navigable;
- evidence for threshold-coupled epistasis near the fold transition.

## Linked experiments/runs

- Experiments: `EXP-051`
- Runs: `RUN-127`

## Notebook record

- Primary notebook: `notebooks/analysis.ipynb`
- Paired text file: `notebooks/analysis.py`
- Derivation script: `notebooks/derive_tables.py`
- Kernel: `python3`

## Produced artifacts

- `source_tables/run-127_summary.json`
- `source_tables/run-127_preflight.json`
- `source_tables/run-127_parameter_grid.csv`
- `source_tables/run-127_graphlet_metrics.csv`
- `source_tables/run-127_knob_summary.csv`
- `source_tables/run-127_monotonicity_checks.csv`
- `source_tables/run-127_support_catalog.csv`
- `source_tables/run-127_hamming_layer_metrics.csv`
- `source_tables/run-127_threshold_epistasis_metrics.csv`
- `source_tables/run-127_core_readout_guardrail_checks.csv`
- `tables/metrics.json`
- `tables/run-127_output_collection_audit.md`
- `tables/run-127_monotonicity_summary.md`
- `tables/run-127_epi_weight_metric_summary.md`
- `tables/run-127_support_mode_knob_summary.md`
- `tables/run-127_joint_regime_summary.md`
- `tables/run-127_core_readout_guardrail_summary.md`
- `tables/run-127_epi0_hamming_beneficial_summary.md`
- `tables/run-127_threshold_epistasis_summary.md`
- `tables/run-127_support_order_summary.md`

## Collection policy

`RUN-127` completed on `lab-slurm` as Slurm job `1935` with inferred scheduler
exit code `0`. The full runtime output remains at:

`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-127/outputs`

The complete compact output bundle was small enough to copy into this analysis:
`721,916` bytes across `10` source files. SHA-256 checksums are recorded in
`tables/run-127_output_collection_audit.md`.

## Output audit

The run evaluated the enumerable binary `L = 8` validation setting specified
for the forward generator. It used three support modes
(`structured_block`, `shuffled_block`, `random`), three seeds
(`501`, `502`, `503`), five epistasis weights (`0.0`, `0.15`, `0.3`, `0.6`,
`1.0`), and `96` hashed GNK terms per parameter set. The output bundle
contains:

- `45` parameter sets;
- `90` graphlet metric rows across latent-margin and observable-fitness
  layers;
- `18` monotonicity checks;
- `4,320` support-term rows;
- `810` hamming-layer rows;
- `954` threshold-epistasis rows;
- `9` core/readout guardrail rows.

## Core/readout guardrails

The core/readout guardrail mostly passed at `epi_weight = 0`. All `9/9`
core-only rows had a declining beneficial fraction toward the peak and both
low and high readout plateaus. The mean Spearman correlation between
beneficial fraction and fraction-of-peak was `-1.0`, and the observable
hamming-layer table shows beneficial fraction moving from `1.0` at the
low-fitness anti-peak to `0.0` at the peak.

This confirms the scratch plan's key implementation guardrail: the generator
was not collapsed to a bare GNK/RMF field. The additive core and saturating
readout still impose the expected peak-proximal depletion of beneficial moves
and plateau structure before any GNK texture is added.

The guardrail is not completely clean. The recorded core-only
`transition_epistasis_enriched_over_deep_regimes` flag passed in `0/9` rows.
In other words, `RUN-127` supports the core/readout plateau and diminishing
returns checks, but not a strong fold-transition enrichment claim at
`epi_weight = 0`.

## Knob-map validation

The primary WS-B result is positive. All `18/18` monotonicity checks passed,
with minimum absolute Spearman `rho = 0.8821` across the preregistered
metric/support-mode checks.

Averaged across support modes and seeds, increasing `epi_weight` moved the
observable-fitness metrics in the expected directions:

- local maxima increased from `1.0` at `epi_weight = 0` to `8.67` at
  `epi_weight = 1.0`;
- accessible path fraction fell from `1.0` to `0.0158`;
- Walsh high-order fraction increased from `0.0229` to `0.3442`;
- sign epistasis fraction increased from `0.0` to `0.4463`;
- reciprocal sign epistasis fraction increased from `0.0` to `0.1286`;
- random-walk lag-1 autocorrelation fell from `0.7341` to `0.4954`.

The small nonzero observable high-order fraction at `epi_weight = 0` comes
from the nonlinear readout, not from GNK texture: the latent-margin layer is
additive at `epi_weight = 0`, while the observable layer has readout-induced
global epistasis.

## Joint multi-peak/navigable regime

The run found a useful joint regime at intermediate `epi_weight`. At
`epi_weight = 0.6`, all `9/9` parameter sets were multi-peaked and still had
nonzero accessible paths. Mean local maxima were `4.89`, mean accessible path
fraction was `0.158`, mean Walsh high-order fraction was `0.193`, and mean
sign epistasis fraction was `0.340`.

At `epi_weight = 1.0`, ruggedness increased further (`8.67` mean local
maxima and `0.344` mean high-order fraction), but accessibility was nearly
collapsed (`0.0158`). Thus the adopted knob is useful because it exposes a
continuous tradeoff rather than a single all-or-nothing ruggedness setting.

## Threshold epistasis

The threshold analysis is mixed. Mean absolute pairwise epistasis was largest
in the fold-transition band at `epi_weight = 0.3`, `0.6`, and `1.0`
(`0.0912`, `0.1638`, and `0.2463`, respectively), which supports a
threshold-coupled epistasis-magnitude interpretation.

However, sign-epistasis and reciprocal-sign enrichment were not consistently
largest at the fold transition. At `epi_weight = 0.6`, fold-transition sign
epistasis (`0.3058`) was similar to deep-stable (`0.3046`) and below
deep-unstable (`0.3227`). At `epi_weight = 1.0`, fold-transition sign
epistasis (`0.4248`) was below both deep-stable (`0.4541`) and deep-unstable
(`0.4384`). This should be recorded as a limitation of the current diagnostic,
not as a clean pass of the scratch plan's threshold-sign-epistasis prediction.

## Interpretation

`RUN-127` validates the first forward-generator step that the PHOT inverse-fit
series lacked: a directly tunable epistatic texture layer whose knob maps are
measured on the generator's own enumerable coordinates. The result supports
`HYP-005` for the hashed GNK `epi_weight` knob in this diagnostic setting.

The result does not yet validate empirical realism. It does not compare PHOT
or a second assay through a certified measurement operator, and it does not
test optimizer or PLM-agent performance. It should therefore motivate later
`HYP-004` and `HYP-009` work rather than be cited as direct support for
biological adaptive-walk realism or benchmark-agent discrimination.

## Decision

Use the hashed GNK composite as the current forward generator candidate for
the next calibrated benchmark-construction stage. The adopted evidence is the
directional knob map and the `epi_weight = 0.6` joint regime. Do not cite
`RUN-127` as empirical PHOT overlap, a full 20-letter protein-scale test, or
evidence that biologically routed supports affect optimizer behaviour.
