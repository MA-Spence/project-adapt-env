# RES-049: RUN-127 hashed GNK knobs validate forward ruggedness-accessibility map

## Summary

RUN-127 supports the EXP-051 WS-B knob-map claim: on enumerable L=8 coordinates, hashed GNK epistasis weight moved all preregistered ruggedness, spectrum, accessibility, sign-epistasis, and autocorrelation metrics in the expected directions, while the retained core/readout guardrail preserved plateaus and declining beneficial fraction at epi_weight=0. The result is a forward-generator validation, not empirical PHOT overlap or optimizer evidence.

## Generated from

- Analyses: ANA-049

## Relevant hypotheses

- Supports: HYP-005
- Motivates: HYP-004, HYP-009
- Weakens: None
- Refutes: None

## Evidence

- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/source_tables
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/metrics.json
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_output_collection_audit.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_monotonicity_summary.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_epi_weight_metric_summary.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_support_mode_knob_summary.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_joint_regime_summary.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_core_readout_guardrail_summary.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_epi0_hamming_beneficial_summary.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_threshold_epistasis_summary.csv
- analyses/2026-06-18_ANA-049_run-127-hashed-gnk-benchmark-generator-knob-map-review/tables/run-127_support_order_summary.csv
- /tank/data/labproj/runtime/project-adapt-env/runs/RUN-127/outputs

## Interpretation

`RUN-127` completed the EXP-051 WS-B forward-generator validation on enumerable
binary `L = 8` coordinates. The run evaluated three support modes
(`structured_block`, `shuffled_block`, `random`), three seeds, five epistasis
weights, and `96` hashed GNK terms per parameter set. `ANA-049` copied the full
compact output bundle into the analysis record and derived the decision tables
listed above.

The primary result is positive for the knob-map question. All `18/18`
monotonicity checks passed, with minimum absolute Spearman `rho = 0.8821`.
Increasing `epi_weight` moved every preregistered observable-fitness summary in
the expected direction across support modes: local maxima and high-order/sign
epistasis increased, while accessible path fraction and random-walk lag-1
autocorrelation decreased.

The quantitative trend is strong enough to be useful as a forward benchmark
knob. Averaged across support modes and seeds, local maxima increased from
`1.0` at `epi_weight = 0` to `8.67` at `epi_weight = 1.0`; Walsh high-order
fraction increased from `0.0229` to `0.3442`; sign epistasis fraction increased
from `0.0` to `0.4463`; and accessible path fraction fell from `1.0` to
`0.0158`.

The useful joint regime appears at the intermediate setting `epi_weight = 0.6`.
All `9/9` parameter sets at this setting were multi-peaked and still had
nonzero accessible paths. Mean local maxima were `4.89`, mean accessible path
fraction was `0.158`, mean Walsh high-order fraction was `0.193`, and mean sign
epistasis fraction was `0.340`. This directly addresses the failure mode from
`RES-047`, where a sparse Walsh/Potts route could not jointly preserve
accessibility and ruggedness under the forward gate.

The retained core/readout also passed its most important guardrail. At
`epi_weight = 0`, all `9/9` core-only rows had declining beneficial fraction
toward the peak and both low and high readout plateaus. The beneficial-fraction
versus fraction-of-peak Spearman value was `-1.0`. This means the implementation
did not silently collapse to a bare GNK/RMF landscape: the stability/Fuji core
and saturating readout remain active and produce the intended global
epistasis/plateau structure.

The threshold-epistasis result is mixed. Mean absolute pairwise epistasis was
largest in the fold-transition band at `epi_weight = 0.3`, `0.6`, and `1.0`,
which supports threshold-coupled epistasis magnitude. But sign-epistasis and
reciprocal-sign enrichment were not consistently largest at the fold
transition. At `epi_weight = 0.6`, fold-transition sign epistasis (`0.3058`)
was similar to deep-stable (`0.3046`) and below deep-unstable (`0.3227`); at
`epi_weight = 1.0`, fold-transition sign epistasis (`0.4248`) was below both
deep-stable (`0.4541`) and deep-unstable (`0.4384`). The sign-enrichment
subclaim should therefore be treated as not cleanly passed.

## Effect on hypothesis

This supports `HYP-005` in the narrow but important sense tested here: the
hashed GNK `epi_weight` knob is directional over the sampled seed/support
grid, with effect sizes larger than seed-to-seed variation for the registered
summary metrics.

For `HYP-004`, the result is motivation rather than direct support. The
generator shows model-native beneficial-fraction decline, readout plateaus,
accessibility loss, and increasing local maxima as ruggedness rises, but it
does not compare those behaviours to empirical evolutionary walks or a real
protein assay.

For `HYP-009`, the result is also motivation rather than direct support. It
shows that structured, shuffled, and random routed regimes can be generated
under the same global knob map, but it does not run optimizers, PLM-guided
agents, observation-noise regimes, or matched routed benchmark comparisons.

## Limitations

- The run uses enumerable binary `L = 8` coordinates, not the full 20-letter
  protein-scale generator.
- It is model-native. It does not compare against PHOT, GB1, TrpB, or another
  empirical assay, and it does not use a certified measurement operator.
- The support modes are structured/shuffled/random diagnostics, not yet real
  DCA or contact-routed protein supports.
- The threshold-coupled sign-epistasis prediction is mixed: magnitude is
  fold-transition enriched, but sign and reciprocal-sign fractions are not
  consistently enriched there.
- The small nonzero observable high-order fraction at `epi_weight = 0` is a
  readout-induced global-epistasis effect, not GNK texture.

## Downstream use

- Cite `RUN-127` as support that the hashed GNK forward generator has a usable
  `epi_weight` knob for ruggedness/accessibility/spectrum tuning.
- Use `epi_weight = 0.6` as the current small-L joint multi-peak/navigable
  regime when motivating the next benchmark-construction stage.
- Do not cite this result as empirical PHOT overlap, optimizer discrimination,
  or biological routing support.
- Preserve the mixed threshold-sign-epistasis caveat when using this result to
  motivate later `HYP-004` or `HYP-009` experiments.
