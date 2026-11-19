# RES-023: RUN-061 DCA-specific epistasis creates sign epistasis but not multi-peak ruggedness

## Summary

RUN-061 shows that disabling the old generic epistasis and adding DCA-prior state-specific pair terms creates strong selected pair epistasis and systematic one-sided sign epistasis in a dense TEM-1 local Hamming cube, with lower random-walk autocorrelation than prior smooth TEM-1 random-walk diagnostics. It does not yet produce reciprocal sign epistasis, multiple local maxima, accessible improving paths from the selected anchor to the all-alternative genotype, or a second-order-dominated rugged topology.

## Generated from

- Analyses: ANA-023

## Relevant hypotheses

- Supports: HYP-005
- Weakens: HYP-004
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-023_run-061-tem1-dca-specific-epistasis-hamming-graph-review/metrics.json`
- `analyses/2026-05-15_ANA-023_run-061-tem1-dca-specific-epistasis-hamming-graph-review/tables/run-061_key_metrics.md`
- `analyses/2026-05-15_ANA-023_run-061-tem1-dca-specific-epistasis-hamming-graph-review/tables/run-061_validation_tests.md`
- `analyses/2026-05-15_ANA-023_run-061-tem1-dca-specific-epistasis-hamming-graph-review/tables/run-061_selected_dca_pairs.md`
- `analyses/2026-05-15_ANA-023_run-061-tem1-dca-specific-epistasis-hamming-graph-review/tables/run-061_pair_epistasis_summary.md`
- `analyses/2026-05-15_ANA-023_run-061-tem1-dca-specific-epistasis-hamming-graph-review/tables/run-061_walsh_summary.md`
- `results/RES-023_run-061-dca-specific-epistasis-creates-sign-epistasis-but-not-multi-peak-ruggedness/metrics.json`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/summary.json`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_graph_metrics.json`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/selected_dca_pairs.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_pair_epistasis.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_walsh_terms.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/landscapy_hamming_landscape/manifest.json`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061.dvc`
- `experiments/2026-05-15_EXP-023_tem1-dca-specific-epistasis-hamming-graph-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-023_tem1-dca-specific-epistasis-hamming-graph-diagnostic/runs/RUN-061.yaml`

## Interpretation

- `RUN-061` provides a clean architecture diagnostic because the old generic epistasis layer is disabled: `n_epistatic_terms = 0`, `epistasis_strength = 0`, `epistasis_covariation_bias = 0`, and `n_higher_order_epistatic_terms = 0`.
- The DCA-prior state-specific override successfully creates the intended local interaction structure. The run selected `5` DCA-enriched state pairs (`10` binary variables) and enumerated the complete `2^10` Hamming cube with `1024` nodes and `5120` edges.
- All selected pairs have positive alignment enrichment, with log2 enrichment from `3.8124` to `8.3239`.
- Mean absolute selected-pair epistasis is `57.0791`, and selected-pair sign epistasis fraction is `1.0000` across `1280` pair/background contrasts.
- The sign epistasis is one-sided compensatory epistasis, not reciprocal sign epistasis. Reciprocal sign epistasis fraction is `0.0000`.
- Random-walk lag-1 autocorrelation drops to `0.7795`, which is much lower than the near-`0.999` autocorrelation observed in the broad `RUN-059` TEM-1 random-walk sweep. This supports the interpretation that the DCA-specific local graph is less smooth under the sampled graph diagnostic.
- The topology remains a single-peak funnel rather than a multi-peak rugged landscape. The graph has `1` local maximum, `1` local minimum, `0` greedy steps from the selected anchor, and `0/3628800` accessible shortest paths from the selected anchor to the all-selected-alternative genotype.
- Walsh decomposition confirms that the pair terms are localized but not variance-dominant: `5/45` second-order terms are nonzero, first-order terms explain `0.9911` of variance, and second-order terms explain `0.0089`.

## Effect on hypothesis

- `HYP-005` is supported for the specific-epistasis architecture subcomponent. DCA-prior state-specific pair overrides behave as a meaningful mechanistic lever and move pairwise epistasis, sign epistasis, and random-walk autocorrelation in the expected direction when isolated from the old generic epistasis layer.
- `HYP-004` is weakened as a full rugged-landscape claim. `RUN-061` supports the narrower claim that targeted state-specific interactions can create constraining sign epistasis, but it does not show reciprocal sign epistasis, multiple local optima, or adaptive accessibility through a rugged topology.

## Limitations

- This is a targeted local dense Hamming graph, not a whole-landscape adaptive-walk experiment.
- The all-ones genotype is the state carrying all selected DCA alternatives, not an independently established peak or biological endpoint.
- The selected anchor is already the global maximum inside the enumerated cube, so zero accessible improving paths from the anchor mostly records that endpoint choice.
- The DCA/covariation signal is used as a contact-like prior, not as a calibrated Potts model with inferred family-specific energies.
- Pairwise terms are intentionally sparse and state-specific. This run does not test whether the same mechanism produces ruggedness under a broader distribution of backgrounds or pair strengths.

## Downstream use

- Use `RES-023` as evidence that the new DCA-specific epistasis pathway is functioning and can produce strong one-sided sign epistasis.
- Do not cite `RUN-061` as evidence that the current model produces reciprocal sign epistasis or multi-peak biological ruggedness.
- Do not compare the absolute autocorrelation value directly to `RUN-059` without noting that `RUN-061` is a targeted dense local cube, while `RUN-059` used broad sampled random walks over TEM-1 homologue starts.
