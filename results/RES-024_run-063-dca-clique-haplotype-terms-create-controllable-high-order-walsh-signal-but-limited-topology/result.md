# RES-024: RUN-063 DCA clique haplotype terms create controllable high-order Walsh signal but limited topology

## Summary

DCA/contact-graph clique exact-match MSA haplotype terms create controllable high-order Walsh coefficients, including a nonzero fifth-order signal at the upper strength, but lower-order/additive variance remains substantial and sign epistasis, RSE, and multi-maximum topology remain limited.

## Generated from

- Analyses: ANA-024

## Relevant hypotheses

- Supports: HYP-005
- Weakens: HYP-004
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-024_run-063-tem1-dca-clique-haplotype-high-order-epistasis-review/metrics.json`
- `analyses/2026-05-15_ANA-024_run-063-tem1-dca-clique-haplotype-high-order-epistasis-review/tables/run-063_key_metrics.md`
- `analyses/2026-05-15_ANA-024_run-063-tem1-dca-clique-haplotype-high-order-epistasis-review/tables/run-063_validation_tests.md`
- `analyses/2026-05-15_ANA-024_run-063-tem1-dca-clique-haplotype-high-order-epistasis-review/tables/run-063_order_strength_summary.md`
- `analyses/2026-05-15_ANA-024_run-063-tem1-dca-clique-haplotype-high-order-epistasis-review/tables/run-063_upper_strength_topology.md`
- `analyses/2026-05-15_ANA-024_run-063-tem1-dca-clique-haplotype-high-order-epistasis-review/tables/run-063_order5_upper_walsh_variance.md`
- `analyses/2026-05-15_ANA-024_run-063-tem1-dca-clique-haplotype-high-order-epistasis-review/tables/run-063_dca_haplotype_terms_summary.md`
- `results/RES-024_run-063-dca-clique-haplotype-terms-create-controllable-high-order-walsh-signal-but-limited-topology/metrics.json`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/summary.json`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/dca_haplotype_terms.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/hamming_map_summary.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/hamming_walsh_terms.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/hamming_sign_epistasis_pairs.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/order_strength_summary.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063.dvc`
- `experiments/2026-05-15_EXP-024_tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-024_tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/runs/RUN-063.yaml`

## Interpretation

- `RUN-063` provides a clean architecture diagnostic because the old generic epistasis layers are disabled: `n_epistatic_terms = 0`, `epistasis_strength = 0`, `epistasis_covariation_bias = 0`, and `n_higher_order_epistatic_terms = 0`.
- The DCA/contact-clique exact-match haplotype overlay behaves as a direct high-order Walsh control. The run selected all requested terms (`48` total; `12` per expected order from `2` through `5`) and evaluated `240` dense `2^(n+1)` maps with `7200` total Hamming nodes.
- Zero-strength maps have no meaningful target-order Walsh signal: maximum zero-strength target coefficient is `2.8422e-13`, and maximum zero-strength target-order variance fraction is `1.6604e-29`.
- Target coefficients match the expected exact-match scaling. Across all tested orders and strengths, the maximum absolute error from `strength / 2^order` is `2.8422e-13`.
- The preregistered fifth-order signal is present at the upper setting. Order-5 maps at strength `3000` have mean fifth-order target coefficient `93.7500` and mean fifth-order variance fraction `0.0294`.
- The high-order term is measurable but not variance-dominant. In the upper order-5 maps, mean Walsh variance fractions are `0.2344` at order 1, `0.2945` at order 2, `0.2945` at order 3, `0.1472` at order 4, and `0.0294` at order 5.
- Sign epistasis and RSE are present but weak side effects. In upper order-5 maps, mean sign epistasis fraction is `0.0611` and mean reciprocal sign epistasis fraction is `0.0194`.
- Local topology remains only weakly rugged. Upper order-5 maps average `1.5` local maxima, have at most `2` local maxima, and `6/12` order-5 maps have multiple local maxima.

## Effect on hypothesis

- `HYP-005` is supported for the DCA clique/haplotype exact-match architecture. Coefficient strength is a highly directional and reproducible knob for target high-order Walsh coefficients, including controllable fifth-order signal.
- `HYP-004` is weakened as a full rugged/evolutionary-landscape claim. The run shows that high-order terms can be inserted by construction, but it does not show strong reciprocal sign epistasis, strong multi-peak topology, adaptive-walk accessibility changes, or peak-proximal evolutionary behavior.

## Limitations

- This is a targeted local dense-Hamming-map diagnostic, not a whole-landscape adaptive-walk or random-walk experiment.
- The exact-match haplotype overlay creates high-order Walsh signal mathematically; this does not by itself establish that the induced topology is biologically realistic.
- The indicator form spreads variance across lower Walsh orders, so high-order coefficients increase without becoming the dominant variance component.
- DCA/covariation is used as a contact-like prior for clique selection, not as a calibrated Potts or structural energy model.
- RSE and local maxima were recorded as side effects, not optimized directly in this design.

## Downstream use

- Use `RES-024` as evidence that DCA clique exact-match haplotype terms are a working model-family lever for controllable high-order Walsh signal.
- Do not cite `RUN-063` as evidence that high-order terms alone solve the broader biological ruggedness or reciprocal-sign-epistasis problem.
- Compare future RSE/topology-focused designs against the `RUN-063` high-order coefficient baseline rather than against the earlier pair-only `RUN-061` baseline.
