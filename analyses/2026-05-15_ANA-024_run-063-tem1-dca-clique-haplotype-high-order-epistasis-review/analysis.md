# ANA-024: RUN-063 TEM-1 DCA clique haplotype high-order epistasis review

## Purpose

Review whether DCA/contact-graph clique exact-match haplotype terms create controllable third-, fourth-, and fifth-order Walsh components in dense local TEM-1 Hamming maps, and record sign epistasis/RSE side effects.

## Linked experiments/runs

- Experiments: EXP-024
- Runs: RUN-063

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-063_key_metrics.md`
- `tables/run-063_validation_tests.md`
- `tables/run-063_order_strength_summary.md`
- `tables/run-063_upper_strength_topology.md`
- `tables/run-063_order5_upper_walsh_variance.md`
- `tables/run-063_dca_haplotype_terms_summary.md`

## Inputs

- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/summary.json`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/alignment_profile_summary.json`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/dca_haplotype_terms.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/hamming_map_specs.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/hamming_map_summary.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/hamming_walsh_terms.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/hamming_sign_epistasis_pairs.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063/order_strength_summary.csv`
- `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063.dvc`
- `experiments/2026-05-15_EXP-024_tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-024_tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/runs/RUN-063.yaml`
- `results/RES-023_run-061-dca-specific-epistasis-creates-sign-epistasis-but-not-multi-peak-ruggedness/result.md`

## Analysis performed

- Reconciled `RUN-063` as the completed live `EXP-024` execution on `lab-slurm` with scheduler job `223` and scheduler exit code `0`.
- Tracked the collected durable output directory with DVC at `data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063.dvc`.
- Confirmed that the run used the reviewed TEM-1 beta-lactamase reference with a 768-sequence MMSeqs/ColabFold alignment profile and `604.1116` effective sequences.
- Confirmed the isolation condition: old generic epistasis layers remained disabled (`n_epistatic_terms = 0`, `epistasis_strength = 0`, `epistasis_covariation_bias = 0`, `n_higher_order_epistatic_terms = 0`).
- Reviewed the DCA/contact clique exact-match haplotype terms, dense `2^(n+1)` Hamming map summaries, Walsh decompositions, sign-epistasis counts, reciprocal sign-epistasis counts, and local maxima counts.
- Recomputed compact review metrics and Markdown tables from the exported JSON/CSV files using the paired text analysis script at `notebooks/analysis.py`.

## Main observations

- `RUN-063` selected all requested DCA haplotype terms: `48` terms total, with `12` terms each for expected orders `2`, `3`, `4`, and `5`.
- The terms are supported by observed alignment haplotypes, with log2 enrichment from `0.6012` to `13.6527` and haplotype frequencies from `0.001655` to `0.003344`.
- The run evaluated `240` dense maps and `7200` Hamming nodes across three seeds, four expected orders, four maps per order, and five coefficient strengths.
- Zero-strength controls show no meaningful target-order signal: the maximum zero-strength target-order Walsh coefficient is `2.8422e-13`, and the maximum zero-strength target-order variance fraction is `1.6604e-29`.
- Target Walsh coefficients follow the exact-match scaling expected from the construction. Across all rows, the maximum absolute error from `strength / 2^order` is `2.8422e-13`.
- Target Walsh coefficients increase monotonically with strength for all tested orders.
- Target-order variance fraction is mostly but not perfectly monotonic. Orders `3`, `4`, and `5` are non-decreasing, while order `2` saturates and drops slightly from `0.333339` at strength `1500` to `0.333158` at strength `3000`.
- The preregistered upper fifth-order signal is present: order-5 maps at strength `3000` have mean target fifth-order coefficient `93.7500` and mean fifth-order variance fraction `0.0294`.
- Fifth-order variance is not dominant. In the upper order-5 maps, mean Walsh variance fractions are `0.2344` at order 1, `0.2945` at order 2, `0.2945` at order 3, `0.1472` at order 4, and `0.0294` at order 5.
- Sign epistasis and reciprocal sign epistasis remain side effects rather than the primary result. In the upper order-5 maps, mean sign epistasis fraction is `0.0611` and mean RSE fraction is `0.0194`.
- Topological ruggedness remains limited. At upper strength, order-5 maps average `1.5` local maxima, have at most `2` local maxima, and only `6/12` order-5 maps have multiple local maxima.

## Interpretation

RUN-063 answers the immediate high-order architecture question raised after `RES-023`: sparse exact-match DCA clique/haplotype terms are a direct and controllable way to introduce higher-order Walsh components into Adapt-Env local maps. The construction behaves almost exactly as expected mathematically: an order-n exact-match term contributes a target Walsh coefficient of `strength / 2^n`, and the zero-strength controls have no high-order signal.

The result should not be over-interpreted as solving biological ruggedness. The exact-match term creates measurable higher-order coefficients, including fifth-order coefficients at the upper strength, but its indicator structure also spreads variance across lower Walsh orders. In the order-5 upper maps, first-, second-, and third-order variance together still dominate the fifth-order component. The local maps show some extra maxima and weak RSE, but not a strongly rugged, multiply peaked landscape.

## Result records created

- `RES-024`

## Hypothesis updates

- `HYP-005` is supported for this model-family lever: DCA clique exact-match haplotype coefficient strength is a highly directional control on target high-order Walsh coefficients.
- `HYP-004` is weakened as a full biological ruggedness claim. The model family can now create high-order terms by construction, but this run does not show strong reciprocal sign epistasis, strong multi-peak topology, or adaptive-walk accessibility changes.
