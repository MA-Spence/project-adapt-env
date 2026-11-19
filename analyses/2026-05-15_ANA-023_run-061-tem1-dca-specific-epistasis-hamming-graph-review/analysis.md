# ANA-023: RUN-061 TEM-1 DCA-specific epistasis Hamming graph review

## Purpose

Review whether disabling old epistasis and using DCA-prior state-specific pair terms creates sign epistasis, lower random-walk autocorrelation, constrained accessibility, and topological ruggedness in the dense TEM-1 local Hamming cube.

## Linked experiments/runs

- Experiments: EXP-023
- Runs: RUN-061

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-061_key_metrics.md`
- `tables/run-061_validation_tests.md`
- `tables/run-061_selected_dca_pairs.md`
- `tables/run-061_pair_epistasis_summary.md`
- `tables/run-061_walsh_summary.md`

## Inputs

- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/summary.json`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/alignment_profile_summary.json`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/selected_dca_pairs.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_variables.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_nodes.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_pair_epistasis.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_walsh_terms.csv`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/hamming_graph_metrics.json`
- `data/processed/tem1-dca-specific-epistasis-hamming-graph-diagnostic/RUN-061/landscapy_hamming_landscape/manifest.json`
- `experiments/2026-05-15_EXP-023_tem1-dca-specific-epistasis-hamming-graph-diagnostic/config.yaml`
- `experiments/2026-05-15_EXP-023_tem1-dca-specific-epistasis-hamming-graph-diagnostic/runs/RUN-061.yaml`
- `results/RES-022_run-059-strong-tem1-epistasis-sweep-creates-pairwise-epistasis-but-leaves-random-walks-smooth/result.md`

## Analysis performed

- Confirmed from `RUN-061.yaml` and `labproj status` that `RUN-061` completed successfully as the live `EXP-023` execution on `lab-slurm` with scheduler job `209` and inferred exit code `0`.
- Checked that the run used the reviewed TEM-1 beta-lactamase reference `P62593`, a 768-sequence MMSeqs/ColabFold alignment profile with `604.1116` effective sequences, and a selected MMSeqs homologue anchor.
- Confirmed the isolation condition: old generic low-rank pairwise epistasis, covariation-biased random placement, higher-order exact-match terms, and legacy epistasis strength were all disabled.
- Reviewed the selected DCA-prior state pairs, their alignment enrichment, and their single-mutant effects.
- Recomputed compact review tables from the exported JSON/CSV files using only the Python standard library in `notebooks/analysis.py`.
- Interpreted the complete `2^10` landscapy Hamming graph using local pair epistasis, Walsh decomposition, local extrema, random-walk autocorrelation, and accessible shortest paths from the selected anchor to the all-selected-alternative genotype.

## Main observations

- The run enumerated the complete local binary Hamming cube: `1024` nodes, `5120` edges, and `10` selected binary mutation variables arranged as `5` DCA-prior pairs.
- The DCA selector found `5/5` required pairs. All selected pairs had positive joint enrichment in the TEM-1 alignment, with log2 enrichment from `3.8124` to `8.3239`.
- The new architecture clearly creates the intended state-specific pairwise epistasis. Mean absolute selected-pair epistasis is `57.0791` across `1280` pair/background contrasts.
- Sign epistasis is systematic in the selected pairs: sign-epistasis fraction is `1.0000`. Reciprocal sign epistasis remains absent (`0.0000`).
- Random-walk lag-1 autocorrelation in this local dense cube is `0.7795`, much lower than the near-`0.999` random-walk autocorrelation observed in `RUN-059`. This is evidence that the DCA-specific local cube is less smooth under the sampled graph diagnostic.
- The graph is not multi-peak rugged. It has exactly `1` local maximum and `1` local minimum; the anchor is the global maximum, greedy ascent from the anchor takes `0` steps, and the accessible shortest-path fraction from the anchor to the all-selected-alternative genotype is `0.0000`.
- Walsh analysis shows that the selected pair terms are exactly localized: `5/45` second-order terms are nonzero, matching the five DCA pairs. However, first-order terms still explain `0.9911` of variance and second-order terms explain only `0.0089`.

## Interpretation

RUN-061 answers the immediate architecture question from `RES-022`: the absence of sign epistasis in the previous broad TEM-1 sweeps was at least partly architectural. Replacing smooth low-rank/random epistasis with sparse state-specific DCA-prior pair terms can produce strong pairwise epistasis and systematic one-sided sign epistasis when the dense graph is built around the targeted local states.

The result should not be over-read as full biological ruggedness. The selected cube is still dominated by large additive single-mutant effects, has only one local maximum, and shows no reciprocal sign epistasis. The all-selected-alternative endpoint is not a peak; it is a lower-fitness genotype carrying all selected DCA alternatives. Therefore, zero accessible shortest paths from the anchor are evidence of a high-fitness anchor and constrained local descent, not evidence of realistic adaptive path accessibility.

## Result records created

- `RES-023`

## Hypothesis updates

- `HYP-005` is supported for the specific-epistasis architecture subcomponent: the DCA-prior override is a meaningful mechanistic lever that changes pairwise and sign-epistasis metrics when the old epistasis layer is disabled.
- `HYP-004` is weakened or, more precisely, only partially supported. The result supports the possibility of constraining sign epistasis in a targeted local graph, but does not show the richer biological landscape pattern of reciprocal sign epistasis, multiple local optima, or accessible adaptive paths through a rugged topology.
