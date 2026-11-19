# EXP-023: TEM-1 DCA-specific epistasis dense Hamming graph diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005

## Question

Can DCA-prior, state-specific pairwise epistasis alone create local topological ruggedness, sign epistasis, and constrained adaptive paths in a TEM-1 stability-margin landscape when the old generic low-rank epistasis layer is disabled?

## Pre-experiment prediction

If smooth low-rank pairwise terms were the limiting architecture, then sparse DCA-selected state-specific compensatory and incompatibility terms should produce measurable reciprocal/sign epistasis, local optima, reduced random-walk autocorrelation, and restricted monotonic paths on dense 2^k local Hamming graphs.

## Rationale

RUN-059 showed that increasing generic epistasis magnitude can create pairwise nonadditivity but leaves random walks highly autocorrelated and sign epistasis nearly absent. This experiment isolates a more biologically specific architecture by routing DCA-prior state-specific pair terms into latent stability and disabling the old generic pairwise and higher-order epistasis layers.

## Experimental design

- Build a TEM-1 MMSeqs alignment profile and use its mean-field DCA-like
  coupling scores as the only source of epistatic pair placement.
- Disable the old generic low-rank pairwise epistasis, covariation-biased random
  placement, higher-order exact-match terms, and noise in the AdaptEnv
  configuration.
- Select a concrete local background by realized effective stability margin, then
  choose DCA-coupled residue pairs whose alternative amino-acid states are
  observed/enriched in the alignment and preferentially have deleterious single
  effects on that background.
- Convert the selected `k` mutations into binary variables and enumerate the full
  `2^k` local Hamming graph with landscapy.
- Analyse local optima, accessible shortest paths, greedy walk endpoint, random
  walk autocorrelation, Dirichlet energy, Walsh epistasis, and pairwise
  sign/reciprocal-sign epistasis across the dense graph.

## Inputs

- TEM-1 beta-lactamase reference sequence from the experiment config.
- MMSeqs/ColabFold alignment fetched or restored through the configured cache.
- AdaptEnv alignment-derived stability and DCA/covariation profile.
- Landscapy dense binary Hamming graph construction.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- tem1_mmseqs_alignment.fasta
- alignment_profile_summary.json
- anchor_candidates.csv
- selected_dca_pairs.csv
- dca_pairwise_override.npz
- hamming_variables.csv
- hamming_nodes.csv
- hamming_pair_epistasis.csv
- hamming_walsh_terms.csv
- hamming_graph_metrics.json
- summary.json
- landscapy_hamming_landscape/manifest.json

## Analysis plan

- Create an ANA record after execution.

## Completion criteria

- [ ] Config committed
- [ ] Run script committed
- [ ] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
