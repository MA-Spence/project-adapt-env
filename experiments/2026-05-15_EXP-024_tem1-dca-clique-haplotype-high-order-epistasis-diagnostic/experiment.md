# EXP-024: TEM-1 DCA clique haplotype high-order epistasis diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-005
- Previous experiment: EXP-023
- Motivating result: RES-023

## Question

Can DCA/contact-graph cliques with observed MSA haplotype exact-match terms create controllable third-, fourth-, and fifth-order Walsh-Hadamard components in local TEM-1 Adapt-Env Hamming maps when the old generic epistasis layer is disabled?

## Pre-experiment prediction

Increasing the DCA haplotype exact-match coefficient should increase the target-order Walsh coefficient and target-order variance fraction. The upper strength setting should produce a clear fifth-order Walsh term in order-5 maps.

Sign epistasis and reciprocal sign epistasis are recorded in every combinatorial map, but reciprocal sign epistasis is not the primary success criterion here. The immediate purpose is to solve the missing higher-order Walsh signal identified in RES-023; RSE will be treated as a separate follow-up.

## Rationale

RES-023 showed that DCA-prior pair terms create strong pairwise/sign epistasis, but the dense TEM-1 Hamming cube remained dominated by first-order Walsh terms and had no higher-order signal. That outcome is expected from a linear stability-margin readout plus sparse pair terms.

This experiment tests a structurally different mechanism aligned with protein-family covariation literature: choose DCA/contact-graph cliques, select observed MSA haplotypes at those clique positions, and add sparse exact-match high-order terms on those haplotypes. An exact-match order-n haplotype term should induce a measurable order-n Walsh component on a dense local binary grid.

## Experimental design

- Use the reviewed TEM-1 beta-lactamase sequence `P62593` and the same MMSeqs/ColabFold alignment source used by EXP-023.
- For each seed, build an Adapt-Env base landscape with the old generic pairwise, higher-order, and noise layers disabled.
- Select an anchor sequence by realized effective stability margin, targeting the folding/viability cliff rather than a nominal config value.
- Build a DCA graph from top alignment couplings.
- For expected orders `2`, `3`, `4`, and `5`, identify DCA cliques of that size and choose observed all-alternative MSA haplotypes relative to the anchor.
- Add exact-match haplotype terms with alternating bonus/penalty signs and strengths `0`, `250`, `750`, `1500`, and `3000` in final fitness units.
- For each order-n term, construct a dense local binary Hamming grid with `n + 1` variables: the `n` clique haplotype states plus one additional context mutation.
- Evaluate every node in each map as Adapt-Env base fitness plus the active DCA haplotype overlay.

## Primary readouts

- Target-order Walsh coefficient and target-order variance fraction by order and strength.
- Fifth-order Walsh coefficient and variance fraction in order-5 maps at the upper strength setting.
- Lift in target-order Walsh coefficient relative to the `strength = 0` base map.
- Number of nonzero Walsh terms by order.
- Local maxima count per map.
- Pairwise sign epistasis and reciprocal sign epistasis counts across all pairs/backgrounds in each dense map.

## Inputs

- `config.yaml`
- `scripts/tem1_dca_clique_haplotype_epistasis.py`
- `external/Adapt-Env`
- TEM-1 MMSeqs/ColabFold alignment fetched or cached through the runner.

## Expected outputs

- `tem1_mmseqs_alignment.fasta`
- `alignment_profile_summary.json`
- `anchor_summary.csv`
- `anchor_candidates.csv`
- `dca_haplotype_terms.csv`
- `dca_clique_rejections.csv`
- `hamming_map_specs.csv`
- `hamming_map_summary.csv`
- `hamming_nodes.csv`
- `hamming_walsh_terms.csv`
- `hamming_sign_epistasis_pairs.csv`
- `order_strength_summary.csv`
- `summary.json`

## Analysis plan

After completion, create an ANA record that checks:

- whether target-order Walsh terms increase monotonically with haplotype strength;
- whether fifth-order maps have substantial fifth-order coefficients at the upper strength;
- whether high-order variance is still swamped by first-order stability effects;
- whether sign epistasis and RSE counts move incidentally with high-order strength.

## Completion criteria

- [x] Config created
- [x] Run script created
- [x] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
