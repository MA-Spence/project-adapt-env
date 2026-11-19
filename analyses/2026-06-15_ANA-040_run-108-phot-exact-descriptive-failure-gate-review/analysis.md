# ANA-040: RUN-108 PHOT exact descriptive failure gate review

## Purpose

Review RUN-108 exact complete-graphlet PHOT descriptive failure tables without full output collection, identify which RUN-039 full-ABC summaries fail, and decide whether the descriptive gate justifies later ABC target construction.

## Linked experiments/runs

- Experiments: EXP-041
- Runs: RUN-108

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `source_tables/run-108_summary.json`
- `source_tables/run-108_abc_descriptive_gate.csv`
- `source_tables/run-108_exact_model_distances.csv`
- `source_tables/run-108_exact_model_failures.csv`
- `source_tables/run-108_graphlet_catalog.csv`
- `tables/metrics.json`
- `tables/run-108_exact_model_distances.md`
- `tables/run-108_abc_descriptive_gate.md`
- `tables/run-108_component_failure_summary.md`
- `tables/run-108_top_component_failures.md`
- `tables/run-108_graphlet_catalog_summary.md`

## Collection policy

`RUN-108` completed on `lab-slurm`, but this analysis deliberately did not run
full `labproj collect`. The remote output directory was inspected directly and
only compact result tables were copied into this analysis record:

| remote output | size | copied |
|:--|--:|:--|
| `summary.json` | 6.4 KB | yes |
| `abc_descriptive_gate.csv` | 1.4 KB | yes |
| `exact_model_distances.csv` | 1.2 KB | yes |
| `exact_model_failures.csv` | 13 KB | yes |
| `graphlet_catalog.csv` | 4.1 KB | yes |
| `exact_graphlet_metrics.csv` | 268 KB | no |
| `exact_graphlet_values.csv` | 1.2 MB | no |

This preserves the relevant result and avoids moving the larger per-node and
per-graphlet value tables over the network. The authoritative full output
remains on the runtime volume under
`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-108/outputs`.

## Output audit

`RUN-108` evaluated `48` complete PHOT graphlets: `16` graphlets each at
dimensions `4`, `5`, and `6`. The empirical layer was the `RUN-039` PHOT
MAVE-NN `latent_phi` landscape, built from `167,529` prepared assay rows and
`167,530` mutation-set entries including the wild type.

The compared fixed model branches were:

- `run039_smc_abc_best_stability_margin`
- `run039_smc_abc_posterior_mean_stability_margin`
- `run039_baseline_shared_stability_readout`
- `run039_predictive_richpair_shared_stability_readout`

The analysis is therefore an exact descriptive benchmark of existing `RUN-039`
branches. It is not a new ABC fit.

## Main result

The exact descriptive gate does not show broad failure across the full summary
package. For all four `RUN-039` branches, the large failures are concentrated
in local Dirichlet-energy concentration summaries:

- `edge_energy_top10_fraction`
- `edge_energy_effective_fraction`
- `site_energy_effective_fraction`
- `site_energy_top1_fraction`

The `RUN-039` SMC-ABC best particle is still the best of the four fixed
branches by mean absolute scaled error, but it fails the same four local-energy
components:

| candidate | mean abs scaled error | max abs scaled error | failed components >= 2 |
|:--|--:|--:|--:|
| SMC-ABC best | `0.5307` | `2.3369` | `4/19` |
| SMC-ABC posterior mean | `0.5394` | `2.3667` | `4/19` |
| predictive rich-pair stability readout | `0.5507` | `2.4097` | `4/19` |
| baseline stability readout | `0.5507` | `2.4097` | `4/19` |

By contrast, the exact gate does not identify severe failure for sign epistasis,
reciprocal sign epistasis, local maxima/minima, basin count, WHT high-order
fraction, GFT low/mid/high energy fractions, random-walk lag-1
autocorrelation, node-incident effective fraction, or edge-energy Moran I.

Accessibility is directionally weak but does not cross the preregistered severe
threshold: all fixed `RUN-039` branches have scaled accessibility error
`-1.183`, below the `2.0` gate.

## Interpretation

`RUN-108` supports the methodological rule behind `EXP-041`: ABC should not be
used as the first diagnostic. The exact analysis gives a much narrower target
than another broad ABC search would have used. The current `RUN-039` model is
not globally failing the complete-graphlet package; it is specifically too
concentrated in edge/site Dirichlet energy relative to empirical PHOT.

This weakens any claim that the `RUN-039` full validation-objective fit already
matches local ruggedness structure. It also avoids overstating the failure:
the exact spectral, sign-epistasis, local-optima, basin, and random-walk
summaries are not the components selected by the descriptive gate.

The result is most relevant to `HYP-004` and `HYP-009`. For `HYP-004`, it shows
that the current PHOT model can match many exact dense-graphlet summaries but
still misses local energy localization. For `HYP-009`, it sharpens the
benchmark-facing target: a future matched benchmark should control edge/site
energy concentration before comparing optimizer behavior, rather than
re-running ABC over the entire old summary vector.

## Limitations

- The empirical target is processed `latent_phi`, not direct raw assay score.
- The analysis uses complete graphlets selected from frequent observed PHOT
  mutation features; it is exact on those subcubes, not on all possible PHOT
  variant combinations.
- Only fixed `RUN-039` branches were benchmarked. `RUN-091`, `RUN-096`, and
  other later coupled-family fits are contextual evidence, not direct
  candidates in this exact gate.
- The source full outputs were not collected into `data/processed`; only the
  compact tables required for this analysis were copied.

## Decision

The descriptive gate permits a later ABC target only for the four local-energy
concentration summaries listed above. A new ABC run that reuses the entire old
validation or graphlet objective would be less disciplined than the evidence
supports.
