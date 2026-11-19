# RES-040: RUN-108 exact PHOT gate isolates local-energy concentration failures

## Summary

RUN-108 exact PHOT complete-graphlet analysis shows the RUN-039 full ABC fit does not broadly fail the graphlet package; severe failures are concentrated in four local Dirichlet-energy concentration summaries, so later ABC is justified only for those summaries rather than the full old objective.

## Generated from

- Analyses: ANA-040

## Relevant hypotheses

- Supports: None
- Weakens: HYP-004
- Motivates: HYP-009
- Refutes: None

## Evidence

- analyses/2026-06-15_ANA-040_run-108-phot-exact-descriptive-failure-gate-review/tables/metrics.json
- analyses/2026-06-15_ANA-040_run-108-phot-exact-descriptive-failure-gate-review/tables/run-108_abc_descriptive_gate.md
- analyses/2026-06-15_ANA-040_run-108-phot-exact-descriptive-failure-gate-review/tables/run-108_component_failure_summary.md
- analyses/2026-06-15_ANA-040_run-108-phot-exact-descriptive-failure-gate-review/tables/run-108_exact_model_distances.md
- analyses/2026-06-15_ANA-040_run-108-phot-exact-descriptive-failure-gate-review/tables/run-108_graphlet_catalog_summary.md
- analyses/2026-06-15_ANA-040_run-108-phot-exact-descriptive-failure-gate-review/tables/run-108_top_component_failures.md

## Interpretation

`RUN-108` completed the `EXP-041` gate exactly as intended: it evaluated
complete observed PHOT graphlets without running another ABC search. The
analysis copied only compact tables into `ANA-040`; the larger exact per-node
and per-graphlet outputs remain on the runtime volume at
`/tank/data/labproj/runtime/project-adapt-env/runs/RUN-108/outputs`.

The empirical layer was the `RUN-039` PHOT MAVE-NN `latent_phi` landscape. The
run evaluated `48` complete graphlets across dimensions `4`, `5`, and `6`, from
`167,529` prepared assay rows and `167,530` mutation-set entries including wild
type.

The fixed `RUN-039` SMC-ABC best particle was the closest candidate by mean
absolute scaled error (`0.5307`), but it still crossed the severe descriptive
gate for `4/19` summary components. The same four components failed for all
four compared branches:

| component | direction | mean absolute scaled error |
|:--|:--|--:|
| `edge_energy_top10_fraction` | model too concentrated in top edge-energy sites | `2.3807` |
| `edge_energy_effective_fraction` | model effective fraction too low | `2.2622` |
| `site_energy_effective_fraction` | model effective fraction too low | `2.2622` |
| `site_energy_top1_fraction` | model too concentrated in top site-energy site | `2.2270` |

The descriptive gate did not identify severe failure for sign epistasis,
reciprocal sign epistasis, local maxima/minima, basin count, WHT high-order
fraction, GFT low/mid/high energy fractions, random-walk lag-1 autocorrelation,
node-incident effective fraction, or edge-energy Moran I. Accessibility was
directionally weak but below the preregistered severe threshold
(`-1.183` scaled error versus an absolute threshold of `2.0`).

## Effect on hypothesis

This weakens the strong adequacy form of `HYP-004`: the current PHOT full-ABC
fit is not a complete local-ruggedness match, even though many exact
graphlet-level summaries are already close. The failure is specific rather
than global, so the result does not justify discarding the current model family
or re-running a broad ABC objective.

For `HYP-009`, the result is methodological evidence rather than a direct
optimizer test. It motivates using edge/site Dirichlet-energy concentration as
a matched benchmark constraint before comparing optimizer classes. It does not
yet show that routing identity affects optimizer performance.

## Limitations

- The empirical target is processed `latent_phi`, not direct raw assay score.
- The graphlets are complete observed PHOT subcubes selected from frequent
  mutation features; the result is exact on those subcubes, not on all possible
  PHOT variant combinations.
- Only fixed `RUN-039` branches were benchmarked. Later coupled-family fits are
  contextual evidence, not direct candidates in this gate.
- The full exact output tables were intentionally not collected locally; this
  record relies on compact copied summaries plus the remote runtime output path.

## Downstream use

Any later PHOT ABC refit should target only the four failed local-energy
concentration summaries unless a separate exact diagnostic identifies
additional failures. Reusing the entire old validation or graphlet summary
vector would be broader than the evidence supports.
