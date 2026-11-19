# ANA-039: RUN-098 PHOT nested soft-min ensemble capacity review

## Purpose

Review RUN-098 nested soft-min ensemble ABC outputs, held-out graphlet/local-ruggedness performance, and empirical-module specificity controls.

## Linked experiments/runs

- Experiments: `EXP-039`
- Runs: `RUN-098`
- Source data:
  `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098`
- DVC pointer:
  `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`

## Scientific context

`EXP-039` tests whether the positive capacity signal from the earlier PHOT
coupled-landscape experiments can be extended by adding a second independently
parameterised soft-min/Boltzmann coupling layer. The first coupled layer is
locked to the `RUN-084` rank-1 parameterisation, so the live fit only estimates
a third latent residual landscape and the second coupling layer:
`c(x)=g(a(x),b(x)); f(x)=h(c(x),d(x))`.

The scientific decision is based on calibration-selected parameters evaluated
on held-out graphlets. Test-oracle posterior particles are reported only as a
capacity diagnostic, not as support, because held-out graphlets were not
available to the selector.

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-098_preflight_and_collection_audit.md`
- `tables/run-098_selected_test_distances.md`
- `tables/run-098_localization_components_test.md`
- `tables/run-098_module_locality_test.md`
- `tables/run-098_test_oracle_capacity_diagnostic.md`
- `tables/run-098_hyp008_decision_table.md`

## Run and data audit

`RUN-098` completed on `lab-slurm` as Slurm job `247` with inferred exit code
`0`. The collected output directory contains 30 files and is DVC-tracked and
backed up at
`data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
(`md5 e1e92177585cf2d09ed03bb3495b6ae4.dir`, 206,813,520 bytes).

The preflight reports 24 PHOT graphlets split into 12 calibration and 12
held-out test graphlets. The objective contains 16 metrics: the original nine
SE/RSE, topology, random-walk, WHT, and GFT components plus seven local
Dirichlet-energy localization components.

## Key observations

- The calibration-selected empirical-module nested model does not beat the
  locked first-stage model or the non-nested fitted baselines on held-out
  full-objective distance: empirical-module is `0.543180`, locked first-stage
  is `0.529590`, fitted current is `0.521007`, and smooth core is `0.526819`.
- The best calibration-selected fitted nested family is the random-module
  negative control (`0.482630`), followed by shuffled-module (`0.510961`). The
  empirical-module family loses to no-module (`0.537199`), random-module, and
  shuffled-module, and only beats uncoupled-module (`0.640141`).
- The local-ruggedness component package is closer to empirical PHOT for some
  nested candidates than for fitted current/smooth core (`0.674567` Local7),
  but not in an empirical-module-specific way. Shuffled-module has Local7
  distance `0.419862`, no-module `0.476565`, empirical-module `0.525045`, and
  random-module `0.553193`.
- The empirical-module selected candidate remains diffuse across empirical
  modules. On held-out graphlets it has edge effective module count
  `11.727783/12` with normalized entropy `0.995298`, and plaquette effective
  module count `10.560660/12` with normalized entropy `0.972360`.
- The empirical readout residual diagnostic remains much closer to empirical
  PHOT (`0.211455` full distance; `0.135199` Local7), reinforcing that residual
  structure exists but is not captured specifically by the tested nested
  empirical-module family.
- Test-oracle posterior particles show unused capacity: the best held-out
  no-module particle reaches `0.395321`, and the best empirical-module particle
  reaches `0.409176`. These particles were not selected by calibration and
  therefore do not count as held-out support.

## Interpretation

`RUN-098` is negative for the nested empirical-module implementation of
`HYP-008`. Adding a third latent landscape and a second soft-min coupling layer
does not produce a calibration-selected empirical-module candidate that beats
the locked first-stage model, fitted current, smooth core, or matched
no/random/shuffled controls on the primary held-out objective.

The run preserves only a weak generic-capacity signal: a random-module nested
control can improve the held-out full-objective distance relative to fitted
current/smooth core and locked first stage, and several posterior particles
would have scored better under test-oracle selection. That is not a mechanism
support result. The empirical-module energy remains diffuse, and the best
selected improvement is carried by a negative control.

## Limitations

- PHOT is the only empirical system in this run.
- The decision uses one fixed calibration/test graphlet split.
- Test-oracle posterior particles are not valid support because they use
  held-out graphlets for selection.
- The empirical target is the processed `latent_phi` layer.
- The local ruggedness metrics are graphlet-local summaries, not full
  landscape-wide localization estimates.
