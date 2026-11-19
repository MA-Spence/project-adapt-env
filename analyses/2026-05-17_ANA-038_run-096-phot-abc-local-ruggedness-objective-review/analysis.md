# ANA-038: RUN-096 PHOT ABC local ruggedness objective review

## Purpose

Review RUN-096 SMC-ABC results after adding local Dirichlet-energy localization statistics to the PHOT calibration objective, testing held-out fit and empirical-module specificity.

## Linked experiments/runs

- Experiments: `EXP-038`
- Runs: `RUN-096`
- Source data:
  `data/processed/phot-abc-local-ruggedness-objective/RUN-096`
- DVC pointer:
  `data/processed/phot-abc-local-ruggedness-objective/RUN-096.dvc`

## Scientific context

`EXP-038` is the direct follow-up to `RES-036` and `RES-037`. `RUN-091`
showed that ABC fitting can improve PHOT held-out graphlet distance but failed
random/uncoupled specificity. `RUN-093` then showed that PHOT has residual
localized ruggedness relative to the fitted current/smooth core, but the
localization was not captured specifically by empirical calibration modules.

`RUN-096` therefore keeps the same PHOT base, graphlet split, and candidate
families as `RUN-091`, but adds local Dirichlet-energy diffuseness statistics
to the calibration objective: adjacent-edge Moran I, edge top-10/effective
fractions, site top-1/effective fractions, and plaquette top-1/effective
fractions. Held-out graphlets remain unused during fitting and are the
scientific comparison.

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-096_preflight_and_collection_audit.md`
- `tables/run-096_selected_test_distances.md`
- `tables/run-096_localization_components_test.md`
- `tables/run-096_hyp008_decision_table.md`

## Run and data audit

`RUN-096` completed with Slurm job `246` and inferred exit code `0`. The
collected output directory contains the planned ABC outputs plus SMC
checkpoints and is DVC-tracked/backed up at
`data/processed/phot-abc-local-ruggedness-objective/RUN-096.dvc`
(`md5 4cd3edc3c50e33ecde7c5c7e32e9ec2a.dir`, 29 files, 115,223,180 bytes).

The preflight reports 24 PHOT graphlets split into 12 calibration and 12
held-out test graphlets. The objective contains 16 metrics: the original nine
SE/RSE, topology, random-walk, WHT, and GFT components plus seven local
Dirichlet-energy localization components.

## Key observations

- The added local objective terms do move fitted-family candidates closer to
  empirical PHOT than fitted current/smooth core on held-out graphlets. Local7
  distance is `0.289976` for uncoupled-module, `0.326809` for no-module, and
  `0.462290` for empirical-module, compared with `0.674567` for fitted
  current and smooth core.
- The calibration-selected held-out full-objective winner among fitted
  families is no-module (`0.277056`), not empirical-module (`0.365272`).
  No-module beats fitted current (`0.521007`), smooth core (`0.526819`),
  empirical-module, random-module (`0.366768`), shuffled-module (`0.379458`),
  and uncoupled-module (`0.291161`).
- The empirical-module family improves over fitted current and smooth core
  (`0.365272` versus `0.521007` and `0.526819`), and narrowly beats random and
  shuffled controls. It still loses to no-module by `0.088216` and uncoupled by
  `0.074111`, which fails the mechanism-specific support rule.
- The empirical readout residual remains the strongest overall selected
  comparator (`0.211455` full distance; `0.135199` Local7), reinforcing that
  the PHOT residual has structure that the tested procedural family is not
  specifically explaining.
- The selected empirical-module energy remains diffuse across modules: held-out
  edge effective module count is `11.010808/12` with normalized entropy
  `0.982059`, and plaquette effective module count is `9.341583/12` with
  normalized entropy `0.944129`.
- A test-oracle empirical-module posterior particle has lower held-out distance
  (`0.259136`) than the selected no-module particle, but it was not selected by
  calibration. This is capacity evidence only, not held-out support.

## Interpretation

`RUN-096` confirms that putting local ruggedness diffuseness statistics into
the ABC objective can select candidates whose held-out localization components
are closer to empirical PHOT than the fitted current/smooth core. That is a
real improvement over the baseline landscape under the expanded metric set.

The result still does not support the tested empirical local-module mechanism.
The calibration-selected empirical-module candidate is not the best held-out
candidate, is not best on the local ruggedness component package, and remains
diffuse across empirical modules. The best selected fitted family is
no-module; the best Local7 selected fitted control is uncoupled-module; and the
empirical readout residual remains best overall. This weakens `HYP-008` as
currently operationalized while preserving only a narrower capacity claim:
ABC-fitted auxiliary landscapes can improve PHOT graphlet/localization metrics,
but the current empirical module construction is not validated as the source of
localized ruggedness.

## Limitations

- PHOT is the only empirical system in this run.
- The result uses one fixed calibration/test graphlet split.
- The support decision is based on calibration-selected particles; test-oracle
  posterior particles are reported only as capacity diagnostics.
- The local ruggedness metrics are graphlet-local summaries, not full
  landscape-wide localization estimates.
- The empirical target is the processed `latent_phi` layer.
