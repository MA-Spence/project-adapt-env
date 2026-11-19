# ANA-037: RUN-093 PHOT localized ruggedness energy diagnostic review

## Purpose

Review RUN-093 localized Dirichlet-energy diagnostics for PHOT empirical, fitted current a(x), smooth-core controls, RUN-091 ABC fitted families, and random/shuffled/uncoupled controls, focusing on whether ruggedness is localized or diffuse.

## Linked experiments/runs

- Experiments: `EXP-037`
- Runs: `RUN-093`
- Source data:
  `data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093`
- DVC pointer:
  `data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093.dvc`

## Scientific context

`EXP-037` was registered after `RUN-091` showed that ABC fitting can improve
fixed held-out PHOT graphlet distance, but without empirical-module specificity:
random-module was fractionally better, uncoupled-module was clearly better, and
module-energy summaries remained diffuse. `RUN-093` therefore does not test a
new operator. It reuses the `RUN-091` fitted ABC parameter sets and asks whether
ruggedness itself is locally concentrated in PHOT empirical graphlets, the
frozen fitted current `a(x)`, and the fitted/control families.

The primary diagnostic is edge Dirichlet energy on each dense graphlet edge,
`(f(v)-f(u))^2`. Localization is read through adjacent-edge energy
autocorrelation, edge/node/site/plaquette energy concentration, and enrichment
inside the empirical, random, shuffled, and uncoupled module sets.

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-093_preflight_and_collection_audit.md`
- `tables/run-093_test_total_localization_summary.md`
- `tables/run-093_test_innovation_localization_summary.md`
- `tables/run-093_test_module_enrichment_summary.md`
- `tables/run-093_hyp008_decision_table.md`

## Run and data audit

`RUN-093` completed with scheduler job `244` and inferred Slurm exit code `0`.
The collected output directory contains 12 files and is now tracked and backed
up with DVC at
`data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093.dvc`
(`md5 8747ed319b4f1a5fed2230d543cc5d5b.dir`, 12 files, 175,474,957 bytes).

The run reused `RUN-091` ABC output tables and `RUN-077` PHOT source graphlet
context. It evaluated 24 complete PHOT graphlets, split into 12 calibration and
12 held-out test graphlets, with 80 reused ABC parameter sets plus static
empirical/current/smooth/residual/GNK controls. Primary comparisons use the
calibration-selected `RUN-091` candidates; posterior rows are retained as fitted
particle spread.

## Key observations

- The empirical PHOT `latent_phi` landscape is more locally concentrated than
  the fitted current/smooth core on held-out graphlets. Empirical test edge
  top-10 energy fraction is `0.3904` versus `0.3085` for fitted current/smooth,
  and empirical site top-1 energy fraction is `0.3863` versus `0.2645`. The
  edge effective fraction is lower for empirical (`0.4085`) than fitted
  current/smooth (`0.4558`), consistent with more concentrated edge energy.
- The empirical readout residual is also locally structured on several readouts:
  held-out adjacent-edge Moran I is `0.0505` versus empirical `0.0126`; edge
  effective fraction is `0.3903` versus empirical `0.4085`; and plaquette
  effective fraction is `0.8664` versus empirical `0.8437`. This supports the
  idea that PHOT has residual local ruggedness beyond the fitted current core.
- The selected empirical-module ABC candidate does not specifically match the
  empirical localization package. On held-out total landscapes, its edge Moran I
  is `0.0045`, edge top-10 fraction `0.3734`, site top-1 fraction `0.3024`, and
  plaquette effective fraction `0.7048`, all worse than either the empirical
  target or simpler controls for key comparisons.
- No-module and random-module controls are at least as competitive as the
  empirical-module candidate. No-module is closest to empirical for edge Moran I
  and edge effective fraction; random-module is closest for edge top-10 fraction
  and site effective fraction. The empirical-module family is not the uniquely
  localized model under these metrics.
- The module-enrichment test is negative. PHOT empirical held-out energy is not
  enriched in the empirical calibration module set: median empirical-module
  site enrichment is `0.8624` and plaquette-pair enrichment is `0.9495`. The
  selected empirical-module model shows only weak empirical-module enrichment
  (`1.0000` edge-site, `1.1119` plaquette-pair), comparable to fitted-current,
  no-module, random, shuffled, and uncoupled controls.
- The innovation fields do not rescue the empirical-module mechanism. The
  empirical-module innovation has negative held-out edge Moran I (`-0.0045`),
  while no-module, random, shuffled, and uncoupled innovations are also
  localized or overconcentrated on some concentration metrics. Empirical-module
  innovation plaquette enrichment in empirical modules is `1.1295`, but
  uncoupled innovation is higher (`1.1697`) and shuffled-module innovation has
  stronger enrichment in its own shuffled pairs (`1.6543`).

## Interpretation

`RUN-093` gives evidence that PHOT ruggedness is not purely diffuse in the
processed graphlet layer: empirical PHOT and its empirical readout residual
show stronger edge/site concentration than the fitted current/smooth core. That
part is directionally consistent with the motivation for `HYP-008`.

It does not support the tested empirical-module mechanism. The localized
ruggedness is not preferentially captured by the empirical-module ABC candidate
or by the empirical calibration module set. The best-matching localization
metrics are split across empirical readout residual, no-module, random-module,
and GNK controls. Because the empirical modules do not isolate the empirical
energy, the evidence weakens `HYP-008` as currently operationalized while
preserving only a narrower statement: PHOT has residual localized ruggedness
that the current fitted core under-represents.

## Limitations

- PHOT is the only empirical system in this diagnostic.
- The diagnostic operates on dense graphlets rather than full-sequence-space
  landscapes.
- The module enrichment measures use RUN-091 empirical/random/shuffled module
  sets, which cover large fractions of graphlet sites and pairs; this limits
  their power to distinguish very local mechanisms.
- The result tests localization structure only. It does not refit the model or
  repeat held-out graphlet distance selection.
