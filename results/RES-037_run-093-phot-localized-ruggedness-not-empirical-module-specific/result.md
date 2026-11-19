# RES-037: RUN-093 PHOT localized ruggedness is real but not empirical-module specific

## Summary

RUN-093 finds that empirical PHOT graphlets and the empirical readout residual show stronger local Dirichlet-energy concentration than fitted current/smooth-core controls, but the RUN-091 empirical-module ABC candidate and empirical calibration modules do not uniquely explain that localization. No-module, random-module, uncoupled, residual, and GNK controls are competitive or better on key localization metrics, so the result weakens the current empirical-module implementation of HYP-008 while preserving a narrower residual-localization diagnosis.

## Generated from

- Analyses: ANA-037

## Relevant hypotheses

- Supports: None
- Weakens: HYP-008
- Refutes: None

## Evidence

- analyses/2026-05-17_ANA-037_run-093-phot-localized-ruggedness-energy-diagnostic-review/metrics.json
- analyses/2026-05-17_ANA-037_run-093-phot-localized-ruggedness-energy-diagnostic-review/tables/run-093_hyp008_decision_table.md
- analyses/2026-05-17_ANA-037_run-093-phot-localized-ruggedness-energy-diagnostic-review/tables/run-093_preflight_and_collection_audit.md
- analyses/2026-05-17_ANA-037_run-093-phot-localized-ruggedness-energy-diagnostic-review/tables/run-093_test_innovation_localization_summary.md
- analyses/2026-05-17_ANA-037_run-093-phot-localized-ruggedness-energy-diagnostic-review/tables/run-093_test_module_enrichment_summary.md
- analyses/2026-05-17_ANA-037_run-093-phot-localized-ruggedness-energy-diagnostic-review/tables/run-093_test_total_localization_summary.md
- results/RES-037_run-093-phot-localized-ruggedness-not-empirical-module-specific/metrics.json
- data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093.dvc
- data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093/summary.json
- data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093/model_localization_summary.csv
- data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093/module_enrichment_summary.csv

## Interpretation

`RUN-093` supports the narrower diagnostic statement that PHOT ruggedness is
partly localized in the processed graphlet layer and that the fitted current
core under-represents this edge/site concentration. On held-out graphlets, the
empirical target has edge top-10 energy fraction `0.3904` and site top-1 energy
fraction `0.3863`, compared with `0.3085` and `0.2645` for fitted
current/smooth core. Empirical edge effective fraction is lower (`0.4085`
versus `0.4558`), also indicating more concentrated edge energy. The empirical
readout residual shows a local signal too: adjacent-edge Moran I is `0.0505`
and edge effective fraction is `0.3903`.

The empirical-module mechanism does not pass the test. The calibration-selected
RUN-091 empirical-module candidate does not uniquely match empirical
localization: held-out edge Moran I is `0.0045` versus empirical `0.0126`, edge
top-10 fraction is `0.3734` versus `0.3904`, and site top-1 fraction is
`0.3024` versus `0.3863`. Simpler or control families are closer on key
metrics: no-module is closest for edge Moran I and edge effective fraction,
random-module is closest for edge top-10 and site effective fraction, and the
empirical readout residual is closest for site top-1 and plaquette-effective
fractions.

The module-enrichment evidence is also negative. Empirical PHOT energy is not
enriched in the empirical calibration module set: median enrichment is `0.8624`
for module-covered edge sites and `0.9495` for module-covered plaquette pairs.
The selected empirical-module candidate has only weak enrichment in those same
modules (`1.0000` edge-site, `1.1119` plaquette-pair), similar to fitted-current
and control families. On innovation fields, empirical-module plaquette
enrichment in empirical modules is `1.1295`, but the uncoupled-module
innovation is higher (`1.1697`) and shuffled-module innovation is stronger in
its own shuffled module set (`1.6543`).

The correct interpretation is therefore mixed. PHOT does show localized
ruggedness relative to the fitted current/smooth core, but `RUN-093` does not
identify the empirical calibration modules, or the RUN-091 empirical-module ABC
candidate, as the mechanism generating that localization.

## Effect on hypothesis

`RES-037` weakens `HYP-008` as currently operationalized. It supports the
motivation that the PHOT empirical landscape contains residual localized
ruggedness not captured by the fitted current core, but it weakens the stronger
claim that the current empirical-module nonlinear operator explains that
localized ruggedness. A mechanism-specific result would require empirical
modules to concentrate empirical energy and the empirical-module fitted family
to beat no-module, random, shuffled, uncoupled, residual-only, and GNK controls
on localization metrics. `RUN-093` does not meet that standard.

## Limitations

- PHOT is the only empirical system in this diagnostic.
- The diagnostic is graphlet-local, not a full landscape-wide localization
  analysis.
- The empirical module sets cover large fractions of graphlet sites and pairs,
  which reduces enrichment-test resolution.
- The result reuses `RUN-091` fitted parameter sets and does not refit an
  operator against localization metrics.
- The empirical target is the processed `latent_phi` layer.

## Downstream use

- Cite `RES-037` as evidence that PHOT empirical graphlets and empirical
  readout residuals contain localized Dirichlet-energy structure relative to
  fitted current/smooth core.
- Do not cite `RUN-093` as support for the tested empirical-module mechanism.
  The empirical modules do not enrich empirical energy, and controls match or
  exceed the empirical-module family on key localization summaries.
- Treat future PHOT local-ruggedness work as a localization-discovery problem,
  not as validation of the current empirical calibration module set.
