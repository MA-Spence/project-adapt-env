| criterion | RUN-093 observation | decision |
|---|---|---|
| Empirical PHOT shows localized ruggedness beyond fitted current/smooth core | Empirical edge top10 fraction 0.3904 vs current/smooth 0.3085; empirical site top1 fraction 0.3863 vs current/smooth 0.2645; empirical edge effective fraction 0.4085 vs current/smooth 0.4558 | supports motivation |
| Empirical readout residual contains local signal | Residual edge Moran I 0.0505, edge effective fraction 0.3903, plaquette effective fraction 0.8664 | supports residual-localization diagnosis |
| Empirical-module ABC candidate uniquely matches empirical localization | No-module is closest for edge Moran I and edge effective fraction; random-module is closest for edge top10 and site effective fraction; empirical-module is not uniquely best | fails |
| Empirical calibration modules enrich empirical energy | Empirical target enrichment is below one for empirical module sites and pairs: 0.8624 edge-site and 0.9495 plaquette-pair | fails |
| ABC empirical-module innovation shows specific empirical-module localization | Empirical-module innovation empirical-module plaquette enrichment is 1.1295, but uncoupled innovation is higher at 1.1697 and shuffled innovation is stronger in its own module set at 1.6543 | fails specificity |
| Overall effect on HYP-008 | PHOT ruggedness appears locally structured, but the tested empirical-module mechanism and module set do not explain it better than controls | weakens current implementation |

