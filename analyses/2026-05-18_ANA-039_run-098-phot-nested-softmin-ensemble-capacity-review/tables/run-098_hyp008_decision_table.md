# RUN-098 HYP-008 decision table

| test | heldout_result | interpretation |
| --- | --- | --- |
| Empirical module vs locked first stage | 0.543180 vs 0.529590; delta +0.013590 | Fails. The calibration-selected empirical module is worse on the primary held-out full objective. |
| Empirical module vs fitted current/smooth core | 0.543180 vs 0.521007/0.526819 | Fails. The empirical module does not beat either non-nested PHOT baseline. |
| Empirical module vs no/random/shuffled controls | loses by +0.005981/+0.060549/+0.032219 | Fails empirical-module specificity. |
| Empirical module vs uncoupled module | wins by -0.096961 | One negative control is beaten, but this is insufficient because no/random/shuffled controls beat empirical-module. |
| Best calibration-selected fitted family | random_module = 0.482630 | Supports generic latent capacity only; the winner is a random-module negative control. |
| Local7 held-out package | empirical_module = 0.525045; no_module = 0.476565; shuffled = 0.419862; fitted current/smooth = 0.674567 | Nested capacity helps some local-ruggedness distances, but not specifically through empirical modules. |
| Empirical readout residual diagnostic | 0.211455 full; 0.135199 Local7 | Shows strong residual structure exists, but it is an empirical diagnostic/oracle comparator, not a generative model. |
| Test-oracle posterior particles | no_module = 0.395321; empirical_module = 0.409176 | There is posterior capacity, but calibration did not select those particles; this is not held-out support. |
