# RUN-096 Preflight And Collection Audit

| Item | Value |
|---|---:|
| Experiment | EXP-038 |
| Run | RUN-096 |
| Scheduler job | 246 |
| Inferred exit code | 0 |
| Assay | PHOT_CHLRE_Chen_2023 |
| Graphlets | 24 |
| Calibration graphlets | 12 |
| Held-out test graphlets | 12 |
| Objective metrics | 16 |
| Local ruggedness objective metrics | 7 |
| Fitted families | no_module, empirical_module, random_module, shuffled_module, uncoupled_module |
| Local output directory | data/processed/phot-abc-local-ruggedness-objective/RUN-096 |
| DVC pointer | data/processed/phot-abc-local-ruggedness-objective/RUN-096.dvc |
| DVC md5 | 4cd3edc3c50e33ecde7c5c7e32e9ec2a.dir |
| DVC files | 29 |
| DVC bytes | 115223180 |

The run extends the RUN-091 ABC objective by adding adjacent-edge
Dirichlet-energy Moran I and edge/site/plaquette concentration metrics. The
outputs were collected with `labproj collect` and tracked/backed up with
`labproj data track`.
