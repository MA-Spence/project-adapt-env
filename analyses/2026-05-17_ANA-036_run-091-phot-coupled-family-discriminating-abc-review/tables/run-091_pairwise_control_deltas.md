# RUN-091 Held-Out Pairwise Control Deltas

| candidate | reference | candidate test distance | reference test distance | delta | candidate better |
| --- | --- | --- | --- | --- | --- |
| empirical module | fitted current | 0.2553 | 0.4016 | -0.1463 | yes |
| empirical module | smooth core | 0.2553 | 0.4119 | -0.1566 | yes |
| empirical module | no module | 0.2553 | 0.3489 | -0.0937 | yes |
| empirical module | random module | 0.2553 | 0.2552 | 9.722e-05 | no |
| empirical module | shuffled module | 0.2553 | 0.3535 | -0.0982 | yes |
| empirical module | uncoupled module | 0.2553 | 0.2237 | 0.0316 | no |
| no module | fitted current | 0.3489 | 0.4016 | -0.0526 | yes |
| no module | smooth core | 0.3489 | 0.4119 | -0.0630 | yes |
| no module | empirical module | 0.3489 | 0.2553 | 0.0937 | no |
| no module | random module | 0.3489 | 0.2552 | 0.0938 | no |
| no module | shuffled module | 0.3489 | 0.3535 | -0.0045 | yes |
| no module | uncoupled module | 0.3489 | 0.2237 | 0.1253 | no |
| uncoupled module | fitted current | 0.2237 | 0.4016 | -0.1779 | yes |
| uncoupled module | smooth core | 0.2237 | 0.4119 | -0.1882 | yes |
| uncoupled module | empirical module | 0.2237 | 0.2553 | -0.0316 | yes |
| uncoupled module | random module | 0.2237 | 0.2552 | -0.0315 | yes |
| uncoupled module | no module | 0.2237 | 0.3489 | -0.1253 | yes |
