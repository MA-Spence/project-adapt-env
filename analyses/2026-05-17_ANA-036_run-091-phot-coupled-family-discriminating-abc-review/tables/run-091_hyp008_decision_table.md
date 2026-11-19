# RUN-091 HYP-008 Decision Table

| criterion | observation | passes |
| --- | --- | --- |
| ABC empirical module beats fitted-current and smooth-core controls | empirical module test distance 0.2553; fitted current 0.4016; smooth core 0.4119 | yes |
| Empirical module beats no-module capacity control | empirical module 0.2553; no module 0.3489; delta -0.0937 | yes |
| Empirical module beats random-module specificity control | empirical module 0.2553; random module 0.2552; delta 0.000097 | no |
| Empirical module beats shuffled-module specificity control | empirical module 0.2553; shuffled module 0.3535; delta -0.0982 | yes |
| Empirical module beats uncoupled-landscape control | empirical module 0.2553; uncoupled module 0.2237; delta 0.0316 | no |
| Calibration-selected empirical module satisfies component guardrails | 0 violations; guardrail penalty 0.0000 | yes |
| Selected empirical-module energy is locally concentrated | plaquette effective modules 9.11/12; plaquette normalized entropy 0.937 | no |
