# RUN-096 HYP-008 Decision Table

| Criterion | Observation | Decision |
|---|---|---|
| Fit improves over fitted current/smooth core | No-module, uncoupled, empirical-module, random, and shuffled selected families all beat fitted current (`0.521007`) and smooth core (`0.526819`) on held-out full objective distance. | Capacity signal present |
| Local ruggedness terms improve over fitted current/smooth core | Local7 distances for selected fitted families are `0.289976` to `0.462290`, below fitted current/smooth core `0.674567`. | Directional local-metric improvement |
| Empirical-module family wins held-out selection | Empirical-module held-out full distance is `0.365272`, worse than no-module `0.277056` and uncoupled-module `0.291161`. | Fails |
| Empirical-module beats random/shuffled controls | Empirical-module narrowly beats random-module by `0.001496` and shuffled-module by `0.014186`, but this is not enough because no-module and uncoupled controls are better. | Partial, not sufficient |
| Empirical-module local metrics are best | Empirical-module Local7 distance is `0.462290`, worse than uncoupled `0.289976`, no-module `0.326809`, random `0.439467`, and shuffled `0.452016`. | Fails |
| Module energy is non-diffuse in empirical modules | Selected empirical-module test energy is diffuse: edge effective module count `11.010808/12`, edge normalized entropy `0.982059`, plaquette effective module count `9.341583/12`, plaquette normalized entropy `0.944129`. | Fails |
| Calibration-selected evidence only | A test-oracle empirical-module particle has test distance `0.259136`, but it was not calibration-selected. | Capacity only, not support |

Overall: RUN-096 supports the narrower statement that ABC-fitted coupled
families can move PHOT held-out graphlet and localization distances away from
the fitted current/smooth core. It does not support the tested empirical local
module mechanism required by HYP-008.
