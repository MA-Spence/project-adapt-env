# RUN-086 HYP-008 Decision Table

| Criterion | Result | Decision |
|---|---|---|
| Calibration graphlet distance improves over controls | ABC best `0.7240` beats fitted current `0.7318` and smooth core `0.7454` only slightly | Weak technical pass |
| Calibration-selected model beats fitted current on held-out graphlets | ABC best `0.4060` versus fitted current `0.4016` | Fail |
| Module-specific effect beats controls | Shuffled-module held-out distance `0.4029` is better than empirical-module ABC best `0.4060` | Fail |
| Empirical graphlet package recovered | Selected model keeps low accessibility and low local maxima relative to empirical graphlets | Fail |
| Innovation is local and non-diffuse | Innovation ratios to empirical residual are near zero | Fail |
| Assay guardrail improves | Best control objective `0.4631`; best budgeted objective `0.5350`; selected particle `0.5560` | Fail |
| Overall effect on `HYP-008` | Budgeted module-tied operator does not pass held-out graphlet, locality, or assay controls | Weakens |

