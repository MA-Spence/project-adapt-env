# RUN-084 HYP-008 decision table

| Criterion | Result | Interpretation |
|---|---|---|
| Procedural validity | Pass | Run completed, expected outputs present, data DVC-tracked. |
| Calibration graphlet improvement | Pass | `abc_best_particle` calibration distance `0.6100` beats fitted current `0.7318` and smooth core `0.7454`. |
| Calibration-selected held-out improvement | Fail | `abc_best_particle` test distance `0.4434` is worse than fitted current `0.4016` and smooth core `0.4119`. |
| Posterior mean held-out improvement | Fail | Posterior mean test distance `0.4543` is worse than fitted current and smooth core. |
| Test-oracle capacity | Weak positive, not support | `abc_particle_rank_023` test distance `0.3964` narrowly beats fitted current, but its calibration distance is poor (`0.7126`). |
| Assay guardrail | Weak positive | Posterior mean assay objective `0.4570` is slightly better than controls, but validation correlations and epistasis prediction do not improve as a package. |
| Spectral/topology package | Fail | Calibration-selected soft-min has held-out WHT high-order fraction `0.5268`, accessibility `0.0037`, and lag-1 autocorrelation `0.0405`, far from empirical `0.3788`, `0.1975`, and `0.1905`. |
| Locality | Fail | Innovation energy is highly concentrated but too large and weakly related to shared features; concentration is not empirical-like locality. |
| HYP-008 effect | Weakens | This implementation does not pass the held-out graphlet and locality controls required for support. |
