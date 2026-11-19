| test | status | evidence |
| --- | --- | --- |
| RUN-065 completed and outputs were collected | PASS | scheduler exit code 0; durable output directory present |
| Durable RUN-065 outputs are DVC-tracked | PASS | data/processed/tem1-aggregate-graphlet-residual-ruggedness-generalization/RUN-065.dvc |
| Graphlet panel spans dramatically different starting classes | FAIL | starts and backgrounds span five classes, but graphlet_catalog contains {'reference_or_consensus_like': 240} |
| Primary residual improves held-out topology over smooth core | PASS | test SE 0.3213 vs 0, RSE 0.1030 vs 0, local maxima 3.54 vs 1, accessible paths 0.0350 vs 0.3097, lag1 autocorrelation 0.4976 vs 0.6455 |
| Primary residual calibration transfers to held-out graphlets | PASS | test metrics closely match calibration: SE +0.0011, RSE +0.0074, maxima -0.12, high-order Walsh +0.0049 |
| Primary residual improves over RES-024 exact-match topology | MIXED | primary has higher test SE/RSE/maxima and lower high-order dominance, but RES-024 exact-match is similarly short-correlated and has lower accessible paths |
| DCA/contact placement contributes beyond matched controls | FAIL | shuffled and random-position controls closely match primary residual on test SE/RSE/maxima and spectral fractions |
| GNK local null is not a stronger comparator | FAIL | GNK has much stronger topology but non-procedural construction and Walsh high-order fraction 0.4532 on test graphlets |
| Pair-only DCA control does not solve the same topology package | FAIL | DCA pair override beats primary topology on test SE/RSE/maxima and accessibility while retaining essentially zero high-order Walsh tail |
| Residual-only behavior survives embedding in the smooth core | MIXED | residual-only pset-000 has about 50.6 test local maxima without the core, while embedded weak residual psets remain near one maximum; the stronger selected primary pset does survive partially |
| Primary residual spectra are low-order dominated with a high-order tail | PASS | test Walsh high-order fraction 0.0737; GFT low/mid/high energy 0.7719/0.2197/0.0084 |
