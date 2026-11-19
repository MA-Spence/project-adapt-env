| test | status | detail |
| --- | --- | --- |
| live run completed | pass | status=completed, scheduler_job_id=223, exit=0 |
| durable outputs DVC-tracked | pass | data/processed/tem1-dca-clique-haplotype-high-order-epistasis-diagnostic/RUN-063.dvc |
| old generic epistasis disabled | pass | {"epistasis_covariation_bias": 0.0, "epistasis_strength": 0.0, "n_epistatic_terms": 0, "n_higher_order_epistatic_terms": 0} |
| DCA haplotype terms selected | pass | 48 terms; counts={2: 12, 3: 12, 4: 12, 5: 12} |
| zero-strength target-order signal absent | pass | max coeff=2.8422e-13; max variance fraction=1.6604e-29 |
| target Walsh coefficient follows exact-match scaling | pass | max \|observed - strength/2^order\|=2.8422e-13 |
| target Walsh coefficient monotonic with strength | pass | {"2": true, "3": true, "4": true, "5": true} |
| target-order variance fraction non-decreasing | partial | {"2": false, "3": true, "4": true, "5": true} |
| fifth-order upper signal present | pass | order-5 upper coeff=93.7500; variance fraction=0.0294 |
| RSE side effect remains weak | observed | order-5 upper sign=0.0611; RSE=0.0194 |
| topological ruggedness remains limited | limited | order-5 upper mean local maxima=1.5000; max=2; multi-maximum maps=6/12 |
