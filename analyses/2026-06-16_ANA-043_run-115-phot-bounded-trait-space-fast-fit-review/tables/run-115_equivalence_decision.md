# RUN-115 Operational Equivalence Decision

`RUN-115` selected `12` bounded trait-space
candidates and produced `12` equivalence rows
against the best available `RUN-098` full-ABC comparator. Operationally
equivalent rows: `0`.

The bounded candidates passed the component-distance guardrail (`12/12` rows),
but failed both the absolute and relative aggregate-distance margins against
the best full-ABC comparator (`12/12` failures for each margin). This means the
minimal bounded route reduced the gross component pathologies enough to satisfy
component guardrails, but did not become equivalent to the previous full-ABC
PHOT comparator.

| fit_method               |   n_rows |   operationally_equivalent |   absolute_margin_pass |   relative_margin_pass |   component_guardrail_pass |   best_candidate_test_distance |   median_candidate_test_distance |   reference_test_distance |   max_component_distance |   best_delta_vs_reference |   best_relative_delta_vs_reference |
|:-------------------------|---------:|---------------------------:|-----------------------:|-----------------------:|---------------------------:|-------------------------------:|---------------------------------:|--------------------------:|-------------------------:|--------------------------:|-----------------------------------:|
| local_refinement_round_1 |       11 |                          0 |                      0 |                      0 |                         11 |                       0.426882 |                         0.462198 |                  0.211455 |                  1.36614 |                  0.215426 |                            1.01878 |
| random_prior             |        1 |                          0 |                      0 |                      0 |                          1 |                       0.533371 |                         0.533371 |                  0.211455 |                  1.3668  |                  0.321915 |                            1.52238 |
