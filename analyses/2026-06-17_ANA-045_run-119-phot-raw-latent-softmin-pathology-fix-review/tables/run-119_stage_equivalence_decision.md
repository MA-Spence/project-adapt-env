# RUN-119 stage equivalence decision

_Operational equivalence is the same RUN-098 best-any-full-ABC decision used by the EXP-046 runner._

| fix_stage   | fixes_applied   |   n_candidate_evaluations |   n_selected_candidates |   best_objective_distance |   best_test_distance |   median_test_distance |   n_operationally_equivalent |   n_component_guardrail_pass |   n_absolute_margin_pass |   n_relative_margin_pass |   best_max_component_distance |
|:------------|:----------------|--------------------------:|------------------------:|--------------------------:|---------------------:|-----------------------:|-----------------------------:|-----------------------------:|-------------------------:|-------------------------:|------------------------------:|
| 1           | 1               |                        73 |                       8 |                     1.943 |                1.224 |                  1.312 |                            0 |                            1 |                        0 |                        0 |                         2.433 |
| 1+2         | 1+2             |                        73 |                       8 |                    21.21  |                1.687 |                  2.012 |                            0 |                            0 |                        0 |                        0 |                         4.486 |
| 1+2+3       | 1+2+3           |                        73 |                       8 |                    24.97  |                1.687 |                  2.012 |                            0 |                            0 |                        0 |                        0 |                         4.486 |
