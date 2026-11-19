# RUN-059 Walk Summary By Mode And Stratum

| Mode | Start stratum | n | Mean steps | Mean total gain | Median decel | Random SD | Lag1 autocorr | Lag5 autocorr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| greedy_adaptive | far_consensus | 60 | 16.15 | 3290 | 0.6421 |  |  |  |
| greedy_adaptive | intermediate_consensus | 60 | 5.4333 | 1046 | 0.9430 |  |  |  |
| greedy_adaptive | near_consensus | 60 | 2.2667 | 518.6 | 1.0967 |  |  |  |
| random_beneficial_adaptive | far_consensus | 60 | 18.48 | 2160 | 1.3860 |  |  |  |
| random_beneficial_adaptive | intermediate_consensus | 60 | 6.3333 | 1120 | 1.2342 |  |  |  |
| random_beneficial_adaptive | near_consensus | 60 | 2.5333 | 548.1 | 1.3379 |  |  |  |
| unconditional_random | far_consensus | 60 | 120.0 | -21814 |  | 231.5 | 0.9993 | 0.9967 |
| unconditional_random | intermediate_consensus | 60 | 120.0 | -22549 |  | 231.7 | 0.9994 | 0.9969 |
| unconditional_random | near_consensus | 60 | 120.0 | -23199 |  | 230.9 | 0.9994 | 0.9969 |

Adaptive total gains are largest from far-consensus starts and greedy walks retain deceleration. Unconditional random walks remain highly autocorrelated in every stratum.
