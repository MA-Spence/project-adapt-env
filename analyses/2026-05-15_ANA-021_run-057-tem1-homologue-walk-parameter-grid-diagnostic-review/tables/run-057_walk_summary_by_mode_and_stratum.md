# RUN-057 Walk Summary By Mode And Start Stratum

| Mode | Stratum | Walks | Mean start distance | Mean start fitness | Mean steps | Mean total gain | Median deceleration | Random step SD | Lag-1 autocorr | Random beneficial | Random deleterious |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| greedy_adaptive | near_consensus | 54 | 151.2563 | 6983 | 2.6667 | 576.8939 | 1.0512 |  |  |  |  |
| greedy_adaptive | intermediate_consensus | 54 | 236.8735 | 5493 | 6.0000 | 1078 | 0.9430 |  |  |  |  |
| greedy_adaptive | far_consensus | 54 | 354.9412 | 2595 | 15.6667 | 3137 | 0.6240 |  |  |  |  |
| random_beneficial_adaptive | near_consensus | 54 | 151.2563 | 6983 | 3.1111 | 626.0298 | 1.2230 |  |  |  |  |
| random_beneficial_adaptive | intermediate_consensus | 54 | 236.8735 | 5493 | 7.4815 | 1194 | 1.1482 |  |  |  |  |
| random_beneficial_adaptive | far_consensus | 54 | 354.9412 | 2595 | 18.9444 | 2098 | 1.2686 |  |  |  |  |
| unconditional_random | near_consensus | 54 | 151.2563 | 6983 | 80.0000 | -16490 |  | 226.0147 | 0.9989 | 0.03657 | 0.6236 |
| unconditional_random | intermediate_consensus | 54 | 236.8735 | 5493 | 80.0000 | -16010 |  | 227.0448 | 0.9988 | 0.04167 | 0.5926 |
| unconditional_random | far_consensus | 54 | 354.9412 | 2595 | 80.0000 | -15460 |  | 223.0877 | 0.9987 | 0.03657 | 0.5563 |

Greedy adaptive walks from far-consensus starts are longer and gain more total
fitness, while near-consensus starts exhaust beneficial single mutants quickly.
Random walks remain highly autocorrelated across all strata, so the ruggedness
signal is directionally correct but small.
