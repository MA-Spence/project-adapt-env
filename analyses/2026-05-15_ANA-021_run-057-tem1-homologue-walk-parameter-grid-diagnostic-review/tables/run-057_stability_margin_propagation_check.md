# RUN-057 Stability-Margin Propagation Check

This check compares identical homologue starts across the nominal
`stability_margin` sweep at `epistasis_strength = 0.0`.

| Start | Parameter set | Nominal stability margin | Background fitness | stability_margin_used | Neutral fraction | Deleterious fraction |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| near_consensus__rep-0 | pset-000 | 4.0 | 8358.9749 | 835.8975 | 0.4818 | 0.5180 |
| near_consensus__rep-0 | pset-001 | 8.0 | 8358.9749 | 835.8975 | 0.4818 | 0.5180 |
| near_consensus__rep-0 | pset-002 | 12.0 | 8358.9749 | 835.8975 | 0.4818 | 0.5180 |
| far_consensus__rep-2 | pset-000 | 4.0 | 10.0000 | 1.0000 | 0.00147 | 0.8966 |
| far_consensus__rep-2 | pset-001 | 8.0 | 10.0000 | 1.0000 | 0.00147 | 0.8966 |
| far_consensus__rep-2 | pset-002 | 12.0 | 10.0000 | 1.0000 | 0.00147 | 0.8966 |

The nominal parameter value is recorded differently, but the realized background
fitness and `stability_margin_used` are identical for these matched starts. The
failed stability-margin validation checks therefore expose a parameterization or
propagation problem in this experiment, not a reliable biological negative.
