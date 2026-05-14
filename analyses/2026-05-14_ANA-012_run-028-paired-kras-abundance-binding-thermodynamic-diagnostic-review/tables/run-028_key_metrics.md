# RUN-028 Key Metrics

## Assay Panel And Observation Layer

| assay | role | variants | singles | multiples | MAVE-NN test Spearman | MAVE-NN test NRMSE |
|---|---:|---:|---:|---:|---:|---:|
| `RASK_HUMAN_Weng_2022_abundance` | abundance | 26,012 | 3,066 | 22,946 | 0.818 | 0.614 |
| `RASK_HUMAN_Weng_2022_binding-DARPin_K55` | binding | 24,873 | 3,084 | 21,789 | 0.933 | 0.362 |

## Branch Validation Metrics

Lower KS and NRMSE are better; higher Spearman is better.

| branch | paired abundance? | thermodynamic readout? | abundance KS | binding KS | single holdout Spearman | double holdout Spearman | epistasis Spearman | epistasis KS | reference fraction of peak | reference distance to peak |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `binding_only_raw` | no | no | NA | 0.506 | 0.308 | 0.195 | 0.106 | 0.673 | 1.000 | 0 |
| `binding_only_biophysical_binding` | no | yes | NA | 0.384 | 0.293 | 0.298 | 0.172 | 0.481 | 0.999 | 1 |
| `paired_abundance_binding_raw` | yes | no | 0.582 | 0.381 | 0.265 | 0.252 | 0.104 | 0.488 | 0.999 | 6 |
| `paired_abundance_binding_biophysical_binding` | yes | yes | 0.832 | 0.258 | 0.105 | 0.144 | 0.287 | 0.428 | 0.999 | 6 |

## Main Contrasts

| contrast | favorable changes | unfavorable changes |
|---|---|---|
| `binding_only_biophysical_binding` vs `binding_only_raw` | binding KS improves 0.506 to 0.384; double holdout Spearman improves 0.195 to 0.298; epistasis Spearman improves 0.106 to 0.172; epistasis KS improves 0.673 to 0.481 | single holdout Spearman falls 0.308 to 0.293; reference remains essentially at peak |
| `paired_abundance_binding_raw` vs `binding_only_raw` | binding KS improves 0.506 to 0.381; double holdout Spearman improves 0.195 to 0.252 | single holdout Spearman falls 0.308 to 0.265; epistasis Spearman is unchanged to slightly worse; reference remains near peak |
| `paired_abundance_binding_biophysical_binding` vs `paired_abundance_binding_raw` | binding KS improves 0.381 to 0.258; epistasis Spearman improves 0.104 to 0.287; epistasis KS improves 0.488 to 0.428 | abundance KS worsens 0.582 to 0.832; single holdout Spearman falls 0.265 to 0.105; double holdout Spearman falls 0.252 to 0.144; reference remains near peak |
