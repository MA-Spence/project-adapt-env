# RUN-051 key metrics

## Assay and paired-readout setup

| Quantity | Value |
| --- | ---: |
| State assay | `OXDA_RHOTO_Vanella_2023_expression` |
| Function assay | `OXDA_RHOTO_Vanella_2023_activity` |
| Sequence length | 364 |
| State single mutants | 6769 |
| Function single mutants | 6396 |
| Multiple mutants | 0 |
| Matched variants for readout correlation | 6387 |
| State/function Spearman | 0.245668548592404 |
| State MAVE-NN test Spearman | 0.5667081580961366 |
| Function MAVE-NN test Spearman | 0.542045485915721 |
| State MAVE-NN test NRMSE | 0.8647862904455272 |
| Function MAVE-NN test NRMSE | 0.8931811662990445 |

## Branch validation metrics

| Branch | State used | Function holdout Spearman | Function holdout NRMSE | Function KS | Function fraction of peak | Reference distance to peak | State unfolding KS | Fitted dims | Fitted sigma |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `function_only_single_trait` | no | 0.18332410997777743 | 0.9809534666770041 | 0.5198561601000625 | 0.9995343920625102 | 1 |  | 4 | 28.0 |
| `paired_state_function_single_trait` | yes | 0.06092087621657774 | 0.9990106501404143 | 0.5231512176313203 | 0.9995139950916463 | 6 | 0.4110663461613225 | 1 | 28.0 |
| `paired_state_function_explicit_two_trait` | yes | 0.18001343719972057 | 0.984857405282575 | 0.5261100687929956 | 0.9999185560270659 | 3 | 0.42424906231851206 | 1 | 10.0 |

## Fitted parameter checks

| Branch | Epistasis strength | Empirical pairwise strength | Noise amplitude | Peak distance setting |
| --- | ---: | ---: | ---: | ---: |
| `function_only_single_trait` | 0.0 | 0.0 | 0.0 | 1 |
| `paired_state_function_single_trait` | 0.0 | 0.0 | 0.0 | 2 |
| `paired_state_function_explicit_two_trait` | 0.0 | 0.0 | 0.0 | 2 |

## Interpretation

- The paired readouts are weakly correlated, so expression could in principle
  add nonredundant state information.
- The paired branches did not translate that information into better activity
  recovery. The explicit two-trait branch nearly matched the function-only
  ranking metric but did not beat it on holdout Spearman, NRMSE, or KS.
- The result is limited to single-mutant recovery and paired-readout structure;
  no double-mutant or epistasis evidence is available for this assay pair.
