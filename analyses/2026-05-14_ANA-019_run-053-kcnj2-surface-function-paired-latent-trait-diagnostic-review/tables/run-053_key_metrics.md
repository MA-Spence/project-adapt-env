# RUN-053 key metrics

## Assay and paired-readout setup

| Quantity | Value |
| --- | ---: |
| State assay | `KCNJ2_MOUSE_Coyote-Maestas_2022_surface` |
| Function assay | `KCNJ2_MOUSE_Coyote-Maestas_2022_function` |
| Sequence length | 428 |
| State single mutants | 6917 |
| Function single mutants | 6963 |
| Multiple mutants | 0 |
| Matched variants for readout correlation | 6789 |
| State/function Spearman | 0.289606296146865 |
| State MAVE-NN test Spearman | 0.5375325146698333 |
| Function MAVE-NN test Spearman | 0.3889332346846113 |
| State MAVE-NN test NRMSE | 0.8899418510969036 |
| Function MAVE-NN test NRMSE | 0.9446255200804048 |

## Branch validation metrics

| Branch | State used | Function holdout Spearman | Function holdout NRMSE | Function KS | Function fraction of peak | Reference distance to peak | State unfolding KS | Fitted dims | Fitted sigma |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `function_only_single_trait` | no | 0.0664910504549503 | 1.0004281797796404 | 0.4596797534067182 | 0.9992867218706184 | 4 |  | 1 | 10.0 |
| `paired_state_function_single_trait` | yes | 0.06523696551447947 | 1.0021137445548625 | 0.5022009738920055 | 0.9997684956545109 | 4 | 0.4925828072740223 | 2 | 28.0 |
| `paired_state_function_explicit_two_trait` | yes | 0.0213506017358205 | 0.9984805162737627 | 0.5105557949159845 | 0.9931934002597993 | 2 | 0.49136806662883015 | 1 | 10.0 |

## Fitted parameter checks

| Branch | Epistasis strength | Empirical pairwise strength | Noise amplitude | Peak distance setting |
| --- | ---: | ---: | ---: | ---: |
| `function_only_single_trait` | 0.0 | 0.0 | 0.0 | 4 |
| `paired_state_function_single_trait` | 0.0 | 0.0 | 0.0 | 2 |
| `paired_state_function_explicit_two_trait` | 0.0 | 0.0 | 0.0 | 1 |

## Interpretation

- The paired readouts are weakly correlated, so surface trafficking could in
  principle add nonredundant state information.
- The function observation layer is weak, with MAVE-NN test Spearman only
  `0.389`, which limits the strength of the negative inference.
- The paired branches did not translate the surface information into better
  function recovery. The explicit two-trait branch had a marginal NRMSE gain
  but worse holdout Spearman and functional KS than the function-only branch.
- The result is limited to single-mutant recovery and paired-readout structure;
  no double-mutant or epistasis evidence is available for this assay pair.
