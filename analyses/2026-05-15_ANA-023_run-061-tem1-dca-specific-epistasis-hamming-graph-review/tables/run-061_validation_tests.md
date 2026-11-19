# RUN-061 Validation Tests

| Test | Criterion | Pass | Value | Interpretation |
| --- | --- | --- | --- | --- |
| old_generic_epistasis_disabled | all legacy epistasis knobs are zero | true | true | Confirms the experiment isolates the DCA-prior override path. |
| complete_dense_hamming_cube | n_nodes == 2^k | true | 1024 == 2^10 | Confirms exhaustive enumeration of the local binary graph. |
| minimum_dca_pairs_selected | selected pairs >= preregistered minimum | true | 5 >= 5 | Confirms the DCA selector found the intended number of pairs. |
| selected_pairs_alignment_enriched | min log2 enrichment > 0 | true | 3.8124 | Selected state pairs have enriched joint support in the TEM-1 alignment. |
| dca_specific_pair_epistasis_present | mean abs pair epistasis > 1 | true | 57.0791 | The new architecture creates strong local pairwise epistasis. |
| selected_pair_sign_epistasis_present | sign epistasis fraction > 0.2 | true | 1.0000 | Selected DCA pairs show systematic sign epistasis. |
| reciprocal_sign_epistasis_present | reciprocal sign epistasis fraction > 0.05 | false | 0.0000 | Fails: this run creates one-sided compensatory sign epistasis, not reciprocal sign epistasis. |
| random_walk_autocorrelation_reduced | lag-1 autocorrelation < 0.95 | true | 0.7795 | The local cube is substantially less autocorrelated than the near-0.999 RUN-059 random-walk regime. |
| multiple_local_maxima_present | local maxima > 1 | false | 1 | Fails: the dense cube is a single-peak funnel, not a multi-peak rugged topology. |
| accessible_shortest_paths_to_all_alternatives | accessible path fraction > 0 | false | 0.0000 | Fails from the anchor: the all-selected-alternative genotype is lower fitness, so no monotone improving shortest path exists. |
| pairwise_walsh_variance_substantial | second-order variance fraction > 0.05 | false | 0.0089 | Fails: the selected pair terms are strong locally but additive single effects dominate total cube variance. |

These are review diagnostics for the scientific record, not a separate preregistered test suite. They distinguish specific pairwise/sign epistasis from broader topological ruggedness.
