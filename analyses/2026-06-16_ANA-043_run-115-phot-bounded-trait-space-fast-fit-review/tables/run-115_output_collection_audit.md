# RUN-115 Output Collection Audit

Runtime output directory: `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs`

Copied compact payload: `940677` bytes.

Large or unnecessary files were left on the runtime volume and are referenced by
absolute runtime path below. In particular, `walsh_terms.csv` and
`gft_coefficients.csv` were not needed for this ANA because their aggregate
spectral summaries were already present in `graphlet_summary.csv` and
`graphlet_metrics.csv`.

| file                                         |   bytes | copied_to_analysis   | analysis_copy                                                      | runtime_path                                                                                                   |
|:---------------------------------------------|--------:|:---------------------|:-------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|
| PHOT_CHLRE_Chen_2023.alignment.fasta         |   74752 | False                |                                                                    | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/PHOT_CHLRE_Chen_2023.alignment.fasta         |
| abc_comparator_selected_family_distances.csv |   18173 | True                 | source_tables/run-115_abc_comparator_selected_family_distances.csv | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/abc_comparator_selected_family_distances.csv |
| bounded_gate_candidates.csv                  |  116295 | True                 | source_tables/run-115_bounded_gate_candidates.csv                  | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/bounded_gate_candidates.csv                  |
| bounded_gate_distances.csv                   |   25846 | True                 | source_tables/run-115_bounded_gate_distances.csv                   | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/bounded_gate_distances.csv                   |
| bounded_gate_equivalence.csv                 |    4444 | True                 | source_tables/run-115_bounded_gate_equivalence.csv                 | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/bounded_gate_equivalence.csv                 |
| gft_coefficients.csv                         | 6404642 | False                |                                                                    | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/gft_coefficients.csv                         |
| graphlet_catalog.csv                         |    2209 | True                 | source_tables/run-115_graphlet_catalog.csv                         | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/graphlet_catalog.csv                         |
| graphlet_metrics.csv                         |  713761 | True                 | source_tables/run-115_graphlet_metrics.csv                         | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/graphlet_metrics.csv                         |
| graphlet_summary.csv                         |   36759 | True                 | source_tables/run-115_graphlet_summary.csv                         | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/graphlet_summary.csv                         |
| module_catalog.csv                           |    4004 | True                 | source_tables/run-115_module_catalog.csv                           | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/module_catalog.csv                           |
| preflight.json                               |    2209 | True                 | source_tables/run-115_preflight.json                               | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/preflight.json                               |
| selected_bounded_gate_candidates.csv         |   15213 | True                 | source_tables/run-115_selected_bounded_gate_candidates.csv         | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/selected_bounded_gate_candidates.csv         |
| summary.json                                 |    1764 | True                 | source_tables/run-115_summary.json                                 | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/summary.json                                 |
| walsh_terms.csv                              | 6167423 | False                |                                                                    | /tank/data/labproj/runtime/project-adapt-env/runs/RUN-115/outputs/walsh_terms.csv                              |
