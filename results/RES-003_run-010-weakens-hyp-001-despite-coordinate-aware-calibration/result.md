# RES-003: RUN-010 weakens HYP-001 despite coordinate-aware calibration

## Summary

RUN-010 used a homogeneous multi-mutant ProteinGym stability panel, raw assay scores, assay-specific MAVE-NN observation models, and real MMseqs alignments, but the shared Adapt-Env calibration still showed weak single- and double-mutant holdout recovery. The measurement-layer confound was reduced, yet the downstream landscape fit remained poor, so the evidence further weakens rather than supports HYP-001.

## Generated from

- Analyses: ANA-003

## Relevant hypotheses

- Supports: None
- Weakens: HYP-001
- Refutes: None

## Evidence

- `analyses/2026-05-08_ANA-003_run-010-raw-scale-latent-observation-calibration-review/tables/run-010_key_metrics.md`
- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/summary.json`
- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/selected_panel.csv`
- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/mavenn_assay_metrics.csv`
- `data/processed/proteingym-raw-scale-latent-observation-calibration-panel/RUN-010/per_assay_latent_fits.csv`
- `experiments/2026-05-08_EXP-003_proteingym-raw-scale-latent-observation-calibration-panel/config.yaml`
- `experiments/2026-05-08_EXP-003_proteingym-raw-scale-latent-observation-calibration-panel/runs/RUN-010.yaml`

## Interpretation

- `RUN-010` addressed the main confound left by `RES-001`: it restricted the panel to one assay class, used raw ProteinGym scores instead of assay-wise z-scores, learned assay-specific nonlinear observation models with `mavenn`, and included real multi-mutant measurements.
- That change worked at the measurement layer. The assay-specific `mavenn` models fit the raw assays well, with mean test Spearman `0.894` and mean test NRMSE `0.491`, so the coordinate transform itself was not the dominant failure mode in this run.
- The downstream shared Adapt-Env calibration remained weak. The better of the two shared branches, `shared_raw`, reached only `0.232` single-mutant holdout Spearman, `1.075` single-mutant holdout NRMSE, `-0.017` double-mutant holdout Spearman, and retained `epistasis_strength = 0.0`.
- The `shared_latent` branch was worse on most ranking metrics and collapsed to a near-peak reference solution with `functional_reference_fraction_of_peak = 0.897`, `empirical_pairwise_strength = 0.0`, and `epistasis_strength = 0.0`.
- The per-assay latent fits also remained poor. Mean single-mutant holdout Spearman was `0.136`, mean double-mutant holdout Spearman was `0.116`, and all six per-assay fits still set `epistasis_strength = 0.0`. This means the failure cannot be assigned mainly to cross-assay pooling or to the previous z-score and wild-type-zeroing convention.

## Effect on hypothesis

- This result further weakens `HYP-001`.
- `RES-001` could still be challenged on the grounds that assay-coordinate choices and mixed assay families confounded the test. `RUN-010` removed much of that confound and still failed to show convincing recovery of local single- and double-mutant statistics.
- The result does not refute `HYP-001` outright because it is still limited to six short stability assays and the current calibration objective, but it materially reduces the plausibility that coordinate-system mismatch alone explains the earlier negative evidence.

## Limitations

- The panel is homogeneous and scientifically cleaner than `RUN-006`, but it is still small: `6` short assays from one study and one experimental readout.
- The experiment did not rerun the original z-scored baseline on the exact same six-assay subset, so comparisons to `RES-001` are mechanistic and qualitative rather than perfectly paired.
- The `mavenn` observation models are learned approximations to the assay measurement process rather than ground-truth latent phenotypes.
- Several per-assay ranking metrics are weak or null, and one assay produced null holdout Spearman values in the stored summary, which limits how much can be inferred from the per-assay ceiling comparison.

## Downstream use

- Use `RES-003` as the current best-controlled test of the claim that coordinate-aware calibration rescues `HYP-001`; it does not.
- Use this result together with `RES-002` when interpreting `HYP-001`: the model family appears locally expressive, but the present empirical calibration and reconstruction strategy remains inadequate even after the main coordinate confound is reduced.
