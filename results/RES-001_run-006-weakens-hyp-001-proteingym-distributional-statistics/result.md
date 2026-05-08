# RES-001: RUN-006 weakens HYP-001 through incomplete recovery of ProteinGym distributional statistics

## Summary

RUN-006 produced a balanced eight-assay ProteinGym panel with real MMseqs alignments, but the best shared parameter regime collapsed epistasis, empirical pairwise structure, and noise to zero and showed weak holdout performance. The synthetic landscapes reproduced some quantiles and conservation trends in a subset of assays, but broadly under-represented beneficial tails and over-produced negative skew, so the evidence weakens rather than supports HYP-001.

## Generated from

- Analyses: ANA-001

## Relevant hypotheses

- Supports: None
- Weakens: HYP-001
- Refutes: None

## Evidence

- `analyses/2026-05-08_ANA-001_run-006-proteingym-distributional-realism-review/tables/run-006_key_metrics.md`
- `data/processed/proteingym-dms-distributional-realism-panel/RUN-006/summary.json`
- `data/processed/proteingym-dms-distributional-realism-panel/RUN-006/selected_panel.csv`
- `experiments/2026-05-07_EXP-001_proteingym-dms-distributional-realism-panel/config.yaml`
- `experiments/2026-05-07_EXP-001_proteingym-dms-distributional-realism-panel/runs/RUN-006.yaml`

## Interpretation

- The run answers the right external-realism question for `HYP-001`: one shared Adapt-Env regime was fit against real ProteinGym single-mutant assays built from real MMseqs alignments rather than internal synthetic alignments.
- The best-fit solution did not preserve the intended richer model structure. It drove `epistasis_strength`, `empirical_pairwise_strength`, and `noise_amplitude` to zero, so the fit is effectively a deterministic low-dimensional functional model rather than the hybrid model described in the hypothesis.
- Predictive adequacy remained weak. Holdout Spearman was only `0.103`, holdout NRMSE was `0.998`, and the functional KS distance was `0.537`, which is too poor to treat the shared regime as a convincing recovery of empirical local landscape statistics.
- The synthetic landscapes captured some coarse local features in a subset of assays, but the dominant error pattern was systematic: beneficial upper tails were compressed, skewness was too negative, and one assay (`OXDA_RHOTO_Vanella_2023_activity`) even reversed the sign of the conservation-sensitivity relationship.

## Effect on hypothesis

- This result weakens `HYP-001`.
- The evidence is not strong enough to refute the hypothesis outright because the current test only covers single-mutant z-scored assay summaries and a limited panel of 8 assays.
- It does show that the current calibration setup does not yet support the stronger claim that one shared synthetic regime can recover empirical local summary statistics across diverse ProteinGym assays without assay-specific fitting.

## Limitations

- The analysis is restricted to single-mutant assays and does not test double-mutant or higher-order epistasis recovery.
- Assay scores were standardized within assay and the wild type was anchored at zero, so the comparison is about distributional shape and relative structure, not absolute assay scale.
- The conclusion is specific to the current calibration grid and objective. A different calibration target could behave differently, but that has not been shown here.
- Durable run outputs are present under `data/processed`, but they are not yet DVC-tracked in the current repository state.

## Downstream use

- Use this result to constrain claims about external realism: the current shared regime should not be described as validated against ProteinGym-like local landscape statistics.
- Use `RES-001` as the current record of the hypothesis impact for `RUN-006` when discussing `HYP-001` in later analyses or manuscript-facing claims.
