# ANA-013: RUN-039 PHOT stability-readout validation-objective control review

## Purpose

Assess whether the EXP-013 PHOT stability-targeted control retains its predictive advantage when fitted with the held-out validation objective used in EXP-012, and record the implication for HYP-007.

## Linked experiments/runs

- Experiments: EXP-013
- Runs: RUN-039

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `tables/run-039_key_metrics.md`

## Inputs Reviewed

- `experiments/2026-05-14_EXP-013_proteingym-phot-stability-readout-validation-objective-control/metadata.yaml`
- `experiments/2026-05-14_EXP-013_proteingym-phot-stability-readout-validation-objective-control/config.yaml`
- `experiments/2026-05-14_EXP-013_proteingym-phot-stability-readout-validation-objective-control/runs/RUN-039.yaml`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/summary.json`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/mavenn_assay_metrics.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/branch_validations.csv`
- `data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039/posterior_parameter_summary.csv`
- `results/RES-010_run-036-explicit-two-latent-trait-validation-objective-does-not-rescue-phot-chlre-recovery/metrics.json`
- `results/RES-011_run-038-stability-readout-control-improves-phot-chlre-recovery-but-does-not-fully-rescue-the-assay/metrics.json`

## Analysis Performed

`RUN-039` was checked as the live `EXP-013` execution. The run completed with
Slurm job `172`, scheduler exit code `0`, and output under
`data/processed/proteingym-phot-stability-readout-validation-objective-control/RUN-039`.
No existing ANA or RES record linked to `RUN-039` was present before this
analysis.

The analysis compared the within-run deterministic raw and stability-readout
branches, the `validation_objective` SMC best particle and posterior mean, and
the two immediately motivating PHOT records: `RES-010` for the explicit
two-latent validation-objective run and `RES-011` for the earlier
stability-readout summary-vector control.

## Observations

- The observation layer remained comparable to the recent PHOT runs. The
  assay-specific MAVE-NN model for `PHOT_CHLRE_Chen_2023` reached test
  Spearman `0.689` and test NRMSE `0.966` on `310` held-out variants.
- The deterministic stability-readout branches again outperformed the raw
  branches on single-mutant holdout ranking and functional KS, but they still
  showed the peak-geometry problem and weak epistasis prediction. The stronger
  deterministic stability branch had single-mutant holdout Spearman `0.505`,
  double-mutant holdout Spearman `0.703`, epistasis-prediction Spearman
  `0.032`, KS `0.190`, and reference fraction of peak `0.987`.
- The Bayesian `validation_objective` stability-targeted fit materially changed
  the picture. The best particle retained single-mutant holdout Spearman
  `0.506`, improved double-mutant holdout Spearman to `0.783`, improved
  epistasis-prediction Spearman to `0.634`, kept KS low at `0.176`, and moved
  the reference away from the peak (`fraction_of_peak 0.323`, distance `45`).
- Relative to the explicit two-latent validation-objective run from `RUN-036`,
  the `RUN-039` Bayesian best particle improved single-mutant holdout Spearman
  by `0.308`, double-mutant holdout Spearman by `0.531`, epistasis-prediction
  Spearman by `0.376`, and KS by `0.320` in the favorable direction, while
  lowering reference fraction of peak by `0.663`.
- Relative to the earlier stability-readout summary-vector control from
  `RUN-038`, the `RUN-039` validation-objective fit preserved single-mutant
  ranking, improved double-mutant holdout Spearman by `0.062`, improved
  epistasis-prediction Spearman by `0.568`, and greatly improved reference
  geometry, at the cost of a small KS regression (`0.166` to `0.176`).
- `validation_objective` mode again disables the matched synthetic-truth
  recovery scaffold, so the run cannot independently quantify inverse-problem
  identifiability under synthetic truth.

## Interpretation

`RUN-039` answers the main control question from `EXP-013`: the
stability-targeted PHOT control does retain and strengthen its predictive
advantage when fit under the same held-out validation objective used by the
explicit two-latent run. The result is not merely an artifact of the older
summary-vector objective used in `RUN-038`.

This is scientifically important because it separates two explanations for the
recent PHOT negatives. The poor `RUN-036` two-latent result cannot be explained
only by the held-out validation objective being too difficult; the same
objective, applied to a stability-targeted causal mapping, produces the
strongest PHOT package recorded so far. Conversely, the result does not prove
that PHOT fluorescence is purely a stability assay. It shows that the current
implemented multi-latent/readout formulations are still structurally mismatched
to this scalar activity readout.

## Hypothesis Update

`RUN-039` weakens the current `HYP-007` implementation on `PHOT_CHLRE_Chen_2023`.
The richer explicit two-trait PHOT model no longer has the objective-mismatch
defense: a simpler stability-targeted mapping fits the same assay substantially
better under the same validation-objective family. The result does not directly
update `HYP-001`, because this is still a single-assay diagnostic rather than a
shared-regime multi-assay realism test.
