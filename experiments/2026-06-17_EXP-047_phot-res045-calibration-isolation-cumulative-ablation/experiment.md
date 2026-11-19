# EXP-047: PHOT RES-045 calibration isolation cumulative ablation

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Do the RES-045 calibration-isolation fixes, applied cumulatively, separate PHOT raw-scale, geometry, texture, sequence-map, and search failures enough to recover a navigable but rugged latent graphlet package?

## Pre-experiment prediction

Anchoring raw latent scale and decoupling scale-invariant geometry should reduce the RES-045 distance inflation. If the remaining failure is texture, the k-trait and hybrid residual stages should raise raw local maxima and Walsh high-order fraction while preserving accessibility. If the fast search was limiting, the expanded-search stage should improve selected candidates beyond the earlier cumulative stages.

## Rationale

RES-045 showed that a two-trait soft-min latent moved raw accessibility for the first time, but it left raw scale hundreds-fold too large, produced too few local maxima, had a weak Walsh high-order tail, fit empirical node ranks only modestly, and returned identical later-stage selected rows. The scratch isolation plan argues that these are separable defects rather than one undifferentiated model failure.

This experiment implements that plan as a single PHOT-only diagnostic while keeping the frozen `RUN-077` PHOT source, the recent fast-fit reporting path, and the `RUN-098` full-ABC comparator fixed.

## Experimental Design

- Run cumulative stages in the plan's recommended order:
  - `A`: z-score anchor the soft-min generator output so raw scale is controlled by one output-scale knob.
  - `A+F`: keep `A` and decouple the objective into scale-invariant geometry plus a separate raw single-mutant-scale term.
  - `A+F+B`: run a model-free k-trait additive soft-min capability probe on empirical graphlet node values; this stage intentionally skips duplicate sequence-conditioned fitting.
  - `A+F+B+C`: re-enable sequence-conditioned fitting and sweep k-trait dimensionality and auxiliary Walsh texture.
  - `A+F+B+C+D`: add a small sparse structured residual to the navigable k-trait soft-min core and restrict the observation readout to a two-parameter monotone sigmoid.
  - `A+F+B+C+D+E`: increase random and local-refinement search on the hybrid form to test whether recent failures were partly search-limited.
- Select top sequence-conditioned candidates within each fitted stage, not globally across stages.
- Report raw and observable metrics together so readout-scale improvement cannot be mistaken for latent-geometry calibration.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
- `scratch/RES-045_epistasis_calibration_isolation_plan.md`
- `scripts/phot_trait_space_bounded_gate_fast_fit.py`
- `src/project_adapt_env/coupled_landscape.py`

## Configuration

- `config.yaml` defines the frozen PHOT source, RUN-098 comparator, graphlet summaries, cumulative stage overrides, capability probe, fast-fit budgets, scale-invariant objective, and quick smoke-test override.

## Execution

- Run through `labproj submit EXP-047 --target lab-slurm --execute`.
- `run.sh` pulls RUN-077/RUN-098 DVC targets if missing, then runs `scripts/phot_trait_space_bounded_gate_fast_fit.py`.

## Expected Outputs

- summary.json
- preflight.json
- bounded_gate_candidates.csv
- selected_bounded_gate_candidates.csv
- bounded_gate_distances.csv
- bounded_gate_equivalence.csv
- graphlet_metrics.csv
- graphlet_summary.csv
- softmin_capability_probe.csv
- walsh_terms.csv
- gft_coefficients.csv
- abc_comparator_selected_family_distances.csv

## Analysis Plan

- Create an ANA record after execution comparing the cumulative stages against RES-045/RUN-119 and RUN-098.
- Primary interpretation hinges on whether raw single-mutant scale becomes controlled without losing raw accessibility, raw local maxima, Walsh high-order fraction, and held-out graphlet distance.
- Use the model-free probe only to judge whether k-trait soft-min surfaces have obvious direct capacity on PHOT graphlets. It does not prove the sequence-conditioned generator can find those surfaces.
- Keep any positive interpretation PHOT-local until the same calibration behavior is shown on another assay.

## Completion Criteria

- [x] Config committed
- [x] Run script committed
- [x] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment Notes

- This experiment follows `scratch/RES-045_epistasis_calibration_isolation_plan.md` as a diagnostic plan, not as prior evidence.
