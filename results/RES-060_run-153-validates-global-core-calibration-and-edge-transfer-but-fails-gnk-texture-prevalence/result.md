# RES-060: RUN-153 validates global core calibration and edge transfer but fails GNK texture prevalence

## Summary

RUN-153 formally passed `6/9` required EXP-062 checks and issued
`CALIBRATED_SIGMOID_TRANSMISSION_FAILED`. The globally frozen mutation-unit
calibration passed all 60 held-out core cells and all 360 complete-h1 edge
cells, and GNK-specific texture remained aligned and dose-directional. However,
only `24/360` actual-core cells and `0/360` frozen-core cells met the
preregistered texture practical-prevalence floor; exact per-landscape
normalization did not rescue it.

The recorded directionality failure is a floating-point equality artifact and
becomes a pass under `1e-12` tolerance, changing the independent audit to `7/9`
without changing the overall decision. This weakens the full observable
protein-scale form of HYP-005 while preserving narrow calibration, edge, dose,
and length-stability evidence.

## Generated from

- Analysis: ANA-060
- Experiment: EXP-062
- Run: RUN-153

## Relevant hypotheses

- Supports: None
- Weakens: HYP-005
- Refutes: None

## Evidence

- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/metrics.json`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_acceptance_audit.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_collection_audit.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_core_calibration_by_length.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_heldout_core_by_length.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_transmission_by_length.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_directionality_by_length.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_dose_curve.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_scale_law_summary.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_ablation_by_length.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_protocol_output_audit.csv`
- `analyses/2026-07-23_ANA-060_run-153-dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-review/tables/run-153_source_hash_audit.csv`
- `data/processed/dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-validation/RUN-153.dvc`

## Interpretation

The additive-core calibration works on the tested panel. The core-only phase
estimated one global mutation scale, `sigma_A = 2.744864`, with coefficient of
variation `0.02274`, then selected `beta=4` and
`mu=logit(0.99)=4.595120` without access to GNK outcomes. The frozen choice
passed every held-out core cell in every length bin. Global realized beta ranged
only from `3.8487` to `4.2042`, explaining why exact per-landscape normalization
adds little.

Complete-h1 edge transfer is also repaired. At `a=1`, all `360/360` cells pass;
median practical fraction is `0.8828`, median correlation `0.8717`, and median
residual RMS `1.5342`. This is a large improvement over the sparse one-step
transfer in RES-059.

GNK-specific two-site texture is transmitted, but not with the preregistered
prevalence. Actual-core texture has median RMS `0.4175`, correlation `0.7322`,
projection `0.9182`, positive projection and permutation `p=0.002` in every
cell, and a crossed-bootstrap projection lower bound of `0.8731`. Its minimum
actual/frozen RMS ratio is `1.263`. Nevertheless, median practical fraction is
only `0.1927`, and just `24/360` cells meet the required `0.25` fraction. The
failure is therefore insufficient prevalence above the heuristic `0.10`
per-contrast floor, not disappearance or reversal of the GNK texture.

The frozen-core positive control also fails practical prevalence: median RMS is
`0.1317`, practical fraction `0.0807`, and `0/360` cells meet the floor, despite
positive projection and permutation support throughout. Because it freezes the
reference at the selected occupancy of `0.99`, it is intrinsically saturated
and behaves more like an upper-shoulder saturation control. Its formal failure
is retained, but RES-059's occupancy-0.5 frozen-core result remains the relevant
positive control for sigmoid texture transmission.

Exact normalization does not rescue the primary result: its median texture RMS
is `0.4153`, practical fraction `0.1901`, and only `23/360` cells meet the
prevalence floor. Per-landscape normalization is therefore not required by
these data.

The dose magnitude is monotone. Global median texture RMS rises
`0.1085 -> 0.2163 -> 0.4175 -> 0.7578` over `a=0.25,0.5,1,2`, all cell RMS
Spearmans are one, and all 12 protein-median curves pass. The runner's global
directionality failure occurs because mathematical `rho=0.80` values are stored
as `0.7999999999999999`. A `1e-12` comparison tolerance changes cell passes from
`305/360` to `354/360`, with every bin above the registered `0.80` fraction.
The formal output remains part of the record, but it is not evidence of
biological nonmonotonicity.

Length stability passes in both arms: the global slope is `-0.10221`
(`[-0.11388,-0.09012]`) and the exact slope is `-0.11358`
(`[-0.12531,-0.10255]`). This is bounded stability with a modest negative
trend, not strict invariance.

## Effect on hypothesis

RES-060 weakens the full observable protein-scale form of HYP-005 because the
registered GNK-specific texture-prevalence and frozen-control requirements fail.
It does not erase the narrower positive evidence: one global mutation-unit
calibration generalizes on the held-out panel, one-step edge transfer is broad,
the intended dose response is monotone under a numerically sound comparison,
and texture RMS remains within the registered length-stability bounds.

The result does not refute all additive-core plus GNK models, require
per-landscape calibration, or show that GNK texture is absent. It also does not
validate the current public generator or service.

## Limitations

- The 12 confirmation proteins are outcome-blind and held out, but form a fixed
  curated panel rather than a probability sample.
- The `0.10` practical floor and `0.25` prevalence criterion are noiseless
  engineering heuristics, not empirical detection limits.
- The frozen positive control is confounded by its reference occupancy of
  `0.99`.
- The formal directionality failure is a floating-point boundary artifact.
- Checkpoint files contain non-standard JSON `NaN` literals.
- Steepness-only, threshold-only, and raw-margin ablations report edges rather
  than texture, so they do not isolate a sufficient texture-repair component.
- No empirical outcomes, observation noise, accessibility, adaptive walks,
  active learning, or optimizer performance were tested.

## Downstream use

- Cite RES-060 for held-out support of the global additive-core calibration and
  complete-h1 edge transfer.
- Cite RES-060 for failure of the preregistered GNK texture-prevalence package
  and for the finding that exact per-landscape normalization does not rescue it.
- Do not cite the formal directionality flag as evidence of nonmonotonicity;
  cite the tolerance-aware audit.
- Do not cite RUN-153 as validation of the public AdaptEnv generator, empirical
  realism, or an optimizer benchmark.
- Keep the Stage 4C/4D/4E hard stop in force.
