# EXP-011: ProteinGym Single-Assay Stability-Readout Control SMC-ABC Diagnostic for PHOT_CHLRE

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001, HYP-007

## Question

Does the simpler stability-targeted single-latent calibration used in RES-007 recover the PHOT_CHLRE activity assay better than the recent activity-specific formulations, or are the PHOT failures robust even under that simpler control?

## Pre-experiment prediction

If the negative PHOT results mainly reflect over-parameterization or confounding in the newer activity-specific formulations, then reverting to the simpler RES-007 stability-targeted calibration may improve holdout ranking and epistasis metrics on PHOT_CHLRE despite its causal mismatch. If it does not, that would strengthen the case that the assay-model mismatch is more fundamental.

## Rationale

RES-007 showed that a stability-targeted readout materially improved recovery on a matched unfolding assay, while RES-008 and RES-009 showed that two richer PHOT_CHLRE activity formulations still did not rescue the same activity landscape. Before attributing that failure entirely to missing latent causal structure, this control tests whether the simpler single-latent stability-targeted fitter from RES-007 can recover PHOT_CHLRE better than the newer activity-oriented versions. A positive control result would suggest that the current multi-latent activity parameterization is harder to fit than necessary; a negative control result would imply that the PHOT assay itself is poorly matched to the stability-only causal chain.

## Experimental design

- Template: generic
- Use the exact single-assay `PHOT_CHLRE_Chen_2023` activity panel from
  `EXP-008` and `EXP-010`.
- Revert the generator and calibration semantics to the successful
  `RES-007` stability-targeted single-latent formulation:
  - default collapsed generator from `EXP-007`
  - `synthetic_readout_mode: stability_margin`
  - `empirical_pairwise_target: stability`
- Keep the paired deterministic controls from `EXP-007`:
  - `baseline_shared_raw`
  - `predictive_richpair_shared_raw`
  - `baseline_shared_stability_readout`
  - `predictive_richpair_shared_stability_readout`
- Run one Bayesian shared fit under the same stability-readout calibration path
  used in `RUN-024`, but on the `PHOT_CHLRE` activity assay.
- Compare holdout ranking, double-mutant recovery, epistasis prediction, KS, and
  reference-to-peak behavior against both the within-run raw controls and the
  historical `EXP-008` / `EXP-010` PHOT results.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) paired branch comparisons between raw-readout and stability-readout controls on PHOT_CHLRE; 2) assay-specific mavenn diagnostics; 3) SMC-ABC posterior particles and round diagnostics under the RES-007-style stability-readout path; 4) synthetic-truth recovery metrics for that path; and 5) a structured summary.json

## Analysis plan

- Create an ANA record after execution.

## Completion criteria

- [ ] Config committed
- [ ] Run script committed
- [ ] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
