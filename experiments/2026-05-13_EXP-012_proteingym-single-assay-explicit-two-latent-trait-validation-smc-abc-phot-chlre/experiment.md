# EXP-012: ProteinGym Single-Assay Explicit Two-Latent-Trait Validation-Objective SMC-ABC Diagnostic for PHOT_CHLRE

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001, HYP-007

## Question

Does an explicit two-latent-trait model, with stability plus one fitted readout trait and only empirical pairwise epistasis routed into that readout trait, recover the `PHOT_CHLRE_Chen_2023` activity assay better than the earlier PHOT implementations once the Bayesian fit is optimized directly on held-out predictive metrics rather than only on summary-vector matching?

## Pre-experiment prediction

If the negative `PHOT_CHLRE` results in `RES-008` and `RES-009` were driven partly by latent-parameter confounding and objective mismatch, then this cleaner explicit two-trait formulation should improve the full predictive package, especially double-mutant holdout, epistasis-prediction Spearman, and reference-to-peak behavior, relative to those earlier PHOT runs.

## Rationale

The later code review after `RES-008` and `RES-009` found that those experiments did not yet fit an explicit second latent trait block. They changed readout semantics, but still relied on the legacy single functional block and also let multiple epistasis levers compete while the Bayesian target remained a bootstrap summary-vector distance rather than the held-out predictive metrics used to judge success. This experiment is the direct replacement test: it uses one explicit readout trait block, keeps stability as the other latent trait, disables generic score-space epistasis, routes only empirical pairwise structure into the readout trait, and makes the SMC distance the same held-out predictive objective that the analysis will later inspect.

## Experimental design

- Template: generic
- Use the exact same single-assay `PHOT_CHLRE_Chen_2023` activity panel used in `EXP-008`, `EXP-010`, and `EXP-011`.
- Replace the legacy single functional block with one explicit latent readout block:
  - built-in latent stability layer
  - one explicit named latent trait block, `readout`
  - public scalar fitness formed as `stability_gate * trait:readout:capacity`
- Remove the confounding generic epistasis lever from the fit:
  - keep `epistasis_strength = 0`
  - keep `n_higher_order_epistatic_terms = 0`
  - fit only `empirical_pairwise_strength`, routed to `trait:readout`
- Drop the earlier free FGM-complexity search:
  - fix the explicit readout block to `1` latent dimension
  - fit only its width through `functional_sigma_base`, which now aliases to the primary explicit trait block
  - fix `peak_distance_from_consensus` instead of fitting it in this run
- Keep paired deterministic controls under the same two-trait generator:
  - a no-pairwise baseline
  - a predictive rich-pair branch
- Run one Bayesian SMC-ABC fit on the same assay, but change the particle distance to the direct validation objective built from held-out single-mutant, held-out double-mutant, epistasis-prediction, KS, and reference-geometry terms.
- Compare the resulting predictive package against `RES-008`, `RES-009`, and the simpler control from `EXP-011`.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) paired deterministic branch comparisons under the explicit two-trait generator; 2) assay-specific `mavenn` diagnostics for `PHOT_CHLRE`; 3) posterior particles, round summaries, and branch validations from the validation-objective SMC run; and 4) a structured `summary.json`

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
