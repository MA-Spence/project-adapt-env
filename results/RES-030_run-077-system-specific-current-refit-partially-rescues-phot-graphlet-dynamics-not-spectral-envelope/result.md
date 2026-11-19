# RES-030: RUN-077 system-specific current refit partially rescues PHOT graphlet dynamics but not the spectral envelope

## Summary

RUN-077 fixes the RUN-073 pooling and stale-cache comparison problem by fitting PHOT and SPTN1 as independent current-model systems and emitting dense graphlet, Walsh, and GFT outputs. The PHOT current model is now the closest held-out candidate for the joint graphlet package and beats GNK and residual-only controls, so the strongest RES-028 evolutionary-dynamics failure is not reproduced under current system-specific refitting. However, PHOT still violates the preregistered spectral envelope, high-order signal is not removed by nonlinear readout correction, and the current model is only marginally better than smooth core. SPTN1 is a smooth 2D negative control and cannot test high-order spectral claims. The result weakens the strong HYP-004 benchmark-adequacy claim while motivating HYP-008 local non-diffuse ruggedness.

## Generated from

- Analyses: `ANA-030`
- Experiments: `EXP-030`
- Runs: `RUN-077`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-004`
- Refutes: None
- Motivates: `HYP-008`

## Evidence

- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/metrics.json`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_preflight_and_output_audit.md`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_graphlet_catalog_summary.md`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_fit_summary.md`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_branch_validation.md`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_empirical_primary_graphlet_metrics.md`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_test_model_graphlet_metrics.md`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_model_distances.md`
- `analyses/2026-05-16_ANA-030_run-077-system-specific-current-model-graphlet-falsifier-review/tables/run-077_falsifier_summary.md`
- `results/RES-030_run-077-system-specific-current-refit-partially-rescues-phot-graphlet-dynamics-not-spectral-envelope/metrics.json`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077/summary.json`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077/falsification_tests.json`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077/system_graphlet_metrics.csv`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077/system_model_distances.csv`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077/system_walsh_terms.csv`
- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077/system_gft_coefficients.csv`

## Interpretation

- `RUN-077` corrects the main procedural defect in `RUN-073`: PHOT and SPTN1
  were fit as two independent current-model systems, the legacy-parameter guard
  passed, and the run emitted dense graphlet, WHT, and GFT statistics.
- The PHOT held-out graphlets partially reverse the strongest negative
  `RES-028` dynamics conclusion. The fitted current model is the closest test
  candidate for the joint graphlet package (`0.4016` distance), ahead of smooth
  core (`0.4119`), posterior mean (`0.4408`), GNK (`0.7594`), and residual-only
  diagnostic (`2.6756`).
- The PHOT spectral-envelope falsifier still fires. Held-out empirical PHOT
  graphlets have WHT high-order fraction `0.3788`, GFT low-energy fraction
  `0.6212`, and GFT high-energy fraction `0.1552`; this exceeds the
  preregistered high-order-tail ceiling and is not removed by nonlinear
  global-readout correction.
- The PHOT current-vs-smooth-core separation is too small to claim a clean
  mechanistic validation of the residual epistasis layer. It is evidence that
  system-specific current refitting improves the graphlet package, not evidence
  that the current diffuse residual structure is the right scientific model.
- SPTN1 functions mainly as a smooth negative control in this run. The retained
  dense graphlets are only 2D, with SE/RSE both `0`, one local maximum, and no
  testable high-order WHT/GFT tail. The runner flags wrong evolutionary
  dynamics because the posterior-mean branch beats the best particle on the
  held-out test distance, but this is branch-selection instability in a trivial
  graphlet regime rather than a meaningful high-order ruggedness falsifier.

## Effect on hypothesis

`RES-030` weakens `HYP-004` in its strong benchmark-adequacy form. The current
model is not simply failing the full graphlet package under current,
system-specific refitting, but it still does not reproduce the empirical PHOT
spectral envelope and is barely separated from smooth core. The appropriate
scientific conclusion is mixed rather than rescued.

The result motivates `HYP-008`. The residual/current model can now recover more
of the PHOT topological package than stale or pooled comparisons suggested, but
the remaining spectral mismatch and weak smooth-core separation are consistent
with the possibility that empirical ruggedness is local and non-diffuse rather
than a broadly distributed residual field.

## Limitations

- SPTN1 graphlets are only 2D, so they cannot evaluate high-order WHT/GFT
  falsifiers.
- This run is not a residual-scale grid and should not replace `RES-027` or the
  scale-control component of `RES-028`.
- The PHOT empirical target still depends on the fitted MAVE-NN/latent assay
  layer; it is not direct exhaustive empirical measurement of all sampled
  graphlets.
- The PHOT current model only marginally beats smooth core, so residual
  mechanism claims should remain provisional.

## Downstream use

- Use `RUN-077`/`RES-030` as the current system-specific PHOT/SPTN1 graphlet
  retest of the `RES-028` falsifiers.
- Do not treat stale pooled `RUN-073` outputs as valid graphlet evidence.
- Keep the PHOT wrong-spectral-envelope result active when designing the
  `HYP-008` local nonlinear trait-operator experiments.
