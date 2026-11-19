# RES-031: RUN-079 no-module coupling improves PHOT graphlet distance but local modules fail the locality test

## Summary

`RUN-079` is a valid deterministic PHOT-only frozen-base grid test of
`HYP-008`. It shows that a correlated nonlinear no-module operator can improve
held-out PHOT graphlet distance over the frozen current model while preserving
base correlation. However, the local empirical plaquette-module coupled
candidates do not beat no-module, random, shuffled, or uncoupled controls under
calibration-selected held-out evaluation, and module-energy diagnostics show
diffuse rather than concentrated innovation. The result weakens the current
local-module coupled-landscape mechanism while leaving broader nonlinear
trait-operator ideas unvalidated rather than refuted.

## Generated from

- Analyses: `ANA-031`
- Experiments: `EXP-031`
- Runs: `RUN-079`

## Relevant hypotheses

- Supports: None
- Weakens: `HYP-008`
- Refutes: None

## Evidence

- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/metrics.json`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_preflight_and_output_audit.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_empirical_target_graphlet_metrics.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_calibration_selected_distances.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_calibration_selected_by_module_source.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_static_control_distances.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_test_oracle_distances_capacity_only.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_relevant_test_graphlet_summaries.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_base_correlation_heldout.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_module_energy_locality_heldout.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_module_catalog.md`
- `analyses/2026-05-16_ANA-031_run-079-phot-frozen-base-coupled-landscape-grid-review/tables/run-079_hyp008_decision_table.md`
- `results/RES-031_run-079-no-module-coupling-improves-phot-graphlet-distance-but-local-modules-fail-locality/metrics.json`
- `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079.dvc`
- `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079/summary.json`
- `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079/model_distances.csv`
- `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079/graphlet_summary.csv`
- `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079/module_energy_summary.csv`
- `data/processed/phot-frozen-base-coupled-landscape-deterministic-grid/RUN-079/base_correlation.csv`

## Interpretation

The procedural execution is sound. `RUN-079` completed with scheduler exit code
`0`, used the `RUN-077` PHOT frozen base, did not rerun ABC, and produced the
planned WHT, GFT, graphlet-topology, model-distance, base-correlation, and
module-locality outputs. The collected directory is DVC-tracked and pushed.

The valid held-out grid-search comparison is calibration-selected, not
test-oracle. Under that comparison, the best no-module coupled operator beats
the frozen current model on held-out PHOT graphlets (`0.3137` versus `0.4016`
distance) and preserves base structure well (`0.982` phenotype correlation and
`0.987` edge-gradient correlation). This is a useful technical signal: the
frozen-base coupled machinery can move graphlet statistics without destroying
the base landscape.

The same result is negative for the local-module mechanism that motivated
`EXP-031`. The calibration-selected empirical plaquette-module coupled
candidate has held-out distance `0.3982`, worse than the no-module operator
(`0.3137`), worse than random modules (`0.3546`), and comparable to shuffled
modules (`0.3946`) and the frozen current model (`0.4016`). The globally best
calibration-selected parameter is an uncoupled local-module setting, but its
held-out distance is only `0.3951`, again barely separated from frozen current.

The locality metric also fails. The calibration-selected empirical-module
candidate has plaquette effective module count `9.33` out of `12` and
normalized entropy `0.942`; the empirical readout residual has effective count
`9.72` and entropy `0.952`. That is diffuse innovation across selected modules,
not concentrated local ruggedness. A random-module control is more concentrated
than the empirical modules without being biologically meaningful, so locality
cannot be inferred from this module-selection procedure.

Relative to `RES-030`, the graphlet distance can be improved, but the PHOT
spectral problem is not rescued as a local-module mechanism. The empirical test
WHT high-order fraction is `0.3788`; the frozen current model is `0.4984`, the
calibration-selected empirical-module coupled candidate is `0.4966`, and the
calibration-selected no-module candidate is `0.4713`. The best empirical
readout residual diagnostic has the closest distance (`0.2708`), but it is an
empirical subtraction and not a procedural generator.

## Effect on hypothesis

`RES-031` weakens `HYP-008` in its current operational form. The run does not
support the specific claim that calibration-derived local empirical modules
provide the missing non-diffuse ruggedness needed to repair PHOT graphlet
spectra/topology.

The result does not fully refute `HYP-008` because the tested implementation is
limited: PHOT-only, deterministic, no ABC refit, and one crude module-selection
procedure. It leaves open the broader possibility that another nonlinear
trait-operator formulation could be useful, but that claim remains unvalidated
and now carries a stronger control burden.

## Limitations

- PHOT is the only empirical system tested.
- The run is a deterministic grid over a frozen `RES-030` base, not a joint
  probabilistic fit.
- The empirical target is still the processed `latent_phi` graphlet target, not
  direct exhaustive empirical measurement.
- Test-oracle rankings are capacity diagnostics only and should not be used as
  held-out evidence.
- The selected modules share site structure and may not represent independent
  molecular patches; the module-energy diagnostics show that this matters.

## Downstream use

- Treat `RUN-079` as evidence that the coupled-landscape machinery executes and
  can preserve base correlation while changing graphlet statistics.
- Do not use the current empirical plaquette-module selection as validated
  evidence for localized ruggedness.
- Keep `RES-031` as a negative control record for future `HYP-008` variants:
  future models must beat no-module, random/shuffled-module, uncoupled-module,
  empirical-readout, and frozen-current controls under calibration-selected
  held-out evaluation and must show genuinely concentrated innovation energy.
