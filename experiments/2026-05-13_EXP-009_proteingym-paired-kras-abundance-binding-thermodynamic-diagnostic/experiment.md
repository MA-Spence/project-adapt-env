# EXP-009: ProteinGym Paired KRAS Abundance-Binding Thermodynamic Diagnostic

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-007

## Question

Does replacing collapsed binding fitness with a paired abundance-plus-binding
thermodynamic readout improve empirical recovery on the matched KRAS ProteinGym
assays `RASK_HUMAN_Weng_2022_abundance` and
`RASK_HUMAN_Weng_2022_binding-DARPin_K55`?

## Pre-experiment prediction

If `RES-008` failed mainly because `PHOT_CHLRE` fluorescence is not a clean
two-latent assay and because one readout leaves stability and function weakly
identifiable, then a paired KRAS abundance-plus-binding fit with a three-state
folded/unfolded/bound readout should improve binding holdout ranking,
double-mutant recovery, epistasis prediction, and reference-to-peak behavior
relative to both binding-only raw controls and paired raw-readout controls.

## Rationale

`RES-008` weakens the specific `PHOT_CHLRE` implementation of `HYP-007`, but it
does not isolate the more plausible technical failure mode. The assay reports
cellular fluorescence from a flavin-binding fluorescent protein, so the observed
score is likely influenced by more than folding stability and one generic
function latent. Published FbFP characterization and the CreiLOV DMS paper both
point to fluorescence depending on folded fraction, FMN occupancy, and
photophysical brightness, while the experiment exposed only one measured
phenotype.

Primary literature suggests a better test. Otwinowski's GB1 reanalysis showed
that a three-state folding-plus-binding model explains binding-assay epistasis
better than simpler collapsed treatments, and Weng et al. showed on KRAS that
paired abundance and binding measurements plus many double-mutant backgrounds
support causal inference of folding and binding free-energy effects. The local
ProteinGym panel contains exactly such a matched pair for KRAS, with
`3,066/22,946` single/multiple mutants for abundance and `3,084/21,789` for
DARPin K55 binding.

This experiment therefore shifts from a single fluorescence assay to one protein
with orthogonal abundance and binding readouts, and from a minimal
stability-times-function map to an explicit thermodynamic binding readout. The
goal is not to broaden scope beyond `HYP-007`, but to test it under a causal
chain that is closer to what the literature says these assays actually measure.

## Experimental design

- Template: generic
- Prepare the matched KRAS ProteinGym assays `RASK_HUMAN_Weng_2022_abundance`
  and `RASK_HUMAN_Weng_2022_binding-DARPin_K55`.
- Treat the abundance assay as the stability-facing readout and the DARPin K55
  assay as the binding-facing readout.
- Compare four branches:
  - `binding_only_raw`
  - `binding_only_biophysical_binding`
  - `paired_abundance_binding_raw`
  - `paired_abundance_binding_biophysical_binding`
- Use the same predictive deterministic calibration machinery across branches,
  with pairwise holdout, iterative functional fitting, and explicit
  reference-to-peak penalties.
- In the thermodynamic branches, fit a folded/unfolded/bound readout using
  `synthetic_readout_mode: stability_binding` and route empirical pairwise
  corrections through the function-facing latent layer.
- Evaluate abundance KS, binding single-mutant holdout, double-mutant holdout,
  epistasis-prediction metrics, and reference-to-peak diagnostics.

## Inputs

- Record declared inputs in metadata.yaml.

## Configuration

- Track configuration in config.yaml.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- 1) branch comparisons between binding-only and paired abundance-plus-binding
  fits on the KRAS assays; 2) assay-specific `mavenn` diagnostics for abundance
  and DARPin K55 binding; 3) fitted thermodynamic-readout parameters and
  predictive validation metrics for each branch; and 4) a structured
  `summary.json`

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
