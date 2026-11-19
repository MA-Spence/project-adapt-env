# EXP-048: PHOT sparse Walsh Potts spectrum-calibrated generator probe

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Can a spectrum-calibrated sparse hashed Potts/Walsh field, routed through the existing additive stability core, pass the PHOT forward generative capability gate and then recover held-out RUN-098 graphlet equivalence under the recent numerical fast-fit protocol?

## Pre-experiment prediction

If RES-046 failed because soft-min is structurally concave/unimodal, the sparse Walsh/Potts field should clear the forward multi-peak/high-order capability probe before sequence-conditioned fitting; if it clears H1 but fails held-out equivalence, the remaining limitation is likely the sequence-to-support map or graphlet/support alignment rather than the operator itself.

## Rationale

The scratch design note RES-046_scalable_multipeak_generator_design.md argues that the next PHOT repair should test Generator A: additive core plus sparse hashed high-order Walsh/Potts terms with empirical spectrum calibration, structured support clusters, and a forward generative capability probe rather than an inverse graphlet fit. This experiment implements A only and keeps the recent A/F/J calibration infrastructure, restricted monotone readout, frozen RUN-077 PHOT source, and RUN-098 comparator.

## Experimental design

- Implement Generator A from `scratch/RES-046_scalable_multipeak_generator_design.md` only: additive stability core plus sparse hashed high-order Potts/Walsh terms.
- Use structured support clusters aligned to the PHOT graphlet/module pool so local graphlets can actually observe the routed epistasis.
- Run staged protocol:
  - `G`: measure empirical PHOT Walsh variance-by-order on calibration graphlets, bootstrap a high-order noise floor, estimate support-hit probability by Monte Carlo over the aligned cluster pool, and freeze the global order weights.
  - `G+H1`: run the forward generative `sparse_walsh_field` capability probe over a small grid around the frozen spectrum and `epi_weight`; empirical summaries must fall inside the generated central interval jointly.
  - `G+H1+FIT`: if H1 passes, promote the sparse-Walsh form to sequence-conditioned deterministic/random fast fitting with the anchored scale, scale-invariant raw objective, restricted monotone readout, and standard RUN-098 held-out equivalence check.
- Do not run the FGM/H2 arm in this experiment. The scratch note treats FGM as a harder complementary setup; this record tests A first.
- Report raw and observable metrics together.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
- `scratch/RES-046_scalable_multipeak_generator_design.md`
- `scripts/phot_trait_space_bounded_gate_fast_fit.py`
- `external/Adapt-Env/adaptenv/landscape.py`

## Configuration

- `config.yaml` defines the frozen PHOT source, RUN-098 comparator, G/H1/FIT stages, sparse-Walsh spectrum calibration, forward capability gate, fast-fit budget, and quick smoke-test override.

## Execution

- Run through `labproj submit EXP-048 --target lab-slurm --execute`.
- `run.sh` pulls RUN-077/RUN-098 DVC targets if missing, then runs `scripts/phot_trait_space_bounded_gate_fast_fit.py`.

## Expected outputs

- sparse_walsh_spectrum_calibration.csv
- sparse_walsh_forward_probe.csv
- bounded_gate_candidates.csv
- selected_bounded_gate_candidates.csv
- bounded_gate_distances.csv
- bounded_gate_equivalence.csv
- graphlet_metrics.csv
- graphlet_summary.csv
- walsh_terms.csv
- gft_coefficients.csv
- summary.json

## Analysis plan

- Create an ANA record after execution.
- First interpret H1. If H1 fails, do not overread missing sequence-conditioned equivalence rows: the operator failed the forward generative gate or the empirical high-order target/noise floor is suspect.
- If H1 passes but FIT fails held-out equivalence, attribute the limit to sequence-to-support mapping, support/graphlet alignment, or readout/fast-fit constraints rather than to sparse-Walsh operator capacity alone.
- If FIT passes held-out equivalence, treat it as a PHOT-local candidate generator only; confirm on a second assay before benchmark use.

## Completion criteria

- [x] Config committed
- [x] Run script committed
- [ ] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- This experiment follows the scratch design note as a prospective protocol, not as evidence.
