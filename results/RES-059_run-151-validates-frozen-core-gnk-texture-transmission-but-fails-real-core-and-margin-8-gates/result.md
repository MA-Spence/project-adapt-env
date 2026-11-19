# RES-059: RUN-151 validates frozen-core GNK texture transmission but fails real-core and margin-8 gates

## Summary

RUN-151 completed the full EXP-061 assay and passed all 12 technical checks but
only `13/16` required checks overall, issuing `SIGMOID_TRANSFER_FAILED`. Pointwise
scalar transfer and frozen-core candidate texture pass; common-threshold real-core
edge coverage, actual-core texture, and margin-8 deployment fail. The sigmoid can
carry normalized GNK texture, but additive-core motion and the declared operating
points leave the signal confined to rare transition-window cases. This weakens
HYP-005 at the observable protein-scale boundary without testing the current
public generator or service end to end.

## Generated from

- Analysis: ANA-059
- Experiment: EXP-061
- Run: RUN-151

## Relevant hypotheses

- Supports: None
- Weakens: HYP-005
- Refutes: None

## Evidence

- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/metrics.json`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_acceptance_audit.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_collection_audit.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_source_hash_audit.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_pointwise_transfer_by_length.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_reference_anchored_edge_by_length.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_texture_transmission_by_length.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_deployment_margin_comparison_by_length.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_scale_law_summary.csv`
- `analyses/2026-07-23_ANA-059_run-151-operator-normalized-gnk-sigmoid-transmission-paired-residual-review/tables/run-151_protocol_output_audit.csv`
- `data/processed/operator-normalized-gnk-sigmoid-transmission-paired-residual-validation/RUN-151.dvc`
- `/tank/data/labproj/runtime/project-adapt-env/runs/RUN-151/outputs`

## Interpretation

RUN-151 executed all `180/180` reference/seed tasks over six fixed proteins at
lengths `48-448` and 30 independent GNK seeds. All 193 collected files match both
their remote SHA-256 digests and DVC object hashes. Runtime configuration and
implementation hashes match the present registered files.

The positive mechanistic result is real but narrow. Per-genotype pointwise
transfer at occupancy `0.5` has minimum median gain `0.99065`, practical fraction
`0.93465`, and log-length slope `-0.00205` with 95% interval
`[-0.01280, 0.00703]`. More importantly, the common-threshold frozen-additive-core
control transmits candidate-additivized texture at every length: all seed
projections are positive, every required permutation test passes, practical
fractions are `0.262-0.393`, and the texture-RMS slope is `-0.10100` with interval
`[-0.11819, -0.08858]`. The sigmoid itself can therefore transmit the normalized
GNK texture, including at all tested `L >= 100`.

The real additive core breaks that conclusion. Anchoring the reference core at
occupancy `0.5` leaves only about `1.4-1.5%` of one-step edges both-dynamic;
approximately `98.5%` are one-crossing. Eligible `n >= 32` seed counts are
`0/30, 4/30, 17/30, 30/30, 30/30, 30/30` with increasing length, failing the
registered 80% per-length coverage gate. The fixed count threshold is mechanically
harsher for short `19L` panels, so this is not evidence for a short-protein-specific
biological effect. Still, across all one-step edges, practical transmission is
only `1.48-1.67%` and correlation only `0.113-0.125`; the failure is not merely a
coverage artifact.

Actual-core `tau_texture` is practically collapsed at all lengths. Its RMS is
only `2.32e-12-1.21e-10`, practical fraction is zero throughout, permutation-pass
fractions are `0.289-0.467`, and its slope is `-0.33194` with interval
`[-1.31142, 0.57766]`. Near-unit correlation and positive formal projections do
not rescue a signal this small. Raw observable epistasis remains substantial, but
the additive-surrogate control also produces observable epistasis with numerically
zero latent epistasis. That establishes why raw `epsilon_R` cannot be interpreted
as transmitted GNK texture and why the additivized contrast is required.

The margin-8 base-template deployment point also fails. Direct one-step practical
fractions are `0-0.000995` and about `98.5-98.7%` of edges are exact float64 ties.
Stable analytic residuals have the same practical fractions, so numerical
subtraction is not the limiting explanation. Supported plaquettes retain positive
projection and `0.844-0.944` permutation-pass fractions, but only
`1.56-2.54%` exceed the `0.10` heuristic floor. Margin 10 is worse: one-step
practical fractions are zero, ties are `99.95-100%`, and texture practical
fraction is `0.003906` at every length. These are rare threshold-crossing effects,
not a broadly usable observable knob.

## Effect on hypothesis

RES-059 weakens HYP-005 specifically at the observable, protein-scale boundary.
The EXP-060 mutation-normalized latent scale law remains supported, and the
frozen-core result adds evidence that the sigmoid can transmit it under controlled
transition-window conditions. What fails is robust common-threshold transmission
when the real additive core moves, together with practical margin-8/10 deployment.

This is not a refutation of all additive-core plus GNK models or of sigmoid
transmission in principle. It is a failure of the tested observable construction
to make GNK strength a broadly directional and practically resolvable knob over
the declared protein panel.

## Limitations

- One protein represents each length, so protein identity and length remain
  confounded; seed uncertainty is conditional on this fixed panel.
- The fixed `n >= 32` edge rule is length-confounded through the complete `19L`
  panel size, although the all-edge result independently shows sparse transfer.
- The square table omits the full set of endpoint core occupancies promised by
  the prose, so double-endpoint saturation is inferred from the edge occupancy and
  actual-versus-frozen controls rather than directly tabulated.
- The runner aggregates seed texture projection as the median of per-order
  projections. A pooled sufficient-statistic calculation leaves pass states
  unchanged.
- The pointwise NRMSE gate is implemented only at occupancy `0.5`; shoulder NRMSE
  would exceed the same bound if the prose were read as applying it there.
- Observation noise is disabled and the `0.10` floor is a heuristic, not an
  empirical detection limit.
- The candidate is experiment-local and replaces the public fixed-support
  generator. Margin 8 is a base-template formula, not a universal service/API
  setting. The public model was not changed.
- No empirical realism, biological routing, accessibility, adaptive-walk,
  active-learning, or optimizer claim was tested.

## Downstream use

- Cite RES-058 for the latent mutation-normalized scale law.
- Cite RES-059 for the boundary that frozen-core sigmoid transmission works but
  real-core/common-threshold and margin-8 deployment do not.
- Do not cite RUN-151 as validation of the current public generator or service.
- Keep the Stage 4C/4D/4E hard stop in force.
