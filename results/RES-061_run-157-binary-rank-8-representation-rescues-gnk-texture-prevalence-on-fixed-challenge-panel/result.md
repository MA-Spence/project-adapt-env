# RES-061: RUN-157 binary rank-8 representation rescues GNK texture prevalence on the fixed challenge panel

## Summary

RUN-157 reproduced the incumbent binary-rank-1 prevalence failure and issued
the preregistered `BINARY_RANK8_SUFFICIENT` decision. At nearly unchanged
median texture RMS, binary rank 8 raised median practical texture prevalence
from `0.198` to `0.461` and passed all `360/360` paired cells, with a
protein-level median gain of `0.262` and bootstrap lower bound of `0.248`.

Dense q20 coding alone worsened prevalence, and dense q20 rank 8 remained below
the per-length package gate. This narrowly supports an experiment-local,
dose-directional GNK knob on the reused challenge panel. It is not independent
protein generalization or public-model validation, and the Stage 4C/4D/4E hard
stop remains.

## Generated from

- Analysis: ANA-061
- Experiment: EXP-063
- Run: RUN-157

## Relevant hypotheses

- Supports: HYP-005
- Weakens: None
- Refutes: None

## Evidence

- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/metrics.json`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_acceptance_audit.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_arm_summary.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_arm_by_length.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_paired_effects.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_factorial_effects.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_dose_curve.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_scale_law_summary.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_representation_audit.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_protocol_output_audit.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_source_hash_audit.csv`
- `analyses/2026-07-24_ANA-061_run-157-paired-categorical-code-and-representation-rank-gnk-texture-rescue-review/tables/run-157_collection_audit.csv`
- `data/processed/paired-categorical-code-and-representation-rank-gnk-texture-rescue/RUN-157.dvc`

## Interpretation

The incumbent failure is procedurally reproducible. With fresh core and GNK
seeds, `binary_r1` retains strong alignment, positive permutation support, and
complete-h1 edge transfer, but only `45/360` actual-core cells pass the practical
texture gate. Its median practical fraction is `0.19792` and median texture RMS
is `0.42520`. Its midpoint control, dose response, and registered length bound
pass, localizing the failure to practical supported-plaquette prevalence.

The binary rank-8 repair is decisive on this panel. `binary_r8` passes every
actual-core cell in every length bin. Median practical fraction rises to
`0.46094`, and the paired protein-level gain over `binary_r1` is `+0.26237`
with crossed-bootstrap lower bound `+0.24805`. Complete-h1 edge, midpoint
control, dose-directionality, projection, permutation, and length gates all
pass.

The rescue is redistribution rather than a larger dose. Median observable
texture RMS is almost identical for binary rank 1 and rank 8
(`0.42520` versus `0.42439`). Rank 8 instead lowers the median latent near-zero
fraction from `0.67057` to `0.05729`, raises the median typical absolute
contrast from numerical zero to `0.30679`, and lowers excess kurtosis from
`6.51` to `1.21`. Similar second-moment energy is spread over many more
plaquettes, so more effects cross the fixed output threshold.

Removing binary ties is not the remedy. `dense_q20_r1` has no sampled nonself
ties but performs worse: median practical fraction is `0.10677`, `0/360`
actual-core cells pass, and the midpoint control fails. Its median typical
absolute contrast is only `0.00694` and excess kurtosis is about `67`, showing
that exact zeros were replaced by many tiny effects plus rare large tails.

`dense_q20_r8` is a partial improvement, not a full rescue. Its median practical
fraction is `0.26563`, `255/360` cells pass the texture conjunction, and its
paired gain over `binary_r1` is positive (`+0.06771`, lower bound `+0.05859`).
But per-length cell-pass fractions fall to `0.617`, `0.583`, and `0.40` at
`L220`, `L238`, and `L448`, below the required `0.80`.

The factorial analysis agrees: the rank-8 main effect is `+0.21126`
(`[+0.20736,+0.21842]`), the dense-code main effect is `-0.14323`
(`[-0.14551,-0.14160]`), and the interaction is `-0.09701`
(`[-0.11719,-0.07943]`). The supported diagnosis is therefore low-rank energy
concentration, not binary allele ties, within this GNK construction.

All arms pass dose directionality and the registered length bounds. For
`binary_r8`, median texture RMS rises
`0.110 -> 0.219 -> 0.424 -> 0.773` over `a=0.25,0.5,1,2`, while practical
fraction rises `0.208 -> 0.349 -> 0.461 -> 0.505`. Its length slope is
`-0.12347` with interval `[-0.13356,-0.11330]`: within the registered bound,
but not invariant with length.

## Effect on hypothesis

RES-061 supports HYP-005 narrowly. An experiment-local binary rank-8
representation makes the GNK amplitude dose-directional and practically
resolvable after the frozen additive core and sigmoid on this fixed challenge
panel. It also demonstrates that the RES-060 rank-one prevalence failure is
repairable.

This does not erase RES-060: the incumbent rank-one arm fails again. It does not
show that rank 8 is minimal, unique, biologically necessary, or sufficient on
new proteins. The public method is unchanged and remains unvalidated at the
full Stage 4 boundary.

## Limitations

- The 12 proteins were reused from EXP-062. Fresh seeds provide paired
  mechanistic replication, not independent protein generalization.
- Rank was tested only at 1 and 8; rank is confounded with the number and
  averaging of independent product components.
- Expected uniform-q20 nonself mutation-operator energy is matched, but realized
  reference-shell dose and higher moments vary, especially for dense rank 1.
- The practical thresholds are noiseless engineering criteria, not empirical
  assay detection limits.
- No empirical outcome, observation noise, DCA/Potts inference, structural
  mechanism, accessibility, adaptive walk, active learning, or optimizer was
  tested.
- The repository and AdaptEnv checkout were dirty at submission, although the
  relevant recorded source hashes all match.
- Local DVC objects and the shared-runtime copy are complete and verified, but
  the configured redundant DVC backup is partial and cannot reconstruct the
  bundle.

## Downstream use

- Cite RES-061 for fixed-panel causal evidence that binary rank 8 repairs the
  practical texture-prevalence failure at nearly unchanged RMS dose.
- Cite RES-061 for evidence against binary-tie removal as a sufficient repair.
- Do not cite RUN-157 as independent protein generalization, biological
  validation, DCA/Potts evidence, or validation of the public AdaptEnv method.
- Keep the Stage 4C/4D/4E hard stop in force.
