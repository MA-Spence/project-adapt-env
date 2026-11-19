# ANA-061: RUN-157 paired categorical-code and representation-rank GNK texture rescue review

## Purpose

Collect and independently audit completed RUN-157 against the preregistered
EXP-063 pairing, manipulation, core, edge, texture, midpoint-control,
dose-directionality, length-stability, paired-superiority, decision, and
provenance gates without post hoc retuning.

## Linked experiments/runs

- Experiment: EXP-063
- Run: RUN-157
- Direct hypothesis: HYP-005
- Motivating result: RES-060

## Collection and reproducibility audit

`RUN-157` completed on `lab-slurm` as Slurm job `2044` with scheduler exit code
`0`. It is the full registered configuration (`quick = false`): the 12
EXP-062 confirmation proteins, six length bins, five fresh additive-core seeds,
six fresh GNK seeds, and four paired representation arms. All `360/360`
protein/core/GNK contexts completed. `stderr.txt` is empty, and the aggregate
tables contain the declared `8,640` edge rows and `104,040` texture rows.

The collected bundle is:

- `data/processed/paired-categorical-code-and-representation-rank-gnk-texture-rescue/RUN-157`;
- `378` files and `456,516,694` bytes: 15 registered outputs, 360 strict-JSON
  checkpoints, and three logs;
- DVC pointer
  `data/processed/paired-categorical-code-and-representation-rank-gnk-texture-rescue/RUN-157.dvc`;
- DVC directory hash `7f2e06c1c337fa49b0e9c4006be443ff.dir`.

All workspace files match the DVC directory manifest and local cache objects by
MD5. All 378 files also match the intact shared-runtime copy by SHA-256 and
size. All 362 JSON files parse as strict JSON, with no `NaN` or infinity tokens.
Checkpoint identities exactly cover the registered `12 x 5 x 6` crossed
design, and their row totals reproduce the aggregate tables. The recorded
EXP-063 configuration and all ten implementation-component hashes match the
present files.

The configured redundant DVC backup is incomplete. The collection and local DVC
tracking succeeded, but the push to `/Volumes/DutchZFS` failed with permission
errors. Only `155/378` content objects are present there, `223` are absent, and
the required directory object is absent, so that remote cannot reconstruct the
bundle. This does not change the scientific decision because the local
content-addressed cache and the byte-identical shared-runtime copy are intact,
but it is a durability warning that remains in the record.

The repository and AdaptEnv checkout were dirty at submission. Exact hashes of
the relevant configuration and source components mitigate, but do not erase,
that broader provenance limitation.

## Registered and independently audited decision

The runner reports `3/3` required contract checks passing and issues:

`BINARY_RANK8_SUFFICIENT`

Independent reaggregation of the raw CSVs reproduces all 30 recorded pass
states, all paired and factorial estimates, all four length fits, and the
terminal decision. The shared additive core passes `60/60` cells, every one of
the 360 cross-arm contexts passes the support/coefficient/probe pairing and
manipulation checks, and every arm passes all `360/360` complete-h1 edge cells.

Only `binary_r8` satisfies the complete observable package plus paired
superiority:

| Arm | Actual texture | Midpoint control | Dose | Length | Paired superiority | Full package |
|---|---:|---:|---:|---:|---:|---:|
| `binary_r1` | fail | pass | pass | pass | baseline | fail |
| `dense_q20_r1` | fail | fail | pass | pass | fail | fail |
| `binary_r8` | pass | pass | pass | pass | pass | pass |
| `dense_q20_r8` | fail | pass | pass | pass | pass | fail |

All four arms pass their technical and complete-h1 edge gates.

## Reproduction of the incumbent failure

The fresh `binary_r1` arm reproduces RES-060 rather than escaping by a new seed
draw. Its actual-core supported texture has:

- median practical fraction `0.19792`;
- median texture RMS `0.42520`;
- median latent/observable correlation `0.74033`;
- median projection `0.92874`;
- positive projection-bootstrap lower bound `0.89004`;
- only `45/360` conjunctive cell passes.

The per-length pass fractions are `0.0167`, `0.0500`, `0`, `0`, `0`, and
`0.6833` for `L093`, `L118`, `L220`, `L238`, `L448`, and the short bin,
respectively, far below the required `0.80` in every bin. Its transition-window
midpoint control, complete-h1 edge transfer, dose directionality, and length
gate all pass. The failure is therefore localized again to practical
plaquette-texture prevalence under the actual additive core.

## Binary rank-8 rescue

`binary_r8` passes the actual-core texture gate in every `360/360` context and
in every length stratum. Its median practical fraction is `0.46094`, compared
with `0.19792` for `binary_r1`. At the protein level, the paired median gain is
`+0.26237`, with crossed-bootstrap lower bound `+0.24805`, comfortably beyond
the preregistered `+0.05` and positive-lower-bound requirements.

This rescue is not an increase in total observable texture energy. Median
texture RMS is nearly unchanged (`0.42520` for rank 1 versus `0.42439` for rank
8), as are correlation and projection. Rank 8 instead redistributes the fixed
second-moment dose:

- median latent near-zero fraction falls from `0.67057` to `0.05729`;
- median exact-zero fraction falls from `0.04427` to `0`;
- the median of per-cell latent absolute medians rises from numerical zero
  (`9.33e-15`) to `0.30679`;
- median latent excess kurtosis falls from `6.5137` to `1.2132`;
- median latent RMS remains nearly unchanged (`0.53981` versus `0.53628`).

Thus eight independent binary product components spread similar RMS energy over
many more supported plaquettes, allowing more contrasts to cross the fixed
`0.10` output threshold. The supported causal statement is about the complete
rank-8 binary representation used here; the experiment does not isolate
algebraic rank from the associated eight-component averaging operation.

## Dense-code arms

Removing binary allele ties is not sufficient and is detrimental in this
construction. `dense_q20_r1` eliminates sampled nonself categorical ties, yet
its median practical fraction falls to `0.10677`, `0/360` actual-core cells
pass, and only `23/360` midpoint-control cells pass. Its typical latent effects
are extremely small despite finite RMS: the median of per-cell absolute medians
is `0.00694`, while median excess kurtosis is `67.09`. Dense rank 1 therefore
replaces exact binary ties with many tiny effects and rare large tails rather
than distributing energy practically.

`dense_q20_r8` is a partial improvement, not a registered rescue. Its median
practical fraction is `0.26563`; `255/360` cells pass the full texture
conjunction, and its paired protein-level gain over `binary_r1` is `+0.06771`
with lower bound `+0.05859`. However, per-length pass fractions decline from
`0.90/0.80` at `L093/L118` to `0.617`, `0.583`, and `0.40` at
`L220/L238/L448`; the minimum is far below `0.80`. Its full package therefore
fails even though its paired gain is positive.

The registered two-by-two factorial summaries make the pattern explicit:

- rank-8 main effect on practical fraction: `+0.21126`,
  95% protein-bootstrap interval `[+0.20736,+0.21842]`;
- dense-code main effect: `-0.14323`,
  interval `[-0.14551,-0.14160]`;
- dense-by-rank interaction: `-0.09701`,
  interval `[-0.11719,-0.07943]`.

These are effects of the complete tested representations. Exact uniform-q20
nonself mutation-operator energy and analytic field variance are matched, but
realized reference-shell dose and higher moments are not forced to be identical
within each protein. The complete-h1 RMS ratio is especially variable for
`dense_q20_r1` (`0.591-2.934`) and less so for `dense_q20_r8`
(`0.763-1.494`). The binary rank comparison is cleaner
(`0.882-1.100` versus `0.902-1.094`) and retains nearly identical observable
RMS.

## Edge, dose, and length behavior

All four arms pass every complete-h1 edge cell. At `a=1`, median practical edge
fractions are `0.883`, `0.736`, `0.885`, and `0.859` for `binary_r1`,
`dense_q20_r1`, `binary_r8`, and `dense_q20_r8`, respectively.

All arms also pass the preregistered `a=0.25,0.5,1,2` directionality package.
For `binary_r8`, median texture RMS increases
`0.110 -> 0.219 -> 0.424 -> 0.773`, while median practical fraction increases
`0.208 -> 0.349 -> 0.461 -> 0.505`. Every `binary_r8` cell and protein passes.

All registered length bounds pass, but this is bounded stability rather than
invariance. The `binary_r8` protein-level log-length slope is `-0.12347`, with
interval `[-0.13356,-0.11330]`; the interval excludes zero while remaining
inside the preregistered limits. The other slopes are `-0.09684`,
`-0.03659`, and `-0.07700` for `binary_r1`, `dense_q20_r1`, and
`dense_q20_r8`.

## Interpretation for HYP-005

RUN-157 gives narrow positive evidence for HYP-005. On the fixed noiseless
challenge panel, an experiment-local rank-8 binary representation turns GNK
amplitude into a dose-directional, practically prevalent, and
registered-length-bounded observable texture knob. It also shows that the
specific RES-060 failure is repairable and that low-rank energy concentration,
not binary ties alone, is the operative bottleneck in this implementation.

This does not invalidate RES-060: the incumbent rank-one construction again
fails. Nor does it establish that rank 8 is minimal, unique, biologically
necessary, or sufficient beyond these 12 previously examined proteins.
Rank was tested only at endpoints 1 and 8, and the population-level unit
remains the protein, not the 360 crossed seed cells.

The candidate remains experiment-local. No public AdaptEnv method was changed,
no empirical fitness outcome or observation noise was used, and no structural
contact mechanism, DCA/Potts inference, accessibility, adaptive walk,
active-learning, or optimizer behavior was tested. The public-model validation
state and the Stage 4C/4D/4E hard stop therefore remain unchanged.

## Limitations

- The 12 proteins are the already examined EXP-062 confirmation panel. Fresh
  core and GNK seeds strengthen paired causal diagnosis and procedural
  replication, but do not provide independent protein generalization.
- Protein-level uncertainty is based on only 12 fixed proteins, two per length
  bin.
- Rank 1 versus rank 8 does not establish the minimum sufficient rank or
  separate rank from component count and averaging.
- Expected mutation-operator energy is controlled exactly; realized local dose
  and higher moments are diagnostic and vary most strongly in the dense arms.
- The `0.10` contrast floor, `0.25` prevalence floor, and `0.80` stratum rule
  are preregistered noiseless engineering criteria, not empirical detection
  limits.
- The repository and external AdaptEnv checkout were dirty at submission,
  although every relevant source hash matches.
- The configured DVC backup is partial and non-reconstructable; the local
  content-addressed cache and shared-runtime copy are intact and verified.

## Produced artifacts

- `tables/metrics.json`
- `tables/run-157_acceptance_audit.csv`
- `tables/run-157_arm_summary.csv`
- `tables/run-157_arm_by_length.csv`
- `tables/run-157_paired_effects.csv`
- `tables/run-157_factorial_effects.csv`
- `tables/run-157_dose_curve.csv`
- `tables/run-157_scale_law_summary.csv`
- `tables/run-157_representation_audit.csv`
- `tables/run-157_protocol_output_audit.csv`
- `tables/run-157_source_hash_audit.csv`
- `tables/run-157_collection_audit.csv`
