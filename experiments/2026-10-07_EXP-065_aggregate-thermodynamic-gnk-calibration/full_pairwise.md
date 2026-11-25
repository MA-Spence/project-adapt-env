# EXP-065: matched TrpB calibration of the full-pairwise model

On 20 November 2026 the scientist authorized a new calibration RUN under the
existing EXP-065, followed by local MacBook execution. This protocol repeats
RUN-178 with the refined 0.6.0rc1 interaction construction. RUN-178 is the direct
matched comparator; RUN-175 remains the recorded fully synthetic control.
Neither control is rerun. No new experiment, ANA, RES or CLAIM is created here.

## One model intervention

Use the locked `thermodynamic_pairwise_v3` implementation: the same profile
additive matrix, full centered 20-by-20 pair tables, pairwise latent terms only,
and the existing sigmoid. The four sites contain all six pairs. Exact mutation-RMS
normalization and the original seed-selected reference policy are retained.
The interaction generator changes from rank-eight mixed-order products to full
pair tables. This comparison tests that combined revision; it cannot separately
attribute a change to rank, binary contrasts, or removal of higher-order terms.

The model uses the Park–Metzger–Thornton pairwise functional architecture, with
synthetic interactions rather than fitted empirical coefficients. Calibration
still fits only the global readout. This run does not infer biological interaction
identities or establish biological validity outside four-site TrpB.

## Matched procedure

`full_pairwise.yaml` preserves RUN-178's complete non-model scientific design:

- The same observed 159,129/160,000 genotype mask, empirical fitness, calibration
  and confirmation cubes, geometry, profile matrix and probability smoothing.
- The same eight calibration seeds, 16 reused and eight fresh confirmation seeds,
  cube assignments, reference sequences and observation RNG seeds.
- Interaction ratios 0, 0.5, 1, 2, 4 and 8; observations `none`,
  `independent_difference` and `shared_stop_sensitivity`; unchanged donor pools,
  donor exclusions, observation uniforms and fixed noise fitting procedure.
- Refit the six r=0/r=4 endpoint maps with RUN-172's original optimizer procedure
  and seeds; then use those anchors in RUN-175's 12 intermediate searches.
  Preserve the quantile/local-SD/floor objective, bounds, population size eight,
  maximum 60 iterations, optimizer seed rules, tolerance and polishing policy.
  Actual evaluation counts can differ through the unchanged stopping rule.
- The same complete-library graph mask, uniform-uphill walk kernel, top-one-percent
  definitions, resolutions 0/0.001/0.01 of the calibration range, WHT/DR summaries,
  intrinsic and observed maps. These diagnostics are not fit targets.

Expected full outputs are 18 selected settings, 24 seed audits, 2,163 graph rows,
2,163 cube-summary rows and 144 latent graph rows, with complete latent fields,
cube values, fitting traces and provenance. There are no scientific adequacy
thresholds or favourable-seed selection. Protocol identity means the same
procedure and inputs, not identical nonzero-r fitted results.

## Preflight and execution

Use the same local environment as RUN-178: Python 3.11.15, NumPy 2.4.4, pandas
2.3.0, SciPy 1.17.1, PyArrow 20.0.0 and PyYAML 6.0.3. Record the native MacBook
platform, Apple Clang/SDK, numerical libraries and installed packages. Use four
workers with one BLAS/OMP thread each; local resource metadata is advisory.
Do not use the historical project dependency overlay. A fresh compiled graph
kernel must pass its original validation. No shortened calibration is substituted.

The local helper verifies the new source lock, the original two input bundles,
and the complete RUN-178 output manifest. It preserves a prelaunch snapshot of
source, complete inputs and environment, then registers a nonexecuting local
preview through labproj. A separate authorized `--execute` invocation registers
and executes the live RUN. The local helper collects outputs, captured logs and
RUN metadata because the local target lacks remote collection-path settings.
Temporary submission settings are restored. Outputs are never overwritten.

Original numerical calibration, observation and diagnostic sources are retained.
`landscape_graph.py` has a pre-existing trailing-whitespace cleanup relative to
the recorded source; the new guard verifies both trailing-whitespace-normalized
bytes and Python AST identity and records the differing hashes. It does not
treat arbitrary source edits as equivalent. Package source is checked against
the frozen 0.6.0rc1 candidate; none is changed by this calibration preparation.

## Completion and independent checks

Check unchanged profile additive values and references across calibration and
all 24 confirmation fields, unchanged empirical arrays/masks and cube assignments,
and expected new interaction fields. The r=0 branch is a negative control: its
parameters/objectives and graph/cube outcomes should replay RUN-178 within
declared floating-point tolerance. Both models use the same local numerical
environment, but bit identity is assessed rather than assumed.

An independent verification script checks manifest hashes, coverage, matched
configuration, inputs, seed cohorts, additive/reference invariance, empirical
outputs and the r=0 control. Technical completion requires return code zero,
full output counts, intact manifests and passing engineering checks. Record
failures without changing optimization or selecting seeds to improve agreement.
Archive the collected RUN with DVC and verify the local cache; remote backup is
separate. Technical completion is not a conclusion that the refined model fits
TrpB adequately and does not change the locked model's default readout.
