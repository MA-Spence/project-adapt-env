# Unpublished 0.7.0rc1 implementation candidate

The scientist requested Gamma concentration and Optuna indirect inference using
the published Park–Metzger–Thornton genetic-architecture estimator. This candidate
implements that request; it is not a published release or a calibrated TrpB preset.

The online model retains the profile additive field, full pairwise tables and
sigmoid. Finite `interaction_concentration=alpha` draws per-entry variances from
Gamma(alpha, rate=alpha), multiplies the existing independent standard normal
coefficients by their square roots, then applies the existing centring and exact
mutation-RMS normalization. `None` preserves the exact Gaussian limit, including
old map fingerprints. No unsupported finite alpha is made the default before
calibration. `global_intercept` optionally places the sigmoid relative to the
uniform-background mean; the old reference-occupancy convention remains exact
when this option is absent. The prior candidates and RUN-175/178/180 are frozen.

Offline inference follows the pinned authors' sigmoid SSE+L1 objective and
cross-validation penalty selection. Sparse one-hot computation, a Python
split-variable solver and algebraic gauge conversion are documented implementation
choices, not an assertion of identical R solver output. The authors' L1 also
penalizes the lower/width readout terms, which is retained here. Objective and
gradient values are checked against an independently executed upstream R fixture.

Optuna fits r, alpha, steepness, global position, lower and span. It matches
inferred-versus-inferred architecture with fixed calibration-bootstrap weights,
identical genotype grouping, masks, observation rules and generator seed lists.
It saves every trial, failure diagnostics and a descriptive near-best region.
Fresh-seed confirmation is separate; optional confirmation genotypes must be
disjoint. Sign epistasis, oriented diminishing returns and induced-panel walks
remain outside the loss. No global graph or genotype identity enters the loss.

The Gamma draw law, summary vector and fixed diagonal distance are declared
generative/inference assumptions; they are not conclusions of Thornton's paper.
The method is indirect inference, not the paper's original outer objective.
The [protocol](../../calibration/README.md) and package's `docs/RFA_IMPLEMENTATION.md`
state these limits and the exact numerical settings. Optuna's trials are not
posterior samples. No new scientific EXP/RUN/ANA/RES/CLAIM is registered here.

Online storage and full-query cost remain O(L) at fixed q and pair degree.
Offline all-pairs inference is quadratic in sites and repeated over CV, seeds,
bootstrap panels and trials; a timing pilot is required before selecting a
biological-run budget. Implementation tests do not establish empirical adequacy,
long-sequence biological realism or identifiability of fitted knobs.

The separate validation environment does not modify the pinned RUN-180 Python
environment. Historical ranked fingerprints depend on Python floating-point sum
semantics: the original fixture uses Python 3.13. Same-interpreter comparison to
frozen RUN-180 source preserves all 28 checked old maps in Python 3.11, and the
historical fixtures are tested unchanged in Python 3.13. The archived provenance
must accompany exact replay. No legacy numerical formula or fixture was changed
to force cross-interpreter fingerprint equality.
