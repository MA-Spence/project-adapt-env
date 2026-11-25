# Unpublished candidate 0.6.0rc1: full pairwise interactions

On 19 November 2026 the scientist requested full 20-by-20 pair matrices and
pairwise-only latent interactions as the default, retaining the profile additive
prior and the sigmoid. This supersedes 0.5.0rc1 as the working default. It is an
implementation revision, not a package publication or a recalibration.

The model is

\[
f(x)=f_{\min}+(f_{\max}-f_{\min})\,
\sigma\!\left[\operatorname{logit}(p_0)+\frac{\beta}{\sigma_A}
\left(A(x)-A(x_0)+r\,[G(x)-G(x_0)]\right)\right],\quad
G(x)=\sum_{(i,j)\in E}J_{ij}(x_i,x_j).
\]

Each selected pair has a full Gaussian table, centered along both amino-acid
axes and scaled with the other tables so that its field has the same exact
uniform-background single-substitution RMS as the additive component.
There is no rank truncation or binary-contrast restriction. A centered q20 table
has 361 independent pure-pair degrees of freedom; storing all 400 entries avoids
reconstructing them during queries. Explicit third- and fourth-order latent
terms are absent; the sigmoid can still induce higher-order observable effects.

This is the pairwise interaction-expansion-plus-sigmoid architecture of
[Park, Metzger and Thornton (2024)](https://doi.org/10.1038/s41467-024-51895-5)
on the selected pair supports. Their fitted coefficients are not reproduced:
AdaptEnv still supplies profile-based additive effects and synthetic pair
coefficients. Default sparsity restricts the available site pairs. Setting
`target_degrees=(L-1,)` includes all pairs; at four sites the default already does.

## Tractability and controls

The default target degree is six, capped at L-1. With M selected pairs, persistent
numeric storage is O(Lq+Mq²), and each full fitness query costs O(L+M).
After one baseline score, each mutation update touches its incident pairs.
At bounded degree and fixed q=20, these are linear in length L. For n=20^L possible
genotypes, this is O(log n) storage and query cost. No complete-sequence lookup
table is built. At L=1,000 and degree six, M=3,000 and pair payload is 9,600,000
bytes. This excludes Python objects, profile/additive arrays, queries and scratch.
The existing support builder includes sorting and possible balancing retries;
we do not assert a strict linear worst-case construction bound.

Ultra-diverse input sequences have the same scoring cost as nearby sequences.
The profile construction itself may have alignment-processing costs separate
from the online oracle. Choosing all pairs changes storage and full-query cost
to quadratic in L, but never requires enumerating 20^L complete sequences.

Computability does not establish informative fitness variation at large L.
With uniform backgrounds and fixed mutation RMS, latent variance grows with L;
the retained generic steepness and r can therefore produce strong sigmoid
saturation. The smoke check establishes execution and storage bounds, not
calibration, learnability, or biological fidelity for distant orthologues.

## Evidence boundary and historical replay

The new default is **uncalibrated**. It inherits no TrpB goodness-of-fit claim.
All six calibrated TrpB presets explicitly retain `ranked_gnk`, their original
orders, weights, degrees and rank, and the `thermodynamic_gnk_v2` model identity.
Their numeric maps and fingerprints are regression-tested against the unchanged
0.5.0rc1 source at commit 745c00b3a369da33e7d7c7b9850ccbc212fac714.
The bundled RUN-178 profile is unchanged. The generic `synthetic_control` now
uses full pair tables with Gaussian additive effects; it is not a historical
calibrated control.

The earlier candidate's manifest, archive, analyses and evidence stay unchanged.
Its lock applies to its archived source, not this revised checkout. The current
manifest records this candidate's source and engineering checks. EXP-067 still
specifies the older ranked construction and is not silently retargeted here.

`tractability.py` reproduces a bounded local smoke check using profile-conditioned
models and 256 uniformly sampled q20 sequences at each of three lengths.
`tractability.json` records timings, payload sizes and process peak RSS; those
measurements are not a comparative benchmark or evidence for HYP-010.
`verification.json` records source and installed-wheel tests and provenance.
No new RUN, empirical analysis, release tag, deployment or publication is created.
