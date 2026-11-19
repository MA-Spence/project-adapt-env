# RES-044 / EXP-045 follow-up: did the A–D fixes resolve the epistasis calibration?

Addendum to `docs/RES-043_epistasis_calibration_review.md`. Reviews RUN-117, the
four-stage cumulative ablation (A, A+B, A+B+C, A+B+C+D), ANA-044, and the staged
fast-fit code in `scripts/phot_trait_space_bounded_gate_fast_fit.py`.

## Verdict: partly confirmed the diagnosis, did not resolve the problem, and exposed the real pathology

The experiment is a clean, well-designed ablation and it earns its negative result.
It confirmed the half of my diagnosis that concerned the **observable scale**, and
it falsified the implicit hope that fixing the readout would also fix accessibility.
That falsification is the valuable part.

What the stages show against empirical (test graphlets; empirical magnitude 0.2926,
SD 0.2704, ratio 1.046, accessibility 0.1975):

| stage | obs. magnitude | raw magnitude | mag/SD ratio | accessibility | guardrails | distance |
|---|---|---|---|---|---|---|
| A (monotone link) | 0.2289 | **178.6** | 1.10 | **0.0164** | 8/8 | 0.4235 |
| A+B (+global sigmoid) | 0 | 0 | 0 | 0 | 0/8 | 2.921 |
| A+B+C (+ratio knob) | 0 | 0 | 0 | 0 | 0/8 | 2.921 |
| A+B+C+D (+conditioned obj.) | 0.2426 | **1.16e-6** | 1.12 | **0.0190** | 8/8 | 0.4122 |

Three things are now established:

1. **The readout-class argument was right.** Swapping affine→monotone-piecewise
   (A) lifted observable magnitude from 0.1708 (RUN-115) to 0.2289, and the ratio
   knob + scale-invariant objective (C+D) pinned magnitude/SD to 1.12 against
   empirical 1.046. The ratio is no longer the failure. Good.
2. **None of A–D touched accessibility.** Every non-collapsed stage sits at
   0.016–0.019 against empirical 0.1975 — an order of magnitude too low, exactly
   as in RES-043. This is now the dominant, unrepaired failure.
3. **The failure mode flipped rather than closed.** RUN-115 had raw latent
   epistasis far too large (~104); RUN-117's best stage has raw latent variation
   essentially **zero** (1.16e-6) while the readout manufactures the observable
   epistasis. Both extremes fail. The model can be pushed off one cliff onto the
   other but not onto the ridge between them.

## Why accessibility could never have been fixed by A–D (the key structural point)

Fixes A and B are both **monotone, single-input transforms** of one additive latent
(A is the monotone-piecewise observation map; B, as implemented in
`candidate_updates`, is a single `sigmoid` `observed_fitness_terms` term on
`stability_margin`). A monotone transform of a scalar is **rank-preserving**, so it
cannot create or destroy sign epistasis, local maxima, or blocked paths. I verified
the bookkeeping consequence in the runner: `record_candidate_metrics` computes all
metrics on `apply_observation_readout(raw_values, readout)`, and that readout is
monotone, so **accessibility of the observable equals accessibility of the raw
latent**. Accessibility is therefore a property of the latent surface alone, and no
amount of readout engineering (A, B, C, D) can move it. My RES-043 memo praised
monotone maps for preserving accessibility while generating magnitude epistasis —
true, but the flip side is that they cannot *create* accessibility either. The
experiment cleanly separated these two halves and landed on the "sourcing is the
defect" branch I flagged in the §6 test.

So the only thing in the current model that can produce sign epistasis / maxima /
reduced accessibility is the **specific pairwise J**, and it is sampled essentially
at random (low-rank random tables + covariation bias). Random specific epistasis at
the amplitude needed to register on sign/maxima shatters paths far more aggressively
than the empirical landscape: empirical PHOT holds sign fraction 0.61 **and** ~5
maxima **and still** 20% accessible paths. Random J cannot occupy that corner — at
matched sign fraction it leaves ~2% accessible. Empirical ruggedness is **structured
and navigable**; the model's is **random and shattering**.

## The two new pathologies, named

**P1 — Accessibility is a latent-geometry deficit, not a readout deficit.** The
model lacks a generator of *correlated, low-dimensional* ruggedness that is rugged
yet navigable. A 1-D additive core + monotone readout is fully accessible (~100%);
adding random J overshoots into shattered (~2%). There is no setting in between
because the model has no structured intermediate. (Gavrilets "holey" landscapes;
the navigability-vs-ruggedness literature — random NK/House-of-Cards gives low
accessibility at given ruggedness, whereas correlated/biophysically-structured
landscapes give much higher accessibility at the same ruggedness.)

**P2 — Readout flexibility made the latent non-identifiable.** Because the
observation readout is re-fit per candidate to map model→`latent_phi`, the objective
rewards *any* latent whose rank order can be monotonically mapped onto empirical.
A near-flat latent with infinitesimal structure (1.16e-6) can be blown up by a steep
piecewise segment to reproduce the observable summaries on calibration — a
readout-induced artifact that does not generalise (held-out distance still 0.41,
accessibility still 0.019). The monotone-piecewise link is too flexible and the
objective places no constraint on the *raw* latent geometry, so stage D's pressure
to match observable scale was satisfied by collapsing the latent rather than by
shaping it.

**P0 (already half-noted) — the global sigmoid (B) collapse.** A+B and A+B+C went to
a flat landscape (distance 2.921, everything 0). The selected B parameters had
`global_readout_midpoint = 0`, `slope = 2.575` while stability margins sit near
+9.6, i.e. the entire genotype distribution lands in the saturated tail of the
sigmoid → constant output. The pre-D objective did not penalise this because a flat
landscape trivially satisfies the near-zero-mean signed-epistasis terms. D repaired
the *conditioning* of that collapse but, lacking a latent constraint (P2), repaired
it by sliding to the opposite degeneracy.

## What this means for the path forward

Keep A, C, D — they correctly solved the observable scale and ratio, and D fixed the
worst objective conditioning. The remaining work is entirely in the **latent
generator** and in **constraining it during the fit**:

1. **Source ruggedness from a low-dimensional, multi-input, non-monotone latent, not
   from random J or a 1-D readout.** The natural candidate is already in the
   codebase: `coupled_landscape.soft_min_values` (weakest-link of two smooth additive
   traits, e.g. stability ⊗ function). Two additive traits combined by soft-min
   create sign epistasis and reduced accessibility *near the trait-crossover
   boundary* — structured, navigable ruggedness — rather than the path-shattering
   noise of random J. This is the generator that can occupy the "0.61 sign / 5 maxima
   / 0.20 accessible" corner. It is also what RES-032/033/034 were probing; the
   missing step is to combine it with the now-working A/C/D observable calibration
   instead of testing it readout-naked.
2. **Constrain the raw latent in the objective so it cannot collapse or explode.**
   Add explicit targets on *raw* quantities: raw single-mutant SD, raw accessibility,
   raw local-maxima count — and fit the observation readout once on calibration and
   **freeze** it (or restrict it to a 2-parameter monotone form, e.g. a Hill/sigmoid,
   not a free piecewise) so the latent must carry the structure. This directly closes
   P2.
3. **Make accessibility a first-class calibration target, measured on the raw
   latent.** It is currently slipping through because it is readout-invariant and the
   objective is dominated by observable summaries the readout can fake.

## Honest correction to my prior memo

My §6 predicted that the monotone link "should" let magnitude rise *and*
accessibility recover together. Only the first happened. I under-weighted that a
monotone readout is structurally incapable of repairing accessibility — accessibility
was always a latent-geometry problem. The experiment's main contribution is to prove
that cleanly: the readout fixes are necessary for scale but provably insufficient for
topology, and the next iteration must change the latent generator (multi-trait
soft-min), not the readout.

## Caveat unchanged

PHOT is still the only empirical system in the loop; treat any latent-generator fix
as validated only when it holds the magnitude/ratio *and* accessibility *and*
local-energy package on a second assay (the GB1/TrpB target).
