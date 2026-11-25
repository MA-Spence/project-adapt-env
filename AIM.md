# Project aims

2026-11-25 execution authorization: register and submit a matched TrpB-only
interaction ablation under EXP-065. Independently refit r=0 and full r in [0,8]
with the frozen RUN-186 shared-stop observation law and restored direct fitness
loss. Evaluate repeatable interaction powers and sign reversals without fitting
those diagnostics. No noise-strength/Gamma fitting or package promotion.

2026-11-24 execution authorization: jointly optimize continuous r in [0,8] and
the four readout parameters using the restored RUN-180 direct fitness loss,
then register and submit a TrpB-only RUN under EXP-065 to Slurm. Keep Gaussian
interactions; Gamma shape fitting is explicitly deferred. Preserve data, seeds,
observation assumptions and post-fit diagnostics; this authorization does not
assert model adequacy or change historical claims.

2026-11-21 execution amendment: stop local RUN-182 at the scientist's request
and submit the unchanged TrpB architecture calibration to the remote Slurm
queue under EXP-065. Preserve the interrupted outputs and restart in a new RUN
from the original frozen source, inputs, seeds and fitting configuration.
Only the execution platform changes; this does not report calibration results.

2026-11-21 execution authorization: register a new TrpB-only RUN under existing
EXP-065 and execute the 0.7.0rc1 Gamma-normal / Optuna architecture calibration
locally on the MacBook. The new inferred-architecture objective and bounded,
disjoint genotype panels are explicitly different from RUN-180's matched-repeat
protocol. Historical controls remain intact; this authorization does not imply
successful calibration, adequacy, publication or recalibrated public presets.

2026-11-20 implementation authorization: add Gamma-normal pair-coefficient
concentration and Optuna indirect inference against the reference-free genetic
architecture inferred by Park, Metzger and Thornton (2024). The same published
auxiliary objective and CV rule are used on empirical and synthetic observations;
the outer minimum-distance objective, Gamma sampling assumption and solver
differences are explicit. This advances the existing synthetic-environment aim,
without fitting empirical interaction identities or asserting biological adequacy.
It does not authorize or report a new empirical calibration RUN. The new method
and tests are documented in docs/calibration and unpublished candidate 0.7.0rc1.


## AIM-001: Create-Model-Landscape

- Status: active
- Owner: matthew-spence

Based on FGM theory, build a procedural and synthetic fitness landscape model that resembles real, biological landscapes from the perspective of population genetics, while also accounting for molecular phenomena such as stability neutral networks and specific epistasis.

2026-10-07 scope clarification from the scientist: the immediate implementation
uses a stability-like latent coordinate and nonlinear readout with synthetic
interactions. The target is conditional recovery of aggregate empirical
statistics, including WHT variance and diminishing returns, with scalable
procedural evaluation. Recovering actual empirical interaction identities is
outside this aim's current implementation scope. One universal parameter setting
across assay classes is not required. EXP-065 tests this narrower target directly.

2026-10-08: EXP-066 extends this aggregate target across existing protein assay
classes with suitable measured combinatorial coverage, for production figure
data and online-environment justification. It preserves mismatches, dimension
limits and observation assumptions; no scientific acceptance gates are imposed.

2026-10-27 contribution clarification from the scientist: the intended novelty
is an efficient implementation for on-demand exploration of large, sparsely
sampled sequence spaces, including intermediate and phylogenetic scales of
adaptation. Novelty in the underlying model class is not an objective beyond
the restrictions needed for scalability. The current implementation remains a
procedural, sparse, low-rank member of the established interaction-expansion plus
sigmoid family. Sparse sampling of sequence space is distinct from sparsity of
the interaction representation. This intended application does not establish
biological validity at large evolutionary distances or a computational advantage
over existing bounded-neighbourhood NK/GNK implementations. Existing empirical
claims and their limits remain unchanged; no new experiment or hypothesis is
created by this clarification.

2026-11-19 authorized intervention: restore an optional alignment-profile
additive prior and repeat the TrpB branch of EXP-065 / RUN-175 with the prior as
the sole scientific change. Reuse the recorded synthetic control and retain
synthetic interactions, original calibration settings and unfitted diagnostics.
This introduces sequence-specific marginal evolutionary preferences, without
claiming to recover empirical interactions or testing PLM performance. Existing
TrpB results do not automatically validate the changed prior.

2026-11-19 release decision: following the mixed profile-versus-synthetic
comparison in ANA-070, the scientist selected alignment/profile-conditioned
additive effects as the default model for unpublished release candidate
0.5.0rc1. This preference is grounded in explicit biological marginal information;
it is not a claim of consistent empirical superiority or statistical equivalence.
Synthetic interactions and the numerical readout remain the same model family.
The synthetic additive prior remains an explicit control. No release, deployment,
new experiment or new empirical claim is authorized by this decision.

2026-11-19 benchmark scope extension authorized by the scientist: retain the
profile-conditioned release candidate as a working research environment and
test (HYP-010) practical computability and comparative efficiency, (HYP-011)
separation of function learnability from fitness searchability, and (HYP-012)
the contribution of biological information to joint and discordant difficulty.
RES-063 synthesizes the existing evidence and the boundary of these new questions.
These are prospective hypotheses; comparative advantage, four-quadrant coverage,
and biological explanation are not conclusions of the existing TrpB experiments.
Synthetic interventions can identify mechanisms within the constructed model;
claims about real biological phenomena require independently measured fitness
data and explicit tests of the resulting predictions. This extension does not
reinstate retired implementations or overturn their historical failed tests.
The scientist authorized Git updates to the source and project repositories;
package release, tagging and deployment remain outside this action.

2026-11-19 backend revision authorized by the scientist: use unrestricted centered
q20 pair matrices and pairwise-only latent interactions as the default, retaining
the profile additive prior and sigmoid. This anchors the functional architecture
to Park, Metzger and Thornton (2024), without claiming fitted empirical
coefficients. Bounded site degree retains tractability for hundreds of variable
sites; sparsity is an explicit modelling restriction. Candidate 0.6.0rc1 is
uncalibrated. Historical TrpB presets and their ranked-backend evidence remain
separate and numerically preserved. Calibration and comparative benchmarking
are not executed by this implementation change.

2026-11-20 calibration authorization: prepare and locally execute a new TrpB
RUN under EXP-065 for the refined full-pairwise model. Repeat RUN-178's protocol,
profile, data, calibration and confirmation seeds, observation assumptions and
unfitted diagnostics, changing only the interaction construction. Reuse the
recorded RUN-178 and RUN-175 controls. This authorizes execution and preservation
in the scientific record, without pre-judging model adequacy or publishing a
release.
