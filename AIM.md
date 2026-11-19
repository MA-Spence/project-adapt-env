# Project aims


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
