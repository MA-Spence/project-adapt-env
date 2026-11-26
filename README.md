# AdaptEnv scientific project

This repository records the development and empirical evaluation of AdaptEnv,
a procedural synthetic protein-fitness environment for online ML testing.
The supported implementation lives in the independently installable
[Adapt-Env package](external/Adapt-Env/README.md), now the unpublished **0.7.0rc1** implementation candidate.

## Current model and usage

The default model uses alignment/profile-derived additive amino-acid preferences
with full centered 20-by-20 pair tables, followed by a sigmoid. The latent
interactions are pairwise only, on a sparse site graph with six partners per site
by default. A supplied profile
is required; fully synthetic additive effects are an explicit control. It evaluates arbitrary sequences
of a configured length without enumerating their complete sequence space.

- [How the model works](external/Adapt-Env/docs/THERMODYNAMIC_V2.md)
- [Python and HTTP APIs](external/Adapt-Env/docs/API_AND_SERVICE.md)
- [Configuration and presets](external/Adapt-Env/docs/CONFIGURATION_AND_PRESETS.md)
- [Full-pairwise migration](external/Adapt-Env/docs/MIGRATION_0_6.md)
- [TrpB calibration: restored RUN-180 targets with Optuna](docs/calibration/README.md)
- [Release decision and source lock](docs/model_release/0.7.0rc1/decision.md)

Install from the project root with:

```bash
python -m pip install -e 'external/Adapt-Env[api]'
```

```python
from adaptenv import ThermodynamicLandscape
from adaptenv.alignment import AlignmentProfile

# Illustrative alignment; supply the mapped alignment for your protein family.
profile = AlignmentProfile.from_sequences(["ACDE", "ACDF", "YCDE"], pseudocount=0.4)
model = ThermodynamicLandscape(alignment_profile=profile)
fitness = model.evaluate(model.random_sequence())
```

The project-side `project_adapt_env` namespace delegates its model API to the
same package when `src` is on the Python path. Existing statistical helpers remain available through their
explicit submodules. Historical experiment scripts can require an earlier
package revision; reproduce them with their frozen RUN source and environment.

## Current empirical evidence

[CLAIM-005](claims/CLAIM-005_scalable-generative-genetic-architecture.md) states the
current generative-model contribution, published foundation, architectural
controls, observation boundary and computational scaling.

[RES-064](results/RES-064_fixed-observation-noise-epistasis-and-trpb-identifiability/result.md)
synthesizes ANA-072/073/074: fixed observation errors add epistasis outside
latent r control but do not account for all intrinsic interaction power. Current
TrpB replicates and stop controls cannot uniquely distinguish shared errors from
biological epistasis. Gamma concentration gave no clear advantage over Gaussian.
The direct joint Optuna and fixed-observation methods are now
[public package APIs](external/Adapt-Env/docs/CALIBRATION.md).

**The project calibration default restores RUN-180's fitness-distribution
objective, using Optuna in place of differential evolution.** It preserves the
profile, Gaussian full-pairwise model, all three
observation modes, native-unit bounds, cubes, seeds and unfitted diagnostics.
The current protocol jointly optimizes continuous r in [0,8] and four readout
parameters in three Optuna studies. Gamma concentration is disabled. A new
TrpB fit completed as RUN-186; RUN-188 tested r=0 and RUN-190 tested Gamma
concentration. No public preset is recalibrated here.
See the [protocol and objective](docs/calibration/README.md).

The alternative Gamma/architecture calibration completed as RUN-184, but
[ANA-071](analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/analysis.md)
found a much too narrow fitness distribution and substantial local-landscape
mismatches. Its achieved fit is not ready as a TrpB surrogate. The architecture
inference API remains available explicitly; this rollback does not remove model
capabilities. The exact Gaussian generator remains the generic default until a
finite concentration is supplied. The established claims below concern
the historical ranked backend, retained explicitly by every named TrpB preset.
Its numerical behaviour and fingerprints are preserved. The new architecture
follows the additive-plus-pairwise expansion through a sigmoid of
[Park, Metzger and Thornton (2024)](https://doi.org/10.1038/s41467-024-51895-5);
its synthetic coefficients and sparse site graph are modelling assumptions.

[CLAIM-004](claims/CLAIM-004_thermodynamic-gnk-trpb-aggregate-and-adaptive-search-resemblance.md)
records useful approximate joint resemblance to TrpB global statistics and
adaptive-search behaviour under a specified observation model. Its evidence is
[ANA-069](analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/analysis.md)
and the linked RUN-175 data. The claim preserves observation dependence, seed
variation, pooled-WHT discrepancies and limits on cross-assay and length transfer.

[ANA-070](analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/analysis.md)
compares profile-conditioned RUN-178 with the recorded synthetic control. It
finds mixed descriptor trade-offs, without consistent superiority or established
statistical equivalence. The scientist selected the profile-conditioned model
as the release default for its explicit biological marginal preferences; ML/pLM
benefits remain untested. CLAIM-004 retains its original synthetic-prior scope.

[RES-063](results/RES-063_profile-model-readiness-and-learning-search-research-boundary/result.md)
synthesizes that evidence and distinguishes current model controls from untested
benchmark advantages. New prospective HYP-010, HYP-011 and HYP-012 address
computational efficiency, learning versus search difficulty, and the role of
biological information. Their definitions, controls and falsifiers are in the
[hypothesis record](docs/hypotheses.md#hyp-010-procedural-landscapes-provide-competitive-computation-and-efficiency).

The public generator is intrinsically noiseless. The TrpB research observation
model is separate; the API's optional Gaussian draws do not reproduce it.
The new `trpb_run178_r05_profile_shared_stop_intrinsic` preset packages the
matched profile and intrinsic readout. The RUN-175 preset remains an explicit
synthetic control. Neither applies an assay-specific observation layer.

## Scientific record and reproduction

- [Aims](AIM.md), [hypotheses](docs/hypotheses.md), and [current state](PROJECT_STATE.md)
- `experiments/`: EXP protocols and RUN metadata
- `analyses/`: executed notebooks, figures and numerical audits
- `results/` and `claims/`: interpreted evidence and manuscript-facing claims
- `registry/`: machine-readable record links
- `data/processed/`: durable datasets tracked with DVC
- `external/`: independently versioned source packages
- `src/` and `scripts/`: project analysis and experiment helpers

Scientific records preserve earlier evidence even when runtime implementations
are removed. RUN bundles retain exact source snapshots and input hashes. Data
restoration and queue operations use [labproj](.agents/skills/labproj/SKILL.md).
Keep `.env` and machine credentials untracked. Agent guidance is in [AGENTS.md](AGENTS.md).

## Current checks versus historical tests

The maintained package suite is `python -m pytest external/Adapt-Env/tests`.
The accumulated project overlay has a documented stale dependency lock and base
stack mismatch; see [environment validity](environment/README.md). This is separate
from the tested standalone model installation above.
The release replay in `docs/model_release/0.5.0rc1/replay.py` compares both
profile and synthetic presets against frozen RUN-178/RUN-175 source and outputs.
Its numerical checks do not rerun calibration or replace the scientific analyses.

Historical experiment scripts and project experiment tests preserve their
original synthetic-default assumptions and frozen source identities. Run them
with their recorded package/environment; they are not a compatibility target
for the new profile-required default. The current public project API delegates
to the release candidate, including `build_landscape_from_preset`.
