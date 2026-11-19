# Model release candidate 0.5.0rc1 — profile-conditioned additive effects

The scientist selected alignment/profile-conditioned additive preferences as the
default on 19 November 2026, following ANA-070. This is an **unpublished release
candidate**, not a published package or deployment. No tag or release was created.

The reason for this choice is explicit biological marginal information in the
additive component. ANA-070 finds mixed descriptor trade-offs against the recorded
synthetic control; it establishes neither consistent superiority nor statistical
equivalence. ML/pLM benefits remain untested. CLAIM-004 keeps its synthetic-prior
scope, and HYP-003's specific biological-conditioning predictions remain open.

## Locked behavior

- Package version: **0.5.0rc1**. Equation/model-family identifier:
  **thermodynamic_gnk_v2**. Numerical normalization, interactions and sigmoid
  are unchanged.
- `ThermodynamicConfig.additive_prior` defaults to `profile`. Ordinary Python,
  service and HTTP initialization requires alignment/profile data. There is no
  automatic protein selection, remote alignment retrieval or Gaussian fallback.
- `additive_prior="synthetic"` and the `synthetic_control` preset explicitly
  select the original Gaussian control. Explicit matrices and effect samples
  remain advanced supplied-data overrides, accurately labelled by `core_prior`;
  they are not automatically treated as alignment-conditioned biology.
- `trpb_run178_r05_profile_shared_stop_intrinsic` bundles the original four
  profile rows and the matching RUN-178 readout at the inherited r=0.5 setting.
  `build_landscape_from_preset` loads them together. The config-only builder
  refuses to drop the bundled profile.
- Calibrated presets permit generator-seed changes, but reject other parameter,
  profile or reference substitutions. Users intentionally changing those inputs
  must create an explicitly custom regime without the calibrated-preset label.
- Runtime observations remain noiseless. Shared-stop describes the research
  calibration sensitivity; the optional API Gaussian noise is a different model.
  Reference sequences remain seed-specific computational anchors.

The source and packaged-data identities are in `manifest.json`. `lock.py --check`
verifies them; subsequent scientific or API changes require a consciously revised
candidate lock. Source remains in the independent `external/Adapt-Env` package.

## Evidence and verification

The comparison is [ANA-070](../../../analyses/2026-11-19_ANA-070_run-178-trpb-profile-prior-versus-synthetic-control/analysis.md),
using frozen RUN-178 and RUN-175. This decision does not change either run or its
analysis. Historical experiment scripts retain their frozen source assumptions;
reproduction uses their archived package/environment rather than the new default.

The maintained package checks cover required profile inputs, prior conflicts,
length inference, direct/service/HTTP agreement, preset substitution guards,
archived synthetic scores and the packaged TrpB prior/readout. Installed-wheel
checks also verify that the bundled profile is available outside the checkout.
Exact results, environment and wheel hash are recorded in `verification.json`.

`replay.py` independently compares both prior-matched presets with frozen source
and outputs on all 24 recorded seeds per prior, using 1,020 sampled latent
sequences and 128 recorded cube sequences per seed. `replay.json` records exact
current-versus-archived numerical agreement, and differences no larger than
2.23e-15 against saved output values. This is an engineering replay of intrinsic
maps, not new calibration, an observation-model replay or whole-landscape proof.
Fingerprints change because configuration now records the prior-selection policy.

The DVC archive and local-cache verification are recorded in `storage.json`.
Remote backup remains pending while Dutch is unavailable. The archive is a local
candidate snapshot and does not constitute publication.
