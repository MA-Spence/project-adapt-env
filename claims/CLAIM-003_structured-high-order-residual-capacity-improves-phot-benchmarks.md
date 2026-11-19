# CLAIM-003: Coupled Auxiliary Capacity Improves PHOT Ruggedness Benchmarks Without Empirical Module Specificity

- Status: completed
- Created: 2026-05-17
- Supersedes: HYP-008
- Aims: AIM-001
- Linked results: RES-035, RES-036, RES-037, RES-038, RES-039

## Claim

For PHOT, adding fitted coupled/auxiliary landscape capacity improves empirical
global and local ruggedness summaries over the current fitted model and smooth
core. The supported claim is about recoverable high-order/local ruggedness
capacity for benchmark construction. It is not a claim that the current
empirical alignment, contact, or module locations are the true amino-acid
interaction structures responsible for PHOT epistasis.

The DCA-, alignment-, covariation-, module-, or mutation-feature-routed terms
should be interpreted as an explicit, biologically legible scaffold that later
optimizer agents may be able to exploit. They should not be interpreted as
evidence that empirical routing is required to fit the landscape statistics.
Random, shuffled, no-module, and uncoupled controls are therefore
interpretation controls: they show that the fitted improvement is generic
capacity rather than validated empirical module specificity.

## Interpretation

The PHOT fitted current model under-represents the empirical graphlet package,
especially high-order and local ruggedness structure. Across the recent PHOT
experiments, giving the model ABC-fittable coupled/auxiliary capacity improves
held-out graphlet and localization summaries relative to fitted current and
smooth-core controls.

The random/shuffled/no-module results should not be read as simple negative
results for the benchmark program. They are negative for the stronger
mechanistic statement that empirical module locations are privileged
landscape-fit features. They are positive or neutral for the benchmark
statement that comparable ruggedness can be routed through several bases while
holding global/local difficulty summaries under experimental control.

This distinction matters for ML benchmarking. Alignment conditioning is kept
because it can make part of the ruggedness signal biologically legible to
pretrained protein language models or evolution-aware active-learning agents.
The benchmark claim does not require that biologically routed terms fit PHOT
better than random or shuffled terms. It requires that biologically routed and
non-biologically routed landscapes can be matched on landscape statistics, then
used to test whether agent priors change optimization behaviour.

## Evidence

- `RES-035` shows that the no-module coupled-operator improvement over fitted
  current and smooth core is stable across replicated calibration/test graphlet
  splits. No-module median held-out graphlet distance is `0.4839`, compared
  with `0.5232` for fitted current and `0.5358` for smooth core. The same
  result weakens empirical module specificity because uncoupled, random,
  shuffled, and empirical module families also improve against current.
- `RES-036` shows that SMC-ABC fitting improves the fixed PHOT held-out graphlet
  package for coupled/auxiliary families. The empirical-module family reaches
  held-out distance `0.2553`, compared with `0.4016` for fitted current and
  `0.4119` for smooth core. Random-module is fractionally better and
  uncoupled-module is clearly better, so this is capacity evidence rather than
  empirical routing evidence.
- `RES-037` shows that the PHOT empirical graphlet layer contains local
  ruggedness signal relative to fitted current/smooth core. Held-out empirical
  edge top-10 energy fraction is `0.3904` versus `0.3085` for fitted
  current/smooth core, and empirical site top-1 energy fraction is `0.3863`
  versus `0.2645`. The empirical calibration modules do not enrich that energy,
  which supports local-structure benchmarking but not module causality.
- `RES-038` shows that adding local Dirichlet-energy ruggedness terms to the
  ABC objective improves the held-out aggregate package over fitted current and
  smooth core. For the empirical-module family, full held-out distance improves
  from `0.5210` fitted current to `0.3653`; original graphlet metrics improve
  from `0.4016` to `0.2898`; and local ruggedness metrics improve from
  `0.6746` to `0.4623`. The no-module fitted family improves still further
  (`0.2771` full distance, `0.2384` original metrics, `0.3268` local metrics).
- `RES-039` adds an important boundary: simply nesting an additional
  soft-min-coupled latent landscape does not monotonically improve the
  calibration-selected empirical-module result. The best selected nested family
  is a random-module control (`0.4826`), while empirical-module (`0.5432`)
  loses to fitted current (`0.5210`), smooth core (`0.5268`), and the locked
  first-stage model (`0.5296`). This does not overturn the narrower capacity
  claim from `RES-035` to `RES-038`, but it rules out citing nested ensemble
  depth alone as support.

## Boundaries

This claim does not say that the fitted model identifies which residues are
truly epistatically linked in PHOT. The shuffled, random, no-module, and
uncoupled controls show that the current empirical module locations are not
validated as specific biological interaction terms.

This claim supersedes the stronger mechanism-specific version of `HYP-008`,
where empirical local modules themselves were the privileged source of PHOT
ruggedness. The supported statement is narrower and more useful for
benchmarking: PHOT contains recoverable high-order/local ruggedness signal, and
coupled/auxiliary residual capacity can improve held-out benchmark statistics
over the current fitted model without relying on specific local amino-acid
module identities.

`RES-039` adds that this capacity claim is not monotonic in architectural depth:
the nested third-landscape soft-min ensemble failed the calibration-selected
held-out empirical-module comparison. Use this claim for the scoped fitted
capacity signal, not as an argument that more latent landscapes necessarily
improve PHOT benchmarks.

## Implications

For benchmark construction, DCA-, alignment-, covariation-, or
mutation-feature-routed high-order terms remain scientifically defensible as an
interpretable scaffold. They should be reported as biologically legible routing
terms for agent-facing benchmarks, not as residue-level mechanistic
discoveries and not as necessary ingredients for landscape-statistic fit.

For downstream supervised models, this is the useful result: the benchmark can
include explicit high-order/local ruggedness while preserving ablation regimes
that separate landscape difficulty from biological side-information leverage.
Biologically routed, random-routed, shuffled-routed, no-module, uncoupled, and
diffuse-residual regimes should be compared only after matching the relevant
global and local landscape summaries.

## Data provenance

- `data/processed/phot-replicated-heldout-split-graphlet-diagnostic/RUN-089.dvc`
- `data/processed/phot-coupled-family-discriminating-abc/RUN-091.dvc`
- `data/processed/phot-localized-ruggedness-energy-diagnostic/RUN-093.dvc`
- `data/processed/phot-abc-local-ruggedness-objective/RUN-096.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`

## Required Before Manuscript Use

- Frame the claim as benchmark capacity and high-order/local-ruggedness
  representation, not residue-level interaction discovery.
- Include the negative-control boundary: empirical module locations are not
  privileged over shuffled/random/no-module/uncoupled controls under the
  current tests.
- Cite the held-out metrics from `RES-035`, `RES-036`, and `RES-038`, the
  empirical localization diagnosis from `RES-037`, and the nested-depth
  boundary from `RES-039`.

_Owner: Codex_
