# CLAIM-002: Empirical Biophysical Fits Are Supported Only As A Single Latent Thermodynamic Trait

- Status: completed
- Created: 2026-05-14
- Supersedes: HYP-007
- Aims: AIM-001
- Linked results: RES-007, RES-008, RES-009, RES-010, RES-011, RES-012, RES-013, RES-014, RES-015, RES-016, RES-017, RES-018, RES-019

## Claim

The current evidence supports a reduced biophysical model in which empirical
fits are driven by one latent thermodynamic trait, interpreted as stability or
stability margin. The evidence does not support `HYP-007` as originally stated:
additional fitted latent molecular traits have not improved empirical recovery
as a complete landscape model.

## Interpretation

The project has repeatedly tested the idea that poor empirical recovery is
mainly caused by missing extra latent phenotypes such as activity, binding,
abundance, or photophysical readout structure. Those additions sometimes
improve one predictive slice, especially double-mutant ranking in diagnostic
models, but they do not improve the full target package of single-mutant
ranking, double-mutant ranking, epistasis recovery, distributional fit, and
non-pathological reference geometry.

The more coherent interpretation is that the empirically recoverable signal in
the current assays is dominated by one thermodynamic coordinate. For these
fits, the biologically realistic abstraction is not a generic functional latent
space plus extra fitted traits, but a single stability-like latent trait with
structured pairwise effects expressed on that thermodynamic coordinate.

## Evidence

- `RES-007` shows that stability-targeted calibration materially improves SPTN1
  recovery, supporting the thermodynamic single-trait direction.
- `RES-008`, `RES-009`, `RES-010`, and `RES-014` show that explicit PHOT
  function, latent-activity, and two-trait formulations do not rescue empirical
  recovery. `RES-014` makes this result reproducible rather than a one-off
  optimisation failure.
- `RES-011` and `RES-013` are the strongest comparative evidence: PHOT
  stability-readout controls outperform explicit two-latent alternatives, and
  the advantage remains under the cleaner validation-objective comparison.
- `RES-012` shows that a paired KRAS abundance/binding thermodynamic readout
  does not rescue predictive recovery.
- `RES-015` and `RES-016` show useful but incomplete diagnostic gains from more
  flexible two-latent or residual-trait formulations. These runs improve some
  ranking metrics but fail to recover epistasis as a package, so they constrain
  rather than support the original multi-latent claim.
- `RES-017` shows that an oracle global-epistasis baseline does not beat the
  PHOT stability control, weakening the argument that the missing ingredient is
  merely a more flexible latent or observation layer.
- `RES-018` and `RES-019` show that paired state/function diagnostics on OXDA
  and KCNJ2 do not improve activity or function recovery over function-only
  branches, within the limitation that those assays are single-mutant-only and
  do not test double-mutant epistasis.

## Implications

`HYP-007` is superseded by a narrower claim: biophysical modelling is useful
for empirical fits, but the supported form is a single thermodynamic latent
trait. Additional fitted latent traits should not be treated as part of the
public empirical core unless future evidence shows that they improve the full
landscape-recovery package, not just isolated prediction metrics.

For benchmarking, this is acceptable if the benchmark is presented honestly as
a thermodynamic synthetic environment. It is not an all-purpose model of every
molecular determinant of protein function. It is a controlled setting where
fitness is generated from a stability-like latent coordinate, with local
statistics calibrated to an empirical system.

## Boundaries

This claim does not say that real proteins have only one molecular phenotype.
It says that, in the current empirical tests and available readouts, fitting
extra latent traits has not produced a more accurate or more biologically
defensible model than the single thermodynamic trait. Future orthogonal data
could justify a richer mechanistic model, but that is not supported by the
current record.
