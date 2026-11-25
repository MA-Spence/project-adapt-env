# RES-064: Fixed observation noise adds epistasis and limits TrpB architectural identifiability

## Summary

**A fixed genotype-level observation-error map adds epistasis outside AdaptEnv's
latent interaction control, r. It does not account for all epistasis in the
synthetic landscape. The existing four-site TrpB calibration data do not uniquely
distinguish shared genotype-level error from biological epistasis.**

This is an identifiability limit of the current observations and assumptions,
not a claim that experimental noise and biological effects can never be separated.
The synthetic components are known and separable exactly. The biological
interpretation of the empirical components is not identified by this calibration.

This result synthesizes EXP-065 RUN-186/ANA-072, RUN-188/ANA-073 and
RUN-190/ANA-074. It introduces no new experiment or fit. The accepted use of
observation noise is an explicit assay-modeling assumption; acceptance does not
turn stop covariance into an identified estimate of biological noise.

## Evidence

The strongest direct test is RUN-188: the same profile-conditioned full-pairwise
model and fixed shared-stop observation law were fitted with r fixed to zero
or allowed to vary. Both arms used 11,712 Optuna trials and the same aggregate
fitness loss. The selected full arm retained r=0.5. Primary confirmation summaries
below use eight originally fresh generator seeds, eight cubes per seed, and 32
empirical cubes. These seeds are not independent biological replicates.

| Quantity | r=0 | Full, r=0.5 | Empirical TrpB |
|---|---:|---:|---:|
| Achieved calibration loss | 0.04207450 | 0.03099517 | — |
| Intrinsic sign-epistatic face fraction, seed median | 0.0% | 9.9% | Not observable directly |
| Observed mean-fitness sign fraction, seed median | 62.5% | 64.1% | 65.1% |
| Repeatable directed-sign fraction, seed median | 43.0% | 39.8% | 38.7% |

Sign classifications use the primary resolution, 0.001 times the frozen empirical
calibration scale (epsilon=0.00019396876295). Full-model loss is 26.3% lower than
the r=0 fit. Neither search established a global optimum; both retained queued
historical candidates. Advantages in unfitted diagnostics were mixed.

At this resolution, a median **90.4% of observed sign-epistatic faces in the full
model were absent in its intrinsic map**. This changes to 72.3% at zero resolution
and 98.3% at the coarser tested resolution. It is a synthetic face-classification
comparison, not an estimate of the fraction of real TrpB epistasis caused by noise.

Noise dominance depends on the metric. ANA-073's exact pooled power accounting
across the 64 fresh-seed cubes gives:

| Spectral order | Intrinsic after-sigmoid fitness | Shared error | Mean independent error | Signed cross terms |
|---|---:|---:|---:|---:|
| 2 | 91.1% | 6.6% | 0.9% | 1.4% |
| 3 | 61.3% | 30.6% | 3.8% | 4.3% |
| 4 | 16.3% | 72.3% | 15.7% | −4.2% |

These are ratios of pooled powers to total observed power, with signed cross
terms retained. They are not nonnegative causal variance fractions. In particular,
the intrinsic landscape still supplies most order-2 power, despite noise's large
influence on near-neutral sign classifications and order-4 power.

### Reproducible analysis and figures

- [ANA-072: fitted landscape, architecture and static difficulty](../../analyses/2026-11-25_ANA-072_run-186-trpb-landscape-behaviour-architecture-and-difficulty-proxies/analysis.md).
- [ANA-073: exact observation/interaction ablation](../../analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/analysis.md) and its [executed figure notebook](../../analyses/2026-11-25_ANA-073_run-188-structured-interactions-versus-observation-randomness/notebooks/analysis.ipynb).
- [ANA-074: Gamma versus Gaussian](../../analyses/2026-11-25_ANA-074_run-190-gamma-concentration-versus-gaussian-trpb-calibration/analysis.md) and its [executed figure notebook](../../analyses/2026-11-25_ANA-074_run-190-gamma-concentration-versus-gaussian-trpb-calibration/notebooks/analysis.ipynb).

Canonical numeric evidence is in ANA-073's `tables/power_budget_fresh_pooled.csv`,
`tables/seed_layer_stratum_summary.csv`, `tables/paired_model_contrasts.csv` and
ANA-074's `tables/calibration.csv` and `tables/key_numbers.json`. The linked
analyses contain their archive hashes, validation reports and reproduction commands.

## Interpretation

There are three distinct sources of observed nonadditivity:

1. **Latent interactions:** the centered full 20×20 pair tables, scaled by r
   relative to additive single-mutation RMS. Their strength and concentration
   are explicit architectural controls.
2. **Nonlinear readout:** the sigmoid transforms the additive-plus-pairwise
   latent score. It can generate spectral epistasis even when r=0. An r=0
   monotone additive model nevertheless has no intrinsic sign reversals.
3. **Fixed observation errors:** genotype-specific residual and shared terms
   are added after the sigmoid. These introduce further interaction contrasts,
   including sign reversals, without changing the latent pair tables or r.

The observation law has no fitted noise-amplitude parameter in these runs.
Its magnitude/distribution is fixed from calibration-only donors and stop controls.
However, the readout changes which fitness-conditioned donor bins are visited.
Consequently, observed ruggedness is not determined by r alone, even with a
fixed observation law. Noise is not an additional fitted architectural interaction
coefficient, but it is an additional source of apparent epistasis in the oracle.

The replicate model is

\[
Y_1(x)=F(x)+C(x)+e_1(x),\qquad Y_2(x)=F(x)+C(x)+e_2(x).
\]

Here F is intrinsic after-sigmoid fitness, C is shared error, and e1/e2 are
independent donor draws. All are frozen for a genotype in the recorded maps.
Repeating an oracle query returns the same error. Averaging replicates reduces
independent error but retains C. Across-replicate repeatability therefore does
not establish that a feature is biological.

For any genotype-dependent function h, replacing F by F+h and C by C−h leaves
both observed replicates unchanged. Replicate differences remove both F and C;
they inform independent error, not their separate contributions. Stop controls
constrain C only under assumptions about stop biology, normalization and transfer
of their error law to non-stop variants. Current data do not establish those
assumptions sufficiently to identify this decomposition uniquely. Negative
stop-adjusted powers near the floor emphasize the limitation; they are not
negative biological variances.

Thus **the model's noise contribution can be computed exactly, but cannot simply
be transferred as an empirical TrpB noise fraction**. The ablation shows that
structured interactions contribute and improve the achieved aggregate fit; it
does not prove that its fitted r is the unique biological interaction strength.

## Gamma concentration follow-up

RUN-190 jointly fitted alpha, r and the four readout coordinates with the same
loss and fixed observation law. The selected alpha=100, r=0.5 candidate had loss
0.03103213 versus Gaussian 0.03099517 (0.12% higher). Fresh-seed latent pair-table
correlation with Gaussian was 0.99875; top-1% coefficient energy was 8.53% versus
8.40%. Confirmation changes were small and mixed. This achieved fit gives no
clear advantage for stronger effect concentration and does not resolve the
observation-layer identifiability limit. It is not proof of equivalence or global
optimality. Gaussian remains the default.

## Effect on hypothesis

HYP-004/HYP-005's aggregate TrpB resemblance remains conditional on readout and
observation assumptions. This result narrows interpretation of apparent genetic
architecture; it does not refute the latent model family or establish a biological
mechanism. RES-063's broader research boundary remains: these results do not
establish computational superiority, independently controlled learning/search
partitions (HYP-011), or biological transfer (HYP-012). No learner or optimizer
benchmark was conducted in the analyses synthesized here.

## Limitations and downstream use

Scope is the assayed four-site TrpB system, the recorded panel and fixed assay
model. Activity is highly sparse, and sign metrics are resolution dependent.
Previously inspected confirmation panels are not new independent validation in
later runs. The selected models retain discrepancies in activity, powers and
seed variation; “all epistasis reproduced” and unrestricted biological graduation
are unsupported.

The model is usable as a conditional research environment if intrinsic fitness,
observed fitness, latent r and observation semantics are kept distinct. A benchmark
using the frozen observed map studies learning/search on that composite map;
it must not attribute all of its difficulty to latent r.

## Public method implementation

AdaptEnv exposes full-pairwise Gaussian/Gamma fields, the direct joint Optuna
calibration, fixed-r ablation, and explicit replicate-error/bank operations as
public Python APIs. [Package calibration documentation](../../external/Adapt-Env/docs/CALIBRATION.md)
describes their inputs, evidence and distinction from the separate experimental
architecture-inference fit. Historical RUN snapshots remain unchanged.

[Public API parity validation](public_api_parity.json) is reproduced by
[validate_public_calibration_res064.py](../../scripts/validate_public_calibration_res064.py).
It checks the frozen observation draws/bank, latent fields and all loss components
against the actual RUN records. This is implementation verification, not new
scientific evidence. The source update does not publish a package release or
promote a new calibrated preset.
