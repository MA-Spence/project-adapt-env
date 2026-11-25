# TrpB calibration protocols

## Public methods and completed evidence (26 November 2026)

The direct loss, joint Optuna optimizer, optional Gamma fit, and fixed replicate
observation operations are now installed AdaptEnv Python APIs. See
[package documentation](../../external/Adapt-Env/docs/CALIBRATION.md).
RUN-186/188/190 are completed; ANA-072/073/074 and
[RES-064](../../results/RES-064_fixed-observation-noise-epistasis-and-trpb-identifiability/result.md)
record the achieved fits and observation-layer identifiability boundary.
The Gamma follow-up did not establish an advantage; Gaussian remains the default.
Frozen study runners remain unchanged. The descriptions below retain the
historical protocol distinction; preparation statements refer to their original
implementation, not the current completion status.

## Current default: joint continuous-r Optuna fitting

The scientist subsequently requested joint optimization of **r continuously on
[0,8] and all four readout coordinates**, with the same restored RUN-180 loss.
The default is now `joint-readout`; see the
[joint protocol](../../experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/joint_readout.md).
There are three studies, one per existing observation assumption, with 11,712
trials each and six historical warm starts each. Historical r values are candidate
starts, not a search grid. Gamma concentration remains disabled; alpha is not
optimized. Model, data, masks, seeds, observation draws and diagnostics retain
RUN-180's meaning. The five-parameter search produces three selected models.

## Superseded fixed-r implementation: RUN-180 targets with Optuna

The initial 2026-11-24 rollback restored RUN-180's fitness-distribution calibration while
replacing differential evolution with Optuna's TPE. It changes the
**project's EXP-065 TrpB workflow**, not the signature or behavior of the separate
`adaptenv.calibration` architecture-inference API. This is a prepared protocol;
no new calibration RUN or improved fit is reported by the implementation.

That superseded implementation remains explicitly available as `ADAPTENV_PROTOCOL=optuna-readout`, using
[`optuna_readout.yaml`](../../experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/optuna_readout.yaml),
[`optuna_pairwise_calibration_exp065.py`](../../scripts/optuna_pairwise_calibration_exp065.py)
and [`optuna_readout_exp065.py`](../../scripts/optuna_readout_exp065.py).
Historical RUN sources and configurations remain unchanged. The earlier
`full-pairwise` dispatch retains the original differential-evolution procedure.

The restored design keeps:

- Profile-derived additive preferences, Gaussian full centered 20 × 20 pair
  matrices, pairwise-only latent terms, and the reference-anchored sigmoid used
  in RUN-180. Gamma concentration and the optional uniform-background intercept
  remain available model features; neither is fitted in this protocol.
- The six fixed interaction ratios `r = [0, 0.5, 1, 2, 4, 8]` and all three
  observation modes: `none`, `independent_difference` and
  `shared_stop_sensitivity`. Each of the 18 settings gets its own readout fit;
  `r` is not an additional continuous Optuna variable.
- The original profile, calibration and confirmation cubes, eight calibration
  generator seeds, 16 reused plus eight fresh confirmation seeds, measured
  genotype mask, observation donor pools and fixed noise draws.
- The original calibration-only affine units, four readout parameter bounds,
  and unfitted WHT, epistasis, diminishing-returns and whole-library graph
  diagnostics. Confirmation data never select readout parameters or trials.

Optuna fits the same four parameters: standardized output offset, positive output
span, reference sigmoid position and positive sigmoid steepness. Their bounds
are respectively `[-2, 2]`, `[0.001, 10]`, `[-14, 8]` and `[0.05, 8]`; span and
steepness are sampled on logarithmic scales. The original conversion to native
assay units is retained, so the old fitted readouts remain inside the search
space. This reference position differs from RUN-184's uniform-background
position and must not be copied between the two parameterizations directly.

For the observed, standardized simulated cube values \(y\) and calibration
cube values \(t\), the loss is exactly

\[
\mathcal L =
\frac{1}{101\operatorname{Var}(t)}
\sum_{p\in\{0,.01,\ldots,1\}}[Q_y(p)-Q_t(p)]^2
+\left(\frac{\overline{\mathrm{SD}}_{\rm cube}(y)}
{\overline{\mathrm{SD}}_{\rm cube}(t)}-1\right)^2
+[F_y(c)-F_t(c)]^2.
\]

Quantiles pool the calibration cube values and include both endpoints. The SD
term uses the mean within-cube population SD. The final term compares fractions
at or below the same standardized near-floor threshold \(c\). Observation modes
are applied before computing these terms, including averaging the same paired
synthetic observations where applicable. All three terms have unit weight;
there is no auxiliary coefficient inference or bootstrap weighting.

The budget is **1,952 Optuna trials per setting**, including queued historical
warm starts; preliminary replay checks are additional objective evaluations.
RUN-180's same-setting readout and endpoint readouts are re-evaluated on the frozen calibration inputs before new proposals can replace
them. They are candidates, not copied results or held-out selection. Retaining
the best evaluated candidate prevents optimization from discarding a better
warm start on this objective; it does not guarantee improved confirmation
behavior or global optimality. Trial histories and all objective components are
saved.

This restores the procedure associated with earlier useful TrpB regimes, not an
unconditional biological adequacy claim. In particular, shared-stop remains a
sensitivity assumption, not an identified assay-error mechanism. WHT and search
outcomes remain diagnostic rather than fitting targets. Historical controls are
reused, not rerun.

## Experimental alternative: architecture-based calibration

Implemented in unpublished candidate **0.7.0rc1**, this remains explicitly
available through `scripts/calibrate_architecture.py` and
`adaptenv.calibration`. It is no longer the selected project TrpB calibration
workflow. RUN-184 completed, but [ANA-071](../../analyses/2026-11-24_ANA-071_run-184-trpb-architecture-calibration-and-landscape-adequacy/analysis.md)
found severe fitness-range compression and local-landscape mismatches; the
achieved fit does not support graduation. Its auxiliary architecture and
bootstrap stability limitations remain recorded, without rejecting the model
family. The method below is preserved for reproducibility and explicit research
use, not silently substituted for the restored protocol.

This workflow fits a **synthetic ensemble to empirical genetic architecture**.
It uses the same reference-free additive-plus-pairwise sigmoid estimator on the
empirical observations and on generated observations, then uses Optuna to match
summaries of those inferred coefficients. It does not match individual empirical
interaction identities, and it does not inherit the published model's empirical
accuracy merely by using its architecture.

This changed the fitting objective rather than exactly repeating the earlier
EXP-065 protocol. Historical RUN-175, RUN-178 and RUN-180 inputs, sources and
outputs remain unchanged. The example configuration is not a calibrated TrpB
preset; the registered RUN-184 settings and results are preserved separately in
EXP-065.

### Architecture model and controls

The generator retains profile-conditioned additive preferences, full centered
20 × 20 pair matrices, exact mutation-RMS scaling and a bounded sigmoid:

\[
f(x)=\ell+s\,\sigma\!\left[b+\beta\{A(x)+rG_\alpha(x)\}/\sigma_A\right].
\]

Here the additive field has uniform-background mean zero and the pure pair field
has zero row and column means. The two fields have equal uniform-background
single-substitution RMS before multiplication by \(r\).

| Parameter | Meaning |
|---|---|
| `r` | Pairwise/additive mutation-RMS ratio. Changes relative interaction strength and can change signs, peaks and paths. |
| `concentration` (\(\alpha\)) | Concentration of pair-table effects. Smaller values yield more uneven effect sizes before centering and normalization. |
| `steepness` (\(\beta\)) | How strongly latent differences pass through the sigmoid. |
| `position` (\(b\)) | Sigmoid position relative to the uniform-background mean latent score. |
| `lower`, `span` (\(\ell,s\)) | Assay-scale nuisance parameters; span is positive. |

The concentration construction is
\(V\sim\mathrm{Gamma}(\alpha,\mathrm{rate}=\alpha)\),
\(Z\sim N(0,1)\), \(J=\sqrt VZ\), independently for each raw table entry.
It has \(E[V]=1\), so alpha changes the spread of effect sizes rather than their
ensemble variance. Large alpha approaches the Gaussian generator. The exact
Gaussian generator is also retained as an explicit backend limit. Table centering
changes the marginal distribution and induces dependence; alpha is a sampling
parameter, not an exact count of nonzero or strong coefficients.

Degree stays fixed during calibration. At four sites degree six gives all six
pairs. Neither degree nor alpha creates latent higher-order terms. The sigmoid
can produce higher-order structure on the observed-fitness scale.

Positive steepness and span preserve exact intrinsic fitness ordering at fixed
\(r\) and alpha. They change saturation, curvature and effect magnitudes, not
which exact noiseless mutations improve fitness. Finite resolution and observation
noise can change apparent ties and paths. The latent uniform-product variance
ratio is \(V_2/V_1=r^2/2\); it is not the observed epistatic fraction.

The example bounds are initial search bounds, not established biological ranges.
`lower_bounds` and `span_bounds` apply after one calibration-only affine fitness
standardization. The result records the conversion back to native assay units.
The profile is fixed input; a global positive profile multiplier would cancel
under mutation-RMS normalization and is therefore not another independent knob.

### Architecture fitting rule and interpretation

1. Freeze the empirical calibration rows, profile, fold assignment, inference
   settings, observation assumption, optimizer budget and seed lists before
   fitting. Genotype replicates stay together in cross-validation. The same
   estimator and fold rule are applied to real and synthetic observations.
2. Infer additive and pairwise effects through a sigmoid using the Park–Metzger–
   Thornton SSE-plus-L1 objective. Cross-validation selects the penalty. Transform
   effects to the uniform-product zero-sum gauge before comparing architecture.
   This Python adaptation follows the authors' published code; it is not their
   R executable, and it does not guarantee the same local optimum.
3. Match the fixed 18-summary vector defined below. Compare inferred
   synthetic coefficients with inferred empirical coefficients, so regularization
   acts on both sides. Raw generated coefficients are not substituted for fitted
   synthetic coefficients.
4. Bootstrap whole genotype groups in the calibration observations once to define
   fixed summary scales. Keep all their replicate measurements together and retain
   their original CV fold identities in every bootstrap refit.
   Optuna's multivariate TPE minimizes the squared summary discrepancy after
   scaling. Every parameter proposal uses the same calibration generator seeds;
   the optimizer never selects a favorable seed. The uncertainty weights remain
   fixed across proposals and do not reward simulated overdispersion.
5. Re-evaluate the selected setting with fresh generator seeds. If separate
   confirmation observations are supplied, evaluate them without refitting the
   generator parameters or uncertainty weights. Inspect diagnostics excluded from
   the objective separately. Fresh generator seeds test realization variability;
   they are not new biological replicates or a substitute for held-out data.

The 18 summaries are:

- `log1p` additive and pairwise latent variances, and the pairwise share of their
  total. Zero-sum effects use the correct \(q^{-1}\) and \(q^{-2}\) variance
  weights for the respective orders.
- Absolute coefficient quantiles (50th, 90th and 99th percentiles), divided by
  the RMS for that order, separately for additive and pairwise effects.
- Fraction of squared coefficient energy held by the strongest 10% of entries,
  separately for each order. This is **not** the paper's T90 predictive-retention
  and refitting statistic.
- Observed fitness 10th, 50th and 90th percentiles; fitted lower and upper sigmoid
  limits; and fractions of sampled sequences with inferred sigmoid occupancy
  below 0.05 and above 0.95.

Each auxiliary fit uses every finite replicate measurement individually; values
are not reduced to genotype means for fitting. The outer loss is

\[
\mathcal L(\theta)=\frac1{18}\sum_j
\frac{\left(\bar s_{\mathrm{sim},j}(\theta)-s_{\mathrm{emp},j}\right)^2}
{\widehat{\operatorname{Var}}_{\mathrm{bootstrap}}(s_{\mathrm{emp},j})+\epsilon^2}.
\]

The fixed floor \(\epsilon=0.05\) is in each summary's normalized units. This is a
declared diagonal distance metric, not a likelihood, calibrated uncertainty model
or posterior. Summary correlations are not modeled. Trial histories include a
descriptive near-best region defined by loss within one unit of the best trial;
this is not a confidence or credible region.

This is **indirect inference**, using a published inference model as the auxiliary
estimator. The Gamma mixing distribution is an additional declared generative
assumption. L1 regularization does not establish that biological coefficients
follow this, or any other, coefficient sampling distribution. Matching a compact
summary vector also does not establish coefficient identity, parameter
identifiability, global optimality or biological validity at long sequence lengths.

The default requires convergence of the selected empirical and synthetic auxiliary
fits. `require_inference_convergence=false` permits an explicitly budget-limited
engineering or sensitivity calculation; its output must not be reported as a
converged calibration. Optuna's finite search budget never proves global optimality.
Keep the trial history and boundary hits when interpreting an apparently good fit.

`noise_mode="none"` adds no measurement error. `"replicate_residual"` is a
replicate-derived residual sensitivity and requires replicate measurements. It is
not the historical shared-stop observation model from EXP-065; that assumption
is not silently carried into this new objective. Residuals are centered within
each genotype, corrected by \(\sqrt{R/(R-1)}\), pooled, and sampled independently
with fixed draws within each generator seed. The pool uses calibration data only;
its magnitude is not an extra optimized parameter. This assumes a common residual
distribution and does not identify correlated or genotype-dependent assay errors.

### Architecture inputs and execution

Install the maintained package's calibration dependencies in the intended
scientific environment:

```bash
python -m pip install './external/Adapt-Env[calibration]'
```

The calibration input is an NPZ containing:

- `sequences`: an integer array of shape N × L, with residues encoded 0–19 in
  `ACDEFGHIKLMNPQRSTVWY` order. The column/site mapping must be fixed and archived.
- `fitness`: an N-vector, or N × R matrix of replicate measurements. Missing
  replicate entries may be NaN; infinite measurements are invalid.
- `additive_matrix`: a finite L × 20 profile-derived effect matrix in the same
  alphabet and site order. It must have nonzero within-site variation overall.
- `folds` (optional): one fold identifier per sequence, to preserve a predefined
  genotype-level split.

A separate confirmation NPZ uses the same schema, but may omit `additive_matrix`.
If supplied, its profile must equal the calibration profile exactly. Extra NPZ
keys and unknown JSON settings are rejected; pickled input arrays are never loaded.
Confirmation genotypes must be disjoint from calibration genotypes. Replicate
holdout on the same genotypes is not implemented as independent confirmation by
this interface. Without a confirmation NPZ, the output explicitly identifies the
assessment as fresh seeds on calibration genotypes only.

```bash
python scripts/calibrate_architecture.py \
  --input /absolute/path/calibration.npz \
  --config docs/calibration/example.json \
  --confirmation /absolute/path/confirmation.npz \
  --validate-only

python scripts/calibrate_architecture.py \
  --input /absolute/path/calibration.npz \
  --config docs/calibration/example.json \
  --confirmation /absolute/path/confirmation.npz \
  --output /absolute/path/new-output-directory
```

Validation checks structure and configuration, not scientific adequacy. It does
not fit or create output files. The execution command refuses to overwrite a
nonempty directory and does not resume old studies. Input hashes, resolved
settings, fitted summaries, trial records and confirmation information accompany
the output; `cli_inputs.json` additionally records the exact supplied file bytes
and command arguments. `source_snapshot/` captures the package source and data
files alongside their SHA-256 hashes in `provenance.json`, including changes in
an uncommitted checkout. Failure artifacts retain the stage and available
inference diagnostics; a failed synthetic inference marks that trial failed and
allows subsequent Optuna trials to continue. A failed empirical inference stops
calibration rather than constructing a target from an unresolved fit.

The unfitted diagnostics are sign and reciprocal-sign epistasis on complete
observed faces, a pooled fixed-orientation diminishing-returns slope, expected
uniform-uphill steps and local peak fraction **on the observed panel's induced
graph**. Missing genotypes and edges are not invented. Sparse-panel paths and
peaks are not whole-landscape searchability estimates; the historical EXP-065
whole-library diagnostic pipeline remains separate.

For TrpB, the archived profile is
`data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/outputs/additive_matrix.npy`.
Existing measured panels and geometry are under
`data/processed/interaction-strength-exp065-inputs/TrpB/`.
These are provenance sources, not a prepared input for this new protocol.
The old calibration panels contain cubes and repeated appearances of genotypes;
they must not be flattened into purported independent observations or a new
held-out split without recording that choice. The registered architecture input preparation and frozen RUN-184 bundle document
its separate genotype partition. These inputs are not substituted into the restored
RUN-180 cube protocol.

### Architecture computational scope

Online oracle evaluation remains sparse in site pairs and does not enumerate
\(20^L\) sequences. Auxiliary inference is a separate offline problem: fitting
all pairwise coefficients uses \(O(L^2q^2)\) coefficients, and its sparse design
uses \(O(N(1+L+P)+Lq+Pq^2)\) storage for \(P\) candidate pairs. Cross-validation,
bootstrap refits and repeated synthetic inference can dominate calibration cost.
The new objective is more expensive than the old quantile/local-SD objective.
The example's 100 trials × 4 calibration seeds × (4 penalties × 10 folds + one
selected-penalty refit) entail **16,400 auxiliary coefficient fits** in the outer
search alone, if every trial completes. Empirical, bootstrap and confirmation
fits are additional. This count is of fitted models, not solver iterations or
potential curvature resets within a fit.
No five-minute runtime or practical all-pairs inference at hundreds of sites is
implied by tractable online scoring. Small, declared calibration panels and a
timing pilot are needed to set a scientific execution budget.

## Sources

- [Park, Metzger & Thornton (2024), *The simplicity of protein sequence–function
  relationships*](https://doi.org/10.1038/s41467-024-51895-5).
- [Authors' RFA implementation, pinned source](https://github.com/JoeThorntonLab/RFA/blob/7a1a8725315a027d9563f6d9b12749174cffefc0/Scripts/RFA-functions.R).
- [Gouriéroux, Monfort & Renault (1993), *Indirect inference*](https://doi.org/10.1002/jae.3950080507).
- [Optuna multivariate TPE sampler](https://optuna.readthedocs.io/en/stable/reference/samplers/generated/optuna.samplers.TPESampler.html).
