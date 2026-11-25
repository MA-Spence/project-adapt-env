# RUN-180 protocol rollback with Optuna

Authorized 2026-11-24. This is the current project TrpB calibration workflow,
selected by `ADAPTENV_PROTOCOL=optuna-readout` (the `run.sh` default). Implementation
and engineering checks do not register, execute or validate a new scientific RUN.
The previous RUN-180 and RUN-184 sources, configurations and results remain intact
in their archives. The architecture fitter remains an explicit experimental method.

## What is restored

Use exactly RUN-180's profile, Gaussian full 20 × 20 pair tables, six site pairs,
reference stream, normalized latent fields, calibration cubes, donor exclusions,
eight frozen observation draws, affine units, near-floor threshold and bounds.
`interaction_concentration=None` and `global_intercept=None` restore the Gaussian
and reference-anchored behavior; the optional Gamma feature remains available.

Keep the fixed ratio grid **r = 0, 0.5, 1, 2, 4, 8**, and separately calibrate
`none`, `independent_difference`, and `shared_stop_sensitivity`. Observation errors
are applied in native assay units before standardizing. The shared-stop arm remains
a sensitivity assumption, not an identified biological mechanism.

For each of these 18 settings, optimize the same four readout coordinates:
standardized offset, log standardized span, reference sigmoid position and log
steepness. Their bounds remain [-2,2], log([0.001,10]), [-14,8], log([0.05,8]).
The frozen calibration origin is -0.03522340992 and scale is 0.19396876295
(these rounded descriptions are not used in computation). The old large native
spans are therefore admissible again. No auxiliary coefficient estimator runs.

The objective is unchanged:

\[
L=\frac{\operatorname{mean}[(Q_{sim}-Q_{emp})^2]}{\operatorname{Var}(y_{emp})}
+\left(\frac{\overline{SD}_{sim}}{\overline{SD}_{emp}}-1\right)^2
+[F_{sim}(c)-F_{emp}(c)]^2.
\]

Q contains all 101 equally spaced quantiles including the minimum and maximum.
SD is calculated within each 16-vertex cube and then averaged over cubes and
observation draws. F(c) is the fraction at or below the frozen near-floor value.
There are 32 empirical calibration cubes and 64 synthetic cubes from eight seeds.
Topology, WHT spectra, sign epistasis and navigation never enter this objective.

## What changes

Replace differential evolution with serial Optuna 4.9.0 TPE, with 1,952 trials
per setting (35,136 total), including warm starts and 20 startup trials. Each
setting first evaluates its previous RUN-180 fitted coordinates, followed by
distinct old r=0/r=4 endpoint coordinates for the same observation mode. All
coordinates are re-evaluated at the current r and checked against the bounds.
The same-r loss must replay before optimization; the selected loss cannot be
worse than that candidate. This guarantee concerns calibration loss only.

Every r now receives a TPE search, including r=0 and r=4. Previously those two
received the endpoint DE searches; intermediate settings also considered the
endpoints. We reuse the archived endpoint candidates instead of performing new
DE searches. This is an optimizer replacement with warm starts, not a matched
test of optimizer efficiency. The old DE budget was 1,952 evaluations plus two
endpoint candidate evaluations for intermediate r.

Seeds are fixed and recorded; no generator seed is optimized. SQLite and CSV
histories retain all trials, components and candidate provenance. Completed
budget is recorded separately from convergence; global optimality is not claimed.
Output directories must be empty. Resuming a database alone is unsupported.

## Confirmation and provenance

Reuse all 16 historical confirmation seeds and eight fresh-to-RUN-180 seeds,
the 32 withheld empirical cubes, genotype mask, relative resolutions and complete
graph/cube diagnostics. The 24 seeds and withheld cubes have already been
inspected in prior work: this is a protocol rollback, not new independent
biological validation. Whole-library graph summaries are post-fit diagnostics;
the whole library is not a held-out panel.

The runner authenticates RUN-180's manifest and payloads and the original input
bundle. Before fitting it requires exact equality of regenerated calibration
fields and replay of every old loss component. It retains the baseline settings,
manifest, replay report, current source hashes, config and dependency versions.
Numerical dependencies are RUN-180's Python 3.11.15, NumPy 2.4.4, SciPy 1.17.1,
pandas 2.3.0, PyArrow 20.0.0 and PyYAML 6.0.3, plus Optuna 4.9.0.

Configuration: [optuna_readout.yaml](optuna_readout.yaml).
Engineering checks: [optuna_readout_validation.json](optuna_readout_validation.json).
Source manifest: [optuna_readout_source.json](optuna_readout_source.json).

The CLI's `--preflight` option only replays the old calibration and checks its
inputs. `--quick` is an eight-trial/two-confirmation-seed engineering smoke mode,
labelled `SMOKE_ONLY`; it is not the scientific protocol. Scientific execution
must be registered under this EXP through labproj when requested.
