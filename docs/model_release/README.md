# Unpublished model candidates

The current working candidate is [0.7.0rc1](0.7.0rc1/decision.md): the same profile
additive/full-pairwise/sigmoid architecture, with Gamma concentration and an
Optuna indirect-inference fitter. The new family/objective is not empirically
calibrated. The exact Gaussian limit remains the generic default until a finite
alpha is supplied or fitted. No biological finite-alpha default is invented.

[0.6.0rc1](0.6.0rc1/decision.md) retains its frozen Gaussian full-pairwise source
and RUN-180's matched historical calibration protocol.

[0.5.0rc1](0.5.0rc1/decision.md) is the prior ranked candidate and retains its
original manifest, replay and archive. Its named TrpB presets remain available
unchanged. An old candidate's lock verifies its archived source, not a newer
checkout; use `python docs/model_release/0.7.0rc1/lock.py --check` for this checkout.

None of these candidates is a published package release or a deployed service.
