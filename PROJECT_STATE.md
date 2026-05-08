# Project state

- Status: active aim, hypotheses recorded
- Active aim: `AIM-001 Create-Model-Landscape`

## Current scientific state

- The parent project still has one active aim and had no project-local hypotheses before this review.
- The practical scientific state currently lives in `external/Adapt-Env`, which already contains a synthetic landscape implementation, a population-genetics layer, and a validation suite aimed at qualitative protein-like behavior.
- The current model claim that is best supported is narrow: AdaptEnv is a biologically informed synthetic environment for benchmarking iterative search methods, not an exact reconstruction of a real protein fitness landscape.

## What the current quick validations support

- `experiment_dfe_characterisation --quick` supports protein-like local DFE structure.
  - unconditioned consensus neutral fraction: `0.677`
  - conditioned consensus neutral fraction: `0.906`
  - beneficial fraction increases with distance from peak: supported
  - DFE left skew: supported
- `experiment_alignment_conditioning --quick` supports a real effect of conditioning, though not yet external biological validity.
  - conservation-sensitivity Pearson: conditioned `0.385`, unconditioned `-0.067`
  - strong-covariation pairs show higher epistasis than random pairs: `0.593` versus `0.241`
  - sampled natural-like sequences remain viable: `1.0`
- `experiment_global_epistasis --quick` supports negative global epistasis across tested FGM dimensionalities.
  - mean Pearson correlations by dimensions `2/4/8`: `-0.557`, `-0.579`, `-0.536`
  - the intended monotonic strengthening with dimensionality is not supported
- `experiment_wf_dynamics --quick` supports qualitative adaptation and diversity turnover but not a strong long-term dynamics match.
  - early diversity increase: supported
  - baseline deceleration target: failed
  - truncation selection outperforms proportional selection in the tested setup
- `experiment_synthetic_config_realism --quick` supports most qualitative preset checks but not all regime claims.
  - alignment sequences viable across presets: supported
  - negative global epistasis majority: supported
  - smooth-bottom / rugged-top majority: supported
  - high-throughput beating low-throughput: failed in the quick run
  - viral earlier improvement than bacterial: failed in the quick run

## Main scientific weaknesses at the current state

- Peak interpretation is still the most important weakness for benchmark use.
  - `experiment_peak_sensitivity_noise_regimes --quick`: noise-free peak retention failed for `0/4` presets, implying the issue is mostly topology rather than noise.
  - `experiment_peak_sensitivity_pairwise_modes --quick`: an internal peak-preserving pairwise mode improved no-noise peak distance in `4/4` presets, but that mode is not the public default.
- Many validation experiments depend on the stored `peak_sequence`, which is local-search derived rather than globally certified.
- Most experiment helpers validate alignment conditioning using internally generated synthetic alignments in `external/Adapt-Env/experiments/common.py`, so those tests establish internal consistency more strongly than external biological realism.
- Several key knobs are not yet behaving as fully interpretable scientific controls.
  - FGM dimensionality does not cleanly strengthen diminishing returns
  - Wright-Fisher deceleration is weaker than intended
  - some preset-level ranking claims are configuration-sensitive
- The current validation suite is largely qualitative and internal.
  - it lacks preregistered external targets drawn from empirical benchmark panels
  - it does not yet cleanly separate "biological realism" from "benchmark usefulness"

## Record update from this review

- Project-local hypotheses have now been recorded in `docs/hypotheses.md` and the `labproj` registry.
- `EXP-001` has been created under `HYP-001` to assemble a real ProteinGym DMS panel, fetch MMseqs alignments from assay wild-type sequences, and fit a shared synthetic regime against external distributional statistics.
- The next scientific work should test those hypotheses against both empirical landscape summaries and benchmark-specific criteria for active learning and optimization.
