# EXP-030: System-specific current-model refit and graphlet falsifier repeat

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-001, HYP-004, HYP-005

## Question

When PHOT and SPTN1 are fit as two separate current-model landscapes rather than pooled, do the RES-028 fired falsifiers still hold: wrong empirical spectral envelope and wrong embedded evolutionary/topological dynamics?

## Pre-experiment prediction

If the RES-028 failure was mainly caused by stale cached preprocessing and the RUN-073 pooling error, then assay-specific current-model ABC fits should produce posterior predictive graphlets whose WHT/GFT spectra, SE/RSE, local maxima, accessibility and random-walk autocorrelation match each system better than smooth-core, residual-only and GNK-style controls. If the fired falsifiers persist, the current residual-ruggedness mechanism remains structurally inadequate even after system-specific recalibration.

## Rationale

CLAIM-001 supersedes the original pooled/shared-regime HYP-001 interpretation: empirical-like statistics are recoverable, if at all, system by system. RUN-073 incorrectly pooled PHOT and SPTN1 into one shared ABC target and did not emit graphlet WHT/GFT statistics. This experiment repeats the current-model refit correctly as two independent landscape fits in one run, then tests the RES-028 fired graphlet breakpoints on each fitted system.

## Experimental design

- Fit `PHOT_CHLRE_Chen_2023` and `SPTN1_CHICK_Tsuboyama_2023_1TUD` as two separate `panel_size=1` ABC targets. There is no pooled posterior and no shared empirical target.
- Use the current causal model only: stability-margin readout, stability-routed pairwise/residual epistasis, and no legacy FGM dimensions, designed peak-distance terms, deterministic roughness, or exact-match high-order haplotype terms.
- Recompute the empirical preprocessing from raw ProteinGym for each assay, then evaluate dense complete Hamming graphlets from each freshly prepared assay.
- Record SE/RSE, local maxima, accessible paths, basin/topology summaries, random-walk autocorrelation, Walsh spectra, and graph Fourier spectra for empirical graphlets, fitted current models, smooth-core controls, residual-only diagnostics, GNK-style local controls, and nonlinear global-readout residual controls.
- Test the RES-028 falsifiers that fired: wrong spectral envelope and wrong embedded evolutionary/topological dynamics. Also retain the scale-confounding and smooth-core-only screens because they determine whether an apparent improvement is actually caused by nonlinear readout alone.

The SPTN1 assay is expected to be graphlet-limited. A preflight check against the RUN-073 prepared assay found complete 2D graphlets but no complete 3D or 4D graphlets even when the top feature pool was expanded to 120. Therefore high-order spectral falsification is expected to be testable for PHOT and may be explicitly recorded as non-testable for SPTN1 under the complete-graphlet criterion.

## Sampling budget

RUN-073 used a pooled two-assay empirical SMC run with 64 particles, four rounds, pool multipliers 6/8, and two simulation replicates per particle. Its empirical round best distances were approximately 10.7975, 10.9120, 10.7028, and 10.8228, with no improvement in the fourth round. The final posterior also remained broad with high effective sample size, so the extra round did not indicate useful contraction.

EXP-030 therefore uses 48 particles, three rounds, pool multipliers 5/6, and one simulation replicate per particle. This gives 816 particle attempts per assay and 1632 total attempts across the two independent systems, while removing the apparently wasteful fourth SMC round and duplicate simulation replicate.

## Inputs

- Raw ProteinGym v1.3 reference and substitutions data.
- MMSeqs/ColabFold alignments for the two assay targets.
- Current AdaptEnv source and project ABC calibration code.
- RES-028 and RES-029 as motivating result records only; legacy prepared outputs are not parameter sources.

## Configuration

- See `config.yaml`. The runner performs a preflight guard to reject pooled panels, legacy fitted parameters, stale deterministic branches, or non-current epistasis/readout semantics.

## Execution

- Use labproj submit for RUN generation.

## Expected outputs

- `current_model_preflight.json`
- `system_fit_summary.csv`
- per-system posterior and branch-validation files
- `system_graphlet_metrics.csv`
- `system_walsh_terms.csv`
- `system_gft_coefficients.csv`
- `system_model_distances.csv`
- `falsification_tests.json`
- `summary.json`

## Analysis plan

- Create an ANA record after execution that reads the per-system posterior fits and graphlet outputs directly. The analysis must state whether each RES-028 fired falsifier still fires for PHOT and SPTN1 separately, and must not average the two systems into a single pooled conclusion.

## Completion criteria

- [ ] Config committed
- [ ] Run script committed
- [ ] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
