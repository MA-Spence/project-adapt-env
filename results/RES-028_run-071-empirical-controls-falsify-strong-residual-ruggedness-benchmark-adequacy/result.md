# RES-028: RUN-071 empirical controls falsify strong residual ruggedness benchmark adequacy but preserve the scale knob

## Summary

RUN-071 rejects smooth-core-only and global-readout-only explanations and confirms that fixed residual scale is directional across the synthetic graphlet grid, but the primary empirical PHOT graphlets violate the preregistered modest high-order-tail spectral envelope and no embedded procedural candidate matches empirical spectra, SE/RSE, local maxima, accessibility, and random-walk statistics as a package. The result falsifies the strong current residual-ruggedness benchmark-adequacy claim while preserving only the narrow HYP-005 residual-scale knob.

## Generated from

- Analyses: ANA-028

## Relevant hypotheses

- Supports: HYP-005
- Weakens: HYP-004
- Refutes: None

## Evidence

- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/metrics.json`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_key_metrics.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_falsification_tests.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_sampling_audit.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_empirical_assay_audit.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_empirical_observed_by_assay_dimension.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_primary_empirical_readout_components.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_global_readout_confounding.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_synthetic_empirical_best_by_model.md`
- `analyses/2026-05-15_ANA-028_run-071-empirical-graphlet-falsification-review/tables/run-071_fixed_scale_directionality.md`
- `results/RES-028_run-071-empirical-controls-falsify-strong-residual-ruggedness-benchmark-adequacy/metrics.json`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/summary.json`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/falsification_tests.json`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_empirical_best_by_model.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_empirical_model_distances.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_graphlet_suite/module_parameter_sweep.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_graphlet_suite/test_graphlet_metrics.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/empirical_graphlet_suite/empirical_graphlet_metrics.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/empirical_graphlet_suite/empirical_global_readout_confounding.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071.dvc`
- `experiments/2026-05-15_EXP-028_empirical-graphlet-falsification-residual-ruggedness/config.yaml`
- `experiments/2026-05-15_EXP-028_empirical-graphlet-falsification-residual-ruggedness/runs/RUN-071.yaml`

## Interpretation

- `RUN-071` completed successfully with scheduler exit code `0`. The collected
  output is about `11.29` MiB, is DVC-tracked, and contains `60` synthetic
  graphlets, `1200` synthetic model-graphlet metric rows, `56` empirical
  graphlets, and `336` empirical metric rows.
- The experiment directly tests the six registered falsifiers. The outcome is
  mixed: four falsifiers do not fire, but two are scientifically material and
  break the strong current model-adequacy claim.
- The smooth-core-only falsifier does not fire. The smooth core is not a good
  empirical graphlet fit: its distance to the primary empirical PHOT target is
  `1.9990`, versus `0.8432` for the selected fixed contact residual. Smooth
  synthetic graphlets have SE/RSE `0` and one local maximum, while the primary
  empirical PHOT graphlets have SE `0.5979`, RSE `0.2524`, and mean local
  maxima `4.9583`.
- The no-residual-scale-effect falsifier does not fire. Across fixed contact
  residual scale `1 -> 4 -> 16`, test SE increases
  `0.1875 -> 0.4176 -> 0.4867`, RSE increases
  `0.0341 -> 0.1259 -> 0.1572`, local maxima increase
  `1.25 -> 1.82 -> 2.02`, WHT/GFT high-order fraction increases
  `0.00089 -> 0.00790 -> 0.04421`, and lag-1 random-walk autocorrelation falls
  `0.4397 -> 0.4250 -> 0.2961`. This reproduces and strengthens the narrow
  RUN-069 fixed-scale knob result across two landscape seeds.
- The wrong-spectral-envelope falsifier does fire. The primary PHOT empirical
  graphlets are not cleanly described as low-order dominated with only a
  modest high-order tail. Their GFT low/mid/high fractions are `0.6056`,
  `0.2472`, and `0.1473`, but their WHT high-order fraction is `0.3944`, above
  the preregistered non-dominant tail ceiling of `0.35`. Dimension-6 PHOT
  graphlets are especially high-order heavy, with WHT high-order fraction
  `0.6782`.
- The matched-control-equivalence falsifier does not fire as written. The
  contact residual distance to the empirical target is `0.8432`, better than
  GNK `1.2414`, random-position residual `1.2853`, and shuffled-contact
  residual `1.3002`. This means the registered GNK/random/shuffled controls do
  not match the empirical target as well as the contact residual in this run.
- The wrong-evolutionary-dynamics falsifier does fire in the broader procedural
  sense. The best overall match is the residual-only diagnostic
  (`distance = 0.7387`), not an embedded procedural model. The selected contact
  residual underestimates empirical local maxima (`2.02` versus `4.96`) and
  WHT high-order fraction (`0.044` versus `0.394`). The residual-only
  diagnostic better matches some ruggedness magnitudes but has zero
  accessibility and negative lag-1 autocorrelation (`-0.197` versus empirical
  `0.174`). No tested embedded candidate matches empirical spectra, SE/RSE,
  local maxima, accessibility, and random-walk statistics as a package.
- The scale-confounding falsifier does not fire. Additive-plus-polynomial
  global-readout correction does not explain away the high-order signal. In
  primary PHOT graphlets, the readout residual has WHT high-order fraction
  `0.9479`, and the mean retained high-order fraction after readout correction
  is `4.8206`.
- The core scientific interpretation is therefore not that the model is
  useless. The fixed residual-scale mechanism is real and directional, and the
  empirical graphlets are not explained by the smooth core or global readout
  alone. The falsified claim is stronger: the current smooth-core plus
  contact-conditioned residual implementation is not yet an empirically
  validated general benchmark landscape model.

## Effect on hypothesis

- `HYP-005` is supported narrowly. RUN-071 confirms that the fixed residual
  scale knob moves the intended local graphlet statistics directionally across
  the sampled synthetic grid.
- `HYP-004` is weakened. RUN-071 falsifies the strong current implementation
  and benchmark-adequacy form of the hypothesis: the empirical spectral
  envelope is wrong, and no embedded procedural candidate matches empirical
  spectra plus topology plus evolutionary dynamics as a package.
- This result does not formally refute every statement in `HYP-004`, because
  HYP-004 also contains broader adaptive-walk and peak-proximal claims that are
  only indirectly tested here. The honest conclusion is that the current
  implementation fails the direct empirical graphlet falsification suite for
  strong benchmark adequacy, while preserving a narrower residual-scale knob.

## Limitations

- The empirical primary target is dominated by complete PHOT graphlets at
  dimensions `4`, `5`, and `6`. SPTN1 and KRAS contribute only dimension-2
  control graphlets, which are useful controls but do not test the high-order
  empirical spectrum.
- The empirical graphlets are extracted from cached MAVE-NN/prepared assay
  outputs, not refit in this run. This is appropriate for reproducibility but
  means RUN-071 inherits those upstream preprocessing and latent-scale choices.
- The synthetic comparison uses dense local graphlets sampled from TEM-1
  procedural landscapes, while the empirical graphlets are from PHOT, SPTN1,
  and KRAS assays. This is a stringent generality test, not a one-protein
  calibrated empirical reconstruction.
- The synthetic graphlet budget is intentionally compact: `60` graphlets over
  two landscape seeds. This is enough for the falsification suite but not a
  broad robustness estimate over many proteins or many seeds.
- WHT and GFT summaries are equivalent for dense binary Hamming graphlets up to
  the order/eigenvalue grouping used here. The result does not test arbitrary
  non-Hamming graph topology.
- The residual-only diagnostic being closest to empirical graphlets is a
  structural warning, not a viable model. It bypasses the smooth thermodynamic
  core and therefore cannot be cited as procedural validation.

## Downstream use

- Use `RES-028` as evidence that fixed residual scale is a real procedural
  knob, consistent with `RES-027`, and that smooth/readout-only explanations
  are insufficient for the empirical graphlets tested here.
- Use `RES-028` as falsifying evidence against the strong claim that the
  current contact-conditioned residual module is already an empirically
  validated general benchmark landscape model.
- Do not cite `RUN-071` as proof of contact-specific biological causality. The
  registered matched controls do not beat contact residuals, but the
  residual-only diagnostic is closer overall and the embedded procedural model
  does not match the full empirical package.
- Do not describe the empirical target as a simple low-order-dominated spectrum
  with only a small high-order tail. The primary PHOT graphlets have a
  measurable and sometimes dominant WHT high-order component.
