# ANA-028: RUN-071 empirical graphlet falsification review

## Purpose

Review the EXP-028 empirical graphlet falsification suite, including the six preregistered falsifiers, empirical spectral envelope, synthetic-control comparisons, nonlinear readout confounding, and implications for HYP-004 and HYP-005.

## Linked experiments/runs

- Experiments: EXP-028
- Runs: RUN-071

## Notebook record

- Primary notebook: notebooks/analysis.ipynb
- Paired text file: notebooks/analysis.py
- Kernel: python3

## Produced artifacts

- `metrics.json`
- `tables/run-071_key_metrics.md`
- `tables/run-071_falsification_tests.md`
- `tables/run-071_sampling_audit.md`
- `tables/run-071_empirical_assay_audit.md`
- `tables/run-071_empirical_observed_by_assay_dimension.md`
- `tables/run-071_primary_empirical_readout_components.md`
- `tables/run-071_global_readout_confounding.md`
- `tables/run-071_synthetic_empirical_best_by_model.md`
- `tables/run-071_fixed_scale_directionality.md`

## Inputs

- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/summary.json`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/falsification_tests.json`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_empirical_best_by_model.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_empirical_model_distances.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_graphlet_suite/module_parameter_sweep.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/synthetic_graphlet_suite/test_graphlet_metrics.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/empirical_graphlet_suite/empirical_graphlet_metrics.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/empirical_graphlet_suite/empirical_graphlet_summary.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/empirical_graphlet_suite/empirical_global_readout_confounding.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071/empirical_graphlet_suite/empirical_global_readout_summary.csv`
- `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071.dvc`
- `experiments/2026-05-15_EXP-028_empirical-graphlet-falsification-residual-ruggedness/config.yaml`
- `experiments/2026-05-15_EXP-028_empirical-graphlet-falsification-residual-ruggedness/runs/RUN-071.yaml`
- `results/RES-027_run-069-fixed-residual-scale-sweep-validates-scale-knob-but-not-peak-stage-unmasking/result.md`

## Analysis performed

- Audited RUN-071 completion, output volume, DVC tracking, synthetic graphlet
  sampling, empirical graphlet extraction, and the six preregistered falsifier
  summaries.
- Recomputed the primary empirical target from the latent `PHOT_CHLRE_Chen_2023`
  graphlets with dimension at least 4, matching the comparison target used by
  `synthetic_empirical_best_by_model.csv`.
- Compared the smooth core, fixed contact residual, shuffled-contact residual,
  random-position residual, GNK local null, RES-024 exact-match haplotype
  control, DCA pair override, and residual-only diagnostic against the primary
  empirical target.
- Checked whether the empirical graphlet spectrum satisfied the registered
  low-order-dominant plus measurable high-order-tail envelope.
- Checked whether additive-plus-polynomial global-readout correction removed
  high-order WHT/GFT signal from empirical graphlets.
- Interpreted each falsifier separately so that a positive scale-knob result is
  not mistaken for validation of the full procedural benchmark model.

## Main observations

- RUN-071 completed with scheduler exit code `0`, generated `60` synthetic
  graphlets and `56` empirical graphlets, and the collected output is
  DVC-tracked at
  `data/processed/empirical-graphlet-falsification-residual-ruggedness/RUN-071.dvc`.
- The primary empirical PHOT graphlets are strongly rugged: SE `0.5979`, RSE
  `0.2524`, mean local maxima `4.9583`, accessible-path fraction `0.0398`, and
  random-walk lag-1 autocorrelation `0.1744`.
- The smooth core alone is not an empirical fit. Its distance to the primary
  empirical target is `1.9990`, versus `0.8432` for the selected contact
  residual. In synthetic test graphlets it has zero SE/RSE and one local
  maximum.
- The fixed contact residual scale grid remains directional across RUN-071's
  two landscape seeds. As scale increases from `1` to `4` to `16`, test SE
  rises `0.1875 -> 0.4176 -> 0.4867`, RSE rises
  `0.0341 -> 0.1259 -> 0.1572`, local maxima rise
  `1.25 -> 1.82 -> 2.02`, WHT/GFT high-order fraction rises
  `0.00089 -> 0.00790 -> 0.04421`, and lag-1 autocorrelation falls
  `0.4397 -> 0.4250 -> 0.2961`.
- The primary empirical spectral envelope is not the preregistered target. PHOT
  graphlets have GFT low/mid/high fractions `0.6056`, `0.2472`, and `0.1473`,
  but WHT high-order fraction is `0.3944`, above the preregistered
  non-dominant high-order-tail ceiling of `0.35`. Dimension-6 PHOT graphlets
  are especially high-order heavy, with WHT high-order fraction `0.6782`.
- Matched controls do not equal the contact residual under the registered
  distance criterion: contact residual distance is `0.8432`, versus GNK
  `1.2414`, random-position `1.2853`, and shuffled-contact `1.3002`.
- The best overall match is nevertheless the residual-only diagnostic
  (`0.7387`), not an embedded procedural model. This is a structural warning:
  the empirical topology looks closer to an unembedded residual field than to
  the current smooth-core plus contact-residual generator.
- No embedded procedural candidate matches the empirical package. The selected
  contact residual underestimates empirical local maxima (`2.02` versus
  `4.96`) and WHT high-order fraction (`0.044` versus `0.394`). The
  residual-only diagnostic gets closer on SE/local maxima but has zero
  accessibility and negative lag-1 autocorrelation (`-0.197` versus empirical
  `0.174`).
- Global-readout correction does not explain away the empirical high-order
  signal. In primary PHOT graphlets, the additive-plus-polynomial fit is
  smoother, but the residual has WHT high-order fraction `0.9479`; the mean
  retained high-order fraction after readout correction is `4.8206`.

## Interpretation

RUN-071 is not a clean validation of the current model as a general empirical
benchmark generator. It rejects two simple alternative explanations that would
have made the residual module unnecessary: the smooth core alone does not match
empirical graphlets, and nonlinear/global-readout confounding does not remove
the high-order residual signal. It also strengthens the narrow HYP-005 result
from RUN-069 by showing that the fixed residual scale grid is directional across
two landscape seeds.

However, the falsification suite also breaks the strong procedural adequacy
claim. The primary empirical graphlets do not have the clean low-order-dominant
WHT spectrum with only a modest high-order tail that the model expected. More
importantly, no embedded procedural candidate matches empirical SE/RSE, local
maxima, accessibility, random-walk autocorrelation, and spectra as a package.
The closest model is a residual-only diagnostic, which is not an acceptable
procedural benchmark mechanism because it bypasses the smooth thermodynamic
core and fails accessibility/autocorrelation.

The scientifically accurate conclusion is therefore mixed but falsifying for
the strong claim: RUN-071 preserves the fixed residual-scale knob and rejects
smooth/readout-only explanations, but falsifies the claim that the current
contact-conditioned residual implementation is already an empirically validated
general benchmark landscape model.

## Result records created

- `RES-028`

## Hypothesis updates

- `HYP-005` is supported narrowly for fixed residual-scale directionality.
- `HYP-004` is weakened. RUN-071 falsifies the strong current implementation
  and benchmark-adequacy form of the hypothesis, but it does not by itself
  refute every adaptive-walk prediction in the broader HYP-004 text.
