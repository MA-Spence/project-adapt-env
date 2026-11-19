# EXP-042: PHOT nested-min cheaper fitting protocol equivalence benchmark

## Status

- planned

## Scientific lineage

- Aims: AIM-001
- Hypotheses: HYP-004, HYP-009

## Question

Can deterministic summary-loss optimization, surrogate-assisted fitting, exact-descriptive-gated fitting, and local ABC initialized from cheap fits recover PHOT nested min-coupled residual landscape behavior similar to the previously computed full ABC fit, well enough to justify scaling the protocol to complete GB1, TrpB, and related combinatorial landscapes?

## Pre-experiment prediction

At least one cheaper protocol will approach the full RUN-098 nested-min ABC comparator on held-out PHOT graphlet distance and component summaries, but parameter equality is not required and failures will be judged on graphlet spectral, epistasis, local ruggedness, accessibility, and diminishing-return summaries.

## Rationale

Full SMC-ABC is too expensive to use as the default fitting protocol for complete 160,000-variant landscapes. This experiment keeps the latest PHOT nested min-coupled sequential residual model family fixed and tests cheaper search protocols against the already computed full ABC comparator before using those protocols on larger complete empirical landscapes.

## Experimental design

- Keep the latest PHOT nested min-coupled, sequential-residual model family from EXP-039/RUN-098 fixed. This experiment does not rerun the full ABC comparator.
- Use RUN-077 as the frozen PHOT prepared-assay and graphlet source, and RUN-098 selected-family distances as the full ABC benchmark.
- Fit the same second-stage residual landscape parameterization with four cheaper protocols:
  - deterministic summary-loss search over the full graphlet summary vector;
  - surrogate-assisted search using a polynomial ridge surrogate trained on real simulator evaluations;
  - exact-descriptive-gated search using the local energy summaries identified by EXP-041/RES-040, followed by full-summary promotion;
  - local ABC initialized around the best cheap-fit candidate with narrowed priors and fewer particles.
- Select the best candidate per protocol and family on calibration graphlets, then evaluate selected candidates on calibration, test, and all graphlets.
- Compare each selected cheaper fit against the full RUN-098 ABC selected fit from the same family and against the best RUN-098 family overall.
- Primary readouts are held-out graphlet summary distance and component guardrails. Parameter equality is not required.
- Secondary readouts include graphlet spectral allocation, sign and reciprocal sign epistasis, magnitude epistasis, diminishing-return epistasis, accessibility, local optima, random-walk autocorrelation, and module-energy localization.
- Runtime and candidate-evaluation counts are recorded to decide whether the protocol is plausibly scalable to complete GB1 and TrpB landscapes.
- Interpretation boundary: operational equivalence on PHOT would support trying the cheaper protocol at GB1/TrpB scale; it would not prove that the protocol is universally equivalent to full ABC.

## Inputs

- `data/processed/system-specific-current-model-refit-graphlet-falsifier-repeat/RUN-077.dvc`
- `data/processed/phot-nested-softmin-ensemble-capacity-diagnostic/RUN-098.dvc`
- `scripts/phot_nested_fit_protocol_benchmark.py`
- `scripts/phot_nested_softmin_ensemble_abc.py`
- `src/project_adapt_env/smc_abc.py`

## Configuration

- `config.yaml` defines the fixed nested-min model family, fitting protocol budgets, exact-gate summaries, equivalence margins, graphlet metrics, and quick-mode smoke-test overrides.

## Execution

- Run through `labproj submit EXP-042 --target lab-slurm --execute`.
- `run.sh` pulls only the required RUN-077 and RUN-098 DVC targets if missing, then runs `scripts/phot_nested_fit_protocol_benchmark.py`.

## Expected outputs

- summary.json
- preflight.json
- fit_protocol_candidates.csv
- selected_fit_protocol_candidates.csv
- fit_protocol_distances.csv
- fit_protocol_equivalence.csv
- graphlet_metrics.csv
- graphlet_summary.csv
- walsh_terms.csv
- gft_coefficients.csv
- module_energy_summary.csv
- abc_comparator_selected_family_distances.csv
- full_abc_comparator_selected_particles.csv

## Analysis plan

- Create an ANA record after execution that reports whether any cheaper protocol is operationally equivalent to the full ABC comparator under the predeclared margins.
- Treat mixed outcomes by protocol and family explicitly; do not collapse a single success into universal equivalence.
- Report whether exact-descriptive-gated fitting solves the RES-040 local energy failures without worsening broader graphlet, epistasis, accessibility, or spectral summaries.

## Completion criteria

- [x] Config defined
- [x] Run script defined
- [x] Environment recorded
- [ ] Outputs generated
- [ ] Important outputs DVC-tracked
- [ ] Analysis record created
- [ ] Result record created
- [ ] Hypothesis updated
- [ ] PROJECT_STATE.md updated

## Post-experiment notes

- Add notes here after execution.
