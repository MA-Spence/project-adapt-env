# CLAIM-004: Thermodynamic GNK landscapes support conditional aggregate and adaptive-search resemblance to TrpB

- Status: completed
- Created: 2026-10-08
- Owner: matthew-spence
- Aims: AIM-001
- Relevant hypotheses: HYP-004, HYP-005
- Direct evidence: ANA-069 / EXP-065 / RUN-175; context from ANA-066, ANA-067 and ANA-068.
- Linked result: RES-063 (2026-11-19 synthesis; original synthetic-prior claim unchanged).
- Supersedes: none. This is a scoped new claim about the current model, not a retrospective change to older claims.

## Claim

The current thermodynamic GNK model admits parameter regimes that, together with
an explicitly specified observation model, produce useful approximate joint
resemblance to the measured four-site TrpB landscape in aggregate epistasis,
diminishing returns and adaptive-search behaviour. This supports its use as an
empirically informed synthetic benchmark environment. “Works / behaves well” here
means recovering several relevant global properties approximately and producing
nontrivial search, not reproducing TrpB's sequence-specific interactions or
matching every statistic.

The strongest joint example highlighted by RUN-175 uses interaction ratio r=0.5,
representation rank 8, order-energy weights (0.35, 0.45, 0.20), reference rank 0.5,
and the fixed shared-stop observation sensitivity. It is a post-sweep illustrative
regime, not a demonstrated optimum or uniquely identified biological mechanism.

## Quantitative evidence

Values below are from ANA-069's primary resolution: epsilon = 0.001 times the
frozen empirical calibration range. Model medians and 10–90% ranges use all eight
fresh seeds; no seeds are discarded. These ranges describe generation variability,
not confidence intervals or total measurement/parameter uncertainty.

| Statistic | Empirical | Model median | Model 10–90% seed range |
|---|---:|---:|---:|
| Expected uniform-uphill substitutions | 8.5887 | 7.4591 | 6.2828–8.2901 |
| Upper-1% walk success probability | 0.6498 | 0.7144 | 0.5194–0.8540 |
| Equal-cube WHT epistatic fraction | 0.6876 | 0.6295 | 0.5772–0.7093 |
| Equal-cube WHT order ≥3 fraction | 0.3156 | 0.2733 | 0.2154–0.3717 |
| Sign-epistatic face fraction | 0.6510 | 0.6562 | 0.5781–0.6771 |
| Median fixed-direction DR slope | -0.9306 | -0.8985 | -1.0623–-0.6985 |
| Mean local fitness SD | 0.0379 | 0.0304 | 0.0203–0.0418 |
| Pooled WHT epistatic fraction | 0.3297 | 0.4838 | 0.3352–0.7014 |
| Pooled WHT order ≥3 fraction | 0.0862 | 0.1760 | 0.0851–0.2850 |

The whole-library graph contains 159,129 observed TrpB genotypes of 160,000
possible four-site sequences. Synthetic maps use the same mask. Mean cube WHT
fractions weight eligible local cubes equally; pooled fractions weight the
combined power and reveal a larger mismatch. Neither weighting is substituted
for the other. The mean local SD and graph outcomes also retain visible differences.

Uniform-uphill walks choose among improving neighbors uniformly. Walks start
uniformly over retained nodes, stop when no improvement exceeds epsilon, and
are evaluated by exact dynamic programming. Top 1% is defined separately for
each landscape. Thus this compares global search properties, not arrival at the
same empirical genotype.

## What the experiment tests

RUN-175 varied r over 0, 0.5, 1, 2, 4 and 8, pairing the same synthetic fields
within each seed. It retained 16 reused and eight fresh seeds per protein, fitted
only the global readout to the original fitness quantile, local-SD and floor-mass
objective, and held the observation donor pools fixed. WHT, DR and graph outcomes
were not in that objective. Nevertheless, choosing a regime after inspecting the
sweep makes this an exploratory demonstration, not independent validation of
a prospectively selected final model.

Calibration/confirmation cube splits were preserved. Whole-library graphs
include calibration nodes and previously inspected empirical data; fresh seeds
test generation variability, not independent biological replication. The finite
optimizer budget did not establish convergence or global optimality.

## Interpretation and constraints

1. The claim concerns the **generator plus a stated observation model**. The bare
   noiseless model does not show the same TrpB agreement. At noiseless r=1, walk
   length is close (8.96 versus 8.59) but top-1% success is 99.5% versus 65.0%.
2. Shared-stop noise already produces substantial apparent epistasis at r=0
   (equal-cube fraction 0.670, versus empirical 0.688). The comparison does not
   establish that latent interactions, or a thermodynamic interpretation, are
   uniquely necessary to explain the observed statistics. Shared-stop is a
   sensitivity assumption, not an identified assay-error mechanism.
3. Seed variability is material. At r=4 under shared-stop observations, top-1%
   success is 70.1% on the 16 reused seeds but 87.5% on the eight fresh seeds.
   A selected favourable realization is not evidence for the entire ensemble.
4. Pooled high-order WHT power remains about twice empirical at the highlighted
   r=0.5 setting (0.176 versus 0.086). The claim is approximate joint resemblance,
   not comprehensive statistical equivalence. Neutral networks are resolution-
   and floor-dependent and are not validated as conserved functional networks.
5. ANA-067 finds appreciable discrepancies across other assays. TrpB supports
   a focused case study, not universal cross-protein realism. Four-site evidence
   does not establish biological resemblance at arbitrary length. Procedural
   scalability and empirical realism are distinct properties.
6. The evidence does not establish transfer of ML algorithm rankings to real
   experiments, physical free-energy inference, or recovery of specific epistatic
   interaction identities. No numerical success gates are used for this claim.

## Implementation covered

The intrinsic model is `thermodynamic_gnk_v2`, retained in Adapt-Env 0.4.0 after
removal of the archived model families (local dev commit `37d23c64b15a26cafa41232940feabc218505cd5`). The new preset
`trpb_run175_r05_shared_stop_intrinsic` records the exact RUN-175 intrinsic readout.
It does not silently add the observation layer. The generic default and a preset
with a different length do not inherit this TrpB evidence automatically.

The public service's optional Gaussian draws are not the research shared-stop
model. Reproducing the observed comparison requires the frozen observation
procedure and inputs in RUN-175. This claim is not evidence that an unconfigured
public API call generates the observed TrpB statistics.

## Figures and executable analysis

- [ANA-069 analysis](../analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/analysis.md).
- [Executed notebook](../analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/notebooks/analysis.ipynb).
- [Figure 1: joint walk outcomes](../analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/figures/01_joint_search_tradeoffs.pdf).
- [Figure 3: TrpB observation sensitivity](../analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/figures/03_trpb_observation.pdf).
- Figures 2, 4 and 5 in that analysis retain aggregate, resolution, neutrality,
  calibration and denominator-availability comparisons. Pooled summaries are
  retained in the numerical tables above and in `aggregate_estimates.csv`.

## Data and source provenance

- [RUN-175 metadata](../experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/runs/RUN-175.yaml).
- [RUN-175 DVC pointer](../data/processed/aggregate-thermodynamic-gnk-calibration/RUN-175.dvc): directory MD5 `b33c8fdb785160b4df3f30d0a26c4711.dir`, 186 files, 151,681,986 bytes.
- Frozen source, input hashes, empirical data, fitted parameters, seed lists,
  generated fields and observation inputs are in that bundle. The ANA independently
  checks 24 objectives and 1,732 WHT rows, and reproduces all 390 prior graph anchors.
- [ANA-069 numerical table](../analyses/2026-10-08_ANA-069_run-175-interaction-strength-global-statistics-and-adaptive-walks/tables/aggregate_estimates.csv) is the direct source of this claim's table.
- [Consolidation verification](CLAIM-004_evidence/verification.json) records removal,
  current-package tests, documentation checks and replay against frozen RUN-175.

This claim is supported directly by the existing ANA/RUN record. The scientist
subsequently requested RES-063 as an evidence synthesis on 2026-11-19. That result
now supplies the RES-level link without a new experiment or extension of this
claim to the profile prior, ML behavior or biological transfer.

## Manuscript wording

“Using TrpB as an empirical reference, we identify synthetic parameter regimes
that approximately reproduce multiple aggregate landscape statistics and
adaptive-search behaviours under an explicit observation model.”

Report the observation assumption, the parameter selection procedure, all fresh
seeds, pooled as well as equal-cube spectra, and the cross-assay limitations.
TrpB is the principal empirical case study, not evidence of universal recovery.

## Literature context

[Johnston et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11317637/)
provides the combinatorial TrpB landscape. [Otwinowski et al. (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6094095/)
supports nonlinear latent-trait mapping as a global-epistasis framework.
[Reddy and Desai (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8057814/)
shows that global epistasis can also emerge from microscopic interactions.
These are conceptual foundations and interpretation constraints; the numerical
claim above rests on this project's recorded comparisons.
