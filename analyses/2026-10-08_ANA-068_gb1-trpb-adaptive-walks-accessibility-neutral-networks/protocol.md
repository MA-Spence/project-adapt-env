# ANA-068 analysis plan — frozen before whole-library calculation

Reuse EXP-066 / RUN-172 parameters without fitting any new parameter. Historical
ANA-020 concerns an older model; ANA-063/064 provide small-cube accessibility and
neutral-edge summaries, not the present whole-library comparison.

## Scope and data

GB1 and TrpB have four variable sites and a 20-letter alphabet: 160,000 possible
nodes, 76 single-substitution neighbours and 6,080,000 undirected edges before
missing measurements. Use finite measured genotypes only and apply that exact
mask to each synthetic landscape. Never impute empirical missing values. Retain
all 16 RUN-172 confirmation landscape seeds. Assess the r=4 public sigmoid and
r=0 additive control for both proteins. For TrpB also assess independent-difference
and shared-stop calibrated r=4 maps, each intrinsically and after observation.
An observation realization is drawn once per genotype and held fixed throughout
walks; this describes a frozen assay map, not noise resampled at every query.
Noise realizations are paired across modes. This whole-library extension includes
previous calibration genotypes and noise donors; it is a descriptive post-fit
analysis, not a wholly independent confirmation dataset.

## Adaptive walks and accessibility

Analyse a directed edge when fitness gain exceeds epsilon. Compare uniform
choice among beneficial neighbours (equal-fixation walk) and greedy ascent
(best available fitness, uniformly breaking exact ties). Dynamic programming on
the fitness-ordered directed acyclic graph gives exact expected substitutions,
expected terminal fitness percentile and probabilities of reaching either the
upper 1% fitness set or the best measured-fitness plateau, from every valid start.
Targets include fitness ties; report their realized sizes. Accessibility is the
fraction of starting genotypes with any uphill path to the target, allowing
indirect paths and reversions. Report accessibility and actual walk hitting
probabilities separately. Uniform starts over all retained genotypes are primary.
Sample 2,048 walks per landscape and policy for trajectory and endpoint figures;
these are not independent landscape replications. Exact nondecreasing reachability
at epsilon=0 separately measures the effect of allowing truly equal-fitness steps.
No population sizes, mutation rates or calendar-time interpretation are assumed.

## Near-neutral connectivity

Compute connected components of the undirected graph with abs(delta fitness) <=
epsilon. Report neutral-edge fraction, component sizes, largest-component fraction
and its fitness span. Near-neutral edges may chain across a broad fitness range;
these components are not necessarily equal-phenotype networks. Repeat component
analysis on the induced subset above the frozen RUN-172 near-floor cutoff (GB1:
0; TrpB: 0.05383169008372466). 'Above floor' is an operational assay distinction,
not proof of biological viability. Isolate exact ties from tolerance-dependent
neutrality, and separate inactive-floor percolation from above-floor connectivity.

## Resolution, uncertainty and missingness

Use epsilon/calibration range = 0, 0.001 (primary) and 0.01. These are numerical
sensitivities, not measured selection coefficients or error bars. Empirical graph
statistics are exact conditional on the published measured map. Synthetic results
show variation across 16 seeds (median and 10th–90th percentiles), not a confidence
interval for a biological population. Existing uncertainty in measurements and
fitted parameters remains. Compare full versus measured-mask model graphs at
primary resolution for the common r=4 fit to quantify mask sensitivity. Missing
empirical states can hide routes or neighbours; graph peaks are conditional.

## Computation and verification

Implicit base-20 adjacency avoids dense NxN matrices. A small C++ kernel, compiled
into a temporary directory, provides O(N log N + E) graph calculations and O(N)
working storage. The notebook compiles and executes it, then generates all figures.
Independent Python tests on exhaustively tractable small graphs check dynamic
programming against recursive path enumeration and components against SciPy.
Compare public model readouts and selected saved RUN-172 cube values. Hash inputs
and frozen source. Generated durable data are stored under data/processed; ANA
contains source, tables, figures and an executed ipynb. No new RUN, RES, CLAIM,
scientific acceptance gate or model change is introduced.

## Scientific grounding

Wu et al. (2016), https://doi.org/10.7554/eLife.16965, motivates multiple uphill
walk rules and accessibility including indirect routes. Johnston et al. (2024),
https://doi.org/10.1073/pnas.2400439121, supplies the TrpB combinatorial landscape.
Neutral graph connectivity is operational here; fitness similarity alone does
not establish the same molecular phenotype or experimentally neutral selection.
