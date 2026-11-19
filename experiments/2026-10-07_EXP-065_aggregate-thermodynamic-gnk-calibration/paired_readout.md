# EXP-065 follow-up: paired readouts and embedded mutational neighbourhoods

Authorized by the scientist on 2026-10-07 as **one additional production RUN
under EXP-065**, implementing the three comparisons proposed after ANA-063.
This is a prospective follow-up informed by RUN-164, not a reinterpretation of
its original confirmation results. The original protocol and config remain
available; no historical numerical acceptance gates apply.

## Scientific question

Does nonlinear transmission of a synthetic additive-plus-interaction trait
explain aggregate empirical patterns better than linear transmission when
measurement amplitude and sampling are comparable? How do the relationships
between WHT power, diminishing returns, and background fitness change with
operating position and sequence length?

TrpB is the primary assay because ANA-063 found consistently resolved WHT spectra
and a clearer preference for bare GNK than the composite. This assay choice is
informed by the earlier result; no across-assay generality is claimed. The
prepared input has no independent replicate errors. Resolution sensitivity is
therefore retained and is not represented as an inferred measurement-noise model.

## 1. Paired linear and sigmoid readouts

Generate the same normalized additive field A and rank-8 GNK field G for both
readouts. Fix latent order weights to (0.35, 0.45, 0.20), the mixture selected
for the earlier TrpB composite. For each ratio r in (0, 0.25, 1, 4), compare

`y_linear = b + c (A + rG)`

`y_sigmoid = b + c sigmoid(logit(p0) + beta (A + rG))`.

Each candidate uses the same analytic calibration of b and positive c: match the
mean within-cube SD and overall mean of the measured calibration cubes. These
two coefficients are subsequently frozen. No clipping or per-cube adjustment
is applied. This isolates shape beyond a shared amplitude-calibration rule;
neither b nor c is interpreted as an identified molecular parameter.

The sigmoid grid is beta=(0.125, 0.5, 2) and p0=(0.2, 0.5, 0.8). Retain the
best calibration sigmoid at **each** r, alongside its linear counterpart. The
primary evidence is the four paired contrasts at fixed latent r. The sigmoid
has additional curvature parameters; held-out comparisons assess their predictive
value. Forty candidates are evaluated, yielding eight readout/ratio arms.

Use 24 measured calibration cubes and 24 confirmation cubes, with no genotype
shared between them and no genotype previously used in either RUN-164 panel.
Draw 144 disjoint candidate cubes and retain the first 48 that avoid the prior
768 TrpB genotypes. Sampling uses observed anchors, four variable sites and
uniform alternative observed alleles, conditional on complete measurements.
The 128 measured genotypes exposed by the follow-up engineering smoke test are
also excluded from production sampling and saved in the exclusion input.
No activity or epistasis filtering is used. Fail explicitly if the draw cannot
provide the declared panel; do not silently alter sampling.

Use four calibration generator seeds and eight fresh confirmation seeds, with
eight cubes per seed. The empirical affine normalization, descriptor embedding,
and readout parameters are fitted on calibration data only. The ten-descriptor
energy distance is retained from RUN-164 for comparability, alongside all raw
descriptors. It is not a scientific success score. Resolution is 0.001 of the
empirical calibration range, with sensitivity at 1e-12 and 0.01. The same frozen
calibration embedding is used throughout that sensitivity comparison.

## 2. Relationships among statistics and operating position

For measured and generated cubes retain absolute WHT power at every order,
resolved WHT fractions, fitness mean and SD, same-mutation beneficial and
deleterious slopes, neutrality, sign epistasis, peaks and accessibility. Save
every oriented 0-to-1 edge with its background fitness and effect, allowing
mutation-specific relationships to be analysed without pooling identities.
Undefined slope/spectrum values remain unavailable, accompanied by counts.

Tabulate descriptor means, medians and availability within background-fitness
bins whose boundaries are quartiles of **empirical calibration cube means**.
These conditional tables and saved edges express relationships; they do not
replace full distributions or support inference about specific interactions.

For the calibration-selected sigmoid, vary p0 over (0.1, 0.5, 0.9) on the same
confirmation latent cubes. Keep beta, r, b and c fixed. This is an intervention
on operating position, not refitting or a new empirical confirmation model.
Interpret coupled changes in WHT, slopes and neutrality across this sweep.
Near-zero spectra and shared measurement-error bias limit biological inference.

## 3. Embed four-site libraries in longer landscapes

Freeze the sigmoid selected on calibration data and its linear partner at the
**same r**, including both arms' four-site affine coefficients. Embed four-site
libraries in length 32, 128 and 512 generators, using three new seeds per length.
Evaluate reference-allele rank 0.5 (matching the fit) and rank 1 (the earlier
scaling stress-test convention) separately. No parameter is renormalized or
selected using the length-transfer observations.

For each seed use two libraries under each of two site policies: four uniformly
sampled sites, and a uniformly sampled existing order-four GNK support. Sites are
fixed within each library. The second policy is an interaction-enriched sampling
control, not a claim to recover TrpB's site-selection process or contact map.

For each library draw four anchors at each Hamming distance 0, 4 and 16 from
the reference, plus four uniformly random anchors. Alternate alleles are nonself
and uniform. Pair anchors and alternatives across site policies, and pair all
genotypes exactly across readouts. Hamming distance describes the anchor; other
cube nodes may be closer or farther from the reference. The reference has the
prescribed occupancy under the unscaled sigmoid; it is not a certified biological
wild type, optimum, or independently selected functional protein.

Record shared-support counts to expose the dilution of direct interactions among
four uniformly chosen sites as total length grows. Retain latent cubes, model
fingerprints, reference rank, sites, anchors, alternate alleles and build times.
Evaluate 1,152 embedded latent cubes and both frozen readouts, reporting the same
statistics and resolution sensitivity as above. Stable edge RMS avoids squaring
underflow; exact float64 zero effects are reported explicitly. This compares
sampling ensembles and model behaviour, not unobserved whole-protein WHT.

## Evidence and interpretation

The main outputs are paired distribution comparisons and conditional curves,
with seed/library/cube identifiers retained for dependence-aware analysis.
The effective replication includes eight confirmation landscapes and three
landscapes per length/rank, not thousands of independent cube observations.
No single comparison is converted to a pass/fail decision.

An advantage for nonlinear readout after amplitude calibration would strengthen
its explanatory role for this assay. Similar performance would indicate weak
identification of readout curvature by these statistics. Residual discrepancies
in absolute and fractional WHT power or its relationship with mutation slopes
would help distinguish insufficient latent structure from an observation-scale
problem. Divergence between local and random backgrounds, or between site
policies, would delimit the meaning of empirical resemblance at large length.

## Literature

- Otwinowski, McCandlish & Plotkin (2018), [Inferring the shape of global epistasis](https://doi.org/10.1073/pnas.1804015115): nonlinear latent maps and separation of variance components.
- Reddy & Desai (2021), [Global epistasis emerges from a generic model of a complex trait](https://doi.org/10.7554/eLife.64740): diminishing returns need not uniquely identify a saturating physical readout.
- Sailer & Harms (2017), [Detecting High-Order Epistasis in Nonlinear Genotype-Phenotype Maps](https://doi.org/10.1534/genetics.116.195214): WHT interpretation depends on phenotype scale.
- Bloom et al. (2006), [Protein stability promotes evolvability](https://doi.org/10.1073/pnas.0510098103): stability margin affects local mutational tolerance.

## Execution

`run.sh` now defaults to this follow-up. Set `ADAPTENV_PROTOCOL=initial` to run
the preserved original `config.yaml` and original calibration script. The source
snapshot collected for RUN-164 also preserves its exact original launcher.
`LABPROJ_QUICK=1` selects separate engineering seeds and reduced dimensions;
quick outputs are labelled SMOKE_ONLY and are not empirical evidence for this run.
The full job remains within the EXP-065 CPU allocation: 2 CPUs, 8 GB, two hours.

The job writes all cube/edge/latent tables, calibration candidates, selected
parameters, input/source hashes and three figures in PNG/PDF. Figures show the
paired comparison, conditional relationships and embedded neighbourhoods. They
are descriptive previews; raw tables retain all resolution and rank-1 controls.
