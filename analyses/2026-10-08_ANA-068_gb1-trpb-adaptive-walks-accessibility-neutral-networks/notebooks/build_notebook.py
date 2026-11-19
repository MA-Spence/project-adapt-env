from pathlib import Path
import os,nbformat
P=Path(__file__).resolve().parent
cells=[]
def md(s):cells.append(nbformat.v4.new_markdown_cell(s))
def code(s):cells.append(nbformat.v4.new_code_cell(s))
md('''# ANA-068 — Adaptive walks, accessibility and neutral connectivity

**GB1 and TrpB · EXP-066 / RUN-172.** Reuse the frozen production fits to compare
whole measured single-amino-acid-substitution graphs with their synthetic counterparts.
This analysis adds dynamics and network structure to ANA-067's aggregate statistics.
It does not fit new parameters or test ML policies. Historical ANA-020 analysed an
older model; ANA-063/064 measured paths and neutral edges on small cubes.

The four-site spaces contain at most 160,000 genotypes and 6,080,000 undirected
edges, making exact graph calculations tractable. The model is always masked to
the same finite measured genotypes as its empirical comparison.

Run all cells to recompute all graph calculations and figures. A C++17 compiler
is required (macOS Command Line Tools or Linux g++); the small graph kernel is
compiled into a temporary directory. Python dependencies are pinned alongside
this notebook. See `../protocol.md` for the analysis plan.''')
code('''from pathlib import Path
import sys, json, time
from IPython.display import display
%matplotlib inline
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p/'labproject.yaml').exists())
ANA = ROOT/'analyses/2026-10-08_ANA-068_gb1-trpb-adaptive-walks-accessibility-neutral-networks'
sys.path.insert(0,str(ANA/'notebooks'))
import analysis_tools as tools
import plot_tools as plots
print('Data:', tools.DATA)
print('Independent kernel validation:', tools.validate_kernel())''')
md('''## Recompute the frozen comparisons

Sixteen existing confirmation landscape seeds are retained. GB1 uses r=4 and r=0
sigmoid fits. TrpB additionally uses r=4 independent-difference and shared-stop
calibrations, with intrinsic and observed maps kept separate. Observations are
**drawn once per genotype**, not redrawn whenever a walk revisits a sequence.

Dynamic programming gives exact expectations under uniform uphill and greedy
walks from every valid start. Reachability includes indirect paths and reversions.
The upper 1% target includes all fitness ties at its boundary; the best target is
the whole best measured-fitness plateau. Their realized sizes are reported.

This whole-library assessment includes calibration genotypes and noise donors.
Parameters are frozen, but the results are descriptive post-fit comparisons, not
an entirely new independent confirmation test. Computational assertions concern
mathematical identities and provenance, not scientific acceptance gates.''')
code('''started = time.perf_counter()
frames = tools.calculate()
print(f'Whole analysis completed in {time.perf_counter()-started:.1f} seconds')
display(frames['coverage'])
audit = json.loads((tools.DATA/'audit.json').read_text())
display({k:v for k,v in audit.items() if k not in ['manifest','source','inputs']})''')
md('''## Figure 1 — Adaptive trajectories

Uniform uphill walks choose uniformly among beneficial neighbours; greedy walks
choose the largest available fitness, breaking exact ties uniformly. Neither
rule asserts a specific effective population size or mutation rate. The x-axis
counts accepted substitutions, not time. Walks start uniformly over measured
nodes and remain at their terminal point after absorption.

Fitness uses the empirical calibration origin and range, so the plot preserves
absolute readout-scale differences. Shading is variation across 16 synthetic
seeds (10th–90th percentiles), not a confidence interval for biological truth.''')
code('plots.trajectories(frames)')
md('''## Figure 2 — Where do walks finish?

Ranks are measured within each landscape, with midranks for ties. Similar endpoint
ranks need not imply similar absolute fitness. The empirical curve is conditional
on the measured map, and each map contributes 2,048 walks per policy.''')
code('plots.endpoints(frames)')
md('''## Figure 3 — Accessibility is not the same as probability of success

A genotype can have a beneficial route to a target even if most adaptive walks
terminate elsewhere. These are exact reachability fractions and walk probabilities
on the measured graph. No path-length cap is applied to the dynamic programming.
The displayed target sets include ties and are specific to each landscape.''')
code('plots.accessibility(frames)')
md('''## Figure 4 — Neutral networks require a stated definition

Here an edge is near-neutral when |Δfitness| ≤ epsilon. We examine exact ties and
tolerances of 0.001 and 0.01 times the empirical calibration range. These choices
are sensitivity scales, not experimentally identified selection thresholds.

Connected near-neutral edges can **chain across a broad fitness range**. A large
component is therefore not automatically one molecular phenotype or one fitness
plateau. We separately analyse the induced above-floor subset: GB1 >0, TrpB >0.05383.
These frozen assay cutoffs do not establish biological viability.''')
code('plots.neutrality(frames)')
md('''## Figure 5 — Observation can change the apparent search landscape

Compare matched intrinsic and observed maps under each TrpB noise assumption.
The noise model is reused from RUN-172. Its effect on topology is itself an
analysis result; it is not proof of the true intrinsic biological topology.''')
code('plots.observation(frames)')
md('''## Figure 6 — Missingness and neutrality diagnostics

Removing unmeasured genotypes can eliminate routes or create apparent optima.
The mask sensitivity compares model graphs with identical target thresholds and
identical measured starting nodes, differing only in availability of intermediate
nodes. It does not impute empirical values or reveal the actual missing routes.

Allowing exact equal-fitness moves is distinct from strict uphill dynamics.
Finally, the span of a component relative to one allowed edge tolerance detects
chaining across fitness values.''')
code('plots.diagnostics(frames)')
md('''## Exact summaries and reproducibility

The empirical graph is one measured map. Synthetic medians and seed ranges
represent variation in this generator with fixed fitted parameters; they omit
measurement and parameter-estimation uncertainty. Adaptive walk simulations are
checked against exact dynamic-programming expectations. Every sampled walk
terminated before the step cap. Raw landscapes, endpoint samples, summaries and
input hashes are retained in the linked DVC bundle.''')
code('''primary = frames['graph_summary'].query('relative_resolution == 0.001')
metrics = ['uniform_steps','greedy_steps','accessible_top','uniform_hit_top','greedy_hit_top',
           'accessible_best','uniform_hit_best','no_beneficial_fraction','neutral_edge_fraction',
           'all_largest_fraction','above_floor_largest_fraction']
display(primary.groupby(['protein','case'])[metrics].median())
display(frames['mask_sensitivity'].groupby(['protein','metric']).difference.agg(['min','median','max']))
for relative, expected in audit['manifest'].items():
    assert tools.sha(tools.DATA/relative) == expected, relative
print('Verified',len(audit['manifest']),'data artifact hashes; six figures exported as PNG and PDF.')''')
md('''## Interpretation and literature

See `../analysis.md` for the figure-led interpretation of this execution. The
important scientific distinctions are route existence versus walk success,
relative versus absolute fitness, and operational near-neutral connectivity
versus a truly conserved biological phenotype. This is a four-site empirical
comparison, not evidence for statistical invariance with arbitrary length.

- [Wu et al. (2016), *Adaptation in protein fitness landscapes is facilitated by indirect paths*](https://doi.org/10.7554/eLife.16965).
- [Johnston et al. (2024), *A combinatorially complete epistatic fitness landscape in an enzyme active site*](https://doi.org/10.1073/pnas.2400439121).

No new RUN, RES or CLAIM is created and no scientific success gate is imposed.''')
nb=nbformat.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3 (ANA-068)','language':'python','name':'python3'},'language_info':{'name':'python'}})
tmp=P/'analysis.tmp.ipynb';nbformat.write(nb,tmp);os.replace(tmp,P/'analysis.ipynb')
paired=[]
for c in cells:paired.append('# %% [markdown]\n'+'\n'.join('# '+s for s in c.source.splitlines()) if c.cell_type=='markdown' else '# %%\n'+c.source.replace('%matplotlib inline',"get_ipython().run_line_magic('matplotlib','inline')"))
(P/'analysis.py').write_text('\n\n'.join(paired)+'\n')
print('Built',sum(c.cell_type=='code' for c in cells),'code cells')
