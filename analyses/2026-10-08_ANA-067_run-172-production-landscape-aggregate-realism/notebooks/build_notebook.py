"""Build the executable figure-led ANA-067 notebook and its paired text source."""
from pathlib import Path
import os,nbformat
HERE=Path(__file__).resolve().parent
cells=[]
def md(s):cells.append(nbformat.v4.new_markdown_cell(s))
def code(s):cells.append(nbformat.v4.new_code_cell(s))
md('''# ANA-067 — How empirical do the production landscapes look?

**EXP-066 / RUN-172.** A comparison of 63 measured assay landscapes with the public
Adapt-Env thermodynamic generator. The scientific target is useful aggregate
resemblance for a synthetic online ML environment—not recovery of specific
empirical interactions or a universal model of protein biology.

One global readout is fitted per assay and arm. Calibration uses only fitness
quantiles, mean within-cube standard deviation, and near-floor mass. **WHT,
diminishing returns, and topology are withheld.** Confirmation uses disjoint
measured genotypes and 16 fresh synthetic seeds. There are no scientific
acceptance gates. The model family has seen some of these assay classes during
development; this is not a wholly unseen-biology discovery test.

Run this notebook top to bottom, or run `python notebooks/execute.py` from this
ANA directory with `notebooks/requirements.txt` installed. Restore the linked
RUN-172 DVC directory first if absent. Every figure is generated below and
exported as both PNG and PDF.''')
code('''from pathlib import Path
import sys, json
from IPython.display import display
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / 'labproject.yaml').exists())
ANA = ROOT / 'analyses/2026-10-08_ANA-067_run-172-production-landscape-aggregate-realism'
sys.path.insert(0, str(ANA / 'notebooks'))
import analysis_tools as tools
import numpy as np
import pandas as pd
print('Raw data:', tools.DATA)
print('Figures:', tools.FIGURES)''')
md('''## 1. Audit the frozen experiment before interpreting it

The audit verifies all output and source hashes, source/API readout identities,
all selected objective values, disjoint calibration/confirmation genotypes and
fresh seeds. It independently recomputes the WHT with a centered butterfly
transform, fixed-direction slopes and sign epistasis from raw cube values.
Very small powers can make fractions sensitive to floating-point cancellation;
these discrepancies receive explicit propagated roundoff bounds and are saved.

Three final excluded FLIP inventory rows were omitted by the RUN's last CSV
flush. The analysis restores those rows in a **derived** inventory table. The
frozen RUN remains unchanged. Statistical availability is not evidence of fit.''')
code('''stats, fits, suitability, audit = tools.audit()
display(pd.Series(audit, name='audit'))
display(suitability.groupby('status').size().rename('inventory entries'))''')
md('''## 2. Estimands and uncertainty

The main numerical resolution is **0.001 × each assay's calibration fitness
range**, with 1e-12 and 0.01 as sensitivities. This is a numerical comparison
scale, not a measured experimental detection limit.

- Equal-cube epistasis averages fractions only where cube variance is resolved.
- Pooled epistasis is the ratio of summed nonadditive power to summed total power.
- DR slopes regress a fixed-direction mutation effect on background fitness;
  means can be dominated by a small denominator. Beneficial-mutation and median
  diagnostics remain separate.
- 95% percentile bootstrap intervals resample empirical **cubes**, or synthetic
  **seeds**. Repeated observation draws are not independent landscapes.

The bootstrap is conditional on selected cubes and fitted parameters; it does
not refit calibration or account for dependence among related assays. The 50
stability assays dominate the count, so family summaries are kept separate.
Model length means variable library coordinates, not whole-protein length.''')
code('''estimates, comparisons, marginal, bootstrap = tools.summarize(stats, fits, draws=1500)
primary = comparisons[comparisons.observation.eq('none')]
summary = primary.assign(absolute_error=lambda d: abs(d.difference)).groupby(['family','arm','metric']).absolute_error.median().unstack('metric')
display(summary[['epistatic_fraction','ho_fraction','slope','sign_epistasis_fraction']])
print('Full numerical tables are saved in', tools.TABLES)''')
md('''## Figure 1 — Coverage and what calibration actually recovered

The scope panel includes every inventory entry. Calibration objective and
held-out quantile error are different quantities and should not be read as a
train–test equality plot. Local-spread ratios expose failures that a plausible
marginal fitness histogram could hide.''')
code('tools.scope_calibration(suitability, marginal, fits)')
md('''## Figure 2 — Statistics deliberately withheld from calibration

Each point is one assay/arm, with dimension-matched empirical and synthetic
sampling. Faint uncertainty bars accompany the r=4 sigmoid. The diagonal is a
reference for resemblance, not a decision boundary. Two-site assays contribute
no evidence about order-three or order-four epistasis.''')
code('tools.global_statistics(comparisons)')
md('''## Figure 3 — All assays, including the mismatches

This atlas shows signed deviations for the primary r=4 sigmoid with no added
observation model. Rows are grouped by family and sorted by label, never by
fit quality. Grey higher-order cells are **unavailable**, not zero. The colour
ranges differ by metric; clipped values are marked and retained exactly in the
CSV. This is a descriptive view, not a composite biological-realism score.''')
code('tools.residual_atlas(comparisons)')
md('''## Figure 4 — What the fitness distributions and local spread look like

Five examples are the protocol's pre-production smoke systems; GB1 is added
explicitly to illustrate seed/split sensitivity. Every confirmation cube and
seed enters each curve. No visually favorable sequence or seed is chosen.
The synthetic comparison targets the sampled library geometry, not an
unmeasured full protein sequence space.''')
code('tools.examples(fits)')
md('''## Figure 5 — TrpB: distinguish intrinsic structure from the observation process

Observed predictions are filled points; the corresponding intrinsic maps are
hollow points. TrpB is the only production assay with a matched replicate/control
observation sensitivity here. A good observed match under shared-stop covariance
does not uniquely identify intrinsic epistasis or a physical thermodynamic
mechanism. No equivalent noise decomposition is available for the other assays.''')
code('tools.trpb_observation(estimates)')
md('''## Figure 6 — Diminishing returns needs a denominator and a sampling context

A two-site cube has only two backgrounds per mutation. Large slope means can
therefore reflect unstable ratios rather than a reproducible biological rate.
Show resolution dependence and median diagnostics alongside the requested mean;
neither replaces the other retrospectively.''')
code('tools.resolution_diagnostics(estimates, stats)')
md('''## Numerical interpretation tables

These tables retain continuous differences rather than introducing a pass/fail
criterion. Optimizer convergence and proximity to bounds qualify the finite
search; they do not certify an optimum or imply that no better model parameters
exist. The public r=0 and r=4 settings are both supplied in portable form.''')
code('''m = marginal[marginal.observation.eq('none') & marginal.split.eq('confirmation')]
display(m.groupby(['family','arm'])[['quantile_nrmse','sd_ratio','calibration_loss']].median())
t = estimates[estimates.dataset_id.eq('trpb_johnston2024_four_site_graphfla') & estimates.split.eq('confirmation') & estimates.relative_resolution.eq(.001)]
display(t[t.metric.isin(['epistatic_fraction','ho_fraction','slope','sd_fitness','near_floor_fraction'])].pivot_table(index=['observation','arm'],columns='metric',values='estimate'))
print('Finite-budget fits converged:', sum(s['optimizer_converged'] for s in fits), '/', len(fits))
print('Fits within 1% of at least one bound:', sum(min(s['relative_boundary_distance']) < .01 for s in fits))
backup = json.loads((tools.TABLES / 'data_bundle_audit.json').read_text())
assert backup['all_md5_match']
display(pd.Series(backup, name='DVC backup'))''')
md('''## Interpretation and limits

The results support a **parameterized synthetic family with useful empirical
regimes**, especially when assay observations are modeled explicitly. They do
not support treating the fixed r=4 production recipe as broadly interchangeable
with all empirical landscapes. Similar marginal fitness distributions are
insufficient when local variation, sign epistasis or diminishing returns differ.

The additive sigmoid control can express magnitude epistasis but cannot produce
intrinsic sign epistasis; the GNK component supplies real synthetic ruggedness.
The need for that capacity is distinct from demonstrating the correct frequency
and strength of ruggedness in every assay. The parameter search is deliberately
bounded and tests r=0 and r=4, so it cannot rule out useful intermediate regimes.

TrpB supplies a concrete observed-fit example, with a material observation-model
qualification. PHOT/HIS7 and the common shortfall in stability-panel local spread
are substantive boundaries on broad empirical resemblance. No claim of perfect
landscapes, fitted empirical interaction identities, or statistical invariance
at arbitrary length follows from these data.

The single-latent nonlinear description is consistent with [Otwinowski,
McCandlish & Plotkin (2018)](https://doi.org/10.1073/pnas.1804015115). Global
fitness-dependent effects can also arise from complex interactions, as shown by
[Reddy & Desai (2021)](https://doi.org/10.7554/eLife.64740); aggregate agreement
alone therefore does not identify a unique physical mechanism.

See `../analysis.md` for the concise figure-led interpretation, and `../tables/`
for exact assay-level estimates and uncertainty. These remain an ANA analysis;
no new RES or CLAIM is created.''')
nb=nbformat.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3 (ANA-067)','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.13'}})
tmp=HERE/'analysis.tmp.ipynb';nbformat.write(nb,tmp);os.replace(tmp,HERE/'analysis.ipynb')
text=[]
for c in cells:
 text.append('# %% [markdown]\n'+'\n'.join('# '+line if line else '#' for line in c.source.splitlines()) if c.cell_type=='markdown' else '# %%\n'+c.source)
(HERE/'analysis.py').write_text('\n\n'.join(text)+'\n')
print('Built notebook with',sum(c.cell_type=='code' for c in cells),'code cells')
