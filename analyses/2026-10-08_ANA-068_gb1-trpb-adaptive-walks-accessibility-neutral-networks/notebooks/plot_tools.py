"""Publication-style figures for ANA-068; estimates are descriptive across seeds."""
import os,tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import analysis_tools as t
C=t.COLORS;L=t.LABELS
PRIMARY={'GB1':['empirical','r4_none','r0_none'],'TrpB':['empirical','r4_none','r4_independent_observed','r4_shared_observed']}

def style():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'figure.titlesize':14,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'savefig.dpi':240,'legend.frameon':False,'legend.fontsize':8})
def save(fig,stem,caption):
 fig.text(.015,.01,caption,fontsize=8,va='bottom',color='#555555')
 for suffix in ['.png','.pdf']:
  with tempfile.NamedTemporaryFile(suffix=suffix,delete=False) as f:p=Path(f.name)
  try:fig.savefig(p,facecolor='white',bbox_inches='tight');os.replace(p,t.FIGURES/(stem+suffix))
  finally:p.unlink(missing_ok=True)
 plt.show()
def envelope(ax,z,x,y,case):
 p=z.groupby(x)[y].agg(median='median',low=lambda a:a.quantile(.1),high=lambda a:a.quantile(.9))
 ax.plot(p.index,p['median'],color=C[case],lw=1.8,label=L[case]);
 if case!='empirical':ax.fill_between(p.index,p.low,p.high,color=C[case],alpha=.14,lw=0)

def trajectories(frames):
 style();fig,ax=plt.subplots(2,2,figsize=(11.5,8));fig.subplots_adjust(left=.085,right=.985,bottom=.15,top=.90,hspace=.4,wspace=.28)
 z=frames['trajectories']
 for r,protein in enumerate(['GB1','TrpB']):
  for col,policy in enumerate(['uniform','greedy']):
   a=ax[r,col]
   for case in PRIMARY[protein]:envelope(a,z[z.protein.eq(protein)&z.policy.eq(policy)&z.case.eq(case)],'step','mean_fitness',case)
   a.set(xlabel='Accepted substitutions (absorbed walks held at endpoint)',ylabel='Mean fitness / calibration affine scale',title=f'{chr(97+2*r+col)}  {protein} · '+('uniform uphill' if policy=='uniform' else 'greedy ascent'),xlim=(0,20));a.legend(loc='center right' if r else 'best')
 fig.suptitle('Whole-library adaptive trajectories on measured genotype graphs')
 save(fig,'01_adaptive_trajectories','2,048 uniformly sampled starts per map/policy. Curves: empirical walk mean or median synthetic walk mean over 16 landscape seeds; shading: 10–90% seed range.\nFitness uses the frozen empirical calibration origin and scale. TrpB observed maps contain one frozen noise realization per genotype; noise is not redrawn on visits.')

def endpoints(frames):
 style();fig,ax=plt.subplots(2,2,figsize=(11.5,8));fig.subplots_adjust(left=.085,right=.985,bottom=.14,top=.90,hspace=.4,wspace=.28)
 z=frames['endpoint_cdf']
 for r,protein in enumerate(['GB1','TrpB']):
  for col,policy in enumerate(['uniform','greedy']):
   a=ax[r,col]
   for case in PRIMARY[protein]:envelope(a,z[z.protein.eq(protein)&z.policy.eq(policy)&z.case.eq(case)],'rank','cdf',case)
   a.set(xlabel='Terminal fitness percentile within its own landscape',ylabel='Cumulative fraction of walks',xlim=(.5,1.005),ylim=(-.015,1.015),title=f'{chr(97+2*r+col)}  {protein} · {policy}');a.legend(loc='upper left')
 fig.suptitle('How high do adaptive walks finish?')
 save(fig,'02_walk_endpoints','Percentiles use midranks, preserving ties; all endpoints remain in the cumulative calculation (the plot focuses on the upper half).\nCurves summarize 2,048 walks per map; shading is the 10–90% range across 16 seeds. These are relative endpoint ranks, not recovery of absolute empirical fitness.')

def accessibility(frames):
 style();fig,ax=plt.subplots(2,2,figsize=(12.5,8));fig.subplots_adjust(left=.20,right=.985,bottom=.15,top=.86,hspace=.55,wspace=.30)
 z=frames['graph_summary'];z=z[z.relative_resolution.eq(.001)]
 markers=[('accessible_','D',-.18,'Any uphill path'),('uniform_hit_','o',0,'Uniform uphill walk'),('greedy_hit_','^',.18,'Greedy walk')]
 for r,protein in enumerate(['GB1','TrpB']):
  for j,target in enumerate(['top','best']):
   a=ax[r,j];cases=PRIMARY[protein]
   for y,case in enumerate(cases):
    w=z[z.protein.eq(protein)&z.case.eq(case)]
    for prefix,marker,offset,_ in markers:
     v=w[prefix+target];m=v.median();lo,hi=v.quantile([.1,.9]);a.errorbar(m,y+offset,xerr=[[m-lo],[hi-m]],fmt=marker,color=C[case],ms=5,capsize=2,lw=1)
   a.set(yticks=range(len(cases)),yticklabels=[L[c] for c in cases],xlim=(-.02,1.02),xlabel='Fraction of all measured starting genotypes',title=f'{chr(97+2*r+j)}  {protein} · '+('upper 1% fitness set' if target=='top' else 'best measured-fitness plateau'));a.invert_yaxis();a.grid(axis='x',alpha=.15)
   if j==1:a.tick_params(axis='y',labelleft=False)
 fig.legend(handles=[Line2D([],[],marker=m,ls='',color='#333333',label=lab) for _,m,_,lab in markers],loc='upper center',bbox_to_anchor=(.55,.945),ncol=3)
 fig.suptitle('A route can exist even when an adaptive walk rarely follows it',y=.99)
 save(fig,'03_accessibility','Exact graph reachability and dynamic-programming walk probabilities. No path-length restriction: indirect substitutions and reversions are allowed.\nBars: 10–90% seed range. Fitness thresholds include ties; targets are defined within each map. Missing empirical genotypes are excluded identically from model graphs.')

def neutrality(frames):
 style();fig,ax=plt.subplots(2,3,figsize=(14,8));fig.subplots_adjust(left=.065,right=.985,bottom=.16,top=.89,hspace=.40,wspace=.33)
 z=frames['graph_summary'].copy();z['resolution_index']=z.relative_resolution.map({0.:0,.001:1,.01:2})
 specs=[('neutral_edge_fraction','Near-neutral edges / measured edges'),('all_largest_fraction','Largest component / all genotypes'),('above_floor_largest_fraction','Largest component / above-floor genotypes')]
 for r,protein in enumerate(['GB1','TrpB']):
  for j,(metric,ylabel) in enumerate(specs):
   a=ax[r,j]
   for case in PRIMARY[protein]:envelope(a,z[z.protein.eq(protein)&z.case.eq(case)],'resolution_index',metric,case)
   a.set(xticks=[0,1,2],xticklabels=['Exact ties','0.001','0.01'],xlabel='Tolerance / calibration fitness range',ylabel=ylabel,ylim=(-.02,1.02),title=f'{chr(97+3*r+j)}  {protein}');a.grid(axis='y',alpha=.15)
   if j==0:a.legend(loc='best')
 fig.suptitle('Near-neutral connectivity depends on tolerance and assay-floor membership')
 save(fig,'04_neutral_connectivity','Components connect single substitutions with |Δfitness| ≤ tolerance. Near-neutral chains can span much more than one tolerance; these are not automatically equal-phenotype networks.\nAbove-floor cutoffs are frozen from RUN-172: GB1 > 0; TrpB > 0.05383. Shading: 10–90% across 16 seeds. An assay floor is not a viability threshold.')

def observation(frames):
 style();fig,ax=plt.subplots(2,2,figsize=(10.5,7.7));fig.subplots_adjust(left=.09,right=.985,bottom=.15,top=.84,hspace=.44,wspace=.3)
 z=frames['graph_summary'];z=z[z.protein.eq('TrpB')&z.relative_resolution.eq(.001)]
 specs=[('uniform_steps','Expected accepted substitutions'),('uniform_hit_top','Probability of reaching upper 1%'),('neutral_edge_fraction','Near-neutral edge fraction'),('above_floor_largest_fraction','Largest above-floor component fraction')]
 for panel,(a,(metric,label)) in enumerate(zip(ax.ravel(),specs)):
  emp=z[z.case.eq('empirical')][metric].iloc[0];a.axhline(emp,color=C['empirical'],lw=1.2,label='Empirical')
  for x,prefix,col in [(0,'r4_independent',C['r4_independent_observed']),(1,'r4_shared',C['r4_shared_observed'])]:
   vals=[]
   for delta,kind in [(-.10,'intrinsic'),(.10,'observed')]:
    v=z[z.case.eq(prefix+'_'+kind)][metric];m=v.median();lo,hi=v.quantile([.1,.9]);vals.append(m)
    a.errorbar(x+delta,m,yerr=[[m-lo],[hi-m]],fmt='o',color=col,mfc='white' if kind=='intrinsic' else col,ms=5,capsize=2)
   a.plot([x-.10,x+.10],vals,color=col,alpha=.6,lw=1)
  a.set(xticks=[0,1],xticklabels=['Independent difference','Shared-stop sensitivity'],ylabel=label,xlim=(-.5,1.5),title=f'{chr(97+panel)}  TrpB');a.grid(axis='y',alpha=.15)
  if metric!='uniform_steps':a.set_ylim(-.02,1.02)
 fig.legend(handles=[Line2D([],[],color='#222222',label='Empirical'),Line2D([],[],marker='o',ls='',mfc='white',color='#555555',label='Intrinsic map'),Line2D([],[],marker='o',ls='',color='#555555',label='Frozen observed map')],loc='upper center',bbox_to_anchor=(.5,.94),ncol=3)
 fig.suptitle('Observation assumptions change search behaviour and network connectivity',y=.99)
 save(fig,'05_observation_sensitivity','Pairs use identical latent seeds and paired noise uniforms. Error bars are 10–90% ranges across 16 seeds, not confidence intervals or independent noise-only uncertainty.\nThe model of observations is transferred from the RUN-172 calibration; this analysis does not identify the true assay-noise process.')

def diagnostics(frames):
 style();fig,ax=plt.subplots(2,2,figsize=(11.5,8));fig.subplots_adjust(left=.10,right=.985,bottom=.15,top=.90,hspace=.5,wspace=.35)
 masks=frames['mask_sensitivity']
 for i,protein in enumerate(['GB1','TrpB']):
  a=ax[0,i];z=masks[masks.protein.eq(protein)]
  for x,metric in enumerate(['accessible_top','uniform_hit_top','no_beneficial_fraction']):
   v=z[z.metric.eq(metric)].difference.to_numpy();a.scatter(x+np.linspace(-.09,.09,len(v)),v,color=C['r4_none'],s=18)
  a.axhline(0,color='#777777',lw=.8);a.set(xticks=[0,1,2],xticklabels=['Any path\nto upper 1%','Uniform-walk\nsuccess','No beneficial\nneighbour'],ylabel='Measured-mask minus full-model fraction',title=f'{chr(97+i)}  {protein} · missing-genotype sensitivity')
 z=frames['graph_summary'];exact=z[z.relative_resolution.eq(0)&z.case.isin(['empirical','r4_none'])]
 a=ax[1,0]
 for i,protein in enumerate(['GB1','TrpB']):
  for j,case in enumerate(['empirical','r4_none']):
   w=exact[exact.protein.eq(protein)&exact.case.eq(case)];x=i*3+j
   for offset,metric,marker in [(-.12,'accessible_best','o'),(.12,'nondecreasing_best','s')]:
    v=w[metric];m=v.median();lo,hi=v.quantile([.1,.9]);a.errorbar(x+offset,m,yerr=[[m-lo],[hi-m]],fmt=marker,color=C[case],capsize=2,ms=5)
 a.set(xticks=[0,1,3,4],xticklabels=['GB1\nempirical','GB1\nmodel','TrpB\nempirical','TrpB\nmodel'],ylim=(-.02,1.02),ylabel='Fraction reaching best plateau',title='c  Equal-fitness steps can connect routes')
 a.legend(handles=[Line2D([],[],marker='o',color='#555555',ls='',label='Strictly uphill'),Line2D([],[],marker='s',color='#555555',ls='',label='Nondecreasing')],loc='best')
 a=ax[1,1];primary=z[z.relative_resolution.eq(.001)]
 for protein,marker in [('GB1','o'),('TrpB','^')]:
  for case in PRIMARY[protein]:
   w=primary[primary.protein.eq(protein)&primary.case.eq(case)]
   a.scatter(w.all_largest_fraction,w.all_largest_span/w.epsilon,color=C[case],marker=marker,s=25,alpha=.7)
 a.set(yscale='symlog',xlabel='Largest near-neutral component fraction',ylabel='Component fitness span / edge tolerance',title='d  Near-neutral chains need not stay isofit',xlim=(-.02,1.02));a.axhline(1,color='#888888',ls='--',lw=.8)
 fig.suptitle('Interpretation depends on missingness, plateaus and the definition of neutrality')
 save(fig,'06_graph_diagnostics','a,b  Same model, target threshold and measured starting nodes; only availability of intermediate genotypes changes. Points: 16 seeds, r=4 no added noise.\nc  Exact equality permits neutral steps. d  Circles: GB1; triangles: TrpB; colours follow the other figures. Large spans warn against interpreting an entire component as one phenotype.')
