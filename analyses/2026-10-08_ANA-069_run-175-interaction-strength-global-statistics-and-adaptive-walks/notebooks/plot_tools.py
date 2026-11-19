"""Figure functions for the RUN-175 paired interaction-strength sweep."""
from pathlib import Path
import tempfile,os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
RATIOS=[0,.5,1,2,4,8]
RCOL=['#440154','#414487','#2A788E','#22A884','#7AD151','#D8BD28']
MCOL={'none':'#0072B2','independent_difference':'#009E73','shared_stop_sensitivity':'#7B3294'}
MLAB={'none':'No added noise','independent_difference':'Independent difference','shared_stop_sensitivity':'Shared-stop sensitivity'}

def style():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'figure.titlesize':13,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'savefig.dpi':240,'legend.frameon':False,'legend.fontsize':8})

def save(fig,folder,stem,caption):
 fig.text(.015,.012,caption,fontsize=8,color='#555555',va='bottom')
 for suffix in ['.png','.pdf']:
  with tempfile.NamedTemporaryFile(suffix=suffix,delete=False) as f:p=Path(f.name)
  try:fig.savefig(p,bbox_inches='tight',facecolor='white');os.replace(p,folder/(stem+suffix))
  finally:p.unlink(missing_ok=True)
 plt.show()

def band(a,z,metric,color,label,ls='-',shade=True,marker='o'):
 s=z.groupby('ratio')[metric].agg(median='median',lo=lambda x:x.quantile(.1),hi=lambda x:x.quantile(.9)).reindex(RATIOS)
 a.plot(range(6),s['median'],color=color,label=label,ls=ls,marker=marker,ms=3.5,lw=1.6)
 if shade:a.fill_between(range(6),s.lo,s.hi,color=color,alpha=.12,lw=0)
 a.set(xticks=range(6),xticklabels=['0','0.5','1','2','4','8'],xlabel='Interaction strength r');a.grid(axis='y',alpha=.12)

def tradeoff(g,c,fit,latent,folder):
 style();fig,axes=plt.subplots(2,2,figsize=(11.6,8));fig.subplots_adjust(left=.085,right=.985,bottom=.15,top=.85,hspace=.42,wspace=.27)
 z=g[g.relative_resolution.eq(.001)]
 for i,protein in enumerate(['GB1','TrpB']):
  e=z[(z.protein==protein)&(z.layer=='empirical')].iloc[0]
  w=z[(z.protein==protein)&(z.observation=='none')&(z.layer=='intrinsic')]
  for j,policy in enumerate(['uniform','greedy']):
   a=axes[i,j];xm=policy+'_steps';ym=policy+'_hit_top';med=[]
   for r,col in zip(RATIOS,RCOL):
    fresh=w[(w.ratio==r)&(w.cohort=='fresh')];old=w[(w.ratio==r)&(w.cohort=='reused')]
    a.scatter(fresh[xm],fresh[ym],s=15,color=col,alpha=.32,zorder=2)
    x=fresh[xm].median();y=fresh[ym].median();med.append((x,y))
    a.plot(x,y,'o',color=col,ms=7,zorder=4)
    a.plot(old[xm].median(),old[ym].median(),'s',mfc='white',mec=col,ms=5,zorder=3)
   a.plot(np.array(med)[:,0],np.array(med)[:,1],color='#999999',lw=.8,zorder=1)
   a.scatter([e[xm]],[e[ym]],marker='*',s=130,c='#202020',zorder=5)
   a.set(title=f'{chr(97+2*i+j)}  {protein} · '+('uniform uphill' if policy=='uniform' else 'greedy'),xlabel='Expected accepted substitutions',ylabel='Probability of reaching upper 1%',ylim=(-.015,1.035));a.grid(alpha=.12)
 handles=[Line2D([],[],marker='o',ls='',color=col,label=f'r={r:g}') for r,col in zip(RATIOS,RCOL)]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.943),ncol=6)
 fig.suptitle('Joint search behaviour across interaction strengths',y=.993)
 save(fig,folder,'01_joint_search_tradeoffs','No added noise; epsilon = 0.001 × empirical calibration range. Filled circles: medians of eight fresh seeds; faint dots: each fresh seed.\nOpen squares: medians of 16 reused seeds. Stars: empirical maps. Connecting lines follow increasing r; they do not interpolate measured settings.')

def epistasis(g,c,fit,latent,folder):
 style();fig,axes=plt.subplots(2,5,figsize=(16,7.2));fig.subplots_adjust(left=.06,right=.985,bottom=.17,top=.86,hspace=.47,wspace=.36)
 z=c[c.relative_resolution.eq(.001)]
 specs=[('epistatic_fraction','Mean cube epistatic fraction'),('ho_fraction','Mean cube order ≥3 fraction'),('sign_epistasis_fraction','Sign-epistatic face fraction'),('slope_median','Median fixed-mutation DR slope'),('sd_fitness','Mean local fitness SD')]
 for i,protein in enumerate(['GB1','TrpB']):
  e=z[(z.protein==protein)&(z.layer=='empirical')].iloc[0];w=z[(z.protein==protein)&(z.observation=='none')&(z.layer=='intrinsic')]
  for j,(metric,label) in enumerate(specs):
   a=axes[i,j];a.axhline(e[metric],color='#222222',lw=1.5,label='Empirical')
   band(a,w[w.cohort=='reused'],metric,'#999999','16 reused seeds',ls='--',shade=False,marker='s')
   band(a,w[w.cohort=='fresh'],metric,'#0072B2','8 fresh seeds')
   a.set(ylabel=label,title=f'{chr(97+5*i+j)}  {protein}')
   if j<3:a.set_ylim(-.02,1.02)
 handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.52,.95),ncol=3)
 fig.suptitle('Global statistics after calibration of the fitness readout',y=.995)
 save(fig,folder,'02_global_statistics','No added noise; confirmation cubes and primary resolution. WHT, sign epistasis and DR were not fit. Shading: 10–90% across fresh seeds, not confidence intervals.\nWHT fractions weight eligible cubes equally; pooled spectra and denominator availability are retained separately. DR is the median across eligible fixed-direction mutation slopes.')

def observation(g,c,fit,latent,folder):
 style();fig,axes=plt.subplots(2,3,figsize=(12.5,7.8));fig.subplots_adjust(left=.08,right=.985,bottom=.16,top=.83,hspace=.42,wspace=.3)
 specs=[(g,'uniform_steps','Expected uniform-uphill substitutions'),(g,'uniform_hit_top','Probability of reaching upper 1%'),(g,'greedy_hit_top','Greedy probability of reaching upper 1%'),(c,'epistatic_fraction','Mean cube epistatic fraction'),(c,'sign_epistasis_fraction','Sign-epistatic face fraction'),(c,'slope_median','Median fixed-mutation DR slope')]
 for i,(data,metric,label) in enumerate(specs):
  a=axes.ravel()[i];z=data[(data.protein=='TrpB')&data.relative_resolution.eq(.001)];e=z[z.layer=='empirical'].iloc[0];a.axhline(e[metric],color='#222222',lw=1.5)
  for mode in MCOL:
   w=z[(z.observation==mode)&(z.cohort=='fresh')]
   band(a,w[w.layer=='intrinsic'],metric,MCOL[mode],MLAB[mode],ls='-' if mode=='none' else '--',shade=mode=='none',marker='o' if mode=='none' else None)
   if mode!='none':band(a,w[w.layer=='observed'],metric,MCOL[mode],MLAB[mode],shade=True)
  a.set(ylabel=label,title=f'{chr(97+i)}  TrpB')
  if metric not in ['uniform_steps','slope_median']:a.set_ylim(-.02,1.02)
 handles=[Line2D([],[],color='#222222',label='Empirical')]+[Line2D([],[],color=col,label=MLAB[mode]) for mode,col in MCOL.items()]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.936),ncol=4)
 fig.suptitle('Observation assumptions reshape the joint comparison',y=.99)
 save(fig,folder,'03_trpb_observation','Fresh seeds: medians and 10–90% ranges. For each observation-aware fit, dashed curves show the intrinsic map and solid curves show its frozen observed map.\nNoise pools are fixed from RUN-172; only the global readout is recalibrated. A noise label describes an explicit sensitivity, not an identified measurement mechanism.')

def sensitivity(g,c,fit,latent,folder):
 style();fig,axes=plt.subplots(2,3,figsize=(12.5,7.8));fig.subplots_adjust(left=.075,right=.985,bottom=.17,top=.86,hspace=.43,wspace=.3)
 colors=['#D55E00','#0072B2','#7B3294'];specs=[('uniform_steps','Expected uniform-uphill substitutions'),('uniform_hit_top','Probability of reaching upper 1%'),('above_floor_largest_fraction','Largest component / above-floor nodes')]
 for i,protein in enumerate(['GB1','TrpB']):
  for j,(metric,label) in enumerate(specs):
   a=axes[i,j]
   for rel,col in zip([0,.001,.01],colors):
    z=g[g.protein.eq(protein)&g.relative_resolution.eq(rel)]
    w=z[z.observation.eq('none')&z.layer.eq('intrinsic')&z.cohort.eq('fresh')]
    band(a,w,metric,col,f'epsilon/range = {rel:g}')
    a.axhline(z[z.layer.eq('empirical')][metric].iloc[0],color=col,lw=1,ls=':')
   if j<2:
    w=latent[latent.protein.eq(protein)&latent.cohort.eq('fresh')]
    band(a,w,metric,'#222222','Latent coordinate, exact',ls='--',shade=False,marker=None)
   a.set(title=f'{chr(97+3*i+j)}  {protein}',ylabel=label)
   if j>0:a.set_ylim(-.02,1.02)
 handles=[Line2D([],[],color=col,label=f'epsilon/range = {rel:g}') for rel,col in zip([0,.001,.01],colors)]+[Line2D([],[],color='#222222',ls='--',label='Latent coordinate, exact')]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.95),ncol=4)
 fig.suptitle('Resolution and neutral connectivity bound the interpretation',y=.997)
 save(fig,folder,'04_resolution_and_neutrality','No added noise; fresh-seed medians and 10–90% ranges. Dotted horizontal lines: empirical values at the matching resolution.\nLatent exact-order curves isolate topology from finite-precision sigmoid ties. Above-floor components exclude the assay floor; connected small changes need not preserve one phenotype.')

def calibration(g,c,fit,latent,folder):
 style();fig,axes=plt.subplots(2,3,figsize=(12.5,7.4));fig.subplots_adjust(left=.075,right=.985,bottom=.17,top=.85,hspace=.43,wspace=.3)
 for i,protein in enumerate(['GB1','TrpB']):
  a=axes[i,0];z=fit[fit.protein.eq(protein)]
  for mode in MCOL:
   w=z[z.observation.eq(mode)].set_index('ratio').reindex(RATIOS)
   if w.objective_total.notna().any():a.plot(range(6),w.objective_total,color=MCOL[mode],marker='o',ms=4,label=MLAB[mode])
  a.set(xticks=range(6),xticklabels=['0','0.5','1','2','4','8'],xlabel='Interaction strength r',ylabel='Calibration objective',title=f'{chr(97+3*i)}  {protein}');a.grid(axis='y',alpha=.12)
  for j,(metric,label) in enumerate([('power_availability','Fraction of cubes with resolvable WHT power'),('slope_availability','Fraction of eligible fixed-mutation slopes')],start=1):
   a=axes[i,j];z=c[c.protein.eq(protein)&c.relative_resolution.eq(.001)];a.axhline(z[z.layer.eq('empirical')][metric].iloc[0],color='#222222',lw=1.5)
   for mode in MCOL:
    w=z[z.cohort.eq('fresh')&z.observation.eq(mode)&z.layer.eq('intrinsic' if mode=='none' else 'observed')]
    if len(w):band(a,w,metric,MCOL[mode],MLAB[mode])
   a.set(title=f'{chr(97+3*i+j)}  {protein}',ylabel=label,ylim=(-.02,1.02))
 handles=[Line2D([],[],color='#222222',label='Empirical availability')]+[Line2D([],[],color=col,label=MLAB[mode]) for mode,col in MCOL.items()]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.946),ncol=4)
 fig.suptitle('Calibration fit and the availability of statistical comparisons',y=.995)
 save(fig,folder,'05_fit_and_availability','Calibration uses only quantiles, mean local SD and floor mass. r=0 and r=4 retain RUN-172 parameters; other ratios use the same finite search budget plus endpoint candidates.\nAvailability uses the primary resolution. Bands: 10–90% across eight fresh seeds. Missing or near-flat measurements do not become zero-valued epistasis or DR estimates.')
