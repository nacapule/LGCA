import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter
out=Path(__file__).resolve().parent.parent
analysis=json.loads((out/'analysis.json').read_text());runs=json.loads((out/'runs.json').read_text());manifest=analysis['manifest']
figdir=out/'figures';figdir.mkdir(exist_ok=True)
alphas=manifest['alphas'];densities=manifest['densities'];colors=['#18698e','#ca6d28','#7357a1','#3e8156']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#19252d','text.color':'#19252d','xtick.color':'#19252d','ytick.color':'#19252d','figure.facecolor':'#fcfcfa','axes.facecolor':'#fcfcfa','savefig.facecolor':'#fcfcfa'})
fmt=FuncFormatter(lambda v,_:'0' if v==0 else f'{v/1000:g}k')
for a in alphas:
 series=[]
 for d,c in zip(densities,colors):
  rs=[r for r in runs if r['job']['alpha']==a and r['job']['dens']==d]
  ts=np.array([m['t'] for m in rs[0]['series']]);v=np.array([[m['clusterDensity'] for m in r['series']] for r in rs]);N=np.mean([r['series'][0]['N'] for r in rs]);mean=v.mean(axis=0);sd=v.std(axis=0,ddof=1)
  series.append((ts,mean,sd,c,f'Mean density {N/manifest["size"]**2:.2f} · N ≈ {N:,.0f}'))
 for t in [10000,20000,40000,80000]:
  for focus in [False,True]:
   if focus:fig,axes=plt.subplots(2,1,figsize=(8.3,5),layout='constrained',sharey=True,gridspec_kw={'height_ratios':[2.4,1]});stops=[5000,t]
   else:fig,ax=plt.subplots(figsize=(8.3,4.4),layout='constrained');axes=[ax];stops=[t]
   for ts,mean,sd,c,label in series:
    for ax,stop in zip(axes,stops):
     use=ts<=stop;ax.plot(ts[use],mean[use],color=c,label=label,lw=1.5);ax.fill_between(ts[use],np.maximum(0,mean[use]-sd[use]),mean[use]+sd[use],color=c,alpha=.10)
   for ax,stop in zip(axes,stops):
    ax.set(xlim=(0,stop),ylim=(0,None));ax.grid(alpha=.18);ax.xaxis.set_major_locator(MaxNLocator(5));ax.xaxis.set_major_formatter(fmt)
   axes[0].set_title(f'α = {a:.2f}  ·  '+('first 5,000 steps  ·  ' if focus else '')+'mean ± SD, 5 seeds',loc='left',fontsize=10,color='#5a676c',pad=12)
   axes[0].set_ylabel('Cluster density (particles/site)');axes[0].legend(frameon=False,fontsize=8.5,loc='upper right',ncol=2)
   axes[-1].set_xlabel('Simulation step')
   if focus:
    axes[1].set_title(f'Full run · {t:,} steps',loc='left',fontsize=10,color='#5a676c',pad=10);axes[1].axvspan(0,5000,color='#5a676c',alpha=.06,zorder=0)
   for ext in ['png','pdf']:fig.savefig(figdir/f'density-time-{"focus-" if focus else ""}a{a:g}-t{t}.{ext}',dpi=170)
   plt.close(fig)
for t in [10000,20000,40000,80000]:
 fig,ax=plt.subplots(figsize=(8.3,4.4),layout='constrained')
 for a,c in zip(alphas,colors):
  gs=[next(g for g in analysis['groups'] if g['alpha']==a and g['dens']==d) for d in densities]
  stats=[next(m for m in g['checkpoints'] if m['t']==t)['lateWindow']['clusterDensity'] for g in gs]
  ax.errorbar([g['N']['mean']/manifest['size']**2 for g in gs],[s['mean'] for s in stats],yerr=[s['sd'] for s in stats],label=f'α = {a:.2f}',color=c,marker='o',capsize=4,lw=1.8)
 ax.set(xlabel='Mean system density (particles/site)',ylabel='Cluster density (particles/site)',xticks=[6*d for d in densities],ylim=(0,None));ax.grid(alpha=.18);ax.legend(frameon=False)
 ax.set_title(f'Steps {int(t*.75)+100:,}–{t:,} · mean ± SD · 5 seeds',loc='left',fontsize=10,color='#5a676c',pad=14)
 for ext in ['png','pdf']:fig.savefig(figdir/f'density-size-t{t}.{ext}',dpi=170)
 plt.close(fig)
print('Density sweep plots rendered.')
