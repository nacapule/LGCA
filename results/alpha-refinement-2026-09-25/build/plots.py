import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter
out=Path(__file__).resolve().parent.parent
analysis=json.loads((out/'analysis.json').read_text());runs=json.loads((out/'runs.json').read_text())
figdir=out/'figures';figdir.mkdir(exist_ok=True)
alphas=[.75,.8,.85,.9];sizes=[90,120,180];colors=['#18698e','#ca6d28','#7357a1','#3e8156']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#19252d','text.color':'#19252d','xtick.color':'#19252d','ytick.color':'#19252d','figure.facecolor':'#fcfcfa','axes.facecolor':'#fcfcfa','savefig.facecolor':'#fcfcfa'})
end=min(r['snaps'][-1]['t'] for r in runs)
for t in [10000,20000,40000,80000]:
 if t>end:continue
 fig,ax=plt.subplots(figsize=(8.3,4.4),layout='constrained')
 for a,c in zip(alphas,colors):
  groups=[next(g for g in analysis['groups'] if g['alpha']==a and g['size']==L) for L in sizes]
  stats=[next(m for m in g['checkpoints'] if m['t']==t)['lateWindow']['clusterDensity'] for g in groups]
  ax.errorbar(sizes,[s['mean'] for s in stats],yerr=[s['sd'] for s in stats],label=f'α = {a:.2f}',color=c,marker='o',capsize=4,lw=1.8)
 ax.set(xlabel='System size L × L',ylabel='Cluster density (particles/site)',xticks=sizes,xticklabels=[f'{L} × {L}' for L in sizes])
 ax.set_ylim(bottom=0);ax.grid(alpha=.18);ax.legend(ncol=2,frameon=False)
 ax.set_title(f'Steps {int(t*.75)+100:,}–{t:,}  ·  mean ± SD  ·  5 seeds',loc='left',fontsize=10,color='#5a676c',pad=14)
 for ext in ['png','pdf']:fig.savefig(figdir/f'density-size-t{t}.{ext}',dpi=170)
 plt.close(fig)
 fig,axes=plt.subplots(2,2,figsize=(9,5.9),layout='constrained')
 for ax,a in zip(axes.flat,alphas):
  for L,c in zip(sizes,colors):
   rs=[r for r in runs if r['job']['alpha']==a and r['job']['size']==L]
   ts=np.array([m['t'] for m in rs[0]['series']]);vals=np.array([[m['clusterDensity'] for m in r['series']] for r in rs]);mean=vals.mean(axis=0);sd=vals.std(axis=0,ddof=1)
   use=ts<=t;ax.plot(ts[use],mean[use],color=c,label=f'L = {L}',lw=1.5);ax.fill_between(ts[use],np.maximum(0,mean[use]-sd[use]),mean[use]+sd[use],color=c,alpha=.10)
  ax.set_title(f'α = {a:.2f}',loc='left',fontweight='bold',fontsize=11);ax.set_ylim(bottom=0);ax.set_xlim(0,t);ax.xaxis.set_major_locator(MaxNLocator(4));ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v/1000:g}k'));ax.grid(alpha=.18)
 axes[0,0].legend(frameon=False,ncol=3,fontsize=8,loc='upper right')
 for ax in axes[-1,:]:ax.set_xlabel('Simulation step')
 for ax in axes[:,0]:ax.set_ylabel('Cluster density (particles/site)')
 fig.supxlabel('Mean ± SD · 5 seeds · independent y-scales',fontsize=9,color='#5a676c')
 for ext in ['png','pdf']:
  fig.savefig(figdir/f'density-time-t{t}.{ext}',dpi=170)
  if t==end:fig.savefig(figdir/f'density-time.{ext}',dpi=170)
 plt.close(fig)
# One alpha per large panel: Manik's density-time/particle-total comparison.
for a in alphas:
 for t in [10000,20000,40000,80000]:
  if t>end:continue
  fig,ax=plt.subplots(figsize=(8.3,4.4),layout='constrained')
  for L,c in zip(sizes,colors):
   rs=[r for r in runs if r['job']['alpha']==a and r['job']['size']==L]
   ts=np.array([m['t'] for m in rs[0]['series']])
   vals=np.array([[m['clusterDensity'] for m in r['series']] for r in rs])
   mean=vals.mean(axis=0);sd=vals.std(axis=0,ddof=1);use=ts<=t
   total=np.mean([r['series'][0]['N'] for r in rs])
   ax.plot(ts[use],mean[use],color=c,label=f'N ≈ {total:,.0f}  ·  {L}×{L}',lw=1.7)
   ax.fill_between(ts[use],np.maximum(0,mean[use]-sd[use]),mean[use]+sd[use],color=c,alpha=.10)
  ax.set(xlabel='Simulation step',ylabel='Cluster density (particles/site)',xlim=(0,t),ylim=(0,None))
  ax.xaxis.set_major_locator(MaxNLocator(5));ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:'0' if v==0 else f'{v/1000:g}k'))
  ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=9,loc='upper right')
  ax.set_title(f'α = {a:.2f}  ·  mean ± SD  ·  5 seeds',loc='left',fontsize=10,color='#5a676c',pad=14)
  for ext in ['png','pdf']:fig.savefig(figdir/f'density-time-a{a:g}-t{t}.{ext}',dpi=170)
  plt.close(fig)
  fig,axes=plt.subplots(2,1,figsize=(8.3,5),layout='constrained',sharey=True,gridspec_kw={'height_ratios':[2.4,1]})
  for L,c in zip(sizes,colors):
   rs=[r for r in runs if r['job']['alpha']==a and r['job']['size']==L]
   ts=np.array([m['t'] for m in rs[0]['series']]);vals=np.array([[m['clusterDensity'] for m in r['series']] for r in rs])
   mean=vals.mean(axis=0);sd=vals.std(axis=0,ddof=1);total=np.mean([r['series'][0]['N'] for r in rs])
   for ax,stop in zip(axes,[5000,t]):
    use=ts<=stop;ax.plot(ts[use],mean[use],color=c,label=f'N ≈ {total:,.0f}  ·  {L}×{L}',lw=1.5)
    ax.fill_between(ts[use],np.maximum(0,mean[use]-sd[use]),mean[use]+sd[use],color=c,alpha=.10)
  for ax,stop in zip(axes,[5000,t]):
   ax.set(xlim=(0,stop),ylim=(0,None));ax.grid(alpha=.18);ax.xaxis.set_major_locator(MaxNLocator(5));ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:'0' if v==0 else f'{v/1000:g}k'))
  axes[0].set_title(f'α = {a:.2f}  ·  first 5,000 steps  ·  mean ± SD, 5 seeds',loc='left',fontsize=10,color='#5a676c',pad=12)
  axes[0].set_ylabel('Cluster density (particles/site)');axes[0].legend(frameon=False,fontsize=9,loc='upper right')
  axes[1].set_title(f'Full run · {t:,} steps',loc='left',fontsize=10,color='#5a676c',pad=10);axes[1].set_xlabel('Simulation step')
  axes[1].axvspan(0,5000,color='#5a676c',alpha=.06,zorder=0)
  for ext in ['png','pdf']:fig.savefig(figdir/f'density-time-focus-a{a:g}-t{t}.{ext}',dpi=170)
  plt.close(fig)
print('Scientific figures rendered through',end)
