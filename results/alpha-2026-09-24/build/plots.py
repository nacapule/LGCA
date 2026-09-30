from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
runs=[json.loads(p.read_text()) for p in (ROOT/'data').glob('*/run.json')]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
colors=['#146b91','#b15e22','#7255a1']
def get(size,alpha,field='neigh',control='fixed-sensitivity',t=10000):
 return [next(m for m in r['snaps'] if m['t']==t) for r in runs if r['job']['size']==size and r['job']['alpha']==alpha and r['job']['field']==field and r['job']['control']==control]
def values(ms,k):
 v=np.array([m[k] for m in ms]);return v.mean(),v.std(ddof=1)
def finish(fig,name):
 fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(12,3.9))
for a,col in zip([0,.5,1],colors):
 for ax,key,title in zip(axes,['nmax','siteFraction','clusterDensity'],['Peak site count','Peak site / total particles','Largest dense cluster: mass / area']):
  vals=[values(get(L,a),key) for L in [60,120,180]]
  ax.errorbar([60,120,180],[v[0] for v in vals],yerr=[v[1] for v in vals],color=col,marker='o',capsize=4,label=f'α={a}')
  ax.set_title(title,fontsize=11);ax.set_xlabel('Lattice side L (sites)');ax.set_xticks([60,120,180]);ax.grid(alpha=.2)
axes[0].set_ylabel('Particles');axes[1].set_ylabel('Fraction');axes[2].set_ylabel('Particles / selected site');axes[0].legend()
fig.suptitle('Neighbours, polar bosons · sensitivity 6 · dens 0.4 · t=10,000',fontsize=13)
fig.text(.5,-.03,'Points: mean of 3 seeds; bars: sample SD. Fixed observation time does not establish the large-system limit.',ha='center',fontsize=10)
fig.tight_layout();finish(fig,'size-scaling')
fig,axes=plt.subplots(1,3,figsize=(12,3.9));alphas=[0,.25,.5,.75,.9,1]
for ax,key,title in zip(axes,['nmax','clusterDensity','polar'],['Peak site count','Largest dense cluster: mass / area','Global polar order']):
 vals=[values(get(90,a),key) for a in alphas];ax.errorbar(alphas,[v[0] for v in vals],yerr=[v[1] for v in vals],color=colors[0],marker='o',capsize=4)
 ax.set_title(title,fontsize=11);ax.set_xlabel('Normalization exponent α');ax.set_xticks(alphas);ax.grid(alpha=.2)
axes[0].set_ylabel('Particles');axes[1].set_ylabel('Particles / selected site');axes[2].set_ylabel('P')
fig.suptitle('90×90 · neighbours · sensitivity 6 · dens 0.4 · t=10,000',fontsize=13)
fig.text(.5,-.03,'Mean ± sample SD, 3 seeds. A finite-alpha scan shows a crossover; it does not locate a critical exponent.',ha='center',fontsize=10)
fig.tight_layout();finish(fig,'alpha-grid')
fig,axes=plt.subplots(1,2,figsize=(12,4.1))
for a,col in zip([0,.5,1],colors):
 ms=get(180,a);distances=sorted(set(b['distance'] for m in ms for b in m.get('edgeProfile',[])))
 for ax,key in zip(axes,['density','localPolar']):
  xx=[];yy=[]
  for d in distances:
   vals=[b[key] for m in ms for b in m.get('edgeProfile',[]) if b['distance']==d and b[key] is not None]
   if len(vals)==3:xx.append(d);yy.append(np.mean(vals) if key!='density' or np.mean(vals)>0 else np.nan)
  ax.plot(xx,yy,'o-',ms=4,color=col,label=f'α={a}')
for ax in axes:
 ax.axvline(0,color='#7b858c',lw=1,ls='--');ax.set_xlabel('Signed graph distance from boundary (lattice hops)');ax.grid(alpha=.2);ax.set_xlim(-8,12)
axes[0].set_title('Site-averaged density');axes[0].set_ylabel('Particles / site');axes[0].set_yscale('log');axes[0].legend()
axes[1].set_title('Site-averaged local polar order (occupied sites)');axes[1].set_ylabel('Local P');axes[1].set_ylim(0,1.05)
fig.suptitle('Largest dense component · 180×180 · neighbours · t=10,000 · threshold ≥5',fontsize=13)
fig.text(.5,-.06,'Negative: interior; 0: boundary sites; positive: exterior. Mean of 3 seeds where all 3 contribute.\nExploratory: exterior bins can include other clusters; singleton local order equals 1.',ha='center',fontsize=10)
fig.tight_layout();finish(fig,'edge-profiles')
# Physical feedback illustration; exact one-step probability for a perfectly aligned own-site stack.
fig,axes=plt.subplots(1,2,figsize=(11,4));n=np.logspace(0,5,300)
for a,col in zip([0,.5,.9,1],['#146b91','#b15e22','#3c7f56','#7255a1']):
 u=6*n**(1-a);rest=2*np.exp(-u/2)+2*np.exp(-1.5*u)+np.exp(-2*u);escape=rest/(1+rest)
 axes[0].loglog(n,n**(1-a),color=col,label=f'α={a}')
 axes[1].loglog(n,np.maximum(1e-14,n*escape),color=col,label=f'α={a}')
axes[0].set_ylabel('Crowd multiplier M^(1−α)');axes[0].set_xlabel('Contributing particle count M');axes[0].legend();axes[0].grid(alpha=.2)
axes[1].set_ylabel('Expected departures per step');axes[1].set_xlabel('Aligned own-site stack size n');axes[1].set_ylim(1e-12,1e5);axes[1].grid(alpha=.2)
fig.suptitle('Why normalizing changes retention · sensitivity 6',fontsize=13)
fig.text(.5,-.045,'Own-site, perfectly aligned configuration only. Departures = n(1−p_forward). This is not a phase diagram.',ha='center')
fig.tight_layout();finish(fig,'feedback')
print('Saved four scientific figure pairs (PNG/PDF).')
