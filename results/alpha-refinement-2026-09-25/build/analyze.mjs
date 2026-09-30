import {readFileSync,writeFileSync} from 'node:fs';
import {resolve,dirname,join} from 'node:path';
import {fileURLToPath} from 'node:url';
const out=resolve(dirname(fileURLToPath(import.meta.url)),'..');
const manifest=JSON.parse(readFileSync(join(out,'data/manifest.json')));
const runs=manifest.jobs.map(j=>JSON.parse(readFileSync(join(out,'data',j.id,'run.json'))));
for(const r of runs)if(r.snaps.at(-1).t<manifest.jobs.find(j=>j.id===r.job.id).end)throw Error('Incomplete: '+r.job.id);
const keys=['nmax','siteFraction','channelFraction','effectiveArea','clusterDensity','clusterMass','clusterArea','clusterFraction','polar'];
const mean=v=>v.reduce((a,b)=>a+b,0)/v.length;
const stat=v=>({mean:mean(v),sd:Math.sqrt(v.reduce((a,b)=>a+(b-mean(v))**2,0)/(v.length-1)),min:Math.min(...v),max:Math.max(...v),values:v});
const groups=[];
for(const size of [90,120,180])for(const alpha of [.75,.8,.85,.9]){
 const rs=runs.filter(r=>r.job.size===size&&r.job.alpha===alpha);if(rs.length!==5)throw Error('Expected five seeds');
 const end=Math.min(...rs.map(r=>r.snaps.at(-1).t)),times=[1000,5000,10000,20000,40000,80000].filter(t=>t<=end);
 const checkpoints=times.map(t=>({t,...Object.fromEntries(keys.map(k=>[k,stat(rs.map(r=>r.snaps.find(m=>m.t===t)[k]))])),lateWindow:{lo:t*.75,hi:t,...Object.fromEntries(keys.map(k=>[k,stat(rs.map(r=>mean(r.series.filter(m=>m.t>t*.75&&m.t<=t).map(m=>m[k]))))]))}}));
 const windows=[[end*.5,end*.75],[end*.75,end]].map(([lo,hi])=>({lo,hi,...Object.fromEntries(keys.map(k=>[k,stat(rs.map(r=>mean(r.series.filter(m=>m.t>lo&&m.t<=hi).map(m=>m[k]))))]))}));
 const drift={};for(const k of keys){const a=windows[0][k].values,b=windows[1][k].values;drift[k]={ensembleRelative:Math.abs(mean(b)-mean(a))/Math.max(Math.abs(mean(a)),1e-12),perSeedRelative:a.map((x,i)=>Math.abs(b[i]-x)/Math.max(Math.abs(x),1e-12))};}
 const thresholdRobustness=[1.5,2,4].map(f=>({factor:f,density:stat(rs.map(r=>{const m=r.snaps.at(-1);return f===2?m.clusterDensity:m.thresholdChecks.find(x=>x.factor===f).density;}))}));
 groups.push({size,alpha,end,N:stat(rs.map(r=>r.snaps[0].N)),checkpoints,windows,drift,thresholdRobustness});
}
const summary={manifest,groups};writeFileSync(join(out,'analysis.json'),JSON.stringify(summary,null,2));
writeFileSync(join(out,'runs.json'),JSON.stringify(runs));
const cols=['size','alpha','seed','dens','sens','t','N',...keys];
writeFileSync(join(out,'measurements.csv'),cols.join(',')+'\n'+runs.flatMap(r=>r.series.map(m=>cols.map(k=>r.job[k]??m[k]).join(','))).join('\n')+'\n');
console.log('L alpha rho_cluster nmax area mass/N drift_rho(last-half quarters)');
for(const g of groups){const m=g.checkpoints.at(-1);console.log(g.size,g.alpha,m.clusterDensity.mean.toFixed(2),m.nmax.mean.toFixed(1),m.clusterArea.mean.toFixed(1),m.clusterFraction.mean.toFixed(3),(100*g.drift.clusterDensity.ensembleRelative).toFixed(1)+'%');}
