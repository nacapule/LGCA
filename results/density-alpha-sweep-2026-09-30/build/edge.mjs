// Edge profiles (Zoom Sep 23 task) from saved exact states: density and local polar order versus hop distance
// to the largest cluster's boundary (same definition as results/alpha-2026-09-24, METHODS "Edge profiles").
// No simulation; reads state-t{T}.bin.gz, checks its hash, averages over seeds at distances all seeds reach.
import {readFileSync,writeFileSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {fileURLToPath} from 'node:url';
import {dirname,resolve,join} from 'node:path';
import {measurements} from '../../alpha-2026-09-24/build/measure.mjs';
const out=resolve(dirname(fileURLToPath(import.meta.url)),'..'),T=+(process.argv[2]||20000);
const alphas=[.5,.8,.85,.9,1],dens=[.1,.4,1.6],seeds=[12345,777,424242,20260925,314159],result=[];
for(const a of alphas)for(const d of dens){
 const per=[];
 for(const s of seeds){
  const dir=join(out,'data',`both-L90-a${a}-d${d}-seed${s}`);let run;
  try{run=JSON.parse(readFileSync(join(dir,'run.json')));}catch{continue;}
  const snap=run.snaps.find(x=>x.t===T);if(!snap)continue;
  const buf=gunzipSync(readFileSync(join(dir,snap.file)));
  if(createHash('sha256').update(buf).digest('hex')!==snap.sha256)throw Error('hash '+dir);
  const m=measurements(new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4),90,90,2,true);
  per.push({seed:s,mean:m.mean,profile:m.edgeProfile||[]});
 }
 if(!per.length)continue;
 const ds=[...new Set(per.flatMap(p=>p.profile.map(b=>b.distance)))].sort((x,y)=>x-y);
 const prof=ds.map(dd=>{const bs=per.map(p=>p.profile.find(b=>b.distance===dd));
  if(bs.some(b=>!b))return null;
  const ord=bs.filter(b=>b.localPolar!==null);
  return {distance:dd,density:bs.reduce((s,b)=>s+b.density,0)/bs.length,localPolar:ord.length===bs.length?ord.reduce((s,b)=>s+b.localPolar,0)/ord.length:null};
 }).filter(Boolean);
 result.push({alpha:a,dens:d,mean:per.reduce((s,p)=>s+p.mean,0)/per.length,seeds:per.length,profile:prof});
}
writeFileSync(join(out,`edge-profiles-t${T}.json`),JSON.stringify(result));
for(const r of result)console.log(r.alpha,r.dens,r.seeds,'seeds',r.profile.length,'bins');
