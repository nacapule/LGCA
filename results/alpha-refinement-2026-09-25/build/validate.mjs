import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {gunzipSync} from 'node:zlib';
import {resolve,dirname,join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {measurements} from '../../alpha-2026-09-24/build/measure.mjs';
const out=resolve(dirname(fileURLToPath(import.meta.url)),'..'),root=resolve(out,'../..');
const src=readFileSync(join(root,'results/alpha-2026-09-24/build/engine-source.js'),'utf8');
const {LGCA}=new Function(src+';return {LGCA};')();const hash=createHash('sha256').update(src).digest('hex');
const old=JSON.parse(readFileSync(join(root,'results/alpha-2026-09-24/validation.json')));assert.equal(hash,old.engineHash,'Engine must match previously parity-verified engine');
const manifest=JSON.parse(readFileSync(join(out,'data/manifest.json')));assert.equal(manifest.jobs.length,60);assert.equal(hash,manifest.engineHash);
let states=0,pilotMatches=0;const initial=new Map();
const keys=['N','nmax','kmax','siteFraction','channelFraction','clusterMass','clusterArea','clusterDensity','clusterFraction','polar','localPolar','effectiveArea'];
const runs=manifest.jobs.map(job=>JSON.parse(readFileSync(join(out,'data',job.id,'run.json'))));
for(const r of runs){
 const j=r.job;assert.equal(r.engineHash,hash);assert.equal(r.snaps.at(-1).t,manifest.jobs.find(x=>x.id===j.id).end);assert.equal(r.series.at(-1).t,r.snaps.at(-1).t);
 const first=r.snaps[0],tag=`${j.size}-${j.seed}`;
 if(initial.has(tag)){assert.equal(first.sha256,initial.get(tag).sha256);assert.deepEqual(first.rng,initial.get(tag).rng);}else initial.set(tag,first);
 for(const m of r.series){assert.equal(m.N,first.N);for(const k of keys)assert(Number.isFinite(m[k]),`${j.id} invalid ${k}`);}
 for(const s of r.snaps){
  const buf=gunzipSync(readFileSync(join(out,'data',j.id,s.file)));assert.equal(createHash('sha256').update(buf).digest('hex'),s.sha256);
  const occ=new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4);assert.equal(occ.length,j.size*j.size*6);assert(occ.every(v=>v>=0));assert.equal(s.rng.length,33);assert(s.rng[32]>=0&&s.rng[32]<32);
  const m=measurements(occ,j.size,j.size);for(const k of keys)assert.equal(m[k],s[k]);assert.equal(s.N,first.N);
  for(const z of s.thresholdChecks){const a=measurements(occ,j.size,j.size,z.factor);assert.equal(a.threshold,z.threshold);assert.equal(a.clusterMass,z.mass);assert.equal(a.clusterArea,z.area);assert.equal(a.clusterDensity,z.density);}
  states++;
 }
 if(j.size===90&&[.75,.9].includes(j.alpha)&&[12345,777,424242].includes(j.seed)){
  const prev=JSON.parse(readFileSync(join(root,`results/alpha-2026-09-24/data/neigh-fixed-sensitivity-L90-a${j.alpha}-seed${j.seed}/run.json`))).snaps.find(s=>s.t===10000);
  const now=r.snaps.find(s=>s.t===10000);assert.equal(now.sha256,prev.sha256);assert.deepEqual(now.rng,prev.rng);pilotMatches++;
 }
}
const r=runs.find(r=>r.job.size===120&&r.job.alpha===.85&&r.job.seed===12345),j=r.job,s=r.snaps.find(s=>s.t===10000),buf=gunzipSync(readFileSync(join(out,'data',j.id,s.file)));
const sim=new LGCA({model:'boson',W:j.size,H:j.size,dens:j.dens,seed:j.seed,kernel:'power',alpha:j.alpha,bosonField:'neigh',bosonAlign:'polar'});
sim.occ.set(new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4));sim.rng.setState(new Uint32Array(s.rng));sim.t=s.t;
for(let t=s.t;t<j.end;t++)sim.step(j.sens);
assert.equal(createHash('sha256').update(Buffer.from(sim.occ.buffer)).digest('hex'),r.snaps.at(-1).sha256);assert.deepEqual(Array.from(sim.rng.getState()),r.snaps.at(-1).rng);
const result={runs:runs.length,states,pilotMatches,engineHash:hash,unchangedFromParityVerifiedEngine:true,nonnegativeAndConserved:true,recomputedClusterStatistics:true,pairedInitialization:true,resumeCheck:{id:j.id,from:s.t,to:j.end,stateAndRngExact:true}};
writeFileSync(join(out,'validation.json'),JSON.stringify(result,null,2));console.log('PASS '+JSON.stringify(result));
