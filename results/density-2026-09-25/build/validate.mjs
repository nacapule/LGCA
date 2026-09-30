import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {gunzipSync} from 'node:zlib';
import {resolve,dirname,join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {measurements} from '../../alpha-2026-09-24/build/measure.mjs';
const out=resolve(dirname(fileURLToPath(import.meta.url)),'..');
const source=readFileSync(join(out,'build/engine-source.js'),'utf8'),hash=createHash('sha256').update(source).digest('hex');
const {LGCA}=new Function(source+';return {LGCA};')();
const manifest=JSON.parse(readFileSync(join(out,'data/manifest.json')));assert.equal(hash,manifest.engineHash);assert.equal(manifest.jobs.length,60);
const provenance=JSON.parse(readFileSync(join(out,'provenance.json')));assert.equal(hash,provenance.engineHash);
const keys=['N','mean','nmax','kmax','siteFraction','channelFraction','effectiveArea','occupied','polar','localPolar','weightedLocalPolar','threshold','clusters','clusterMass','clusterArea','clusterDensity','clusterFraction','clusterPolar'];
const runs=manifest.jobs.map(j=>JSON.parse(readFileSync(join(out,'data',j.id,'run.json'))));
const init=new Map();let states=0;
const expectedTimes=[...Array.from({length:201},(_,i)=>25*i),...Array.from({length:750},(_,i)=>5100+100*i)];
for(const r of runs){
 assert.equal(r.engineHash,hash);assert.deepEqual(r.job,manifest.jobs.find(j=>j.id===r.job.id));assert.equal(r.series.at(-1).t,80000);assert.equal(r.snaps.at(-1).t,80000);
 const tag=`${r.job.dens}-${r.job.seed}`,first=r.snaps[0];
 if(init.has(tag)){assert.equal(first.sha256,init.get(tag).sha256);assert.deepEqual(first.rng,init.get(tag).rng);}else{
  init.set(tag,first);const j=r.job;const fresh=new LGCA({model:'boson',W:j.size,H:j.size,dens:j.dens,seed:j.seed,kernel:'power',alpha:j.alpha,bosonField:j.field,bosonAlign:'polar'});
  assert.equal(createHash('sha256').update(Buffer.from(fresh.occ.buffer)).digest('hex'),first.sha256);assert.deepEqual(Array.from(fresh.rng.getState()),first.rng);assert.equal(fresh.nuparts,first.N);
 }
 let previous=-1;
 for(const m of r.series){assert(m.t>previous);previous=m.t;assert.equal(m.N,first.N);for(const k of keys)assert(Number.isFinite(m[k]),`${r.job.id}: ${k}`);}
 assert.deepEqual(r.series.map(m=>m.t),expectedTimes);
 for(const m of r.series){assert(m.clusterArea>=0&&m.clusterArea<=r.job.size**2);assert(m.clusterMass>=0&&m.clusterMass<=m.N);assert.equal(m.clusterDensity,m.clusterArea?m.clusterMass/m.clusterArea:0);assert.equal(m.clusterFraction,m.clusterMass/m.N);}
 assert.equal(r.series.length,951);assert.deepEqual(r.snaps.map(s=>s.t),manifest.checkpoints);
 for(const s of r.snaps){
  const b=gunzipSync(readFileSync(join(out,'data',r.job.id,s.file)));assert.equal(createHash('sha256').update(b).digest('hex'),s.sha256);
  const occ=new Int32Array(b.buffer,b.byteOffset,b.byteLength/4);assert.equal(occ.length,r.job.size**2*6);assert(occ.every(v=>v>=0));
  const m=measurements(occ,r.job.size,r.job.size);for(const k of keys)assert.equal(s[k],m[k],`${r.job.id} t${s.t}: ${k}`);assert.equal(m.N,first.N);
  assert.equal(s.rng.length,33);assert(s.rng[32]>=0&&s.rng[32]<32);
  for(const z of s.thresholdChecks){const q=measurements(occ,r.job.size,r.job.size,z.factor);assert.equal(z.threshold,q.threshold);assert.equal(z.mass,q.clusterMass);assert.equal(z.area,q.clusterArea);assert.equal(z.density,q.clusterDensity);}
  states++;
 }
}
const resumeChecks=[];
for(const [alpha,dens,seed,from,to] of [[.8,.8,12345,40000,80000],[1,.1,777,1000,2000]]){
 const r=runs.find(r=>r.job.alpha===alpha&&r.job.dens===dens&&r.job.seed===seed),j=r.job,s=r.snaps.find(s=>s.t===from),target=r.snaps.find(s=>s.t===to);
 const b=gunzipSync(readFileSync(join(out,'data',j.id,s.file)));
 const sim=new LGCA({model:'boson',W:j.size,H:j.size,dens:j.dens,seed:j.seed,kernel:'power',alpha:j.alpha,bosonField:j.field,bosonAlign:'polar'});
 sim.occ.set(new Int32Array(b.buffer,b.byteOffset,b.byteLength/4));sim.rng.setState(new Uint32Array(s.rng));sim.t=from;
 for(let t=from;t<to;t++)sim.step(j.sens);
 assert.equal(createHash('sha256').update(Buffer.from(sim.occ.buffer)).digest('hex'),target.sha256);assert.deepEqual(Array.from(sim.rng.getState()),target.rng);
 resumeChecks.push({id:j.id,from,to,stateAndRngExact:true});
}
const analysisSha256=createHash('sha256').update(readFileSync(join(out,'analysis.json'))).digest('hex');
const result={runs:runs.length,states,engineHash:hash,analysisSha256,horizon:80000,nonnegativeAndConserved:true,recomputedClusterStatistics:true,pairedInitialization:true,resumeChecks};
writeFileSync(join(out,'validation.json'),JSON.stringify(result,null,2));console.log('PASS '+JSON.stringify(result));
