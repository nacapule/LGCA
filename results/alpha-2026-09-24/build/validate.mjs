import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {gunzipSync} from 'node:zlib';
import {resolve,dirname,join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {measurements} from './measure.mjs';
const root=resolve(dirname(fileURLToPath(import.meta.url)),'../../..'),out=join(root,'results/alpha-2026-09-24'),data=join(out,'data');
const src=readFileSync(join(root,'results/alpha-2026-09-24/build/engine-source.js'),'utf8');
const {LGCA}=new Function(src+';return {LGCA};')();const hash=createHash('sha256').update(src).digest('hex');
// Hand-constructed geometry tests: periodic wrap, adjacency, mass and thresholds.
const lattice=new Int32Array(8*6*6);lattice[0]=10;lattice[7*6]=20;lattice[(3*8+3)*6]=7;
const m=measurements(lattice,8,6,2,true);assert.equal(m.N,37);assert.equal(m.nmax,20);assert.equal(m.clusterArea,2);assert.equal(m.clusterMass,30);assert.equal(m.clusters,2);assert.equal(m.clusterDensity,15);assert.equal(m.localPolar,1);
assert.equal(m.edgeProfile.find(p=>p.distance===0).density,15);
const jobs=JSON.parse(readFileSync(join(data,'manifest.json'))).jobs;let states=0;const cycles=[];
for(const job of jobs){
 const r=JSON.parse(readFileSync(join(data,job.id,'run.json')));assert.equal(r.engineHash,hash);assert.equal(r.snaps.at(-1).t,job.end);
 for(const s of r.snaps){const buf=gunzipSync(readFileSync(join(data,job.id,s.file)));assert.equal(createHash('sha256').update(buf).digest('hex'),s.sha256);const occ=new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4);assert.equal(occ.length,job.size*job.size*6);assert(occ.every(x=>x>=0));const z=measurements(occ,job.size,job.size,2,true);for(const k of ['N','nmax','kmax','clusterMass','clusterArea','clusterDensity','polar','localPolar','effectiveArea'])assert.equal(z[k],s[k]);assert.equal(z.N,r.snaps[0].N);states++;}
 if(job.field==='neigh'&&job.control==='fixed-sensitivity'){
  const s=r.snaps.at(-1),buf=gunzipSync(readFileSync(join(data,job.id,s.file))),occ=new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4);
  const sim=new LGCA({model:'boson',W:job.size,H:job.size,dens:.4,seed:job.seed,kernel:'power',alpha:job.alpha,bosonField:job.field,bosonAlign:'polar'});sim.occ.set(occ);sim.rng.setState(new Uint32Array(s.rng));sim.t=s.t;
  for(let t=0;t<job.size;t++)sim.step(job.sens);
  cycles.push({id:job.id,alpha:job.alpha,size:job.size,seed:job.seed,t:s.t,lag:job.size,exactStateRecurrence:sim.occ.every((v,i)=>v===occ[i])});
 }
}
// Resume check against an independent uninterrupted run from the seed.
const sample=jobs.find(j=>j.size===120&&j.alpha===.5&&j.seed===12345);const sim=new LGCA({model:'boson',W:120,H:120,dens:.4,seed:sample.seed,kernel:'power',alpha:.5,bosonField:'neigh',bosonAlign:'polar'});
for(let t=0;t<sample.end;t++)sim.step(sample.sens);
const saved=JSON.parse(readFileSync(join(data,sample.id,'run.json'))).snaps.at(-1);
assert.equal(createHash('sha256').update(Buffer.from(sim.occ.buffer)).digest('hex'),saved.sha256);assert.deepEqual(Array.from(sim.rng.getState()),saved.rng);
const result={states,completedRuns:jobs.length,engineHash:hash,geometryChecks:'passed',resumeCheck:'120x120 alpha 0.5 seed12345: state + RNG identical after 10000 steps',cycles};writeFileSync(join(out,'validation.json'),JSON.stringify(result,null,2));
console.log(`PASS ${jobs.length} runs / ${states} states: hashes, nonnegative counts, conservation, recomputed metrics; periodic-cluster fixtures; resumed/uninterrupted state and RNG.`);
console.log('Exact recurrence after one box-crossing time at t=10000: '+cycles.filter(c=>c.exactStateRecurrence).length+'/'+cycles.length+' neighbour runs.');
