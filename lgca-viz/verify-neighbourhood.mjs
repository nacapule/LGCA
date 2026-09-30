#!/usr/bin/env node
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {dirname,resolve,join} from 'node:path';
const root=resolve(dirname(fileURLToPath(import.meta.url)),'..');
const load=src=>new Function(src+';return {LGCA,JP,JN,makeWell};')();
const source=readFileSync(join(root,'lgca-viz/lgca-lab.html'),'utf8').match(/<script id="engineSrc">([\s\S]*?)<\/script>/)[1];
const {LGCA,JP,JN,makeWell}=load(source);
const {LGCA:Before}=load(readFileSync(join(root,'results/alpha-2026-09-24/build/engine-source.js'),'utf8'));
let legacy=0;
for(const model of ['boson','fermion'])for(const field of ['site','neigh'])for(const align of ['polar','nematic'])for(const kernel of ['sum','avg','power']){
 const o={model,W:27,H:19,dens:model==='boson'?1.25:.85,seed:4294967295,bosonField:field,bosonAlign:align,kernel,alpha:.85};
 const a=new LGCA(o),b=new Before(o);
 for(let t=0;t<150;t++){a.step(6);b.step(6);assert.deepEqual(a.occ,b.occ);assert.deepEqual(a.rng.getState(),b.rng.getState());}
 legacy++;
}
console.log(`PASS ${legacy} unchanged field/kernel/alignment trajectories against pinned study engine, state + RNG`);
// Directly populated asymmetric seven-site neighbourhood, including periodic
// corners; a populated second-ring site must not enter the field or denominator.
let cdfs=0;
for(const [i,j] of [[0,0],[10,6],[5,3]])for(const align of ['polar','nematic'])for(const alpha of [0,.5,.85,1,1.5,2]){
 const o={model:'boson',W:11,H:7,dens:0,seed:17,kernel:'power',alpha,bosonField:'both',bosonAlign:align};
 const sim=new LGCA(o),put=(x,y,ns)=>sim.occ.set(ns,(((y+7)%7)*11+(x+11)%11)*6);
 const centre=[4,2,0,1,0,0],ring=[[1,0,[2,0,0,1,0,0]],[0,1,[0,3,0,0,0,1]],[-1,0,[0,0,2,0,0,0]],[0,-1,[0,1,0,4,0,0]],[1,-1,[0,0,0,0,3,0]],[-1,1,[1,0,0,0,0,2]]];
 put(i,j,centre);for(const [dx,dy,ns] of ring)put(i+dx,j+dy,ns);put(i+2,j,[100,0,0,0,0,0]);
 const m=centre.map((n,k)=>n+ring.reduce((s,r)=>s+r[2][k],0)),M=m.reduce((s,n)=>s+n,0),J=align==='polar'?JP:JN;
 const w=Array.from({length:6},(_,c)=>Math.exp(.3*m.reduce((s,n,k)=>s+n*J[k][c],0)/M**alpha)),z=w.reduce((s,n)=>s+n,0);
 const rng=makeWell(0);rng.setState(sim.rng.getState());const expected=new Int32Array(6);let cdf=0;const cdfsExpected=w.map(x=>(cdf+=x/z));
 for(let p=0;p<7;p++){const r=rng();let c=5;for(let k=0;k<6;k++)if(r<=cdfsExpected[k]){c=k;break;}expected[c]++;}
 sim.collideBoson(i,j,.3,7);
 for(let c=0;c<6;c++)assert(Math.abs(sim.prob[c]-cdfsExpected[c])<1e-13);
 assert.deepEqual(sim.tmp.slice((j*11+i)*6,(j*11+i+1)*6),expected);assert.deepEqual(sim.rng.getState(),rng.getState());cdfs++;
}
console.log(`PASS ${cdfs} independent seven-site CDF/sample checks, including periodic corners and denominator`);
// With empty neighbours, centre-inclusive and own-site collision coincide.
for(const align of ['polar','nematic'])for(const alpha of [0,.85,1,2]){
 const o={model:'boson',W:11,H:7,dens:0,seed:777,kernel:'power',alpha,bosonAlign:align};
 const own=new LGCA({...o,bosonField:'site'}),both=new LGCA({...o,bosonField:'both'});
 for(const sim of [own,both]){sim.occ[0]=100;sim.nuparts=100;sim.collideBoson(0,0,6,100);}
 assert.deepEqual(own.prob,both.prob);assert.deepEqual(own.tmp,both.tmp);assert.deepEqual(own.rng.getState(),both.rng.getState());
}
// The new option conserves N and consumes exactly one draw per particle per
// step; alpha endpoints use the existing sum/average paths; replay is exact.
for(const alpha of [0,.85,1,2])for(const align of ['polar','nematic']){
 const o={model:'boson',W:31,H:17,dens:2.25,seed:424242,kernel:'power',alpha,bosonField:'both',bosonAlign:align};
 const sim=new LGCA(o),rng=makeWell(0);rng.setState(sim.rng.getState());
 const endpoint=alpha===0||alpha===1?new LGCA({...o,kernel:alpha===0?'sum':'avg'}):null;
 let saved;
 for(let t=1;t<=200;t++){
  sim.step(8);for(let p=0;p<sim.nuparts;p++)rng();
  assert.equal(sim.occ.reduce((s,n)=>s+n,0),sim.nuparts);assert(sim.occ.every(n=>n>=0));assert.deepEqual(sim.rng.getState(),rng.getState());
  if(endpoint){endpoint.step(8);assert.deepEqual(sim.occ,endpoint.occ);assert.deepEqual(sim.rng.getState(),endpoint.rng.getState());}
  if(t===100)saved={occ:sim.occ.slice(),rng:sim.rng.getState()};
 }
 const resumed=new LGCA(o);resumed.occ.set(saved.occ);resumed.rng.setState(saved.rng);resumed.t=100;
 for(let t=100;t<200;t++)resumed.step(8);
 assert.deepEqual(resumed.occ,sim.occ);assert.deepEqual(resumed.rng.getState(),sim.rng.getState());
}
console.log('PASS centre-inclusive isolated-site limit, conservation, RNG count, endpoints, and resumed rectangular trajectories');
