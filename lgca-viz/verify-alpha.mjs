#!/usr/bin/env node
import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,mkdtempSync,rmSync} from 'node:fs';
import {execFileSync} from 'node:child_process';
import {tmpdir} from 'node:os';
import {resolve,join,dirname} from 'node:path';
import {fileURLToPath} from 'node:url';
const root=resolve(dirname(fileURLToPath(import.meta.url)),'..');
const load=p=>new Function(readFileSync(p,'utf8').match(/<script id="engineSrc">([\s\S]*?)<\/script>/)[1]+';return {LGCA,JP,JN};')();
const {LGCA,JP}=load(join(root,'lgca-viz/lgca-lab.html'));
const {LGCA:Before}=load(join(root,'results/alpha-2026-09-24/build/lgca-lab.before.html'));
let endpoints=0;
for(const model of ['boson','fermion']) for(const field of ['site','neigh']) for(const align of ['polar','nematic']) for(const a of [0,1]){
  const o={model,W:24,H:18,dens:model==='boson'?2.25:0.85,seed:54321,bosonField:field,bosonAlign:align};
  const old=new Before({...o,kernel:a?'avg':'sum'}), sim=new LGCA({...o,kernel:'power',alpha:a});
  for(let t=0;t<120;t++) { old.step(6);sim.step(6);assert.deepEqual(sim.occ,old.occ);assert.deepEqual(sim.rng.getState(),old.rng.getState()); }
  endpoints++;
}
console.log(`PASS ${endpoints} endpoint/legacy trajectories, 120 steps each, rectangular lattice`);
for(const a of [0.25,0.5,0.75,1.5,2]) for(const model of ['boson','fermion']) for(const field of ['site','neigh']){
 const o={model,W:24,H:18,dens:model==='boson'?1.2:0.7,seed:42,kernel:'power',alpha:a,bosonField:field,bosonAlign:'polar'};
 const sim=new LGCA(o), zero=new LGCA(o), sum=new LGCA({...o,kernel:'sum'});
 for(let t=0;t<100;t++){
  sim.step(8);zero.step(0);sum.step(0);
  assert.deepEqual(zero.occ,sum.occ);assert.deepEqual(zero.rng.getState(),sum.rng.getState());
  assert.equal(sim.occ.reduce((s,n)=>s+n,0),sim.nuparts);assert(sim.occ.every(n=>n>=0&&(model!=='fermion'||n<=1)));
 }
}
// Independent CDF check: asymmetric local field, explicit exp(score)/Z.
for(const a of [0,0.5,1,1.5]){
 const sim=new LGCA({model:'boson',W:3,H:3,dens:0,seed:1,kernel:'power',alpha:a,bosonAlign:'polar'});
 const ns=[4,2,0,1,0,0]; sim.occ.set(ns,0);sim.collideBoson(0,0,0.7,7);
 const w=JP.map((_,c)=>Math.exp(0.7*ns.reduce((s,n,k)=>s+n*JP[k][c],0)/Math.pow(7,a)));
 const z=w.reduce((s,n)=>s+n,0);let cdf=0;
 for(let c=0;c<6;c++){cdf+=w[c]/z;assert(Math.abs(sim.prob[c]-cdf)<1e-14);}
}
// Isolated site with neighbours-only has empty interaction field: always uniform.
const iso=new LGCA({model:'boson',W:7,H:5,dens:0,seed:7,kernel:'power',alpha:0.5,bosonField:'neigh',bosonAlign:'polar'});
iso.occ[0]=100;iso.collideBoson(0,0,8,100);
for(let c=0;c<6;c++) assert(Math.abs(iso.prob[c]-(c+1)/6)<1e-14);
console.log('PASS intermediate-alpha conservation, exclusion, zero-sensitivity equality, analytic CDF, empty-neighbour field');
// Real C++: every 50th streamed state and all 33 WELL words, stress + rectangular.
const temp=mkdtempSync(join(tmpdir(),'lgca-adversarial-'));
try{
 const rng=readFileSync(join(root,'rng/WELL1024a.c'),'utf8');
 const cleanup=['W','R','M1','M2','M3','MAT0POS','MAT0NEG','Identity','V0','VM1','VM2','VM3','VRm1','newV0','newV1','FACT'].map(x=>'#undef '+x).join('\n');
 let checkpoints=0;
 for(const model of ['boson','fermion']){
  let cpp=readFileSync(join(root,'lgca/lgca_clean-1.cpp'),'utf8').replace('XDIM  = 120','XDIM  = 96').replace('YDIM  = 120','YDIM  = 60');
  if(model==='fermion')cpp=cpp.replace(/^#define MODEL\s+BOSON.*$/m,'#define MODEL FERMION');
  const probe=`if((t+1)%50==0){ std::cout << "CHECK " << t+1 << " " << state_i; for(int q=0;q<32;q++)std::cout<<" "<<STATE[q]; for(int i=0;i<XDIM;i++)for(int j=0;j<YDIM;j++)for(int k=0;k<6;k++)std::cout<<" "<<sim.lattice[i][j][k];std::cout<<"\\n";}\n`;
  cpp=cpp.replace('            // (3) OBSERVABLES.',probe+'            // (3) OBSERVABLES.');
  const path=join(temp,model); writeFileSync(path+'.cpp',rng+'\n'+cleanup+'\n'+cpp);
  execFileSync('c++',['-std=c++17','-O2','-I',join(root,'lgca'),path+'.cpp','-o',path]);
  for(const dens of model==='boson'?[0.4,2.25]:[0.2,0.95]){
   const lines=execFileSync(path,['300','1','8',String(dens),'4294967295'],{encoding:'utf8',maxBuffer:30e6}).split('\n').filter(x=>x.startsWith('CHECK '));
   const sim=new LGCA({model,W:96,H:60,dens,seed:4294967295});
   for(let t=1;t<=300;t++){
    sim.step(8);
    if(t%50)continue;
    const a=lines[t/50-1].split(' ').slice(1).map(Number); assert.equal(a.shift(),t);
    const idx=a.shift(),state=a.splice(0,32); assert.deepEqual(Array.from(sim.rng.getState()),[...state,idx]);
    let q=0;for(let i=0;i<96;i++)for(let j=0;j<60;j++)for(let k=0;k<6;k++)assert.equal(sim.occ[(j*96+i)*6+k],a[q++]);
    checkpoints++;
   }
  }
 }
 console.log(`PASS ${checkpoints} real C++ state/RNG checkpoints over 1200 steps, 96×60 stress scenarios`);
} finally {rmSync(temp,{recursive:true,force:true});}
