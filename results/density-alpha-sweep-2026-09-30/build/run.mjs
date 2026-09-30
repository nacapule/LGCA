// Expanded density x alpha sweep (Manik, Sep 25 request): many starting densities, same grid for alpha 0..1.
// Same physics, measurement schedule and checkpoints as results/density-2026-09-25 (pinned engine, centre+neighbours).
// Staged: `node run.mjs --end=20000 --seeds=12345` runs/extends every job to `end`; later stages resume exactly
// from the last saved state+RNG. Matching completed runs from density-2026-09-25 are copied, not re-simulated.
import {readFileSync,writeFileSync,mkdirSync,existsSync,renameSync,cpSync} from 'node:fs';
import {gzipSync,gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {dirname,resolve,join} from 'node:path';
import {measurements} from '../../alpha-2026-09-24/build/measure.mjs';
const here=dirname(fileURLToPath(import.meta.url)),out=resolve(here,'..'),data=join(out,'data');
const oldData=resolve(out,'../density-2026-09-25/data');
const PINNED='34edbb46d96ee9aaea500e8c5e0927a5858366a2ed52fc3db1284d37b7d5ff6d';
const source=readFileSync(join(here,'engine-source.js'),'utf8'),engineHash=createHash('sha256').update(source).digest('hex');
if(engineHash!==PINNED)throw Error('Engine differs from pinned source');
const {LGCA}=new Function(source+';return {LGCA};')();
// Per-channel density inputs; mean particles/site = 6x. Existing .1/.2/.4/.8 retained.
const densities=[.025,.05,.1,.15,.2,.3,.4,.6,.8,1.2,1.6,2.4];
const baseAlphas=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1];
// Refinement where the coarse pass (Sep 30, seed 12345, 20k) changes behaviour: between .8 and 1.
const refineAlphas=[.75,.81,.82,.83,.84,.85,.86,.87,.88,.89,.92,.94,.96,.98];
const allSeeds=[12345,777,424242,20260925,314159],size=90,sens=6,field='both';
const ticks=[0,100,500,1000,2000,5000,10000,20000,40000,60000,80000];
const isSample=t=>t<=5000?t%25===0:t%100===0;
const arg=k=>{const a=process.argv.find(x=>x.startsWith(`--${k}=`));return a&&a.slice(k.length+3);};
const idOf=(alpha,dens,seed)=>`both-L${size}-a${alpha}-d${dens}-seed${seed}`;
const writeJSON=(path,value)=>{writeFileSync(path+'.tmp',JSON.stringify(value));renameSync(path+'.tmp',path);};
const reached=dir=>{const p=join(dir,'run.json');if(!existsSync(p))return -1;const r=JSON.parse(readFileSync(p));if(r.engineHash!==engineHash)throw Error('Engine changed: '+dir);return r.snaps.at(-1).t;};

if(process.argv[2]==='job'){
 const [,, ,id,endArg]=process.argv,end=+endArg;
 const m=id.match(/^both-L90-a([\d.]+)-d([\d.]+)-seed(\d+)$/);if(!m)throw Error('Bad id');
 const job={id,size,alpha:+m[1],dens:+m[2],seed:+m[3],sens,field};
 const dir=join(data,id);mkdirSync(dir,{recursive:true});
 const sim=new LGCA({model:'boson',W:size,H:size,dens:job.dens,seed:job.seed,kernel:'power',alpha:job.alpha,bosonField:field,bosonAlign:'polar'});
 const path=join(dir,'run.json');let series=[],snaps=[],start=0,previousSeconds=0;
 if(existsSync(path)){
  const old=JSON.parse(readFileSync(path));if(old.engineHash!==engineHash)throw Error('Engine changed');
  for(const k of ['size','alpha','dens','seed','sens','field'])if(old.job[k]!==job[k])throw Error('Settings changed');
  const last=old.snaps.at(-1),buf=gunzipSync(readFileSync(join(dir,last.file)));
  if(createHash('sha256').update(buf).digest('hex')!==last.sha256)throw Error('State hash mismatch');
  sim.occ.set(new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4));sim.rng.setState(new Uint32Array(last.rng));sim.t=last.t;
  if(last.N!==sim.nuparts)throw Error('Particle count mismatch');
  start=last.t+1;series=old.series.filter(x=>x.t<=last.t);snaps=old.snaps;previousSeconds=old.seconds||0;
 }
 const begin=Date.now();
 for(let t=start;t<=end;t++){
  if(t)sim.step(sens);
  if(isSample(t)){
   const q=measurements(sim.occ,size,size,2,false);
   if(q.N!==sim.nuparts||sim.occ.some(n=>n<0))throw Error('Invalid particle count');
   delete q.histogram;
   if(ticks.includes(t)){
    const buf=Buffer.from(sim.occ.buffer),file=`state-t${t}.bin.gz`;writeFileSync(join(dir,file),gzipSync(buf));
    const checks=[1.5,4].map(f=>{const c=measurements(sim.occ,size,size,f);return {factor:f,threshold:c.threshold,mass:c.clusterMass,area:c.clusterArea,density:c.clusterDensity};});
    snaps.push({t,file,sha256:createHash('sha256').update(buf).digest('hex'),rng:Array.from(sim.rng.getState()),...q,thresholdChecks:checks});
   }
   series.push({t,...q});
   if(ticks.includes(t))writeJSON(path,{job:{...job,end},engineHash,model:'boson',alignment:'polar',stage:'after streaming',seconds:previousSeconds+(Date.now()-begin)/1000,series,snaps});
  }
 }
 console.log(`${id}: t=${end}, ${(previousSeconds+(Date.now()-begin)/1000).toFixed(1)}s`);
}else{
 const end=+(arg('end')||20000),workers=+(arg('workers')||6);
 const seeds=(arg('seeds')||'12345').split(',').map(Number),set=arg('set')||'base',alphas=set==='refine'?refineAlphas:baseAlphas;
 if(!ticks.includes(end))throw Error('end must be a checkpoint');
 mkdirSync(data,{recursive:true});
 writeJSON(join(data,'manifest.json'),{created:'2026-09-30',engineHash,densities,baseAlphas,refineAlphas,allSeeds,size,sens,field,model:'boson',alignment:'polar',kernel:'power',checkpoints:ticks,measurementSchedule:'every 25 steps through 5000, every 100 thereafter',question:'Cluster density versus time for many starting densities, repeated across alpha 0..1 (centre plus neighbours, 90x90)'});
 const jobs=[];
 for(const seed of seeds)for(const dens of [...densities].reverse())for(const alpha of alphas)jobs.push(idOf(alpha,dens,seed));
 let next=0,done=0,copied=0;const t0=Date.now();
 async function worker(){while(next<jobs.length){const id=jobs[next++],dir=join(data,id),old=join(oldData,id);
  if(reached(dir)<0&&existsSync(join(old,'run.json'))&&reached(old)>=0){cpSync(old,dir,{recursive:true});copied++;}
  if(reached(dir)>=end){done++;continue;}
  await new Promise((res,rej)=>{const p=spawn(process.execPath,[fileURLToPath(import.meta.url),'job',id,String(end)],{stdio:'inherit'});p.on('error',rej);p.on('exit',c=>c===0?res():rej(Error(id+' failed')));});
  done++;console.log(`${done}/${jobs.length} complete (${((Date.now()-t0)/60000).toFixed(1)} min)`);
 }}
 await Promise.all(Array.from({length:workers},worker));
 console.log(`Stage ${set} end=${end}, seeds ${seeds}: all ${jobs.length} jobs reached ${end} (${copied} copied from density-2026-09-25).`);
}
