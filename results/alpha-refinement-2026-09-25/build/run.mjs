import {readFileSync,writeFileSync,mkdirSync,existsSync,renameSync} from 'node:fs';
import {gzipSync,gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {dirname,resolve,join} from 'node:path';
import {measurements} from '../../alpha-2026-09-24/build/measure.mjs';
const here=dirname(fileURLToPath(import.meta.url)),root=resolve(here,'../../..'),out=resolve(here,'..'),data=join(out,'data');
const source=readFileSync(join(root,'results/alpha-2026-09-24/build/engine-source.js'),'utf8');
const {LGCA}=new Function(source+';return {LGCA};')();
const engineHash=createHash('sha256').update(source).digest('hex');
const horizon=Number(process.env.LGCA_REFINEMENT_END||80000);
if(![20000,40000,80000].includes(horizon))throw Error('unsupported horizon');
const seeds=[12345,777,424242,20260925,314159],jobs=[];
for(const size of [90,120,180])for(const alpha of [.75,.8,.85,.9])for(const seed of seeds)
 jobs.push({id:`neigh-L${size}-a${alpha}-d0.4-seed${seed}`,size,alpha,seed,dens:.4,field:'neigh',sens:6,end:horizon});
const writeJSON=(path,value)=>{writeFileSync(path+'.tmp',JSON.stringify(value));renameSync(path+'.tmp',path);};
if(process.argv[2]==='job'){
 const job=jobs.find(j=>j.id===process.argv[3]);if(!job)throw Error('job missing');
 const dir=join(data,job.id);mkdirSync(dir,{recursive:true});
 const sim=new LGCA({model:'boson',W:job.size,H:job.size,dens:job.dens,seed:job.seed,kernel:'power',alpha:job.alpha,bosonField:job.field,bosonAlign:'polar'});
 const ticks=[0,1000,5000,10000,20000,40000,80000].filter(t=>t<=horizon);let series=[],snaps=[],start=0,previousSeconds=0;
 const path=join(dir,'run.json');
 if(existsSync(path)){
  const previous=JSON.parse(readFileSync(path));if(previous.engineHash!==engineHash)throw Error('engine changed: refuse resume');
  for(const k of ['size','alpha','seed','dens','field','sens'])if(previous.job[k]!==job[k])throw Error('settings changed');
  const last=previous.snaps.at(-1),buf=gunzipSync(readFileSync(join(dir,last.file)));
  if(createHash('sha256').update(buf).digest('hex')!==last.sha256)throw Error('checkpoint hash mismatch');
  sim.occ.set(new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4));sim.rng.setState(new Uint32Array(last.rng));sim.t=last.t;
  if(last.N!==sim.nuparts)throw Error('checkpoint particle count mismatch');
  start=last.t+1;series=previous.series.filter(m=>m.t<=last.t);snaps=previous.snaps;previousSeconds=previous.seconds||0;
 }
 const t0=Date.now();
 for(let t=start;t<=horizon;t++){
  if(t)sim.step(job.sens);
  if(t%100===0){
   const full=ticks.includes(t),m=measurements(sim.occ,job.size,job.size,2,false);
   if(m.N!==sim.nuparts||sim.occ.some(n=>n<0))throw Error('invalid particle count');
   if(full){
    const state=Buffer.from(sim.occ.buffer),file=`state-t${t}.bin.gz`;writeFileSync(join(dir,file),gzipSync(state));
    snaps.push({t,file,sha256:createHash('sha256').update(state).digest('hex'),rng:Array.from(sim.rng.getState()),...m,thresholdChecks:[1.5,4].map(f=>{const z=measurements(sim.occ,job.size,job.size,f);return {factor:f,threshold:z.threshold,mass:z.clusterMass,area:z.clusterArea,density:z.clusterDensity};})});
   }
   delete m.histogram;series.push({t,...m});
   if(full)writeJSON(path,{job,engineHash,model:'boson',alignment:'polar',stage:'after streaming',seconds:previousSeconds+(Date.now()-t0)/1000,series,snaps});
  }
 }
 console.log(`${job.id}: t=${horizon}, ${(previousSeconds+(Date.now()-t0)/1000).toFixed(1)}s`);
}else{
 mkdirSync(data,{recursive:true});
 for(const job of jobs){const path=join(data,job.id,'run.json');if(existsSync(path)&&JSON.parse(readFileSync(path)).snaps.at(-1).t>horizon)throw Error('Refuse to lower the recorded horizon; use a separate output directory for a shorter study');}
 writeJSON(join(data,'manifest.json'),{created:'2026-09-25',engineHash,seeds,jobs,measureEvery:100,question:'Does the alpha .75–.90 concentration crossover persist with size and observation time at density .4?',densityPolicy:'Hold density fixed for finite-size comparison; a density scan is a separate parameter axis.'});
 let next=0,done=0;
 async function worker(){while(next<jobs.length){const job=jobs[next++],path=join(data,job.id,'run.json');
  if(existsSync(path)){const old=JSON.parse(readFileSync(path));if(old.engineHash!==engineHash)throw Error('engine changed');if(old.snaps.at(-1).t>=horizon){done++;continue;}}
  await new Promise((res,rej)=>{const p=spawn(process.execPath,[fileURLToPath(import.meta.url),'job',job.id],{stdio:'inherit'});p.on('error',rej);p.on('exit',code=>code===0?res():rej(Error(job.id+' failed')));});
  done++;console.log(`${done}/${jobs.length} complete`);
 }}
 await Promise.all([worker(),worker(),worker()]);console.log(`All ${jobs.length} refinement runs complete through ${horizon}.`);
}
