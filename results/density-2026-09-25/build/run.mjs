import {readFileSync,writeFileSync,mkdirSync,existsSync,renameSync} from 'node:fs';
import {gzipSync,gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {dirname,resolve,join} from 'node:path';
import {measurements} from '../../alpha-2026-09-24/build/measure.mjs';
const here=dirname(fileURLToPath(import.meta.url)),out=resolve(here,'..'),data=join(out,'data');
const source=readFileSync(join(here,'engine-source.js'),'utf8'),engineHash=createHash('sha256').update(source).digest('hex');
const {LGCA}=new Function(source+';return {LGCA};')();
const end=80000,seeds=[12345,777,424242,20260925,314159],densities=[.1,.2,.4,.8],alphas=[.8,.9,1],size=90;
const jobs=[];
// Interleave conditions so saved early progress covers the parameter grid.
for(const seed of seeds)for(const alpha of alphas)for(const dens of densities)jobs.push({id:`both-L${size}-a${alpha}-d${dens}-seed${seed}`,size,alpha,dens,seed,sens:6,field:'both',end});
const writeJSON=(path,value)=>{writeFileSync(path+'.tmp',JSON.stringify(value));renameSync(path+'.tmp',path);};
const ticks=[0,100,500,1000,2000,5000,10000,20000,40000,60000,80000];
const isSample=t=>t<=5000?t%25===0:t%100===0;
if(process.argv[2]==='job'){
 const job=jobs.find(j=>j.id===process.argv[3]);if(!job)throw Error('Job missing');
 const dir=join(data,job.id);mkdirSync(dir,{recursive:true});
 const sim=new LGCA({model:'boson',W:job.size,H:job.size,dens:job.dens,seed:job.seed,kernel:'power',alpha:job.alpha,bosonField:job.field,bosonAlign:'polar'});
 const path=join(dir,'run.json');let series=[],snaps=[],start=0,previousSeconds=0;
 if(existsSync(path)){
  const old=JSON.parse(readFileSync(path));if(old.engineHash!==engineHash)throw Error('Engine changed');
  for(const k of ['size','alpha','dens','seed','sens','field','end'])if(old.job[k]!==job[k])throw Error('Settings changed');
  const last=old.snaps.at(-1),buf=gunzipSync(readFileSync(join(dir,last.file)));
  if(createHash('sha256').update(buf).digest('hex')!==last.sha256)throw Error('State hash mismatch');
  sim.occ.set(new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4));sim.rng.setState(new Uint32Array(last.rng));sim.t=last.t;
  if(last.N!==sim.nuparts)throw Error('Particle count mismatch');
  start=last.t+1;series=old.series.filter(x=>x.t<=last.t);snaps=old.snaps;previousSeconds=old.seconds||0;
 }
 const begin=Date.now();
 for(let t=start;t<=end;t++){
  if(t)sim.step(job.sens);
  if(isSample(t)){
   const m=measurements(sim.occ,job.size,job.size,2,false);
   if(m.N!==sim.nuparts||sim.occ.some(n=>n<0))throw Error('Invalid particle count');
   if(ticks.includes(t)){
    const buf=Buffer.from(sim.occ.buffer),file=`state-t${t}.bin.gz`;writeFileSync(join(dir,file),gzipSync(buf));
    const checks=[1.5,4].map(f=>{const q=measurements(sim.occ,job.size,job.size,f);return {factor:f,threshold:q.threshold,mass:q.clusterMass,area:q.clusterArea,density:q.clusterDensity};});
    delete m.histogram;snaps.push({t,file,sha256:createHash('sha256').update(buf).digest('hex'),rng:Array.from(sim.rng.getState()),...m,thresholdChecks:checks});
   }
   delete m.histogram;series.push({t,...m});
   if(ticks.includes(t))writeJSON(path,{job,engineHash,model:'boson',alignment:'polar',stage:'after streaming',seconds:previousSeconds+(Date.now()-begin)/1000,series,snaps});
  }
 }
 console.log(`${job.id}: t=${end}, ${(previousSeconds+(Date.now()-begin)/1000).toFixed(1)}s`);
}else{
 mkdirSync(data,{recursive:true});
 const manifest={created:'2026-09-25',engineHash,seeds,densities,alphas,size,jobs,checkpoints:ticks,measurementSchedule:'every 25 steps through 5000, every 100 thereafter',question:'Cluster density rise and late-time behaviour versus mean density, at fixed box size and centre-plus-neighbours interactions'};
 const manifestPath=join(data,'manifest.json');
 if(existsSync(manifestPath)&&JSON.stringify(JSON.parse(readFileSync(manifestPath)))!==JSON.stringify(manifest))throw Error('Existing manifest differs');
 writeJSON(manifestPath,manifest);
 let next=0,done=0;
 async function worker(){while(next<jobs.length){const j=jobs[next++],p=join(data,j.id,'run.json');if(existsSync(p)){const r=JSON.parse(readFileSync(p));if(r.engineHash!==engineHash)throw Error('Engine changed');if(r.snaps.at(-1).t===end){done++;continue;}}
  await new Promise((resolve,reject)=>{const p=spawn(process.execPath,[fileURLToPath(import.meta.url),'job',j.id],{stdio:'inherit'});p.on('error',reject);p.on('exit',c=>c===0?resolve():reject(Error(j.id+' failed')));});done++;console.log(`${done}/${jobs.length} complete`);
 }}
 await Promise.all([worker(),worker(),worker()]);console.log(`All ${jobs.length} density runs complete through ${end}.`);
}
