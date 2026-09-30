import {readFileSync,writeFileSync,mkdirSync,existsSync} from 'node:fs';
import {gzipSync,gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {dirname,resolve,join} from 'node:path';
import {measurements} from './measure.mjs';
const here=dirname(fileURLToPath(import.meta.url)),root=resolve(here,'../../..'),data=resolve(here,'../data');
const source=readFileSync(join(root,'results/alpha-2026-09-24/build/engine-source.js'),'utf8');
const {LGCA}=new Function(source+';return {LGCA};')();
const engineHash=createHash('sha256').update(source).digest('hex');
const jobs=[];
function add(size,alpha,seed,field='neigh',control='fixed-sensitivity'){
 const sens=control==='matched-field'?6/Math.pow(14.4,1-alpha):6;
 const end=control==='fixed-sensitivity'?10000:2000;
 jobs.push({id:`${field}-${control}-L${size}-a${alpha}-seed${seed}`,size,alpha,seed,field,control,sens,end});
}
for(const size of [60,120,180])for(const alpha of [0,.5,1])for(const seed of [12345,777,424242])add(size,alpha,seed);
for(const alpha of [0,.25,.5,.75,.9,1])for(const seed of [12345,777,424242])add(90,alpha,seed);
for(const alpha of [0,.5,1])for(const seed of [12345,777,424242])add(60,alpha,seed,'site');
for(const alpha of [0,.5])for(const seed of [12345,777,424242])add(60,alpha,seed,'neigh','matched-field');
if(process.argv[2]==='job'){
 const job=jobs.find(j=>j.id===process.argv[3]);if(!job)throw Error('job missing');
 const dir=join(data,job.id);mkdirSync(dir,{recursive:true});
 const sim=new LGCA({model:'boson',W:job.size,H:job.size,dens:.4,seed:job.seed,kernel:'power',alpha:job.alpha,bosonField:job.field,bosonAlign:'polar'});
 const ticks=[0,250,1000,2000,...(job.end>2000?[5000,10000]:[])]; let series=[],snaps=[],start=0;
 if(existsSync(join(dir,'run.json'))){
  const previous=JSON.parse(readFileSync(join(dir,'run.json'),'utf8'));
  if(previous.engineHash!==engineHash)throw Error('engine changed: refuse resume');
  series=previous.series;snaps=previous.snaps;const last=snaps.at(-1);
  const buf=gunzipSync(readFileSync(join(dir,last.file)));
  sim.occ.set(new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4));sim.rng.setState(new Uint32Array(last.rng));sim.t=last.t;start=last.t+1;
 }
 const t0=Date.now();
 for(let t=start;t<=job.end;t++){
  if(t)sim.step(job.sens);
  if(t%50===0){
   const m=measurements(sim.occ,job.size,job.size,2,ticks.includes(t));
   if(m.N!==sim.nuparts)throw Error('mass not conserved');
   if(ticks.includes(t)){
    const state=Buffer.from(sim.occ.buffer);const filename=`state-t${t}.bin.gz`;writeFileSync(join(dir,filename),gzipSync(state));
    snaps.push({t,file:filename,sha256:createHash('sha256').update(state).digest('hex'),rng:Array.from(sim.rng.getState()),...m,thresholdChecks:[1.5,4].map(f=>{const z=measurements(sim.occ,job.size,job.size,f);return {factor:f,threshold:z.threshold,mass:z.clusterMass,area:z.clusterArea,density:z.clusterDensity};})});
   }
   delete m.histogram;delete m.edgeProfile;series.push({t,...m});
  }
 }
 writeFileSync(join(dir,'run.json'),JSON.stringify({job,engineHash,model:'boson',dens:.4,alignment:'polar',stage:'after streaming',stateEncoding:'Int32 little-endian, (j*W+i)*6+k',seconds:(Date.now()-t0)/1000,series,snaps}));
 console.log(job.id+' finished '+((Date.now()-t0)/1000).toFixed(1)+'s');
}else{
 writeFileSync(join(data,'manifest.json'),JSON.stringify({created:'2026-09-24',engineHash,jobs},null,2));
 let next=0,done=0;
 async function worker(){while(next<jobs.length){const job=jobs[next++];if(existsSync(join(data,job.id,'run.json'))&&JSON.parse(readFileSync(join(data,job.id,'run.json'),'utf8')).snaps.at(-1).t>=job.end){done++;continue;}
 await new Promise((res,rej)=>{const p=spawn(process.execPath,[fileURLToPath(import.meta.url),'job',job.id],{stdio:'inherit'});p.on('exit',code=>code===0?res():rej(Error(job.id+' failed')));});done++;console.log(`${done}/${jobs.length}`);}}
 await Promise.all([worker(),worker(),worker()]);
 console.log('All '+jobs.length+' pilot runs complete.');
}
