import {readFileSync,writeFileSync,mkdirSync,existsSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {dirname,resolve,join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {encodePNG} from '../../shared/png.mjs';
import {DLUT} from '../../shared/render.mjs';
const here=dirname(fileURLToPath(import.meta.url)),out=resolve(here,'..');
const analysis=JSON.parse(readFileSync(join(out,'analysis.json'))),runs=JSON.parse(readFileSync(join(out,'runs.json'))),findings=JSON.parse(readFileSync(join(out,'findings.json')));
const validation=JSON.parse(readFileSync(join(out,'validation.json')));
if(validation.resumeCheck.to!==findings.horizon||createHash('sha256').update(readFileSync(join(out,'analysis.json'))).digest('hex')!==findings.analysisSha256||validation.runs!==runs.length||runs.some(r=>r.snaps.at(-1).t!==findings.horizon))throw Error('Findings or validation do not match the completed run horizon');
mkdirSync(join(out,'images'),{recursive:true});
for(const r of runs)for(const m of r.snaps.filter(m=>m.t>=10000)){
 const L=r.job.size,buf=gunzipSync(readFileSync(join(out,'data',r.job.id,m.file))),occ=new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4),w=L*2,rgba=new Uint8Array(w*w*4);
 for(let q=0;q<w*w;q++)rgba.set([11,14,19,255],q*4);
 for(let j=0;j<L;j++)for(let i=0;i<L;i++){
  let n=0;for(let k=0;k<6;k++)n+=occ[(j*L+i)*6+k];if(!n)continue;
  const col=DLUT[Math.min(255,Math.floor(255*Math.log1p(n)/Math.log1p(1000)))];
  for(let dy=0;dy<2;dy++)for(let dx=0;dx<2;dx++)rgba.set([...col,255],((2*j+dy)*w+(2*i+j+dx)%w)*4);
 }
 writeFileSync(join(out,'images',`${r.job.id}-t${m.t}.png`),encodePNG(w,w,rgba));
}
const data={groups:analysis.groups,runs:runs.map(r=>({job:r.job,snaps:r.snaps.map(({t,nmax,clusterDensity,clusterMass,clusterArea})=>({t,nmax,clusterDensity,clusterMass,clusterArea}))}))};
const gradient=Array.from({length:256},(_,q)=>{
 const logCount=Math.log(2)+(Math.log(1001)-Math.log(2))*q/255;
 const c=DLUT[Math.min(255,Math.floor(255*logCount/Math.log(1001)))];
 return `rgb(${c.join(',')}) ${100*q/255}%`;
}).join(',');
const densityRoot=resolve(out,'../density-2026-09-25');
let density=null;
if(existsSync(join(densityRoot,'display.json'))){
 const v=JSON.parse(readFileSync(join(densityRoot,'validation.json')));
 density=JSON.parse(readFileSync(join(densityRoot,'display.json')));
 if(v.runs!==60||v.horizon!==80000||density.analysisSha256!==v.analysisSha256||createHash('sha256').update(readFileSync(join(densityRoot,'analysis.json'))).digest('hex')!==v.analysisSha256)throw Error('Stale density display');
}
let html=readFileSync(join(here,'report.html.in'),'utf8');html=html.replace('/*DATA*/',JSON.stringify(data)).replace('/*DENSITY_GRADIENT*/',gradient).replace('/*DENSITY_DATA*/',JSON.stringify(density));writeFileSync(join(out,'index.html'),html);console.log('Report and density snapshots packaged.');
