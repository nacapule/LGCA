import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {dirname,resolve,join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {encodePNG} from '../../shared/png.mjs';
import {DLUT} from '../../shared/render.mjs';
const out=resolve(dirname(fileURLToPath(import.meta.url)),'..');
const analysis=JSON.parse(readFileSync(join(out,'analysis.json'))),runs=JSON.parse(readFileSync(join(out,'runs.json'))),validation=JSON.parse(readFileSync(join(out,'validation.json')));
const analysisSha256=createHash('sha256').update(readFileSync(join(out,'analysis.json'))).digest('hex');
if(validation.runs!==60||validation.horizon!==80000||validation.analysisSha256!==analysisSha256||runs.some(r=>r.snaps.at(-1).t!==80000||r.engineHash!==validation.engineHash))throw Error('Unvalidated or incomplete density data');
mkdirSync(join(out,'images'),{recursive:true});
for(const r of runs)for(const m of r.snaps.filter(m=>[10000,20000,40000,80000].includes(m.t))){
 const L=r.job.size,b=gunzipSync(readFileSync(join(out,'data',r.job.id,m.file))),occ=new Int32Array(b.buffer,b.byteOffset,b.byteLength/4),w=L*2,rgba=new Uint8Array(w*w*4);
 for(let q=0;q<w*w;q++)rgba.set([11,14,19,255],q*4);
 for(let j=0;j<L;j++)for(let i=0;i<L;i++){
  let n=0;for(let k=0;k<6;k++)n+=occ[(j*L+i)*6+k];if(!n)continue;
  const c=DLUT[Math.min(255,Math.floor(255*Math.log1p(n)/Math.log1p(1000)))];
  for(let dy=0;dy<2;dy++)for(let dx=0;dx<2;dx++)rgba.set([...c,255],((2*j+dy)*w+(2*i+j+dx)%w)*4);
 }
 writeFileSync(join(out,'images',`${r.job.id}-t${m.t}.png`),encodePNG(w,w,rgba));
}
writeFileSync(join(out,'display.json'),JSON.stringify({analysisSha256,groups:analysis.groups,runs:runs.map(r=>({job:r.job,snaps:r.snaps.map(({t,nmax,clusterDensity,clusterMass,clusterArea})=>({t,nmax,clusterDensity,clusterMass,clusterArea}))}))}));
console.log('Validated density data and snapshots packaged.');
