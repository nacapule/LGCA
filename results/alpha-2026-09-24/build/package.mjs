import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {dirname,resolve,join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {encodePNG} from '../../shared/png.mjs';
import {DLUT,renderDirection} from '../../shared/render.mjs';
const here=dirname(fileURLToPath(import.meta.url)),out=resolve(here,'..'),data=join(out,'data');
const manifest=JSON.parse(readFileSync(join(data,'manifest.json')));
mkdirSync(join(out,'images'),{recursive:true});
const runs=manifest.jobs.map(j=>JSON.parse(readFileSync(join(data,j.id,'run.json'))));
for(const r of runs){
 if(r.snaps.at(-1).t!==r.job.end)throw Error('incomplete '+r.job.id);
 const W=r.job.size,H=W;
 for(const m of r.snaps){
  const buf=gunzipSync(readFileSync(join(data,r.job.id,m.file))),occ=new Int32Array(buf.buffer,buf.byteOffset,buf.byteLength/4);
  const w=2*W,h=2*H,rgba=new Uint8Array(w*h*4);for(let q=0;q<w*h;q++)rgba.set([11,14,19,255],q*4);
  const local=new Uint8Array(rgba);
  for(let j=0;j<H;j++)for(let i=0;i<W;i++){
   let n=0,x=0,y=0;const cx=[1,.5,-.5,-1,-.5,.5],cy=[0,Math.sqrt(3)/2,Math.sqrt(3)/2,0,-Math.sqrt(3)/2,-Math.sqrt(3)/2];
   for(let k=0;k<6;k++){const z=occ[(j*W+i)*6+k];n+=z;x+=z*cx[k];y+=z*cy[k];}
   if(!n)continue;
   const c=DLUT[Math.min(255,Math.floor(255*Math.log1p(n)/Math.log1p(1000)))],p=DLUT[Math.min(255,Math.floor(255*Math.hypot(x,y)/n))];
   for(let dy=0;dy<2;dy++)for(let dx=0;dx<2;dx++){const q=((2*j+dy)*w+(2*i+j+dx)%w)*4;rgba.set([...c,255],q);local.set([...p,255],q);}
  }
  const base=`${r.job.id}-t${m.t}`;
  writeFileSync(join(out,'images',base+'-density.png'),encodePNG(w,h,rgba));
  writeFileSync(join(out,'images',base+'-local.png'),encodePNG(w,h,local));
  const d=renderDirection(occ,W,H);writeFileSync(join(out,'images',base+'-direction.png'),encodePNG(d.w,d.h,d.rgba));
 }
}
writeFileSync(join(out,'summary.json'),JSON.stringify({manifest,runs}));
const recipes=runs.map(r=>({model:'boson',size:r.job.size,dens:.4,sens:r.job.sens,seed:r.job.seed,kernel:'power',alpha:r.job.alpha,bosonField:r.job.field,bosonAlign:'polar',name:r.job.id,desc:'Independent alpha pilot; reset, then advance to a recorded step.',when:'2026-09-24'}));
writeFileSync(join(out,'tool-settings.json'),JSON.stringify(recipes,null,2));
const csvkeys=['t','N','mean','nmax','kmax','siteFraction','channelFraction','effectiveArea','polar','localPolar','weightedLocalPolar','threshold','clusterMass','clusterArea','clusterDensity','clusterFraction','clusterPolar'];
writeFileSync(join(out,'measurements.csv'),'id,field,control,size,alpha,sensitivity,seed,'+csvkeys.join(',')+'\n'+runs.flatMap(r=>r.series.map(m=>[r.job.id,r.job.field,r.job.control,r.job.size,r.job.alpha,r.job.sens,r.job.seed,...csvkeys.map(k=>m[k])].join(','))).join('\n'));
let html=readFileSync(join(here,'report-template.html'),'utf8');html=html.replace('/*DATA*/',JSON.stringify(runs));writeFileSync(join(out,'index.html'),html);
console.log(`Packaged ${runs.length} runs, ${runs.reduce((s,r)=>s+r.snaps.length,0)} exact checkpoints, density/direction/local-order images.`);
