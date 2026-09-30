// Descriptive statistics; no simulation RNG. Axial periodic six-neighbour graph.
export const NB=[[1,0],[0,1],[-1,0],[0,-1],[1,-1],[-1,1]];
export function measurements(occ,W,H,thresholdMultiplier=2,profiles=false){
 const A=W*H,n=new Int32Array(A),px=new Float64Array(A),py=new Float64Array(A);
 const cx=[1,.5,-.5,-1,-.5,.5],cy=[0,Math.sqrt(3)/2,Math.sqrt(3)/2,0,-Math.sqrt(3)/2,-Math.sqrt(3)/2];
 let N=0,nmax=0,kmax=0,sum2=0,occupied=0,lpol=0,weighted=0,gx=0,gy=0;
 const hist={};
 for(let x=0;x<A;x++){
  for(let k=0;k<6;k++){const v=occ[6*x+k];n[x]+=v;px[x]+=v*cx[k];py[x]+=v*cy[k];kmax=Math.max(kmax,v);}
  N+=n[x];sum2+=n[x]*n[x];nmax=Math.max(nmax,n[x]);hist[n[x]]=(hist[n[x]]||0)+1;
  const p=Math.hypot(px[x],py[x]);gx+=px[x];gy+=py[x];weighted+=p;if(n[x]){occupied++;lpol+=p/n[x];}
 }
 const threshold=Math.max(1,Math.ceil(thresholdMultiplier*N/A));
 const neigh=(x,d)=>((Math.floor(x/W)+NB[d][1]+H)%H)*W+(x%W+NB[d][0]+W)%W;
 const seen=new Uint8Array(A),queue=new Int32Array(A);let best=[],mass=0,clusters=0;
 for(let x=0;x<A;x++)if(!seen[x]&&n[x]>=threshold){
  let h=0,t=1,m=0;queue[0]=x;seen[x]=1;
  while(h<t){const v=queue[h++];m+=n[v];for(let d=0;d<6;d++){const z=neigh(v,d);if(!seen[z]&&n[z]>=threshold){seen[z]=1;queue[t++]=z;}}}
  clusters++;if(m>mass){mass=m;best=Array.from(queue.subarray(0,t));}
 }
 let bx=0,by=0;for(const x of best){bx+=px[x];by+=py[x];}
 const result={N,mean:N/A,nmax,kmax,siteFraction:N?nmax/N:0,channelFraction:N?kmax/N:0,effectiveArea:sum2?N*N/sum2:0,occupied,polar:N?Math.hypot(gx,gy)/N:0,localPolar:occupied?lpol/occupied:0,weightedLocalPolar:N?weighted/N:0,threshold,clusters,clusterMass:mass,clusterArea:best.length,clusterDensity:best.length?mass/best.length:0,clusterFraction:N?mass/N:0,clusterPolar:mass?Math.hypot(bx,by)/mass:0,histogram:hist};
 if(profiles&&best.length&&best.length<A){
  const mask=new Uint8Array(A);best.forEach(x=>mask[x]=1);
  const dist=new Int32Array(A);dist.fill(-1);let h=0,t=0;
  for(const x of best)if(NB.some((_,d)=>!mask[neigh(x,d)])){dist[x]=0;queue[t++]=x;}
  while(h<t){const x=queue[h++];for(let d=0;d<6;d++){const z=neigh(x,d);if(dist[z]<0){dist[z]=dist[x]+1;queue[t++]=z;}}}
  const bins={};for(let x=0;x<A;x++){
   const d=mask[x]?-dist[x]:dist[x];if(Math.abs(d)>12)continue;
   const b=bins[d]||(bins[d]={distance:d,sites:0,mass:0,occupied:0,localSum:0,weightedSum:0});
   b.sites++;b.mass+=n[x];if(n[x]){b.occupied++;const norm=Math.hypot(px[x],py[x]);b.localSum+=norm/n[x];b.weightedSum+=norm;}
  }
  result.edgeProfile=Object.values(bins).sort((a,b)=>a.distance-b.distance).map(b=>({...b,density:b.mass/b.sites,localPolar:b.occupied?b.localSum/b.occupied:null,weightedPolar:b.mass?b.weightedSum/b.mass:null}));
 }
 return result;
}
