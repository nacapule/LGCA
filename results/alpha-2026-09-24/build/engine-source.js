
"use strict";
function makeWell(seed){
  const S = new Uint32Array(32);
  let s = seed >>> 0;
  for(let j=0;j<32;j++){ S[j]=s; s=(Math.imul(s,1664525)+1013904223)>>>0; }
  let si = 0;
  const FACT = 2.32830643653869628906e-10;
  const gen = function(){
    const z0 = S[(si+31)&31];
    const v0 = S[si], vm1 = S[(si+3)&31], vm2 = S[(si+24)&31], vm3 = S[(si+10)&31];
    const z1 = (v0 ^ (vm1 ^ (vm1>>>8))) >>> 0;
    const z2 = ((vm2 ^ ((vm2<<19)>>>0)) ^ (vm3 ^ ((vm3<<14)>>>0))) >>> 0;
    S[si] = (z1 ^ z2) >>> 0;
    S[(si+31)&31] = ((z0 ^ ((z0<<11)>>>0)) ^ (z1 ^ ((z1<<7)>>>0)) ^ (z2 ^ ((z2<<13)>>>0))) >>> 0;
    si = (si+31)&31;
    return S[si] * FACT;
  };
  gen.getState = function(){ const c=new Uint32Array(33); c.set(S); c[32]=si; return c; };
  gen.setState = function(c){ S.set(c.subarray(0,32)); si=c[32]; };
  return gen;
}

/* Channel geometry exactly as the C++ builds it (sqrt(3)/2 computed, not a literal,
   so the J tables match the C++ doubles bit for bit). */
const EH = Math.sqrt(3.0)/2.0;
const CVEC = [[1.0,0.0],[0.5,EH],[-0.5,EH],[-1.0,0.0],[-0.5,-EH],[0.5,-EH]];
const CX = CVEC.map(v=>v[0]), CY = CVEC.map(v=>v[1]);
const JP = [], JN = [];
for(let i=0;i<6;i++){
  JP.push(new Float64Array(6)); JN.push(new Float64Array(6));
  for(let j=0;j<6;j++){
    const dot = CVEC[i][0]*CVEC[j][0] + CVEC[i][1]*CVEC[j][1];
    JP[i][j] = dot; JN[i][j] = dot*dot;
  }
}
const DST = [[1,0],[1,-1],[0,-1],[-1,0],[-1,1],[0,1]];
const NB  = [[1,0],[0,1],[-1,0],[0,-1],[1,-1],[-1,1]];
const METROPOLIS_STEPS = 12;

class LGCA {
  constructor(o){
    this.model = o.model;
    this.W = o.W; this.H = o.H;
    this.dens = o.dens;
    this.kernel = o.kernel || "sum";
    this.alpha = Number.isFinite(o.alpha) ? o.alpha : 0;
    if(this.alpha<0 || this.alpha>2) throw new RangeError("alpha must be in [0,2]");
    this.bosonField = o.bosonField || "site";
    this.bosonAlign = o.bosonAlign || "nematic";
    this.initialSeed=o.seed;
    this.rng = makeWell(o.seed);
    const N6 = this.W*this.H*6;
    /* C++ stores boson channel populations in int.  Int32Array preserves that
       range; Uint16Array would silently wrap a clustered channel at 65,536. */
    this.occ = new Int32Array(N6);
    this.src = new Int32Array(N6);
    this.tmp = new Int32Array(N6);
    this.prob = new Float64Array(6);
    this.m = new Float64Array(6);
    this.h = new Float64Array(6);
    this.cur = new Uint8Array(6);
    this.tri = new Uint8Array(6);
    this.perm = new Uint8Array(6);
    this.initLattice();
  }
  /* Mirrors the C++ realization init exactly, including its RNG draw order:
     one draw per (i,j,k) for the fractional fill, i outer, j inner, k innermost. */
  initLattice(){
    const {W,H,occ,rng} = this;
    const base = Math.trunc(this.dens);
    const frac = this.dens - base;
    let n = base*6*W*H;
    if(this.model==="boson") occ.fill(base); else occ.fill(0);
    for(let i=0;i<W;i++)
      for(let j=0;j<H;j++)
        for(let k=0;k<6;k++)
          if(rng()<frac){
            const b=(j*W+i)*6+k;
            if(this.model==="boson") occ[b]++; else occ[b]=1;
            n++;
          }
    this.nuparts = n;
    this.t = 0;
    this.src.fill(0);
  }
  neighbourField(i,j){
    const {W,H,occ,m} = this;
    m.fill(0);
    for(let d=0;d<6;d++){
      const X=(i+NB[d][0]+W)%W, Y=(j+NB[d][1]+H)%H, bb=(Y*W+X)*6;
      for(let k=0;k<6;k++) m[k]+=occ[bb+k];
    }
  }
  collideFermion(i,j,sens){
    const {occ,cur,tri,perm,m,h,rng} = this;
    const b=(j*this.W+i)*6;
    this.neighbourField(i,j);
    for(let a=0;a<6;a++){
      let s=0.0;
      for(let k=0;k<6;k++) s += JP[a][k]*m[k];
      h[a]=s;
    }
    if(this.kernel==="avg" || (this.kernel==="power" && this.alpha===1)){
      let M=0; for(let k=0;k<6;k++) M+=m[k];
      if(M!==0) for(let a=0;a<6;a++) h[a]/=M;
    } else if(this.kernel==="power" && this.alpha!==0){
      let M=0; for(let k=0;k<6;k++) M+=m[k];
      if(M!==0){ const divisor=Math.pow(M,this.alpha); for(let a=0;a<6;a++) h[a]/=divisor; }
    }
    for(let k=0;k<6;k++) cur[k]=occ[b+k];
    let s=0.0;
    for(let k=0;k<6;k++) if(cur[k]) s+=h[k];
    let E = sens*s;
    for(let it=0;it<METROPOLIS_STEPS;it++){
      for(let k=0;k<6;k++) perm[k]=k;
      for(let k=5;k>0;k--){
        const r = (rng()*(k+1))|0;
        const tp=perm[k]; perm[k]=perm[r]; perm[r]=tp;
      }
      tri.fill(0);
      for(let k=0;k<6;k++) if(cur[k]) tri[perm[k]]=1;
      let s2=0.0;
      for(let k=0;k<6;k++) if(tri[k]) s2+=h[k];
      const E2 = sens*s2;
      if(Math.log(rng()) < E2-E){ cur.set(tri); E=E2; }
    }
    const tb=b;
    for(let k=0;k<6;k++) this.tmp[tb+k]=cur[k];
  }
  collideBoson(i,j,sens,n){
    const {occ,prob,m,rng} = this;
    const b=(j*this.W+i)*6;
    const JT = this.bosonAlign==="nematic" ? JN : JP;
    let srcArr, numparts;
    if(this.bosonField==="site"){
      srcArr = null; numparts = n;
      for(let c=0;c<6;c++){
        let s=0.0;
        for(let k=0;k<6;k++) s += occ[b+k]*JT[k][c];
        prob[c] = sens*s;
      }
    } else {
      this.neighbourField(i,j);
      numparts = 0; for(let k=0;k<6;k++) numparts += m[k];
      for(let c=0;c<6;c++){
        let s=0.0;
        for(let k=0;k<6;k++) s += m[k]*JT[k][c];
        prob[c] = sens*s;
      }
    }
    let mx=prob[0];
    for(let c=1;c<6;c++) if(prob[c]>mx) mx=prob[c];
    let z=0.0;
    if(this.kernel==="avg" || (this.kernel==="power" && this.alpha===1)){
      for(let c=0;c<6;c++){ let v=prob[c]-mx; if(numparts!==0) v/=numparts; prob[c]=Math.exp(v); z+=prob[c]; }
    } else if(this.kernel==="power" && this.alpha!==0){
      const divisor=numparts===0?1:Math.pow(numparts,this.alpha);
      for(let c=0;c<6;c++){ prob[c]=Math.exp((prob[c]-mx)/divisor); z+=prob[c]; }
    } else {
      for(let c=0;c<6;c++){ prob[c]=Math.exp(prob[c]-mx); z+=prob[c]; }
    }
    for(let c=0;c<6;c++) prob[c]/=z;
    for(let c=1;c<6;c++) prob[c]+=prob[c-1];
    const tb=b;
    for(let p=0;p<n;p++){
      const r=rng(); let c=5;
      for(let q=0;q<6;q++) if(r<=prob[q]){ c=q; break; }
      this.tmp[tb+c]++;
    }
  }
  step(sens){
    const {W,H,occ,tmp} = this;
    tmp.fill(0);
    if(this.model==="boson"){
      for(let i=0;i<W;i++)
        for(let j=0;j<H;j++){
          const b=(j*W+i)*6;
          const n=occ[b]+occ[b+1]+occ[b+2]+occ[b+3]+occ[b+4]+occ[b+5];
          if(n!==0) this.collideBoson(i,j,sens,n);
        }
    } else {
      for(let i=0;i<W;i++)
        for(let j=0;j<H;j++){
          const b=(j*W+i)*6;
          const n=occ[b]+occ[b+1]+occ[b+2]+occ[b+3]+occ[b+4]+occ[b+5];
          if(n!==0) this.collideFermion(i,j,sens);
        }
    }
    this.src.set(tmp);
    for(let i=0;i<W;i++)
      for(let j=0;j<H;j++){
        const b=(j*W+i)*6;
        for(let k=0;k<6;k++){
          const X=(i+DST[k][0]+W)%W, Y=(j+DST[k][1]+H)%H;
          occ[(Y*W+X)*6+k]=tmp[b+k];
        }
      }
    this.t++;
  }
  measure(){
    const {W,H,occ,nuparts} = this;
    let x1=0,y1=0,x2=0,y2=0,S=0,occup=0;
    /* Keep i -> j -> channel accumulation order identical to C++.  The order
       is immaterial mathematically, but not to last-bit floating-point output. */
    for(let i=0;i<W;i++) for(let j=0;j<H;j++){
      const b=(j*W+i)*6; let nl=0;
      for(let k=0;k<6;k++){
        const n=occ[b+k]; if(!n) continue;
        nl+=n;
        x1+=n*CX[k]; y1+=n*CY[k];
        x2+=n*(2*CX[k]*CX[k]-1); y2+=n*(2*CX[k]*CY[k]);
      }
      if(nl===0) continue;
      occup++;
      const p=nl/nuparts; S+=p*Math.log(p);
    }
    return {
      polar:Math.sqrt(x1*x1+y1*y1)/nuparts,
      nematic:Math.sqrt(x2*x2+y2*y2)/nuparts,
      entropy:-S,
      spatial:1-(-S)/Math.log(W*H),
      occFrac:occup/(W*H),
      invOcc:occup>0?1/occup:0
    };
  }
  band(){
    const {W,H,occ,nuparts} = this;
    let pre=0;
    for(let i=0;i<W;i++) for(let j=0;j<H;j++){
      const b=(j*W+i)*6;
      let any=0; for(let k=0;k<6;k++) any+=occ[b+k];
      if(!any) continue;
      let x3=0,y3=0,nband=0,nsite=0;
      for(let k=1;k<5;k++){
        const idx=[[(i+k)%W,j],[(i+k)%W,(j-k+H)%H],[i,(j-k+H)%H],
                   [(i-k+W)%W,j],[(i-k+W)%W,(j+k)%H],[i,(j+k)%H]];
        for(let m=1;m<6;m++){
          for(let d=0;d<6;d++){
            const o=occ[(idx[d][1]*W+idx[d][0])*6+m];
            x3+=o*(2*CX[d]*CX[d]-1); y3+=o*(2*CX[d]*CY[d]); nband+=o;
          }
          nsite+=occ[b+m];
        }
      }
      if(nband>0) pre+=nsite*Math.sqrt(x3*x3+y3*y3)/(4*nband);
    }
    return pre/nuparts;
  }
}

/* Mirrors the C++ program end to end: one RNG stream across realizations. */
function runBatch(o){
  const opt = Object.assign({model:"boson",tsteps:100,iters:1,sens:2,dens:0.4,seed:12345,W:120,H:120,kernel:"sum",bosonField:"site",bosonAlign:"nematic"},o);
  const sim = new LGCA(opt);
  const acc = {p:0,p2:0,n:0,n2:0,e:0,e2:0,b:0,b2:0,i:0,i2:0};
  for(let r=0;r<opt.iters;r++){
    if(r>0) sim.initLattice();
    let m=null,bv=0;
    for(let t=0;t<opt.tsteps;t++) sim.step(opt.sens);
    m = sim.measure(); bv = sim.band();
    acc.p+=m.polar; acc.p2+=m.polar*m.polar;
    acc.n+=m.nematic; acc.n2+=m.nematic*m.nematic;
    acc.e+=m.entropy; acc.e2+=m.entropy*m.entropy;
    acc.b+=bv; acc.b2+=bv*bv;
    acc.i+=m.invOcc; acc.i2+=m.invOcc*m.invOcc;
  }
  const I=opt.iters;
  const se=(s,s2)=>{ const mean=s/I; return Math.sqrt(Math.abs(mean*mean-s2/I)/I); };
  const logA=Math.log(opt.W*opt.H);
  const line=[opt.sens*opt.dens, acc.p/I, se(acc.p,acc.p2), acc.n/I, se(acc.n,acc.n2),
    1-(acc.e/I)/logA, se(acc.e,acc.e2)/logA, acc.b/I, se(acc.b,acc.b2), acc.i/I, se(acc.i,acc.i2)].join(" ");
  console.log(line);
  return {line, polar:acc.p/I, nematic:acc.n/I};
}

function dumpState(o){
  const opt = Object.assign({model:"boson",tsteps:30,sens:2,dens:0.4,seed:424242,W:120,H:120,kernel:"sum",bosonField:"site",bosonAlign:"nematic"},o);
  const sim = new LGCA(opt);
  for(let t=0;t<opt.tsteps;t++) sim.step(opt.sens);
  const out=["DUMP t="+sim.t+" N="+sim.nuparts];
  for(let i=0;i<opt.W;i++)
    for(let j=0;j<opt.H;j++)
      for(let k=0;k<6;k++)
        out.push(String(sim.occ[(j*opt.W+i)*6+k]));
  return out.join("\n");
}
if(typeof window!=="undefined"){ window.runBatch=runBatch; window.dumpState=dumpState; }
