#!/usr/bin/env node

/*
 * State-for-state regression test for lgca-lab.html.
 *
 * The test compiles an instrumented copy of the actual lgca_clean-1.cpp source
 * in a temporary directory, then compares the C++ and JavaScript lattice,
 * WELL1024a state, particle conservation, exclusion invariant, and observables.
 * No generated source or binary is left in the repository.
 */

import assert from "node:assert/strict";
import {execFileSync} from "node:child_process";
import {mkdtempSync, readFileSync, rmSync, writeFileSync} from "node:fs";
import {tmpdir} from "node:os";
import {dirname, join, resolve} from "node:path";
import {fileURLToPath} from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(here, "..");
const htmlPath = join(here, "lgca-lab.html");
const cppPath = join(root, "lgca", "lgca_clean-1.cpp");
const rngPath = join(root, "rng", "WELL1024a.c");
const W = 120;
const H = 120;

const html = readFileSync(htmlPath, "utf8");
const engineMatch = html.match(/<script id="engineSrc">([\s\S]*?)<\/script>/);
assert(engineMatch, "could not find the v2 engine script");
const engineLogs=[];
const loadEngine = new Function("console",`${engineMatch[1]}\nreturn {LGCA, makeWell, runBatch, dumpState};`);
const {LGCA,runBatch,dumpState} = loadEngine({log:value=>engineLogs.push(String(value))});

function verifyHtmlStructureAndDefaults() {
  const scripts=[...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)].map(m=>m[1]);
  for (const script of scripts) new Function(script);

  const ids=[...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]);
  assert.equal(new Set(ids).size,ids.length,"HTML contains duplicate ids");
  const refs=[...scripts.join("\n").matchAll(/\$\("([^"]+)"\)/g)].map(m=>m[1]);
  for (const id of refs) assert(ids.includes(id),`UI script references missing id ${id}`);

  const selectedValue=id=>{
    const group=html.match(new RegExp(`<span[^>]*id="${id}"[^>]*>([\\s\\S]*?)<\\/span>`));
    assert(group,`missing segmented control ${id}`);
    const active=group[1].match(/<button[^>]*data-v="([^"]+)"[^>]*class="on"/);
    assert(active,`segmented control ${id} has no active default`);
    return active[1];
  };
  assert.equal(selectedValue("segModel"),"boson");
  assert.equal(selectedValue("segKernel"),"sum");
  assert.equal(selectedValue("segField"),"site");
  assert.equal(selectedValue("segAlign"),"nematic");
  assert.match(html,/<option selected>120<\/option>/,"C++ lattice size 120 must be the default");
  assert.match(html,/id="inSens"[^>]*value="2"/,"C++ sensitivity preset changed unexpectedly");
  assert.match(html,/id="inDens"[^>]*value="0\.4"/,"C++ density preset changed unexpectedly");
}

function verifyReplayRoundTrip() {
  const sens=2;
  const sim=new LGCA({model:"boson",W:60,H:60,dens:0.4,seed:12345,
    kernel:"sum",bosonField:"site",bosonAlign:"nematic"});
  const snapshots=[];
  for (let t=1;t<=8;t++) {
    sim.step(sens);
    snapshots[t]={occ:new Int32Array(sim.occ),src:new Int32Array(sim.src),rng:sim.rng.getState()};
  }

  const restored=snapshots[4];
  sim.occ.fill(0);
  for(let i=0;i<sim.W;i++) for(let j=0;j<sim.H;j++) for(let k=0;k<6;k++){
    const b=(j*sim.W+i)*6;
    const X=(i+DST_FOR_REPLAY[k][0]+sim.W)%sim.W;
    const Y=(j+DST_FOR_REPLAY[k][1]+sim.H)%sim.H;
    sim.occ[(Y*sim.W+X)*6+k]=restored.src[b+k];
  }
  assert.deepEqual(sim.occ,restored.occ,"replay source state does not reconstruct streamed state");
  sim.src.set(restored.src);
  sim.rng.setState(restored.rng);
  sim.t=4;
  sim.step(sens);
  assert.deepEqual(sim.src,snapshots[5].src,"replay did not reproduce the next collision state");
  assert.deepEqual(sim.occ,snapshots[5].occ,"replay did not reproduce the next streamed state");
  assert.deepEqual(sim.rng.getState(),snapshots[5].rng,"replay did not reproduce the RNG state");
}

const DST_FOR_REPLAY=[[1,0],[1,-1],[0,-1],[-1,0],[-1,1],[0,1]];

const scenarios = [
  {model:"boson",   tsteps:1, iters:1, sens:0,   dens:0.4,  seed:0},
  {model:"boson",   tsteps:7, iters:1, sens:2,   dens:0.4,  seed:12345},
  {model:"boson",   tsteps:3, iters:1, sens:8,   dens:1,    seed:987654321},
  {model:"boson",   tsteps:17,iters:1, sens:4.25,dens:0.05, seed:42},
  {model:"boson",   tsteps:5, iters:2, sens:7.5, dens:2.25, seed:4294967295},
  {model:"fermion", tsteps:1, iters:1, sens:0,   dens:0.2,  seed:0},
  {model:"fermion", tsteps:6, iters:1, sens:0.8, dens:0.2,  seed:12345},
  {model:"fermion", tsteps:11,iters:1, sens:3.25,dens:0.05, seed:42},
  {model:"fermion", tsteps:7, iters:1, sens:8,   dens:0.6,  seed:987654321},
  {model:"fermion", tsteps:4, iters:2, sens:8,   dens:0.95, seed:4294967295},
];

const macroCleanup = `
#undef W
#undef R
#undef M1
#undef M2
#undef M3
#undef MAT0POS
#undef MAT0NEG
#undef Identity
#undef V0
#undef VM1
#undef VM2
#undef VM3
#undef VRm1
#undef newV0
#undef newV1
#undef FACT
`;

const probe = `
    std::cout << "\\n__PARITY_RNG__ " << state_i;
    for (int q = 0; q < 32; ++q) std::cout << " " << STATE[q];
    std::cout << "\\n__PARITY_STATE__\\n";
    for (int i = 0; i < XDIM; ++i)
        for (int j = 0; j < YDIM; ++j)
            for (int k = 0; k < NODES; ++k)
                std::cout << static_cast<unsigned long long>(sim.lattice[i][j][k]) << " ";
    std::cout << "\\n";
`;

function buildReference(model, tempDir) {
  let cpp = readFileSync(cppPath, "utf8");
  if (model === "fermion") {
    cpp = cpp.replace(/^#define MODEL\s+BOSON.*$/m, "#define MODEL   FERMION");
  }
  cpp = cpp.replace(
    "    std::cout << sens * dens << \" \"",
    "    std::cout << std::setprecision(17);\n    std::cout << sens * dens << \" \"",
  );
  cpp = cpp.replace("\n    return 0;\n}", `${probe}\n    return 0;\n}`);
  assert(cpp.includes("__PARITY_STATE__"), "failed to instrument C++ reference");

  const combined = `${readFileSync(rngPath, "utf8")}\n${macroCleanup}\n#include <iomanip>\n${cpp}`;
  const sourcePath = join(tempDir, `reference-${model}.cpp`);
  const binaryPath = join(tempDir, `reference-${model}`);
  writeFileSync(sourcePath, combined);
  execFileSync(process.env.CXX || "c++", [
    "-std=c++17", "-O2", "-I", join(root, "lgca"), sourcePath, "-o", binaryPath,
  ], {stdio:"inherit"});
  return binaryPath;
}

function verifyInvariants(sim, model, label) {
  let total = 0;
  for (const value of sim.occ) {
    assert(value >= 0, `${label}: negative channel population`);
    if (model === "fermion") {
      assert(value === 0 || value === 1, `${label}: fermion exclusion violated`);
    }
    total += value;
  }
  assert.equal(total, sim.nuparts, `${label}: particle number not conserved`);
}

function runJavaScript(o) {
  const sim = new LGCA({
    model:o.model, W, H, dens:o.dens, seed:o.seed,
    kernel:"sum", bosonField:"site", bosonAlign:"nematic",
  });
  assert(sim.occ instanceof Int32Array, "channel populations must have C++ int range");

  const acc = {p:0,p2:0,n:0,n2:0,e:0,e2:0,b:0,b2:0,i:0,i2:0};
  for (let r=0; r<o.iters; r++) {
    if (r>0) sim.initLattice();
    verifyInvariants(sim, o.model, `${o.model} realization ${r} init`);
    for (let t=0; t<o.tsteps; t++) {
      sim.step(o.sens);
      verifyInvariants(sim, o.model, `${o.model} realization ${r} step ${t+1}`);
    }
    const m=sim.measure(), b=sim.band();
    acc.p+=m.polar;   acc.p2+=m.polar*m.polar;
    acc.n+=m.nematic; acc.n2+=m.nematic*m.nematic;
    acc.e+=m.entropy; acc.e2+=m.entropy*m.entropy;
    acc.b+=b;         acc.b2+=b*b;
    acc.i+=m.invOcc;  acc.i2+=m.invOcc*m.invOcc;
  }
  const I=o.iters;
  const se=(sum,sum2)=>Math.sqrt(Math.abs((sum/I)*(sum/I)-sum2/I)/I);
  const logA=Math.log(W*H);
  const metrics=[o.sens*o.dens,
    acc.p/I, se(acc.p,acc.p2), acc.n/I, se(acc.n,acc.n2),
    1-(acc.e/I)/logA, se(acc.e,acc.e2)/logA,
    acc.b/I, se(acc.b,acc.b2), acc.i/I, se(acc.i,acc.i2)];
  return {sim, metrics};
}

function runCpp(binary, o) {
  const stdout = execFileSync(binary,
    [o.tsteps, o.iters, o.sens, o.dens, o.seed].map(String),
    {encoding:"utf8", maxBuffer:64*1024*1024});
  const [beforeRng, afterRng] = stdout.split("__PARITY_RNG__ ");
  const [rngText, stateText] = afterRng.split("\n__PARITY_STATE__\n");
  assert(beforeRng && rngText && stateText, "malformed C++ probe output");
  const metricLine = beforeRng.trim().split("\n").at(-1);
  const metrics = metricLine.trim().split(/\s+/).map(Number);
  const rng = rngText.trim().split(/\s+/).map(Number);
  const state = stateText.trim().split(/\s+/).map(Number);
  return {metrics, rng, state};
}

function compareMetrics(actual, expected, label) {
  assert.equal(actual.length, expected.length, `${label}: metric count differs`);
  for (let q=0; q<actual.length; q++) {
    const scale=Math.max(1,Math.abs(expected[q]));
    const error=Math.abs(actual[q]-expected[q]);
    assert(error <= 5e-14*scale,
      `${label}: metric ${q} differs: JS=${actual[q]} C++=${expected[q]} error=${error}`);
  }
}

function compareRng(sim, cppRng, label) {
  const jsRng=sim.rng.getState();
  assert.equal(cppRng.length, 33, `${label}: malformed C++ RNG state`);
  assert.equal(jsRng[32], cppRng[0], `${label}: WELL state index differs`);
  for (let q=0; q<32; q++) {
    assert.equal(jsRng[q], cppRng[q+1], `${label}: WELL state word ${q} differs`);
  }
}

function compareState(sim, cppState, label) {
  assert.equal(cppState.length, W*H*6, `${label}: lattice size differs`);
  let q=0;
  for (let i=0; i<W; i++) for (let j=0; j<H; j++) for (let k=0; k<6; k++,q++) {
    const js=sim.occ[(j*W+i)*6+k];
    assert.equal(js, cppState[q],
      `${label}: first lattice difference at site (${i},${j}), channel ${k}: JS=${js}, C++=${cppState[q]}`);
  }
}

verifyHtmlStructureAndDefaults();
verifyReplayRoundTrip();
console.log("PASS HTML structure, C++ defaults, and exact replay round trip");

const tempDir=mkdtempSync(join(tmpdir(), "lgca-parity-"));
try {
  const binaries={
    boson:buildReference("boson",tempDir),
    fermion:buildReference("fermion",tempDir),
  };
  for (const scenario of scenarios) {
    const label=`${scenario.model} t=${scenario.tsteps} r=${scenario.iters} s=${scenario.sens} d=${scenario.dens} seed=${scenario.seed}`;
    const js=runJavaScript(scenario);
    const cpp=runCpp(binaries[scenario.model],scenario);
    compareState(js.sim,cpp.state,label);
    compareRng(js.sim,cpp.rng,label);
    compareMetrics(js.metrics,cpp.metrics,label);
    console.log(`PASS ${label}`);
  }

  const apiScenario=scenarios.find(s=>s.model==="boson"&&s.iters===2);
  engineLogs.length=0;
  const batch=runBatch({...apiScenario,W,H,kernel:"sum",bosonField:"site",bosonAlign:"nematic"});
  const batchCpp=runCpp(binaries.boson,apiScenario);
  compareMetrics(batch.line.split(/\s+/).map(Number),batchCpp.metrics,"runBatch console API");
  assert.deepEqual(engineLogs,[batch.line],"runBatch did not log exactly its returned output line");

  const dumpScenario=scenarios.find(s=>s.model==="fermion"&&s.tsteps===6);
  const dump=dumpState({...dumpScenario,W,H,kernel:"sum",bosonField:"site",bosonAlign:"nematic"});
  const dumpLines=dump.split("\n");
  const dumpInit=new LGCA({...dumpScenario,W,H,kernel:"sum",bosonField:"site",bosonAlign:"nematic"});
  assert.equal(dumpLines.shift(),`DUMP t=${dumpScenario.tsteps} N=${dumpInit.nuparts}`,
    "dumpState header differs");
  const dumpCpp=runCpp(binaries.fermion,dumpScenario);
  assert.deepEqual(dumpLines.map(Number),dumpCpp.state,"dumpState console API differs from C++");
  console.log("PASS runBatch and dumpState console APIs");
  console.log(`\nAll ${scenarios.length} C++/JavaScript parity scenarios passed.`);
} finally {
  rmSync(tempDir,{recursive:true,force:true});
}
