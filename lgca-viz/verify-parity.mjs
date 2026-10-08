#!/usr/bin/env node

// Compares lgca-lab.html with the C++: compiles a copy of lgca_clean-1.cpp with
// prints added (in a temp folder), runs both, and checks lattice, WELL1024a state,
// conservation, fermion exclusion and observables. The reference is the one committed
// at git HEAD (lgca_clean-1.cpp and rng/WELL1024a.c/.h), so an edit in the working tree
// cannot change what the Lab is compared with; outside a git checkout that has those
// files it is read from the working tree, with a note. Also checked: the empty lattice
// against the reference's NaN conventions, and boson sampling at forced draws exactly on
// the boundaries of the channel distribution.

import assert from "node:assert/strict";
import {execFileSync} from "node:child_process";
import {mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync} from "node:fs";
import {tmpdir} from "node:os";
import {dirname, join, resolve} from "node:path";
import {fileURLToPath} from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(here, "..");
const htmlPath = join(here, "lgca-lab.html");
const W = 120;
const H = 120;

const html = readFileSync(htmlPath, "utf8");
const engineMatch = html.match(/<script id="engineSrc">([\s\S]*?)<\/script>/);
assert(engineMatch, "could not find the v2 engine script");
const engineLogs=[];
const loadEngine = new Function("console",`${engineMatch[1]}\nreturn {LGCA, makeWell, runBatch, dumpState};`);
const {LGCA,makeWell,runBatch,dumpState} = loadEngine({log:value=>engineLogs.push(String(value))});

// The reference as committed: the three files at git HEAD when this is a git checkout
// whose HEAD has them, otherwise the working tree (with a note).
const REFERENCE_FILES = ["lgca/lgca_clean-1.cpp", "rng/WELL1024a.c", "rng/WELL1024a.h"];
function referenceSources() {
  const git = args => execFileSync("git", ["-C", root, ...args],
    {encoding:"utf8", stdio:["ignore", "pipe", "ignore"], maxBuffer:64 * 1024 * 1024});
  const working = file => { try { return readFileSync(join(root, file), "utf8"); } catch { return null; } };
  try {
    const head = git(["rev-parse", "--short", "HEAD"]).trim();
    const text = Object.fromEntries(REFERENCE_FILES.map(file => [file, git(["show", `HEAD:./${file}`])]));
    const differs = REFERENCE_FILES.filter(file => working(file) !== text[file]);
    return {text, source:`git HEAD ${head}` +
      (differs.length ? ` (the working tree's ${differs.join(", ")} differs from it; HEAD is used)` : "")};
  } catch {
    return {text:Object.fromEntries(REFERENCE_FILES.map(file => [file, readFileSync(join(root, file), "utf8")])),
            source:"the working tree (note: not a git checkout whose HEAD has these files, so nothing pins them)"};
  }
}
const reference = referenceSources();

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
  // The engine control was removed from the page (2026-10-01): WebAssembly is simply the
  // engine, the JavaScript engine (the one this suite tests) is the silent fallback when it
  // cannot load, and no choice is stored. ?engine= in the URL and labEngine() on the
  // console stay as developer paths.
  assert(!/id="segEngine"|id="engineNote"/.test(html),"the page must have no engine control");
  assert.match(html,/const ENGINE_DEFAULT="wasm"/,"the default engine must be WebAssembly");
  assert(!/localStorage\.(getItem|setItem)\([^)]*(lgcaEngine|ENGINE_KEY)/.test(html)&&!/ENGINE_KEY/.test(html),
    "nothing may read or write a stored engine choice (localStorage lgcaEngine)");
  assert.match(html,/localStorage\.removeItem\("lgcaEngine"\)/,"an old stored engine choice must be removed on load");
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
  let cpp = reference.text["lgca/lgca_clean-1.cpp"];
  if (model === "fermion") {
    cpp = cpp.replace(/^#define MODEL\s+BOSON.*$/m, "#define MODEL   FERMION");
  }
  cpp = cpp.replace(
    "    std::cout << sens * dens << \" \"",
    "    std::cout << std::setprecision(17);\n    std::cout << sens * dens << \" \"",
  );
  // Our one deliberate change to the reference (docs/FIDELITY.md, "The band correction"):
  // the band counts channel 0 too, so the reference copy compiled here gets the same line.
  const BAND_FIX = ["for (int m = 1; m < NODES; ++m) {         // channels (skips 0)",
                    "for (int m = 0; m < NODES; ++m) {         // channels (all six: band fix)"];
  assert.equal(cpp.split(BAND_FIX[0]).length, 2, "the reference's band loop is not there once");
  cpp = cpp.replace(...BAND_FIX);
  cpp = cpp.replace("\n    return 0;\n}", `${probe}\n    return 0;\n}`);
  assert(cpp.includes("__PARITY_STATE__"), "failed to instrument C++ reference");

  const combined = `${reference.text["rng/WELL1024a.c"]}\n${macroCleanup}\n#include <iomanip>\n${cpp}`;
  // laid out as in the repository: the reference includes "../rng/WELL1024a.h"
  mkdirSync(join(tempDir, "lgca"), {recursive:true});
  mkdirSync(join(tempDir, "rng"), {recursive:true});
  writeFileSync(join(tempDir, "rng", "WELL1024a.h"), reference.text["rng/WELL1024a.h"]);
  const sourcePath = join(tempDir, "lgca", `reference-${model}.cpp`);
  const binaryPath = join(tempDir, `reference-${model}`);
  writeFileSync(sourcePath, combined);
  // -ffp-contract=off: clang would otherwise fuse a*b+c into one FMA instruction (one
  // rounding instead of two); JavaScript never fuses, so the reference must not either.
  execFileSync(process.env.CXX || "c++", [
    "-std=c++17", "-O2", "-ffp-contract=off", sourcePath, "-o", binaryPath,
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

// The empty lattice (density 0), which the reference runs: its polar, nematic and band
// divide by N = 0 and print NaN, and so do their errors; spatial order is 1 with error 0,
// 1/occupied and its error are 0. The Lab keeps those conventions (docs/PHYSICS.md, "The
// printed output and empty lattices"): NaN exactly where the reference prints it.
const emptyScenarios = [
  {model:"boson",   tsteps:3, iters:2, sens:2,   dens:0, seed:12345},
  {model:"fermion", tsteps:3, iters:2, sens:0.8, dens:0, seed:12345},
];
const EMPTY_NAN = [1, 2, 3, 4, 7, 8];       // polar, error, nematic, error, band, error
function compareEmptyMetrics(actual, expected, label) {
  assert.equal(actual.length, expected.length, `${label}: metric count differs`);
  for (let q = 0; q < actual.length; q++) {
    if (EMPTY_NAN.includes(q)) {
      assert(Number.isNaN(expected[q]), `${label}: the reference prints ${expected[q]} for metric ${q}, not NaN`);
      assert(Number.isNaN(actual[q]), `${label}: metric ${q} is ${actual[q]}; the reference gives NaN for an empty lattice`);
    } else {
      assert.equal(actual[q], expected[q], `${label}: metric ${q} differs: JS=${actual[q]} C++=${expected[q]}`);
    }
  }
  assert.deepEqual(expected.filter((_, q) => !EMPTY_NAN.includes(q)), [0, 1, 0, 0, 0],
    `${label}: the reference's finite values (sens*dens, spatial, its error, 1/occupied, its error)`);
}

// Boson sampling on the boundaries of the channel distribution. The reference takes the
// first channel c with r <= cdf[c] (sample_channel), so a draw exactly on a boundary, or on
// a boundary repeated by zero-weight channels, belongs to the lower channel. One particle
// at one site with sensitivity 1000: every weight but the largest underflows to 0, so the
// distribution is exact. The draw is forced through the generator state: with every word
// 0 but S[i] = x, the next output is x ^ (x << 7) (mod 2^32), and
// x = w ^ w<<7 ^ w<<14 ^ w<<21 ^ w<<28 makes it the word w.
function stateForDraw(word) {
  const x = (word ^ (word << 7) ^ (word << 14) ^ (word << 21) ^ (word << 28)) >>> 0;
  const state = new Uint32Array(33);
  state[5] = x; state[32] = 5;
  return state;
}
const HALF = 0x80000000;
const CDF_CASES = [
  // occupied channel, the exact distribution, then [forced word, channel the rule picks]
  {channel:0, cdf:[0.5, 0.5, 0.5, 1, 1, 1], draws:[[HALF, 0], [HALF + 1, 3], [0, 0], [0xffffffff, 3]]},
  {channel:1, cdf:[0, 0.5, 0.5, 0.5, 1, 1], draws:[[0, 0], [1, 1], [HALF, 1], [HALF + 1, 4]]},
  {channel:3, cdf:[0.5, 0.5, 0.5, 1, 1, 1], draws:[[HALF, 0], [HALF + 1, 3]]},
];
const ruleChannel = (r, cdf) => { for (let c = 0; c < 6; c++) if (r <= cdf[c]) return c; return 5; };
// One collision of the lone particle with the forced draw, on any engine with LGCA's
// interface; returns the channel it took, the distribution, and the generator state after.
function forcedCollision(make, channel, word) {
  const sim = make({model:"boson", W:6, H:5, dens:0, seed:99, kernel:"sum", bosonField:"site", bosonAlign:"nematic"});
  try {
    const site = (3 * 6 + 2) * 6;
    sim.occ.fill(0); sim.occ[site + channel] = 1; sim.nuparts = 1;
    sim.rng.setState(stateForDraw(word));
    sim.step(1000);
    const src = sim.src;
    const taken = [0, 1, 2, 3, 4, 5].filter(c => src[site + c] === 1);
    let total = 0; for (const v of src) total += v;
    assert(taken.length === 1 && total === 1, "the forced collision did not move exactly one particle");
    return {channel:taken[0], prob:sim.prob ? Array.from(sim.prob) : null, rng:sim.rng.getState()};
  } finally { if (sim.destroy) sim.destroy(); }
}
function verifyCdfBoundaries(make, engine) {
  let checks = 0;
  for (const {channel, cdf, draws} of CDF_CASES)
    for (const [word, expected] of draws) {
      const well = makeWell(0); well.setState(stateForDraw(word));
      const r = well();
      assert.equal(r, word * 2.32830643653869628906e-10, "the forced generator state does not give the wanted draw");
      assert.equal(ruleChannel(r, cdf), expected, "the test's own expectation disagrees with the reference rule");
      const got = forcedCollision(make, channel, word);
      if (got.prob) assert.deepEqual(got.prob, cdf, `${engine}: channel ${channel}: the distribution differs from ${cdf}`);
      assert.equal(got.channel, expected,
        `${engine}: particle in channel ${channel}, draw ${r} (cdf ${cdf.join(" ")}): took channel ${got.channel}, the reference rule gives ${expected}`);
      assert.deepEqual(got.rng, well.getState(), `${engine}: the collision did not consume exactly the one forced draw`);
      checks++;
    }
  return checks;
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
const cdfChecks = verifyCdfBoundaries(o => new LGCA(o), "JavaScript");
console.log(`PASS boson sampling at ${cdfChecks} forced draws on channel-distribution boundaries ` +
  "(ties and zero-weight channels) picks the reference's channel");
console.log(`Reference C++: ${reference.source}`);

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
  for (const scenario of emptyScenarios) {
    const label=`empty lattice: ${scenario.model} t=${scenario.tsteps} r=${scenario.iters} s=${scenario.sens} d=0 seed=${scenario.seed}`;
    const js=runJavaScript(scenario);
    const cpp=runCpp(binaries[scenario.model],scenario);
    compareState(js.sim,cpp.state,label);
    compareRng(js.sim,cpp.rng,label);
    compareEmptyMetrics(js.metrics,cpp.metrics,label);
    console.log(`PASS ${label}: lattice, generator and the reference's NaN conventions`);
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
