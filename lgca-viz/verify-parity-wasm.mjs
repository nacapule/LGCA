#!/usr/bin/env node

// Compares the Lab's WebAssembly engine with the C++ reference. The sibling of
// verify-parity.mjs, which does the same for the Lab's JavaScript engine; run both after
// any change to lgca-lab.html:
//
//   node verify-parity.mjs && node verify-parity-wasm.mjs
//
// The engine under test is the one the page carries: the text of its
// <script id="wasmEngine" type="text/plain"> block (written there by
// native/inline-wasm.mjs), loaded and driven by the page's own code (loadWasmEngine,
// class WasmSim and simWorkerMain in the workerSrc block). Five parts:
//   1. the block: before the page's script, its text exactly the file its comment names
//      (size and sha256), nothing in it the HTML parser could misread. Only a note if
//      native/build/lgca_wasm.js is missing or differs (then the block is out of date,
//      or the build is newer than the page);
//   2. the 10 scenarios of verify-parity.mjs: lattice, WELL1024a state and output line
//      against the reference C++ (compiled the same way, from the copy committed at git
//      HEAD), and every step against the JavaScript engine (lattices, generator state, N
//      and t identical); and the empty lattice against the reference's NaN conventions;
//   3. a replay round trip as the worker does it: restore a saved state, step on; the
//      research options (every kernel, boson field and alignment, also changed between
//      steps), step for step against the JavaScript engine; boson sampling at draws forced
//      onto the boundaries of the channel distribution, on both engines; a malformed
//      generator state refused by both engines, which stay as they were; and the density
//      diagnostics of the worker against numbers worked out by hand;
//   4. the worker protocol: a worker running each engine receives a session of requests
//      (steps, option changes, rewinds, looks at buffered ticks) and every reply must be
//      the one worked out without the worker: fresh JavaScript engines stepped from the
//      seed, a rewind or a look at tick k answered by running the seed again to k (so the
//      generator words after a rewind are checked against an independent run), also
//      when requests arrive while the module loads. Further sessions run past the replay
//      buffer's capacity (eviction, a rewind to the oldest tick, a branch after the ring
//      has wrapped), across step 2^31, with batches cut short by their time budget, and
//      with the engine switched while the module loads. Deliberately broken copies of the
//      worker (the generator not restored on rewind; a buffered tick's look given the
//      current generator; the ring's head advanced twice once full) must fail that check.
//      A worker falls back to the JavaScript engine, and says so, when the module cannot
//      load or refuses the settings; after a failed load, an init with retryWasm tries
//      again (a module that fails once, then loads), while an init without it does not; a
//      module that never loads is tried once per such request, and the requests after it
//      are still answered in order;
//   5. the page: its three scripts run here against a small stand-in for the document,
//      with the main-thread host and the real module: the engine switched with
//      labEngine("js") while the module loads, a rewind whose target leaves the buffer,
//      the CSV and the history kept for it through repeated thinning, the study link's
//      defaults, and the smaller controls and displays.
// Entropy (and the spatial order computed from it) may differ by one unit in the last
// place between the engines: the WebAssembly build has its own log(), as the C++ has
// (docs/FIDELITY.md). Those two are compared to 1e-12 relative; everything else exactly.

import assert from "node:assert/strict";
import {execFileSync} from "node:child_process";
import {resolveObjectURL} from "node:buffer";
import {createHash} from "node:crypto";
import {existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync} from "node:fs";
import {createRequire} from "node:module";
import {tmpdir} from "node:os";
import {dirname, join, resolve} from "node:path";
import {fileURLToPath} from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(here, "..");
const htmlPath = join(here, "lgca-lab.html");
const gluePath = join(root, "native", "build", "lgca_wasm.js");
const W = 120;
const H = 120;

const html = readFileSync(htmlPath, "utf8");
const scriptBlock = id => {
  const m = html.match(new RegExp(`<script id="${id}">([\\s\\S]*?)<\\/script>`));
  assert(m, `could not find the ${id} script`);
  return m[1];
};

// The page's engine and worker code, joined as the worker joins them (workerSource()).
// `workerText` is the worker block, or a deliberately broken copy of it (part 4).
const loadLab = workerText => new Function("console",
  `${scriptBlock("engineSrc")}\n${workerText}\n` +
  "return {LGCA, WasmSim, loadWasmEngine, simWorkerMain, streamInto, densityStats};")({log() {}});
const {LGCA, WasmSim, loadWasmEngine, simWorkerMain, streamInto, densityStats} = loadLab(scriptBlock("workerSrc"));

// The loader text runs its Node branch here: it looks for these CommonJS names, which a
// module (this file) does not have. In a browser it takes its web branch instead.
globalThis.require = createRequire(import.meta.url);
globalThis.__dirname = here;
globalThis.__filename = htmlPath;

const sha256 = text => createHash("sha256").update(text).digest("hex");
const relativeClose = (a, b) => a === b || Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(b));

// ------------------------------------------------------------------------------------
//  1. The wasmEngine block
// ------------------------------------------------------------------------------------

const BLOCK = new RegExp(
  "<!-- wasmEngine: generated by native/inline-wasm\\.mjs from (\\S+)\\n" +
  "\\s+\\((\\d+) bytes, sha256 ([0-9a-f]{64})\\)\\.\\n[\\s\\S]*?-->\\n" +
  '<script id="wasmEngine" type="text/plain">\\n([\\s\\S]*?)\\n</script>');

function verifyBlock() {
  const found = html.match(BLOCK);
  assert(found, "no wasmEngine block with its generated comment (run node native/inline-wasm.mjs)");
  const [whole, source, bytes, hash, glue] = found;
  assert.equal(html.split('id="wasmEngine"').length, 2, "more than one wasmEngine block");

  // The text is exactly the file the comment names: the file minus its last newline.
  const fileText = `${glue}\n`;
  assert.equal(Buffer.byteLength(fileText), Number(bytes), "the block's size differs from its comment");
  assert.equal(sha256(fileText), hash, "the block's sha256 differs from its comment (edited by hand?)");

  // The HTML parser must hand the text over unchanged (see native/inline-wasm.mjs).
  assert(!/<\/script|<!--|<script/i.test(glue), 'the block holds "</script", "<!--" or "<script"');
  assert(/^[\x09\x0a\x20-\x7e]*$/.test(glue), "the block holds characters other than printable ASCII");
  assert(glue.includes("createLgcaWasm"), "the block does not define createLgcaWasm");

  // The page's script reads the block when it starts the simulation, so the block must
  // come before that script, and after the worker block (where inline-wasm.mjs puts it).
  const at = found.index;
  assert(at > html.indexOf('<script id="workerSrc">'), "the block is not after the worker block");
  const pageScript = html.indexOf("<script>");
  assert(pageScript > at + whole.length, "the block is not before the page's script");

  let note = "";
  if (!existsSync(gluePath)) note = `note: ${source} is not there (not built here); the block was not compared with it`;
  else if (readFileSync(gluePath, "utf8") !== fileText)
    note = `note: ${source} differs from the block: the page holds an older build, or the build is newer ` +
           "than the page (node native/inline-wasm.mjs copies it in)";
  return {glue, source, hash, note};
}

// ------------------------------------------------------------------------------------
//  2. The 10 scenarios against the C++ (as in verify-parity.mjs)
// ------------------------------------------------------------------------------------

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

// The reference C++, with prints of the final generator state and lattice added: the
// same instrumentation and compiler flags as verify-parity.mjs, and the same sources: the
// three files committed at git HEAD when this is a git checkout whose HEAD has them (an
// edit in the working tree cannot change the oracle), otherwise the working tree.
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
const macroCleanup = ["W", "R", "M1", "M2", "M3", "MAT0POS", "MAT0NEG", "Identity", "V0", "VM1",
  "VM2", "VM3", "VRm1", "newV0", "newV1", "FACT"].map(name => `#undef ${name}`).join("\n");

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
  if (model === "fermion") cpp = cpp.replace(/^#define MODEL\s+BOSON.*$/m, "#define MODEL   FERMION");
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
  assert(cpp.includes("__PARITY_STATE__"), "failed to instrument the C++ reference");

  // laid out as in the repository: the reference includes "../rng/WELL1024a.h"
  mkdirSync(join(tempDir, "lgca"), {recursive:true});
  mkdirSync(join(tempDir, "rng"), {recursive:true});
  writeFileSync(join(tempDir, "rng", "WELL1024a.h"), reference.text["rng/WELL1024a.h"]);
  const sourcePath = join(tempDir, "lgca", `reference-${model}.cpp`);
  const binaryPath = join(tempDir, `reference-${model}`);
  writeFileSync(sourcePath, `${reference.text["rng/WELL1024a.c"]}\n${macroCleanup}\n#include <iomanip>\n${cpp}`);
  // -ffp-contract=off: no fused multiply-adds, which neither JavaScript nor this
  // WebAssembly build uses.
  execFileSync(process.env.CXX || "c++", [
    "-std=c++17", "-O2", "-ffp-contract=off", sourcePath, "-o", binaryPath,
  ], {stdio:"inherit"});
  return binaryPath;
}

function runCpp(binary, o) {
  const stdout = execFileSync(binary, [o.tsteps, o.iters, o.sens, o.dens, o.seed].map(String),
    {encoding:"utf8", maxBuffer:64 * 1024 * 1024});
  const [beforeRng, afterRng] = stdout.split("__PARITY_RNG__ ");
  const [rngText, stateText] = afterRng.split("\n__PARITY_STATE__\n");
  assert(beforeRng && rngText && stateText, "malformed C++ probe output");
  return {
    metrics: beforeRng.trim().split("\n").at(-1).trim().split(/\s+/).map(Number),
    rng: rngText.trim().split(/\s+/).map(Number),
    state: stateText.trim().split(/\s+/).map(Number),
  };
}

function verifyInvariants(sim, model, label) {
  let total = 0;
  for (const value of sim.occ) {
    assert(value >= 0, `${label}: negative channel population`);
    if (model === "fermion") assert(value === 0 || value === 1, `${label}: fermion exclusion violated`);
    total += value;
  }
  assert.equal(total, sim.nuparts, `${label}: particle number not conserved`);
}

// The two engines after the same steps: lattices, generator state, N and t identical.
function compareEngines(wasm, js, label) {
  const first = (a, b) => { for (let q = 0; q < a.length; q++) if (a[q] !== b[q]) return q; return -1; };
  for (const name of ["occ", "src"]) {
    const a = wasm[name], b = js[name];
    assert.equal(a.length, b.length, `${label}: ${name} sizes differ`);
    const q = first(a, b);
    assert.equal(q, -1, `${label}: ${name} differs from the JavaScript engine at index ${q}: ` +
      `WebAssembly ${a[q]}, JavaScript ${b[q]}`);
  }
  assert.deepEqual(wasm.rng.getState(), js.rng.getState(), `${label}: generator state differs from the JavaScript engine`);
  assert.equal(wasm.nuparts, js.nuparts, `${label}: N differs from the JavaScript engine`);
  assert.equal(wasm.t, js.t, `${label}: t differs from the JavaScript engine`);
}

// One scenario as the reference runs it, on both engines in step; returns the
// WebAssembly engine and its output line.
function runScenario(engine, o, label) {
  const options = {model:o.model, W, H, dens:o.dens, seed:o.seed, kernel:"sum", bosonField:"site", bosonAlign:"nematic"};
  const sim = new WasmSim(engine, options);
  const js = new LGCA(options);
  const acc = {p:0, p2:0, n:0, n2:0, e:0, e2:0, b:0, b2:0, i:0, i2:0};
  for (let r = 0; r < o.iters; r++) {
    if (r > 0) { sim.initLattice(); js.initLattice(); }
    verifyInvariants(sim, o.model, `${label} realization ${r} init`);
    compareEngines(sim, js, `${label} realization ${r} init`);
    for (let t = 0; t < o.tsteps; t++) {
      sim.step(o.sens); js.step(o.sens);
      verifyInvariants(sim, o.model, `${label} realization ${r} step ${t + 1}`);
      compareEngines(sim, js, `${label} realization ${r} step ${t + 1}`);
    }
    const m = sim.measure(), b = sim.band();
    const mj = js.measure();
    for (const key of ["polar", "nematic", "occFrac", "invOcc"])
      assert.equal(m[key], mj[key], `${label}: ${key} differs from the JavaScript engine`);
    for (const key of ["entropy", "spatial"])
      assert(relativeClose(m[key], mj[key]), `${label}: ${key} differs from the JavaScript engine: ${m[key]} vs ${mj[key]}`);
    assert.equal(b, js.band(), `${label}: band differs from the JavaScript engine`);
    acc.p += m.polar;   acc.p2 += m.polar * m.polar;
    acc.n += m.nematic; acc.n2 += m.nematic * m.nematic;
    acc.e += m.entropy; acc.e2 += m.entropy * m.entropy;
    acc.b += b;         acc.b2 += b * b;
    acc.i += m.invOcc;  acc.i2 += m.invOcc * m.invOcc;
  }
  // The reference's output line, with verify-parity.mjs's expressions.
  const I = o.iters;
  const se = (sum, sum2) => Math.sqrt(Math.abs((sum / I) * (sum / I) - sum2 / I) / I);
  const logA = Math.log(W * H);
  const metrics = [o.sens * o.dens,
    acc.p / I, se(acc.p, acc.p2), acc.n / I, se(acc.n, acc.n2),
    1 - (acc.e / I) / logA, se(acc.e, acc.e2) / logA,
    acc.b / I, se(acc.b, acc.b2), acc.i / I, se(acc.i, acc.i2)];
  return {sim, metrics};
}

// The empty lattice (density 0), which the reference runs: its polar, nematic and band
// divide by N = 0 and print NaN, and so do their errors; spatial order is 1 with error 0,
// 1/occupied and its error are 0. Lattice and generator are compared as for the other
// scenarios.
const emptyScenarios = [
  {model:"boson",   tsteps:3, iters:2, sens:2,   dens:0, seed:12345},
  {model:"fermion", tsteps:3, iters:2, sens:0.8, dens:0, seed:12345},
];
const EMPTY_NAN = [1, 2, 3, 4, 7, 8];       // polar, error, nematic, error, band, error

function compareWithCpp(sim, metrics, cpp, label, {empty = false} = {}) {
  assert.equal(cpp.state.length, W * H * 6, `${label}: lattice size differs`);
  const occ = sim.occ;
  let q = 0;
  for (let i = 0; i < W; i++) for (let j = 0; j < H; j++) for (let k = 0; k < 6; k++, q++) {
    const value = occ[(j * W + i) * 6 + k];
    assert.equal(value, cpp.state[q],
      `${label}: first lattice difference at site (${i},${j}), channel ${k}: WebAssembly=${value}, C++=${cpp.state[q]}`);
  }
  const rng = sim.rng.getState();
  assert.equal(cpp.rng.length, 33, `${label}: malformed C++ generator state`);
  assert.equal(rng[32], cpp.rng[0], `${label}: WELL state index differs`);
  for (let w = 0; w < 32; w++) assert.equal(rng[w], cpp.rng[w + 1], `${label}: WELL state word ${w} differs`);
  assert.equal(metrics.length, cpp.metrics.length, `${label}: metric count differs`);
  if (empty) {
    for (let m = 0; m < metrics.length; m++) {
      if (EMPTY_NAN.includes(m)) {
        assert(Number.isNaN(cpp.metrics[m]), `${label}: the reference prints ${cpp.metrics[m]} for metric ${m}, not NaN`);
        assert(Number.isNaN(metrics[m]), `${label}: metric ${m} is ${metrics[m]}; the reference gives NaN for an empty lattice`);
      } else {
        assert.equal(metrics[m], cpp.metrics[m], `${label}: metric ${m} differs: WebAssembly=${metrics[m]} C++=${cpp.metrics[m]}`);
      }
    }
    assert.deepEqual(cpp.metrics.filter((_, m) => !EMPTY_NAN.includes(m)), [0, 1, 0, 0, 0],
      `${label}: the reference's finite values (sens*dens, spatial, its error, 1/occupied, its error)`);
    return;
  }
  for (let m = 0; m < metrics.length; m++) {
    const error = Math.abs(metrics[m] - cpp.metrics[m]);
    assert(error <= 5e-14 * Math.max(1, Math.abs(cpp.metrics[m])),
      `${label}: metric ${m} differs: WebAssembly=${metrics[m]} C++=${cpp.metrics[m]} error=${error}`);
  }
}

// ------------------------------------------------------------------------------------
//  3. Replay round trip (verify-parity.mjs's, on the WebAssembly engine)
// ------------------------------------------------------------------------------------

// The worker rewinds by streaming a saved post-collision lattice into occ, copying it
// into src, and restoring the generator and t (simWorkerMain, request "jump").
function verifyReplay(engine) {
  const sens = 2;
  const sim = new WasmSim(engine, {model:"boson", W:60, H:60, dens:0.4, seed:12345,
    kernel:"sum", bosonField:"site", bosonAlign:"nematic"});
  const saved = [];
  for (let t = 1; t <= 8; t++) {
    sim.step(sens);
    saved[t] = {occ:sim.occ.slice(), src:sim.src.slice(), rng:sim.rng.getState(), N:sim.nuparts};
  }
  sim.occ.fill(0);
  streamInto(saved[4].src, sim.occ, sim.W, sim.H);
  assert.deepEqual(sim.occ, saved[4].occ, "replay: streaming the saved state does not give the saved lattice");
  sim.src.set(saved[4].src);
  sim.rng.setState(saved[4].rng);
  sim.t = 4;
  sim.step(sens);
  assert.deepEqual(sim.src, saved[5].src, "replay did not reproduce the next collision state");
  assert.deepEqual(sim.occ, saved[5].occ, "replay did not reproduce the next streamed state");
  assert.deepEqual(sim.rng.getState(), saved[5].rng, "replay did not reproduce the generator state");
  assert.equal(sim.t, 5, "replay: t after the step");
  assert.equal(sim.nuparts, saved[5].N, "replay: N changed");
}

// ------------------------------------------------------------------------------------
//  3b. The research options, through the page's WasmSim
// ------------------------------------------------------------------------------------

// The 10 scenarios use the default options only. Here every combination of kernel
// (power at alpha 0, 0.85 and 1), boson field and alignment, for both models, runs on both engines and is compared
// after every step; then one boson run changes the options between steps as the worker
// does (applyCfg sets them one by one), through every field. This tests WasmSim's
// translation of the option names, which native/verify-wasm.mjs (its own wrapper) cannot.
function optionRun(engine, model, opt, label) {
  const options = {model, W:30, H:24, dens:model === "boson" ? 0.6 : 0.3, seed:2718, ...opt};
  const sim = new WasmSim(engine, options);
  const js = new LGCA(options);
  try {
    compareEngines(sim, js, `${label} init`);
    return {sim, js};
  } catch (err) { sim.destroy(); throw err; }
}

function compareMeasures(sim, js, label) {
  const m = sim.measure(), mj = js.measure();
  for (const key of ["polar", "nematic", "occFrac", "invOcc"])
    assert.equal(m[key], mj[key], `${label}: ${key} differs from the JavaScript engine`);
  for (const key of ["entropy", "spatial"])
    assert(relativeClose(m[key], mj[key]), `${label}: ${key} differs from the JavaScript engine: ${m[key]} vs ${mj[key]}`);
  assert.equal(sim.band(), js.band(), `${label}: band differs from the JavaScript engine`);
}

function verifyOptions(engine) {
  let runs = 0;
  for (const model of ["boson", "fermion"])
    // power at alpha 0 and 1 are special cases (alpha 0 is the "polar neighbours" preset)
    for (const [kernel, alpha] of [["sum", 0.5], ["avg", 0.5], ["power", 0], ["power", 0.85], ["power", 1]])
      for (const bosonField of ["site", "neigh", "both"])
        for (const bosonAlign of ["nematic", "polar"]) {
          const label = `options ${model} ${kernel} alpha ${alpha} ${bosonField} ${bosonAlign}`;
          const {sim, js} = optionRun(engine, model, {kernel, alpha, bosonField, bosonAlign}, label);
          try {
            for (let t = 1; t <= 6; t++) {
              sim.step(3); js.step(3);
              compareEngines(sim, js, `${label} step ${t}`);
            }
            compareMeasures(sim, js, label);
          } finally { sim.destroy(); }
          runs++;
        }

  // Options changed mid-run, each set for 3 steps.
  const phases = [
    {kernel:"sum",   alpha:0.5,  bosonField:"site",  bosonAlign:"nematic"},
    {kernel:"power", alpha:1.3,  bosonField:"neigh", bosonAlign:"polar"},
    {kernel:"avg",   alpha:1.3,  bosonField:"both",  bosonAlign:"nematic"},
    {kernel:"sum",   alpha:1.3,  bosonField:"neigh", bosonAlign:"nematic"},
    {kernel:"power", alpha:0.2,  bosonField:"site",  bosonAlign:"polar"},
    {kernel:"power", alpha:2,    bosonField:"neigh", bosonAlign:"polar"},
    {kernel:"power", alpha:0,    bosonField:"neigh", bosonAlign:"polar"},
    {kernel:"power", alpha:1,    bosonField:"both",  bosonAlign:"nematic"},
    {kernel:"avg",   alpha:0,    bosonField:"site",  bosonAlign:"nematic"},
  ];
  const {sim, js} = optionRun(engine, "boson", phases[0], "option switches");
  try {
    phases.forEach((phase, p) => {
      for (const [name, value] of Object.entries(phase)) { sim[name] = value; js[name] = value; }
      for (let t = 1; t <= 3; t++) {
        sim.step(4); js.step(4);
        compareEngines(sim, js, `option switches, phase ${p} (${phase.bosonField} ${phase.bosonAlign} ${phase.kernel}) step ${t}`);
      }
    });
    compareMeasures(sim, js, "option switches");
  } finally { sim.destroy(); }
  return runs;
}

// ------------------------------------------------------------------------------------
//  3c. Boson sampling on the boundaries of the channel distribution
// ------------------------------------------------------------------------------------

// The reference takes the first channel c with r <= cdf[c] (sample_channel), so a draw
// exactly on a boundary, or on a boundary repeated by zero-weight channels, belongs to the
// lower channel. One particle at one site with sensitivity 1000: every weight but the
// largest underflows to 0, so the distribution is exact. The draw is forced through the
// generator state, the same way on both engines: with every word 0 but S[i] = x, the next
// output is x ^ (x << 7) (mod 2^32), and x = w ^ w<<7 ^ w<<14 ^ w<<21 ^ w<<28 makes it w.
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
function verifyCdfBoundaries(engine) {
  const {makeWell} = new Function(`${scriptBlock("engineSrc")}\nreturn {makeWell};`)();
  let checks = 0;
  for (const [name, make] of [["WebAssembly", o => new WasmSim(engine, o)], ["JavaScript", o => new LGCA(o)]])
    for (const {channel, cdf, draws} of CDF_CASES)
      for (const [word, expected] of draws) {
        const well = makeWell(0); well.setState(stateForDraw(word));
        const r = well();
        assert.equal(r, word * 2.32830643653869628906e-10, "the forced generator state does not give the wanted draw");
        assert.equal(ruleChannel(r, cdf), expected, "the test's own expectation disagrees with the reference rule");
        const got = forcedCollision(make, channel, word);
        if (got.prob) assert.deepEqual(got.prob, cdf, `${name}: channel ${channel}: the distribution differs from ${cdf}`);
        assert.equal(got.channel, expected, `${name}: particle in channel ${channel}, draw ${r} (cdf ${cdf.join(" ")}): ` +
          `took channel ${got.channel}, the reference rule gives ${expected}`);
        assert.deepEqual(got.rng, well.getState(), `${name}: the collision did not consume exactly the one forced draw`);
        checks++;
      }
  return checks;
}

// The generator state contract (rng.setState, both engines): exactly 33 integers from 0 to 2^32 - 1, the
// last (the index) below 32. Anything else throws a RangeError and leaves the generator as
// it was; a valid state, typed or plain array, is taken word for word.
function badRngStates(good) {
  const edit = (q, v) => Array.from(good, (x, i) => i === q ? v : x);
  const hole = Array.from(good); delete hole[3];
  // every word different from the current state, so a setter that writes before it has
  // checked everything shows (the index stays valid until the last word)
  const other = Uint32Array.from(good, (x, i) => i < 32 ? ~x >>> 0 : (x + 7) % 32);
  const late = (q, v) => Array.from(other, (x, i) => i === q ? v : x);
  return [["32 words", good.slice(0, 32)], ["34 words", Uint32Array.from([...good, 0])],
    ["index 32", edit(32, 32)], ["index -1", edit(32, -1)], ["a word of 2^32", edit(3, 2 ** 32)],
    ["a negative word", edit(3, -1)], ["a fractional word", edit(3, 1.5)], ["a NaN word", edit(3, NaN)],
    ["a missing word", edit(3, undefined)], ["a hole", hole], ["null", null], ["no words", []],
    ["new words, then index 32", late(32, 32)], ["new words, then a NaN index", late(32, NaN)],
    ["new words, the last 2^32", late(31, 2 ** 32)], ["32 new words", other.slice(0, 32)],
    ["34 words, new ones", Uint32Array.from([...other, 0])]];
}
function verifyRngStateContract(make, engine) {
  const options = {model:"boson", W:12, H:10, dens:0.4, seed:4242, kernel:"sum", bosonField:"site", bosonAlign:"nematic"};
  const sim = make(options), twin = new LGCA(options);
  try {
    for (let t = 0; t < 3; t++) { sim.step(2); twin.step(2); }
    const good = sim.rng.getState();
    const cases = badRngStates(good);
    for (const [what, bad] of cases) {
      assert.throws(() => sim.rng.setState(bad), RangeError, `${engine}: rng.setState accepted ${what}`);
      assert.deepEqual(Array.from(sim.rng.getState()), Array.from(good), `${engine}: a refused state (${what}) changed the generator`);
    }
    sim.rng.setState(Array.from(good));
    for (let t = 0; t < 2; t++) { sim.step(2); twin.step(2); }
    assert.deepEqual(Array.from(sim.rng.getState()), Array.from(twin.rng.getState()), `${engine}: the generator moved after the refusals`);
    assert.deepEqual(Array.from(sim.occ), Array.from(twin.occ), `${engine}: the lattice moved after the refusals`);
    const odd = Uint32Array.from(good, (x, i) => i === 0 ? 0 : i === 1 ? 0xffffffff : i === 32 ? 31 : x);
    sim.rng.setState(odd);
    assert.deepEqual(Array.from(sim.rng.getState()), Array.from(odd), `${engine}: a valid state was not taken word for word`);
    return cases.length;
  } finally { if (sim.destroy) sim.destroy(); }
}

// ------------------------------------------------------------------------------------
//  3d. Density diagnostics (densityStats), against numbers worked out by hand
// ------------------------------------------------------------------------------------

// The worker reports these with every state (the density readout, the CSV): N, the most
// particles at a site (max) and in one channel of a site (maxChannel, the studies' kmax),
// occupied sites, mean per site, max/N and the participation area N^2 / sum of n_x^2.
// Four 2 x 2 lattices (four sites, channels listed per site), each with its expected values.
function verifyDensityStats() {
  const lattice = sites => Int32Array.from(sites.flat());
  const empty = [0, 0, 0, 0, 0, 0];
  const cases = [
    ["empty", lattice([empty, empty, empty, empty]),
     {N:0, max:0, maxChannel:0, occupied:0, mean:0, fraction:0, effectiveArea:0}],
    ["one pile", lattice([empty, empty, [3, 0, 0, 2, 0, 0], empty]),
     {N:5, max:5, maxChannel:3, occupied:1, mean:1.25, fraction:1, effectiveArea:1}],
    ["unequal piles", lattice([[1, 0, 0, 0, 0, 0], [2, 2, 0, 0, 0, 0], empty, [0, 0, 0, 0, 0, 3]]),
     {N:8, max:4, maxChannel:3, occupied:3, mean:2, fraction:0.5, effectiveArea:64 / 26}],
    ["uniform", lattice([[1, 1, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0], [1, 0, 0, 0, 0, 1]]),
     {N:8, max:2, maxChannel:1, occupied:4, mean:2, fraction:0.25, effectiveArea:4}],
  ];
  for (const [what, arr, expected] of cases)
    assert.deepEqual(densityStats(arr), expected, `density diagnostics of ${what} lattice`);
  return cases.length;
}

// ------------------------------------------------------------------------------------
//  4. The worker protocol
// ------------------------------------------------------------------------------------

// simWorkerMain on a stand-in for the worker's port: requests go in through onmessage,
// replies are collected. Nothing is transferred, as on a MessagePort without transfers.
function startWorker(glue, main = simWorkerMain) {
  const replies = [];
  const port = {onmessage:null, postMessage(msg) { replies.push(msg); }};
  main(port, glue);
  return {replies, send(msg) { port.onmessage({data:msg}); }};
}

async function until(done, what, ms = 120000) {
  const start = Date.now();
  while (!done()) {
    if (Date.now() - start > ms) throw new Error(`timed out waiting for ${what}`);
    await new Promise(resolve => setTimeout(resolve, 1));
  }
}

const settings = extra => ({sens:2, kernel:"sum", alpha:0.5, bosonField:"site", bosonAlign:"nematic", ...extra});
const buffers = (w, h) => ({src:new ArrayBuffer(w * h * 6 * 4), occ:new ArrayBuffer(w * h * 6 * 4)});

// A session of requests as the page sends them: start, run, switch options mid-run (a
// deliberate branch), rewind, run on with other options, look at a buffered tick; then a second init with
// the other model on the same worker (which frees the first engine), and the same again.
function requests(engine) {
  const boson = {model:"boson", W:48, H:40, dens:0.4, seed:777, kernel:"sum", alpha:0.5,
                 bosonField:"site", bosonAlign:"nematic", engine};
  const fermion = {model:"fermion", W:36, H:30, dens:0.3, seed:4242, kernel:"avg", alpha:0.5,
                   bosonField:"site", bosonAlign:"nematic", engine};
  const bosonCfg = settings();
  const branch = settings({sens:6, kernel:"power", alpha:0.85, bosonField:"both", bosonAlign:"polar"});
  const neighbours = settings({sens:5, kernel:"avg", bosonField:"neigh", bosonAlign:"polar"});
  const fermionCfg = settings({sens:1.5, kernel:"avg"});
  const fermionBranch = settings({sens:3, kernel:"power", alpha:1});
  const list = [
    {type:"init", params:boson, cfg:bosonCfg, band:true, bufs:buffers(48, 40)},
    {type:"step", n:25, ms:1e9, cfg:bosonCfg, band:true, bufs:buffers(48, 40)},
    {type:"step", n:12, ms:1e9, cfg:branch, band:true, bufs:buffers(48, 40)},
    {type:"jump", k:20, bufs:buffers(48, 40)},
    {type:"step", n:9, ms:1e9, cfg:neighbours, band:true, bufs:buffers(48, 40)},
    {type:"frame", bufs:buffers(48, 40)},
    {type:"slot", k:22, band:true, buf:new ArrayBuffer(48 * 40 * 6 * 4)},
    {type:"jump", k:5},
    {type:"jump", k:1000, bufs:buffers(48, 40)},
    {type:"step", n:3, ms:1e9, cfg:bosonCfg, band:false},
    {type:"init", params:fermion, cfg:fermionCfg, band:true, bufs:buffers(36, 30)},
    {type:"step", n:15, ms:1e9, cfg:fermionCfg, band:true, bufs:buffers(36, 30)},
    {type:"jump", k:7, bufs:buffers(36, 30)},
    {type:"step", n:4, ms:1e9, cfg:fermionBranch, band:true, bufs:buffers(36, 30)},
    {type:"slot", k:9},
    {type:"frame"},
  ];
  return list.map((msg, seq) => ({...msg, seq, epoch:1}));
}

// Same replies: every field equal, typed arrays element for element; entropy and the
// spatial order (computed from it) to 1e-12 relative. Returns the first difference.
function difference(a, b, path = "") {
  const loose = /(^|\.)(entropy|spatial)$/.test(path);
  if (typeof a === "number" && typeof b === "number")
    return (loose ? relativeClose(a, b) : a === b || (Number.isNaN(a) && Number.isNaN(b))) ? null
      : `${path}: ${a} vs ${b}`;
  if (ArrayBuffer.isView(a) || ArrayBuffer.isView(b)) {
    if (!ArrayBuffer.isView(a) || !ArrayBuffer.isView(b) || a.constructor !== b.constructor || a.length !== b.length)
      return `${path}: arrays of different kinds or sizes`;
    for (let q = 0; q < a.length; q++) if (a[q] !== b[q]) return `${path}[${q}]: ${a[q]} vs ${b[q]}`;
    return null;
  }
  if (a && b && typeof a === "object" && typeof b === "object") {
    if (Array.isArray(a) !== Array.isArray(b)) return `${path}: array vs object`;
    const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
    for (const key of keys) {
      const d = difference(a[key], b[key], path ? `${path}.${key}` : key);
      if (d) return d;
    }
    return null;
  }
  return a === b ? null : `${path}: ${JSON.stringify(a)} vs ${JSON.stringify(b)}`;
}

// The replies of a worker after the requests, without the engine fields of init replies
// (compared on their own).
function withoutEngine(replies) {
  return replies.map(({engine, engineNote, ...rest}) => rest);
}

async function session(glue, list, {oneAtATime = false, main = simWorkerMain, ms = 120000} = {}) {
  const worker = startWorker(glue, main);
  assert.equal(worker.replies[0]?.type, "hello", "the worker did not say hello first");
  for (const msg of list) {
    worker.send(msg);
    if (oneAtATime) await until(() => worker.replies.length === 2 + list.indexOf(msg), `the reply to ${msg.type}`, ms);
  }
  await until(() => worker.replies.length === 1 + list.length, "all replies", ms);
  await new Promise(resolve => setTimeout(resolve, 20));   // no stray extra replies
  assert.equal(worker.replies.length, 1 + list.length, "the worker sent more replies than requests");
  const replies = worker.replies.slice(1);
  replies.forEach((reply, q) => {
    assert.equal(reply.type, list[q].type, `reply ${q} is not for request ${q} (order kept?)`);
    assert.equal(reply.seq, list[q].seq, `reply ${q} carries the wrong seq`);
  });
  return replies;
}

// The documented size of the worker's replay buffer: lattices of about 24 MB, at least 20
// and at most 300 ticks.
const capacity = (w, h) => Math.min(300, Math.max(20, Math.floor(24e6 / (w * h * 6 * 4))));

// The replies a correct worker sends for `list`, worked out without the worker's code:
// a fresh JavaScript LGCA per init, stepped as the requests say. A rewind to tick k, and
// a look at buffered tick k ("slot"), run the seed again to k with the settings each step
// had, instead of restoring anything the worker saved; so the lattices, the generator
// words and the measurements after a rewind come from an independent run. The buffer
// holds the newest `capacity` ticks: after init only tick 1, each step adds one and drops
// the oldest once full, a rewind to k keeps the ticks up to k. `startT` starts the step
// counter elsewhere than 0 (for a worker made to start there); `did` gives, by seq, the
// steps a time-limited request took (the worker decides that; every other request takes
// all n).
function expectedReplies(list, {startT = 0, did = new Map()} = {}) {
  const replies = [];
  let params = null, perStep = [], sim = null, lastM = null, lastRho = null, oldest = 0, cap = 0;
  const setOptions = (s, cfg) => {
    s.kernel = cfg.kernel; s.alpha = cfg.alpha; s.bosonField = cfg.bosonField; s.bosonAlign = cfg.bosonAlign;
  };
  const fresh = () => { const s = new LGCA(params); s.t = startT; return s; };
  const rerun = k => {
    const s = fresh();
    for (let t = startT + 1; t <= k; t++) { const cfg = perStep[t - startT - 1]; setOptions(s, cfg); s.step(cfg.sens); }
    return s;
  };
  const step = cfg => {
    setOptions(sim, cfg); sim.step(cfg.sens); perStep.push(cfg);
    oldest = Math.max(oldest, sim.t - cap + 1);
  };
  const record = (cfg, band) => {
    const m = sim.measure(), fermion = sim.model === "fermion";
    lastM = m; lastRho = densityStats(sim.occ);
    const rec = {t:sim.t, polar:m.polar, nematic:m.nematic, spatial:m.spatial, rho:lastRho,
      cfg:{sens:cfg.sens, alpha:cfg.kernel === "sum" ? 0 : cfg.kernel === "avg" ? 1 : cfg.alpha, kernel:cfg.kernel,
           field:fermion ? "neigh" : cfg.bosonField, alignment:fermion ? "polar" : cfg.bosonAlign}};
    if (band && sim.t % 10 === 0) rec.band = sim.band();
    return rec;
  };
  const reply = (msg, extra) => {
    const out = {type:msg.type, seq:msg.seq, epoch:msg.epoch, t:sim.t, nuparts:sim.nuparts, rng:sim.rng.getState(),
                 m:lastM, rho:lastRho, range:[oldest, sim.t], ...extra};
    if (msg.bufs) { out.src = sim.src.slice(); out.occ = sim.occ.slice(); out.frameT = sim.t; }
    replies.push(out);
  };
  const buffered = k => Number.isInteger(k) && k >= oldest && k <= sim.t;
  for (const msg of list) {
    if (msg.type === "init") {
      params = msg.params; perStep = []; cap = capacity(params.W, params.H);
      sim = fresh(); oldest = startT + 1;
      setOptions(sim, msg.cfg);
      const records = [record(msg.cfg, msg.band)];
      step(msg.cfg);
      records.push(record(msg.cfg, msg.band));
      reply(msg, {records});
    } else if (msg.type === "step") {
      const records = [], n = did.has(msg.seq) ? did.get(msg.seq) : msg.n;
      for (let q = 0; q < n; q++) { step(msg.cfg); records.push(record(msg.cfg, msg.band)); }
      reply(msg, {records, did:n});
    } else if (msg.type === "jump") {
      if (!buffered(msg.k)) { reply(msg, {ok:false}); continue; }
      perStep.length = msg.k - startT;
      sim = rerun(msg.k);
      lastM = sim.measure(); lastRho = densityStats(sim.occ);
      reply(msg, {ok:true, k:msg.k});
    } else if (msg.type === "slot") {
      const at = buffered(msg.k) ? rerun(msg.k) : null;
      replies.push({type:"slot", seq:msg.seq, epoch:msg.epoch, k:at ? msg.k : null,
                    src:at ? at.src.slice() : new Int32Array(sim.W * sim.H * 6), rng:at ? at.rng.getState() : null,
                    m:at ? at.measure() : null, rho:at ? densityStats(at.occ) : null, band:at&&msg.band ? at.band() : null});
    } else if (msg.type === "frame") {
      reply(msg, {});
    }
  }
  return replies;
}

// Deliberately broken copies of the worker block: each must fail the comparison with the
// independent replies, on both engines, or that comparison could not see the defect.
const BREAKAGES = [
  ["the generator state is not restored on a rewind", "sim.rng.setState(rb.rng[i]);", ""],
  ["a look at a buffered tick gets the current generator state", "rng=new Uint32Array(rb.rng[i]);", "rng=sim.rng.getState();"],
];

// Sessions beyond the main one. Each runs on both engines against expectedReplies; each
// names the deliberately broken worker it must catch (a change that left every earlier
// check passing).
const cfgA = settings(), cfgB = settings({sens:5, kernel:"power", alpha:0.85, bosonField:"neigh", bosonAlign:"polar"});
const steps = (n, cfg, extra = {}) => ({type:"step", n, ms:1e9, cfg, band:true, ...extra});
const numbered = (list, epoch) => list.map((msg, seq) => ({epoch, ...msg, seq}));
// Past the buffer: 300 ticks at 12 x 10. At t = 401 it holds 102..401; 101 is gone. Rewind
// to the oldest tick, branch with other settings, fill the ring again past the rewind point
// (its indices wrap twice), and rewind and look again.
function rolloverRequests(engine) {
  const params = {model:"boson", W:12, H:10, dens:0.8, seed:31337, kernel:"sum", alpha:0.5,
                  bosonField:"site", bosonAlign:"nematic", engine};
  return numbered([
    {type:"init", params, cfg:cfgA, band:true, bufs:buffers(12, 10)},
    steps(150, cfgA), steps(250, cfgA, {bufs:buffers(12, 10)}),
    {type:"slot", k:101, band:true}, {type:"slot", k:102, band:true}, {type:"slot", k:401, band:false},
    {type:"jump", k:101}, {type:"jump", k:102, bufs:buffers(12, 10)},
    steps(30, cfgB), {type:"slot", k:110, band:true}, steps(300, cfgB),
    {type:"jump", k:132}, {type:"jump", k:133}, {type:"slot", k:200, band:true},
    steps(7, cfgA, {bufs:buffers(12, 10)}), {type:"frame", bufs:buffers(12, 10)},
  ], 7);
}
// Across step 2^31 (the engines count to 2^53): a worker started at t = 2^31 - 4.
const BIG_T = 2 ** 31 - 4;
function bigTickRequests(engine) {
  const params = {model:"fermion", W:20, H:16, dens:0.3, seed:2024, kernel:"sum", alpha:0.5,
                  bosonField:"site", bosonAlign:"nematic", engine};
  return numbered([
    {type:"init", params, cfg:settings({sens:1.2}), band:true, bufs:buffers(20, 16)},
    steps(6, settings({sens:1.2})), {type:"slot", k:2 ** 31, band:true}, {type:"slot", k:2 ** 31 + 3},
    {type:"jump", k:2 ** 31 - 1, bufs:buffers(20, 16)}, steps(3, settings({sens:2.5})),
    {type:"jump", k:2 ** 31 + 2}, {type:"frame", bufs:buffers(20, 16)},
  ], 8);
}
// A worker (from the worker block's text, possibly broken) whose engines start counting at
// t0: both engine classes set t right after they are made, as a rewind does.
const workerFrom = (text, t0 = 0) => loadLab(t0 ? `${text}
{ const init = LGCA.prototype.initLattice; LGCA.prototype.initLattice = function () { init.call(this); this.t = ${t0}; }; }
WasmSim = class extends WasmSim { constructor(wasm, o) { super(wasm, o); this.t = ${t0}; } };` : text).simWorkerMain;
// Batches cut short by their time budget: ms 0 takes exactly one step (the worker always
// takes one), a short budget some number between 1 and n that the worker decides.
function partialRequests(engine) {
  const params = {model:"boson", W:48, H:40, dens:0.9, seed:5150, kernel:"sum", alpha:0.5,
                  bosonField:"site", bosonAlign:"nematic", engine};
  return numbered([
    {type:"init", params, cfg:cfgA, band:true, bufs:buffers(48, 40)},
    steps(50, cfgA, {ms:0, bufs:buffers(48, 40)}), steps(400, cfgB, {ms:3}),
    {type:"jump", k:2, bufs:buffers(48, 40)}, steps(9, cfgA, {ms:0}), steps(5, cfgB),
    {type:"slot", k:3, band:true}, {type:"frame", bufs:buffers(48, 40)},
  ], 9);
}
// The engine switched while the module loads (labEngine("js") on the page): an init for
// the WebAssembly engine, then one for the JavaScript engine under a new epoch, all sent
// before the module has loaded. Both are answered, in order, the second on JavaScript.
function switchRequests() {
  const params = engine => ({model:"boson", W:30, H:24, dens:0.5, seed:808, kernel:"sum", alpha:0.5,
                             bosonField:"site", bosonAlign:"nematic", engine});
  return [
    {type:"init", params:params("wasm"), cfg:cfgA, band:true, bufs:buffers(30, 24), epoch:10},
    {...steps(5, cfgA), epoch:10},
    {type:"init", params:params("js"), cfg:cfgA, band:true, bufs:buffers(30, 24), epoch:11},
    {...steps(12, cfgA), epoch:11}, {type:"jump", k:6, epoch:11}, {...steps(4, cfgB), epoch:11},
    {type:"frame", bufs:buffers(30, 24), epoch:11},
  ].map((msg, seq) => ({...msg, seq}));
}
const MORE_BREAKAGES = [
  ["the ring's head advances twice once the buffer is full", rolloverRequests,
   "rb.head=(rb.head+1)%rb.cap;", "rb.head=(rb.head+(rb.len===rb.cap?2:1))%rb.cap;"],
  ["replay ticks are kept as 32-bit integers", bigTickRequests,
   "ticks:new Float64Array(cap)", "ticks:new Int32Array(cap)"],
  ["a time-limited batch reports n steps instead of those it took", partialRequests,
   "reply(msg,{records,did},msg.bufs);", "reply(msg,{records,did:msg.n},msg.bufs);"],
];

async function verifyProtocol(glue) {
  const results = [];
  const expected = expectedReplies(requests("js"));
  assert(expected[3].ok === true && expected[8].ok === false, "the session's rewinds do not succeed and fail as intended");
  const js = await session(glue, requests("js"));
  for (const r of js.filter(r => r.type === "init"))
    assert(r.engine === "js" && r.engineNote === "", "the JavaScript engine's init reply names another engine");
  const dJs = difference(withoutEngine(js), expected);
  assert.equal(dJs, null, `JavaScript worker differs from the independent run: ${dJs}`);
  results.push(`PASS worker protocol, JavaScript engine: every reply equals independently stepped engines, ` +
    `rewinds and buffered ticks run again from the seed (${js.length} requests)`);

  // All requests sent at once: the first init starts loading the module and the others
  // wait for it, in order. Then a second worker, one request at a time.
  for (const [how, options] of [["all requests sent while the module loads", {}],
                                ["one request at a time", {oneAtATime:true}]]) {
    const wasm = await session(glue, requests("wasm"), options);
    for (const r of wasm.filter(r => r.type === "init"))
      assert(r.engine === "wasm" && r.engineNote === "", `init reply: engine ${r.engine} (${r.engineNote})`);
    const d = difference(withoutEngine(wasm), expected);
    assert.equal(d, null, `WebAssembly worker (${how}) differs from the independent run: ${d}`);
    results.push(`PASS worker protocol, WebAssembly engine, ${how}: every reply equals independently stepped ` +
      `engines (${wasm.length} requests)`);
  }

  // The comparison catches a broken worker: each breakage, on each engine.
  for (const [what, text, replacement] of BREAKAGES) {
    const workerText = scriptBlock("workerSrc");
    assert.equal(workerText.split(text).length, 2, `the worker block no longer holds "${text}" once (update BREAKAGES)`);
    const broken = loadLab(workerText.replace(text, replacement)).simWorkerMain;
    for (const engine of ["js", "wasm"]) {
      const replies = await session(glue, requests(engine), {main:broken});
      assert.notEqual(difference(withoutEngine(replies), expected), null,
        `a worker in which ${what} (${engine}) passes the comparison with the independent run`);
    }
    results.push(`PASS a worker in which ${what} fails the comparison, on both engines`);
  }

  // Fallback: no block, a block that fails to load, and settings the engine refuses (a
  // 3x3 lattice; the JavaScript engine runs it). The JavaScript engine then runs, the
  // init reply says so, and the replies are those of the JavaScript worker.
  const short = engine => requests(engine).slice(0, 7);
  const jsShort = expected.slice(0, 7);
  for (const [what, brokenGlue, expect] of [
    ["no wasmEngine block", "", "the page has no wasmEngine block"],
    ["a block that fails to load", 'throw new Error("broken on purpose");', "broken on purpose"]]) {
    const replies = await session(brokenGlue, short("wasm"));
    assert(replies[0].engine === "js" && replies[0].engineNote.includes(expect),
      `${what}: init reply engine ${replies[0].engine}, note "${replies[0].engineNote}"`);
    const d = difference(withoutEngine(replies), jsShort);
    assert.equal(d, null, `${what}: the fallback differs from the independent run: ${d}`);
    results.push(`PASS worker falls back to the JavaScript engine with ${what}, and says so`);
  }
  const tiny = {model:"boson", W:3, H:3, dens:0.4, seed:9, kernel:"sum", alpha:0, bosonField:"site", bosonAlign:"nematic"};
  const list = [
    {type:"init", params:{...tiny, engine:"wasm"}, cfg:settings(), band:true, bufs:buffers(3, 3)},
    {type:"step", n:5, ms:1e9, cfg:settings(), band:true, bufs:buffers(3, 3)},
    requests("wasm")[0],
    requests("wasm")[1],
  ].map((msg, seq) => ({...msg, seq, epoch:2}));
  const refused = await session(glue, list);
  assert(refused[0].engine === "js" && refused[0].engineNote.includes("refused these settings"),
    `refused settings: init reply engine ${refused[0].engine}, note "${refused[0].engineNote}"`);
  assert(refused[2].engine === "wasm", "after a refusal, the next valid init did not get the WebAssembly engine");
  const d = difference(withoutEngine(refused), expectedReplies(list));
  assert.equal(d, null, `refused settings: replies differ from the independent run: ${d}`);
  results.push("PASS worker falls back to the JavaScript engine for settings the WebAssembly engine refuses, then uses it again");

  // A module that fails to load once, then loads: an init without retryWasm does not try
  // again (the failure is kept); one with it does, and gets the WebAssembly engine. Both
  // ways of sending; each worker counts its own attempts.
  const init = extra => ({...requests("wasm")[0], ...extra});
  const steps = () => requests("wasm")[1];
  const retryList = [init(), steps(), init(), init({retryWasm:true}), steps(), {type:"jump", k:7, bufs:buffers(48, 40)},
                     steps(), init({retryWasm:true})].map((msg, seq) => ({...msg, seq, epoch:3}));
  const expectedRetry = expectedReplies(retryList);
  for (const [n, [how, options]] of [["all requests sent at once", {}], ["one request at a time", {oneAtATime:true}]].entries()) {
    const key = `__lgcaLoadAttempts${n}`;
    const failOnce = `globalThis.${key} = (globalThis.${key} || 0) + 1;\n` +
      `if (globalThis.${key} === 1) throw new Error("fails once on purpose");\n${glue}`;
    const replies = await session(failOnce, retryList, options);
    const engines = replies.filter(r => r.type === "init").map(r => r.engine).join(" ");
    assert(engines === "js js wasm wasm" && replies[0].engineNote.includes("fails once on purpose") &&
      replies[2].engineNote.includes("fails once on purpose") && replies[3].engineNote === "" && globalThis[key] === 2,
      `retry (${how}): init engines ${engines}, ${globalThis[key]} load attempts, notes ` +
      JSON.stringify(replies.filter(r => r.type === "init").map(r => r.engineNote)));
    const dRetry = difference(withoutEngine(replies), expectedRetry);
    assert.equal(dRetry, null, `retry (${how}): replies differ from the independent run: ${dRetry}`);
  }
  results.push("PASS after a failed load, an init with retryWasm loads the module (a module that fails once), " +
    "one without it does not try again; replies equal the independent run (both ways of sending)");

  // A module that never loads: the retry is tried once for its request, fails again, and
  // the requests after it are answered in order on the JavaScript engine (a worker that
  // retried the same request again and again would never answer them).
  for (const [n, [how, options]] of [["all requests sent at once", {}], ["one request at a time", {oneAtATime:true}]].entries()) {
    const key = `__lgcaFailedAttempts${n}`;
    const neverLoads = `globalThis.${key} = (globalThis.${key} || 0) + 1;\nthrow new Error("never loads, on purpose");\n`;
    const failList = [init(), steps(), init({retryWasm:true}), steps(), {type:"frame", bufs:buffers(48, 40)}]
      .map((msg, seq) => ({...msg, seq, epoch:4}));
    const replies = await session(neverLoads, failList, {...options, ms:15000});
    const engines = replies.filter(r => r.type === "init").map(r => r.engine).join(" ");
    assert(engines === "js js" && globalThis[key] === 2 &&
      replies.filter(r => r.type === "init").every(r => r.engineNote.includes("never loads, on purpose")),
      `retry of a module that never loads (${how}): init engines ${engines}, ${globalThis[key]} load attempts`);
    const dFail = difference(withoutEngine(replies), expectedReplies(failList));
    assert.equal(dFail, null, `retry of a module that never loads (${how}): replies differ from the independent run: ${dFail}`);
  }
  results.push("PASS a module that never loads: the retry is tried once (2 attempts in all), and the requests after it " +
    "are answered in order on the JavaScript engine (both ways of sending)");

  // Past the buffer's capacity, across step 2^31, and time-limited batches, on both
  // engines; then the broken worker each must catch.
  const extra = [
    ["past the buffer's capacity: evicted ticks refused, a rewind to the oldest tick, branches after the ring wraps",
     rolloverRequests, 0],
    ["across step 2^31: ticks, buffered range, rewinds and looks", bigTickRequests, BIG_T],
    ["with batches cut short by their time budget (did, records, lattice and generator)", partialRequests, 0],
  ];
  for (const [what, make, startT] of extra) {
    for (const engine of ["js", "wasm"]) {
      const list = make(engine);
      const replies = await session(glue, list, {main:workerFrom(scriptBlock("workerSrc"), startT)});
      for (const r of replies.filter(r => r.type === "init"))
        assert(r.engine === engine && r.engineNote === "", `${what}: init reply engine ${r.engine} (${r.engineNote})`);
      const problem = protocolProblem(list, replies, startT);
      assert.equal(problem, null, `worker protocol ${what} (${engine}): ${problem}`);
    }
    results.push(`PASS worker protocol ${what}: every reply equals independently stepped engines, on both engines`);
  }
  for (const [what, make, text, replacement] of MORE_BREAKAGES) {
    const workerText = scriptBlock("workerSrc");
    assert.equal(workerText.split(text).length, 2, `the worker block no longer holds "${text}" once (update MORE_BREAKAGES)`);
    const startT = make === bigTickRequests ? BIG_T : 0;
    const broken = workerFrom(workerText.replace(text, replacement), startT);
    for (const engine of ["js", "wasm"]) {
      const list = make(engine);
      let problem;
      try { problem = protocolProblem(list, await session(glue, list, {main:broken, ms:20000}), startT); }
      catch (err) { problem = err.message; }
      assert.notEqual(problem, null, `a worker in which ${what} (${engine}) passes the comparison with the independent run`);
    }
    results.push(`PASS a worker in which ${what} fails the comparison, on both engines`);
  }

  // The engine switched while the module loads: both inits answered in order, the second
  // on the JavaScript engine, every reply as worked out independently. A worker that let
  // the JavaScript init skip the queue would answer it before the earlier one.
  const switchList = switchRequests();
  const switched = await session(glue, switchList);
  const engines = switched.filter(r => r.type === "init").map(r => r.engine).join(" ");
  assert.equal(engines, "wasm js", `engine switched while the module loads: init replies on ${engines}`);
  const dSwitch = difference(withoutEngine(switched), expectedReplies(switchList));
  assert.equal(dSwitch, null, `engine switched while the module loads: replies differ from the independent run: ${dSwitch}`);
  {
    const text = "if(waiting){ waiting.push(msg); return; }";
    const workerText = scriptBlock("workerSrc");
    assert.equal(workerText.split(text).length, 2, `the worker block no longer holds "${text}" once`);
    const broken = workerFrom(workerText.replace(text,
      'if(waiting&&!(msg.type==="init"&&msg.params.engine!=="wasm")){ waiting.push(msg); return; }'));
    let problem;
    try { problem = difference(withoutEngine(await session(glue, switchList, {main:broken, ms:20000})), expectedReplies(switchList)); }
    catch (err) { problem = err.message; }
    assert.notEqual(problem, null, "a worker in which a JavaScript init skips the queue while the module loads passes");
  }
  results.push("PASS worker protocol, engine switched while the module loads: both inits answered in order " +
    "(WebAssembly, then JavaScript), replies equal the independent run; a worker letting the second init skip the queue fails");
  return results;
}

// The first difference between a session's replies and the independent run, or null.
// Time-limited steps (ms below 1e9) take what the worker decided, at least 1 and at most n,
// and exactly 1 when the budget is 0; those counts drive the independent run.
function protocolProblem(list, replies, startT) {
  const did = new Map();
  for (const [q, msg] of list.entries()) {
    if (msg.type !== "step" || msg.ms >= 1e9) continue;
    const n = replies[q].did;
    if (!(Number.isInteger(n) && n >= 1 && n <= msg.n)) return `request ${q}: did ${n} of ${msg.n}`;
    if (msg.ms === 0 && n !== 1) return `request ${q}: a zero time budget took ${n} steps, not 1`;
    did.set(msg.seq, n);
  }
  return difference(withoutEngine(replies), expectedReplies(list, {startT, did}));
}

// ------------------------------------------------------------------------------------
//  5. The page
// ------------------------------------------------------------------------------------

// The page's three scripts (engine, worker, page) run here as one function, as they share
// the page's global scope, against a small stand-in for the document: the elements of the
// page's own markup, canvases that draw nothing (a check can record their calls), browser
// storage in a Map, a clipboard that keeps what it is given. Node has a MessageChannel and
// no Worker, so the page starts its main-thread host, the fallback it uses where no worker
// can start: the same worker code and the real WebAssembly module. Animation frames run
// when a check asks for them (tick).

class FakeClassList {
  constructor(names) { this.names = new Set(names.split(/\s+/).filter(Boolean)); }
  add(...names) { for (const n of names) this.names.add(n); }
  remove(...names) { for (const n of names) this.names.delete(n); }
  contains(name) { return this.names.has(name); }
  toggle(name, force) {
    const on = force === undefined ? !this.names.has(name) : !!force;
    if (on) this.names.add(name); else this.names.delete(name);
    return on;
  }
}

// A 2D context that draws nothing; while `calls` is an array it records every method call
// and property set, which the drawing checks read.
function fakeContext() {
  const state = {calls:null};
  return new Proxy(state, {
    get(target, key) {
      if (key in target) return target[key];
      if (key === "measureText") return () => ({width:0});
      if (key === "createImageData" || key === "getImageData")
        return (...a) => { const [w, h] = a.length >= 4 ? a.slice(2) : a; return {width:w, height:h, data:new Uint8ClampedArray(w * h * 4)}; };
      return (...args) => { if (target.calls) target.calls.push([key, ...args]); };
    },
    set(target, key, value) {
      target[key] = value;
      if (target.calls && key !== "calls") target.calls.push([`=${String(key)}`, value]);
      return true;
    },
  });
}
class FakePath2D { rect() {} moveTo() {} lineTo() {} arc() {} }

// One selector as the page uses them: tag, .class, [attribute="value"], or a combination.
function matchesSelector(el, selector) {
  return selector.split(",").some(one => {
    const m = one.trim().match(/^([a-z0-9]+)?((?:\.[\w-]+)*)(?:\[([\w-]+)(?:="([^"]*)")?\])?$/i);
    if (!m) throw new Error(`the fake document does not know the selector ${one}`);
    const [, tag, classes, attr, value] = m;
    if (tag && el.localName !== tag.toLowerCase()) return false;
    for (const c of classes.split(".").filter(Boolean)) if (!el.classList.contains(c)) return false;
    if (attr) { const v = el.getAttribute(attr); if (v === null || (value !== undefined && v !== value)) return false; }
    return true;
  });
}

class FakeElement {
  constructor(tag, attrs = {}) {
    this.localName = tag; this.tagName = tag.toUpperCase(); this.attributes = {...attrs};
    this.children = []; this.parentNode = null; this.listeners = {}; this.style = {}; this.dataset = {};
    for (const [k, v] of Object.entries(attrs)) if (k.startsWith("data-")) this.dataset[k.slice(5)] = v;
    for (const decl of (attrs.style || "").split(";")) {
      const at = decl.indexOf(":");
      if (at > 0) this.style[decl.slice(0, at).trim()] = decl.slice(at + 1).trim();
    }
    this.classList = new FakeClassList(attrs.class || "");
    this.id = attrs.id || ""; this.title = attrs.title || "";
    this.checked = "checked" in attrs; this.open = "open" in attrs; this.selected = "selected" in attrs;
    this.width = +(attrs.width ?? 300); this.height = +(attrs.height ?? 150);
    this.ownText = ""; this.html = ""; this.files = [];
    this.storedValue = attrs.value;
    this.context = null;
  }
  get textContent() { return this.ownText + this.children.map(c => c.textContent).join(""); }
  set textContent(v) { this.children = []; this.ownText = String(v); }
  get innerHTML() { return this.html; }
  set innerHTML(v) { this.children = []; this.ownText = ""; this.html = String(v); }
  get value() {
    if (this.localName === "option") return this.attributes.value ?? this.textContent;
    if (this.localName === "select" && this.storedValue === undefined) {
      const options = this.options;
      return (options.find(o => o.selected) || options[0] || {value:""}).value;
    }
    return this.storedValue ?? "";
  }
  set value(v) { this.storedValue = String(v); }
  get max() { return this.attributes.max ?? ""; }
  set max(v) { this.attributes.max = String(v); }
  get min() { return this.attributes.min ?? ""; }
  set min(v) { this.attributes.min = String(v); }
  get options() { return this.children.filter(c => c.localName === "option"); }
  get nextSibling() { const p = this.parentNode; return p ? p.children[p.children.indexOf(this) + 1] || null : null; }
  getAttribute(name) {
    if (name.startsWith("data-") && this.dataset[name.slice(5)] !== undefined) return this.dataset[name.slice(5)];
    return this.attributes[name] ?? null;
  }
  setAttribute(name, v) { this.attributes[name] = String(v); }
  hasAttribute(name) { return name in this.attributes; }
  removeAttribute(name) { delete this.attributes[name]; }
  detach() { const p = this.parentNode; if (p) p.children.splice(p.children.indexOf(this), 1); this.parentNode = null; }
  appendChild(c) { c.detach(); c.parentNode = this; this.children.push(c); return c; }
  append(...cs) { for (const c of cs) if (c instanceof FakeElement) this.appendChild(c); }
  prepend(c) { c.detach(); c.parentNode = this; this.children.unshift(c); }
  insertBefore(c, ref) {
    c.detach(); c.parentNode = this;
    const at = ref ? this.children.indexOf(ref) : -1;
    if (at < 0) this.children.push(c); else this.children.splice(at, 0, c);
    return c;
  }
  descendants() { return this.children.flatMap(c => [c, ...c.descendants()]); }
  querySelectorAll(selector) { return this.descendants().filter(e => matchesSelector(e, selector)); }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
  matches(selector) { return matchesSelector(this, selector); }
  closest(selector) { for (let e = this; e; e = e.parentNode) if (e.localName !== "#root" && e.matches(selector)) return e; return null; }
  addEventListener(type, fn) { (this.listeners[type] ||= []).push(fn); }
  dispatch(type, init = {}) {
    const event = {type, target:this, preventDefault() {}, ...init};
    for (let e = this; e; e = e.parentNode) for (const fn of e.listeners[type] || []) fn(event);
    return event;
  }
  click() { this.dispatch("click"); }
  getContext(kind) { return kind === "2d" ? (this.context ||= fakeContext()) : null; }
  getBoundingClientRect() {
    return {left:0, top:0, width:parseFloat(this.style.width) || 660, height:parseFloat(this.style.height) || 572};
  }
  setPointerCapture() {}
  toDataURL() { return "data:,"; }
}

// The page's markup (the body up to its first script) as elements.
const VOID_ELEMENTS = new Set(["input", "br", "img", "hr", "meta", "link", "source", "wbr"]);
function parseMarkup(markup) {
  const root = new FakeElement("#root"), stack = [root];
  const pattern = /<!--[\s\S]*?-->|<(\/?)([a-zA-Z][\w-]*)((?:\s+[\w:-]+(?:="[^"]*")?)*)\s*\/?>|([^<]+)/g;
  for (const [whole, close, tagName, attrText, text] of markup.matchAll(pattern)) {
    if (whole.startsWith("<!--")) continue;
    const top = stack.at(-1);
    if (text !== undefined) { top.ownText += text.replace(/&amp;/g, "&"); continue; }
    const tag = tagName.toLowerCase();
    if (close) { const at = stack.findLastIndex(e => e.localName === tag); if (at > 0) stack.length = at; continue; }
    const attrs = {};
    for (const [, name, value] of attrText.matchAll(/([\w:-]+)(?:="([^"]*)")?/g)) attrs[name] = (value ?? "").replace(/&amp;/g, "&");
    const el = top.appendChild(new FakeElement(tag, attrs));
    if (!VOID_ELEMENTS.has(tag)) stack.push(el);
  }
  assert.equal(stack.length, 1, "the page's markup did not close every element (the fake document cannot read it)");
  return root;
}

const PAGE_MARKUP = html.slice(html.indexOf("<body>") + 6, html.indexOf('<script id="engineSrc">'));
const PAGE_SCRIPT = (() => { const m = html.match(/<script>([\s\S]*?)<\/script>/); assert(m, "no page script"); return m[1]; })();
// The colour variables of the light theme (the first :root block) and of the dark one.
const themeVars = block => Object.fromEntries([...block.matchAll(/(--[\w-]+):([^;]+);/g)].map(([, name, value]) => [name, value.trim()]));
const CSS_VARS = themeVars(html.match(/:root\{([\s\S]*?)\}/)[1]);
const DARK_VARS = themeVars(html.match(/prefers-color-scheme: dark\)\{\s*:root\{([\s\S]*?)\}/)[1]);
// What the checks read of the page's own state.
const PAGE_EXPORTS = `
;return {get sim(){return sim}, get hist(){return hist}, get replay(){return replay}, get pend(){return pend},
  get queued(){return queued}, get playing(){return playing}, get host(){return host}, get epoch(){return epoch},
  get view(){return view}, set view(v){view=v}, get snaps(){return snaps}, P, $, frame, reset, back, fwdN, jumpTo,
  enterReplay, queueSteps, applyRecords, drawChart, drawParticles, drawOverlay, labEngine:window.labEngine,
  HIST_LIMIT, BAND_LIMIT, ccv, ctx, off};`;

// Node's MessageChannel, but an error thrown by a message handler is recorded (as a
// browser reports it as a page error) instead of ending this process. The page's side
// (port1) also records what the page sends and receives.
function guardedChannel(out) {
  const guard = (port, page) => {
    let handler = null;
    return {
      postMessage:(msg, transfer) => {
        if (page) out.posted.push({type:msg.type, epoch:msg.epoch, k:msg.k, n:msg.n});
        port.postMessage(msg, transfer);
      },
      close:() => port.close(),
      get onmessage() { return handler; },
      set onmessage(fn) {
        handler = fn;
        port.onmessage = e => {
          if (page) out.received.push({type:e.data.type, epoch:e.data.epoch, engine:e.data.engine, ok:e.data.ok});
          try { fn(e); } catch (err) { out.pageErrors.push(String(err && err.message || err)); }
        };
      },
    };
  };
  return class { constructor() { const ch = new MessageChannel(); this.port1 = guard(ch.port1, true); this.port2 = guard(ch.port2, false); } };
}

class FakeStorage {
  constructor(entries = {}, refuse = false) { this.map = new Map(Object.entries(entries)); this.refuse = refuse; }
  getItem(key) { return this.map.has(key) ? this.map.get(key) : null; }
  setItem(key, value) { if (this.refuse) throw new Error("QuotaExceededError (refused on purpose)"); this.map.set(key, String(value)); }
  removeItem(key) { this.map.delete(key); }
}

// Open the page. `search` is its URL query, `storage` what browser storage holds and
// `refuse` makes storing fail, `before(document)` changes the markup before the page's
// scripts run (as a browser restoring form values would), `mutate` is [text, replacement]
// in the page script (a deliberately broken page), `inject` code run after the worker
// block (to make an engine fail), `dpr` the device pixel ratio. Returns the page's state
// (lab), its document, the elements it created, and what it sent (posted), received
// (received), warned and alerted.
function openLab(glue, {search = "", storage = {}, refuse = false, before = null, mutate = null, inject = "", dpr = 1} = {}) {
  const root = parseMarkup(PAGE_MARKUP);
  const byId = new Map(root.descendants().filter(e => e.id).map(e => [e.id, e]));
  const wasmBlock = new FakeElement("script", {id:"wasmEngine"}); wasmBlock.ownText = glue;
  byId.set("wasmEngine", wasmBlock);
  const rose = new FakeElement("canvas", {id:"roseCv"});
  const document = {
    documentElement:new FakeElement("html"),
    getElementById:id => byId.get(id) || (id === "roseCv" ? rose : null),
    querySelector:selector => root.querySelector(selector),
    createElement:tag => { const el = new FakeElement(tag.toLowerCase()); created.push(el); return el; },
    createTextNode:text => { const node = new FakeElement("#text"); node.ownText = String(text); return node; },
    addEventListener() {},
  };
  const created = [];
  if (before) before(document);
  const out = {document, rose, created, posted:[], received:[], warnings:[], alerts:[], clipboard:[], pageErrors:[],
               frameFn:null, storage:new FakeStorage(storage, refuse), theme:CSS_VARS, themeListeners:[]};
  const window = {devicePixelRatio:dpr, addEventListener() {}, matchMedia:() => ({addEventListener:(type, fn) => out.themeListeners.push(fn)})};
  const navigator = {clipboard:{writeText:async text => { out.clipboard.push(text); }}};
  out.navigator = navigator;
  const quiet = {log() {}, info() {}, warn:m => out.warnings.push(String(m)), error:m => out.warnings.push(String(m))};
  let script = PAGE_SCRIPT;
  if (mutate) {
    assert.equal(script.split(mutate[0]).length, 2, `the page script no longer holds "${mutate[0]}" once`);
    script = script.replace(mutate[0], mutate[1]);
  }
  const names = ["document", "window", "location", "localStorage", "navigator", "getComputedStyle",
                 "requestAnimationFrame", "alert", "prompt", "Path2D", "Worker", "MessageChannel", "console"];
  const run = new Function(...names, `${scriptBlock("engineSrc")}\n${scriptBlock("workerSrc")}\n${inject}\n${script}${PAGE_EXPORTS}`);
  out.lab = run(document, window, {search}, out.storage, navigator,
    () => ({getPropertyValue:name => out.theme[name] ?? ""}), fn => { out.frameFn = fn; return 1; },
    m => out.alerts.push(String(m)), () => null, FakePath2D, undefined, guardedChannel(out), quiet);
  out.$ = id => document.getElementById(id);
  out.tick = () => { const fn = out.frameFn; out.frameFn = null; if (fn) fn(performance.now()); };
  // run animation frames until `done`, letting the worker's replies in between
  out.settle = async (done, what, ms = 60000) => {
    const start = Date.now();
    while (!done()) {
      if (Date.now() - start > ms) throw new Error(`the page: timed out waiting for ${what}`);
      out.tick();
      await new Promise(resolve => setTimeout(resolve, 1));
    }
  };
  out.idle = () => { const {pend, queued} = out.lab; return !pend.stepping && pend.jumps === 0 && !pend.frame && queued === 0; };
  out.close = () => { if (out.lab.host) out.lab.host.terminate(); };
  out.setTheme = vars => { out.theme = vars; for (const fn of out.themeListeners) fn({}); };
  return out;
}
const sleep = ms => new Promise(resolve => setTimeout(resolve, ms));

// Run `check(page)` on a freshly opened page: null when it passes, otherwise what failed
// (a page that throws while it starts fails too). A broken page's replacement text must
// be in the page script once; if not, this throws, so a check cannot pass a broken page
// that was never built.
async function onPage(glue, options, check) {
  if (options.mutate)
    assert.equal(PAGE_SCRIPT.split(options.mutate[0]).length, 2, `the page script no longer holds "${options.mutate[0]}" once`);
  let page = null;
  try { page = openLab(glue, options); await check(page); return null; }
  catch (err) { return err.message; }
  finally { if (page) page.close(); }
}

// The JavaScript engine run independently with the page's default settings to step t, and
// its chart samples (a record per step, as the worker makes them).
function defaultRun(t) {
  const sim = new LGCA({model:"boson", W:120, H:120, dens:0.4, seed:12345, kernel:"sum", alpha:0.5,
                        bosonField:"site", bosonAlign:"nematic"});
  const samples = [sim.measure()];
  for (let q = 0; q < t; q++) { sim.step(2); samples.push(sim.measure()); }
  return {sim, samples};
}

async function verifyPage(glue) {
  const results = [];

  // (a) labEngine("js") while the module loads, at the start: the page restarts on the
  // JavaScript engine; the first run's replies, which arrive later from the module, are
  // dropped; the lattice, generator and chart are those of an independent run. A page
  // that applied replies of an earlier run would chart the steps twice.
  const switchCheck = async mutate => {
    return onPage(glue, {mutate}, async page => {
      page.lab.labEngine("js");
      page.lab.queueSteps(20);
      await page.settle(() => page.lab.sim.t === 21 && page.idle() && page.lab.sim.frameT === 21, "21 steps");
      await sleep(30);
      const {sim, hist} = page.lab;
      const first = page.received.find(r => r.type === "init" && r.epoch === 1);
      assert(first && first.engine === "wasm", "the first run's init was not answered on the WebAssembly engine");
      assert.equal(page.lab.labEngine(), "js", "the engine after labEngine(\"js\") during the load");
      const ref = defaultRun(21);
      assert.deepEqual(sim.occ, ref.sim.occ, "the page's lattice differs from an independent run");
      assert.deepEqual(sim.frameRng, ref.sim.rng.getState(), "the page's generator state differs from an independent run");
      assert.deepEqual(hist.t, [...Array(22).keys()], "the chart's steps are not 0..21 once each");
      assert.deepEqual(hist.p, ref.samples.map(m => m.polar), "the chart's polar order differs from an independent run");
    });
  };
  assert.equal(await switchCheck(null), null, "labEngine(\"js\") while the module loads");
  assert.notEqual(await switchCheck(["if(r.epoch!==epoch) return;", ""]), null,
    "a page that keeps replies from before its last reset passes the engine-switch check");
  results.push("PASS the page: labEngine(\"js\") while the module loads restarts on the JavaScript engine; lattice, " +
    "generator and chart equal an independent run (a page keeping an earlier run's replies fails)");

  // (b) A rewind whose target leaves the buffer before the worker gets it (a batch in
  // flight evicts it), with steps queued after it (+N while replaying): the steps must not
  // run on from the live state; play stops and the label says the limit was reached.
  const rewindCheck = async mutate => {
    return onPage(glue, {mutate}, async page => {
      await page.settle(() => page.lab.sim.engine !== "", "the init reply");
      page.lab.$("bPlay").click();                                  // pause
      page.lab.queueSteps(80);                                      // 120^2: 69 buffered ticks
      await page.settle(() => page.lab.sim.t === 81 && page.idle(), "81 steps");
      const [oldest] = page.lab.sim.range;
      page.lab.queueSteps(5);                                       // in flight: evicts the oldest ticks
      assert(page.lab.pend.stepping, "the batch is not in flight");
      page.lab.enterReplay(oldest);
      page.lab.fwdN(25);                                            // rewind there, then 25 steps
      await page.settle(() => page.received.some(r => r.type === "jump"), "the rewind's reply");
      const tAtJump = page.lab.sim.t;           // 81 and the steps the batch in flight took (its time budget decides)
      for (let q = 0; q < 20; q++) { page.tick(); await sleep(2); }
      const jumpReply = page.received.find(r => r.type === "jump");
      const after = page.posted.slice(page.posted.findIndex(p => p.type === "jump") + 1);
      assert.equal(jumpReply.ok, false, "the rewind did not fail (its target was still buffered)");
      assert(!after.some(p => p.type === "step"), "steps were sent after the failed rewind: " + JSON.stringify(after));
      assert(tAtJump > 81 && page.lab.sim.t === tAtJump, `the simulation ran on after the failed rewind (${tAtJump} → ${page.lab.sim.t})`);
      assert(page.lab.queued === 0 && !page.lab.playing, "steps still queued, or playing, after the failed rewind");
      assert.equal(page.lab.$("lblBuf").textContent, "at history limit", "the buffer label after the failed rewind");
    });
  };
  assert.equal(await rewindCheck(null), null, "a rewind whose target left the buffer");
  assert.notEqual(await rewindCheck(["clearQueue(); playing=false; $(\"bPlay\").textContent=\"Play\";\n      exitReplay(); flashBufLimit();",
                                     "exitReplay(); flashBufLimit();"]), null,
    "a page that keeps the queued steps after a failed rewind passes");
  results.push("PASS the page: a rewind whose target left the buffer runs no queued steps from the live state and " +
    "stops play (the page without that fix fails)");

  // (c) The history through repeated thinning, and the CSV made from it. Records as the
  // worker sends them, 480,000 steps: band every 10 steps; the sensitivity 2, then 3 at
  // step 3 only, 2 again, 5 from step 250,001; three density diagnostics in turn. Every
  // band sample kept must be in the CSV with its full row; every change of settings must
  // survive, so the settings of every step can be read from the last row at or before it;
  // the series stay within their limits; the chart draws about one point per pixel; the
  // run's identity is on every row; kmax is the busiest channel.
  const historyCheck = async mutate => {
    return onPage(glue, {mutate}, async page => {
      await page.settle(() => page.lab.sim.engine !== "", "the init reply");
      page.lab.$("bPlay").click();
      await page.settle(() => page.idle(), "a quiet page");
      const {lab} = page;
      const cfg = sens => ({sens, alpha:0, kernel:"sum", field:"site", alignment:"nematic"});
      const cfgs = {2:cfg(2), 3:cfg(3), 5:cfg(5)};
      const sensAt = t => t === 3 ? 3 : t > 250000 ? 5 : 2;
      const rhos = [0, 1, 2].map(q => ({N:5760, max:7 + q, maxChannel:3 + q, occupied:4000 + q, mean:0.4, fraction:0.001, effectiveArea:3000 + q}));
      lab.hist.t.length = lab.hist.p.length = lab.hist.n.length = lab.hist.s.length = lab.hist.rho.length = lab.hist.cfg.length = 0;
      lab.hist.bt.length = lab.hist.bv.length = 0;
      const T = 480000;
      let primaryThinned = 0, bandThinned = 0;
      for (let t0 = 0; t0 <= T; t0 += 1000) {
        const recs = [];
        for (let t = t0; t < Math.min(T + 1, t0 + 1000); t++) {
          const rec = {t, polar:t / T, nematic:0.5, spatial:0.25, rho:rhos[t % 3], cfg:cfgs[sensAt(t)]};
          if (t % 10 === 0) rec.band = t / 1e6;
          recs.push(rec);
        }
        const before = [lab.hist.t.length, lab.hist.bt.length];
        lab.applyRecords(recs);
        if (lab.hist.t.length < before[0]) primaryThinned++;
        if (lab.hist.bt.length < before[1]) bandThinned++;
      }
      assert(primaryThinned >= 2 && bandThinned >= 2, `thinned ${primaryThinned} and ${bandThinned} times, not at least twice each`);
      assert(lab.hist.t.length <= lab.HIST_LIMIT && lab.hist.bt.length <= lab.BAND_LIMIT,
        `the series hold ${lab.hist.t.length} and ${lab.hist.bt.length} samples, past their limits`);
      lab.$("bCsv").click();
      await sleep(5);
      const csv = page.clipboard.at(-1).trim().split("\n");
      const header = csv.shift().split(",");
      assert.deepEqual(header, ["step", "polar", "nematic", "spatial", "band", "N", "mean_count", "nmax", "kmax",
        "peak_fraction", "effective_area", "sensitivity", "alpha", "kernel", "field", "alignment",
        "model", "W", "H", "density", "seed", "band_definition"], "the CSV header");
      const rows = csv.map(line => line.split(","));
      assert(rows.every(r => r.length === header.length), "a CSV row has the wrong number of fields");
      const col = name => header.indexOf(name);
      assert(rows.every(r => r.slice(col("model")).join(",") === "boson,120,120,0.4,12345,all-six-channels"),
        "the run's identity is not on every row");
      const banded = rows.filter(r => r[col("band")] !== "");
      assert.deepEqual(banded.map(r => +r[col("step")]), lab.hist.bt, "the CSV's band steps differ from the band samples kept");
      assert.deepEqual(banded.map(r => +r[col("band")]), lab.hist.bv, "the CSV's band values differ from the band samples kept");
      assert(banded.every(r => r[col("polar")] !== "" && r[col("sensitivity")] !== ""), "a band value without its full row");
      for (const r of rows) assert.equal(+r[col("kmax")], 3 + (+r[col("step")] % 3), `kmax at step ${r[col("step")]}`);
      assert.equal(+rows.at(-1)[col("step")], T, "the newest step is not the CSV's last row");
      for (let t = 0, q = 0; t <= T; t++) {          // the settings of every step, from the last row at or before it
        while (q + 1 < rows.length && +rows[q + 1][col("step")] <= t) q++;
        if (+rows[q][col("sensitivity")] !== sensAt(t)) throw new Error(`the CSV loses the sensitivity of step ${t}`);
      }
      // a band sample whose step has no primary row (thinning keeps those rows; this is the
      // guard) still gets a row: its band, empty measurements and settings, the run's identity
      lab.hist.bt.push(T + 5); lab.hist.bv.push(0.77);
      lab.$("bCsv").click();
      await sleep(5);
      const last = page.clipboard.at(-1).trim().split("\n").at(-1).split(",");
      assert.deepEqual(last, ["480005", "", "", "", "0.77", ...Array(11).fill(""), "boson", "120", "120", "0.4", "12345", "all-six-channels"],
        "a band sample without its primary row");
      lab.hist.bt.pop(); lab.hist.bv.pop();
      const ctx = lab.ccv.getContext("2d");
      ctx.calls = []; lab.drawChart();
      const strokes = ctx.calls.filter(c => c[0] === "lineTo").length;
      ctx.calls = null;
      assert(strokes < 6 * 420, `the chart drew ${strokes} line segments for ${lab.hist.t.length + lab.hist.bt.length} samples`);
    });
  };
  assert.equal(await historyCheck(null), null, "the history through repeated thinning and its CSV");
  for (const [what, mutate] of [
    ["thins without keeping band steps and changes of settings",
     ["const keep=all.length<=room?all:all.slice(all.length-room);", "const keep=[...Array(N).keys()].filter(q=>q%2===0);"]],
    ["never thins band", ["if(hist.bt.length>BAND_LIMIT) thinBand();", ""]],
    ["draws every band sample", ["const step=Math.max(1,Math.floor(nb/iw));", "const step=1;"]],
    ["exports the busiest site as kmax", ["d.max,d.maxChannel,", "d.max,d.max,"]],
  ])
    assert.notEqual(await historyCheck(mutate), null, `a page that ${what} passes the history check`);
  results.push("PASS the page: 480,000 steps of history through repeated thinning (both series at least twice): the CSV " +
    "keeps every retained band value on its full row and every change of settings, the run's identity and kmax; " +
    "memory and chart drawing stay bounded (four broken pages fail)");

  // (c2) Settings that change on every step (sensitivity 2 and 3 in turn) for 330,000
  // steps: the changes alone would fill the history, so its oldest part goes instead. What
  // stays must still give the settings of every step from its first row on, and no band
  // sample may be older than that row; the history stays within its limit.
  const denseCheck = mutate => onPage(glue, {mutate}, async page => {
    await page.settle(() => page.lab.sim.engine !== "", "the init reply");
    page.lab.$("bPlay").click();
    await page.settle(() => page.idle(), "a quiet page");
    const {lab} = page, T = 330000;
    const cfgs = [2, 3].map(sens => ({sens, alpha:0, kernel:"sum", field:"site", alignment:"nematic"}));
    const rho = {N:5760, max:7, maxChannel:3, occupied:4000, mean:0.4, fraction:0.001, effectiveArea:3000};
    for (const k of ["t", "p", "n", "s", "rho", "cfg", "bt", "bv"]) lab.hist[k].length = 0;
    let thinned = 0;
    for (let t0 = 0; t0 <= T; t0 += 1000) {
      const recs = [];
      for (let t = t0; t < Math.min(T + 1, t0 + 1000); t++)
        recs.push({t, polar:0.5, nematic:0.5, spatial:0.25, rho, cfg:cfgs[t % 2], ...(t % 10 ? {} : {band:t / 1e6})});
      const before = lab.hist.t.length;
      lab.applyRecords(recs);
      if (lab.hist.t.length < before) thinned++;
    }
    assert(thinned >= 2, `the history thinned ${thinned} times`);
    assert(lab.hist.t.length <= lab.HIST_LIMIT, `the history holds ${lab.hist.t.length} samples`);
    lab.$("bCsv").click();
    await sleep(5);
    const rows = page.clipboard.at(-1).trim().split("\n").slice(1).map(line => line.split(","));
    const first = +rows[0][0];
    assert(first > 0 && +rows.at(-1)[0] === T, `the CSV runs from step ${first} to ${rows.at(-1)[0]}`);
    assert(rows.every(r => r[1] !== ""), "a band value without its full row (band older than the history kept)");
    for (let t = first, q = 0; t <= T; t++) {
      while (q + 1 < rows.length && +rows[q + 1][0] <= t) q++;
      if (+rows[q][11] !== 2 + t % 2) throw new Error(`the CSV loses the sensitivity of step ${t}`);
    }
  });
  assert.equal(await denseCheck(null), null, "settings changed on every step through repeated thinning");
  for (const [what, mutate] of [
    ["keeps every other sample when the changes fill the history",
     ["const keep=all.length<=room?all:all.slice(all.length-room);", "const keep=all.length<=room?all:all.filter((q,i)=>i%2===0);"]],
    ["keeps band samples older than the history it keeps",
     ["  if(keep!==all){\n    let b=0;", "  if(false){\n    let b=0;"]],
  ])
    assert.notEqual(await denseCheck(mutate), null, `a page that ${what} passes the dense-changes check`);
  results.push("PASS the page: settings changing on every step for 330,000 steps: the oldest history goes, band included, " +
    "and the rest keeps the settings of every step (two broken pages fail)");

  // (d) The smaller fixes, each with the broken page it must catch.
  const small = [
    ["the speed label follows a restored slider at start",
     {before:document => { document.getElementById("inSpeed").value = "6"; }},
     async page => assert.equal(page.$("lblSpeed").textContent, "20 /s"),         // SPEEDS[6]
     ["updateSpeedLabel();\nupdateAlpha();", "updateAlpha();"]],
    ["typed sensitivity ignores negative and non-numbers, keeps values above the slider", {},
     async page => {
       const box = page.$("numSens");
       for (const v of ["-1", "abc", "1e400", ""]) { box.value = v; box.dispatch("input"); assert.equal(page.lab.P.sens, 2, `after "${v}"`); }
       box.value = "40"; box.dispatch("input"); assert.equal(page.lab.P.sens, 40);
       box.value = "-3"; box.dispatch("input"); box.dispatch("change");
       assert(page.lab.P.sens === 40 && box.value === "40", "the box does not show the sensitivity in use");
     },
     ["if(!(Number.isFinite(v)&&v>=0)) return;", ""]],
    ["band shows – while its computation is off", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       page.lab.sim._b = 0.321; page.$("inBand").checked = false; page.$("inBand").dispatch("change"); page.tick();
       assert.equal(page.$("vBan").textContent, "–");
     },
     ["const b=!$(\"inBand\").checked?null:shownReplay?replay.slotBand:sim._b;", "const b=shownReplay?replay.slotBand:sim._b;"]],
    ["boson density has no upper limit (fermion keeps its own)", {},
     async page => {
       const seg = page.$("segModel"), dens = page.$("inDens");
       assert(!dens.hasAttribute("max"), "boson density has a maximum at start");
       seg.querySelector('[data-v="fermion"]').click(); assert.equal(dens.max, "0.95");
       seg.querySelector('[data-v="boson"]').click(); assert(!dens.hasAttribute("max"), "boson density has a maximum after the model switch");
       page.$("presetRow").children.at(-1).click(); assert(!dens.hasAttribute("max"), "boson density has a maximum after a recipe");
     },
     ["else $(\"inDens\").removeAttribute(\"max\");            // boson density has no upper limit", "else $(\"inDens\").max=6;"]],
    ["damaged card preferences are ignored",
     {storage:{lgcaCards:JSON.stringify({order:["chart", "view\"]", 5, "nope"], open:{chart:false, view:"yes", model:true}})}},
     async page => {
       const col = page.document.querySelector(".col");
       assert.deepEqual(col.children.map(c => c.dataset.card), ["transport", "model", "view", "chart"]);
       assert(!col.children[3].open && col.children[2].open, "open states from a damaged preference");
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
     },
     ["const d=typeof k===\"string\"?cards.get(k):null;", "const d=colEl.querySelector(`details[data-card=\"${k}\"]`);"]],
    ["a preference that is not an object is ignored", {storage:{lgcaCards:'{"order":[],"open":true}'}},
     async page => { await page.settle(() => page.lab.sim.engine !== "", "the init reply"); },
     ["if(s.open&&typeof s.open===\"object\"&&!Array.isArray(s.open))\n    for(const [k,d] of cards)\n      if(Object.prototype.hasOwnProperty.call(s.open,k)&&typeof s.open[k]===\"boolean\") d.open=s.open[k];",
      "if(s.open)\n    for(const [k,d] of cards) if(k in s.open) d.open=!!s.open[k];"]],
    ["a recipe saved where storage refuses it says so", {refuse:true},
     async page => {
       page.$("bSnapSave").click();
       assert.equal(page.$("bSnapSave").textContent, "saved for this session");
       assert.equal(page.$("snapStore").style.display, "", "the session-only note is hidden");
       assert.equal(page.lab.snaps.length, 1);
     },
     ["try{ localStorage.setItem(SNAP_KEY,JSON.stringify(snaps)); }catch(e){ ok=false; }", "try{ localStorage.setItem(SNAP_KEY,JSON.stringify(snaps)); }catch(e){}"]],
    ["a recipe saved where storage works says so", {},
     async page => {
       page.$("bSnapSave").click();
       assert.equal(page.$("bSnapSave").textContent, "saved ✓");
       assert.equal(page.$("snapStore").style.display, "none");
       assert.equal(JSON.parse(page.storage.getItem("lgcaSettingsSnapshots")).length, 1);
     }, null],
    ["one recipe, as copy JSON gives it, imports", {},
     async page => {
       const recipe = {model:"fermion", sens:0.8, dens:0.2, size:60, seed:7, kernel:"sum", alpha:0.5, bosonField:"site",
                       bosonAlign:"nematic", name:"shared", desc:"", when:"2026-10-08"};
       page.$("snapFile").files = [{text:async () => JSON.stringify(recipe)}];
       page.$("snapFile").dispatch("change");
       await sleep(5);
       assert.equal(page.lab.snaps.length, 1); assert.equal(page.lab.snaps[0].name, "shared");
       assert.equal(page.$("bSnapImp").textContent, "imported 1");
       assert.equal(page.alerts.length, 0, page.alerts.join(" "));
     },
     ["const list=Array.isArray(a)?a:(a&&typeof a===\"object\")?[a]:null;", "const list=Array.isArray(a)?a:null;"]],
    ["a named preset sets the settings its behaviour was measured with, also after an α preset", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       const {lab} = page, button = name => page.$("presetRow").children.find(b => b.textContent === name);
       page.$("inSeed").value = "777";
       button("polar neighbours · α 0.5").click();
       button("stripes").click();
       const now = () => ({model:lab.P.model, size:+page.$("inSize").value, sens:lab.P.sens, dens:+page.$("inDens").value,
                           seed:+page.$("inSeed").value, kernel:lab.P.kernel, field:lab.P.bosonField, align:lab.P.bosonAlign});
       assert.deepEqual(now(), {model:"boson", size:120, sens:2, dens:0.4, seed:12345, kernel:"sum", field:"site", align:"nematic"});
       assert.equal(lab.sim.W, 120, "the run after the preset");
       page.$("segModel").querySelector('[data-v="fermion"]').click();
       page.$("segKernel").querySelector('[data-v="power"]').click();
       button("flocking").click();
       assert.deepEqual(now(), {model:"fermion", size:120, sens:0.8, dens:0.2, seed:12345, kernel:"sum", field:"site", align:"nematic"});
       for (const model of ["boson", "fermion"]) {
         page.$("segModel").querySelector(`[data-v="${model}"]`).click();
         const named = page.$("presetRow").children.filter(b => !b.textContent.startsWith("polar neighbours"));
         assert.equal(named.length, 4, `${model} presets`);
         assert(named.every(b => /^Sensitivity [\d.]+, density [\d.]+ on 120 × 120 with the C\+\+ physics/.test(b.title)),
           `a ${model} preset tooltip: ${named.map(b => b.title).join(" | ")}`);
       }
     },
     ["applySnap({model:P.model,size:120,dens:+dd,sens:+ps,seed:(+$(\"inSeed\").value)>>>0,kernel:\"sum\",alpha:P.alpha,\n                 bosonField:\"site\",bosonAlign:\"nematic\"});",
      "P.sens=+ps; $(\"inSens\").value=ps; $(\"numSens\").value=ps; $(\"inDens\").value=dd; reset();"]],
    ["the inspector follows new frames under a resting pointer; isotropic sites have no axis; net motion in the text colour", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       page.$("bPlay").click();
       await page.settle(() => page.idle() && page.lab.sim.frameT === page.lab.sim.t, "a quiet page");
       const {lab} = page, cvs = page.$("simCanvas");
       const site = (j, i) => (j * 120 + i) * 6;
       // the pointer over site (0, 0): world (0.5 s0, 0.5 s0 sqrt(3)/2), s0 = 5.5
       const at = {clientX:0.5 * 5.5, clientY:0.5 * 5.5 * Math.sqrt(3) / 2};
       lab.sim.occ.set([1, 1, 1, 1, 1, 1], site(0, 0));
       cvs.dispatch("mousemove", at);
       const tip = page.$("siteTip");
       assert(tip.innerHTML.startsWith("site (0, 0) — <b>6</b>"), tip.innerHTML);
       assert(tip.innerHTML.includes("axis <b>–</b>"), `an isotropic site shows an axis: ${tip.innerHTML}`);
       lab.sim.occ.set([9, 0, 0, 0, 0, 0], site(0, 0));
       const rc = page.rose.getContext("2d"); rc.calls = [];
       cvs.dispatch("mousemove", at);
       assert(rc.calls.some(c => c[0] === "=strokeStyle" && c[1] === CSS_VARS["--text-primary"]), "net motion not in the text colour");
       rc.calls = null;
       assert(tip.innerHTML.startsWith("site (0, 0) — <b>9</b>"), tip.innerHTML);
       lab.queueSteps(1);                                            // a new frame arrives; the pointer rests
       await page.settle(() => page.idle() && lab.sim.frameT === lab.sim.t, "the new frame");
       page.tick();
       const now = Array.from(lab.sim.occ.subarray(site(0, 0), site(0, 0) + 6));
       assert(tip.innerHTML.includes(`channels [${now.join(", ")}]`), `the inspector still shows the old frame: ${tip.innerHTML}`);
     },
     ["if(hoverAt&&!panning&&inspectKey()!==hoverKey) inspect(hoverAt.x,hoverAt.y);", ""]],
    ["the inspector redraws in the new colours when the system theme changes", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       page.$("bPlay").click();
       await page.settle(() => page.idle(), "a quiet page");
       page.lab.sim.occ.set([9, 0, 0, 0, 0, 0], 0);                  // site (0, 0): net motion
       page.$("simCanvas").dispatch("mousemove", {clientX:0.5 * 5.5, clientY:0.5 * 5.5 * Math.sqrt(3) / 2});
       const rc = page.rose.getContext("2d"); rc.calls = [];
       page.setTheme({...CSS_VARS, ...DARK_VARS}); page.tick();         // the system turns dark; the pointer rests
       const stroke = rc.calls.some(c => c[0] === "=strokeStyle" && c[1] === DARK_VARS["--text-primary"]);
       rc.calls = null;
       assert(stroke, "the inspector kept the old theme's colour");
     },
     ["addEventListener(\"change\",()=>{ drawChart(); hoverKey=\"\"; });", "addEventListener(\"change\",()=>{ drawChart(); });"]],
    ["the inspector redraws when the particle colours or the view change under a resting pointer", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       page.$("bPlay").click();
       await page.settle(() => page.idle(), "a quiet page");
       page.lab.sim.occ.set([9, 0, 0, 0, 0, 0], 0);                  // site (0, 0): channel 0 only
       page.$("simCanvas").dispatch("mousemove", {clientX:0.5 * 5.5, clientY:0.5 * 5.5 * Math.sqrt(3) / 2});
       const rc = page.rose.getContext("2d");
       const after = (seg, v) => {                                   // a control changes; the pointer rests
         rc.calls = []; page.$(seg).querySelector(`[data-v="${v}"]`).click(); page.tick();
         const strokes = rc.calls.filter(c => c[0] === "=strokeStyle").map(c => c[1]); rc.calls = null;
         return strokes;
       };
       assert(after("segPCol", "uni").includes(CSS_VARS["--accent"]), "uniform colour: the rose kept the channel colours");
       assert(after("segPCol", "dir").includes("#e5484d"), "direction colour: the rose kept the uniform colour");
       assert(after("segMode", "density").includes(CSS_VARS["--accent"]), "density view: the rose kept the channel colours");
     },
     ['+"|"+P.pcol+"|"+P.mode; }', '; }']],
    ["the CSV file is named by the run it holds when the clipboard refuses after a reset", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "" && page.lab.sim.frameT >= 1, "the first frame");
       const {lab} = page, t = lab.sim.t;
       let refuse = null;
       page.navigator.clipboard.writeText = () => new Promise((resolve, reject) => { refuse = reject; });
       page.$("bCsv").click();
       page.$("inSize").value = "60"; page.$("inSeed").value = "777"; lab.reset();
       await page.settle(() => lab.sim.W === 60 && lab.sim.engine !== "", "the 60 x 60 run");
       refuse(new Error("refused on purpose"));
       await sleep(5);
       const csv = page.created.find(e => e.localName === "a" && String(e.download).endsWith(".csv"));
       assert(csv, "no CSV file was saved");
       assert.equal(csv.download, `lgca-boson-120x120-d0.4-seed12345-t${t}.csv`);
     },
     ["    a.download=name;", "    a.download=`lgca-${sim.model}-${sim.W}x${sim.H}-d${sim.dens}-seed${sim.seed}-t${sim.t}.csv`;"]],
    ["guides follow the periodic copies across the edge", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       const {lab} = page;
       lab.$("inSize").value = "60"; lab.reset();
       await page.settle(() => lab.sim.W === 60 && lab.sim.engine !== "" && page.idle(), "the 60 x 60 run");
       lab.view = {x:650, y:0, zoom:4};
       const ctx = lab.ctx; ctx.calls = [];
       lab.drawParticles(new Int32Array(60 * 60 * 6), 0, true);
       const calls = ctx.calls; ctx.calls = null;
       // site (0, 0) is drawn at world x = 0.5 cells = 5.5 px, here only in its copy one world
       // width (660 px) to the right: on screen (665.5 - 650) * 4 = 62, y = 0.5 sqrt(3)/2 * 11 * 4
       const x = 62, y = 0.5 * Math.sqrt(3) / 2 * 11 * 4, near = (a, b) => Math.abs(a - b) < 1e-9;
       assert(calls.some(c => c[0] === "moveTo" && near(c[1], x) && near(c[2], y)), "no guide line from the wrapped site");
       assert(calls.some(c => c[0] === "fillRect" && near(c[1], x - 1) && near(c[2], y - 1)), "no guide dot at the wrapped site");
     },
     ["        const bwx=bwx0+mshift*worldW;", "        const bwx=bwx0+mshift*worldW; if(mshift) continue;"]],
    ["counts follow the periodic copies across the edge", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       const {lab} = page;
       lab.$("inSize").value = "60"; lab.reset();
       await page.settle(() => lab.sim.W === 60 && lab.sim.engine !== "" && page.idle(), "the 60 x 60 run");
       const arr = new Int32Array(60 * 60 * 6); arr[0] = 100;            // 100 particles in channel 0 at site (0, 0)
       lab.view = {x:650, y:0, zoom:4};
       const ctx = lab.ctx; ctx.calls = [];
       lab.drawParticles(arr, 0, true);
       const labels = ctx.calls.filter(c => c[0] === "fillText" && c[1] === 100 && c[2] >= 0 && c[2] <= 660);
       ctx.calls = null;
       assert.equal(labels.length, 1, "the wrapped site's count label");
     },
     ["const sx=(wx+mshift*worldW-view.x)*z, sy=(wy-view.y)*z;", "const sx=(wx+(mshift?1e9:0)-view.x)*z, sy=(wy-view.y)*z;"]],
    ["overlay arrows sit at the centre of partial blocks", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       const {lab} = page;
       lab.$("inSize").value = "60"; lab.reset();
       await page.settle(() => lab.sim.W === 60 && lab.sim.engine !== "" && page.idle(), "the 60 x 60 run");
       lab.$("inBlock").value = "8";
       const arr = new Int32Array(60 * 60 * 6); arr[(58 * 60 + 58) * 6] = 10;   // the last block: sites 56..59
       lab.view = {x:0, y:0, zoom:1};
       const ctx = lab.ctx; ctx.calls = [];
       lab.drawOverlay("flux", arr);
       const starts = ctx.calls.filter(c => c[0] === "moveTo").map(c => [c[1], c[2]]);
       ctx.calls = null;
       const near = (p, x, y) => Math.abs(p[0] - x) < 1e-9 && Math.abs(p[1] - y) < 1e-9;
       // the block of sites 56..59 in both directions: its drawn site centres, (i + j/2 + 1/2) cells
       // across (mod 60) and (j + 1/2) sqrt(3)/2 cells down, average to x = 57.5 + 28.75 + 0.5
       // = 86.75 = 26.75 (mod 60) and y = 58 sqrt(3)/2; 11 px a cell
       assert(starts.some(p => near(p, 26.75 * 11, 58 * Math.sqrt(3) / 2 * 11)), `arrow starts ${JSON.stringify(starts)}`);
     },
     ["const jm=J+Math.min(bsz,H-J)/2, im=I+Math.min(bsz,W-I)/2;", "const jm=J+bsz/2, im=I+bsz/2;"]],
    ["overlay arrows sit at the centroid of their sites", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       const {lab} = page;
       lab.$("inBlock").value = "4";
       const arr = new Int32Array(120 * 120 * 6); arr[(5 * 120 + 9) * 6] = 10;   // the block of sites 8..11, 4..7
       lab.view = {x:0, y:0, zoom:1};
       const ctx = lab.ctx; ctx.calls = [];
       lab.drawOverlay("flux", arr);
       const starts = ctx.calls.filter(c => c[0] === "moveTo").map(c => [c[1], c[2]]);
       ctx.calls = null;
       // centres (i + j/2 + 1/2) average to 9.5 + 2.75 + 0.5 = 12.75 cells; rows (j + 1/2) to 6; 5.5 px a cell
       const near = (p, x, y) => Math.abs(p[0] - x) < 1e-9 && Math.abs(p[1] - y) < 1e-9;
       assert(starts.some(p => near(p, 12.75 * 5.5, 6 * Math.sqrt(3) / 2 * 5.5)), `arrow starts ${JSON.stringify(starts)}`);
     },
     ["const bwx=(((im+jm*0.5-0.25)%W)+W)%W*s0", "const bwx=(((im+jm*0.5)%W)+W)%W*s0"]],
    ["the chart is drawn at the device pixel ratio, and not again while nothing changes", {dpr:2},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       page.$("bPlay").click();
       await page.settle(() => page.idle(), "a quiet page");
       page.tick();
       const {ccv} = page.lab;
       assert.deepEqual([ccv.width, ccv.height], [840, 372], "the chart's backing store at device pixel ratio 2");
       const c = ccv.getContext("2d"); c.calls = [];
       page.tick(); page.tick();
       const redraws = c.calls.filter(k => k[0] === "clearRect").length;
       c.calls = null;
       assert.equal(redraws, 0, "the chart was redrawn with nothing changed");
     },
     ["if(chartState()!==chartKey) drawChart();", "drawChart();"]],
    ["the chart's backing store follows the device pixel ratio", {dpr:2},
     async page => { page.lab.drawChart(); assert.deepEqual([page.lab.ccv.width, page.lab.ccv.height], [840, 372]); },
     ["const w=ccv.clientWidth||420, h=186, bw=Math.round(w*DPR), bh=Math.round(h*DPR);",
      "const w=ccv.clientWidth||420, h=186, bw=w, bh=h;"]],
    ["a paused field view is not rebuilt every frame", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       page.$("bPlay").click();
       page.$("segMode").querySelector('[data-v="density"]').click();
       await page.settle(() => page.idle() && page.lab.sim.frameT === page.lab.sim.t, "a quiet page");
       page.tick();
       const c = page.lab.off.getContext("2d"); c.calls = [];
       page.tick(); page.tick(); page.tick();
       const builds = c.calls.filter(k => k[0] === "putImageData").length;
       page.$("inRhoMax").value = "50"; page.tick();
       const after = c.calls.filter(k => k[0] === "putImageData").length;
       c.calls = null;
       assert.equal(builds, 0, "the picture was rebuilt with nothing changed");
       assert.equal(after, 1, "the picture was not rebuilt after its ceiling changed");
     },
     ["if(arr!==lastFArr||key!==lastFKey){ drawFields(arr); lastFArr=arr; lastFKey=key; }", "drawFields(arr);"]],
    ["a rewind clamped at the oldest buffered tick shows the history limit", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       page.$("bPlay").click();
       page.lab.queueSteps(80);
       await page.settle(() => page.lab.sim.t === 81 && page.idle(), "81 steps");
       page.lab.back(1000);
       await page.settle(() => page.received.some(r => r.type === "jump"), "the rewind");   // the label shows for 1.4 s
       assert.equal(page.lab.sim.t, page.lab.sim.range[0], "the rewind did not reach the oldest tick");
       assert.equal(page.$("lblBuf").textContent, "at history limit");
     },
     ["if(bufFlashT===null) $(\"lblBuf\").textContent=(b-a)+\" ticks\";", "$(\"lblBuf\").textContent=(b-a)+\" ticks\";"]],
    ["a failure in a run already replaced by Reset is not reported against the new run",
     {inject:"{ const step = WasmSim.prototype.step; WasmSim.prototype.step = function (s) { if (this.initialSeed === 999) throw new Error(\"engine failure on purpose\"); step.call(this, s); }; }"},
     async page => {
       await page.settle(() => page.lab.sim.engine === "wasm", "the init reply");
       page.$("inSeed").value = "999"; page.lab.reset();          // this run's init will fail ...
       page.$("inSeed").value = "12345"; page.lab.reset();        // ... after it has been replaced
       await page.settle(() => page.lab.sim.engine !== "" && page.lab.sim.t === 1, "the new run");
       page.lab.queueSteps(3);
       await page.settle(() => page.idle() && page.lab.sim.t === 4, "steps of the new run");
       assert(!page.lab.sim.failed, "the new run was stopped");
       assert(!page.warnings.some(w => w.includes("the simulation stopped")), `a warning for the replaced run: ${page.warnings.join(" | ")}`);
     },
     ["  if(!sim||msg.epoch!==epoch) return;          // a run already replaced by Reset: nothing to stop\n  console.warn(",
      "  console.warn(\"LGCA Lab: the simulation stopped\");\n  if(!sim||msg.epoch!==epoch) return;\n  console.warn("]],
    ["an engine failure on the main-thread host stops the run and says so",
     {inject:"{ const step = WasmSim.prototype.step; WasmSim.prototype.step = function (s) { if (this.t >= 3) throw new Error(\"engine failure on purpose\"); step.call(this, s); }; }"},
     async page => {
       await page.settle(() => page.lab.sim.engine === "wasm", "the init reply");
       page.$("bPlay").click();
       page.lab.queueSteps(5);
       await page.settle(() => page.warnings.some(w => w.includes("the simulation stopped")), "the failure's warning", 10000);
       await page.settle(() => page.idle(), "a quiet page", 10000);
       assert(!page.lab.playing && page.lab.sim.failed, "the run did not stop");
       const sent = page.posted.length;
       page.$("bStep1").click(); page.$("bPlay").click();
       for (let q = 0; q < 10; q++) { page.tick(); await sleep(2); }
       assert(!page.posted.slice(sent).some(p => p.type === "step"), "steps were sent to the failed run");
       page.$("bReset").click();
       await page.settle(() => page.lab.sim.t === 1 && !page.lab.sim.failed, "the new run");
       assert.deepEqual(page.pageErrors, [], "the failure also surfaced as a page error");
     },
     ["simWorkerMain(ch.port2,wasmGlueText(),hostFailed);", "simWorkerMain(ch.port2,wasmGlueText());"]],
    ["a study link with an empty value uses the default for it", {search:"?study=alpha&alpha=0.9&size=90&dens=&sens=6"},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "" && page.lab.sim.W === 90, "the study run");
       assert.equal(page.lab.sim.dens, 0.4);
       assert.equal(page.posted.filter(p => p.type === "init").length, 1, "a default run was started before the study run");
     },
     ["const num=(k,d)=>{ const v=q.get(k); return v===null||v.trim()===\"\"?d:Number(v); };",
      "const num=(k,d)=>{ const v=q.get(k); return v===null?d:Number(v); };"]],
    ["a study link that cannot be used is reported and the defaults run", {search:"?study=alpha&alpha=0.9&size=100&sens=6"},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       assert(page.warnings.some(w => w.includes("the study link was not used")), "no warning about the unusable link");
       assert(page.lab.sim.W === 120 && page.lab.playing, "the ordinary defaults are not running");
     },
     ["    console.warn(\"LGCA Lab: the study link was not used", "    void(\"LGCA Lab: the study link was not used"]],
    ["the inspector names the site whose drawn centre is nearest the pointer", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "", "the init reply");
       const {lab} = page, cvs = page.$("simCanvas"), tip = page.$("siteTip"), s0 = 5.5, H3 = Math.sqrt(3) / 2;
       let seed = 1; const random = () => (seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 2 ** 32;
       for (const zoom of [1, 3]) {
         lab.view = {x:0, y:0, zoom};
         for (let q = 0; q < 300; q++) {
           const clientX = random() * 600, clientY = random() * 500, wx = clientX / zoom, wy = clientY / zoom;
           let best = null;
           for (let jj = Math.floor(wy / (H3 * s0)) - 2; jj <= Math.floor(wy / (H3 * s0)) + 2; jj++)
             for (let ii = Math.floor(wx / s0 - jj / 2) - 2; ii <= Math.floor(wx / s0 - jj / 2) + 2; ii++) {
               const d = (wx - (ii + jj / 2 + 0.5) * s0) ** 2 + (wy - (jj + 0.5) * H3 * s0) ** 2;
               if (!best || d < best.d) best = {d, i:((ii % 120) + 120) % 120, j:((jj % 120) + 120) % 120};
             }
           cvs.dispatch("mousemove", {clientX, clientY});
           assert(tip.innerHTML.startsWith(`site (${best.i}, ${best.j}) `),
             `pointer (${clientX.toFixed(1)}, ${clientY.toFixed(1)}) at zoom ${zoom}: nearest site (${best.i}, ${best.j}), shown ${tip.innerHTML.slice(0, 16)}`);
         }
       }
     },
     ["for(let row=Math.floor(wy/(H3*s0)-0.5), r2=row+1; row<=r2; row++){",
      "for(let row=Math.round(wy/(H3*s0)-0.5), r2=row; row<=r2; row++){"]],
    ["other fields start from their own settings after the density view", {},
     async page => {
       const seg = page.$("segMode"), fields = () => [page.$("inView").value, page.$("inOver").value, page.$("inSmooth").value];
       const view = v => seg.querySelector(`[data-v="${v}"]`).click();
       view("density"); assert.deepEqual(fields(), ["countlog", "none", "crisp"]);
       view("fields"); assert.deepEqual(fields(), ["dens", "flux", "soft"]);
       page.$("inView").value = "nem"; page.$("inView").dispatch("change");
       view("particles"); view("fields");
       assert.equal(page.$("inView").value, "nem", "a choice made in other fields was not kept");
       view("density"); view("particles"); view("fields");
       assert.deepEqual(fields(), ["dens", "flux", "soft"], "density, particles, then other fields");
       view("density"); page.$("inView").value = "dir"; page.$("inView").dispatch("change"); view("fields");
       assert.equal(page.$("inView").value, "dir", "a choice made in the density view was not kept");
     },
     ["else if(v===\"fields\"&&densitySetFields)", "else if(false)"]],
    ["a typed α applies while typing when it is between 0 and 2", {},
     async page => {
       page.$("segKernel").querySelector('[data-v="power"]').click();
       const box = page.$("numAlpha");
       box.value = "1.25"; box.dispatch("input");
       assert.equal(page.lab.P.alpha, 1.25);
       box.value = "3"; box.dispatch("input");
       assert.equal(page.lab.P.alpha, 1.25, "an α above 2 applied while typing");
       assert(page.$("alphaNote").textContent.startsWith("α = 1.25 · field"), page.$("alphaNote").textContent);
     },
     ["if(a>=0&&a<=2){ P.alpha=a; $(\"inAlpha\").value=a; updateAlpha(); }", ""]],
    ["snapshot files are named by run and record the α in effect", {},
     async page => {
       await page.settle(() => page.lab.sim.engine !== "" && page.lab.sim.frameT >= 1, "the first frame");
       page.$("bCapture").click();
       const png = page.created.find(e => e.localName === "a" && String(e.download).startsWith("lgca-density-"));
       assert.equal(png.download, `lgca-density-boson-120x120-seed12345-alpha0-t${page.lab.sim.frameT}.png`);
       page.$("densitySnaps").children[0].children.find(e => e.localName === "button").click();
       const state = page.created.find(e => e.localName === "a" && String(e.download).startsWith("lgca-state-"));
       assert.equal(state.download, `lgca-state-boson-120x120-seed12345-t${page.lab.sim.frameT}.json`);
       const saved = JSON.parse(await resolveObjectURL(state.href).text());
       assert.equal(saved.settings.alpha, 0, "the α in effect under the sum kernel");
     },
     ["settings={...cur,model:sim.model,size:sim.W,dens:sim.dens,seed:sim.initialSeed,alpha:alphaNow};",
      "settings={...cur,model:sim.model,size:sim.W,dens:sim.dens,seed:sim.initialSeed};"]],
  ];
  for (const [what, options, check, mutate] of small) {
    const attempt = m => onPage(glue, {...options, mutate:m}, check);
    const problem = await attempt(null);
    assert.equal(problem, null, `the page: ${what}: ${problem}`);
    if (mutate) assert.notEqual(await attempt(mutate), null, `the page without the fix passes: ${what}`);
  }
  results.push(`PASS the page: ${small.length} checks of controls and displays (speed label, typed sensitivity and α, ` +
    "band when off, boson density limit, card preferences, recipe storage and import, presets, the site inspector, " +
    "wrapped counts and guides, overlay centres, chart sharpness and redraws, field-view rebuilds, the history-limit label, an " +
    "engine failure on the main thread, study links, other-fields defaults, snapshot names, the inspector after colour and view " +
    "changes, the CSV file name after a late clipboard refusal); each fix's absence fails its check");

  // (e) The study link's defaults.
  {
    const page = openLab(glue, {search:"?study=alpha"});
    try {
      await page.settle(() => page.lab.sim.engine !== "" && page.lab.sim.W === 180, "the study run");
      const {lab} = page;
      assert.deepEqual({model:lab.P.model, size:lab.sim.W, dens:lab.sim.dens, sens:lab.P.sens, seed:lab.sim.seed,
                        kernel:lab.P.kernel, alpha:lab.P.alpha, field:lab.P.bosonField, align:lab.P.bosonAlign,
                        playing:lab.playing, mode:lab.P.mode},
        {model:"boson", size:180, dens:0.4, sens:6, seed:12345, kernel:"power", alpha:0, field:"neigh", align:"polar",
         playing:false, mode:"density"}, "?study=alpha with no other parameters");
    } finally { page.close(); }
  }
  results.push("PASS the page: ?study=alpha alone starts the documented study setup, paused, in the density view");
  return results;
}

// ------------------------------------------------------------------------------------

const block = verifyBlock();
console.log(`PASS wasmEngine block: ${block.source}, sha256 ${block.hash.slice(0, 16)}..., inlined intact`);
if (block.note) console.log(block.note);

const engine = await loadWasmEngine(block.glue);
verifyReplay(engine);
console.log("PASS exact replay round trip on the WebAssembly engine");
const optionRuns = verifyOptions(engine);
console.log(`PASS research options: ${optionRuns} combinations of kernel, boson field and alignment, ` +
  "and options changed between steps, step for step with the JavaScript engine");
const cdfChecks = verifyCdfBoundaries(engine);
console.log(`PASS boson sampling at ${cdfChecks} forced draws on channel-distribution boundaries ` +
  "(ties and zero-weight channels) picks the reference's channel, on both engines");
const rngRefusals = verifyRngStateContract(o => new WasmSim(engine, o), "WebAssembly") +
  verifyRngStateContract(o => new LGCA(o), "JavaScript");
console.log(`PASS the generator refuses ${rngRefusals} malformed saved states (short, long, index out of range, ` +
  "words not 32-bit integers) on both engines, which stay as they were");
const densityCases = verifyDensityStats();
console.log(`PASS density diagnostics on ${densityCases} hand-worked lattices (empty, one pile, unequal piles, uniform)`);
console.log(`Reference C++: ${reference.source}`);

const tempDir = mkdtempSync(join(tmpdir(), "lgca-parity-wasm-"));
try {
  const binaries = {boson:buildReference("boson", tempDir), fermion:buildReference("fermion", tempDir)};
  for (const scenario of scenarios) {
    const label = `${scenario.model} t=${scenario.tsteps} r=${scenario.iters} s=${scenario.sens} d=${scenario.dens} seed=${scenario.seed}`;
    const {sim, metrics} = runScenario(engine, scenario, label);
    compareWithCpp(sim, metrics, runCpp(binaries[scenario.model], scenario), label);
    console.log(`PASS ${label} (C++ and JavaScript)`);
  }
  for (const scenario of emptyScenarios) {
    const label = `empty lattice: ${scenario.model} t=${scenario.tsteps} r=${scenario.iters} s=${scenario.sens} d=0 seed=${scenario.seed}`;
    const {sim, metrics} = runScenario(engine, scenario, label);
    compareWithCpp(sim, metrics, runCpp(binaries[scenario.model], scenario), label, {empty:true});
    const rho = densityStats(sim.occ);
    assert.deepEqual(rho, {N:0, max:0, maxChannel:0, occupied:0, mean:0, fraction:0, effectiveArea:0},
      `${label}: density diagnostics of the empty lattice`);
    console.log(`PASS ${label}: lattice, generator and the reference's NaN conventions (C++ and JavaScript)`);
  }
} finally {
  rmSync(tempDir, {recursive:true, force:true});
}

for (const line of await verifyProtocol(block.glue)) console.log(line);
for (const line of await verifyPage(block.glue)) console.log(line);
console.log(`\nAll ${scenarios.length} C++/WebAssembly parity scenarios passed.`);
