# C++ source and tools

The simulation is written in C++. The same engine runs in two places: compiled to
WebAssembly inside [lgca-lab.html](../lgca-viz/lgca-lab.html), where opening the Lab
needs no build tools or installation, and as a command-line program that runs the
studies' jobs and measures their islands.

| File | Contents |
|---|---|
| [lgca_engine.cpp](lgca_engine.cpp) | Initialization, collisions, streaming and measurements |
| [lgca_engine.hpp](lgca_engine.hpp) | Parameters, research settings and engine interface |
| [well1024a.hpp](well1024a.hpp) | WELL1024a random generator |
| [wasm/lgca_wasm.cpp](wasm/lgca_wasm.cpp) | Functions that connect the browser to the C++ engine |
| [build-wasm.sh](build-wasm.sh) | Compiles the engine for the browser |
| [inline-wasm.mjs](inline-wasm.mjs) | Places the compiled engine in the Lab HTML |
| [build.sh](build.sh) | Builds the command-line program |
| [lgca_cli.cpp](lgca_cli.cpp) | The program's command line, including the reference model's arguments |
| [run_cli.cpp](run_cli.cpp), [trace_cli.cpp](trace_cli.cpp), [measure_cli.cpp](measure_cli.cpp) | Its `run`, `trace` and `measure` modes |
| [measure.cpp](measure.cpp), [measure.hpp](measure.hpp) | Observables and islands |
| [sweep_job.cpp](sweep_job.cpp), [sweep_job.hpp](sweep_job.hpp) | A job's fixed settings, name, sample steps and checkpoints |
| [file_io.cpp](file_io.cpp), [json.cpp](json.cpp), [sha256.cpp](sha256.cpp) | Files, JSON and SHA-256, with their headers |
| [sweep.py](sweep.py), [trace_sweep.py](trace_sweep.py) | Run or trace many jobs in parallel |
| [sweeps/](sweeps/) | The island study's job lists |

The [model guide](../docs/PHYSICS.md) describes the collision rules and
observables. The [verification guide](../docs/FIDELITY.md) explains the
reference comparisons and numerical limits.

## Running the program

The command-line program needs a C++17 compiler and zlib (part of macOS; on Linux,
the zlib development package). The job drivers need Python 3.11 or later; the island
study's scripts also need NumPy.

### Build

From the repository root:

```sh
bash native/build.sh
```

This writes `native/build/lgca_native`. The build turns off fused multiply-add
(`-ffp-contract=off`) so that the arithmetic matches the reference; do not add
`-ffast-math` or `-Ofast`.

### Run one job

A job is the studies' setting: polar bosons on an L × L lattice, the field from a site
and its six neighbours, division by M^alpha, sensitivity 6. Only the side L, alpha, the
density per channel and the seed vary, and the job's name holds all four:

```sh
native/build/lgca_native run --job both-L90-a0.85-d0.4-seed12345 --end 20000 \
    runs/both-L90-a0.85-d0.4-seed12345
```

`--size 90 --alpha 0.85 --dens 0.4 --seed 12345` names the same job. L goes from 4 to
10,000, alpha from 0 to 2, the density from 0 to 357,913,941 and the seed from 0 to
4,294,967,295.

The folder receives `run.json`, with the measurements (every 25 steps up to step 5,000,
then every 100) and a checkpoint record, and the lattice at every checkpoint,
`state-t<T>.bin.gz`. The end must be a checkpoint:

```text
0  100  500  1000  2000  5000  10000  20000  40000  60000  80000
then every 20000 up to 2560000
```

The same command with a later end continues the job from its last checkpoint; with an
end already reached it says so and writes nothing. A job run to the last checkpoint has
135 saved lattices and a `run.json` of 11 to 12 MB, whatever the lattice size.

### Run a sweep

A TOML file lists the jobs; each `[[grid]]` stands for every combination of its lists:

```toml
end = 20000

[[grid]]
sizes = [90]
alphas = [0.8, 0.9, 1]
densities = [0.4]
seeds = [12345, 777]
```

```sh
python3 native/sweep.py my-sweep.toml runs --processes 4
```

Each job gets its own folder in `runs/`. `--processes` sets how many jobs run at once
(by default one per processor), and `--dry-run` lists what would run. Ctrl-C stops the
sweep; running the same command again continues every job from its last checkpoint.

### Trace jobs: the islands at every step

`trace` runs a job from its seed and measures its islands at every sample step:

```sh
native/build/lgca_native trace --job both-L90-a0.85-d0.4-seed12345 --end 20000 \
    --states traces/both-L90-a0.85-d0.4-seed12345 \
    traces/both-L90-a0.85-d0.4-seed12345/trace.json
```

`trace.json` holds one row per sample step; `columns` names the values in each row.
With `--states`, the lattice of every checkpoint is saved beside it. With `--check DIR`,
where DIR is a `run` folder of the same job, every saved row and checkpoint of DIR must
come out the same. `trace_sweep.py` traces every job of a config in parallel.

### The island study

The island study's runs come from the five configs in [sweeps/](sweeps/), run in six
stages. Each stage continues the traces that the stages before it left in the same
folder:

| Stage | Config | Runs |
|---|---|---|
| 1 | `islands-size-2026-10-06.toml` with `--end 20000` | Larger boxes, L = 60 to 480, to step 20,000 |
| 2 | `islands-density-2026-10-06.toml` | The denser 90 × 90 box, to step 20,000 |
| 3 | `islands-size-2026-10-06.toml` | Larger boxes continued to step 80,000 |
| 4 | `islands-size-L720-2026-10-06.toml` | L = 720 at four alphas, to step 80,000 |
| 5 | `islands-size-160k-2026-10-08.toml` | Alpha 0.85 to 0.94 continued to step 160,000 (L = 720 at 0.9 and 0.94) |
| 6 | `islands-size-320k-2026-10-09.toml` | Alpha 0.9 and 0.92 continued to step 320,000 (L = 720 at 0.9) |

From the repository root:

```sh
S=native/sweeps
OUT=results/islands-2026-10-06
python3 native/trace_sweep.py $S/islands-size-2026-10-06.toml $OUT/size --no-check --end 20000
python3 native/trace_sweep.py $S/islands-density-2026-10-06.toml $OUT/density --no-check
python3 native/trace_sweep.py $S/islands-size-2026-10-06.toml $OUT/size --no-check
python3 native/trace_sweep.py $S/islands-size-L720-2026-10-06.toml $OUT/size --no-check
python3 native/trace_sweep.py $S/islands-size-160k-2026-10-08.toml $OUT/size --no-check
python3 native/trace_sweep.py $S/islands-size-320k-2026-10-09.toml $OUT/size --no-check
```

Every config saves the lattice at every checkpoint, which the study's scripts read.
Four of them also compare each job that the density and alpha study ran with its saved
runs in `results/density-alpha-sweep-2026-09-30/data`. That folder is not part of this
repository, so `--no-check` leaves the comparison out. A stage run on its own traces its
jobs from the seed, with the same result. The whole study took roughly 80
processor-hours; `--end 500` stops every job at step 500, and the same command without
it continues them.

When the stages have finished, the study's scripts make the page and the CSV files. They
need every denser-box job at step 20,000, and at least four seeds for every box, density
and alpha, all ending at the same step. From `results/islands-2026-10-06`:

```sh
python3 build/heaviest.py      # the islands of every saved lattice
python3 build/page_data.py     # seed means and late values
python3 build/make_page.py     # island-study.html
python3 build/export_csv.py    # the three CSV files
```

### Other modes

Given the reference program's arguments (steps, realizations, sensitivity, density,
seed), the program runs the reference model and prints the reference's output line.
The band column counts all six channels ([band correction](../docs/FIDELITY.md#the-band-correction)).
Options choose the lattice size, the fermion model and the research settings of the
[model guide](../docs/PHYSICS.md); the top of [lgca_cli.cpp](lgca_cli.cpp) lists them.

```sh
native/build/lgca_native 500 4 6 0.4 12345
native/build/lgca_native 20000 1 6 0.4 12345 --W 90 --H 90 \
    --kernel power --alpha 0.9 --boson-field both --boson-align polar
```

`lgca_native measure` measures saved lattices and prints one line of JSON per file;
with `--island` it gives the island and the heaviest and densest islands
([measure_cli.cpp](measure_cli.cpp)).

### Where the outputs go

| Command | Writes |
|---|---|
| `lgca_native run ... DIR` | `DIR/run.json`, `DIR/state-t<T>.bin.gz` |
| `sweep.py CONFIG DIR` | one such folder per job, `DIR/<job>/` |
| `lgca_native trace ... OUT` | `OUT`; with `--states DIR`, `DIR/state-t<T>.bin.gz` |
| `trace_sweep.py CONFIG OUT` | `OUT/<job>/trace.json` and, with `states = true`, the job's lattices beside it |

A job's folder or trace can be extended later with a later end. Stopping any of these
programs is safe: each file is replaced only by a complete new version.

## Rebuild the Lab engine

For readers changing the C++ source, rebuilding needs Emscripten (`em++`
on PATH), Node.js and a C++ compiler for the reference checks.
The included engine was built with Emscripten 6.0.10; other compiler
versions can produce different module bytes.

From the repo root:

```sh
bash native/build-wasm.sh
node native/inline-wasm.mjs
node lgca-viz/verify-parity-wasm.mjs
node lgca-viz/verify-parity.mjs
```

The build produces `native/build/lgca_wasm.js`, which contains the compiled
engine and its loader. The second command embeds it in the HTML. The last
two commands check the resulting Lab against the reference model.

To compare a build with the embedded copy without changing the HTML:

```sh
node native/inline-wasm.mjs --check
```
