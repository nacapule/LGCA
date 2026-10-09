"""The islands of every saved lattice: heaviest.json.

    python3 build/heaviest.py

The runs saved their lattices at steps 100, 500, 1000, 2000, 5000, 10000 and 20000, and in
the larger boxes also every 20000 steps after that, to the end of each run (80000, or 160000
for the runs at alpha 0.85 to 0.94). This script finds every island of each saved lattice
and records, as [particles, sites]:

- island: the busiest site's island, the same island as the trace's;
- heaviest: the island with the most particles;
- densest: the island with the most particles per site;
- densest2: the densest island of at least two sites.

An empty lattice has no island: island, heaviest and densest are then [0, 0] (as in the
trace) and densest2 is null.

The 90 x 90 lattices are those of the September sweep (../density-alpha-sweep-2026-09-30/data),
or density/<job>/ for the jobs that sweep did not run to step 20,000. Each lattice must have
the SHA-256 that the trace recorded for that checkpoint, and its busiest-site island,
heaviest mass, island count and N must equal the trace row at the same step; any difference
stops the script.
"""
import gzip
import hashlib
import json
import multiprocessing
import re
import sys

import numpy as np

from islands import DENSITY_DIR, SIZE_DIR, STUDY

SEPT_DATA = STUDY.parent / "density-alpha-sweep-2026-09-30" / "data"
STATE = re.compile(r"state-t(\d+)\.bin\.gz$")
# three of the six neighbours; the other three are their reverses
DIRECTIONS = [(1, 0), (0, 1), (1, -1)]


def islands(n, L):
    """Island label of every site (connected occupied sites, six periodic neighbours)."""
    occupied = n > 0
    i = np.tile(np.arange(L), L)
    j = np.repeat(np.arange(L), L)
    parent = np.arange(L * L)
    us, vs = [], []
    for di, dj in DIRECTIONS:
        v = ((j + dj) % L) * L + (i + di) % L
        both = occupied & occupied[v]
        us.append(np.flatnonzero(both))
        vs.append(v[both])
    u, v = np.concatenate(us), np.concatenate(vs)
    while True:                     # hook roots onto the smaller root, then compress
        pu, pv = parent[u], parent[v]
        differ = pu != pv
        if not differ.any():
            break
        np.minimum.at(parent, np.maximum(pu, pv)[differ], np.minimum(pu, pv)[differ])
        while True:
            jumped = parent[parent]
            if np.array_equal(jumped, parent):
                break
            parent = jumped
    return parent, occupied


def read_lattice(path, L, sha256=None):
    """The particles of each channel of a saved lattice, (L*L, 6). With `sha256`, the
    lattice bytes must have that hash (the one the trace recorded for this checkpoint)."""
    raw = gzip.decompress(path.read_bytes())
    if sha256 is not None and hashlib.sha256(raw).hexdigest() != sha256:
        raise ValueError(f"{path}: SHA-256 differs from the trace's checkpoint")
    if len(raw) != L * L * 6 * 4:
        raise ValueError(f"{path}: {len(raw)} bytes, not a {L} x {L} lattice")
    return np.frombuffer(raw, dtype="<i4").reshape(L * L, 6)


def measure(occ, L):
    n = occ.sum(1).astype(np.int64)
    if not n.any():                           # no occupied site: no island
        return {"N": 0, "components": 0, "island": [0, 0], "heaviest": [0, 0],
                "densest": [0, 0], "densest2": None, "densestIsIsland": False}
    label, occupied = islands(n, L)
    roots, inverse = np.unique(label[occupied], return_inverse=True)
    mass = np.bincount(inverse, weights=n[occupied]).astype(np.int64)
    area = np.bincount(inverse)
    busiest = int(np.argmax(n))
    b = int(np.searchsorted(roots, label[busiest]))
    h = int(np.argmax(mass))                  # first heaviest in root order
    d = int(np.argmax(mass / area))
    multi = area >= 2
    d2 = int(np.flatnonzero(multi)[np.argmax((mass / area)[multi])]) if multi.any() else None
    return {"N": int(n.sum()), "components": len(roots),
            "island": [int(mass[b]), int(area[b])],
            "heaviest": [int(mass[h]), int(area[h])],
            "densest": [int(mass[d]), int(area[d])],
            "densest2": [int(mass[d2]), int(area[d2])] if d2 is not None else None,
            "densestIsIsland": bool(d == b)}


def job(task):
    route, name, trace_path, folders = task
    record = json.loads(trace_path.read_text())
    # the trace's checkpoints after step 0; the lattices must cover every one of them
    wanted = {s["t"] for s in record["snaps"] if 0 < s["t"] <= record["reached"]}
    state_dir = next((f for f in folders
                      if wanted <= {int(STATE.search(p.name).group(1))
                                    for p in f.glob("state-t*.bin.gz")}), None)
    if state_dir is None:
        return None, [f"{name}: no folder holds the lattices of every checkpoint {sorted(wanted)}"]
    L = record["job"]["size"]
    col = {c: k for k, c in enumerate(record["columns"])}
    rows = {row[0]: row for row in record["rows"]}
    sha = {s["t"]: s["sha256"] for s in record["snaps"]}
    out = {"route": route, "job": name, **{k: record["job"][k] for k in ("size", "alpha", "dens", "seed")},
           "steps": {}}
    problems = []
    for path in sorted(state_dir.glob("state-t*.bin.gz"), key=lambda p: int(STATE.search(p.name).group(1))):
        t = int(STATE.search(path.name).group(1))
        if t == 0 or t > record["reached"]:
            continue
        if t not in sha:
            problems.append(f"{name} t={t}: the trace has no checkpoint at this step")
            continue
        try:
            m = measure(read_lattice(path, L, sha[t]), L)
        except ValueError as err:
            problems.append(f"{name} t={t}: {err}")
            continue
        row = rows.get(t)
        if row is None:
            problems.append(f"{name} t={t}: no trace row")
            continue
        expected = (row[col["islandMass"]], row[col["islandArea"]], row[col["largestMass"]],
                    row[col["components"]], row[col["N"]])
        got = (m["island"][0], m["island"][1], m["heaviest"][0], m["components"], m["N"])
        if expected != got:
            problems.append(f"{name} t={t}: trace {expected} != lattice {got}")
        out["steps"][t] = m
    if set(out["steps"]) != wanted:
        problems.append(f"{name}: measured steps {sorted(out['steps'])} != checkpoints {sorted(wanted)}")
    return out, problems


def tasks():
    for trace in sorted(SIZE_DIR.glob("*/trace.json")):
        yield "size", trace.parent.name, trace, [trace.parent]
    for trace in sorted(DENSITY_DIR.glob("*/trace.json")):
        yield "density", trace.parent.name, trace, [trace.parent, SEPT_DATA / trace.parent.name]


def main():
    work = list(tasks())
    print(f"{len(work)} runs", flush=True)
    results, problems = [], []
    with multiprocessing.Pool(4) as pool:
        for k, (out, bad) in enumerate(pool.imap_unordered(job, work, chunksize=4), 1):
            if out is not None:
                results.append(out)
            problems.extend(bad)
            if k % 200 == 0:
                print(f"  {k}/{len(work)}", flush=True)
    if problems:
        print("\n".join(problems[:20]), file=sys.stderr)
        raise SystemExit(f"{len(problems)} lattices disagree with their traces")
    results.sort(key=lambda r: (r["route"], r["job"]))
    checked = sum(len(r["steps"]) for r in results)
    path = STUDY / "heaviest.json"
    path.write_text(json.dumps({"runs": results}, separators=(",", ":")))
    print(f"wrote {path}: {len(results)} runs, {checked} lattices, every one matching its trace "
          "(SHA-256 and island numbers)")


if __name__ == "__main__":
    main()
