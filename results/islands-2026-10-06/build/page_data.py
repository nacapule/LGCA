"""The page's data: page-data.json.

    python3 build/page_data.py

Seed means of the island density over time, on one grid of log-spaced sample steps shared
by both boxes (so the two boxes are compared at the same steps), and late values: the mean
over the last quarter of each run, with the smallest and largest seed value. The larger-box
runs go on to step 80,000, and at alpha 0.85 to 0.94 to 160,000. Their late values are also
given at step 20,000, where the denser box ends, and those of the runs to 160,000 also at
step 80,000; the grid stops at step 80,000. At every saved lattice the trace's heaviest and
densest islands must equal heaviest.json (build/heaviest.py); any difference stops the script.

Every curve and late value is a mean over at least four seeds sampled at the same steps;
anything else stops the script. Values are rounded to 4 significant digits; a curve is null
after its runs end.
"""
import json
from datetime import datetime

import numpy as np

from islands import (BOX_DENS, DENSITY_DIR, DENSITY_L, SIZE_DIR, STUDY, check_group, describe,
                     group, in_late_window, late, load)

CURVE_KEYS = ["rho", "rhoOld", "rhoDensest"]    # new island, old cluster, densest island
HEAVIEST = {}                       # (route, job name) -> {step: islands}, filled in main()
LAST_STEP = 80000                   # the grid ends here, also for the runs that go on to 160,000
GRID_POINTS = 180                   # log-spaced points from step 25 to LAST_STEP


def saved_steps(end):
    """The steps at which the runs saved their lattice, up to step `end`: 100, 500, 1,000,
    2,000, 5,000, 10,000, then every 20,000 (the ends of the runs are among them)."""
    return [t for t in [100, 500, 1000, 2000, 5000, 10000, *range(20000, end + 1, 20000)] if t <= end]


def sample_grid():
    """The sample steps the page and the CSV use, the same for every run: step 0, the saved
    steps up to LAST_STEP, and GRID_POINTS log-spaced steps from 25 to LAST_STEP, each rounded
    to a step the runs sample (every 25 steps to 5,000, then every 100). Up to step 20,000
    this keeps 106 steps, and 139 up to LAST_STEP."""
    steps = {0, *saved_steps(LAST_STEP)}
    for t in np.geomspace(25, LAST_STEP, GRID_POINTS):
        spacing = 25 if t <= 5000 else 100
        steps.add(int(round(t / spacing)) * spacing)
    return sorted(steps)


def rounded(values):
    out = []
    for v in values:
        if v is None or not np.isfinite(v):
            out.append(None)
        elif v == 0:
            out.append(0)
        else:
            out.append(float(f"{v:.4g}"))
    return out


def curve(runs, key, steps):
    """Seed mean at each of `steps` up to the runs' common end, null after it. The runs
    must all have a sample at every one of those steps (check_group: same sample steps)."""
    horizon = runs[0]["reached"]
    index = {int(s): k for k, s in enumerate(runs[0]["t"])}
    missing = [s for s in steps if s <= horizon and s not in index]
    if missing:
        raise SystemExit(f"{describe(runs[0])}: no sample at grid steps {missing[:5]}")
    values = np.array([r[key] for r in runs])
    out = []
    for s in steps:
        if s > horizon:
            out.append(None)
            continue
        column = values[:, index[s]]
        if not np.isfinite(column).all():
            raise SystemExit(f"{describe(runs[0])}: undefined {key} at step {s}")
        out.append(column.mean())
    return rounded(out)


def late_summary(runs):
    t = runs[0]["t"]
    window = t[in_late_window(t, runs[0]["reached"])]
    out = {"seeds": len(runs), "N": int(round(np.mean([r["N0"] for r in runs]))),
           "mean": float(f"{np.mean([r['mean'] for r in runs]):.4g}"), "T": runs[0]["reached"],
           "window": [int(window[0]), int(window[-1])]}
    for key in CURVE_KEYS:
        values = [late(r, key, strict=True) for r in runs]
        out[key] = rounded([np.mean(values), np.min(values), np.max(values)])
    return out


def check_saved(route, runs):
    """Stop unless heaviest.json has every run of the route and every saved lattice of each,
    and the run's trace row at each of them has the heaviest island's sites and the densest
    islands that heaviest.json found on that lattice."""
    names = set()
    for r in runs:
        name = f"both-L{r['size']}-a{r['alpha']:g}-d{r['dens']:g}-seed{r['seed']}"
        names.add(name)
        saved = HEAVIEST.get((route, name))
        if not saved:
            raise SystemExit(f"heaviest.json has no lattices for {name}: run build/heaviest.py")
        steps = sorted(int(t) for t in saved if int(t) <= r["reached"])
        if steps != saved_steps(r["reached"]):
            raise SystemExit(f"{route} {describe(r)}: heaviest.json has the lattices of steps {steps}: "
                             "run build/heaviest.py")
        row = {int(t): k for k, t in enumerate(r["t"])}
        for t in steps:
            if t not in row:
                raise SystemExit(f"{route} {describe(r)}: no trace row at saved step {t}")
            m = saved[str(t)]
            densest2 = m["densest2"] or [None, None]
            expected = [m["heaviest"][1], *m["densest"], *densest2]
            got = [r[c][row[t]] for c in ("largestArea", "densestMass", "densestArea",
                                          "densest2Mass", "densest2Area")]
            got = [None if np.isnan(v) else float(v) for v in got]
            if got != expected:
                raise SystemExit(f"{route} {describe(r)} step {t}: trace {got} != heaviest.json {expected}")
    missing = sorted(job for rt, job in HEAVIEST if rt == route and job not in names)
    if missing:
        raise SystemExit(f"{route}: {len(missing)} runs of heaviest.json have no trace, such as {missing[0]}")


def cut(runs, horizon):
    out = []
    for r in runs:
        r = dict(r)
        keep = r["t"] <= horizon
        for key in ["t"] + CURVE_KEYS:
            r[key] = r[key][keep]
        r["reached"] = horizon
        out.append(r)
    return out


def size_route(steps):
    runs = load(SIZE_DIR)
    check_saved("size", runs)
    entries = []
    for alpha in sorted({r["alpha"] for r in runs}):
        for L in sorted({r["size"] for r in runs}):
            g = group(runs, alpha=alpha, size=L)
            if not g:
                continue
            check_group(g, f"larger box L {L} alpha {alpha:g}")   # also: every seed ends at one step
            horizon = g[0]["reached"]
            at20000 = horizon >= 20000   # the same runs at step 20,000, like the denser box
            entry = {
                "alpha": alpha, "L": L, "horizon": horizon,
                "curves": {key: curve(g, key, steps) for key in CURVE_KEYS},
                "late": late_summary(g),
                "late20000": late_summary(cut(g, 20000)) if at20000 else None,
            }
            if horizon > 80000:          # the runs to 160,000 also at 80,000, where the others end
                entry["late80000"] = late_summary(cut(g, 80000))
            entries.append(entry)
    return {"steps": steps, "density": BOX_DENS, "entries": entries}


def density_route(steps):
    runs = load(DENSITY_DIR, 20000)
    check_saved("density", runs)
    entries = []
    for alpha in sorted({r["alpha"] for r in runs}):
        for dens in sorted({r["dens"] for r in runs}):
            g = group(runs, alpha=alpha, dens=dens)
            if g:
                check_group(g, f"denser box dens {dens:g} alpha {alpha:g}")
                entries.append({
                    "alpha": alpha, "dens": dens,
                    "N": int(round(np.mean([r["N0"] for r in g]))),
                    "curves": {key: curve(g, key, steps) for key in CURVE_KEYS},
                    "late": late_summary(g),
                })
    return {"steps": steps, "L": DENSITY_L, "entries": entries}


def main():
    for run in json.loads((STUDY / "heaviest.json").read_text())["runs"]:
        HEAVIEST[(run["route"], run["job"])] = run["steps"]
    grid = sample_grid()
    data = {
        "generated": datetime.now().isoformat(timespec="minutes"),
        "settings": {"model": "boson", "alignment": "polar", "kernel": "power (divide by M^alpha)",
                     "field": "centre + 6 neighbours", "sensitivity": 6,
                     "boundaries": "periodic", "lattice": "hexagonal, 6 velocity channels"},
        "size": size_route(grid),
        "density": density_route([s for s in grid if s <= 20000]),
    }
    path = STUDY / "page-data.json"
    path.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False))
    print("wrote", path, f"{path.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
