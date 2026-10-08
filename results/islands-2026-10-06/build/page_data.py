"""The page's data: page-data.json.

    python3 build/page_data.py

Seed means of the island density over time, on one grid of log-spaced sample steps shared
by both boxes (so the two boxes are compared at the same steps), and late values: the mean
over the last quarter of each run, with the smallest and largest seed value. The larger-box
runs go on to step 80,000 (40,000 at L = 720); their late values are also given at step
20,000, where the denser box ends. The densest island comes from heaviest.json
(build/heaviest.py) and exists only at the saved lattices: its curve has a point at each
saved step, and its late value is the mean over the saved lattices inside the late window
(the one at the end of the run) until the densest island is traced at every step.

Every curve and late value is a mean over at least four seeds sampled at the same steps;
anything else stops the script. Values are rounded to 4 significant digits; a curve is null
after its runs end.
"""
import json
from datetime import datetime

import numpy as np

from islands import (BOX_DENS, DENSITY_DIR, DENSITY_L, SIZE_DIR, STUDY, check_group, describe,
                     group, in_late_window, late, load)

CURVE_KEYS = ["rho", "rhoOld"]      # island density, old cluster density
HEAVIEST = {}                       # (route, job name) -> {step: islands}, filled in main()
# the saved lattices, and the ends of the runs (20,000 in the denser box; 40,000 at L = 720)
SAVED_STEPS = [100, 500, 1000, 2000, 5000, 10000, 20000, 40000, 60000, 80000]
LAST_STEP = 80000
GRID_POINTS = 180                   # log-spaced points from step 25 to LAST_STEP


def sample_grid():
    """The sample steps the page and the CSV use, the same for every run: step 0, the saved
    steps (which include every run's end), and GRID_POINTS log-spaced steps from 25 to
    LAST_STEP, each rounded to a step the runs sample (every 25 steps to 5,000, then every
    100). Up to step 20,000 this keeps 106 steps."""
    steps = {0, *SAVED_STEPS}
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


def island_density(record):
    """Particles per site of a [particles, sites] island record; 0 for no island."""
    mass, sites = record
    return mass / sites if sites else 0.0


def densest_block(route, runs, horizon):
    """The densest island at the saved lattices up to `horizon`: seed means of its density,
    and `late` = [seed mean, min, max] over the saved lattices in the late window."""
    saved_by_run = []
    for r in runs:
        name = f"both-L{r['size']}-a{r['alpha']:g}-d{r['dens']:g}-seed{r['seed']}"
        saved = HEAVIEST.get((route, name))
        if not saved:
            raise SystemExit(f"heaviest.json has no lattices for {name}: run build/heaviest.py")
        saved_by_run.append(saved)
    steps = sorted({int(t) for saved in saved_by_run for t in saved if int(t) <= horizon})
    for r, saved in zip(runs, saved_by_run):    # every seed at every step, so no mean drops a seed
        missing = [t for t in steps if saved.get(str(t)) is None]
        if missing:
            raise SystemExit(f"{route} {describe(r)}: no saved lattice at steps {missing}")
    density = np.array([[island_density(saved[str(t)]["densest"]) for t in steps]
                        for saved in saved_by_run])
    late_steps = [t for t in steps if in_late_window(t, horizon)]
    if not late_steps:
        raise SystemExit(f"{route} {describe(runs[0])}: no saved lattice in the late window")
    per_seed = density[:, [steps.index(t) for t in late_steps]].mean(axis=1)
    return {"steps": steps, "seeds": len(runs), "rhoDensest": rounded(density.mean(axis=0)),
            "late": {"rhoDensest": rounded([per_seed.mean(), per_seed.min(), per_seed.max()])},
            "lateSteps": late_steps}


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
    entries = []
    for alpha in sorted({r["alpha"] for r in runs}):
        for L in sorted({r["size"] for r in runs}):
            g = group(runs, alpha=alpha, size=L)
            if not g:
                continue
            check_group(g, f"larger box L {L} alpha {alpha:g}")   # also: every seed ends at one step
            horizon = g[0]["reached"]
            at20000 = horizon >= 20000   # the same runs at step 20,000, like the denser box
            entries.append({
                "alpha": alpha, "L": L, "horizon": horizon,
                "curves": {key: curve(g, key, steps) for key in CURVE_KEYS},
                "late": late_summary(g),
                "late20000": late_summary(cut(g, 20000)) if at20000 else None,
                "checkpoints": densest_block("size", g, horizon),
                "checkpoints20000": densest_block("size", g, 20000) if at20000 else None,
            })
    return {"steps": steps, "density": BOX_DENS, "entries": entries}


def density_route(steps):
    runs = load(DENSITY_DIR, 20000)
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
                    "checkpoints": densest_block("density", g, 20000),
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
