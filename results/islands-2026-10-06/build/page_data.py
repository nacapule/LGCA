"""The page's data: page-data.json.

    python3 build/page_data.py

Seed means of the island density over time, at 100 to 120 log-spaced sample steps, and
late values: the mean over the last quarter of each run, with the smallest and largest seed
value. The larger-box runs go on to step 80,000 (40,000 at L = 720); their late values are
also given at step 20,000, where the denser box ends. The densest island comes from
heaviest.json (build/heaviest.py) and exists only at the saved lattices. Values are rounded
to 4 significant digits; a curve is null after its runs end.
"""
import json
import warnings
from datetime import datetime

import numpy as np

from islands import BOX_DENS, DENSITY_DIR, DENSITY_L, SIZE_DIR, STUDY, group, late, load

CURVE_KEYS = ["rho", "rhoOld"]      # island density, old cluster density
HEAVIEST = {}                       # (route, job name) -> {step: islands}, filled in main()


def sample_steps(end):
    """The runs sample every 25 steps to 5000, then every 100. This keeps 100 to 120 of
    those steps, log-spaced, and step 0."""
    steps = {0}
    for t in np.geomspace(25, end, 160):
        spacing = 25 if t <= 5000 else 100
        steps.add(int(round(t / spacing)) * spacing)
    return sorted(s for s in steps if s <= end)


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
    """Seed mean at `steps` (null where a run has no sample, or no seed has an old cluster)."""
    index = {int(s): k for k, s in enumerate(runs[0]["t"])}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean([r[key] for r in runs], axis=0)
    return rounded([mean[index[s]] if s in index else None for s in steps])


def late_summary(runs):
    out = {"seeds": len(runs), "N": int(round(np.mean([r["N0"] for r in runs]))),
           "mean": float(f"{np.mean([r['mean'] for r in runs]):.4g}"), "T": runs[0]["reached"]}
    for key in CURVE_KEYS:
        values = [late(r, key) for r in runs]
        out[key] = rounded([np.mean(values), np.min(values), np.max(values)])
    return out


def densest_block(route, runs, horizon):
    """The densest island at the saved lattices up to `horizon`: seed means of its density,
    and `late` = [seed mean, min, max] over the lattices in the last quarter (steps 60,000
    and 80,000 for an 80,000-step run, the last lattice otherwise)."""
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
            raise SystemExit(f"{route} alpha {r['alpha']} size {r['size']} dens {r['dens']} "
                             f"seed {r['seed']}: no saved lattice at steps {missing}")
    density = np.array([[saved[str(t)]["densest"][0] / saved[str(t)]["densest"][1] for t in steps]
                        for saved in saved_by_run])
    late_steps = [t for t in steps if t >= 0.75 * horizon]
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


def size_route():
    runs = load(SIZE_DIR)
    steps = sample_steps(max(r["reached"] for r in runs))
    entries = []
    for alpha in sorted({r["alpha"] for r in runs}):
        for L in sorted({r["size"] for r in runs}):
            g = group(runs, alpha=alpha, size=L)
            if not g:
                continue
            horizon = min(r["reached"] for r in g)
            g = cut(g, horizon)          # every seed at the common horizon of this alpha and L
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


def density_route():
    runs = load(DENSITY_DIR, 20000)
    steps = sample_steps(20000)
    entries = []
    for alpha in sorted({r["alpha"] for r in runs}):
        for dens in sorted({r["dens"] for r in runs}):
            g = group(runs, alpha=alpha, dens=dens)
            if g:
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
    data = {
        "generated": datetime.now().isoformat(timespec="minutes"),
        "settings": {"model": "boson", "alignment": "polar", "kernel": "power (divide by M^alpha)",
                     "field": "centre + 6 neighbours", "sensitivity": 6,
                     "boundaries": "periodic", "lattice": "hexagonal, 6 velocity channels"},
        "size": size_route(),
        "density": density_route(),
    }
    path = STUDY / "page-data.json"
    path.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False))
    print("wrote", path, f"{path.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
