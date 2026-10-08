"""Read the island traces of this study.

Each run's trace.json has one row per sample step: N, the island's mass and number of
sites, the old cluster's mass and area, the densest island's mass and number of sites, and
other columns.

Densities: an island or cluster with no sites has density 0. The old cluster has no sites
when no site holds max(1, ceil(2N/A)) particles: at step 0 in the 90 x 90 boxes with 3.6 or
more particles per site, where no site starts with that many. The native measurement also
reports its density as 0 then.
"""
import json
import math
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent
SIZE_DIR = STUDY / "size"          # more particles from a larger box (dens 0.4)
DENSITY_DIR = STUDY / "density"    # more particles from a denser 90 x 90 box

BOX_DENS = 0.4                     # per channel: 2.4 particles per site
DENSITY_L = 90

# Late values: the samples at steps t with LATE_FRACTION * T < t <= T, where T is the end
# of the run. The same window for every island.
LATE_FRACTION = 0.75
# Every averaged curve and late value needs at least this many distinct seeds.
MIN_SEEDS = 4


def load(folder, horizon=None):
    """Every trace in `folder` that reached `horizon` (all of them if None), as dicts
    of job settings plus one array per column, cut at the horizon."""
    runs = []
    for path in sorted(folder.glob("*/trace.json")):
        record = json.loads(path.read_text())
        if horizon is not None and record["reached"] < horizon:
            continue
        columns = record["columns"]
        rows = np.array([[np.nan if v is None else v for v in row] for row in record["rows"]],
                        dtype=float)
        if horizon is not None:
            rows = rows[rows[:, 0] <= horizon]
        run = {key: record["job"][key] for key in ("size", "alpha", "dens", "seed")}
        run["reached"] = record["reached"] if horizon is None else horizon
        run["checked"] = record["checked"]
        for k, name in enumerate(columns):
            run[name] = rows[:, k]
        run["seedSite"] = run["seed"]          # the trace column: the busiest site's index
        run["seed"] = record["job"]["seed"]    # the job's random seed
        add_derived(run)
        runs.append(run)
    return runs


def density(mass, area, steps=None, what="island"):
    """Particles per site; 0 where there are no sites and no particles (no island, or no old
    cluster). An undefined (null) mass or area, or particles without sites or sites without
    particles, stops."""
    mass, area = np.asarray(mass, float), np.asarray(area, float)
    bad = ~(np.isfinite(mass) & np.isfinite(area)) | ((mass == 0) != (area == 0))
    if bad.any():
        k = np.flatnonzero(bad)[0]
        at = f"step {steps[k]:g}" if steps is not None else f"sample {k}"
        raise SystemExit(f"{what}: {mass.flat[k]:g} particles on {area.flat[k]:g} sites at {at}")
    return np.divide(mass, area, out=np.zeros_like(mass), where=area > 0)


def add_derived(run):
    sites = run["size"] ** 2
    N = run["N"][0]
    run["N0"] = N
    run["sites"] = sites
    run["mean"] = N / sites
    where = describe(run)
    with np.errstate(invalid="ignore", divide="ignore"):
        run["rho"] = density(run["islandMass"], run["islandArea"], run["t"], f"{where}: new island")
        run["frac"] = run["islandMass"] / N                          # island's share of N
        run["ratio"] = run["rho"] / run["mean"]
        run["rhoOld"] = density(run["clusterMass"], run["clusterArea"], run["t"], f"{where}: old cluster")
        run["rhoDensest"] = density(run["densestMass"], run["densestArea"], run["t"],
                                    f"{where}: densest island")
        run["fracOld"] = run["clusterMass"] / N                      # old cluster's share
        run["occFrac"] = run["occupied"] / sites
        run["compsPerSite"] = run["components"] / sites
        run["isHeaviest"] = (run["islandMass"] == run["largestMass"]).astype(float)


def in_late_window(t, horizon):
    """True for the steps of the late window: LATE_FRACTION * horizon < t <= horizon."""
    t = np.asarray(t)
    return (t > LATE_FRACTION * horizon) & (t <= horizon)


def late(run, key, strict=False):
    """Mean of `key` over the late window of the run. Columns that can be undefined at a
    sample (such as rg) skip those samples; with strict=True an undefined sample stops."""
    values = run[key][in_late_window(run["t"], run["reached"])]
    if not len(values):
        raise SystemExit(f"{describe(run)}: no samples in the late window")
    if strict:
        if not np.isfinite(values).all():
            raise SystemExit(f"{describe(run)}: undefined {key} in the late window")
        return float(np.mean(values))
    return float(np.nanmean(values))


def describe(run):
    return f"L {run['size']} alpha {run['alpha']:g} dens {run['dens']:g} seed {run['seed']}"


def sample_steps(end):
    """The steps the runs sample up to step `end`: every 25 steps to step 5,000, then every 100."""
    return np.array([t for t in range(0, int(end) + 1, 25) if t <= 5000 or t % 100 == 0], float)


def check_group(runs, what):
    """Stop unless `runs` are at least MIN_SEEDS runs with distinct seeds, the same
    horizon and every sample step up to it, so that their means are seed means at equal times."""
    seeds = [r["seed"] for r in runs]
    if len(set(seeds)) != len(seeds):
        raise SystemExit(f"{what}: a seed appears twice ({sorted(seeds)})")
    if len(seeds) < MIN_SEEDS:
        raise SystemExit(f"{what}: {len(seeds)} seeds ({sorted(seeds)}), at least {MIN_SEEDS} needed")
    if len({r["reached"] for r in runs}) != 1:
        raise SystemExit(f"{what}: the runs end at different steps "
                         f"({sorted({r['reached'] for r in runs})})")
    steps = sample_steps(runs[0]["reached"])
    for r in runs:
        if not np.array_equal(r["t"], steps):
            diff = sorted(set(r["t"].tolist()) ^ set(steps.tolist()))
            raise SystemExit(f"{what}: seed {r['seed']} does not have exactly the sample steps up to "
                             f"step {r['reached']}" + (f" (step {diff[0]:g})" if diff else ""))
    return len(seeds)


def group(runs, **settings):
    return [r for r in runs if all(math.isclose(r[k], v) for k, v in settings.items())]
