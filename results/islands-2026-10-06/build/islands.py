"""Read the island traces of this study.

Each run's trace.json has one row per sample step: N, the island's mass and number of
sites, the old cluster's mass and area, and other columns.
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


def add_derived(run):
    sites = run["size"] ** 2
    N = run["N"][0]
    run["N0"] = N
    run["sites"] = sites
    run["mean"] = N / sites
    with np.errstate(invalid="ignore", divide="ignore"):
        run["rho"] = run["islandMass"] / run["islandArea"]           # island density
        run["frac"] = run["islandMass"] / N                          # island's share of N
        run["ratio"] = run["rho"] / run["mean"]
        run["rhoOld"] = np.where(run["clusterArea"] > 0,
                                 run["clusterMass"] / run["clusterArea"], np.nan)
        run["fracOld"] = run["clusterMass"] / N                      # old cluster's share
        run["occFrac"] = run["occupied"] / sites
        run["compsPerSite"] = run["components"] / sites
        run["isHeaviest"] = (run["islandMass"] == run["largestMass"]).astype(float)


def late(run, key, start_fraction=0.75):
    """Mean of `key` over the last quarter of the run (t > 3/4 of its horizon)."""
    t = run["t"]
    return float(np.nanmean(run[key][t > start_fraction * run["reached"]]))


def group(runs, **settings):
    return [r for r in runs if all(math.isclose(r[k], v) for k, v in settings.items())]
