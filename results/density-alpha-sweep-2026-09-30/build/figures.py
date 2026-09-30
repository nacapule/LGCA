"""Figures of the density x alpha sweep.

python3 build/figures.py              makes the figures from the CSVs
python3 build/figures.py --data data  first makes the CSVs from the runs (not in repo)
"""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, NullFormatter

STUDY = Path(__file__).resolve().parent.parent
FIGURES = STUDY / "figures"

T = 20000
LATE_START = 15100  # last quarter of the run
MAIN_SIZE = 90
AREA = MAIN_SIZE * MAIN_SIZE
COARSE_ALPHAS = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
FOCUS_ALPHAS = [0.8, 0.85, 0.9, 1]
BOX_SIZES = [90, 120, 180, 240]
BOX_DENSITY = 0.4
BOX_SEEDS = [12345, 777, 424242]

SUMMARY_CSV = STUDY / "summary-t20000-prelim.csv"
TIME_CSV = STUDY / "cluster-density-vs-time-t20000.csv"
BOX_CSV = STUDY / "condensate-vs-box-size-t20000.csv"


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def curve_column(alpha, dens):
    return f"a{alpha:g}_d{dens:g}"


def load_finished_runs(data):
    runs = []
    for path in sorted(data.glob("both-L*/run.json")):
        record = json.loads(path.read_text())
        if record["snaps"][-1]["t"] < T:
            continue
        series = [sample for sample in record["series"] if sample["t"] <= T]
        job = record["job"]
        runs.append({
            "size": job["size"], "alpha": job["alpha"], "dens": job["dens"],
            "seed": job["seed"],
            "t": np.array([sample["t"] for sample in series]),
            "N": series[0]["N"],
            "clusterDensity": np.array([s["clusterDensity"] for s in series]),
            "clusterFraction": np.array([s["clusterFraction"] for s in series]),
            "kmax": np.array([s["kmax"] for s in series]),
        })
    return runs


def late_mean(run, key):
    return run[key][run["t"] >= LATE_START].mean()


def export_main_sweep(runs):
    main = [run for run in runs if run["size"] == MAIN_SIZE]
    densities = sorted({run["dens"] for run in main})
    seeds_by_cell = {}
    for run in main:
        seeds_by_cell.setdefault((run["alpha"], run["dens"]), set()).add(run["seed"])
    # each alpha uses the seeds that finished for all densities (5, or 3 for finer alphas)
    row_seeds = {}
    for alpha in sorted({run["alpha"] for run in main}):
        common = set.intersection(*[seeds_by_cell.get((alpha, dens), set())
                                    for dens in densities])
        if len(common) >= 3:
            row_seeds[alpha] = common

    summary_rows = []
    curves = {}
    times = None
    for alpha, seeds in row_seeds.items():
        for dens in densities:
            group = [run for run in main if run["alpha"] == alpha
                     and run["dens"] == dens and run["seed"] in seeds]
            mean_density = np.mean([run["N"] for run in group]) / AREA
            late = [late_mean(run, "clusterDensity") for run in group]
            late_fraction = [late_mean(run, "clusterFraction") for run in group]
            summary_rows.append({
                "alpha": alpha,
                "dens_input": dens,
                "mean_density": round(mean_density, 4),
                "seeds": len(group),
                "T": T,
                "late_cluster_density": float(np.mean(late)),
                "seed_min": float(min(late)),
                "seed_max": float(max(late)),
                "ratio_to_mean": float(np.mean(late) / mean_density),
                "cutoff": math.ceil(2 * mean_density),
                "late_cluster_mass_fraction": float(np.mean(late_fraction)),
            })
            if times is None:
                times = group[0]["t"]
            for run in group:
                if not np.array_equal(run["t"], times):
                    raise ValueError("sample times differ between runs")
            curves[curve_column(alpha, dens)] = np.mean(
                [run["clusterDensity"] for run in group], axis=0)

    write_csv(SUMMARY_CSV, list(summary_rows[0]), summary_rows)
    time_rows = []
    for i, t in enumerate(times):
        row = {"t": int(t)}
        row.update({name: f"{values[i]:.6g}" for name, values in curves.items()})
        time_rows.append(row)
    write_csv(TIME_CSV, ["t", *curves], time_rows)


def export_box_sizes(runs):
    rows = []
    for alpha in FOCUS_ALPHAS:
        for size in BOX_SIZES:
            group = [run for run in runs if run["size"] == size
                     and run["alpha"] == alpha and run["dens"] == BOX_DENSITY
                     and run["seed"] in BOX_SEEDS]
            if len(group) != len(BOX_SEEDS):
                raise ValueError(f"missing box-size runs for alpha {alpha}, L {size}")
            late = [late_mean(run, "kmax") for run in group]
            particles = np.mean([run["N"] for run in group])
            rows.append({
                "alpha": alpha,
                "L": size,
                "N": round(particles),
                "mean_density": round(particles / size**2, 4),
                "seeds": len(group),
                "T": T,
                "late_kmax": round(float(np.mean(late)), 2),
                "seed_min": round(float(min(late)), 2),
                "seed_max": round(float(max(late)), 2),
            })
    write_csv(BOX_CSV, list(rows[0]), rows)


def number_label(value, _position=None):
    return f"{value:g}" if value < 1000 else f"{value / 1000:g}k"


def seed_label(alphas, seeds):
    counts = sorted({seeds[alpha] for alpha in alphas})
    if len(counts) == 1:
        return f"mean of {counts[0]} seeds"
    return f"mean of {counts[0]}–{counts[-1]} seeds (per panel)"


def density_colours(densities):
    norm = LogNorm(6 * densities[0], 6 * densities[-1])
    cmap = plt.get_cmap("viridis")
    colour = {dens: cmap(norm(6 * dens)) for dens in densities}
    scale = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    return colour, scale


def save(fig, name):
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / name, dpi=150)
    plt.close(fig)


def plot_time_grid(alphas, times, curves, densities, seeds, name):
    colour, scale = density_colours(densities)
    mixed = len({seeds[alpha] for alpha in alphas}) > 1
    ncol = 4
    nrow = math.ceil(len(alphas) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.55 * nrow + 0.6),
                             layout="constrained", squeeze=False)
    for ax, alpha in zip(axes.flat, alphas):
        for dens in densities:
            values = curves[curve_column(alpha, dens)]
            ax.plot(times[1:], values[1:], color=colour[dens], lw=1.1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(25, T)
        ax.grid(alpha=0.18, which="major")
        ax.xaxis.set_major_formatter(FuncFormatter(number_label))
        ax.yaxis.set_major_formatter(FuncFormatter(number_label))
        ax.yaxis.set_minor_formatter(NullFormatter())
        title = f"α = {alpha:g}" + (f" · {seeds[alpha]} seeds" if mixed else "")
        ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    for ax in axes.flat[len(alphas):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("Cluster density\n(particles/site)")
    for ax in axes[-1, :]:
        ax.set_xlabel("Simulation step")
    fig.colorbar(scale, ax=axes, shrink=0.6, pad=0.01).set_label(
        "Starting average density (particles/site)")
    fig.suptitle(f"Cluster density vs time · {len(densities)} starting densities · "
                 f"90×90, centre + neighbours, sens 6 · {seed_label(alphas, seeds)} · "
                 f"to {T:,} steps", x=0.01, ha="left", fontsize=10, color="#5a676c")
    save(fig, name)


def plot_focus(times, curves, densities, seeds):
    colour, scale = density_colours(densities)
    fig, axes = plt.subplots(2, len(FOCUS_ALPHAS), figsize=(16, 7),
                             layout="constrained", squeeze=False, sharey=True,
                             gridspec_kw={"height_ratios": [3, 1]})
    ticks = [0, 100, 200, 300, 500, 1000, 2000, 5000, 10000, 20000]
    for column, alpha in enumerate(FOCUS_ALPHAS):
        top, bottom = axes[:, column]
        for dens in densities:
            values = curves[curve_column(alpha, dens)]
            top.plot(times, values, color=colour[dens], lw=1.1)
            early = times <= 500
            top.plot(times[early], values[early], ls="none", marker="o", ms=2.6,
                     color=colour[dens])
            bottom.plot(times, values, color=colour[dens], lw=1.1)
        top.set_xscale("symlog", linthresh=500, linscale=1.5)
        top.set_xlim(0, T)
        top.set_xticks(ticks)
        top.xaxis.set_minor_formatter(NullFormatter())
        top.axvline(500, color="#9aa3a7", lw=0.8, ls=":")
        top.text(500, 0.015, " linear ← | → log ", transform=top.get_xaxis_transform(),
                 ha="center", va="bottom", fontsize=7.5, color="#7a8488",
                 backgroundcolor="#fcfcfa")
        bottom.set_xlim(0, T)
        bottom.set_xlabel("Simulation step")
        for ax in (top, bottom):
            ax.set_yscale("log")
            ax.grid(alpha=0.18)
            ax.xaxis.set_major_formatter(FuncFormatter(number_label))
            ax.yaxis.set_major_formatter(FuncFormatter(number_label))
            ax.yaxis.set_minor_formatter(NullFormatter())
        top.set_title(f"α = {alpha:g} · {seeds[alpha]} seeds",
                      loc="left", fontsize=11, fontweight="bold")
    axes[0, 0].set_ylabel("Cluster density (particles/site)\n"
                          "time axis linear to 500, log after")
    axes[1, 0].set_ylabel("full run, linear time")
    fig.colorbar(scale, ax=axes, shrink=0.7, pad=0.01).set_label(
        "Starting average density (particles/site)")
    fig.suptitle("Cluster density vs time · 90×90, centre + neighbours, sens 6 · "
                 f"{seed_label(FOCUS_ALPHAS, seeds)}",
                 x=0.01, ha="left", fontsize=10, color="#5a676c")
    save(fig, "focus-transition-t20000-prelim.png")


def plot_box_sizes(rows):
    cmap = plt.get_cmap("plasma")
    fig, ax = plt.subplots(figsize=(7.6, 5), layout="constrained")
    for k, alpha in enumerate(FOCUS_ALPHAS):
        group = sorted((row for row in rows if row["alpha"] == alpha),
                       key=lambda row: row["L"])
        sizes = [row["L"] for row in group]
        means = [row["late_kmax"] for row in group]
        below = [row["late_kmax"] - row["seed_min"] for row in group]
        above = [row["seed_max"] - row["late_kmax"] for row in group]
        ax.errorbar(sizes, means, yerr=[below, above], fmt="-o", ms=4, capsize=3,
                    lw=1.4, color=cmap(0.85 * k / (len(FOCUS_ALPHAS) - 1)),
                    label=f"α = {alpha:g}")
    first = sorted((row for row in rows if row["alpha"] == FOCUS_ALPHAS[0]),
                   key=lambda row: row["L"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(BOX_SIZES, [f"{row['L']:g}\nN ≈ {number_label(round(row['N'], -2))}"
                              for row in first])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(number_label))
    ax.grid(alpha=0.18, which="major")
    ax.set_xlabel("Box side L (particle number N grows with the area)")
    ax.set_ylabel("Condensate size $N_c$\n(most particles in one channel of one site)")
    ax.legend(frameon=False)
    ax.set_title(f"Condensate size vs box size · average density ≈ "
                 f"{first[0]['mean_density']:.1f} particles/site\n"
                 f"late window, steps {LATE_START:,}–{T:,} · mean of {len(BOX_SEEDS)} seeds · "
                 "bars = range over seeds", loc="left", fontsize=9.5, color="#5a676c")
    save(fig, "condensate-vs-box-size-t20000.png")


def draw_figures():
    summary = read_csv(SUMMARY_CSV)
    seeds = {float(row["alpha"]): int(row["seeds"]) for row in summary}
    densities = sorted({float(row["dens_input"]) for row in summary})
    alphas = sorted(seeds)

    table = read_csv(TIME_CSV)
    times = np.array([int(row["t"]) for row in table])
    curves = {name: np.array([float(row[name]) for row in table])
              for name in table[0] if name != "t"}
    # 0 = no cluster yet (t = 0, densest starts), don't plot it on log axis
    curves = {name: np.where(values > 0, values, np.nan)
              for name, values in curves.items()}

    box_rows = [{key: float(value) for key, value in row.items()}
                for row in read_csv(BOX_CSV)]

    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.labelcolor": "#19252d", "text.color": "#19252d",
        "xtick.color": "#19252d", "ytick.color": "#19252d",
        "figure.facecolor": "#fcfcfa", "axes.facecolor": "#fcfcfa",
        "savefig.facecolor": "#fcfcfa",
    })
    coarse = [alpha for alpha in alphas if alpha in COARSE_ALPHAS]
    fine = [alpha for alpha in alphas
            if alpha not in COARSE_ALPHAS or alpha in (0.8, 0.9, 1)]
    plot_time_grid(coarse, times, curves, densities, seeds, "grid-time-t20000-prelim.png")
    plot_time_grid(fine, times, curves, densities, seeds,
                   "grid-time-refine-t20000-prelim.png")
    plot_focus(times, curves, densities, seeds)
    plot_box_sizes(box_rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, help="folder with the saved runs")
    args = parser.parse_args()
    if args.data:
        finished = load_finished_runs(args.data)
        export_main_sweep(finished)
        export_box_sizes(finished)
    draw_figures()
