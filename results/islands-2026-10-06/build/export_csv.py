"""The study's data as CSV.

    python3 build/export_csv.py

Reads page-data.json (build/page_data.py) and heaviest.json (build/heaviest.py) and writes:

  island-density-vs-time.csv   seed-mean island density at each sampled step
  late-island-density.csv      late island density, with the smallest and largest seed value
  islands-at-saved-steps.csv   the islands of every saved lattice, one row per lattice
"""
import csv
import json

from islands import STUDY

BOX = {"density": "denser", "size": "larger"}


def write(name, fields, rows):
    with (STUDY / name).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    print("wrote", name, len(rows), "rows")


def cell(value):
    if value is None:
        return ""
    return f"{value:g}" if isinstance(value, float) else value


def settings(route, e, data):
    side = e["L"] if route == "size" else data["density"]["L"]
    dens = data["size"]["density"] if route == "size" else e["dens"]
    N = e["late"]["N"] if route == "size" else e["N"]
    return [BOX[route], e["alpha"], side, dens, N]


def window(late):
    """The late window: from the first sample after three quarters of the run to its end."""
    return f"{int(0.75 * late['T']) + 100}-{late['T']}"


def main():
    data = json.loads((STUDY / "page-data.json").read_text())
    entries = [(route, e) for route in ("density", "size") for e in data[route]["entries"]]

    rows = []
    for route, e in entries:
        head = settings(route, e, data) + [e["late"]["seeds"]]
        for k, step in enumerate(data[route]["steps"]):
            old, new = e["curves"]["rhoOld"][k], e["curves"]["rho"][k]
            if old is not None or new is not None:
                rows.append(head + [step, cell(old), cell(new)])
    write("island-density-vs-time.csv",
          ["box", "alpha", "side", "dens", "N", "seeds", "step", "old", "new"], rows)

    rows = []
    for route, e in entries:
        head = settings(route, e, data)
        # the larger-box runs go past step 20,000: their late values at 20,000 and at the end
        lates = [x for x in (e.get("late20000"), e["late"]) if x]
        blocks = [x for x in (e.get("checkpoints20000"), e["checkpoints"]) if x]
        for island, key in (("old", "rhoOld"), ("new", "rho")):
            for late in lates:
                rows.append(head + [island, window(late), late["seeds"], *map(cell, late[key])])
        for block in blocks:
            steps = "+".join(str(t) for t in block["lateSteps"])
            rows.append(head + ["densest", steps, block["seeds"], *map(cell, block["late"]["rhoDensest"])])
    write("late-island-density.csv",
          ["box", "alpha", "side", "dens", "N", "island", "steps", "seeds", "mean",
           "seed_min", "seed_max"], rows)

    rows = []
    for run in json.loads((STUDY / "heaviest.json").read_text())["runs"]:
        for step, m in sorted(run["steps"].items(), key=lambda item: int(item[0])):
            densest2 = m["densest2"] or [None, None]
            rows.append([BOX[run["route"]], run["alpha"], run["size"], run["dens"], run["seed"],
                         int(step), m["N"], m["components"], *m["island"], *m["heaviest"],
                         *m["densest"], cell(densest2[0]), cell(densest2[1])])
    rows.sort(key=lambda r: r[:6])
    write("islands-at-saved-steps.csv",
          ["box", "alpha", "side", "dens", "seed", "step", "N", "islands", "new_mass",
           "new_sites", "heaviest_mass", "heaviest_sites", "densest_mass", "densest_sites",
           "densest2_mass", "densest2_sites"], rows)


if __name__ == "__main__":
    main()
