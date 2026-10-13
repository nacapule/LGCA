#!/usr/bin/env python3
"""native/trace_sweep.py: trace many sweep jobs in parallel (`lgca_native trace`).

Usage:

    python3 native/trace_sweep.py CONFIG OUT [--processes N] [--end T] [--dry-run]
                                  [--program PATH] [--no-check]

CONFIG is a TOML file in the format of native/sweeps/*.toml (an `end` and [[grid]]
tables of sizes, alphas, densities and seeds; native/sweep.py describes it and checks it
the same way here), with two optional top-level keys:

    states = true          also save each job's lattice at every checkpoint
                           (OUT/<job>/state-t<T>.bin.gz), which lets a later, longer
                           trace continue from the last one instead of the seed
    check = ["DIR", ...]   folders holding earlier runs of the same jobs
                           (DIR/<job>/run.json, written by `lgca_native run` or the
                           sweep's JavaScript runner). A job found there is checked
                           against it: every saved row and checkpoint must be reproduced
                           exactly.
                           The first folder that has the job is used. A relative path
                           is resolved against the repository root (the folder above
                           native/), so a copy of native/ elsewhere needs absolute
                           paths. A check folder that does not exist is an error;
                           --no-check traces without them.

Each job writes OUT/<job>/trace.json (trace_cli.cpp describes its contents). What
happens to each job, judged from that file:

  done       the trace reached the job's end, with the current columns, checked
             against the job's check folder when there is one (it names that folder)
             and, with states, with a state file for each of its checkpoints. Still
             handed to the native program, which checks it in full (every row and
             checkpoint, the check fingerprint of the check folder's rows and the
             trace itself, every state file's lattice) and writes nothing. Done jobs
             come last.
  continue   the trace stops before the end, or lacks the check or state files: the
             native program continues or completes it (with states it resumes from the
             last saved checkpoint, otherwise it replays from the seed).
  start      no trace yet.
  check      trace.json is unreadable, another job's or engine's, or has other columns
             (written before the densest-island columns were added): handed to the
             native program, which refuses it with the reason and changes nothing; such
             a job fails. Trace those jobs into a new OUT.

The longest jobs start first. Two sweeps cannot work in one OUT at once (a lock on
OUT/.sweep.lock, the lock native/sweep.py takes), and the native program locks each
job's folder while it works there. Stopping is safe at any moment (Ctrl-C, SIGTERM, or
closing the terminal): the running programs are stopped (after five seconds, killed),
and every trace keeps its last complete state. Run the same command again to continue.

Options:

  --processes N   how many native programs run at once (default: the number of CPUs)
  --end T         trace every job to step T instead of the ends in CONFIG
  --dry-run       list each job, the step it is traced to and what would happen to it,
                  then stop; changes nothing
  --program PATH  another build of the native program (default native/build/lgca_native)
  --no-check      leave out the config's check folders (for a checkout without them)

Exit status: 0 when every job has reached its end, 1 when some job failed (or a native
program could not be started), 2 for a wrong command line or config, a missing check
folder or program, or an output folder it cannot use, 130 when stopped.

Needs Python 3.11 or later (for tomllib); standard library only.
"""

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path

# The run sweep's config checks, job costs and process handling.
from sweep import (ENGINE_HASH, ConfigError, Progress, check_end, config_jobs, cost, format_time,
                   load_config, lock_folder, request_stop, run_jobs)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROGRAM = ROOT / "native" / "build" / "lgca_native"
# The columns lgca_native trace writes (COLUMNS in trace_cli.cpp; keep the two in step).
COLUMNS = ["t", "N", "nmax", "kmax", "occupied", "threshold", "clusterMass", "clusterArea",
           "seed", "islandMass", "islandArea", "islandPolar", "wraps", "ci", "cj", "rg",
           "components", "largestMass", "largestArea", "densestMass", "densestArea",
           "densest2Mass", "densest2Area"]


def check_folder(check_roots, name):
    """The first earlier run of the job to check against, or None."""
    for root in check_roots:
        if (root / name / "run.json").exists():
            return root / name
    return None


def job_state(folder, name, end, check, states):
    """What the job in `folder` needs, as (state, step its trace can resume from, note):
    state is "done", "continue", "start" or "check" (see the top of this file). Only
    trace.json and the names of the state files are looked at; the native program checks
    the rest when it gets the job."""
    path = folder / "trace.json"
    if not path.exists():
        return "start", 0, "new"
    try:
        trace = json.loads(path.read_text())
        same = trace["job"]["id"] == name and trace["engineHash"] == ENGINE_HASH
        columns = trace["columns"] == COLUMNS
        reached = trace["reached"]
        snaps = [snap["t"] for snap in trace["snaps"]]
        checked = trace["checked"]
        if isinstance(reached, bool) or not isinstance(reached, int):
            raise TypeError("reached is not a whole number")
    except (OSError, ValueError, KeyError, TypeError) as e:
        return "check", 0, f"trace.json unreadable ({type(e).__name__})"
    if not same:
        return "check", 0, "trace.json is another job's or engine's"
    if not columns:
        return "check", 0, "trace.json has other columns"
    resume = reached if states else 0
    if reached < end:
        return "continue", resume, f"t={reached}"
    if check is not None and not (isinstance(checked, dict) and checked.get("dir") == str(check)):
        return "continue", 0, f"t={reached}, not checked against {check}"
    if states and not all((folder / f"state-t{t}.bin.gz").exists() for t in snaps):
        return "continue", 0, f"t={reached}, state files missing"
    return "done", reached, f"t={reached}"


def main():
    parser = argparse.ArgumentParser(
        description="Trace sweep jobs in parallel with the native program "
                    "(the top of native/trace_sweep.py explains it).")
    parser.add_argument("config", type=Path, help="the TOML file listing the jobs")
    parser.add_argument("out", type=Path, metavar="OUT", help="the folder of job folders")
    parser.add_argument("--processes", type=int, default=os.cpu_count() or 1,
                        help="native programs running at once (default: the CPU count)")
    parser.add_argument("--end", type=int, help="trace every job to this step instead")
    parser.add_argument("--dry-run", action="store_true",
                        help="list each job, its end and what would happen to it; change nothing")
    parser.add_argument("--program", type=Path, default=DEFAULT_PROGRAM,
                        help="the native program (default: native/build/lgca_native)")
    parser.add_argument("--no-check", action="store_true",
                        help="leave out the config's check folders")
    args = parser.parse_args()

    try:
        if args.processes < 1:
            raise ConfigError("--processes must be at least 1")
        if args.end is not None:
            check_end(args.end, "--end")
        config = load_config(args.config)
        jobs = config_jobs(config, args.config, args.end, extra_keys=("states", "check"))
        states = config.get("states", False)
        if not isinstance(states, bool):
            raise ConfigError(f"{args.config}: `states` must be true or false, found {states!r}")
        check = config.get("check", [])
        if not isinstance(check, list) or not all(isinstance(d, str) and d for d in check):
            raise ConfigError(f"{args.config}: `check` must be a list of folder names, "
                              f"found {check!r}")
        check_roots = [] if args.no_check else [ROOT / d for d in check]
        for root in check_roots:
            if not root.is_dir():
                raise ConfigError(f"check folder {root} does not exist (a relative path is "
                                  f"resolved against {ROOT}); give --no-check to trace "
                                  f"without checking")
        # The program by its full path (sweep.py explains why).
        try:
            program = args.program.resolve()
        except (OSError, RuntimeError):
            program = None
        if not args.dry_run and not (program and program.is_file() and os.access(program, os.X_OK)):
            raise ConfigError(f"no native program at {args.program} "
                              f"(build it with: bash native/build.sh)")
    except ConfigError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    # What each job needs.
    found = {name: job_state(args.out / name, name, job["end"],
                             check_folder(check_roots, name), states)
             for name, job in jobs.items()}
    count = {s: sum(state == s for state, _, _ in found.values())
             for s in ("done", "continue", "start", "check")}
    summary = (f"{len(jobs)} jobs: {count['done']} done, {count['continue']} to continue, "
               f"{count['start']} to start" +
               (f", {count['check']} to check" if count["check"] else ""))

    if args.dry_run:
        for name in sorted(jobs):
            state, _, note = found[name]
            print(f"{state:9} {name}  to {jobs[name]['end']}  ({note})")
        print(summary)
        return 0

    # Every job goes to the native program, the longest first; the done ones (cost 0)
    # come last and are only checked.
    todo = sorted(((name, job, 0 if found[name][0] == "done" else cost(job, found[name][1]))
                   for name, job in jobs.items()),
                  key=lambda item: (-item[2], item[0]))

    try:
        args.out.mkdir(parents=True, exist_ok=True)
        lock = lock_folder(args.out)
    except OSError as e:
        reason = "not a folder" if args.out.exists() and not args.out.is_dir() else e.strerror
        print(f"error: cannot use {args.out}: {reason}", file=sys.stderr)
        return 2
    if lock is None:
        print(f"error: another sweep (sweep.py or trace_sweep.py) is working in {args.out} "
              f"(or the jobs it started are still finishing)", file=sys.stderr)
        return 2

    def command(name, job):
        folder = args.out / name
        argv = [str(program), "trace", "--job", name, "--end", str(job["end"])]
        if states:
            argv += ["--states", str(folder)]
        check = check_folder(check_roots, name)
        if check is not None:
            argv += ["--check", str(check)]
        return argv + [str(folder / "trace.json")]

    # From here on a stop is handled by run_jobs (it also stops the native programs).
    for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        signal.signal(signum, request_stop)
    print(f"{summary}; running {min(args.processes, len(todo))} at a time in {args.out}" +
          (" (done traces are checked, not traced again)" if count["done"] else "") +
          (" (no check folders: --no-check)" if args.no_check and check else ""), flush=True)
    progress = Progress(len(jobs), sum(c for _, _, c in todo))
    try:
        failed = run_jobs(command, todo, args.processes, progress, lock)
    except OSError as e:
        progress.close()
        print(f"error: {e}; the native programs started were stopped, and every trace keeps "
              f"its last complete state", file=sys.stderr)
        return 1
    progress.close()
    if failed is None:
        print("stopped; every trace keeps its last complete state. "
              "Run the same command again to continue.")
        return 130

    elapsed = format_time(time.monotonic() - progress.start)
    if failed:
        print(f"{len(failed)} of {len(todo)} jobs FAILED ({elapsed}): {', '.join(failed)}")
        return 1
    traced = len(todo) - count["done"]
    print(f"all {len(jobs)} jobs have reached their end ({traced} traced, " +
          (f"{count['done']} checked, " if count["done"] else "") + f"{elapsed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
