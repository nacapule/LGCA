#!/usr/bin/env python3
"""native/sweep.py: run many sweep jobs in parallel with the native program.

Usage:

    python3 native/sweep.py CONFIG DIR [--processes N] [--end T] [--dry-run]

CONFIG is a TOML file listing the jobs:

    end = 20000                # the step every job runs to (a checkpoint)

    [[grid]]                   # one or more grids; each is every combination
    sizes = [90]               # of its four lists
    alphas = [0, 0.1, 0.5, 1]
    densities = [0.1, 0.4]
    seeds = [12345, 777]
    # end = 40000              # a grid may set its own end

DIR is the folder that holds one folder per job, named as the sweep names them
(both-L90-a0.85-d0.4-seed12345), the layout the density and alpha sweep used. Each job
is one run of the native program,

    native/build/lgca_native run --job ID --end T DIR/ID

which writes run.json and the state-t<T>.bin.gz lattices (native/README.md). The script
itself computes nothing; it decides what each job needs, starts up to N native programs
at a time, and reports progress.

What happens to each job, judged from DIR/ID/run.json:

  done       run.json has reached the job's end step (it belongs to this job and to the
             pinned engine, and its last snap is at or past the end). Not rerun, but
             still handed to the native program, which checks the folder as it would
             before continuing it (the job's settings, the last lattice's size and
             SHA-256, its particle count) and writes nothing. A few milliseconds per job.
  continue   run.json stops at an earlier checkpoint: the native program continues it
             from there. This is how a finished sweep is extended: raise `end` in the
             config, or pass --end, and run again.
  start      no run.json yet: the native program starts the job from its seed.
  check      run.json is unreadable or belongs to another job or engine: handed to the
             native program, which refuses it with the reason (and changes nothing).

A done or check job that the native program refuses is reported as failed.

The longest jobs start first, so that no long job is left running alone at the end;
the checks of done folders come last.

Two sweeps cannot work in one DIR at once (a lock on DIR/.sweep.lock, which
native/trace_sweep.py takes too); and the native program locks each job folder while
it works there, so a job folder reached by another way (a direct `lgca_native run`, a
link from another DIR) is never written by two programs at once: the second is
refused, and its job counts as failed.

Stopping is safe at any moment (Ctrl-C, or closing the terminal): the native program
writes every file atomically, so each job folder always resumes from its last
checkpoint. Run the same command again to continue.

Options:

  --processes N   how many native programs run at once (default: the number of CPUs)
  --end T         run every job to step T instead of the ends in CONFIG
  --dry-run       list each job, the step it runs to and what would happen to it,
                  then stop; changes nothing
  --program PATH  another build of the native program (default native/build/lgca_native)

Exit status: 0 when every job has reached its end, 1 when some job failed (or a
native program could not be started), 2 for a wrong command line or config or an
output folder it cannot use (in use by another sweep.py, or not a folder), 130 when
stopped.

Needs Python 3.11 or later (for tomllib); standard library only.
"""

import argparse
import fcntl
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from decimal import Decimal
from pathlib import Path

if sys.version_info < (3, 11):
    sys.exit("sweep.py needs Python 3.11 or later (it reads the config with tomllib)")
import tomllib

HERE = Path(__file__).resolve().parent                 # native/
DEFAULT_PROGRAM = HERE / "build" / "lgca_native"

# The steps at which a job saves its lattice; a job's end must be one of them. Up to 80000
# they are the sweep's JavaScript runner's; after it they go on every 20000 steps up to
# 2560000. The same list is TICKS in native/sweep_job.hpp.
CHECKPOINTS = [0, 100, 500, 1000, 2000, 5000, 10000] + list(range(20000, 2560001, 20000))
CHECKPOINT_TEXT = "0 100 500 1000 2000 5000 10000, then every 20000 up to 2560000"

# The engine every job folder must come from: the "engineHash" that the JavaScript runner
# and the native program write into run.json (ENGINE_HASH in native/sweep_job.hpp). A
# folder with another value was written by another engine; it is never counted as done.
ENGINE_HASH = "34edbb46d96ee9aaea500e8c5e0927a5858366a2ed52fc3db1284d37b7d5ff6d"

# The settings the native program accepts (native/run_cli.cpp checks them again).
SIZE_RANGE = (4, 10000)
ALPHA_RANGE = (0, 2)
DENS_RANGE = (0, 357913941)
SEED_RANGE = (0, 2**32 - 1)


class ConfigError(Exception):
    """A wrong config file or command line: reported in one line, exit status 2."""


# ---------------------------------------------------------------------------------------
# Job names
# ---------------------------------------------------------------------------------------

def js_number(x):
    """Write the number x as JavaScript writes it (String(x), JSON.stringify(x)).

    Job names hold their alpha and dens written this way: "a0.85", "a1" (not "a1.0"),
    "d0.025". Python's repr differs (1.0, and an exponent from 1e-5 on), so this follows
    the ECMAScript rule (Number::toString) itself. It is json_number in native/measure.cpp.
    """
    x = float(x)
    if x == 0:
        return "0"
    # repr gives the shortest digits that read back as x, as JavaScript requires.
    _, digit_tuple, exponent = Decimal(repr(abs(x))).as_tuple()
    digits = "".join(map(str, digit_tuple))
    n = len(digits) + exponent             # x = 0.digits * 10**n
    digits = digits.rstrip("0")            # repr writes 90.0 as the digits 900
    k = len(digits)
    if k <= n <= 21:                       # an integer: the digits, then n - k zeros
        text = digits + "0" * (n - k)
    elif 0 < n <= 21:                      # a decimal point inside the digits
        text = digits[:n] + "." + digits[n:]
    elif -6 < n <= 0:                      # 0.000ddd
        text = "0." + "0" * -n + digits
    else:                                  # exponent form: d.ddde+X or de-X
        text = digits[0] + ("." + digits[1:] if k > 1 else "")
        text += ("e+" if n - 1 >= 0 else "e-") + str(abs(n - 1))
    return "-" + text if x < 0 else text


def job_name(size, alpha, dens, seed):
    """The sweep's name for a job, e.g. both-L90-a0.85-d0.4-seed12345."""
    return f"both-L{size}-a{js_number(alpha)}-d{js_number(dens)}-seed{seed}"


# ---------------------------------------------------------------------------------------
# The config file
# ---------------------------------------------------------------------------------------

def read_config(path, end_override=None):
    """Read CONFIG and return the jobs as a dict {name: job}.

    A job is a dict with size, alpha, dens, seed and end. A job listed by two grids is
    run once, to the larger of their ends.
    """
    return config_jobs(load_config(path), path, end_override)


def load_config(path):
    """The TOML file at `path`, as a dict."""
    try:
        with open(path, "rb") as f:
            return tomllib.load(f)
    except OSError as e:
        raise ConfigError(f"cannot read {path}: {e.strerror}")
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{path} is not valid TOML: {e}")


def config_jobs(config, path, end_override=None, extra_keys=()):
    """The jobs of a config read by load_config, as read_config returns them. The
    config's top level may hold `end`, [[grid]] and the settings in `extra_keys` (which
    the caller checks)."""
    unknown = set(config) - {"end", "grid", *extra_keys}
    if unknown:
        allowed = ", ".join(["end", "[[grid]]", *extra_keys])
        raise ConfigError(f"{path}: unknown setting(s) {', '.join(sorted(unknown))} "
                          f"(allowed: {allowed})")
    grids = config.get("grid")
    if not isinstance(grids, list) or not grids:
        raise ConfigError(f"{path}: no [[grid]] table")

    jobs = {}
    for number, grid in enumerate(grids, start=1):
        where = f"{path}, grid {number}"
        if not isinstance(grid, dict):
            raise ConfigError(f"{where} is not a table (write each grid as [[grid]])")
        unknown = set(grid) - {"sizes", "alphas", "densities", "seeds", "end"}
        if unknown:
            raise ConfigError(f"{where}: unknown setting(s) {', '.join(sorted(unknown))} "
                              f"(allowed: sizes, alphas, densities, seeds, end)")
        sizes = number_list(grid, "sizes", where, SIZE_RANGE, integers=True)
        alphas = number_list(grid, "alphas", where, ALPHA_RANGE)
        densities = number_list(grid, "densities", where, DENS_RANGE)
        seeds = number_list(grid, "seeds", where, SEED_RANGE, integers=True)
        end = end_override if end_override is not None else grid.get("end", config.get("end"))
        if end is None:
            raise ConfigError(f"{where}: no end step (set `end` at the top or in the grid)")
        check_end(end, where)

        for size in sizes:
            for alpha in alphas:
                for dens in densities:
                    for seed in seeds:
                        if "e" in js_number(alpha) + js_number(dens):
                            # JavaScript writes values below 0.000001 with an exponent;
                            # the sweep's names allow only digits and dots.
                            raise ConfigError(f"{where}: alpha {alpha} or dens {dens} has "
                                              f"no job name (use 0 or at least 0.000001)")
                        name = job_name(size, alpha, dens, seed)
                        old = jobs.get(name)
                        jobs[name] = dict(size=size, alpha=alpha, dens=dens, seed=seed,
                                          end=max(end, old["end"]) if old else end)
    return jobs


def number_list(grid, key, where, allowed, integers=False):
    """grid[key], checked to be a non-empty list of numbers within `allowed`."""
    values = grid.get(key)
    if not isinstance(values, list) or not values:
        raise ConfigError(f"{where}: `{key}` must be a non-empty list")
    low, high = allowed
    for v in values:
        # (bool is a kind of int in Python; `true` is not a number here.)
        is_number = isinstance(v, int) if integers else isinstance(v, (int, float))
        if isinstance(v, bool) or not is_number:
            kind = "whole numbers" if integers else "numbers"
            raise ConfigError(f"{where}: `{key}` must hold {kind}, found {v!r}")
        if not low <= v <= high:
            raise ConfigError(f"{where}: {key} value {v} is outside {low} to {high}")
    return values


def check_end(end, where):
    if isinstance(end, bool) or not isinstance(end, int) or end not in CHECKPOINTS:
        raise ConfigError(f"{where}: end must be one of the checkpoints {CHECKPOINT_TEXT}, "
                          f"found {end!r}")


# ---------------------------------------------------------------------------------------
# What each job needs
# ---------------------------------------------------------------------------------------

def job_state(folder, name, end):
    """What the job in `folder` needs, as (state, last checkpoint reached, note).

    state is "done", "continue", "start" or "check" (see the top of this file). Only
    run.json is read, as the JavaScript runner read it to decide what to skip; the native
    program checks the rest (the settings, the last lattice's SHA-256) when it gets the job.
    """
    path = folder / "run.json"
    if not path.exists():
        return "start", 0, "new"
    try:
        record = json.loads(path.read_text())
        reached = record["snaps"][-1]["t"]
        same_job = record["job"]["id"] == name
        same_engine = record["engineHash"] == ENGINE_HASH
        if isinstance(reached, bool) or not isinstance(reached, int):
            raise TypeError("t is not a whole number")
    except (OSError, ValueError, KeyError, IndexError, TypeError) as e:
        return "check", 0, f"run.json unreadable ({type(e).__name__})"
    if not same_job:
        return "check", 0, "run.json is another job's"
    if not same_engine:
        return "check", 0, "run.json comes from another engine"
    if reached >= end:
        return "done", reached, f"t={reached}"
    return "continue", reached, f"t={reached}"


def cost(job, reached):
    """A rough measure of a job's running time, used to start the long jobs first and
    to estimate the time left: lattice sites x steps still to run x (0.5 + dens). The
    last factor was measured: at L = 90 a job at dens 2.4 took 5.5 to 6.5 times as
    long as one at dens 0.025. A larger alpha also slows a job (alpha 1 about
    twice as slow as alpha 0.8 or less), which this leaves out. A done job costs 0:
    it is only checked."""
    return job["size"] ** 2 * max(job["end"] - reached, 0) * (0.5 + job["dens"])


# ---------------------------------------------------------------------------------------
# Running the jobs
# ---------------------------------------------------------------------------------------

def format_time(seconds):
    seconds = int(round(seconds))
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m{seconds % 60:02d}s"
    return f"{seconds // 3600}h{seconds // 60 % 60:02d}m"


class Progress:
    """The progress report. On a terminal: one line, rewritten in place. In a log file:
    one line per finished job. A job counts as finished when its native program has
    ended: run to its end, checked (a done job), or failed."""

    def __init__(self, total, total_cost):
        self.total = total
        self.finished = self.failed = 0
        self.total_cost, self.finished_cost = total_cost, 0
        self.start = time.monotonic()
        self.live = sys.stdout.isatty()
        self.shown = ""

    def line(self, running):
        """e.g. "[57/132] 8 running, 0 failed, 1m23s, about 2m10s left". The time left
        assumes the jobs still to finish run at the speed of those finished so far."""
        elapsed = time.monotonic() - self.start
        text = (f"[{self.finished}/{self.total}] {running} running, "
                f"{self.failed} failed, {format_time(elapsed)}")
        if 0 < self.finished_cost < self.total_cost:
            left = elapsed * (self.total_cost - self.finished_cost) / self.finished_cost
            text += f", about {format_time(left)} left"
        return text

    def show(self, running):
        text = self.line(running)
        if self.live and text != self.shown:
            print("\r" + text.ljust(len(self.shown)), end="", flush=True)
            self.shown = text

    def job_finished(self, name, ok, job_cost, output, running):
        self.finished += 1
        self.finished_cost += job_cost
        self.failed += not ok
        output = output.strip() or "(no output)"
        if self.live:
            if not ok:
                print("\r" + " " * len(self.shown) + "\r" + f"FAILED {name}:\n{output}")
                self.shown = ""
            self.show(running)
        else:
            print(f"{self.line(running)}  " +
                  (output.splitlines()[-1] if ok else f"FAILED {name}:\n{output}"), flush=True)

    def close(self):
        if self.live and self.shown:
            print()


def run_jobs(command, todo, processes, progress, lock):
    """Run the jobs in `todo`, a list of (name, job, cost), at most `processes` at a
    time; `command(name, job)` is the job's command line (a list). Returns the names of
    the jobs that failed, or None if it was stopped.

    A stop signal (see request_stop) only sets `stop_signal`; the loop sees it within
    0.05 seconds, then stops the running native programs and waits for them. An error
    here (a native program that cannot be started, say) stops them the same way and
    is then passed on.
    """
    waiting = list(todo)
    running = {}            # name -> (process, file holding its output, cost)
    failed = []
    try:
        while (waiting or running) and stop_signal is None:
            # Start jobs while there is room.
            while waiting and len(running) < processes and stop_signal is None:
                name, job, job_cost = waiting.pop(0)
                output = tempfile.TemporaryFile(mode="w+")
                try:
                    process = subprocess.Popen(
                        command(name, job),
                        stdout=output, stderr=subprocess.STDOUT,
                        # The native programs inherit the lock on DIR (see lock_folder),
                        # so it stays held until the last of them has stopped, even if
                        # this script is killed outright.
                        pass_fds=(lock,),
                        # A --program wrapper may launch the native program as a
                        # child. Give the whole job a group of its own so stopping
                        # the sweep stops both the wrapper and its children.
                        start_new_session=True)
                except OSError:
                    output.close()
                    raise
                running[name] = (process, output, job_cost)

            # Collect the jobs that have ended.
            for name in [n for n, (p, _, _) in running.items() if p.poll() is not None]:
                if stop_signal is not None:
                    break   # a requested stop is not a failed job
                process, output, job_cost = running.pop(name)
                # A wrapper must wait for the program it launches. If it exits
                # early, stop what it left running and report a failed job.
                left_children = stop_programs([process])
                output.seek(0)
                text = output.read()
                output.close()
                ok = process.returncode == 0 and not left_children
                if left_children:
                    text += "\nerror: the job program exited while its children were still running\n"
                if not ok:
                    failed.append(name)
                progress.job_finished(name, ok, job_cost, text,
                                      len(running))
            progress.show(len(running))
            time.sleep(0.05)
        return failed if stop_signal is None else None
    finally:
        # Stopped, or an error: end the native programs still running (each job folder
        # keeps its last checkpoint) and wait until they are gone. After a normal end
        # there are none left.
        stop_programs(process for process, _, _ in running.values())
        for process, output, _ in running.values():
            output.close()


def signal_group(pid, signum, deadline):
    """Send `signum` to the process group `pid`: True if it was sent, False if the
    group is gone. macOS can refuse (EPERM) a group that is still being torn down,
    even right after its last process was reaped (about 1 time in 300 for a job that
    has just ended); the refusal is retried until the group is gone. One that lasts
    past `deadline` (time.monotonic()) is passed on: it is not evidence that the group
    stopped."""
    while True:
        try:
            os.killpg(pid, signum)
            return True
        except ProcessLookupError:
            return False
        except PermissionError:
            if time.monotonic() >= deadline:
                raise
            time.sleep(0.005)


def stop_programs(programs):
    """Stop each job's process group, including children of a --program wrapper.

    Returns whether any group received SIGTERM. The native program normally stops
    at once; its atomic files keep the last checkpoint. Give our direct children
    up to five seconds to exit, then send SIGKILL to any remaining original groups
    and reap our children. Descendants are signalled too; an already-exited wrapper
    can cause SIGKILL immediately. We cannot reap a wrapper's children ourselves.
    """
    programs = list(programs)
    groups = set()
    deadline = time.monotonic() + 5
    for process in programs:
        if signal_group(process.pid, signal.SIGTERM, deadline):
            groups.add(process.pid)
    were_present = bool(groups)
    deadline = time.monotonic() + 5
    for process in programs:
        try:
            process.wait(timeout=max(0, deadline - time.monotonic()))
        except subprocess.TimeoutExpired:
            pass
    # A reaped wrapper may still have children, so finish every group that received
    # SIGTERM, even when its original process has exited. Avoid probing with signal 0:
    # macOS can return EPERM while a group is exiting rather than report it gone.
    deadline = time.monotonic() + 5
    for pid in groups:
        signal_group(pid, signal.SIGKILL, deadline)
    for process in programs:
        process.wait()
    return were_present


def lock_folder(out):
    """Make sure only one sweep works in `out` at a time (sweep.py or trace_sweep.py,
    which takes the same lock): two programs writing the same job folder could corrupt
    it. Returns the locked file descriptor, or None if the folder is in use."""
    lock = os.open(out / ".sweep.lock", os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        os.close(lock)
        return None
    return lock


# The signal that asked the sweep to stop, or None.
stop_signal = None


def request_stop(signum, frame):
    """Ctrl-C, `kill` or closing the terminal: ask run_jobs to stop. (It only sets a
    flag, so a stop never falls between starting a job and recording it.)"""
    global stop_signal
    stop_signal = signum


# ---------------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run sweep jobs in parallel with the native program "
                    "(the top of native/sweep.py explains it).")
    parser.add_argument("config", type=Path, help="the TOML file listing the jobs")
    parser.add_argument("out", type=Path, metavar="DIR", help="the folder of job folders")
    parser.add_argument("--processes", type=int, default=os.cpu_count() or 1,
                        help="native programs running at once (default: the CPU count)")
    parser.add_argument("--end", type=int, help="run every job to this step instead")
    parser.add_argument("--dry-run", action="store_true",
                        help="list each job, its end and what would happen to it; change nothing")
    parser.add_argument("--program", type=Path, default=DEFAULT_PROGRAM,
                        help="the native program (default: native/build/lgca_native)")
    args = parser.parse_args()

    try:
        if args.processes < 1:
            raise ConfigError("--processes must be at least 1")
        if args.end is not None:
            check_end(args.end, "--end")
        jobs = read_config(args.config, args.end)
        # The program by its full path: a bare name such as "lgca_native" would
        # otherwise be looked up in PATH instead of the current folder. A path that
        # cannot be resolved is no program either (Python 3.11 and 3.12 raise
        # RuntimeError on a symlink loop, which later versions resolve quietly).
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
    states = {name: job_state(args.out / name, name, job["end"]) for name, job in jobs.items()}
    count = {s: sum(state == s for state, _, _ in states.values())
             for s in ("done", "continue", "start", "check")}
    summary = (f"{len(jobs)} jobs: {count['done']} done, {count['continue']} to continue, "
               f"{count['start']} to start" +
               (f", {count['check']} to check" if count["check"] else ""))

    if args.dry_run:
        for name in sorted(jobs):
            state, _, note = states[name]
            print(f"{state:9} {name}  to {jobs[name]['end']}  ({note})")
        print(summary)
        return 0

    # Every job goes to the native program, the longest first; the done ones (cost 0)
    # come last and are only checked.
    todo = sorted(((name, job, cost(job, states[name][1])) for name, job in jobs.items()),
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

    # From here on a stop is handled by run_jobs (it also stops the native programs).
    for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        signal.signal(signum, request_stop)
    print(f"{summary}; running {min(args.processes, len(todo))} at a time in {args.out}" +
          (" (done folders are checked, not rerun)" if count["done"] else ""), flush=True)
    progress = Progress(len(jobs), sum(c for _, _, c in todo))
    try:
        failed = run_jobs(lambda name, job: [str(program), "run", "--job", name, "--end",
                                             str(job["end"]), str(args.out / name)],
                          todo, args.processes, progress, lock)
    except OSError as e:
        progress.close()
        print(f"error: {e}; the native programs started were stopped, and every job "
              f"folder resumes from its last checkpoint", file=sys.stderr)
        return 1
    progress.close()
    if failed is None:
        print("stopped; every job folder resumes from its last checkpoint. "
              "Run the same command again to continue.")
        return 130

    elapsed = format_time(time.monotonic() - progress.start)
    if failed:
        print(f"{len(failed)} of {len(todo)} jobs FAILED ({elapsed}): {', '.join(failed)}")
        return 1
    run = len(todo) - count["done"]
    print(f"all {len(jobs)} jobs have reached their end ({run} run, " +
          (f"{count['done']} checked, " if count["done"] else "") + f"{elapsed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
