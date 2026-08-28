#!/usr/bin/env python3
"""The compute loop. Deliberately stupid, and never talks to Claude.

Claim a job -> run its cheap precheck -> run the numerical job under a hard
wall-clock cap with the divergence guard polling alongside -> report -> repeat.

Nothing here makes a scientific decision. It only decides whether a run is
still worth electricity. Because it holds no API credentials and needs no
session, a rate limit or a closed laptop does not stop it: work keeps draining
from the queue, and the analysis pass finds finished results waiting.

  python driver.py                  # run until the queue is dry
  python driver.py --idle-forever   # keep waiting for new work instead
  touch autolab/STOP                # halt after the current job

Wall clock: every job is capped at --max-hours (7 h ceiling, enforced here and
again in ledger.py at enqueue time).
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
LEDGER = HERE / "ledger.py"
HEALTH = HERE / "health.py"
STOP = HERE / "STOP"
LOG = HERE / "driver.log"

HARD_CAP_HOURS = 7.0
PRECHECK_TIMEOUT_S = 600  # a precheck that takes >10 min is not a precheck
POLL_S = 120              # how often the divergence guard runs
TICK_S = 5                # how fast the loop notices exit / STOP / wall clock


def log(msg: str) -> None:
    line = f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def ledger(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(LEDGER), *args],
                          capture_output=True, text=True, cwd=REPO)


def claim() -> dict | None:
    r = ledger("claim")
    if r.returncode != 0:
        return None
    return json.loads(r.stdout)


def report(job_id: str, status: str, *, exit_code=None, kind=None, note=None, result=None) -> None:
    args = ["report", job_id, "--status", status]
    if exit_code is not None:
        args += ["--exit-code", str(exit_code)]
    if kind:
        args += ["--failure-kind", kind]
    if note:
        args += ["--failure-note", note[:2000]]
    if result:
        args += ["--result", json.dumps(result)]
    r = ledger(*args)
    log(f"  report -> {r.stdout.strip() or r.stderr.strip()}")


def run_precheck(job: dict) -> tuple[bool, str]:
    cmd = job.get("precheck") or ""
    if not cmd.strip():
        return True, ""
    log(f"  precheck: {cmd}")
    try:
        r = subprocess.run(cmd, shell=True, cwd=REPO, capture_output=True,
                           text=True, timeout=PRECHECK_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return False, f"precheck exceeded {PRECHECK_TIMEOUT_S}s"
    if r.returncode != 0:
        tail = (r.stderr or r.stdout or "").strip().splitlines()[-15:]
        return False, "precheck failed:\n" + "\n".join(tail)
    return True, ""


def health_of(csv_path: Path, started_at: float) -> tuple[bool, str]:
    r = subprocess.run(
        [sys.executable, str(HEALTH), "--csv", str(csv_path), "--started-at", str(started_at)],
        capture_output=True, text=True, cwd=REPO,
    )
    try:
        d = json.loads(r.stdout)
        return bool(d["ok"]), d["reason"]
    except (json.JSONDecodeError, KeyError):
        return True, "health check unreadable; not killing on that alone"


def kill_tree(proc: subprocess.Popen) -> None:
    """TERM the whole process group, then KILL what is left."""
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except ProcessLookupError:
        return
    for _ in range(60):
        if proc.poll() is not None:
            return
        time.sleep(1)
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except ProcessLookupError:
        pass


def kind_from_reason(reason: str) -> str:
    head = reason.split(":", 1)[0].strip()
    return {"nan": "diverged", "diverged": "diverged", "stalled": "stalled"}.get(head, "diverged")


def run_job(job: dict, max_hours_override: float | None) -> None:
    jid = job["id"]
    cap_h = min(float(job.get("max_hours") or HARD_CAP_HOURS),
                max_hours_override or HARD_CAP_HOURS, HARD_CAP_HOURS)
    deadline = time.time() + cap_h * 3600
    started = time.time()

    log(f"{jid} [{job['branch']}] {job['claim']}")
    log(f"  cap {cap_h:.2f} h | predicts: {job['predicts']}")

    ok, note = run_precheck(job)
    if not ok:
        log(f"  PRECHECK FAILED: {note.splitlines()[0] if note else ''}")
        report(jid, "failed", kind="precheck", note=note)
        return

    job_log = HERE / "runs" / f"{jid}.log"
    job_log.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(job["monitor_csv"]) if job.get("monitor_csv") else None
    if csv_path and not csv_path.is_absolute():
        csv_path = REPO / csv_path

    log(f"  exec: {job['cmd']}")
    with open(job_log, "w") as out:
        proc = subprocess.Popen(job["cmd"], shell=True, cwd=REPO, stdout=out,
                                stderr=subprocess.STDOUT, start_new_session=True)

        rc = None
        last_health = 0.0
        while True:
            # Return the instant the job exits; a cheap probe must not be held
            # hostage to the health-poll interval.
            try:
                rc = proc.wait(timeout=TICK_S)
                break
            except subprocess.TimeoutExpired:
                pass

            if time.time() > deadline:
                log(f"  WALL CLOCK: {cap_h:.2f} h cap hit; terminating")
                kill_tree(proc)
                report(jid, "failed", kind="timeout",
                       note=f"hit the {cap_h:.2f} h cap; partial output in {job_log}")
                return

            if STOP.exists():
                log("  STOP sentinel; terminating and re-queueing")
                kill_tree(proc)
                ledger("unblock", jid, "--to", "queued")
                return

            if csv_path is not None and time.time() - last_health > POLL_S:
                last_health = time.time()
                healthy, reason = health_of(csv_path, started)
                if not healthy:
                    log(f"  UNHEALTHY: {reason}; terminating")
                    kill_tree(proc)
                    report(jid, "failed", kind=kind_from_reason(reason), note=reason)
                    return

    elapsed_h = (time.time() - started) / 3600
    if rc != 0:
        tail = "\n".join(job_log.read_text(errors="replace").splitlines()[-25:])
        log(f"  CRASH rc={rc} after {elapsed_h:.2f} h")
        report(jid, "failed", exit_code=rc, kind="crash", note=f"exit {rc}\n{tail}")
        return

    # A clean exit still has to pass the guard: a run can finish "successfully"
    # having quietly melted an hour earlier.
    if csv_path is not None:
        healthy, reason = health_of(csv_path, started)
        if not healthy:
            log(f"  exited 0 but UNHEALTHY: {reason}")
            report(jid, "failed", exit_code=0, kind=kind_from_reason(reason),
                   note=f"exited cleanly but {reason}")
            return

    log(f"  DONE in {elapsed_h:.2f} h -> {job_log}")
    report(jid, "done", exit_code=0,
           result={"elapsed_hours": round(elapsed_h, 3), "log": str(job_log),
                   "monitor_csv": str(csv_path) if csv_path else None})


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--max-hours", type=float, default=None,
                   help=f"override per-job cap downward (ceiling {HARD_CAP_HOURS} h)")
    p.add_argument("--idle-forever", action="store_true",
                   help="wait for new work instead of exiting when the queue is dry")
    p.add_argument("--idle-sleep", type=int, default=300)
    args = p.parse_args()

    log("driver up")
    idle_notified = False
    while True:
        if STOP.exists():
            log("STOP sentinel present; driver exiting")
            return

        job = claim()
        if job is None:
            if not args.idle_forever:
                log("queue dry; driver exiting")
                return
            if not idle_notified:
                log(f"queue dry; idling (poll {args.idle_sleep}s)")
                idle_notified = True
            time.sleep(args.idle_sleep)
            continue

        idle_notified = False
        try:
            run_job(job, args.max_hours)
        except Exception as exc:  # never let the loop die on one bad job
            log(f"  driver exception: {exc!r}")
            report(job["id"], "failed", kind="crash", note=f"driver exception: {exc!r}")


if __name__ == "__main__":
    main()
