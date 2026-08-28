#!/usr/bin/env python3
"""Durable hypothesis queue for the autonomous loop.

The ledger is the only state that survives a rate limit, a dead session or a
reboot.  Everything the loop knows lives in autolab/ledger.json; nothing
important is allowed to exist only in a model's context or a /tmp scratchpad.

Job lifecycle
-------------
  queued  -> running -> done          numerical run finished, result harvested
                     -> failed        crashed, diverged, or blew the wall clock
  queued  -> blocked                  a job it depends on failed; premise is gone
  blocked -> queued | dead            triaged by the analysis pass

`blocked` is the mechanism that stops compute being spent on dead ends: the
moment a job fails, everything downstream of it stops being claimable, while
independent work keeps the GPU busy.

Usage:
  python ledger.py add --claim "..." --branch detuning --cmd "..." --predicts "..."
  python ledger.py claim                 # driver: atomically take the next job
  python ledger.py report <id> --status done --result '{"chi2": 4.1}'
  python ledger.py list [--status queued]
  python ledger.py dead <branch> --why "..."   # retire a whole line of enquiry
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEDGER = HERE / "ledger.json"
LOCK = HERE / ".ledger.lock"
DEAD_ENDS = HERE / "DEAD_ENDS.md"
NEEDS_ANALYSIS = HERE / "NEEDS_ANALYSIS"

MAX_HOURS_CAP = 7.0  # hard ceiling on any single numerical run

TERMINAL = {"done", "failed", "dead"}


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@contextmanager
def locked():
    """Exclusive lock around every read-modify-write of the ledger."""
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    fh = open(LOCK, "w")
    try:
        fcntl.flock(fh, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(fh, fcntl.LOCK_UN)
        fh.close()


def _read() -> list[dict]:
    if not LEDGER.exists():
        return []
    return json.loads(LEDGER.read_text())


def _write(jobs: list[dict]) -> None:
    tmp = LEDGER.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(jobs, indent=2) + "\n")
    os.replace(tmp, LEDGER)  # atomic; a crash mid-write cannot corrupt the ledger


def dead_branches() -> set[str]:
    """Branches retired in DEAD_ENDS.md, as `## branch: <name>` headings."""
    if not DEAD_ENDS.exists():
        return set()
    out = set()
    for line in DEAD_ENDS.read_text().splitlines():
        if line.startswith("## branch:"):
            out.add(line.split(":", 1)[1].strip())
    return out


def next_id(jobs: list[dict]) -> str:
    n = 0
    for j in jobs:
        if j["id"].startswith("h") and j["id"][1:].isdigit():
            n = max(n, int(j["id"][1:]))
    return f"h{n + 1:03d}"


def cmd_add(args) -> None:
    with locked():
        jobs = _read()
        dead = dead_branches()
        if args.branch in dead and not args.force:
            sys.exit(
                f"refusing: branch '{args.branch}' is retired in DEAD_ENDS.md.\n"
                "Re-read why it died before spending compute on it; --force to override."
            )
        if any(j["claim"].strip() == args.claim.strip() for j in jobs):
            sys.exit("refusing: a job with this exact claim already exists.")

        job = {
            "id": next_id(jobs),
            "status": "queued",
            "branch": args.branch,
            "claim": args.claim,
            "predicts": args.predicts,
            "falsifies": args.falsifies,
            "cmd": args.cmd,
            "precheck": args.precheck,
            "max_hours": min(args.max_hours, MAX_HOURS_CAP),
            "monitor_csv": args.monitor_csv,
            "depends_on": args.depends_on or [],
            "created": now(),
            "started": None,
            "finished": None,
            "exit_code": None,
            "failure_kind": None,
            "failure_note": None,
            "result": None,
            "verdict": None,
        }
        jobs.append(job)
        _write(jobs)
    print(job["id"])


def cmd_claim(args) -> None:
    """Atomically hand the driver the next runnable job, or exit 1 if there is none."""
    with locked():
        jobs = _read()
        by_id = {j["id"]: j for j in jobs}
        dead = dead_branches()
        dirty = False
        for j in jobs:
            if j["status"] != "queued":
                continue
            if j["branch"] in dead:
                # Retired while it sat in the queue: never spend the GPU on it.
                j["status"] = "dead"
                dirty = True
                continue
            deps = [by_id.get(d) for d in j["depends_on"]]
            if any(d is None or d["status"] != "done" for d in deps):
                continue
            j["status"] = "running"
            j["started"] = now()
            _write(jobs)
            print(json.dumps(j))
            return
        if dirty:
            _write(jobs)
    sys.exit(1)


def _block_dependents(jobs: list[dict], failed_id: str) -> list[str]:
    """Transitively block everything whose premise just died."""
    blocked, frontier = [], {failed_id}
    while frontier:
        nxt = set()
        for j in jobs:
            if j["status"] in ("queued", "blocked") and set(j["depends_on"]) & frontier:
                if j["status"] == "queued":
                    j["status"] = "blocked"
                    blocked.append(j["id"])
                nxt.add(j["id"])
        frontier = nxt - frontier
    return blocked


def cmd_report(args) -> None:
    with locked():
        jobs = _read()
        job = next((j for j in jobs if j["id"] == args.id), None)
        if job is None:
            sys.exit(f"no such job: {args.id}")
        job["status"] = args.status
        job["finished"] = now()
        job["exit_code"] = args.exit_code
        job["failure_kind"] = args.failure_kind
        job["failure_note"] = args.failure_note
        if args.result:
            job["result"] = json.loads(args.result)

        blocked = []
        if args.status == "failed":
            blocked = _block_dependents(jobs, job["id"])
            # A failure is a decision point: the loop must think before it spends
            # more GPU-hours down the same line.
            NEEDS_ANALYSIS.write_text(
                f"{now()} {job['id']} failed ({args.failure_kind}): "
                f"{args.failure_note}\nblocked: {', '.join(blocked) or 'none'}\n"
            )
        _write(jobs)
    print(json.dumps({"id": args.id, "status": args.status, "blocked": blocked}))


def cmd_verdict(args) -> None:
    """Record what a finished job actually decided about its claim."""
    with locked():
        jobs = _read()
        job = next((j for j in jobs if j["id"] == args.id), None)
        if job is None:
            sys.exit(f"no such job: {args.id}")
        job["verdict"] = {"call": args.call, "note": args.note, "at": now()}
        _write(jobs)
    print(f"{args.id}: {args.call}")


def cmd_unblock(args) -> None:
    with locked():
        jobs = _read()
        job = next((j for j in jobs if j["id"] == args.id), None)
        if job is None:
            sys.exit(f"no such job: {args.id}")
        job["status"] = args.to
        _write(jobs)
    print(f"{args.id} -> {args.to}")


def cmd_dead(args) -> None:
    """Retire a branch so no future pass can re-enqueue it by accident."""
    with locked():
        jobs = _read()
        entry = (
            f"\n## branch: {args.branch}\n\n"
            f"- retired: {now()}\n"
            f"- why: {args.why}\n"
            f"- evidence: {args.evidence or 'see ledger'}\n"
        )
        with open(DEAD_ENDS, "a") as fh:
            fh.write(entry)
        killed = []
        for j in jobs:
            if j["branch"] == args.branch and j["status"] in ("queued", "blocked"):
                j["status"] = "dead"
                killed.append(j["id"])
        _write(jobs)
    print(f"retired branch '{args.branch}'; killed {killed or 'no pending jobs'}")


def cmd_list(args) -> None:
    jobs = _read()
    if args.status:
        jobs = [j for j in jobs if j["status"] == args.status]
    if args.json:
        print(json.dumps(jobs, indent=2))
        return
    if not jobs:
        print("(no jobs)")
        return
    width = max(len(j["branch"]) for j in jobs)
    for j in jobs:
        v = (j.get("verdict") or {}).get("call", "")
        flag = f"  [{v}]" if v else ""
        print(f"{j['id']}  {j['status']:<8} {j['branch']:<{width}}  {j['claim'][:72]}{flag}")


def cmd_status(args) -> None:
    """One-screen summary: what the loop should do next."""
    jobs = _read()
    counts: dict[str, int] = {}
    for j in jobs:
        counts[j["status"]] = counts.get(j["status"], 0) + 1
    print("counts:", json.dumps(counts))
    print("needs_analysis:", NEEDS_ANALYSIS.exists())
    unharvested = [j["id"] for j in jobs if j["status"] == "done" and not j.get("verdict")]
    print("done but unharvested:", unharvested or "none")
    print("blocked:", [j["id"] for j in jobs if j["status"] == "blocked"] or "none")
    print("dead branches:", sorted(dead_branches()) or "none")
    runnable = [j for j in jobs if j["status"] == "queued"]
    print("queued:", [j["id"] for j in runnable] or "none")
    if not runnable and not unharvested and not NEEDS_ANALYSIS.exists():
        print("\nIDLE: queue is empty and nothing needs triage.")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("add", help="enqueue a hypothesis")
    a.add_argument("--claim", required=True, help="the falsifiable statement")
    a.add_argument("--branch", required=True, help="line of enquiry, for dead-end retirement")
    a.add_argument("--cmd", required=True, help="shell command that tests it")
    a.add_argument("--predicts", required=True, help="the number that must come out, with a threshold")
    a.add_argument("--falsifies", default="", help="what result kills the claim")
    a.add_argument("--precheck", default="", help="cheap (<10 min) command that must pass first")
    a.add_argument("--max-hours", type=float, default=MAX_HOURS_CAP)
    a.add_argument("--monitor-csv", default="", help="metrics CSV the health guard watches")
    a.add_argument("--depends-on", nargs="*", default=[])
    a.add_argument("--force", action="store_true")
    a.set_defaults(func=cmd_add)

    c = sub.add_parser("claim", help="take the next runnable job (driver only)")
    c.set_defaults(func=cmd_claim)

    r = sub.add_parser("report", help="record a job outcome")
    r.add_argument("id")
    r.add_argument("--status", required=True, choices=["done", "failed"])
    r.add_argument("--exit-code", type=int, default=None)
    r.add_argument("--failure-kind", default=None,
                   choices=[None, "precheck", "crash", "diverged", "stalled", "timeout"])
    r.add_argument("--failure-note", default=None)
    r.add_argument("--result", default=None, help="JSON blob of harvested numbers")
    r.set_defaults(func=cmd_report)

    v = sub.add_parser("verdict", help="record what a finished job decided")
    v.add_argument("id")
    v.add_argument("--call", required=True, choices=["supported", "refuted", "inconclusive"])
    v.add_argument("--note", required=True)
    v.set_defaults(func=cmd_verdict)

    u = sub.add_parser("unblock", help="re-queue or kill a blocked job after triage")
    u.add_argument("id")
    u.add_argument("--to", required=True, choices=["queued", "dead"])
    u.set_defaults(func=cmd_unblock)

    d = sub.add_parser("dead", help="retire a branch")
    d.add_argument("branch")
    d.add_argument("--why", required=True)
    d.add_argument("--evidence", default="")
    d.set_defaults(func=cmd_dead)

    ls = sub.add_parser("list")
    ls.add_argument("--status", default=None)
    ls.add_argument("--json", action="store_true")
    ls.set_defaults(func=cmd_list)

    st = sub.add_parser("status")
    st.set_defaults(func=cmd_status)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
