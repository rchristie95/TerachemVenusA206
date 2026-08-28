#!/usr/bin/env python3
"""Divergence guard for a running numerical job.

This project's production runs have collapsed silently before, so the kill
decision lives here -- in the driver, polled every couple of minutes -- and not
in the judgement of a model that only wakes up once an hour.

Watches the metrics CSV a run appends to (the format run_overnight.py writes:
time_ps, ..., potential_kj_mol, temperature_k, cos_alpha ...) and fails the run on:

  nan       a NaN or inf appeared in any numeric column
  stalled   the file has not been touched for --stall-min minutes
  diverged  potential energy drifted >--drift-frac from the run's own baseline,
            or temperature left the thermostat by >--temp-tol K

The baseline is the run's own first --baseline-rows samples rather than a
hard-coded energy, so the same guard works for any system size.

Usage:
  python health.py --csv path/to/run_metrics.csv
  python health.py --csv ... --grace-min 15      # before the first row exists
Exit code 0 = healthy, 1 = kill it.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

NUMERIC_HINTS = ("potential", "temperature", "energy", "alpha", "ang", "cos")


def _rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _floats(rows: list[dict], col: str) -> list[float]:
    out = []
    for r in rows:
        try:
            out.append(float(r[col]))
        except (TypeError, ValueError, KeyError):
            out.append(math.nan)
    return out


def _median(xs: list[float]) -> float:
    xs = sorted(x for x in xs if math.isfinite(x))
    if not xs:
        return math.nan
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def check(
    csv_path: Path,
    baseline_rows: int = 20,
    tail_rows: int = 10,
    drift_frac: float = 0.05,
    temp_tol: float = 25.0,
    stall_min: float = 25.0,
    grace_min: float = 20.0,
    started_at: float | None = None,
) -> tuple[bool, str, dict]:
    info: dict = {"csv": str(csv_path)}

    if not csv_path.exists():
        age_min = (time.time() - started_at) / 60 if started_at else 0.0
        if age_min > grace_min:
            return False, f"stalled: no metrics file after {age_min:.0f} min", info
        return True, "warming up (no metrics file yet)", info

    age_min = (time.time() - csv_path.stat().st_mtime) / 60.0
    info["mtime_age_min"] = round(age_min, 1)
    if age_min > stall_min:
        return False, f"stalled: metrics CSV untouched for {age_min:.0f} min", info

    rows = _rows(csv_path)
    info["n_rows"] = len(rows)
    if len(rows) < 2:
        return True, "warming up (<2 rows)", info

    cols = [c for c in rows[0] if any(h in c.lower() for h in NUMERIC_HINTS)]

    # 1. NaN / inf anywhere in the tail is fatal and unambiguous.
    for col in cols:
        tail = _floats(rows[-tail_rows:], col)
        if any(not math.isfinite(x) for x in tail):
            return False, f"nan: non-finite value in '{col}'", info

    baseline = rows[: min(baseline_rows, len(rows) // 2 or 1)]
    tail = rows[-tail_rows:]

    # 2. Potential energy must stay near where the run equilibrated.
    for col in cols:
        if "potential" not in col.lower() and "energy" not in col.lower():
            continue
        b, t = _median(_floats(baseline, col)), _median(_floats(tail, col))
        if not (math.isfinite(b) and math.isfinite(t)) or b == 0:
            continue
        drift = abs(t - b) / abs(b)
        info[f"{col}_baseline"] = round(b, 1)
        info[f"{col}_now"] = round(t, 1)
        info[f"{col}_drift_frac"] = round(drift, 4)
        if drift > drift_frac:
            return False, f"diverged: {col} drifted {drift:.1%} from baseline ({b:.4g} -> {t:.4g})", info

    # 3. Thermostat must still be holding.
    for col in cols:
        if "temperature" not in col.lower():
            continue
        b, t = _median(_floats(baseline, col)), _median(_floats(tail, col))
        if not (math.isfinite(b) and math.isfinite(t)):
            continue
        info[f"{col}_now"] = round(t, 2)
        if abs(t - b) > temp_tol:
            return False, f"diverged: {col} left the thermostat ({b:.1f} -> {t:.1f} K)", info

    return True, "healthy", info


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", type=Path, required=True)
    p.add_argument("--baseline-rows", type=int, default=20)
    p.add_argument("--tail-rows", type=int, default=10)
    p.add_argument("--drift-frac", type=float, default=0.05)
    p.add_argument("--temp-tol", type=float, default=25.0)
    p.add_argument("--stall-min", type=float, default=25.0)
    p.add_argument("--grace-min", type=float, default=20.0)
    p.add_argument("--started-at", type=float, default=None, help="unix time the job started")
    args = p.parse_args()

    ok, reason, info = check(
        args.csv, args.baseline_rows, args.tail_rows, args.drift_frac,
        args.temp_tol, args.stall_min, args.grace_min, args.started_at,
    )
    print(json.dumps({"ok": ok, "reason": reason, **info}))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
