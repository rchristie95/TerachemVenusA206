#!/usr/bin/env python3
"""Bounded, resumable one-process-per-GPU TeraChem launcher."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def complete(output: Path) -> bool:
    return output.is_file() and "Job finished:" in output.read_text(errors="replace")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("frame_dir", type=Path)
    p.add_argument("--terachem", type=Path, required=True)
    p.add_argument("--gpu", default="0")
    p.add_argument("--input", default="tddft.in")
    args = p.parse_args()
    if "," in args.gpu:
        raise ValueError("This smoke launcher accepts one explicit GPU; scale-out uses one invocation per GPU")
    records = []
    for site in ("A", "B"):
        work = args.frame_dir / f"site_{site}"
        inp, out = work / args.input, work / "tddft.out"
        if complete(out):
            records.append({"site": site, "status": "reused_complete", "output_sha256": sha256(out)})
            continue
        if out.exists():
            out.rename(work / f"tddft.incomplete.{int(time.time())}.out")
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = args.gpu
        started = time.time()
        with out.open("w") as log:
            result = subprocess.run([str(args.terachem), inp.name], cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
        record = {
            "site": site, "status": "complete" if result.returncode == 0 and complete(out) else "failed",
            "returncode": result.returncode, "wall_seconds": time.time() - started, "gpu": args.gpu,
            "input_sha256": sha256(inp), "geometry_sha256": sha256(work / "geometry.xyz"),
            "point_charges_sha256": sha256(work / "mm_charges.dat"), "output_sha256": sha256(out),
        }
        records.append(record)
        (work / "job.json").write_text(json.dumps(record, indent=2) + "\n")
        if record["status"] != "complete":
            break
    manifest = {"jobs": records, "completed_sites": sum(r["status"] in {"complete", "reused_complete"} for r in records)}
    (args.frame_dir / "launch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    if manifest["completed_sites"] != 2:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
