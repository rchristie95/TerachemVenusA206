#!/usr/bin/env python3
"""Run the 50-frame rapid TDDFT/TDC campaign with measured GPU concurrency."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
import threading
import time

from run_ensemble import (
    FrameFailure,
    ROOT,
    load_manifest,
    process_frame,
    save_manifest,
)


WINDOW = "long100ns50"
MANIFEST = ROOT / "manifests/long100ns50_frames.csv"
OUTDIR = ROOT / "rapid_100ns50"


def gpu_memory() -> tuple[int, int]:
    output = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=memory.used,memory.total",
            "--format=csv,noheader,nounits",
        ],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.splitlines()[0]
    used, total = (int(value.strip()) for value in output.split(","))
    return used, total


def assert_caps(frame: int) -> None:
    prep = json.loads(
        (ROOT / f"frames/{WINDOW}/frame_{frame:04d}/preparation.json").read_text()
    )
    caps = [
        cap
        for site in ("A", "B")
        for cap in prep["sites"][site]["link_caps"]
    ]
    if len(caps) != 6:
        raise RuntimeError(f"Frame {frame}: expected six link-cap diagnostics")
    for cap in caps:
        if abs(cap["cap_direction_dot_parent_to_partner"] - 1.0) > 1.0e-10:
            raise RuntimeError(f"Frame {frame}: reversed or invalid link cap: {cap}")


def worker(
    frame: int,
    row: dict[str, str],
    gpu: str,
    roots: int,
    profile: str,
    delta: float,
    reference: Path,
) -> dict:
    started = time.monotonic()
    try:
        result = process_frame(
            WINDOW,
            frame,
            row,
            gpu,
            {},
            roots,
            profile,
            delta,
            reference,
        )
        return {
            "frame": frame,
            "status": "complete",
            "row": row,
            "result": result,
            "wall_seconds": time.monotonic() - started,
        }
    except FrameFailure as error:
        row["convergence_status"] = "failed"
        row["failure_class"] = error.failure_class
        row["notes"] = str(error)
        return {
            "frame": frame,
            "status": "failed",
            "row": row,
            "failure_class": error.failure_class,
            "error": str(error),
            "wall_seconds": time.monotonic() - started,
        }
    except Exception as error:
        row["convergence_status"] = "failed"
        row["failure_class"] = "gpu_memory_or_numerical_failure"
        row["notes"] = repr(error)
        return {
            "frame": frame,
            "status": "failed",
            "row": row,
            "failure_class": row["failure_class"],
            "error": repr(error),
            "wall_seconds": time.monotonic() - started,
        }


def run_batch(
    frames: list[int],
    workers: int,
    rows: list[dict[str, str]],
    columns: list[str],
    gpu: str,
    roots: int,
    profile: str,
    delta: float,
    reference: Path,
) -> dict:
    if not frames:
        return {
            "frames": [],
            "workers": workers,
            "wall_seconds": 0.0,
            "peak_gpu_memory_mib": gpu_memory()[0],
            "gpu_memory_total_mib": gpu_memory()[1],
            "complete": 0,
            "failed": 0,
        }
    for frame in frames:
        rows[frame]["convergence_status"] = "running"
        rows[frame]["failure_class"] = ""
        rows[frame]["gpu"] = gpu
    save_manifest(MANIFEST, columns, rows)

    stop_monitor = threading.Event()
    samples: list[tuple[float, int, int]] = []

    def monitor() -> None:
        while not stop_monitor.is_set():
            try:
                used, total = gpu_memory()
                samples.append((time.time(), used, total))
            except Exception:
                pass
            stop_monitor.wait(1.0)

    monitor_thread = threading.Thread(target=monitor, daemon=True)
    monitor_thread.start()
    started = time.monotonic()
    outcomes: list[dict] = []
    try:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    worker,
                    frame,
                    dict(rows[frame]),
                    gpu,
                    roots,
                    profile,
                    delta,
                    reference,
                ): frame
                for frame in frames
            }
            for future in as_completed(futures):
                outcome = future.result()
                frame = outcome["frame"]
                rows[frame] = outcome["row"]
                outcomes.append(outcome)
                save_manifest(MANIFEST, columns, rows)
                if outcome["status"] == "complete":
                    result = outcome["result"]
                    print(
                        f"[done] frame={frame} roots="
                        f"{result['site_A']['root']}/{result['site_B']['root']} "
                        f"EA={result['site_A']['energy_eV']:.6f} "
                        f"EB={result['site_B']['energy_eV']:.6f} "
                        f"J={result['tdc']['J_screened_cm-1']:.3f} cm-1",
                        flush=True,
                    )
                else:
                    print(
                        f"[failed] frame={frame} class={outcome['failure_class']}: "
                        f"{outcome['error']}",
                        flush=True,
                    )
    finally:
        stop_monitor.set()
        monitor_thread.join(timeout=3.0)
    wall = time.monotonic() - started
    if not samples:
        used, total = gpu_memory()
        samples.append((time.time(), used, total))
    return {
        "frames": frames,
        "workers": workers,
        "roots": roots,
        "profile": profile,
        "wall_seconds": wall,
        "throughput_frames_per_hour": len(frames) * 3600.0 / wall,
        "peak_gpu_memory_mib": max(item[1] for item in samples),
        "gpu_memory_total_mib": samples[-1][2],
        "gpu_memory_peak_fraction": max(item[1] for item in samples) / samples[-1][2],
        "complete": sum(item["status"] == "complete" for item in outcomes),
        "failed": sum(item["status"] == "failed" for item in outcomes),
        "outcomes": [
            {
                key: value
                for key, value in item.items()
                if key not in {"row", "result"}
            }
            for item in outcomes
        ],
    }


def archive_for_retry(frame: int, profile: str, roots: int) -> None:
    frame_dir = ROOT / f"frames/{WINDOW}/frame_{frame:04d}"
    archive = ROOT / (
        f"frames/{WINDOW}/frame_{frame:04d}_attempt0_{profile}_{roots}roots"
    )
    if archive.exists():
        return
    if not frame_dir.exists():
        raise RuntimeError(f"Missing failed frame directory: {frame_dir}")
    frame_dir.rename(archive)


def main() -> None:
    global WINDOW, MANIFEST, OUTDIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu", default="0")
    parser.add_argument("--maximum-workers", type=int, default=4, choices=(1, 2, 3, 4))
    parser.add_argument(
        "--window",
        choices=(
            "long100ns50",
            "short100fs50",
            "short100fs50_fixedmm",
            "short100fs100_wb97xd3",
        ),
        default="long100ns50",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Campaign result directory; defaults from --window",
    )
    args = parser.parse_args()
    WINDOW = args.window
    MANIFEST = ROOT / "manifests" / f"{WINDOW}_frames.csv"
    OUTDIR = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else ROOT
        / (
            "rapid_100fs50_fixedmm"
            if WINDOW == "short100fs50_fixedmm"
            else "wb97xd3_100fs100"
            if WINDOW == "short100fs100_wb97xd3"
            else "rapid_100fs50"
            if WINDOW == "short100fs50"
            else "rapid_100ns50"
        )
    )
    OUTDIR.mkdir(parents=True, exist_ok=True)
    lock_handle = (OUTDIR / "campaign.lock").open("w")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        raise RuntimeError("Another adaptive campaign driver already holds the lock") from error
    lock_handle.write(f"pid={__import__('os').getpid()}\n")
    lock_handle.flush()

    calibration = json.loads((OUTDIR / "calibration.json").read_text())
    delta = float(calibration["steom_constant_shift_eV"])
    reference = Path(calibration["state_tracking_reference"])
    primary_profile = str(calibration.get("production_profile", calibration["profile"]))
    primary_roots = int(calibration.get("production_roots", calibration["roots"]))
    columns, rows = load_manifest(MANIFEST)
    expected_frames = 100 if WINDOW == "short100fs100_wb97xd3" else 50
    if len(rows) != expected_frames:
        raise RuntimeError(
            f"Expected {expected_frames} manifest rows, found {len(rows)}"
        )

    benchmark_path = OUTDIR / "adaptive_parallel_benchmark.json"
    benchmark: dict = {
        "status": "running",
        "policy": (
            "Increase concurrency only when the preceding measured level has no "
            "failures, improves frame throughput by at least 5%, and remains below "
            "90% GPU-memory occupancy."
        ),
        "levels": [],
        "selected_workers": 1,
    }
    last = expected_frames - 1
    midpoint = expected_frames // 2
    pilot_groups = {
        1: [0],
        2: [midpoint - 1, last],
        3: [1, midpoint, last - 1],
        4: [2, midpoint - 2, midpoint + 1, last - 2],
    }
    previous_throughput = None
    selected_workers = 1
    for level in range(1, args.maximum_workers + 1):
        group = [
            frame
            for frame in pilot_groups[level]
            if rows[frame]["convergence_status"] != "complete"
        ]
        if not group:
            continue
        print(f"[adaptive-pilot] workers={level} frames={group}", flush=True)
        result = run_batch(
            group,
            level,
            rows,
            columns,
            args.gpu,
            primary_roots,
            primary_profile,
            delta,
            reference,
        )
        if level in (1, 2):
            for frame in group:
                if rows[frame]["convergence_status"] == "complete":
                    assert_caps(frame)
        throughput = result["throughput_frames_per_hour"]
        gain = (
            None
            if previous_throughput is None
            else throughput / previous_throughput - 1.0
        )
        accepted = (
            result["failed"] == 0
            and result["gpu_memory_peak_fraction"] < 0.90
            and (gain is None or gain >= 0.05)
        )
        result["throughput_gain_over_preceding_level"] = gain
        result["accepted"] = accepted
        benchmark["levels"].append(result)
        if accepted:
            selected_workers = level
            previous_throughput = throughput
            benchmark["selected_workers"] = selected_workers
            benchmark_path.write_text(json.dumps(benchmark, indent=2) + "\n")
        else:
            print(
                f"[adaptive-stop] workers={level} accepted=false "
                f"gain={gain} memory={result['gpu_memory_peak_fraction']:.3f} "
                f"failures={result['failed']}",
                flush=True,
            )
            break

    pending = [
        index
        for index, row in enumerate(rows)
        if row["convergence_status"] not in {"complete", "failed"}
    ]
    print(
        f"[production] selected_workers={selected_workers} pending_frames={len(pending)}",
        flush=True,
    )
    if pending:
        production = run_batch(
            pending,
            selected_workers,
            rows,
            columns,
            args.gpu,
            primary_roots,
            primary_profile,
            delta,
            reference,
        )
        benchmark["production"] = production

    failed_before_retry = [
        index for index, row in enumerate(rows) if row["convergence_status"] == "failed"
    ]
    retry_audit = []
    for frame in failed_before_retry:
        failure_class = rows[frame]["failure_class"]
        archive_for_retry(frame, primary_profile, primary_roots)
        rows[frame]["retry_count"] = str(int(rows[frame].get("retry_count") or 0) + 1)
        if failure_class == "missing_target_state_or_root_flipping":
            roots = 9 if primary_roots >= 7 else 7
            profile = primary_profile
        elif failure_class in {
            "ground_state_scf_convergence",
            "excited_state_eigensolver_convergence",
        }:
            roots = primary_roots
            profile = (
                "rapid_cam_retry"
                if primary_profile == "rapid_cam"
                else "rapid_retry"
            )
        else:
            # A serial retry is the minimally stronger response to a concurrent
            # GPU-memory/numerical failure and also regenerates malformed input.
            roots, profile = primary_roots, primary_profile
        print(
            f"[retry] frame={frame} original_class={failure_class} "
            f"profile={profile} roots={roots} workers=1",
            flush=True,
        )
        outcome = run_batch(
            [frame],
            1,
            rows,
            columns,
            args.gpu,
            roots,
            profile,
            delta,
            reference,
        )
        retry_audit.append(
            {
                "frame": frame,
                "original_failure_class": failure_class,
                "retry_profile": profile,
                "retry_roots": roots,
                "result": outcome,
            }
        )

    final_complete = sum(row["convergence_status"] == "complete" for row in rows)
    final_failed = sum(row["convergence_status"] == "failed" for row in rows)
    unresolved = [
        int(row["frame_index"])
        for row in rows
        if row["convergence_status"] not in {"complete", "failed"}
    ]
    if unresolved:
        raise RuntimeError(f"Unclassified frames after campaign: {unresolved}")
    benchmark.update(
        {
            "status": "complete",
            "selected_workers": selected_workers,
            "failed_before_retry": failed_before_retry,
            "retry_audit": retry_audit,
            "final_complete": final_complete,
            "final_failed": final_failed,
            "failed_frame_indices": [
                int(row["frame_index"])
                for row in rows
                if row["convergence_status"] == "failed"
            ],
        }
    )
    benchmark_path.write_text(json.dumps(benchmark, indent=2) + "\n")
    print(
        f"[campaign-complete] complete={final_complete}/{expected_frames} "
        f"failed={final_failed}/{expected_frames}",
        flush=True,
    )


if __name__ == "__main__":
    main()
