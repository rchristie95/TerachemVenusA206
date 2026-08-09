#!/usr/bin/env python3
"""Orchestrate the staged TDDFT/STEOM site-energy CD campaign."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
import numpy as np


def run_command(cmd, env=None):
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, env=env)


def generate_frames(frame_indices, out_prefix, args, tddft_in="tddft.in"):
    for idx in frame_indices:
        out_dir = args.results_dir / f"{out_prefix}_frame_{idx:04d}"
        if not out_dir.exists():
            cmd = [
                args.python,
                str(args.scripts_dir / "prepare_production_frame.py"),
                "--topology", str(args.topology),
                "--trajectory", str(args.trajectory),
                "--frame-index", str(idx),
                "--output-dir", str(out_dir),
                "--embedding-cache", str(args.results_dir / "embedding_charges.npz"),
                "--amber-cr2-prmtop", str(args.amber_cr2_prmtop),
                "--retain-partner-cr2-charges",
                "--conserve-boundary-residue-charge"
            ]
            if out_prefix != "pilot" or idx != 0:
                cmd.extend(["--fixed-embedding-selection", str(args.results_dir / "fixed_embedding_selection.json")])
            else:
                # Frame 0 of pilot creates the fixed embedding selection
                cmd.extend(["--fixed-embedding-selection", str(args.results_dir / "fixed_embedding_selection.json")])
            
            run_command(cmd)
        else:
            print(f"Skipping preparation of {out_dir}, already exists.")


def launch_frames(frame_indices, out_prefix, args):
    env = os.environ.copy()
    if args.gpu:
        env["CUDA_VISIBLE_DEVICES"] = args.gpu

    for idx in frame_indices:
        out_dir = args.results_dir / f"{out_prefix}_frame_{idx:04d}"
        cmd = [
            args.python,
            str(args.scripts_dir / "launch_jobs.py"),
            str(out_dir),
            "--terachem", str(args.terachem),
            "--gpu", "0"
        ]
        try:
            run_command(cmd, env=env)
        except subprocess.CalledProcessError:
            print(f"Job failed for frame {idx}. Check logs in {out_dir}.")
            sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--amber-cr2-prmtop", type=Path, required=True)
    parser.add_argument("--terachem", type=Path, required=True)
    parser.add_argument("--python", type=str, default=sys.executable)
    parser.add_argument("--gpu", type=str, default="0")
    parser.add_argument("--scripts-dir", type=Path, default=Path(__file__).parent)
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).parent.parent / "results")
    
    parser.add_argument("--pilot", action="store_true", help="Run 5-frame state-tracking pilot")
    parser.add_argument("--feasibility", action="store_true", help="Run 50-frame feasibility ensemble")
    parser.add_argument("--calibration", action="store_true", help="Prepare 12-frame STEOM calibration subset")
    parser.add_argument("--production", action="store_true", help="Run 1000-frame production (if authorized)")
    
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Run audit first
    audit_cmd = [
        args.python,
        str(args.scripts_dir / "audit_inputs.py"),
        "--topology", str(args.topology),
        "--trajectory", str(args.trajectory),
        "--amber-cr2-prmtop", str(args.amber_cr2_prmtop),
        "--terachem", str(args.terachem),
        "--output", str(args.results_dir / "run_manifest.json")
    ]
    run_command(audit_cmd)
    
    with open(args.results_dir / "run_manifest.json") as f:
        manifest = json.load(f)
        if not manifest.get("production_join_allowed", False):
            print("Audit failed! Fix hash mismatch before continuing.")
            sys.exit(1)
            
    if args.pilot:
        # 5 stratified frames: 0, 249, 499, 749, 999
        pilot_frames = [0, 249, 499, 749, 999]
        print("=== Step 2: Preparing and launching 5 stratified pilot frames ===")
        # Wait, the TDDFT input template was already modified to include NTO tracking
        generate_frames(pilot_frames, "pilot", args)
        launch_frames(pilot_frames, "pilot", args)
        print("Pilot jobs launched. Validate transition-density/NTO overlap manually before proceeding.")
        
    if args.feasibility:
        # 50 matched frames: every 20 frames
        feasibility_frames = list(range(0, 1000, 20))
        print(f"=== Step 3: Preparing and launching {len(feasibility_frames)} feasibility frames ===")
        generate_frames(feasibility_frames, "feasibility", args)
        launch_frames(feasibility_frames, "feasibility", args)
        
    if args.calibration:
        # 12 representative site/frame points with STEOM
        calib_frames = np.linspace(0, 999, 12, dtype=int).tolist()
        print(f"=== Step 4: Preparing {len(calib_frames)} frames for STEOM calibration ===")
        generate_frames(calib_frames, "calibration", args)
        print("Frames prepared. Run STEOM calculations separately.")
        
    if args.production:
        print("=== Step 5: Preparing and launching 1000-frame production ===")
        prod_frames = list(range(1000))
        generate_frames(prod_frames, "production", args)
        launch_frames(prod_frames, "production", args)


if __name__ == "__main__":
    main()
