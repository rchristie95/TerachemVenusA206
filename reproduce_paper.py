#!/usr/bin/env python3
"""
reproduce_paper.py — one-shot reproduction of the numerical data in the
Venus_A206 excitonic-coupling paper.

Chains the existing tools into a single, reproducible, reuse-aware orchestrator:

  TDDFT  (TeraChem, GPU)      site energy + transition density   [single-excitation reference]
  STEOM-CCSD (ORCA, CPU)      in-protein bright state (~532 nm)
  EOM-CCSD(fT)/ADC(2) (Q-Chem) triples energy + doubles character [the doubles/triples that
                               TDDFT, being single-excitation, structurally cannot show]
  Davydov coupling J (TDC)    static (single geometry) + thermal NVT ensemble, for both
                               the TDDFT and STEOM transition densities

Design:
  * Config block below (SEED, EPS, N_FRAMES, COUPLING_BACKEND, REUSE toggles).
  * Each stage reuses cached outputs by default (REUSE[...]=True) and only recomputes
    when toggled off (or via the matching --run-* flag). Recompute shells out to the
    right environment (TeraChem env for OpenMM/TeraChem/coupling; the openmpi416+ORCA
    env via go_par.sh for STEOM; Q-Chem at $HOME/qchem for EOM/ADC).
  * Coupling runs on the GPU via the OpenCL backend (numba-CUDA is PTX-blocked on this box).
  * Aggregates everything to paper_data_summary.json + a printed table.

Honest-provenance notes baked in:
  * The doubles/triples character is a STEOM/EOM-CCSD property; the script reports the
    method ladder (TDDFT -> bare EOM-CCSD -> EOM-CCSD(fT) ~ STEOM), not a parsed "% doubles".
  * The published TDDFT thermal coupling (65.3) is read from its cached distribution; a
    live recompute of it is a documented open item (see STEOM_COUPLING_FINDINGS, sec 2.5).
    The STEOM thermal coupling is recomputed live (reproduces ~96.4).
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# ============================================================
#  CONFIG  (edit here; every value is overridable on the CLI)
# ============================================================
SEED = 20260618                 # global reproducibility seed (NVT integrator, frame choice)
EPS = 1.77                      # optical screening (matches the published TDDFT ensemble)
N_FRAMES = 200                  # NVT frames for the thermal coupling ensemble
COUPLING_BACKEND = "opencl"     # GPU backend for the TDC kernel (numba-CUDA is PTX-blocked here)

REUSE = {                       # True = reuse cached outputs if present; False = recompute
    "nvt":     True,            # reuse dimer_nvt_restrained_clean.pdb
    "tddft":   True,            # reuse tc_tddft_*/energy.out + transition density
    "steom":   True,            # reuse neo_model/orca_steom/steom_phenol_svpd.out
    "eomft":   True,            # reuse qchem_validation/eomcc_ft_*.out + eomccsd_bare/adc2
    "density": True,            # reuse the built/spec-normalised/matched STEOM density npz
}

# ----- paths (relative to this file) -----
REPO = Path(__file__).resolve().parent
PY = sys.executable                              # TeraChem-env python running this script

TRAJ          = "dimer_nvt_restrained_clean.pdb"     # restrained NVT, coupling-ready
DIMER         = "venus_dimer.pdb"                    # crystal dimer (static-geometry J)
OLD_MONOMER   = "tc_simple_old/classical_relaxed.pdb"  # frame the dimer chains were built in
ANION_MONOMER = "tc_simple_anionic/monomer_relaxed.pdb"

STEOM_DIR      = REPO / "neo_model/orca_steom"
STEOM_SPECNORM = STEOM_DIR / "steom_transdens_specnorm.npz"
STEOM_MATCHED  = STEOM_DIR / "steom_transdens_specnorm_oldframe.npz"
STEOM_OUT      = STEOM_DIR / "steom_phenol_svpd.out"
STEOM_INP      = "steom_phenol_svpd.inp"             # arg to go_par.sh (run from STEOM_DIR)
GO_PAR         = STEOM_DIR / "go_par.sh"

TDDFT_DIRS = ["tc_tddft_old_current_current", "tc_tddft_prod_current",
              "tc_tddft_anionic_current", "tc_tddft_44"]
TDDFT_THERMAL_JSON = "coupling_sampling_out/coupling_distribution.json"  # cached 65.3

QCHEM_DIR   = REPO / "qchem_validation"
QCHEM_EOMFT = "eomcc_ft_631g.out"
QCHEM_BARE  = "eomccsd_bare.out"
QCHEM_ADC2  = "adc2_bare.out"

OUT_DIR_STEOM_STATIC  = "coupling_paper_steom_static"
OUT_DIR_STEOM_THERMAL = "coupling_paper_steom_thermal"
SUMMARY = "paper_data_summary.json"

LOG_DIR = REPO / "pipeline_logs"

# Reference numbers (for provenance/labels only; never overwrite a freshly parsed value)
PUBLISHED_TDDFT_STATIC_J = 74.38      # paper, single minimised geometry


# ============================================================
#  small utilities
# ============================================================
def log(msg):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def run(cmd, log_name, cwd=None, env=None):
    """Run a subprocess, tee to a log file, return (returncode, tail)."""
    LOG_DIR.mkdir(exist_ok=True)
    logf = LOG_DIR / log_name
    log(f"  $ {' '.join(str(c) for c in cmd)}   (log: {logf.name})")
    with open(logf, "w") as fh:
        proc = subprocess.run(cmd, cwd=cwd, env=env, stdout=fh,
                              stderr=subprocess.STDOUT, text=True)
    tail = "\n".join(logf.read_text(errors="replace").splitlines()[-15:])
    return proc.returncode, tail


# ============================================================
#  parsers (reuse path)
# ============================================================
def parse_orca_steom_spectrum(out_path):
    """Bright (max-fosc) row of the ORCA STEOM absorption spectrum -> dict or None."""
    p = Path(out_path)
    if not p.exists():
        return None
    lines = p.read_text(errors="replace").splitlines()
    start = None
    for i, ln in enumerate(lines):
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE" in ln:
            start = i
    if start is None:
        return None
    rows = []
    for ln in lines[start + 4: start + 30]:
        m = re.match(r"\s*\d+-\d+\w+\s+->\s+\d+-\d+\w+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
                     r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", ln)
        if not m:
            if rows:
                break
            continue
        ev, cm, nm, fosc, d2, dx, dy, dz = map(float, m.groups())
        rows.append({"ev": ev, "cm": cm, "nm": nm, "fosc": fosc,
                     "mu_au": (dx**2 + dy**2 + dz**2) ** 0.5})
    if not rows:
        return None
    bright = max(rows, key=lambda r: r["fosc"])
    txt = "\n".join(lines)
    # Success criterion for STEOM here is "the converged spectrum printed", NOT clean
    # termination: this calc reliably error-terminates in MDCI (the DoSTEOMNatTransOrb
    # post-step) *after* emitting the valid spectrum. Treat that as a converged result.
    bright["terminated_normally"] = "TERMINATED NORMALLY" in txt
    bright["mdci_error_after_spectrum"] = "error termination in MDCI" in txt
    bright["spectrum_converged"] = True  # we parsed a STEOM absorption block
    return bright


def parse_qchem(out_path):
    """EOM/ADC excitation energies (eV) and the ground-state CCSD T1^2/T2^2 from a Q-Chem out.

    NB: T2^2 is the GROUND-STATE CC doubles-amplitude norm (always large), NOT the
    excited-state double-excitation character. The excited-state 2p2h weight (the ~12%
    that matters here) comes from ADC(2); these cached files are heterogeneous (6-31G
    exploratory) and one ADC(2) run was truncated by a disk-full event, so we report the
    raw parsed energies with provenance rather than asserting a single bright value.
    """
    p = Path(out_path)
    if not p.exists():
        return {"file": str(p.name), "present": False}
    txt = p.read_text(errors="replace")
    evs = [float(x) for x in re.findall(r"Excitation energy\s*=\s*([\d.]+)\s*eV", txt)]
    t = re.search(r"T1\^2\s*=\s*([\d.]+)\s+T2\^2\s*=\s*([\d.]+)", txt)
    return {
        "file": p.name,
        "present": True,
        "excitation_eV": [round(e, 4) for e in evs],   # all parsed roots, in print order
        "lowest_eV": (min(evs) if evs else None),       # NOT necessarily the bright state
        "gs_ccsd_T1sq": (float(t.group(1)) if t else None),
        "gs_ccsd_T2sq": (float(t.group(2)) if t else None),
        "complete": ("Thank you very much for using Q-Chem" in txt),
    }


def parse_terachem_all(dirs):
    """Brightest (max-fosc) state from EVERY cached tc_tddft_* dir with an energy.out.

    Reports each geometry separately: the old generic-FF geometry (~549 nm, the
    error-cancellation case) and the physically-correct anionic geometry (~420 nm, blue)
    are different data points and both belong in the paper's TDDFT comparison.
    """
    try:
        from coupling_core import parse_excited_state_candidates
    except Exception as e:
        return [{"error": f"coupling_core import failed: {e}"}]
    out = []
    for d in dirs:
        cands = parse_excited_state_candidates(REPO / d / "energy.out")
        if not cands:
            continue
        vis = [c for c in cands if 450 <= c["nm"] <= 650] or cands
        b = max(vis, key=lambda c: c["osc"])
        out.append({"dir": d, "root": b["root"], "ev": round(b["ev"], 4),
                    "nm": round(b["nm"], 1), "fosc": round(b["osc"], 4)})
    return out or [{"error": f"no energy.out with parsed states in {dirs}"}]


def coupling_summary(out_dir):
    """Read a coupling_ensemble.py output distribution.json."""
    j = REPO / out_dir / "coupling_distribution.json"
    if not j.exists():
        return None
    s = json.loads(j.read_text())
    return {"J_mean": s["mean"], "J_std": s.get("std", 0.0),
            "two_J": 2 * abs(s["mean"]), "n": s["n"], "eps": s.get("epsilon")}


# ============================================================
#  stages
# ============================================================
def preflight(args):
    info = {"gpu_opencl": False, "missing": []}
    try:
        from coupling_core import _is_opencl_ready, _is_cuda_ready
        info["gpu_opencl"] = bool(_is_opencl_ready())
        info["gpu_cuda_numba"] = bool(_is_cuda_ready())
    except Exception as e:
        info["coupling_core_error"] = str(e)
    for f in [TRAJ, DIMER, OLD_MONOMER, ANION_MONOMER, STEOM_SPECNORM]:
        if not (REPO / f).exists():
            info["missing"].append(str(f))
    if not info["gpu_opencl"]:
        log("  [!] OpenCL GPU not available — coupling stages need it (numba-CUDA is "
            "PTX-blocked). If the driver is down, rebuild it (see memory: "
            "gpu-driver-rebuild-after-kernel-update).")
    if info["missing"]:
        log(f"  [!] missing required inputs: {info['missing']}")
    return info


def stage_tddft(args):
    if not args.reuse["tddft"]:
        log("  TDDFT recompute requested — invoking qmmm_tddft_pipeline stage2 (GPU).")
        run([PY, "qmmm_tddft_pipeline.py", "--skip-simple", "--skip-coupling",
             "--skip-visualize"], "tddft_stage2.log", cwd=REPO)
    return {"site_energies": parse_terachem_all(TDDFT_DIRS),
            "note": "single-excitation reference; misses the doubles/triples character"}


def stage_steom(args):
    if not args.reuse["steom"]:
        # Protect the authoritative result: go_par.sh truncates the .out at launch, and
        # this calc reliably aborts in MDCI (DoSTEOMNatTransOrb) AFTER printing the valid
        # spectrum. If a re-run dies before the spectrum, restore the cache so the pipeline
        # still finishes end-to-end (critical for a publishable, reproducible run).
        bak = STEOM_DIR / "_authoritative_backup"
        bak.mkdir(exist_ok=True)
        for f in list(STEOM_DIR.glob("steom_phenol_svpd.out")) + \
                 list(STEOM_DIR.glob("steom_phenol_svpd.gbw")) + \
                 list(STEOM_DIR.glob("steom_phenol_svpd.s1.*.cube")):
            shutil.copy2(f, bak / f.name)
        log("  STEOM recompute requested — launching ORCA DLPNO-STEOM-CCSD (CPU, hours).")
        log("  (expected: prints the 532.6 spectrum, then a benign MDCI/NatTransOrb abort)")
        run(["bash", str(GO_PAR), STEOM_INP], "steom_run.log", cwd=STEOM_DIR)
        bright = parse_orca_steom_spectrum(STEOM_OUT)
        if not bright:
            log("  [!] re-run produced no STEOM spectrum — restoring authoritative backup.")
            for f in bak.glob("steom_phenol_svpd.*"):
                shutil.copy2(f, STEOM_DIR / f.name)
            return {"bright": parse_orca_steom_spectrum(STEOM_OUT),
                    "recompute_failed_restored": True,
                    "source": str(STEOM_OUT.relative_to(REPO))}
        return {"bright": bright, "recomputed": True,
                "source": str(STEOM_OUT.relative_to(REPO))}
    return {"bright": parse_orca_steom_spectrum(STEOM_OUT),
            "source": str(STEOM_OUT.relative_to(REPO))}


def stage_eom_triples(args):
    if not args.reuse["eomft"]:
        log("  EOM-CCSD(fT)/ADC(2) recompute requested — Q-Chem (qchem_validation/).")
        log("  [i] regenerate inputs with qchem_validation/make_qchem_inputs.py then run "
            "qchem; skipping automatic launch (engine/queue specific).")
    return {
        "raw_qchem": {
            "eom_ccsd_fT": parse_qchem(QCHEM_DIR / QCHEM_EOMFT),
            "eom_ccsd_bare": parse_qchem(QCHEM_DIR / QCHEM_BARE),
            "adc2": parse_qchem(QCHEM_DIR / QCHEM_ADC2),
        },
        # Established def2-SVP bare-anion comparison (reference; see memory
        # steom-vs-eomccsd-validation). Reproduce live with --run-eomft (Q-Chem, hours).
        "validated_ladder_eV": {
            "EOM_CCSD_bare": 3.72,   # too blue — misses the triples
            "EOM_CCSD_fT": 3.29,     # (fT) triples correction, ~-0.43 eV
            "STEOM_CCSD": 3.335,     # ~ EOM-CCSD(fT) -> STEOM validated
            "ADC2_2p2h_doubles_pct": 12,   # excited-state double-excitation character
        },
        "note": ("Doubles/triples is a STEOM/EOM-CCSD property; TDDFT (single-excitation) "
                 "cannot show it. Evidence: ADC(2) ~12% 2p2h doubles character, and the (fT) "
                 "triples correction (-0.43 eV) brings bare EOM-CCSD (3.72) onto STEOM (3.335). "
                 "raw_qchem lists are the cached (heterogeneous-basis) runs; validated_ladder "
                 "is the established def2-SVP comparison."),
    }


def stage_density(args):
    if args.reuse["density"] and STEOM_MATCHED.exists():
        log(f"  reusing matched STEOM density {STEOM_MATCHED.name}")
        return {"matched_density": str(STEOM_MATCHED.relative_to(REPO)), "rebuilt": False}
    log("  building matched STEOM density (Kabsch into the dimer-chain frame)")
    rc, tail = run([PY, "align_steom_density.py",
                    "--density", str(STEOM_SPECNORM),
                    "--anion-pdb", ANION_MONOMER, "--old-pdb", OLD_MONOMER,
                    "--out", str(STEOM_MATCHED)], "match_density.log", cwd=REPO)
    return {"matched_density": str(STEOM_MATCHED.relative_to(REPO)),
            "rebuilt": True, "rc": rc}


def _sample_coupling(out_dir, traj, monomer, density, n_frames, args, random=False):
    cmd = [PY, "coupling_ensemble.py", "--traj", traj, "--monomer", monomer,
           "--density", str(density), "--mode", "rigid",
           "--backend", args.backend, "--epsilon", str(args.eps),
           "--n-frames", str(n_frames), "--out", out_dir]
    if random:
        cmd += ["--random", "--seed", str(args.seed)]
    rc, tail = run(cmd, f"{out_dir}.log", cwd=REPO)
    return rc, tail


def stage_static_J(args):
    # STEOM: single-geometry coupling on the crystal dimer (live, reproducible).
    _sample_coupling(OUT_DIR_STEOM_STATIC, DIMER, OLD_MONOMER, STEOM_MATCHED, 1, args)
    return {
        "STEOM": coupling_summary(OUT_DIR_STEOM_STATIC),
        "TDDFT_published": {"J": PUBLISHED_TDDFT_STATIC_J, "two_J": 2 * PUBLISHED_TDDFT_STATIC_J,
                            "note": "paper single-min geometry; same kernel on the 44-atom "
                                    "spectroscopy geom gives ~118 (open item, findings 2.5)"},
    }


def stage_thermal_J(args):
    # STEOM thermal ensemble: recomputed live on the restrained NVT trajectory.
    _sample_coupling(OUT_DIR_STEOM_THERMAL, TRAJ, OLD_MONOMER, STEOM_MATCHED,
                     args.n_frames, args, random=(args.n_frames < 200))
    steom = coupling_summary(OUT_DIR_STEOM_THERMAL)
    # TDDFT thermal ensemble: read the cached published distribution (65.3).
    tddft = coupling_summary(Path(TDDFT_THERMAL_JSON).parent.name) if \
        (REPO / TDDFT_THERMAL_JSON).exists() else None
    return {
        "STEOM": steom,
        "TDDFT_cached": tddft,
        "note": ("QM-dipole-normalised (NOT the empirical 7.3 D). TDDFT value is the cached "
                 "published ensemble; live TDDFT recompute is the documented open item."),
    }


# ============================================================
#  main
# ============================================================
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--eps", type=float, default=EPS)
    p.add_argument("--n-frames", type=int, default=N_FRAMES)
    p.add_argument("--backend", default=COUPLING_BACKEND, choices=["opencl", "gpu", "auto"])
    # --run-X flips REUSE[X] to False (force recompute of that stage)
    for k in REUSE:
        p.add_argument(f"--run-{k}", action="store_true", help=f"recompute the {k} stage")
    p.add_argument("--out", default=SUMMARY)
    a = p.parse_args(argv)
    a.reuse = {k: (not getattr(a, f"run_{k}")) for k in REUSE}
    return a


def main(argv=None):
    args = parse_args(argv)
    log("=" * 60)
    log("Venus_A206 paper-data pipeline")
    log(f"  seed={args.seed} eps={args.eps} n_frames={args.n_frames} backend={args.backend}")
    log(f"  reuse={args.reuse}")
    log("=" * 60)

    results = {"config": {"seed": args.seed, "eps": args.eps, "n_frames": args.n_frames,
                          "backend": args.backend, "reuse": args.reuse},
               "timestamp": datetime.now().isoformat()}

    log("[0/7] preflight");                results["preflight"] = preflight(args)
    log("[1/7] TDDFT site energy (GPU)");  results["tddft"] = stage_tddft(args)
    log("[2/7] STEOM-CCSD (CPU)");         results["steom"] = stage_steom(args)
    log("[3/7] EOM-CCSD(fT)/ADC(2) doubles+triples (Q-Chem)"); \
        results["doubles_triples"] = stage_eom_triples(args)
    log("[4/7] STEOM coupling density");   results["density"] = stage_density(args)
    log("[5/7] static Davydov J (TDC)");   results["static_J"] = stage_static_J(args)
    log("[6/7] thermal NVT J ensemble");   results["thermal_J"] = stage_thermal_J(args)

    log("[7/7] aggregate")
    (REPO / args.out).write_text(json.dumps(results, indent=2))
    _print_table(results)
    log(f"wrote {args.out}")
    return 0


def _g(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


def _print_table(r):
    print("\n" + "=" * 64)
    print("  VENUS_A206 PAPER DATA — SUMMARY")
    print("=" * 64)
    for td in (_g(r, "tddft", "site_energies", default=[]) or []):
        if "error" in td:
            print(f"  TDDFT  S1 (GPU)        : {td['error']}")
        else:
            print(f"  TDDFT  S1 (GPU)        : {td['nm']} nm  f={td['fosc']}   [{td['dir']}]")
    st = _g(r, "steom", "bright")
    print(f"  STEOM-CCSD S1 (CPU)    : {_g(st,'nm',default='?')} nm  f={_g(st,'fosc',default='?')}"
          f"  |mu|={_g(st,'mu_au',default='?')} au")
    lad = _g(r, "doubles_triples", "validated_ladder_eV", default={})
    print(f"  Doubles/triples ladder : bare EOM-CCSD {lad.get('EOM_CCSD_bare','?')} -> "
          f"EOM-CCSD(fT) {lad.get('EOM_CCSD_fT','?')} ~ STEOM {lad.get('STEOM_CCSD','?')} eV "
          f"(ADC(2) {lad.get('ADC2_2p2h_doubles_pct','?')}% 2p2h)")
    ftraw = _g(r, "doubles_triples", "raw_qchem", "eom_ccsd_fT", "excitation_eV")
    print(f"  (cached Q-Chem fT roots: {ftraw}  — heterogeneous basis; reference ladder above)")
    sj = _g(r, "static_J", "STEOM")
    print(f"  Static J  STEOM        : {_g(sj,'J_mean',default='?')} cm^-1  (2|J|={_g(sj,'two_J',default='?')})")
    print(f"  Static J  TDDFT (pub)  : {_g(r,'static_J','TDDFT_published','J')} cm^-1")
    tj = _g(r, "thermal_J", "STEOM")
    tjt = _g(r, "thermal_J", "TDDFT_cached")
    print(f"  Thermal J STEOM (NVT)  : {_g(tj,'J_mean',default='?')} +/- {_g(tj,'J_std',default='?')}"
          f"  (2|J|={_g(tj,'two_J',default='?')}, n={_g(tj,'n',default='?')})")
    print(f"  Thermal J TDDFT (cache): {_g(tjt,'J_mean',default='?')} +/- {_g(tjt,'J_std',default='?')}"
          f"  (2|J|={_g(tjt,'two_J',default='?')})")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
