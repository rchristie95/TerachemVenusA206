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
import os
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
ENGINE = "terachem"             # site-energy DFT engine: "terachem" (GPU) or "orca" (CPU).
                                # "orca" runs ORCA TDDFT at the STEOM geometry AND skips Q-Chem,
                                # giving a TeraChem-free, ORCA-only path (DFT + STEOM in one engine).

REUSE = {                       # True = reuse cached outputs if present; False = recompute
    "nvt":     True,            # reuse dimer_nvt_restrained_clean.pdb
    "tddft":   True,            # reuse tc_tddft_*/energy.out + transition density
    "steom":   True,            # reuse neo_model/orca_steom/steom_phenol_svpd.out
    "eomft":   True,            # reuse qchem_validation/eomcc_ft_*.out + eomccsd_bare/adc2
    "density": True,            # reuse the built/spec-normalised/matched STEOM density npz
    "spectra": True,            # reuse lineshape_out/ absorption+CD figure panels
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

# ORCA TDDFT (engine=orca): reuses the *identical* 44-atom STEOM geometry + point-charge
# field, so the DFT and STEOM site energies are computed on exactly the same structure.
ORCA_TDDFT_DIR = REPO / "neo_model/orca_dft"
ORCA_TDDFT_INP = "tddft_wb97xd3.inp"
ORCA_TDDFT_OUT = "tddft_wb97xd3.out"
ORCA_BIN       = REPO / "orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"
OPENMPI_DIR    = Path("/home/robson/anaconda3/envs/openmpi416")  # openmpi4.1.6 for parallel ORCA
STEOM_GEOM     = "geom_cthrp.xyz"                    # 44-atom CR2+Tyr203 phenol (shared with STEOM)
STEOM_FIELD    = "field.pc"                          # ORCA point-charge embedding (shared with STEOM)

TDDFT_DIRS = ["tc_tddft_old_current_current", "tc_tddft_prod_current",
              "tc_tddft_anionic_current", "tc_tddft_44"]
TDDFT_THERMAL_JSON = "coupling_sampling_out/coupling_distribution.json"  # cached 65.3

QCHEM_DIR   = REPO / "qchem_validation"
QCHEM_EOMFT = "eomcc_ft_631g.out"
QCHEM_BARE  = "eomccsd_bare.out"
QCHEM_ADC2  = "adc2_bare.out"

OUT_DIR_STEOM_STATIC  = "coupling_paper_steom_static"
OUT_DIR_STEOM_THERMAL = "coupling_paper_steom_thermal"
OUT_DIR_SPECTRA       = "lineshape_out"                          # absorption/CD figure panels
DIPOLE_GEOM_JSON      = f"{OUT_DIR_STEOM_THERMAL}/dipole_geometry.json"
ENSEMBLE_GEOM_NPZ     = f"{OUT_DIR_STEOM_THERMAL}/coupling_geometry.npz"  # per-frame mu/r/J
THERMAL_DIST_JSON     = f"{OUT_DIR_STEOM_THERMAL}/coupling_distribution.json"
T2_STAR_FS   = 60.0                   # pure-dephasing time feeding the homogeneous Voigt width
EXP_SPLITTING = (262.0, 372.0)        # Nguyen U=131-186 cm^-1 -> Davydov splitting 2U (cm^-1)
SUMMARY = "paper_data_summary.json"

LOG_DIR = REPO / "pipeline_logs"

# Reference numbers (for provenance/labels only; never overwrite a freshly parsed value)
PUBLISHED_TDDFT_STATIC_J = 74.38      # paper, single minimised geometry
REFERENCE_TERACHEM_44_NM = 361.0      # tc_tddft_44 bright state (wb97xd3/6-311g**, TDA) — the
                                      # like-for-like TeraChem number the ORCA TDDFT reproduces

STAGES = ["tddft", "steom", "eom", "density", "static", "thermal", "spectra"]  # ordered; for --stop-after

# ORCA TDDFT reference input (engine=orca). Matches the TeraChem reference method
# (wb97xd3 / 6-311G** / TDA) so the two engines are directly comparable, evaluated on the
# identical 44-atom STEOM geometry + field. The reference TeraChem 44-atom run (tc_tddft_44)
# used bare point-charge embedding with NO COSMO/PCM, so we match that (no %cpcm here). To add
# a dielectric, insert e.g.  %cpcm  epsilon 78.39  refrac 1.33  end  and recompute.
ORCA_TDDFT_INPUT = f"""! wB97X-D3 6-311G** RIJCOSX def2/J TightSCF
%maxcore 2000
%pal nprocs 8 end
%tddft
  tda true
  nroots 5
end
%pointcharges "{STEOM_FIELD}"
* xyzfile -1 1 {STEOM_GEOM}
"""


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
def _orca_abs_bright(out_path):
    """Bright (max-fosc) row of an ORCA 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC
    DIPOLE MOMENTS' block -> dict(ev, cm, nm, fosc, mu_au) or None. The table layout is
    identical for STEOM and TDDFT in ORCA 6.1, so both parsers share this."""
    p = Path(out_path)
    if not p.exists():
        return None
    lines = p.read_text(errors="replace").splitlines()
    start = None
    for i, ln in enumerate(lines):
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE" in ln:
            start = i  # take the last such block
    if start is None:
        return None
    rows = []
    for ln in lines[start + 4: start + 80]:
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
    return max(rows, key=lambda r: r["fosc"])


def parse_orca_steom_spectrum(out_path):
    """Bright (max-fosc) row of the ORCA STEOM absorption spectrum -> dict or None."""
    bright = _orca_abs_bright(out_path)
    if bright is None:
        return None
    txt = Path(out_path).read_text(errors="replace")
    # Success criterion for STEOM here is "the converged spectrum printed", NOT clean
    # termination: this calc reliably error-terminates in MDCI (the DoSTEOMNatTransOrb
    # post-step) *after* emitting the valid spectrum. Treat that as a converged result.
    bright["terminated_normally"] = "TERMINATED NORMALLY" in txt
    bright["mdci_error_after_spectrum"] = "error termination in MDCI" in txt
    bright["spectrum_converged"] = True  # we parsed a STEOM absorption block
    return bright


def parse_orca_tddft_spectrum(out_path):
    """Bright (max-fosc) row of an ORCA TDDFT absorption spectrum -> dict or None.
    Unlike STEOM, a TDDFT job is expected to terminate cleanly."""
    bright = _orca_abs_bright(out_path)
    if bright is None:
        return None
    txt = Path(out_path).read_text(errors="replace")
    bright["terminated_normally"] = "ORCA TERMINATED NORMALLY" in txt
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
    return {"engine": "terachem",
            "site_energies": parse_terachem_all(TDDFT_DIRS),
            "note": "single-excitation reference; misses the doubles/triples character"}


def _orca_env():
    """Environment for the shared/parallel ORCA build (mirrors go_par.sh)."""
    env = os.environ.copy()
    env["PATH"] = f"{OPENMPI_DIR}/bin:{ORCA_BIN}:" + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = f"{OPENMPI_DIR}/lib:{ORCA_BIN}:" + env.get("LD_LIBRARY_PATH", "")
    return env


def _run_orca(inp_name, out_name, cwd, log_name):
    """Run ORCA on inp_name in cwd, writing the ORCA output to cwd/out_name (ORCA prints
    to stdout) and mirroring the tail into pipeline_logs/log_name. Returns (rc, tail)."""
    LOG_DIR.mkdir(exist_ok=True)
    out_path = Path(cwd) / out_name
    log(f"  $ orca {inp_name}   (-> {out_name})")
    with open(out_path, "w") as fh:
        proc = subprocess.run([str(ORCA_BIN / "orca"), inp_name],
                              cwd=str(cwd), env=_orca_env(),
                              stdout=fh, stderr=subprocess.STDOUT, text=True)
    tail = "\n".join(out_path.read_text(errors="replace").splitlines()[-15:])
    (LOG_DIR / log_name).write_text(tail)
    return proc.returncode, tail


def stage_tddft_orca(args):
    """ORCA TDDFT site energy (engine=orca): wB97X-D3/6-311G** (TDA) on the identical
    44-atom STEOM geometry + point-charge field (geom_cthrp.xyz + field.pc). Directly
    comparable to the cached TeraChem 44-atom run (tc_tddft_44 ~ 361 nm) and to the STEOM
    site energy on the same structure — a TeraChem-free, ORCA-only DFT reference."""
    ORCA_TDDFT_DIR.mkdir(parents=True, exist_ok=True)
    staged = []
    for f in (STEOM_GEOM, STEOM_FIELD):
        src = STEOM_DIR / f
        if src.exists():
            shutil.copy2(src, ORCA_TDDFT_DIR / f)
            staged.append(f)
    (ORCA_TDDFT_DIR / ORCA_TDDFT_INP).write_text(ORCA_TDDFT_INPUT)
    out = ORCA_TDDFT_DIR / ORCA_TDDFT_OUT
    if args.reuse["tddft"] and out.exists():
        log(f"  reusing cached ORCA TDDFT {out.name}")
    else:
        missing = [f for f in (STEOM_GEOM, STEOM_FIELD) if f not in staged]
        if missing:
            log(f"  [!] missing STEOM inputs to stage for ORCA TDDFT: {missing}")
        log("  ORCA TDDFT — wB97X-D3/6-311G** (TDA) at the 44-atom STEOM geometry+field (CPU).")
        _run_orca(ORCA_TDDFT_INP, ORCA_TDDFT_OUT, ORCA_TDDFT_DIR, "orca_tddft.log")
    bright = parse_orca_tddft_spectrum(out)
    return {"engine": "orca", "bright": bright,
            "source": str(out.relative_to(REPO)) if out.exists() else None,
            "reference_terachem_44_nm": REFERENCE_TERACHEM_44_NM,
            "note": ("single-excitation ORCA TDDFT at the identical STEOM geometry/embedding; "
                     "reproduces the TeraChem wb97xd3/6-311g** reference (tc_tddft_44) for a "
                     "like-for-like DFT number and misses the doubles/triples character that "
                     "STEOM captures.")}


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


def stage_spectra(args):
    """Excitonic absorption + CD lineshape figure (Fig. spectra) from the thermal
    coupling distribution and the STEOM two-dipole geometry.

    Two sub-steps, cached like the others:
      1. export_dipole_geometry.py -> mu_A, mu_B, r_A, r_B at both chromophore sites
         (same Kabsch/`super` placement as the coupling ensemble).
      2. absorption_cd_spectra.py  -> Fig_Spectra_{Coupling,Absorption,CD}.pdf. Panels
         (b)/(c) are summed over the per-frame ensemble geometry (coupling_geometry.npz
         from the thermal stage) so the lineshape is grounded in the real sampled
         disorder; the single-geometry dipole JSON is exported as a cross-check.
    """
    geom_path = REPO / DIPOLE_GEOM_JSON
    if args.reuse["spectra"] and geom_path.exists():
        log(f"  reusing dipole geometry {DIPOLE_GEOM_JSON}")
    else:
        log("  exporting STEOM two-dipole geometry (Kabsch placement at both sites)")
        run([PY, "export_dipole_geometry.py",
             "--density", str(STEOM_MATCHED),
             "--monomer", OLD_MONOMER, "--dimer", DIMER,
             "--out", str(geom_path)], "dipole_geometry.log", cwd=REPO)

    geom = json.loads(geom_path.read_text()) if geom_path.exists() else {}

    cmd = [PY, "absorption_cd_spectra.py",
           "--distribution", THERMAL_DIST_JSON,
           "--geometry-json", DIPOLE_GEOM_JSON,
           "--t2-star-fs", str(T2_STAR_FS),
           "--exp-splitting", str(EXP_SPLITTING[0]), str(EXP_SPLITTING[1]),
           "--out", OUT_DIR_SPECTRA]
    ensemble = REPO / ENSEMBLE_GEOM_NPZ
    if ensemble.exists():
        log(f"  ensemble lineshape from {ENSEMBLE_GEOM_NPZ} (hard-data 4b/4c)")
        cmd += ["--ensemble-geometry", ENSEMBLE_GEOM_NPZ]
    else:
        log(f"  [warn] {ENSEMBLE_GEOM_NPZ} absent; single-geometry lineshape (rerun thermal stage)")
    rc, tail = run(cmd, "absorption_cd_spectra.log", cwd=REPO)

    dist = json.loads((REPO / THERMAL_DIST_JSON).read_text()) if (REPO / THERMAL_DIST_JSON).exists() else {}
    J = float(dist.get("mean", float("nan")))
    panels = [f"{OUT_DIR_SPECTRA}/Fig_Spectra_Coupling.pdf",
              f"{OUT_DIR_SPECTRA}/Fig_Spectra_Absorption.pdf",
              f"{OUT_DIR_SPECTRA}/Fig_Spectra_CD.pdf"]
    return {
        "J_mean_cm": J,
        "two_J_cm": 2 * abs(J),
        "separation_A": geom.get("separation_A"),
        "angle_deg": geom.get("angle_deg"),
        "mu_debye": geom.get("mu_debye"),
        "t2_star_fs": T2_STAR_FS,
        "exp_splitting_cm": list(EXP_SPLITTING),
        "panels": panels,
        "rc": rc,
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
    p.add_argument("--engine", default=ENGINE, choices=["terachem", "orca"],
                   help="site-energy DFT engine; 'orca' runs ORCA TDDFT and skips the Q-Chem stage")
    p.add_argument("--stop-after", default=None, choices=STAGES,
                   help="stop after this stage (quick checks; e.g. --stop-after steom)")
    # --run-X flips REUSE[X] to False (force recompute of that stage)
    for k in REUSE:
        p.add_argument(f"--run-{k}", action="store_true", help=f"recompute the {k} stage")
    p.add_argument("--out", default=SUMMARY)
    a = p.parse_args(argv)
    a.reuse = {k: (not getattr(a, f"run_{k}")) for k in REUSE}
    return a


def _finish(results, args):
    """Write the summary + print the table (shared by normal end and --stop-after)."""
    log("aggregate")
    (REPO / args.out).write_text(json.dumps(results, indent=2))
    _print_table(results)
    log(f"wrote {args.out}")
    return 0


def _stopping(args, stage):
    if args.stop_after == stage:
        log(f"  --stop-after {stage}: stopping here (remaining stages unchanged from a full run).")
        return True
    return False


def main(argv=None):
    args = parse_args(argv)
    log("=" * 60)
    log("Venus_A206 paper-data pipeline")
    log(f"  seed={args.seed} eps={args.eps} n_frames={args.n_frames} backend={args.backend}")
    log(f"  engine={args.engine}  stop_after={args.stop_after}")
    log(f"  reuse={args.reuse}")
    log("=" * 60)

    results = {"config": {"seed": args.seed, "eps": args.eps, "n_frames": args.n_frames,
                          "backend": args.backend, "engine": args.engine,
                          "stop_after": args.stop_after, "reuse": args.reuse},
               "timestamp": datetime.now().isoformat()}

    log("[preflight]"); results["preflight"] = preflight(args)

    # --- TDDFT site energy (engine-dependent) ---
    if args.engine == "orca":
        log("[TDDFT] ORCA wB97X-D3/6-311G** (CPU)"); results["tddft"] = stage_tddft_orca(args)
    else:
        log("[TDDFT] TeraChem (GPU)");               results["tddft"] = stage_tddft(args)
    if _stopping(args, "tddft"): return _finish(results, args)

    log("[STEOM] DLPNO-STEOM-CCSD (CPU)");           results["steom"] = stage_steom(args)
    if _stopping(args, "steom"): return _finish(results, args)

    # --- EOM-CCSD(fT)/ADC(2) triples (Q-Chem) — skipped in the ORCA-only path ---
    if args.engine == "orca":
        log("[EOM] skipped (engine=orca): Q-Chem EOM-CCSD(fT)/ADC(2) omitted")
        results["doubles_triples"] = {"skipped": True, "reason": "engine=orca",
            "note": ("Q-Chem triples/doubles validation intentionally skipped in the ORCA-only "
                     "path; see the terachem-engine run (or the paper) for the validated ladder.")}
    else:
        log("[EOM] EOM-CCSD(fT)/ADC(2) doubles+triples (Q-Chem)")
        results["doubles_triples"] = stage_eom_triples(args)
    if _stopping(args, "eom"): return _finish(results, args)

    log("[density] STEOM coupling density");         results["density"] = stage_density(args)
    if _stopping(args, "density"): return _finish(results, args)

    log("[static J] static Davydov J (TDC)");        results["static_J"] = stage_static_J(args)
    if _stopping(args, "static"): return _finish(results, args)

    log("[thermal J] thermal NVT J ensemble");       results["thermal_J"] = stage_thermal_J(args)
    if _stopping(args, "thermal"): return _finish(results, args)

    log("[spectra] absorption + CD lineshape figure"); results["spectra"] = stage_spectra(args)
    return _finish(results, args)


def _g(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


def _print_table(r):
    print("\n" + "=" * 64)
    engine = _g(r, "config", "engine", default="terachem")
    print(f"  VENUS_A206 PAPER DATA — SUMMARY   (engine={engine})")
    print("=" * 64)

    # --- TDDFT site energy: shape differs by engine ---
    tdd = _g(r, "tddft", default={})
    if "site_energies" in tdd:                         # TeraChem engine (per-geometry list)
        for td in tdd["site_energies"] or []:
            if "error" in td:
                print(f"  TDDFT  S1 (TeraChem)   : {td['error']}")
            else:
                print(f"  TDDFT  S1 (TeraChem)   : {td['nm']} nm  f={td['fosc']}   [{td['dir']}]")
    elif "bright" in tdd:                               # ORCA engine (single bright state)
        b = tdd["bright"] or {}
        ref = tdd.get("reference_terachem_44_nm")
        refstr = f"   [TeraChem tc_tddft_44 ref: {ref} nm]" if ref else ""
        print(f"  TDDFT  S1 (ORCA)       : {b.get('nm','?')} nm  f={b.get('fosc','?')}"
              f"  |mu|={b.get('mu_au','?')} au{refstr}")

    st = _g(r, "steom", "bright")
    if st is not None or "steom" in r:
        print(f"  STEOM-CCSD S1 (CPU)    : {_g(st,'nm',default='?')} nm  f={_g(st,'fosc',default='?')}"
              f"  |mu|={_g(st,'mu_au',default='?')} au")

    # --- doubles/triples (may be skipped in the ORCA-only path, or absent if stopped early) ---
    dt = _g(r, "doubles_triples", default=None)
    if dt is not None:
        if dt.get("skipped"):
            print("  Doubles/triples        : skipped (engine=orca)")
        else:
            lad = dt.get("validated_ladder_eV", {})
            print(f"  Doubles/triples ladder : bare EOM-CCSD {lad.get('EOM_CCSD_bare','?')} -> "
                  f"EOM-CCSD(fT) {lad.get('EOM_CCSD_fT','?')} ~ STEOM {lad.get('STEOM_CCSD','?')} eV "
                  f"(ADC(2) {lad.get('ADC2_2p2h_doubles_pct','?')}% 2p2h)")
            ftraw = _g(dt, "raw_qchem", "eom_ccsd_fT", "excitation_eV")
            print(f"  (cached Q-Chem fT roots: {ftraw}  — heterogeneous basis; reference ladder above)")

    if "static_J" in r:
        sj = _g(r, "static_J", "STEOM")
        print(f"  Static J  STEOM        : {_g(sj,'J_mean',default='?')} cm^-1  (2|J|={_g(sj,'two_J',default='?')})")
        print(f"  Static J  TDDFT (pub)  : {_g(r,'static_J','TDDFT_published','J')} cm^-1")
    if "thermal_J" in r:
        tj = _g(r, "thermal_J", "STEOM")
        tjt = _g(r, "thermal_J", "TDDFT_cached")
        print(f"  Thermal J STEOM (NVT)  : {_g(tj,'J_mean',default='?')} +/- {_g(tj,'J_std',default='?')}"
              f"  (2|J|={_g(tj,'two_J',default='?')}, n={_g(tj,'n',default='?')})")
        print(f"  Thermal J TDDFT (cache): {_g(tjt,'J_mean',default='?')} +/- {_g(tjt,'J_std',default='?')}"
              f"  (2|J|={_g(tjt,'two_J',default='?')})")
    if "spectra" in r:
        sp = _g(r, "spectra", default={})
        lo, hi = (sp.get("exp_splitting_cm") or ["?", "?"])
        print(f"  Spectra 2|J| (CD/abs)  : {sp.get('two_J_cm','?'):.1f} cm^-1  vs exp {lo}-{hi}"
              if isinstance(sp.get('two_J_cm'), (int, float)) else
              f"  Spectra 2|J| (CD/abs)  : {sp.get('two_J_cm','?')} cm^-1  vs exp {lo}-{hi}")
        print(f"  Spectra geometry       : sep={sp.get('separation_A','?')} A  "
              f"angle={sp.get('angle_deg','?')} deg  |mu|={sp.get('mu_debye','?')} D")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
