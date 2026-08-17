#!/usr/bin/env python3
"""gamma_0: the monomer angle between mu and delta_mu, from ORCA output.

Cusick et al. 2026 measure gamma_0 = 22 deg for the Venus chromophore in
dVenus-TDX -- the angle between the transition dipole mu and the difference
dipole delta_mu = mu_excited - mu_ground. This needs no dimer, no MD and no
coupling: one excited state of one monomer settles it. That makes it the
cleanest available validation of the electronic structure underpinning the whole
transition density, and therefore of the coupling itself.

No parser for ORCA state dipoles existed in this repo -- reproduce_paper.py:191
and weekend_dipole/collect_weekend.py:33 both read only the absorption block.
This is that parser.

TWO TRAPS, both live in these files:

1. There are FOUR "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
   blocks in a STEOM output. The first is the CIS/TDDFT guess (wrong energy
   entirely -- 331.8 nm vs 523.9 nm in steom_phenol_svpd.out). The other three
   are the RIGHT, LEFT and LEFT-RIGHT transition-moment conventions of the
   non-Hermitian coupled-cluster response. Taking the first block silently gives
   a garbage answer. We key every block off the "SPECTRUM FOR <X> TRANSITION
   MOMENTS" header that precedes it.

2. ORCA's "FINAL STEOM-CCSD ABSORPTION SPECTRUM" block is internally mixed: its
   fosc/D2 come from the LEFT-RIGHT geometric mean, but its printed DX/DY/DZ are
   the RIGHT vector. In steom_phenol_svpd.out the FINAL row shows D2 = 15.26429
   (so |mu| = 3.9070 au) alongside components whose norm is 3.8810 au. That 0.7%
   gap is not rounding -- it is the L/R asymmetry of EOM-CC transition moments.
   For gamma_0 only the DIRECTION matters, but the LEFT and RIGHT directions are
   not identical either, so we report gamma_0 under each convention and use the
   spread as an honest error bar.

CAVEAT ON delta_mu. These are UNRELAXED state dipoles -- ORCA's
"UNRELAXED EXCITED STATE DIPOLE MOMENTS" block, computed without orbital
relaxation. The relaxed/unrelaxed distinction matters a great deal for
delta_mu specifically (it is a difference of two large, nearly cancelling
vectors). A disagreement with Cusick's 22 deg is therefore a real, reportable
finding about this level of theory, not automatically an error.

Read-only. Writes results/gamma0.json.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = Path(__file__).resolve().parent / "results"

TARGETS = [
    ("steom_phenol_svpd", ROOT / "neo_model/orca_steom/steom_phenol_svpd.out",
     "production Venus chromophore, DLPNO-STEOM/def2-SVPD"),
    ("weekend_def2_svp", ROOT / "weekend_dipole/steom_def2_svp.out",
     "basis check, DLPNO-STEOM/def2-SVP, 44-atom anion in MM field"),
    ("weekend_def2_svpd", ROOT / "weekend_dipole/steom_def2_svpd.out",
     "basis check, DLPNO-STEOM/def2-SVPD, 44-atom anion in MM field"),
]

CUSICK_GAMMA0_DEG = 22.0

STATE_DIPOLE_HEADER = "UNRELAXED EXCITED STATE DIPOLE MOMENTS"
ABSORPTION_HEADER = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
CONVENTION_RE = re.compile(r"SPECTRUM FOR (RIGHT|LEFT|LEFT-RIGHT) TRANSITION MOMENTS")
FINAL_RE = re.compile(r"FINAL STEOM-CCSD ABSORPTION SPECTRUM")
IROOT_RE = re.compile(
    r"^IROOT=\s*(\d+):\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
)
TRANS_RE = re.compile(
    r"^\s*0-1A\s+->\s+(\d+)-1A\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
    r"([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
)


def parse_state_dipoles(lines):
    """IROOT -> {E_eV, D_au (3,)} from the unrelaxed state-dipole block."""
    out = {}
    for i, line in enumerate(lines):
        if STATE_DIPOLE_HEADER not in line:
            continue
        for probe in lines[i : i + 40]:
            m = IROOT_RE.match(probe.strip())
            if m:
                out[int(m.group(1))] = {
                    "E_eV": float(m.group(2)),
                    "D_au": [float(m.group(3)), float(m.group(4)), float(m.group(5))],
                    "D_debye": float(m.group(6)),
                }
            elif out and probe.strip().startswith("---"):
                break
        if out:
            break
    return out


def parse_absorption_blocks(lines):
    """Every absorption block, tagged by the convention header preceding it.

    Returns {label: {root: {...}}}. Labels are RIGHT / LEFT / LEFT-RIGHT /
    FINAL / CIS_GUESS. Anything before the first "SPECTRUM FOR" header is the
    CIS/TDDFT guess and must not be mistaken for a STEOM result.
    """
    blocks = {}
    pending = "CIS_GUESS"
    for i, line in enumerate(lines):
        m = CONVENTION_RE.search(line)
        if m:
            pending = m.group(1)
            continue
        if FINAL_RE.search(line):
            pending = "FINAL"
            continue
        if ABSORPTION_HEADER not in line:
            continue
        roots = {}
        for probe in lines[i + 1 : i + 40]:
            t = TRANS_RE.match(probe)
            if t:
                roots[int(t.group(1))] = {
                    "E_eV": float(t.group(2)),
                    "E_cm": float(t.group(3)),
                    "wavelength_nm": float(t.group(4)),
                    "fosc": float(t.group(5)),
                    "D2_au2": float(t.group(6)),
                    "mu_au": [float(t.group(7)), float(t.group(8)), float(t.group(9))],
                }
            elif roots and (probe.strip().startswith("---") or not probe.strip()):
                break
        if roots:
            blocks[pending] = roots
        pending = "CIS_GUESS"
    return blocks


def angle_deg(u, v):
    u = np.asarray(u, float)
    v = np.asarray(v, float)
    c = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def fold_0_90(a):
    """Transition-dipole sign is an arbitrary phase, so gamma and 180-gamma are
    the same physical angle. Fold to [0, 90] and report the raw value too."""
    return a if a <= 90.0 else 180.0 - a


def analyse(path):
    lines = path.read_text(errors="replace").splitlines()
    states = parse_state_dipoles(lines)
    blocks = parse_absorption_blocks(lines)
    if not states or not blocks:
        return {"error": "missing state-dipole or absorption block"}

    steom = {k: v for k, v in blocks.items() if k != "CIS_GUESS"}
    if not steom:
        return {"error": "no STEOM absorption block (only the CIS guess)"}

    # Bright state = max fosc in the FINAL block if present, else RIGHT.
    ref_label = "FINAL" if "FINAL" in steom else sorted(steom)[0]
    ref = steom[ref_label]
    bright = max(ref, key=lambda r: ref[r]["fosc"])

    if 0 not in states or bright not in states:
        return {"error": f"state dipoles missing IROOT 0 or {bright}"}

    # Guard the root<->IROOT mapping by energy rather than trusting the index.
    e_abs = ref[bright]["E_eV"]
    e_state = states[bright]["E_eV"]
    if abs(e_abs - e_state) > 0.01:
        return {
            "error": f"root/IROOT energy mismatch: absorption {e_abs} eV vs "
                     f"state dipole {e_state} eV"
        }

    d_gs = np.array(states[0]["D_au"])
    d_ex = np.array(states[bright]["D_au"])
    delta_mu = d_ex - d_gs

    per_convention = {}
    for label, roots in sorted(steom.items()):
        if bright not in roots:
            continue
        mu = np.array(roots[bright]["mu_au"])
        raw = angle_deg(mu, delta_mu)
        d2 = roots[bright]["D2_au2"]
        per_convention[label] = {
            "mu_au": mu.tolist(),
            "mu_norm_from_components_au": float(np.linalg.norm(mu)),
            "mu_norm_from_D2_au": float(np.sqrt(d2)) if d2 > 0 else None,
            "wavelength_nm": roots[bright]["wavelength_nm"],
            "gamma0_raw_deg": raw,
            "gamma0_deg": fold_0_90(raw),
        }

    folded = [v["gamma0_deg"] for v in per_convention.values()]
    return {
        "bright_root": bright,
        "bright_label": ref_label,
        "wavelength_nm": ref[bright]["wavelength_nm"],
        "ground_state_dipole_au": d_gs.tolist(),
        "excited_state_dipole_au": d_ex.tolist(),
        "delta_mu_au": delta_mu.tolist(),
        "delta_mu_norm_au": float(np.linalg.norm(delta_mu)),
        "state_dipoles_are_unrelaxed": True,
        "per_convention": per_convention,
        "gamma0_deg_mean": float(np.mean(folded)),
        "gamma0_deg_spread": float(np.max(folded) - np.min(folded)),
        "cusick_gamma0_deg": CUSICK_GAMMA0_DEG,
        "discrepancy_deg": float(np.mean(folded) - CUSICK_GAMMA0_DEG),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    for key, path, description in TARGETS:
        if not path.exists():
            results[key] = {"error": f"missing file {path}"}
            continue
        r = analyse(path)
        r["source"] = str(path.relative_to(ROOT))
        r["description"] = description
        results[key] = r

    payload = {
        "caveat": (
            "delta_mu is built from ORCA UNRELAXED excited-state dipoles. Orbital "
            "relaxation matters disproportionately for delta_mu because it is a "
            "difference of two large, nearly cancelling vectors. Treat a "
            "disagreement with Cusick's 22 deg as a finding about this level of "
            "theory, not as a settled refutation."
        ),
        "convention_note": (
            "ORCA's FINAL STEOM block mixes conventions: fosc/D2 are the "
            "LEFT-RIGHT geometric mean while DX/DY/DZ are the RIGHT vector. The "
            "|mu| discrepancy that produces is real, not rounding. gamma_0 is "
            "reported per convention; the spread is the honest error bar."
        ),
        "results": results,
    }
    (OUT_DIR / "gamma0.json").write_text(json.dumps(payload, indent=2))

    for key, r in results.items():
        print(f"=== {key} ===")
        if "error" in r:
            print(f"  ERROR: {r['error']}")
            continue
        print(f"  {r['description']}")
        print(f"  bright root {r['bright_root']} at {r['wavelength_nm']} nm")
        print(f"  |delta_mu| = {r['delta_mu_norm_au']:.4f} au (unrelaxed)")
        for label, v in r["per_convention"].items():
            nd2 = v["mu_norm_from_D2_au"]
            nd2s = f"{nd2:.4f}" if nd2 is not None else "  n/a "
            print(
                f"    {label:10s} gamma0 = {v['gamma0_deg']:6.2f} deg "
                f"(raw {v['gamma0_raw_deg']:6.2f})   "
                f"|mu| comp {v['mu_norm_from_components_au']:.4f} vs sqrt(D2) {nd2s}"
            )
        print(
            f"  gamma0 = {r['gamma0_deg_mean']:.2f} deg "
            f"(spread {r['gamma0_deg_spread']:.2f})   "
            f"Cusick {CUSICK_GAMMA0_DEG} deg   "
            f"discrepancy {r['discrepancy_deg']:+.2f} deg"
        )
        print()
    print(f"wrote {OUT_DIR/'gamma0.json'}")


if __name__ == "__main__":
    main()
