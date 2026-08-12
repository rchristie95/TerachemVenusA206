"""Audit the magnetic transition dipoles / rotational strengths TeraChem prints.

TeraChem *does* print, for every root,

    Magnetic transition dipole moments and rotational strengths:
      Root   Lx   Ly   Lz   R(len)   R(vel) (a.u.)

The question this script answers is whether those numbers are *usable* for an
absolute CD prediction. Three things have to hold:

  1. Conventions. What are L, R(len), R(vel) actually? We recover them by
     numerical identity against the length and velocity transition dipoles.
  2. Gauge agreement. R(len) and R(vel) are the same observable computed two
     ways. They agree only to the extent that the velocity dipole satisfies
     mu_vel = -omega * mu_len (the exact off-diagonal hypervirial relation).
  3. Origin independence. R(len) is origin independent *only* through that same
     hypervirial relation: shifting the gauge origin by a changes

         L -> L - a x mu_vel   =>   Delta R(len) = a . (mu_vel x mu_len)

     which vanishes iff mu_vel || mu_len. So (2) and (3) are the same defect
     measured two ways, and (3) is the one with teeth: it sets the length scale
     over which the answer changes by 100%.

Run from the repository root:
    python terachem_site_energy_cd/scripts/audit_magnetic_dipoles.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
BOHR_TO_ANGSTROM = 0.529177210903

_ENERGY_HEAD = re.compile(r"^\s*Root\s+Total Energy \(a\.u\.\)")
_ROW = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+(?:\s+-?\d+\.\d+)*)\s*$")


def _table(lines: list[str], start: int, ncol: int) -> dict[int, np.ndarray]:
    """Read a `Root  v1 v2 ...` table beginning at/after line `start`."""
    out: dict[int, np.ndarray] = {}
    for line in lines[start:]:
        m = _ROW.match(line)
        if m:
            vals = np.array([float(x) for x in m.group(2).split()], float)
            if vals.size < ncol:
                continue
            out[int(m.group(1))] = vals[:ncol]
        elif out:
            break
    return out


def _find(lines: list[str], exact: str) -> int | None:
    for i, line in enumerate(lines):
        if line.strip() == exact:
            return i
    return None


def parse(path: Path) -> dict | None:
    lines = path.read_text(errors="replace").splitlines()

    i = next((k for k, l in enumerate(lines) if _ENERGY_HEAD.match(l)), None)
    if i is None:
        return None
    energies = _table(lines, i + 1, 4)  # Etot, Eex(au), Eex(eV), Eex(nm)

    idx_len = _find(lines, "Transition dipole moments:")
    idx_vel = _find(lines, "Velocity transition dipole moments:")
    idx_mag = _find(
        lines, "Magnetic transition dipole moments and rotational strengths:"
    )
    if None in (idx_len, idx_vel, idx_mag):
        return None

    mu_len = _table(lines, idx_len + 1, 4)  # Tx Ty Tz |T|
    mu_vel = _table(lines, idx_vel + 1, 4)
    mag = _table(lines, idx_mag + 1, 5)  # Lx Ly Lz R(len) R(vel)

    roots = sorted(set(energies) & set(mu_len) & set(mu_vel) & set(mag))
    if not roots:
        return None
    return {
        "path": path,
        "roots": roots,
        "omega": {r: energies[r][1] for r in roots},
        "nm": {r: energies[r][3] for r in roots},
        "mu_len": {r: mu_len[r][:3] for r in roots},
        "mu_vel": {r: mu_vel[r][:3] for r in roots},
        "L": {r: mag[r][:3] for r in roots},
        "r_len": {r: mag[r][3] for r in roots},
        "r_vel": {r: mag[r][4] for r in roots},
    }


def main() -> None:
    files = sorted(RESULTS.glob("v2_camb3lyp_frame_*/site_*/tddft.out"))
    runs = [d for d in (parse(p) for p in files) if d]
    print(f"parsed {len(runs)} / {len(files)} runs\n")

    # ---- 1. conventions, tested on every root of every run -----------------
    res_len, res_vel, rows = [], [], []
    for d in runs:
        for r in d["roots"]:
            mu, mv, L = d["mu_len"][r], d["mu_vel"][r], d["L"][r]
            w = d["omega"][r]
            res_len.append(d["r_len"][r] - (-mu @ L))
            res_vel.append(d["r_vel"][r] - (mv @ L) / w)
            rows.append((d, r))
    res_len = np.array(res_len)
    res_vel = np.array(res_vel)
    print("CONVENTION CHECK  (all roots, all runs; n = %d)" % len(res_len))
    print(f"  R(len) == -mu_len . L        max |residual| = {np.abs(res_len).max():.2e}")
    print(f"  R(vel) == (mu_vel . L)/omega max |residual| = {np.abs(res_vel).max():.2e}")
    print("  -> L is the raw angular-momentum matrix element (no 1/2c),")
    print("     R is quoted in the same convention. Both gauges use the SAME L.\n")

    # ---- 2. bright state, gauge agreement, origin sensitivity ---------------
    rec = []
    for d in runs:
        r = max(d["roots"], key=lambda k: np.linalg.norm(d["mu_len"][k]))
        mu, mv, L = d["mu_len"][r], d["mu_vel"][r], d["L"][r]
        w = d["omega"][r]
        nmu, nmv = np.linalg.norm(mu), np.linalg.norm(mv)
        # exact relation: mu_vel = -omega mu_len  (TeraChem's own sign convention,
        # fixed by R(len) = -mu.L vs R(vel) = +mu_vel.L/omega)
        cos = float(mv @ mu) / (nmv * nmu)
        f_len = (2.0 / 3.0) * w * nmu**2
        f_vel = (2.0 / 3.0) * nmv**2 / w
        # origin sensitivity of R(len): dR/da = mu_vel x mu_len, per bohr
        sens = np.linalg.norm(np.cross(mv, mu))
        rec.append(
            dict(
                path=str(d["path"].relative_to(RESULTS)),
                root=r,
                nm=d["nm"][r],
                r_len=d["r_len"][r],
                r_vel=d["r_vel"][r],
                f_len=f_len,
                f_vel=f_vel,
                ratio_f=f_len / f_vel,
                cos=cos,
                angle_deg=float(np.degrees(np.arccos(np.clip(cos, -1, 1)))),
                # ratio of velocity dipole to its exact value -omega*mu_len
                vel_defect=nmv / (w * nmu),
                sens_per_bohr=sens,
                # origin shift (Angstrom) that changes R(len) by 100 % of itself
                d100_ang=abs(d["r_len"][r]) / sens * BOHR_TO_ANGSTROM,
                L_norm=float(np.linalg.norm(L)),
            )
        )

    def col(k):
        return np.array([x[k] for x in rec], float)

    print("BRIGHT STATE  (largest |mu_len| per run; n = %d)" % len(rec))
    print(f"  lambda            = {col('nm').mean():7.2f} +- {col('nm').std():.2f} nm")
    print(f"  f (length gauge)  = {col('f_len').mean():7.4f} +- {col('f_len').std():.4f}")
    print(f"  f (velocity gauge)= {col('f_vel').mean():7.4f} +- {col('f_vel').std():.4f}")
    print(f"  f_len / f_vel     = {col('ratio_f').mean():7.2f}  "
          f"[{col('ratio_f').min():.1f}, {col('ratio_f').max():.1f}]")
    print(f"  |mu_vel|/(w|mu_len|) = {col('vel_defect').mean():.4f} "
          f"(exact value: 1.0)")
    print(f"  angle(mu_vel, mu_len)= {col('angle_deg').mean():7.1f} deg "
          f"+- {col('angle_deg').std():.1f}   (exact value: 180 deg)")
    print()

    print("SIGN CONSISTENCY OF THE ROTATIONAL STRENGTH")
    same = np.sign(col("r_len")) == np.sign(col("r_vel"))
    print(f"  R(len) and R(vel) agree in sign in {same.sum()} / {len(rec)} runs "
          f"({100*same.mean():.1f} %)")
    print(f"  R(len) = {col('r_len').mean():+.4f} +- {col('r_len').std():.4f}")
    print(f"  R(vel) = {col('r_vel').mean():+.4f} +- {col('r_vel').std():.4f}")
    print()

    print("ORIGIN DEPENDENCE OF R(len)")
    print(f"  |L| about TeraChem's origin = {col('L_norm').mean():.3f} a.u. "
          f"(would be ~20 a.u. if the origin were the coordinate origin at "
          f"~112 A, so TeraChem centres it on the molecule)")
    print(f"  |dR(len)/da| = {col('sens_per_bohr').mean():.4f} +- "
          f"{col('sens_per_bohr').std():.4f} per bohr")
    d100 = col("d100_ang")
    print(f"  origin shift that changes R(len) by 100 %:")
    print(f"      mean {d100.mean():.3f} A, median {np.median(d100):.3f} A, "
          f"max {d100.max():.3f} A")
    print(f"  fraction of runs where < 1 A of origin ambiguity flips R(len) by "
          f"100 %: {100*(d100 < 1.0).mean():.1f} %")
    print()

    out = RESULTS / "magnetic_dipole_audit.json"
    out.write_text(json.dumps(rec, indent=1))
    print(f"wrote {out.relative_to(ROOT.parent)}")


if __name__ == "__main__":
    main()
