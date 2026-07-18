#!/usr/bin/env python3
r"""
build_tandem_dimer.py  --  Build the actual VenusA206 tandem-dimer construct
(Kim 2019 / Nguyen 2025) for MM simulation.

The photophysically-relevant tandem is VenusA206-TD: two Venus barrels joined by
a ~33-residue flexible linker (with a TEV site), associating through the weak
A206 interface. The excitonic coupling is set by the A206 barrel-barrel docking,
which venus_dimer.pdb (built by build_dimer.py from the 1MYW crystal lattice)
already captures and the manuscript validated (25.5 A chromophore separation).
So the build KEEPS that validated interface and adds the covalent linker:

  1. Extract the two A206-docked barrels (chain A = FP1, chain B = FP2), dropping
     crystal waters/ions; keep the CR2 chromophore.
  2. Build the 33-aa linker  SGLRSENLYFQGPREFCRYPAQWRPLESRPRTT  (TEV ENLYFQG at
     6-12) with PyMOL `fab` as a helix (~50 A end-to-end, close to the 54 A
     FP1(C-term)->FP2(N-term) gap), and place it bridging the two termini.
  3. Close the C-terminal junction with an OpenMM minimisation of the LINKER
     ALONE (all standard residues -> amber ff14SB, no chromophore FF needed):
     freeze the N-terminal anchor, restrain the C-terminal C onto FP2's N.
  4. Assemble FP1 + closed linker + FP2 into ONE continuous chain, renumbered, so
     PDBFixer/OpenMM infer both peptide bonds. Output: tandem_dimer.pdb.

Then feed tandem_dimer.pdb to run_nvt.py (WITHOUT --restrain-interface): the
covalent linker is now the only tether, so the run tests whether the A206
interface (and the coupling) survives without the artificial restraint.

Env: run under the TeraChem conda env (pymol + openmm + pdbfixer).
"""
import argparse
from pathlib import Path

import numpy as np

LINKER_SEQ = "SGLRSENLYFQGPREFCRYPAQWRPLESRPRTT"   # 33 aa, TEV ENLYFQG at 6-12
PEPTIDE_CN = 1.33                                   # C(i)-N(i+1) bond length (A)


# --------------------------------------------------------------------------- #
# 1. barrels + anchors (PyMOL)
# --------------------------------------------------------------------------- #
def extract_barrels(dimer_pdb, outdir):
    from pymol import cmd
    cmd.reinitialize()
    cmd.load(str(dimer_pdb), "dim")
    # keep protein + the CR2 chromophore (a HETATM, so NOT `polymer`); drop only
    # crystal waters/ions.
    keep = "not (solvent or resn NA+CL+SO4+GOL+EDO+HOH+WAT+PO4+MG+ZN+K+ACT)"
    cmd.remove("solvent")
    fp1 = outdir / "fp1.pdb"
    fp2 = outdir / "fp2.pdb"
    cmd.save(str(fp1), f"dim and chain A and ({keep})")
    cmd.save(str(fp2), f"dim and chain B and ({keep})")

    def anchor(sel):
        xyz = cmd.get_coords(sel)
        return None if xyz is None else np.asarray(xyz[0], float)

    a_res = cmd.get_model("dim and chain A and polymer and name CA")
    b_res = cmd.get_model("dim and chain B and polymer and name CA")
    a_max = max(int(at.resi) for at in a_res.atom)
    b_min = min(int(at.resi) for at in b_res.atom)
    c_anchor = anchor(f"dim and chain A and resi {a_max} and name C")   # FP1 C-term carbonyl C
    n_anchor = anchor(f"dim and chain B and resi {b_min} and name N")   # FP2 N-term amide N
    com = np.asarray(cmd.centerofmass("dim and polymer"), float)        # dimer centre
    # backbone steric shell (CA/CB) the linker must route around
    shell = cmd.get_coords("dim and polymer and name CA+CB")
    return fp1, fp2, a_max, b_min, c_anchor, n_anchor, com, np.asarray(shell, float)


# --------------------------------------------------------------------------- #
# 2. build + place the linker (PyMOL)
# --------------------------------------------------------------------------- #
def build_place_linker(seq, c_anchor, n_anchor, com, outdir, bow=28.0):
    """Place the helical linker as a shallow bow offset OUTWARD from the dimer
    centre so it runs alongside the barrels (in solvent) rather than through them.
    Both ends start ~`bow` A outward from their junctions; the closure pulls the
    ends in while the barrel repulsion keeps the middle out."""
    from pymol import cmd
    cmd.delete("lnk")
    cmd.fab(seq, "lnk", ss=1)          # helix ~1.5 A/residue -> ~50 A for 33 aa
    cmd.remove("lnk and hydro")

    nN = np.asarray(cmd.get_coords("lnk and resi 1 and name N")[0], float)
    cC = np.asarray(cmd.get_coords(f"lnk and resi {len(seq)} and name C")[0], float)
    axis_l = cC - nN
    axis_l /= np.linalg.norm(axis_l)

    chord = n_anchor - c_anchor
    L = np.linalg.norm(chord)
    axis_t = chord / L
    mid = 0.5 * (c_anchor + n_anchor)
    outward = mid - com
    outward -= np.dot(outward, axis_t) * axis_t          # component perpendicular to chord
    outward /= (np.linalg.norm(outward) + 1e-9)

    R = _rot_between(axis_l, axis_t)                      # align helix axis to chord
    coords = cmd.get_coords("lnk")
    coords = (R @ (coords - 0.5 * (nN + cC)).T).T         # rotate about helix centre
    coords = coords + mid + bow * outward                 # centre on the outward-bowed midpoint
    cmd.load_coords(coords, "lnk")

    placed = outdir / "linker_placed.pdb"
    cmd.save(str(placed), "lnk")
    nN2 = np.asarray(cmd.get_coords("lnk and resi 1 and name N")[0], float)
    cC2 = np.asarray(cmd.get_coords(f"lnk and resi {len(seq)} and name C")[0], float)
    print(f"[place] chord {L:.1f} A, bow {bow:.0f} A outward; N-end->FP1(C) "
          f"{np.linalg.norm(nN2-c_anchor):.1f} A, C-end->FP2(N) {np.linalg.norm(cC2-n_anchor):.1f} A")
    return placed


def _rot_between(a, b):
    """Rotation matrix taking unit vector a onto unit vector b."""
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


# --------------------------------------------------------------------------- #
# 3. close the C-terminal junction (OpenMM, linker only, standard AAs)
# --------------------------------------------------------------------------- #
def close_linker(placed_pdb, c_anchor, n_anchor, shell, outdir, nsteps=40000):
    """Close BOTH junctions: pull the linker ends onto FP1(C)/FP2(N) while the
    barrel backbone (CA/CB) acts as frozen soft-repulsive centres so the linker
    routes around the barrels. Uses a raw Context (extra frozen particles beyond
    the linker topology)."""
    import openmm as mm
    import openmm.app as app
    from openmm import unit, Vec3
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(placed_pdb))
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    top = fixer.topology
    ff = app.ForceField("amber14-all.xml")
    system = ff.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None)
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    n_lnk = system.getNumParticles()

    pos = list(fixer.positions.value_in_unit(unit.nanometer))
    # add barrel CA/CB as frozen soft-repulsive particles, but EXCLUDE the
    # backbone near the two junctions (within 8 A of either anchor) so the linker
    # ends can reach the surface junction points; the mid-body stays repelled.
    kept = 0
    for xyz in shell:
        if (np.linalg.norm(xyz - c_anchor) < 8.0 or np.linalg.norm(xyz - n_anchor) < 8.0):
            continue
        system.addParticle(0.0)
        nb.addParticle(0.0, 0.40, 0.40)          # q=0, sigma=0.40 nm, eps=0.40 kJ/mol
        pos.append(Vec3(*(xyz / 10.0)))
        kept += 1
    print(f"[close] {kept}/{len(shell)} barrel shell atoms active (junction-proximal excluded)")

    res = list(top.residues())
    nterm_N = next(a.index for a in res[0].atoms() if a.name == "N")
    cterm_C = next(a.index for a in res[-1].atoms() if a.name == "C")
    nN0 = np.asarray(fixer.positions[nterm_N].value_in_unit(unit.angstrom))
    cC0 = np.asarray(fixer.positions[cterm_C].value_in_unit(unit.angstrom))
    n_tgt = c_anchor + PEPTIDE_CN * (nN0 - c_anchor) / np.linalg.norm(nN0 - c_anchor)
    c_tgt = n_anchor + PEPTIDE_CN * (cC0 - n_anchor) / np.linalg.norm(cC0 - n_anchor)

    rest = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    rest.addGlobalParameter("k", 500.0)
    for p in ("x0", "y0", "z0"):
        rest.addPerParticleParameter(p)
    rest.addParticle(nterm_N, [c / 10.0 for c in n_tgt])
    rest.addParticle(cterm_C, [c / 10.0 for c in c_tgt])
    system.addForce(rest)

    integ = mm.LangevinMiddleIntegrator(350 * unit.kelvin, 2.0 / unit.picosecond,
                                        0.001 * unit.picoseconds)
    ctx = mm.Context(system, integ, _platform())
    ctx.setPositions(pos)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 5000)     # relax fab clashes + soft pull
    for k in (1000.0, 4000.0, 12000.0):                   # ramp the end restraints
        ctx.setParameter("k", k)
        ctx.setVelocitiesToTemperature(350 * unit.kelvin)
        integ.step(nsteps // 3)
    ctx.setParameter("k", 30000.0)                        # snap the junctions shut
    mm.LocalEnergyMinimizer.minimize(ctx, 2.0, 8000)

    final = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    n_gap = float(np.linalg.norm(np.asarray(final[nterm_N]) - c_anchor))
    c_gap = float(np.linalg.norm(np.asarray(final[cterm_C]) - n_anchor))
    closed = outdir / "linker_closed.pdb"
    with open(closed, "w") as f:
        app.PDBFile.writeFile(top, (final[:n_lnk]) * unit.angstrom, f)
    print(f"[close] N-junction {n_gap:.2f} A, C-junction {c_gap:.2f} A "
          f"({'both BONDABLE' if max(n_gap, c_gap) < 2.0 else 'still open'})")
    return closed, max(n_gap, c_gap)


def _platform():
    # Prefer CPU for the tiny linker closure so it never contends with a GPU
    # coupling/NVT job running in parallel.
    import openmm as mm
    for name in ("CPU", "Reference", "OpenCL", "CUDA"):
        try:
            return mm.Platform.getPlatformByName(name)
        except Exception:
            continue
    return None


# --------------------------------------------------------------------------- #
# 4. assemble one continuous chain
# --------------------------------------------------------------------------- #
def _read_atoms(pdb):
    return [l for l in Path(pdb).read_text().splitlines()
            if l[:6] in ("ATOM  ", "HETATM")]


def _renumber(atom_lines, chain, start_res, start_serial):
    """Renumber residues consecutively from start_res, chain id, serials."""
    out = []
    serial = start_serial
    resmap = {}
    nxt = start_res
    for l in atom_lines:
        key = (l[21], l[22:27])              # orig chain + resSeq+icode
        if key not in resmap:
            resmap[key] = nxt; nxt += 1
        rn = resmap[key]
        out.append(f"{l[:6]}{serial:5d}{l[11:21]}{chain}{rn:4d} {l[27:54]}"
                   f"{l[54:] if len(l) > 54 else ''}")
        serial += 1
    return out, nxt, serial


def assemble(fp1, linker_closed, fp2, out_pdb):
    # linker_closed has PDBFixer hydrogens AND terminal caps (OXT) from when it
    # was a standalone peptide; strip both so it splices in as a clean INTERNAL
    # segment (run_nvt.py re-protonates). CAPS on the junction residues would
    # otherwise break the mid-chain FF template.
    caps = {"OXT", "OT1", "OT2", "HXT", "H1", "H2", "H3"}
    lnk = [l for l in _read_atoms(linker_closed)
           if l[76:78].strip() != "H" and l[12:16].strip()[:1] != "H"
           and l[12:16].strip() not in caps]
    a1 = _read_atoms(fp1)
    a2 = _read_atoms(fp2)

    all_out = []
    seg1, nxt, ser = _renumber(a1, "A", 1, 1)
    all_out += seg1
    seg2, nxt, ser = _renumber(lnk, "A", nxt, ser)
    all_out += seg2
    seg3, nxt, ser = _renumber(a2, "A", nxt, ser)
    all_out += seg3

    with open(out_pdb, "w") as f:
        f.write("REMARK   TANDEM VenusA206-TD: FP1 - 33aa flexible linker - FP2 "
                "(A206 interface preserved from venus_dimer.pdb)\n")
        f.writelines(l + "\n" for l in all_out)
        f.write("END\n")
    print(f"[assemble] wrote {out_pdb}: {len(all_out)} atoms, single chain A, "
          f"residues 1-{nxt-1}")
    return out_pdb


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dimer", default="venus_dimer.pdb")
    ap.add_argument("--seq", default=LINKER_SEQ)
    ap.add_argument("--out", default="tandem_dimer.pdb")
    ap.add_argument("--workdir", default="tandem_build")
    ap.add_argument("--close-steps", type=int, default=20000)
    args = ap.parse_args(argv)

    outdir = Path(args.workdir); outdir.mkdir(exist_ok=True)
    print(f"[*] Building VenusA206 tandem dimer from {args.dimer}")
    print(f"    linker ({len(args.seq)} aa): {args.seq}")

    fp1, fp2, a_max, b_min, c_anchor, n_anchor, com, shell = extract_barrels(args.dimer, outdir)
    print(f"[extract] FP1 C-term res {a_max}, FP2 N-term res {b_min}, "
          f"barrel shell {len(shell)} CA/CB atoms")

    placed = build_place_linker(args.seq, c_anchor, n_anchor, com, outdir)
    closed, gap = close_linker(placed, c_anchor, n_anchor, shell, outdir, nsteps=args.close_steps)
    if gap >= 2.5:
        print(f"[!] warning: junction still {gap:.2f} A open; peptide bond may not "
              f"be inferred. Consider more --close-steps.")
    assemble(fp1, closed, fp2, args.out)
    print("[done] feed to run_nvt.py WITHOUT --restrain-interface for the "
          "linker-tethered (unrestrained-interface) test.")


if __name__ == "__main__":
    main()
