#!/usr/bin/env python3
"""Register diagnostics for the AlphaFold Server runs in this directory.

What this measures, and why each number is here:

  barrel fold RMSD   CA RMSD of each predicted barrel against the corresponding
                     chain of the crystal dimer (venus_dimer.pdb).  AF3 is
                     excellent at this; it is the control that says the model is
                     not broken.

  register RMSD      Superpose barrel 1 on crystal chain A, carry that same
                     transform onto barrel 2, and measure barrel 2 against
                     crystal chain B.  This is the number the whole project
                     turns on: it says how far AF3 puts the second barrel from
                     where the crystal puts it, once the first is aligned.

  inter-domain PAE   Mean predicted aligned error between the two structured
                     cores, from *_full_data_*.json.  AF3's own statement of how
                     little it knows about the register.  nan when the run
                     directory has no full_data files.

  pLDDT              Mean per-residue confidence over each structured barrel and
                     over the linker between them.  High barrels next to a low
                     linker is the signature of a model that knows each domain
                     and not their arrangement.

  delta, alpha       Chromophore-chromophore angles, defined only when the
                     construct has two chromophores.  Reported as n/a for TDX,
                     which has one.

Run directories are the unpacked AlphaFold Server outputs:

    python3 alphafold_tandem/analyse_af.py
    python3 alphafold_tandem/analyse_af.py af3_tdx_server_run

The construct is read from each run's *_job_request.json rather than assumed,
including where its chromophore-forming triads are, so a 511-aa TD run and a
510-aa TDX run can be compared in one table without editing anything here.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
CRYSTAL = REPO / "venus_dimer.pdb"
DEFAULT_RUNS = ("af3_server_run", "af3_tdx_server_run")
BARREL_START = "MVSKGEELFTG"
CORE = 230              # crystal residues 0..230 are the structured barrel
IMID = ("N2", "CA2", "C2", "N3", "O2")


def angle(u, v):
    return float(np.degrees(np.arccos(np.clip(
        np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def read_crystal(path: Path):
    ca, cr2 = {}, {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain, name = line[21], line[12:16].strip()
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if name == "CA":
            ca.setdefault(chain, {})[int(line[22:26])] = xyz
        if line[17:20].strip() == "CR2":
            cr2.setdefault(chain, {})[name] = xyz
    return ca, cr2


def read_cif(path: Path):
    """auth_seq_id -> {atom: xyz}, plus residue names."""
    atoms: dict[int, dict[str, np.ndarray]] = {}
    names: dict[int, str] = {}
    plddt: dict[int, float] = {}
    inloop = False
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            inloop = True
            continue
        if not inloop:
            continue
        p = line.split()
        if len(p) < 18:
            if atoms:
                break
            continue
        if p[0] != "ATOM":
            continue
        number = int(p[15])
        names[number] = p[5]
        plddt[number] = float(p[14])          # AF3 writes pLDDT in B_iso_or_equiv
        atoms.setdefault(number, {})[p[3]] = np.array(
            [float(p[10]), float(p[11]), float(p[12])])
    return atoms, names, plddt


def construct(sequence: str):
    """Where the barrels start and which chromophore triads survived.

    Returns (barrel_starts, tyrosines, deletions) in 1-based numbering of THIS
    sequence.  `deletions` are positions where this construct is one residue
    shorter than the two-chromophore parent, which is what shifts the mapping
    onto the crystal downstream of a dark barrel's triad.
    """
    starts = [m.start() + 1 for m in re.finditer(BARREL_START, sequence)]
    tyrosines = [m.start() + 2 for m in re.finditer("GYG", sequence)]
    deletions = []
    for start in starts:
        # each barrel should carry a triad at barrel-local 66-68; a barrel whose
        # local 67 is not a tyrosine has lost it, either to a substitution
        # (same length) or to a deletion (one shorter).
        local = start + 65                       # 1-based position of triad G
        window = sequence[local - 1:local + 2]
        if window[:1] == "G" and window[1:2] == "G":
            deletions.append(local + 1)          # the missing residue's parent position
    return starts, tyrosines, deletions


def crystal_index(number: int, offset: int, deletions: list[int]) -> int:
    """AF residue number -> crystal residue number, allowing for deletions."""
    return number - offset + sum(1 for d in deletions if number >= d)


def analyse_run(directory: Path, ca, cr2):
    requests = sorted(directory.glob("*job_request.json"))
    if not requests:
        print(f"  {directory.name}: no *job_request.json, skipping", file=sys.stderr)
        return []
    request = json.load(open(requests[0]))
    request = request[0] if isinstance(request, list) else request
    sequence = request["sequences"][0]["proteinChain"]["sequence"]
    seed = request["modelSeeds"][0]
    starts, tyrosines, deletions = construct(sequence)
    if len(starts) != 2:
        print(f"  {directory.name}: found {len(starts)} barrels, skipping", file=sys.stderr)
        return []
    offsets = [starts[0], starts[1]]      # AF number of each barrel's Met -> crystal 0

    rows = []
    for k in range(5):
        models = sorted(directory.glob(f"*_model_{k}.cif"))
        if not models:
            continue
        atoms, _, plddt = read_cif(models[0])
        cas = {n: a["CA"] for n, a in atoms.items() if "CA" in a}

        fits = []
        for offset, chain in zip(offsets, ("A", "B")):
            pairs = [(n, crystal_index(n, offset, deletions)) for n in cas]
            pairs = [(n, c) for n, c in pairs if 0 <= c <= CORE and c in ca[chain]]
            P = np.array([cas[n] for n, _ in pairs])
            Q = np.array([ca[chain][c] for _, c in pairs])
            R, t = kabsch(Q, P)
            fits.append({
                "rmsd": float(np.sqrt(((Q @ R.T + t - P) ** 2).sum(1).mean())),
                "n": len(pairs), "P": P, "Q": Q, "R": R, "t": t,
            })

        # register: fit barrel 1 onto crystal A, carry the transform to barrel 2
        R1, t1 = kabsch(fits[0]["Q"], fits[0]["P"])
        register = float(np.sqrt(((fits[1]["Q"] @ R1.T + t1 - fits[1]["P"]) ** 2).sum(1).mean()))

        # chromophore angles, only where both barrels have a chromophore
        delta = alpha = float("nan")
        if len(tyrosines) == 2:
            axes = []
            for tyr, fit, chain in zip(tyrosines, fits, ("A", "B")):
                moved = {a: (x @ fit["R"].T + fit["t"]) for a, x in cr2[chain].items()}
                axes.append(np.mean([moved[a] for a in IMID], axis=0) - moved["OH"])
            alpha = angle(axes[0], axes[1])
            u = atoms[tyrosines[0]]["CB"] - atoms[tyrosines[0]]["OH"]
            v = atoms[tyrosines[1]]["CB"] - atoms[tyrosines[1]]["OH"]
            delta = 90 - angle(u, v) / 2

        barrel_plddt = [float(np.mean([plddt[n] for n in range(o, o + CORE + 1)
                                       if n in plddt])) for o in offsets]
        linker = [n for n in plddt if offsets[0] + CORE < n < offsets[1]]
        linker_plddt = float(np.mean([plddt[n] for n in linker])) if linker else float("nan")

        # Inter-domain PAE.  The *_full_data_*.json files are MB-scale per model
        # and .gitignore excludes them, so the derived number is cached in a
        # small tracked file the first time it is computed; after that the run
        # directory reports its PAE with or without the bulk output present.
        pae = float("nan")
        cache = directory / "pae_inter_barrel.json"
        full = sorted(directory.glob(f"*full_data_{k}.json"))
        if full:
            matrix = np.array(json.load(open(full[0]))["pae"], float)
            A = slice(offsets[0] - 1, offsets[0] + CORE)
            B = slice(offsets[1] - 1, offsets[1] + CORE)
            pae = float(np.concatenate([matrix[A, B].ravel(), matrix[B, A].ravel()]).mean())
            stored = json.loads(cache.read_text()) if cache.exists() else {}
            if stored.get(str(k)) != round(pae, 3):
                stored[str(k)] = round(pae, 3)
                stored["_note"] = ("mean PAE between the two structured barrel cores, "
                                   "from *_full_data_*.json, which .gitignore excludes")
                cache.write_text(json.dumps(stored, indent=1, sort_keys=True) + "\n")
        elif cache.exists():
            pae = json.loads(cache.read_text()).get(str(k), float("nan"))

        rows.append(dict(run=directory.name, n=len(sequence), seed=seed, model=k,
                         chromophores=len(tyrosines), fit=np.mean([f["rmsd"] for f in fits]),
                         matched=fits[0]["n"], register=register, pae=pae,
                         delta=delta, alpha=alpha,
                         plddt1=barrel_plddt[0], plddt2=barrel_plddt[1],
                         plddt_linker=linker_plddt))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*", default=list(DEFAULT_RUNS),
                    help="AlphaFold Server run directories (default: both repo runs)")
    ap.add_argument("--crystal", type=Path, default=CRYSTAL)
    args = ap.parse_args()

    ca, cr2 = read_crystal(args.crystal)
    rows = []
    for name in args.runs:
        directory = Path(name)
        if not directory.is_dir():
            directory = HERE / name
        if not directory.is_dir():
            print(f"  {name}: not a directory, skipping", file=sys.stderr)
            continue
        rows.extend(analyse_run(directory, ca, cr2))

    if not rows:
        raise SystemExit("no runs analysed")

    print(f"crystal reference: {args.crystal}")
    for run in dict.fromkeys(r["run"] for r in rows):
        r = [x for x in rows if x["run"] == run]
        F = np.array([x["fit"] for x in r])
        G = np.array([x["register"] for x in r])
        P = np.array([x["pae"] for x in r])
        D = np.array([x["delta"] for x in r])
        A = np.array([x["alpha"] for x in r])
        seeds = sorted({x["seed"] for x in r})
        print(f"\n=== {run}: {r[0]['n']}-aa construct, {r[0]['chromophores']} chromophore(s), "
              f"{len(seeds)} seed(s) x {len(r)} model(s) ===")
        print(f"   CA pairs matched to crystal  : {r[0]['matched']} per barrel")
        print(f"   barrel fold RMSD vs crystal  : {F.mean():.2f} A "
              f"(range {F.min():.2f}-{F.max():.2f})")
        print(f"   inter-domain PAE (cores)     : "
              + ("n/a (no full_data files)" if np.isnan(P).all()
                 else f"{np.nanmean(P):.1f} A (range {np.nanmin(P):.1f}-{np.nanmax(P):.1f})"))
        B1 = np.array([x["plddt1"] for x in r])
        B2 = np.array([x["plddt2"] for x in r])
        BL = np.array([x["plddt_linker"] for x in r])
        print(f"   pLDDT barrel 1 / barrel 2    : {B1.mean():.1f} / {B2.mean():.1f}")
        print(f"   pLDDT linker                 : {BL.mean():.1f}")
        print(f"   register RMSD vs crystal     : {G.mean():.1f} A "
              f"(range {G.min():.1f}-{G.max():.1f})")
        if np.isnan(D).all():
            print("   delta, alpha                 : n/a (one chromophore)")
        else:
            print(f"   delta (Cusick)               : {D.mean():.1f} +- {D.std():.1f}")
            print(f"   grafted alpha                : {A.mean():.1f} +- {A.std():.1f}")


if __name__ == "__main__":
    main()
