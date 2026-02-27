# PyMOL Script: visualise_dimer_old_PCM_current.pml
# Visualize old-PCM-current TDDFT transition/difference densities on Venus dimer sites.

# 1. Environment
reinitialize
set bg_rgb, [1, 1, 1]
set grid_mode, 0
set transparency, 0.5
set cartoon_transparency, 0.6
set cartoon_fancy_helices, 1
set line_width, 1.0
set depth_cue, 0
set auto_zoom, 0

# 2. Load dimer arrangement reference
python
from pathlib import Path
from pymol import cmd

dimer_candidates = [
    Path("venus_dimer.pdb"),
    Path("tc_tddft_old_current_frame0002/venus_dimer.pdb"),
    Path("tc_tddft_old_current/venus_dimer.pdb"),
    Path("tc_tddft_old_PCM_current_frame0002/venus_dimer.pdb"),
]

dimer_path = None
for candidate in dimer_candidates:
    if candidate.exists():
        dimer_path = candidate
        break

if dimer_path is None:
    raise FileNotFoundError("Could not find venus_dimer.pdb in expected locations.")

cmd.load(str(dimer_path), "dimer_ref")
cmd.hide("everything", "dimer_ref")
print(f"Loaded dimer reference: {dimer_path}")
python end

# 3. Build alignment scaffold (TDDFT-first; no hard requirement on final_oldmin_relaxed.pdb)
python
from pathlib import Path
from pymol import cmd

monomer_candidates = [
    Path("tc_simple_old_PCM/final_oldmin_relaxed.pdb"),
    Path("tc_simple_old_PCM/classical_relaxed.pdb"),
    Path("tc_simple_old/final_oldmin_relaxed.pdb"),
    Path("tc_simple_old/classical_relaxed.pdb"),
]

monomer_path = None
for candidate in monomer_candidates:
    if candidate.exists():
        monomer_path = candidate
        break

if monomer_path is None:
    cmd.create("siteA", "dimer_ref and chain A")
    cmd.create("siteB", "dimer_ref and chain A")
    cmd.select("context_structure", "dimer_ref and (chain A or chain B)")
    cmd.select("alignment_scaffold_only", "siteA or siteB")
    print("No local monomer PDB found; using dimer_ref as alignment scaffold/context.")
else:
    cmd.load(str(monomer_path), "siteA")
    cmd.load(str(monomer_path), "siteB")
    if "final_oldmin_relaxed.pdb" in monomer_path.name:
        cmd.select("context_structure", "siteA or siteB")
        cmd.select("alignment_scaffold_only", "none")
        print(f"Loaded QM-minimized monomer source: {monomer_path}")
    else:
        cmd.select("context_structure", "dimer_ref and (chain A or chain B)")
        cmd.select("alignment_scaffold_only", "siteA or siteB")
        print(
            f"Using non-final monomer only as alignment scaffold: {monomer_path}. "
            "Context cartoon will use dimer_ref."
        )
python end

# 4. Align monomers to dimer sites
print("Aligning siteA to dimer_ref chain A...")
super siteA and chain A, dimer_ref and chain A
print("Aligning siteB to dimer_ref chain B...")
super siteB and chain A, dimer_ref and chain B

# 5. Visual appearance
as cartoon, context_structure
color gray85, context_structure
hide everything, alignment_scaffold_only

# Build full QM-region selection from PCM qm_region_atoms.txt when available.
python
from pathlib import Path
from pymol import cmd

qm_table_candidates = [
    Path("tc_simple_old_PCM/qm_region_atoms.txt"),
    Path("tc_simple_old/qm_region_atoms.txt"),
]

qm_table = None
for candidate in qm_table_candidates:
    if candidate.exists():
        qm_table = candidate
        break

qm_residues = []
seen = set()

if qm_table is not None:
    for raw_line in qm_table.read_text(errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("local_index"):
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        resn = parts[2].strip()
        resid = parts[3].strip()
        status = parts[6].strip()
        if status != "qm_atom":
            continue
        if resn == "LNK" or resid == "0":
            continue
        key = (resn, resid)
        if key in seen:
            continue
        seen.add(key)
        qm_residues.append(key)

if qm_residues:
    qm_clause = " or ".join(f"(resn {resn} and resi {resid})" for resn, resid in qm_residues)
    cmd.select("qm_region_siteA", f"siteA and ({qm_clause})")
    cmd.select("qm_region_siteB", f"siteB and ({qm_clause})")
    cmd.select("qm_residue_context", f"context_structure and ({qm_clause})")
    print(f"Loaded full QM region from {qm_table}: {len(qm_residues)} residues")
else:
    cmd.select("qm_region_siteA", "siteA and resn CR2")
    cmd.select("qm_region_siteB", "siteB and resn CR2")
    cmd.select("qm_residue_context", "context_structure and resn CR2")
    print("qm_region_atoms.txt not found/usable; fallback to CR2-only selection")

cmd.select("qm_region", "qm_region_siteA or qm_region_siteB")
python end

select non_qm_nonwater, context_structure and not qm_residue_context and not resn HOH+WAT+SOL
hide cartoon, non_qm_nonwater
show lines, non_qm_nonwater
spectrum count, rainbow, non_qm_nonwater

show sticks, qm_region
show spheres, qm_region
set sphere_scale, 0.25, qm_region
util.cbaw qm_region
color yellow, qm_region and elem H

# 6. Auto-detect and load density maps
python
import re
from pathlib import Path
from pymol import cmd


def first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def unique_dirs(paths):
    out = []
    seen = set()
    for p in paths:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        out.append(p)
    return out


def frame_index(path):
    m = re.search(r"frame(\d+)$", path.name)
    return int(m.group(1)) if m else -1


pcm_dirs = [p for p in sorted(Path(".").glob("tc_tddft_old_PCM_current_frame*")) if p.is_dir()]
legacy_dirs = [p for p in sorted(Path(".").glob("tc_tddft_old_current*")) if p.is_dir()]
analysis_dirs = [p for p in sorted(Path(".").glob("tc_tddft_analysis*")) if p.is_dir()]

candidate_dirs = []
candidate_dirs.extend(sorted(pcm_dirs, key=lambda x: (frame_index(x), x.stat().st_mtime), reverse=True))
candidate_dirs.extend(sorted(legacy_dirs, key=lambda x: x.stat().st_mtime, reverse=True))
candidate_dirs.extend(sorted(analysis_dirs, key=lambda x: x.stat().st_mtime, reverse=True))
candidate_dirs = unique_dirs(candidate_dirs)

latest_trans_file = None
selected_base = None
for wd in candidate_dirs:
    trans_candidates = list(wd.glob("abs_transdens_*.dx"))
    trans_candidates += list(wd.glob("transdens_*.dx"))
    trans_candidates += list((wd / "scr_plot").glob("abs_transdens_*.dx"))
    trans_candidates += list((wd / "scr_plot").glob("transdens_*.dx"))
    if not trans_candidates:
        continue
    latest_trans_file = max(trans_candidates, key=lambda p: p.stat().st_mtime)
    selected_base = wd
    break

if latest_trans_file is None:
    raise FileNotFoundError(
        "No transition-density files found in tc_tddft_old_PCM_current* (or fallback old/current dirs)."
    )

match = re.search(r"(?:abs_)?transdens_(\d+)\.dx$", latest_trans_file.name)
if not match:
    raise RuntimeError(f"Could not parse root index from: {latest_trans_file}")
selected_root = int(match.group(1))
selected_dir = latest_trans_file.parent
if selected_dir.name == "scr_plot":
    selected_base = selected_dir.parent
elif selected_base is None:
    selected_base = selected_dir

trans_path = first_existing([
    selected_base / f"abs_transdens_{selected_root}.dx",
    selected_base / f"transdens_{selected_root}.dx",
    selected_base / "scr_plot" / f"abs_transdens_{selected_root}.dx",
    selected_base / "scr_plot" / f"transdens_{selected_root}.dx",
])

if trans_path is None:
    raise FileNotFoundError(f"Could not resolve transition density for root {selected_root} in {selected_base}")

diff_path = first_existing([
    selected_base / f"abs_diffdens_{selected_root}.dx",
    selected_base / f"diffdens_{selected_root}.dx",
    selected_base / "scr_plot" / f"abs_diffdens_{selected_root}.dx",
    selected_base / "scr_plot" / f"diffdens_{selected_root}.dx",
])

print(f"Using TDDFT directory: {selected_base}")
print(f"Using root: {selected_root}")
print(f"Newest transition-density file: {latest_trans_file}")
print(f"Transition density: {trans_path}")
if diff_path is not None:
    print(f"Difference density: {diff_path}")
else:
    print("Difference density not found; only transition density will be displayed.")

cmd.load(str(trans_path), "map_trans")
if diff_path is not None:
    cmd.load(str(diff_path), "map_diff")

# Prefer the exact optimized QM geometry used to generate the maps.
geom_path = selected_base / "geometry.xyz"
if geom_path.exists():
    hyd_ids = []
    try:
        with geom_path.open() as geom_handle:
            atom_count = int(geom_handle.readline().strip())
            geom_handle.readline()
            for atom_id in range(1, atom_count + 1):
                line = geom_handle.readline()
                if not line:
                    break
                symbol = line.split()[0].upper()
                if symbol == "H":
                    hyd_ids.append(atom_id)
    except Exception:
        hyd_ids = []

    cmd.load(str(geom_path), "qm_opt_A")
    cmd.load(str(geom_path), "qm_opt_B")
    cmd.matrix_copy("siteA", "qm_opt_A")
    cmd.matrix_copy("siteB", "qm_opt_B")
    cmd.select("qm_carve_A", "qm_opt_A")
    cmd.select("qm_carve_B", "qm_opt_B")
    cmd.hide("sticks", "qm_region")
    cmd.hide("spheres", "qm_region")
    cmd.show("sticks", "qm_opt_A or qm_opt_B")
    cmd.show("spheres", "qm_opt_A or qm_opt_B")
    cmd.set("sphere_scale", 0.25, "qm_opt_A or qm_opt_B")
    util.cbaw("qm_opt_A or qm_opt_B")
    if hyd_ids:
        hyd_clause = "+".join(str(i) for i in hyd_ids)
        cmd.select("qm_hyd_A", f"qm_opt_A and id {hyd_clause}")
        cmd.select("qm_hyd_B", f"qm_opt_B and id {hyd_clause}")
        cmd.color("yellow", "qm_hyd_A or qm_hyd_B")
        print(f"Colored hydrogens yellow from geometry indices: {len(hyd_ids)}")
    else:
        cmd.color("yellow", "(qm_opt_A or qm_opt_B) and elem H")
        print("Hydrogen index parse failed; fallback to element-based hydrogen coloring.")
    print(f"Loaded optimized QM geometry for display/carving: {geom_path}")
else:
    cmd.select("qm_carve_A", "qm_region_siteA")
    cmd.select("qm_carve_B", "qm_region_siteB")
    print("No optimized geometry.xyz in selected TDDFT directory; using residue-based QM carve selections.")
python end

# 7. Create and align density isosurfaces
python
from pymol import cmd

ISO_VAL = 0.000005

if "map_trans" not in cmd.get_names("objects"):
    raise RuntimeError("Transition density map was not loaded; check TDDFT folder selection and DX files.")

# Transition density maps per monomer site
cmd.copy("map_trans_A", "map_trans")
cmd.copy("map_trans_B", "map_trans")
cmd.matrix_copy("siteA", "map_trans_A")
cmd.matrix_copy("siteB", "map_trans_B")

cmd.isosurface("trans_A_pos", "map_trans_A", ISO_VAL, "qm_carve_A", carve=6.0)
cmd.isosurface("trans_A_neg", "map_trans_A", -ISO_VAL, "qm_carve_A", carve=6.0)
cmd.isosurface("trans_B_pos", "map_trans_B", ISO_VAL, "qm_carve_B", carve=6.0)
cmd.isosurface("trans_B_neg", "map_trans_B", -ISO_VAL, "qm_carve_B", carve=6.0)

cmd.color("red", "trans_*_pos")
cmd.color("blue", "trans_*_neg")

# Difference density is optional
if "map_diff" in cmd.get_names("objects"):
    cmd.copy("map_diff_A", "map_diff")
    cmd.copy("map_diff_B", "map_diff")
    cmd.matrix_copy("siteA", "map_diff_A")
    cmd.matrix_copy("siteB", "map_diff_B")

    cmd.isosurface("diff_A_pos", "map_diff_A", ISO_VAL, "qm_carve_A", carve=6.0)
    cmd.isosurface("diff_A_neg", "map_diff_A", -ISO_VAL, "qm_carve_A", carve=6.0)
    cmd.isosurface("diff_B_pos", "map_diff_B", ISO_VAL, "qm_carve_B", carve=6.0)
    cmd.isosurface("diff_B_neg", "map_diff_B", -ISO_VAL, "qm_carve_B", carve=6.0)

    cmd.color("green", "diff_*_pos")
    cmd.color("magenta", "diff_*_neg")

# Keep source maps hidden after surfaces are built
for obj_name in [
    "map_trans",
    "map_trans_A",
    "map_trans_B",
    "map_diff",
    "map_diff_A",
    "map_diff_B",
]:
    if obj_name in cmd.get_names("objects"):
        cmd.disable(obj_name)

python end

# 8. Organization
group Transition_Density, trans_*
group Difference_Density, diff_*
disable Difference_Density
enable Transition_Density

# 9. Final view
orient context_structure
zoom context_structure
ray
png dimer_view_old_PCM_current.png, width=1600, height=1200, dpi=300
