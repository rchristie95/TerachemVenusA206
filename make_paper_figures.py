# make_paper_figures.py  — STEOM-CCSD versions of the paper's four density panels.
# Format matched to visualise_dimer.pml / the JPCL figures:
#   transition density  rho_0n         -> red (+) / blue (-)
#   difference density  rho_e - rho_g  -> green (gain) / violet (loss)
#   each shown as (a) full dimer view (rainbow MM backbone) and (b) QM close-up (site B).
from pymol import cmd, util
import os

DIMER = "/home/robson/PetaChem/tandem_dimer.pdb"
MONO  = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"
TRANS = "/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.cube"
DIFF  = "/home/robson/PetaChem/neo_model/orca_steom/steom_diffdens.cube"
OUT   = "/home/robson/PetaChem/plots_steom"; os.makedirs(OUT, exist_ok=True)
ISO, W, H, QM = 0.04, 1600, 1200, "(resi 66 or resi 202)"     # CR2 + Tyr203 phenol (202 in PDB)

cmd.reinitialize()
cmd.bg_color("white")
for k, v in [("ray_opaque_background",1),("ray_shadows",0),("antialias",2),("ambient",0.45),
             ("transparency",0.35),("cartoon_transparency",0.65),("cartoon_fancy_helices",1),
             ("depth_cue",0),("valence",0),("sphere_scale",0.22)]:
    cmd.set(k, v)

# dimer scaffold + a monomer aligned onto each site
cmd.load(DIMER, "dimer_ref")
cmd.load(MONO, "siteA"); cmd.load(MONO, "siteB")
cmd.align("siteA and name CA", "dimer_ref and resi 1-229 and name CA")
cmd.align("siteB and name CA", "dimer_ref and resi 263-491 and name CA")
cmd.matrix_copy("siteA", "siteA")  # ensure matrix is applied
cmd.matrix_copy("siteB", "siteB")

cmd.hide("everything")
cmd.show("sticks", "dimer_ref and backbone")
cmd.set("stick_radius", 0.05, "dimer_ref")
cmd.set("stick_transparency", 0.45, "dimer_ref")
cmd.spectrum("count", "rainbow", "dimer_ref and resi 1-240")
cmd.spectrum("count", "rainbow", "dimer_ref and resi 250-500")

for s in ("siteA","siteB"):                                       # QM chromophore sticks+spheres
    cmd.show("sticks",  f"{s} and {QM}")
    cmd.show("spheres", f"{s} and {QM}")
util.cbaw(f"(siteA or siteB) and {QM}")                          # white carbons
cmd.color("yellow", f"(siteA or siteB) and {QM} and elem H")

# map the monomer densities onto both sites via the alignment matrices
cmd.load(TRANS, "mt"); cmd.load(DIFF, "md")
for src, dst, site in [("mt","mt_A","siteA"),("mt","mt_B","siteB"),
                       ("md","md_A","siteA"),("md","md_B","siteB")]:
    cmd.copy(dst, src); cmd.matrix_copy(site, dst)
cmd.disable("mt"); cmd.disable("md")

for s in ("A","B"):                                              # transition density: red/blue
    cmd.isosurface(f"tr_{s}_p", f"mt_{s}",  ISO)
    cmd.isosurface(f"tr_{s}_n", f"mt_{s}", -ISO)
    cmd.isosurface(f"df_{s}_p", f"md_{s}",  ISO)                 # difference density: green/violet
    cmd.isosurface(f"df_{s}_n", f"md_{s}", -ISO)
cmd.color("red","tr_*_p"); cmd.color("blue","tr_*_n")
cmd.color("green","df_*_p"); cmd.color("violet","df_*_n")
cmd.group("Transition_Density","tr_*"); cmd.group("Difference_Density","df_*")

def render(fn, kind, view):
    cmd.disable("Transition_Density"); cmd.disable("Difference_Density")
    cmd.enable("Transition_Density" if kind=="tr" else "Difference_Density")
    if view == "dimer":
        cmd.enable("siteA"); cmd.enable("siteB"); cmd.show("sticks","dimer_ref and backbone")
        cmd.enable(f"{kind}_A_p"); cmd.enable(f"{kind}_A_n")
        cmd.orient("dimer_ref")
    else:                                                        # site-B QM close-up
        cmd.hide("sticks", "dimer_ref"); cmd.disable("siteA")
        cmd.disable(f"{kind}_A_p"); cmd.disable(f"{kind}_A_n")
        cmd.orient(f"siteB and {QM}"); cmd.zoom(f"siteB and {QM}", 3)
    cmd.ray(W, H); cmd.png(f"{OUT}/{fn}", dpi=300); print("wrote", fn)

render("01_transition_dimer.png", "tr", "dimer")
render("02_transition_qm.png",    "tr", "qm")
render("03_difference_dimer.png", "df", "dimer")
render("04_difference_qm.png",    "df", "qm")

# clean state for interactive viewing
cmd.enable("siteA"); cmd.enable("siteB"); cmd.show("sticks","dimer_ref and backbone")
cmd.enable("Transition_Density"); cmd.disable("Difference_Density"); cmd.orient("dimer_ref")
print("DONE — 4 paper-format figures in", OUT)
