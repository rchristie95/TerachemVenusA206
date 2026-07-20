# make_steom_paper_figs.py — STEOM 44-atom density figures in the paper's EXACT style.
# Dimer view: black bg, whole protein as rainbow lines, white-stick QM region, transparent density.
# QM close-up: white bg, white-stick QM region + density. Matches visualise_dimer.pml conventions.
from pymol import cmd, util
import os

B     = "/home/robson/PetaChem"
DIMER = f"{B}/venus_dimer.pdb"
MONO  = f"{B}/tc_simple_anionic/monomer_relaxed.pdb"
QMXYZ = f"{B}/neo_model/orca_steom/steom_qm.xyz"          # the honest 44 atoms STEOM computed
TRANS = f"{B}/neo_model/orca_steom/steom_transdens.cube"
DIFF  = f"{B}/neo_model/orca_steom/steom_diffdens.cube"
OUT   = f"{B}/plots_steom"; os.makedirs(OUT, exist_ok=True)
ISO, W, H = 0.04, 1600, 1200

cmd.reinitialize()
for k, v in [("ray_shadows",0),("antialias",2),("transparency",0.5),("depth_cue",0),
             ("line_width",1.0),("sphere_scale",0.25),("valence",0)]:
    cmd.set(k, v)

cmd.load(DIMER, "dimer_ref")
cmd.load(MONO, "siteA"); cmd.load(MONO, "siteB")
cmd.super("siteA and chain A", "dimer_ref and chain A")
cmd.super("siteB and chain A", "dimer_ref and chain B")

cmd.load(QMXYZ, "qmA"); cmd.load(QMXYZ, "qmB")            # 44-atom region -> each site
cmd.matrix_copy("siteA", "qmA"); cmd.matrix_copy("siteB", "qmB")

cmd.load(TRANS, "mt"); cmd.load(DIFF, "md")               # densities -> each site
for src, dst, site in [("mt","mt_A","siteA"),("mt","mt_B","siteB"),
                       ("md","md_A","siteA"),("md","md_B","siteB")]:
    cmd.copy(dst, src); cmd.matrix_copy(site, dst)
cmd.disable("mt"); cmd.disable("md")

cmd.hide("everything")
for q in ("qmA","qmB"):
    cmd.show("sticks", q); cmd.show("spheres", q)
util.cbaw("qmA or qmB"); cmd.color("yellow", "(qmA or qmB) and elem H")

for s, q in [("A","qmA"),("B","qmB")]:                    # carve 6 A around the QM region
    cmd.isosurface(f"tr_{s}_p", f"mt_{s}",  ISO, q, carve=6.0)
    cmd.isosurface(f"tr_{s}_n", f"mt_{s}", -ISO, q, carve=6.0)
    cmd.isosurface(f"df_{s}_p", f"md_{s}",  ISO, q, carve=6.0)
    cmd.isosurface(f"df_{s}_n", f"md_{s}", -ISO, q, carve=6.0)
cmd.color("red","tr_*_p");   cmd.color("blue","tr_*_n")
cmd.color("green","df_*_p"); cmd.color("magenta","df_*_n")
cmd.group("Transition_Density","tr_*"); cmd.group("Difference_Density","df_*")

def enable_only(kind):
    cmd.disable("Transition_Density"); cmd.disable("Difference_Density")
    cmd.enable("Transition_Density" if kind == "tr" else "Difference_Density")

def dimer_view(kind, fn):                                 # black bg, rainbow lines
    cmd.bg_color("black"); cmd.set("ray_opaque_background", 1)
    enable_only(kind); cmd.enable(f"{kind}_A_p"); cmd.enable(f"{kind}_A_n")
    cmd.enable("qmA"); cmd.enable("qmB")
    cmd.hide("everything","dimer_ref"); cmd.show("lines","dimer_ref")
    cmd.spectrum("count","rainbow","dimer_ref")
    cmd.orient("dimer_ref"); cmd.ray(W, H); cmd.png(f"{OUT}/{fn}", dpi=300); print("wrote", fn)

def qm_view(kind, fn):                                    # white bg, QM close-up (site B)
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 0)
    enable_only(kind); cmd.hide("everything","dimer_ref"); cmd.disable("qmA")
    cmd.disable(f"{kind}_A_p"); cmd.disable(f"{kind}_A_n"); cmd.enable("qmB")
    cmd.orient("qmB"); cmd.zoom("qmB", 2)
    cmd.ray(W, H); cmd.png(f"{OUT}/{fn}", dpi=300); print("wrote", fn)
    cmd.enable("qmA"); cmd.enable(f"{kind}_A_p"); cmd.enable(f"{kind}_A_n")

dimer_view("tr", "01_transition_dimer.png")
qm_view("tr",   "02_transition_qm.png")
dimer_view("df", "03_difference_dimer.png")
qm_view("df",   "04_difference_qm.png")

cmd.bg_color("black"); cmd.set("ray_opaque_background", 1)
enable_only("tr"); cmd.show("lines","dimer_ref"); cmd.spectrum("count","rainbow","dimer_ref")
cmd.enable("qmA"); cmd.enable("qmB"); cmd.orient("dimer_ref")
print("STEOM PAPER-STYLE FIGS DONE")
