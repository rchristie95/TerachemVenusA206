# make_steom_gui.py — interactive paper-style STEOM scene (no auto-render).
from pymol import cmd, util

B     = "/home/robson/PetaChem"
DIMER = f"{B}/venus_dimer.pdb"
MONO  = f"{B}/tc_simple_anionic/monomer_relaxed.pdb"
QMXYZ = f"{B}/neo_model/orca_steom/steom_qm.xyz"
TRANS = f"{B}/neo_model/orca_steom/steom_transdens.cube"
DIFF  = f"{B}/neo_model/orca_steom/steom_diffdens.cube"
ISO = 0.04

cmd.reinitialize()
for k, v in [("ray_shadows",0),("antialias",1),("transparency",0.5),("depth_cue",0),
             ("line_width",1.0),("sphere_scale",0.25),("valence",0),("surface_quality",0)]:
    cmd.set(k, v)

cmd.load(DIMER, "dimer_ref")
cmd.load(MONO, "siteA"); cmd.load(MONO, "siteB")
cmd.super("siteA and chain A", "dimer_ref and chain A")
cmd.super("siteB and chain A", "dimer_ref and chain B")

cmd.load(QMXYZ, "qmA"); cmd.load(QMXYZ, "qmB")
cmd.matrix_copy("siteA", "qmA"); cmd.matrix_copy("siteB", "qmB")

cmd.load(TRANS, "mt"); cmd.load(DIFF, "md")
for src, dst, site in [("mt","mt_A","siteA"),("mt","mt_B","siteB"),
                       ("md","md_A","siteA"),("md","md_B","siteB")]:
    cmd.copy(dst, src); cmd.matrix_copy(site, dst)
cmd.disable("mt"); cmd.disable("md")

cmd.hide("everything")
for q in ("qmA","qmB"):
    cmd.show("sticks", q); cmd.show("spheres", q)
util.cbaw("qmA or qmB"); cmd.color("yellow", "(qmA or qmB) and elem H")

for s, q in [("A","qmA"),("B","qmB")]:
    cmd.isosurface(f"tr_{s}_p", f"mt_{s}",  ISO, q, carve=6.0)
    cmd.isosurface(f"tr_{s}_n", f"mt_{s}", -ISO, q, carve=6.0)
    cmd.isosurface(f"df_{s}_p", f"md_{s}",  ISO, q, carve=6.0)
    cmd.isosurface(f"df_{s}_n", f"md_{s}", -ISO, q, carve=6.0)
cmd.color("red","tr_*_p");   cmd.color("blue","tr_*_n")
cmd.color("green","df_*_p"); cmd.color("magenta","df_*_n")
cmd.group("Transition_Density","tr_*"); cmd.group("Difference_Density","df_*")

cmd.bg_color("black"); cmd.set("ray_opaque_background", 1)
cmd.show("lines", "dimer_ref"); cmd.spectrum("count", "rainbow", "dimer_ref")
cmd.enable("Transition_Density"); cmd.disable("Difference_Density")
cmd.orient("dimer_ref")
print("STEOM PAPER-STYLE SCENE READY")
