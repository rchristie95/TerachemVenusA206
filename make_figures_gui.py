# make_figures_gui.py — interactive STEOM density scene (paper format), no auto-render.
from pymol import cmd, util

DIMER = "/home/robson/PetaChem/venus_dimer.pdb"
MONO  = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"
TRANS = "/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.cube"
DIFF  = "/home/robson/PetaChem/neo_model/orca_steom/steom_diffdens.cube"
ISO, QM = 0.04, "resi 66"

cmd.reinitialize()
cmd.bg_color("white")
for k, v in [("ray_opaque_background",0),("ray_shadows",0),("antialias",2),("ambient",0.45),
             ("transparency",0.35),("cartoon_transparency",0.65),("cartoon_fancy_helices",1),
             ("depth_cue",0),("valence",0),("sphere_scale",0.22)]:
    cmd.set(k, v)

cmd.load(DIMER, "dimer_ref")
cmd.load(MONO, "siteA"); cmd.load(MONO, "siteB")
cmd.super("siteA and chain A", "dimer_ref and chain A")
cmd.super("siteB and chain A", "dimer_ref and chain B")

cmd.hide("everything")
cmd.show("cartoon", "dimer_ref")
cmd.spectrum("count", "rainbow", "dimer_ref and name CA")

for s in ("siteA", "siteB"):
    cmd.show("sticks",  f"{s} and {QM}")
    cmd.show("spheres", f"{s} and {QM}")
util.cbaw(f"(siteA or siteB) and {QM}")
cmd.color("yellow", f"(siteA or siteB) and {QM} and elem H")

cmd.load(TRANS, "mt"); cmd.load(DIFF, "md")
for src, dst, site in [("mt", "mt_A", "siteA"), ("mt", "mt_B", "siteB"),
                       ("md", "md_A", "siteA"), ("md", "md_B", "siteB")]:
    cmd.copy(dst, src); cmd.matrix_copy(site, dst)
cmd.disable("mt"); cmd.disable("md")

for s in ("A", "B"):
    cmd.isosurface(f"tr_{s}_p", f"mt_{s}",  ISO)
    cmd.isosurface(f"tr_{s}_n", f"mt_{s}", -ISO)
    cmd.isosurface(f"df_{s}_p", f"md_{s}",  ISO)
    cmd.isosurface(f"df_{s}_n", f"md_{s}", -ISO)
cmd.color("red", "tr_*_p");   cmd.color("blue", "tr_*_n")
cmd.color("green", "df_*_p"); cmd.color("violet", "df_*_n")
cmd.group("Transition_Density", "tr_*"); cmd.group("Difference_Density", "df_*")

cmd.enable("Transition_Density"); cmd.disable("Difference_Density")
cmd.orient("dimer_ref")
print("=" * 60)
print("STEOM density scene ready (paper format).")
print(" Transition_Density (red/blue) [on] | Difference_Density (green/violet)")
print(" zoom a monomer:  zoom siteB and resi 66")
print(" 4 PNGs already saved in plots_steom/")
print("=" * 60)
