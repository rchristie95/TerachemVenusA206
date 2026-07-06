# make_toc.py — graphical-abstract (TOC) render for the JPC B revision.
# Landscape ACS TOC aspect (~1.85:1): the two Venus beta-barrels as translucent
# yellow cartoon, each chromophore's 44-atom STEOM QM region as sticks, with the
# STEOM transition density (red/blue) showing the delocalised, coupled excitation.
from pymol import cmd, util
import os

B     = "/home/robson/PetaChem"
DIMER = f"{B}/tandem_dimer.pdb"
MONO  = f"{B}/tc_simple_anionic/monomer_relaxed.pdb"
QMXYZ = f"{B}/neo_model/orca_steom/steom_qm.xyz"
TRANS = f"{B}/neo_model/orca_steom/steom_transdens.cube"
OUT   = f"{B}/plots_steom"; os.makedirs(OUT, exist_ok=True)
ISO, W, H = 0.04, 2000, 1080          # 1.852:1 ~ ACS TOC 3.25 x 1.75 in

cmd.reinitialize()
for k, v in [("ray_shadows", 0), ("antialias", 2), ("transparency", 0.4),
             ("cartoon_transparency", 0.55), ("cartoon_fancy_helices", 1),
             ("depth_cue", 0), ("sphere_scale", 0.22), ("valence", 0),
             ("ray_opaque_background", 1)]:
    cmd.set(k, v)
cmd.bg_color("white")

cmd.load(DIMER, "dimer_ref")
cmd.load(MONO, "siteA"); cmd.load(MONO, "siteB")
cmd.align("siteA and name CA", "dimer_ref and resi 1-229 and name CA")
cmd.align("siteB and name CA", "dimer_ref and resi 263-491 and name CA")

cmd.load(QMXYZ, "qmA"); cmd.load(QMXYZ, "qmB")
cmd.matrix_copy("siteA", "qmA"); cmd.matrix_copy("siteB", "qmB")

cmd.load(TRANS, "mt")
for dst, site in [("mt_A", "siteA"), ("mt_B", "siteB")]:
    cmd.copy(dst, "mt"); cmd.matrix_copy(site, dst)
cmd.disable("mt")

cmd.hide("everything")
cmd.show("lines", "dimer_ref")
cmd.spectrum("count", "rainbow", "dimer_ref and resi 1-240 and name CA")
cmd.spectrum("count", "rainbow", "dimer_ref and resi 250-500 and name CA")
for q in ("qmA", "qmB"):
    cmd.show("sticks", q); cmd.show("spheres", q)
util.cbaw("qmA or qmB")
cmd.color("gray90", "(qmA or qmB) and elem C")
cmd.color("yellow", "(qmA or qmB) and elem H")
cmd.set("stick_radius", 0.13, "qmA or qmB")
cmd.set("sphere_scale", 0.22, "qmA or qmB")

for s in ("A", "B"):                                   # STEOM transition density, carved to QM
    cmd.isosurface(f"tr_{s}_p", f"mt_{s}",  ISO, f"qm{s}", carve=5.0)
    cmd.isosurface(f"tr_{s}_n", f"mt_{s}", -ISO, f"qm{s}", carve=5.0)
cmd.color("red", "tr_*_p"); cmd.color("blue", "tr_*_n")

cmd.orient("dimer_ref")                                # long (barrel-to-barrel) axis -> horizontal
cmd.zoom("dimer_ref", 3)
cmd.ray(W, H)
cmd.png(f"{OUT}/toc_steom_white.png", dpi=300)
print("wrote toc_steom_white.png", W, "x", H)
