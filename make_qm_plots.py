# make_qm_plots.py  — run: pymol -cq make_qm_plots.py
# Figures for the in-protein STEOM/QM-region result:
#   (1) qm_overview.png  : chromophore + Tyr203 pi-stack + His148 in the protein barrel
#   (2) qm_pistack.png   : close-up of the CR2 phenolate / Tyr203 phenol pi-stack + His148 H-bond
# AMBER numbering in monomer_min.pdb: CR2=66, Tyr203=202(pi-stack), His148=147, Ser205=204.
from pymol import cmd

PDB = "/home/robson/PetaChem/anionic_build/monomer_min.pdb"
OUT = "/home/robson/PetaChem/plots"
import os; os.makedirs(OUT, exist_ok=True)

cmd.reinitialize()
cmd.load(PDB, "prot")
cmd.remove("solvent")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_transparency", 0.75)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("stick_radius", 0.16)
cmd.set("valence", 0)

cmd.select("chromo", "resi 66")          # CR2 chromophore (QM)
cmd.select("tyr", "resi 202")            # Tyr203 pi-stack (QM, sidechain)
cmd.select("his", "resi 147")            # His148 (MM)
cmd.select("ser", "resi 204")            # Ser205

cmd.hide("everything")
cmd.show("cartoon", "prot")
cmd.color("gray80", "prot")

for sel, fn in [("chromo", cmd.util.cbac), ("tyr", cmd.util.cbao), ("his", cmd.util.cbag)]:
    cmd.show("sticks", sel)
    fn(sel)                               # cbac=cyan C, cbao=orange C, cbag=green C
cmd.set("stick_radius", 0.22, "chromo")
cmd.hide("everything", "elem H and not (neighbor elem N+O)")   # keep only polar H

# H-bond His148...phenolate-O  and any chromophore-Tyr/Ser polar contacts
cmd.distance("hb_his", "his and elem N", "chromo and elem O", 3.4, mode=2)
cmd.distance("hb_ser", "ser and elem O", "chromo and elem O+N", 3.4, mode=2)
cmd.color("black", "hb_his hb_ser")
cmd.set("dash_gap", 0.3); cmd.set("dash_width", 3)
cmd.hide("labels")

# ---- Figure 1: overview ----
cmd.orient("chromo or tyr or his")
cmd.zoom("chromo or tyr or his", 6)
cmd.turn("y", 20)
cmd.set("ray_trace_mode", 1)
cmd.ray(1600, 1200)
cmd.png(f"{OUT}/qm_overview.png", dpi=300)
print("wrote qm_overview.png")

# ---- Figure 2: pi-stack close-up (drop the cartoon) ----
cmd.hide("cartoon")
cmd.show("sticks", "chromo or tyr or his")
cmd.orient("chromo and (resn CR2 and name CG+CD1+CD2+CE1+CE2+CZ+OH)")  # phenolate ring face-on
cmd.zoom("chromo or tyr", 2.5)
cmd.ray(1600, 1200)
cmd.png(f"{OUT}/qm_pistack.png", dpi=300)
print("wrote qm_pistack.png")
print("DONE")
