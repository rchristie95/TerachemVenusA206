# Fig. S1: the 44-atom QM region at three orientations.
# Answers reviewer 2, minor point 2 ("I could not figure out how to get a QM
# region with 44 atoms"). Rendered from the exact geometry used in every
# site-energy calculation, neo_model/orca_steom/geom_cthrp.xyz.
reinitialize
set bg_rgb, [1,1,1]
set ray_opaque_background, off
set antialias, 4
set ray_trace_mode, 1
set ray_trace_gain, 0.15
set stick_radius, 0.13
set sphere_scale, 0.19
set label_size, 15
set label_color, black
set label_font_id, 7
set depth_cue, 0
set specular, 0.2

load /home/robson/PetaChem/neo_model/orca_steom/geom_cthrp.xyz, qm

hide everything
show sticks, qm
show spheres, qm
util.cbaw("qm")
color grey30, qm and elem C
color firebrick, qm and elem O
color steelblue, qm and elem N
color grey85, qm and elem H

# The three link hydrogens are the last three atoms of the file; mark them.
select links, qm and id 42-44
color limegreen, links
set sphere_scale, 0.30, links

orient qm
python
from pymol import cmd
views = {"a": (0, 0, 0), "b": (0, 90, 0), "c": (90, 0, 0)}
for tag, (rx, ry, rz) in views.items():
    cmd.orient("qm")
    if rx: cmd.turn("x", rx)
    if ry: cmd.turn("y", ry)
    if rz: cmd.turn("z", rz)
    cmd.zoom("qm", buffer=1.2)
    cmd.png(f"/home/robson/PetaChem/manuscript/FigS1_qmregion_{tag}.png",
            width=1400, height=1400, dpi=300, ray=1)
    print(f"rendered {tag}", flush=True)
python end
quit
