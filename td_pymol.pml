load /home/robson/PetaChem/neo_model/orca_steom/steom_qm.xyz, mol
load /home/robson/PetaChem/neo_model/orca_steom/steom_transdens.cube, tdmap
hide everything
show sticks, mol
set stick_radius, 0.14, mol
color gray40, mol and elem C
color red,    mol and elem O
color blue,   mol and elem N
color gray90, mol and elem H
isosurface td_pos, tdmap,  0.04
isosurface td_neg, tdmap, -0.04
set transparency, 0.30
color marine,    td_pos
color firebrick, td_neg
bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set ambient, 0.5
set specular, 0.15
orient mol
ray 1800, 1350
png /home/robson/PetaChem/neo_model/orca_steom/steom_transition_density.png, dpi=300
