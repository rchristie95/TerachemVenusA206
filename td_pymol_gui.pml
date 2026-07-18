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
set ray_shadows, 0
orient mol
