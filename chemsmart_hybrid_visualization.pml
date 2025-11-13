pymol_style all
unset stick_color, all
hide everything, all
show sticks, all
set_color light_C, [0.8, 0.8, 0.9]
set_color light_N, [0.6, 0.8, 1.0]
set_color light_O, [1.0, 0.7, 0.7]
set_color light_P, [1.0, 0.85, 0.6]
color light_C, elem C
color light_P, elem P
color light_O, elem O
color light_N, elem N
select group1,  id 233 or id 468-512
select group2,  id 308 or id 397-414 or id 416-423
util.cbap group1
util.cbac group2
set stick_transparency, 0, all
set stick_radius, 0.25, (group1 or group2)
show surface, all
set surface_color, grey, all
set transparency, 0.7, all
