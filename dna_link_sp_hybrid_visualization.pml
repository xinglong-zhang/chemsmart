unset stick_color, all
hide everything, all
show sticks, all
set_color light_C, [3, 3, 3]
set_color light_N, [0.6, 0.8, 1.0]
set_color light_O, [1.0, 0.7, 0.7]
set_color light_P, [1.0, 0.85, 0.6]
set_color light_S, [1.0, 0.7, 0.7]
color light_C, elem C
color light_P, elem P
color light_O, elem O
color light_N, elem N
color light_S, elem S
select group1, id 1 or id 2 or id 3
util.cbap group1
set stick_transparency, 0, all
set stick_radius, 0.25, (group1)
show surface, all
set surface_color, grey, all
set transparency, 0.8, all
