#####################################
# Viewport & General Settings
#####################################
viewport 1356, 1538

#####################################
# Stick Settings
#####################################
hide everything
show sticks, all
set stick_transparency, 0, all
set stick_radius, 0.25, (sub1 or sub2 or sub3 or sub4)

#####################################
# Custom Colors
#####################################
set_color light_C, [0.8, 0.8, 0.9]   # light gray-blue for carbon
set_color light_N, [0.6, 0.8, 1.0]   # light blue for nitrogen
set_color light_O, [1.0, 0.7, 0.7]   # light red/pink for oxygen
set_color light_P, [1.0, 0.85, 0.6]  # light orange for phosphorus

#####################################
# Apply Colors
#####################################
color light_C, elem_C
color light_P, elem_P
color light_O, elem_O
color light_H, elem_H
color light_N, elem_N
#####################################
# Ligands Display
#####################################
hide everything, ({groups})
show sticks, ({groups})
unset stick_color, ({groups})

# Distinguish ligands with util color schemes
util.cbap sub1
util.cbac sub2
util.cbay sub3
util.cbag sub4
#####################################
# Labels
#####################################
set label_digits,2

#####################################
# Surface Settings (applied last)
#####################################
show surface, all
set surface_color, grey, all
set transparency, 0.7, all
set sphere_scale, 0.2

