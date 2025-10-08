from pymol import cmd


def setup_hybrid_visualization(sub1=None, sub2=None, sub3=None, sub4=None):
    """
    Configure PyMOL viewport, sticks, colors, ligands, labels, and surfaces.

    Args:
        sub1, sub2, sub3, sub4, po: Optional selection strings for ligands/substrates.
    """
    #####################################
    # Viewport & General Settings
    #####################################
    cmd.viewport(1356, 1538)

    #####################################
    # Stick Settings
    #####################################
    # Build dynamic selection string
    subs = [s for s in [sub1, sub2, sub3, sub4] if s]
    selection_str = " or ".join(subs) if subs else "all"

    cmd.set("stick_transparency", 0, "all")
    cmd.set("stick_radius", 0.25, selection_str)

    #####################################
    # Custom Colors
    #####################################
    cmd.set_color("light_C", [0.8, 0.8, 0.9])  # light gray-blue for carbon
    cmd.set_color("light_N", [0.6, 0.8, 1.0])  # light blue for nitrogen
    cmd.set_color("light_O", [1.0, 0.7, 0.7])  # light red/pink for oxygen
    cmd.set_color("light_P", [1.0, 0.85, 0.6])  # light orange for phosphorus

    #####################################
    # Apply Colors
    #####################################
    cmd.color("light_C", "elem_C")
    cmd.color("light_P", "elem_P")
    cmd.color("light_O", "elem_O")
    cmd.color("light_H", "elem_H")
    cmd.color("light_N", "elem_N")

    #####################################
    # Ligands Display
    #####################################
    ligands = [s for s in [sub1, sub2, sub3, sub4] if s]
    lig_str = " or ".join(ligands) if ligands else "none"
    cmd.hide("everything", lig_str)
    cmd.show("sticks", lig_str)
    cmd.unset("stick_color", lig_str)

    # Distinguish ligands with util color schemes
    if sub1:
        cmd.util.cbap(sub1)
    if sub2:
        cmd.util.cbac(sub2)
    if sub3:
        cmd.util.cbay(sub3)
    if sub4:
        cmd.util.cbag(sub4)

    #####################################
    # Labels
    #####################################
    # cmd.set("label_size", 25, "h_bond")  # optional
    cmd.set("label_digits", 2)

    #####################################
    # Surface Settings (applied last)
    #####################################
    cmd.show("surface", "all")
    cmd.set("surface_color", "grey", "all")
    cmd.set("transparency", 0.7, "all")
    cmd.show("sphere", "id 413")
    cmd.set("sphere_scale", 0.2, "id 413")
