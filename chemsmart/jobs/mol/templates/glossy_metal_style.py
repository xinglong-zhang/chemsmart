"""
Glossy semi-metallic PyMOL style for publication and presentation figures.

Load structures, then run for example::

    glossy_complex all, white
    glossy_complex all, dark
"""

from pymol import cmd


def _safe_set(setting, value, selection=None):
    """
    Set a PyMOL setting without crashing on unsupported versions.
    """
    try:
        if selection is None:
            cmd.set(setting, value)
        else:
            cmd.set(setting, value, selection)
    except Exception:
        pass


def _define_colors():
    """Define soft high-contrast colors for glossy rendering."""
    colors = {
        "glossy_carbon": [0.72, 0.72, 0.72],
        "glossy_hydrogen": [0.96, 0.96, 0.96],
        "glossy_nitrogen": [0.20, 0.32, 0.95],
        "glossy_oxygen": [0.95, 0.08, 0.06],
        "glossy_sulfur": [0.95, 0.76, 0.18],
        "glossy_phosphorus": [1.00, 0.55, 0.15],
        "glossy_halogen": [0.20, 0.85, 0.35],
        "metal_mn_gold": [0.92, 0.62, 0.25],
        "metal_mn_pink": [0.86, 0.48, 0.78],
        "metal_silver": [0.82, 0.82, 0.86],
        "metal_copper": [0.95, 0.48, 0.24],
        "label_dark": [0.02, 0.02, 0.02],
        "label_light": [1.00, 1.00, 1.00],
    }

    for name, rgb in colors.items():
        try:
            cmd.set_color(name, rgb)
        except Exception:
            pass


def default_render(background="white", transparent=None, orthoscopic="on"):
    """
    High-quality glossy / semi-metallic rendering style.

    Parameters
    ----------
    background : str
        "white" for publication figures, or "dark" for high-contrast slides.
    transparent : bool or None
        If None, white mode uses transparent background, dark mode opaque.
    orthoscopic : str
        "on" gives a flatter graphical look; "off" gives perspective.
    """
    _define_colors()

    bg = str(background).lower()
    is_dark = bg in ["dark", "black", "presentation", "slide", "slides"]

    if transparent is None:
        transparent = not is_dark

    if is_dark:
        cmd.bg_color("black")
    else:
        cmd.bg_color("white")

    _safe_set("ray_opaque_background", 0 if transparent else 1)
    _safe_set("ray_shadow", 1)
    _safe_set("ray_trace_mode", 0)
    _safe_set("ray_trace_gain", 0.08)
    _safe_set("ray_shadow_decay_factor", 0.25)
    _safe_set("ray_shadow_decay_range", 2.0)
    _safe_set("ambient_occlusion_mode", 1)
    _safe_set("ambient_occlusion_scale", 18)
    _safe_set("ambient_occlusion_smooth", 12)
    _safe_set("ambient", 0.18 if is_dark else 0.25)
    _safe_set("direct", 0.82)
    _safe_set("reflect", 0.55)
    _safe_set("specular", 0.85)
    _safe_set("spec_reflect", 0.70)
    _safe_set("spec_power", 350)
    _safe_set("shininess", 90)
    _safe_set("light_count", 8)
    _safe_set("two_sided_lighting", 1)
    _safe_set("use_shaders", 1)
    _safe_set("antialias", 2)
    _safe_set("ray_trace_antialias", 2)
    _safe_set("sphere_quality", 3)
    _safe_set("stick_quality", 30)
    _safe_set("cartoon_sampling", 14)

    if str(orthoscopic).lower() in ["on", "1", "true", "yes"]:
        _safe_set("orthoscopic", 1)
    else:
        _safe_set("orthoscopic", 0)

    _safe_set("field_of_view", 45)
    _safe_set("depth_cue", 1)
    _safe_set("fog_start", 0.28 if is_dark else 0.55)
    _safe_set("label_font_id", 7)
    _safe_set("label_size", 24)
    _safe_set("label_position", [0, 0, 0])
    _safe_set("label_connector", 0)

    if is_dark:
        _safe_set("label_color", "white")
        _safe_set("label_outline_color", "black")
    else:
        _safe_set("label_color", "black")
        _safe_set("label_outline_color", "white")

    _safe_set("label_shadow_mode", 2)
    cmd.rebuild()


def color_by_element(selection="all", metal_color="metal_mn_gold"):
    """Apply a clean element color scheme."""
    _define_colors()

    sel = f"({selection})"

    cmd.color("glossy_carbon", f"{sel} and elem C")
    cmd.color("glossy_hydrogen", f"{sel} and elem H")
    cmd.color("glossy_nitrogen", f"{sel} and elem N")
    cmd.color("glossy_oxygen", f"{sel} and elem O")
    cmd.color("glossy_sulfur", f"{sel} and elem S")
    cmd.color("glossy_phosphorus", f"{sel} and elem P")
    cmd.color("glossy_halogen", f"{sel} and elem F+Cl+Br+I")
    cmd.color(metal_color, f"{sel} and elem Mn")
    cmd.color("metal_silver", f"{sel} and elem Fe+Co+Ni+Ru+Rh+Pd+Ir+Pt")
    cmd.color("metal_copper", f"{sel} and elem Cu")


def style_metal_complex(
    selection="all",
    metal_sel="elem Mn",
    coord_sel=None,
    coord_cutoff=2.6,
    donor_elements="N+O+S+P+H",
    show_all_hydrogens="off",
    metal_color="metal_mn_gold",
):
    """Clean ligand-stick + central ball/sphere style."""
    _define_colors()

    sel = f"({selection})"
    metal = f"({sel}) and ({metal_sel})"

    if coord_sel is None:
        coord = (
            f"({sel}) and "
            f"(elem {donor_elements}) within {float(coord_cutoff)} of ({metal})"
        )
    else:
        coord = f"({sel}) and ({coord_sel})"

    core = f"({metal}) or ({coord})"

    cmd.hide("everything", sel)
    cmd.show("sticks", sel)
    _safe_set("stick_radius", 0.105, sel)
    _safe_set("stick_quality", 30, sel)

    if str(show_all_hydrogens).lower() not in ["on", "1", "true", "yes"]:
        cmd.hide("sticks", f"{sel} and elem H and not ({coord})")

    cmd.show("spheres", core)
    _safe_set("sphere_scale", 0.62, metal)
    _safe_set("sphere_scale", 0.40, coord)
    _safe_set("sphere_scale", 0.32, f"({coord}) and elem H")
    _safe_set("stick_ball", 0, sel)

    color_by_element(sel, metal_color=metal_color)
    _safe_set("stick_radius", 0.14, f"{sel} and within 2.8 of ({metal})")
    cmd.rebuild()


def make_internal_labels(
    selection,
    label_obj="internal_atom_labels",
    label_mode="element",
    label_color="auto",
    outline_color="auto",
    size=24,
    font_id=7,
):
    """Create flat centered atom labels using pseudoatoms."""
    _define_colors()
    cmd.delete(label_obj)

    model = cmd.get_model(selection)
    if len(model.atom) == 0:
        print(f"No atoms found for label selection: {selection}")
        return

    for i, atom in enumerate(model.atom, start=1):
        elem = atom.symbol.strip()
        name = atom.name.strip()
        text = (
            name
            if label_mode.lower() in ["name", "atomname", "atom_name"]
            else elem
        )
        cmd.pseudoatom(
            object=label_obj,
            name=f"L{i}",
            pos=atom.coord,
            label=text,
        )

    cmd.hide("everything", label_obj)
    cmd.show("labels", label_obj)
    _safe_set("label_font_id", int(font_id), label_obj)
    _safe_set("label_size", float(size), label_obj)
    _safe_set("label_position", [0, 0, 0], label_obj)
    _safe_set("label_connector", 0, label_obj)
    _safe_set("label_shadow_mode", 2, label_obj)

    if label_color == "auto":
        label_color = "white"
    if outline_color == "auto":
        outline_color = "black"

    _safe_set("label_color", label_color, label_obj)
    _safe_set("label_outline_color", outline_color, label_obj)
    cmd.rebuild()


def label_coordination_core(
    selection="all",
    metal_sel="elem Mn",
    coord_sel=None,
    coord_cutoff=2.6,
    donor_elements="N+O+S+P+H",
    label_obj="coordination_core_labels",
    label_color="auto",
    outline_color="auto",
    size=24,
):
    """Label the central metal and coordinating atoms."""
    sel = f"({selection})"
    metal = f"({sel}) and ({metal_sel})"

    if coord_sel is None:
        coord = (
            f"({sel}) and "
            f"(elem {donor_elements}) within {float(coord_cutoff)} of ({metal})"
        )
    else:
        coord = f"({sel}) and ({coord_sel})"

    label_selection = f"({metal}) or ({coord})"
    make_internal_labels(
        label_selection,
        label_obj=label_obj,
        label_mode="element",
        label_color=label_color,
        outline_color=outline_color,
        size=size,
    )


def glossy_complex(
    selection="all",
    background="white",
    metal_sel="elem Mn",
    coord_sel=None,
    coord_cutoff=2.6,
    donor_elements="N+O+S+P+H",
    label_core="on",
    metal_color="metal_mn_gold",
):
    """
    Apply the full glossy semi-metallic molecular illustration style.

    Example
    -------
    glossy_complex all, dark
    glossy_complex all, white
    """
    bg = str(background).lower()
    is_dark = bg in ["dark", "black", "presentation", "slide", "slides"]

    default_render(background=background)

    style_metal_complex(
        selection=selection,
        metal_sel=metal_sel,
        coord_sel=coord_sel,
        coord_cutoff=coord_cutoff,
        donor_elements=donor_elements,
        show_all_hydrogens="off",
        metal_color=metal_color,
    )

    if str(label_core).lower() in ["on", "1", "true", "yes"]:
        label_color = "white" if is_dark else "black"
        outline_color = "black" if is_dark else "white"
        label_coordination_core(
            selection=selection,
            metal_sel=metal_sel,
            coord_sel=coord_sel,
            coord_cutoff=coord_cutoff,
            donor_elements=donor_elements,
            label_color=label_color,
            outline_color=outline_color,
            size=24,
        )

    cmd.zoom(selection, buffer=2.0)
    cmd.orient(selection)
    cmd.rebuild()


def render_png(filename="render.png", width=2400, height=1800, ray="on"):
    """Render a high-resolution PNG."""
    if str(ray).lower() in ["on", "1", "true", "yes"]:
        cmd.ray(int(width), int(height))
    cmd.png(filename, dpi=300)


cmd.extend("default_render", default_render)
cmd.extend("style_metal_complex", style_metal_complex)
cmd.extend("make_internal_labels", make_internal_labels)
cmd.extend("label_coordination_core", label_coordination_core)
cmd.extend("glossy_complex", glossy_complex)
cmd.extend("render_png", render_png)
