"""
Glossy metallic-poster PyMOL style for publication and presentation figures.

ChemSmart applies this template when ``visualize -s glossy`` is used.
In PyMOL directly::

    run glossy_metal_style.py
    metallic_poster_render all
    metallic_poster_render all, elem Mn, None, 2.6, N+O+S+P+H, dark
"""

from pymol import cmd
from pymol_style_imports import load_metallic_element_colors

_metallic_colors = load_metallic_element_colors()
define_metallic_poster_colors = _metallic_colors.define_metallic_poster_colors
apply_metallic_poster_element_colors = (
    _metallic_colors.apply_metallic_poster_element_colors
)
DEFAULT_METAL_COLOR = _metallic_colors.DEFAULT_METAL_COLOR


def _safe_set(setting, value, selection=None):
    """Set PyMOL setting safely across different PyMOL versions."""
    try:
        if selection is None:
            cmd.set(setting, value)
        else:
            cmd.set(setting, value, selection)
    except Exception:
        pass


def _normalize_none(value):
    """Convert PyMOL command-line 'None' strings to Python None."""
    if value is None:
        return None
    if str(value).lower() in ["none", "null", "false", ""]:
        return None
    return value


def _make_centered_element_labels(
    atoms_to_label,
    label_prefix="metallic_poster_labels",
    label_size=24,
    label_font_id=7,
):
    """
    Create centered pseudoatom labels at the exact coordinates of selected atoms.

    This is the most reliable PyMOL method for the 'letter inside sphere'
    or 'circled atom label' effect.
    """
    define_metallic_poster_colors()

    cmd.delete(label_prefix + "*")

    model = cmd.get_model(atoms_to_label)

    if len(model.atom) == 0:
        print("No atoms found for centered labels:", atoms_to_label)
        return

    label_color_by_element = {
        "Mn": "poster_label_black",
        "Fe": "poster_label_black",
        "Co": "poster_label_black",
        "Ni": "poster_label_black",
        "Cu": "poster_label_black",
        "Zn": "poster_label_black",
        "Ru": "poster_label_black",
        "Rh": "poster_label_black",
        "Pd": "poster_label_black",
        "Ir": "poster_label_black",
        "Pt": "poster_label_black",
        "N": "poster_label_white",
        "O": "poster_label_white",
        "S": "poster_label_black",
        "P": "poster_label_black",
        "H": "poster_label_black",
    }

    counters = {}

    for atom in model.atom:
        elem = atom.symbol.strip()
        if not elem:
            elem = atom.name.strip()[0]

        counters[elem] = counters.get(elem, 0) + 1

        obj_name = f"{label_prefix}_{elem}"
        pseudo_name = f"L_{elem}_{counters[elem]}"

        cmd.pseudoatom(
            object=obj_name,
            name=pseudo_name,
            pos=atom.coord,
            label=elem,
        )

    for elem in counters:
        obj_name = f"{label_prefix}_{elem}"

        cmd.hide("everything", obj_name)
        cmd.show("labels", obj_name)

        _safe_set("label_position", [0, 0, 0], obj_name)
        _safe_set("label_font_id", int(label_font_id), obj_name)
        _safe_set("label_size", float(label_size), obj_name)
        _safe_set("label_connector", 0, obj_name)
        _safe_set("label_shadow_mode", 2, obj_name)

        label_color = label_color_by_element.get(elem, "poster_label_white")
        outline_color = (
            "poster_label_black"
            if label_color == "poster_label_white"
            else "poster_label_white"
        )

        _safe_set("label_color", label_color, obj_name)
        _safe_set("label_outline_color", outline_color, obj_name)


def metallic_poster_render(
    selection="all",
    metal_sel="elem Mn",
    coord_sel=None,
    coord_cutoff=2.6,
    donor_elements="N+O+S+P+H",
    background="white",
    label_core="on",
    label_size=24,
    metal_color="poster_mn_gold",
):
    """
    Apply a glossy metallic poster-style rendering to a molecular complex.

    Parameters
    ----------
    selection : str
        Whole molecule or object to style. Default: "all".

    metal_sel : str
        Selection for the central metal. Default: "elem Mn".

    coord_sel : str or None
        Optional explicit selection for coordinating atoms.
        If None, coordinating atoms are detected by distance from the metal.

    coord_cutoff : float
        Distance cutoff in Angstrom for detecting coordinating atoms.

    donor_elements : str
        Elements considered as coordinating atoms when coord_sel is None.

    background : str
        "white", "offwhite", "transparent", or "dark".

    label_core : str
        "on" or "off". If on, labels metal and coordinating atoms.

    label_size : int
        Size of centered element labels.

    metal_color : str
        PyMOL color for the metal center.
        Useful choices: poster_mn_gold, poster_mn_pink, poster_metal_gray.

    Example
    -------
    metallic_poster_render all

    metallic_poster_render all, elem Mn, None, 2.6, N+O+S+P+H, white

    metallic_poster_render all, elem Mn, None, 2.6, N+O+S+P+H, dark
    """
    define_metallic_poster_colors()

    coord_sel = _normalize_none(coord_sel)

    sel = f"({selection})"
    metal = f"({sel}) and ({metal_sel})"

    if cmd.count_atoms(metal) == 0:
        fallback_metals = "elem Mn+Fe+Co+Ni+Cu+Zn+Ru+Rh+Pd+Ag+Ir+Pt+Au"
        metal = f"({sel}) and ({fallback_metals})"

    if coord_sel is None:
        coord = (
            f"({sel}) and "
            f"(elem {donor_elements}) within {float(coord_cutoff)} of ({metal})"
        )
    else:
        coord = f"({sel}) and ({coord_sel})"

    core = f"({metal}) or ({coord})"

    bg = str(background).lower()

    if bg in ["dark", "black", "presentation", "slide"]:
        cmd.bg_color("black")
        _safe_set("ray_opaque_background", 1)
        _safe_set("ambient", 0.12)
        _safe_set("direct", 0.88)
        _safe_set("fog_start", 0.25)
    elif bg in ["offwhite", "cream", "paper"]:
        cmd.bg_color("white")
        _safe_set("ray_opaque_background", 0)
        _safe_set("ambient", 0.22)
        _safe_set("direct", 0.82)
        _safe_set("fog_start", 0.60)
    else:
        cmd.bg_color("white")
        _safe_set("ray_opaque_background", 0)
        _safe_set("ambient", 0.22)
        _safe_set("direct", 0.82)
        _safe_set("fog_start", 0.60)

    _safe_set("orthoscopic", 1)
    _safe_set("field_of_view", 35)
    _safe_set("depth_cue", 1)

    _safe_set("specular", 0.85)
    _safe_set("spec_reflect", 0.60)
    _safe_set("spec_power", 260)
    _safe_set("reflect", 0.45)
    _safe_set("shininess", 90)
    _safe_set("light_count", 8)
    _safe_set("two_sided_lighting", 1)

    _safe_set("ray_shadow", 1)
    _safe_set("ray_trace_mode", 0)
    _safe_set("ray_trace_gain", 0.08)
    _safe_set("ray_trace_disco_factor", 1)

    try:
        cmd.util.ray_shadows("light")
    except Exception:
        pass

    _safe_set("ray_shadow_decay_factor", 0.25)
    _safe_set("ray_shadow_decay_range", 2.0)

    _safe_set("ambient_occlusion_mode", 1)
    _safe_set("ambient_occlusion_scale", 18)
    _safe_set("ambient_occlusion_smooth", 12)

    _safe_set("antialias", 2)
    _safe_set("ray_trace_antialias", 2)

    cmd.hide("everything", sel)

    cmd.show("sticks", sel)
    _safe_set("stick_radius", 0.13, sel)
    _safe_set("stick_quality", 30, sel)
    _safe_set("stick_ball", 0, sel)
    _safe_set("valence", 0, sel)

    cmd.hide("sticks", f"{sel} and elem H and not ({coord})")

    cmd.show("spheres", core)
    _safe_set("sphere_quality", 4)
    _safe_set("sphere_scale", 0.62, metal)
    _safe_set("sphere_scale", 0.50, coord)
    _safe_set("sphere_scale", 0.34, f"({coord}) and elem H")

    _safe_set("stick_radius", 0.16, f"{sel} and within 2.9 of ({metal})")

    apply_metallic_poster_element_colors(
        sel,
        metal,
        metal_color=metal_color,
    )

    if str(label_core).lower() in ["on", "1", "true", "yes"]:
        _safe_set("label_font_id", 7)
        _safe_set("label_size", float(label_size))
        _safe_set("label_position", [0, 0, 0])
        _safe_set("label_connector", 0)

        _make_centered_element_labels(
            core,
            label_prefix="metallic_poster_labels",
            label_size=label_size,
            label_font_id=7,
        )

    cmd.zoom(selection, buffer=2.0)
    cmd.orient(selection)
    cmd.rebuild()


def metallic_poster_png(
    filename="metallic_poster.png", width=2400, height=1800, dpi=300
):
    """Ray-trace and save a high-resolution PNG."""
    cmd.ray(int(width), int(height))
    cmd.png(filename, dpi=int(dpi))


cmd.extend("metallic_poster_render", metallic_poster_render)
cmd.extend("metallic_poster_png", metallic_poster_png)
