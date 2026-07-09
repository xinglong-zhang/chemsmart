"""
PyMOL visualization styles for publication, cover, and presentation figures.

ChemSmart applies this template for ``visualize -s`` choices including
``glossy``, ``comic``, ``soft-cartoon``, ``editorial-minimal``, and
``black-gold-cover``.
In PyMOL directly::

    run scientific_styles.py
    metallic_poster_render all
    render_comic_metallic_labeled_final all
    render_soft_cartoon all
    render_editorial_minimal all
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
    """Set a PyMOL setting safely across different PyMOL versions."""
    try:
        if selection is None:
            cmd.set(setting, value)
        else:
            cmd.set(setting, value, selection)
    except Exception:
        pass


def _safe_ray_shadows(mode="light"):
    try:
        cmd.util.ray_shadows(mode)
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
    """Create centered pseudoatom labels at selected atom coordinates."""
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
    """Apply a glossy metallic poster-style rendering to a molecular complex."""
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
    _safe_ray_shadows("light")

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


def _metal_element_label(metal_selection):
    """Return the element symbol of the first atom in a metal selection."""
    model = cmd.get_model(metal_selection)
    if not model.atom:
        return "Mn"
    symbol = model.atom[0].symbol.strip()
    return symbol or "Mn"


def render_comic_metallic_labeled_final(selection="all", background="white"):
    """Apply comic metallic ball-and-stick rendering with black outlines."""
    sel = f"({selection})"
    metal_name = "comic_metal_atom"

    cmd.delete(metal_name)

    cmd.hide("everything", sel)
    cmd.show("sticks", sel)
    cmd.show("spheres", sel)

    _safe_set("stick_radius", 0.14, sel)
    _safe_set("sphere_scale", 0.25, f"{sel} and elem C+N+O+S+P")
    _safe_set("sphere_scale", 0.15, f"{sel} and elem H")

    cmd.select(
        metal_name,
        f"{sel} and (elem Mn or elem Fe or elem Co or elem Ni or elem Ru or elem Rh)",
    )
    if cmd.count_atoms(metal_name) > 0:
        _safe_set("sphere_scale", 0.42, metal_name)

    apply_metallic_poster_element_colors(
        sel,
        metal_name,
        metal_color=DEFAULT_METAL_COLOR,
    )

    if cmd.count_atoms(metal_name) > 0:
        cmd.bond(metal_name, f"{sel} and elem N")
        cmd.bond(metal_name, f"{sel} and elem P")
        cmd.bond(metal_name, f"{sel} and elem S")

    _safe_set("specular", 0.85)
    _safe_set("spec_power", 300)
    _safe_set("spec_reflect", 0.70)
    _safe_set("ambient", 0.25)
    _safe_set("direct", 0.75)
    _safe_set("shininess", 90)

    _safe_set("ray_trace_mode", 1)
    _safe_set("ray_trace_color", "black")
    _safe_set("ray_trace_gain", 0.6)

    cmd.label(sel, '""')
    if cmd.count_atoms(metal_name) > 0:
        metal_label = _metal_element_label(metal_name)
        cmd.label(metal_name, f'"{metal_label}"')
    cmd.label(f"{sel} and elem N", '"N"')
    cmd.label(f"{sel} and elem P", '"P"')
    cmd.label(f"{sel} and elem S", '"S"')

    _safe_set("label_position", [0.0, 0.0, 2.0])
    _safe_set("label_shadow_mode", 0)
    _safe_set("ray_label_specular", 0)
    _safe_set("label_font_id", 7)
    _safe_set("label_color", "white")
    _safe_set("label_size", 28)

    bg = str(background).lower()
    if bg in ["dark", "black", "presentation", "slide"]:
        cmd.bg_color("black")
        _safe_set("ray_opaque_background", 1)
    else:
        cmd.bg_color("white")
        _safe_set("ray_opaque_background", 1)

    _safe_set("orthoscopic", 1)
    _safe_set("antialias", 2)

    cmd.zoom(selection, buffer=2.0)
    cmd.orient(selection)
    cmd.refresh()
    print(
        "Comic metallic style applied. Run 'ray 1200, 1200' to generate the image."
    )


def comic_render(selection="all", background="white", *_args, **_kwargs):
    """ChemSmart entry point for the comic style."""
    render_comic_metallic_labeled_final(
        selection=selection, background=background
    )


def comic_png(filename="comic.png", width=2400, height=1800, dpi=300):
    """Ray-trace and save a high-resolution PNG."""
    cmd.ray(int(width), int(height))
    cmd.png(filename, dpi=int(dpi))


def render_soft_cartoon(selection="all", background="white"):
    """Apply soft cartoon ball-and-stick rendering with premium cover colors."""
    sel = f"({selection})"
    metal_name = "soft_cartoon_metal_atom"
    coord_core_name = "soft_cartoon_coord_core"

    cmd.delete(metal_name)
    cmd.delete(coord_core_name)

    cmd.hide("everything", sel)
    cmd.show("sticks", sel)
    cmd.show("spheres", sel)

    _safe_set("stick_radius", 0.12, sel)
    _safe_set("sphere_scale", 0.24, f"{sel} and elem C+N+O+S+P+F+Cl+Br+I")
    _safe_set("sphere_scale", 0.14, f"{sel} and elem H")

    cmd.select(
        metal_name,
        f"{sel} and (elem Mn or elem Fe or elem Co or elem Ni or elem Ru "
        f"or elem Rh or elem Pd or elem Pt or elem Ir)",
    )
    if cmd.count_atoms(metal_name) > 0:
        _safe_set("sphere_scale", 0.48, metal_name)

    cmd.select(
        coord_core_name,
        f"{metal_name} or ((elem N+O+S+P+H) within 2.6 of {metal_name})",
    )
    _safe_set("sphere_scale", 0.36, f"{coord_core_name} and not {metal_name}")
    _safe_set(
        "stick_radius",
        0.14,
        f"byres ({coord_core_name} around 2.8)",
    )

    cmd.set_color("cover_carbon", [0.66, 0.66, 0.64])
    cmd.set_color("cover_hydrogen", [0.96, 0.96, 0.93])
    cmd.set_color("cover_nitrogen", [0.16, 0.26, 0.90])
    cmd.set_color("cover_oxygen", [0.92, 0.05, 0.04])
    cmd.set_color("cover_sulfur", [0.95, 0.70, 0.14])
    cmd.set_color("cover_mn_gold", [0.95, 0.66, 0.24])

    cmd.color("cover_carbon", f"{sel} and elem C")
    cmd.color("cover_hydrogen", f"{sel} and elem H")
    cmd.color("cover_nitrogen", f"{sel} and elem N")
    cmd.color("cover_oxygen", f"{sel} and elem O")
    cmd.color("cover_sulfur", f"{sel} and elem S")
    cmd.color("orange", f"{sel} and elem P")
    if cmd.count_atoms(metal_name) > 0:
        cmd.color("cover_mn_gold", metal_name)

    _safe_set("specular", 0.78)
    _safe_set("spec_power", 280)
    _safe_set("spec_reflect", 0.60)
    _safe_set("ambient", 0.22)
    _safe_set("direct", 0.78)
    _safe_set("reflect", 0.38)
    _safe_set("shininess", 85)

    _safe_set("ray_trace_mode", 0)
    _safe_set("ray_shadow", "on")
    _safe_ray_shadows("light")
    _safe_set("ray_trace_gain", 0.10)

    cmd.label(sel, '""')

    bg = str(background).lower()
    if bg in ["dark", "black", "presentation", "slide"]:
        cmd.bg_color("black")
        _safe_set("ray_opaque_background", 1)
    else:
        cmd.bg_color("white")
        _safe_set("ray_opaque_background", "off")

    _safe_set("orthoscopic", 0)
    _safe_set("field_of_view", 35)
    _safe_set("antialias", 2)
    _safe_set("depth_cue", 1)
    _safe_set("fog_start", 0.35)

    cmd.zoom(selection, buffer=2.0)
    cmd.orient(selection)
    cmd.refresh()
    print("Soft cartoon style applied.")


def soft_cartoon_render(
    selection="all", background="white", *_args, **_kwargs
):
    """ChemSmart entry point for the soft cartoon style."""
    render_soft_cartoon(selection=selection, background=background)


def soft_cartoon_png(
    filename="soft_cartoon.png", width=2400, height=1800, dpi=300
):
    """Ray-trace and save a high-resolution PNG."""
    cmd.ray(int(width), int(height))
    cmd.png(filename, dpi=int(dpi))


def _common_select_core(selection="all", cutoff=2.6):
    sel = f"({selection})"

    cmd.delete("metal_atom")
    cmd.delete("coord_core")
    cmd.delete("near_core")

    cmd.select(
        "metal_atom",
        f"{sel} and (elem Mn or elem Fe or elem Co or elem Ni or elem Cu or elem Zn "
        f"or elem Ru or elem Rh or elem Pd or elem Ag or elem Ir or elem Pt or elem Au "
        f"or elem Mg or elem Al or elem Ti or elem Zr)",
    )

    cmd.select(
        "coord_core",
        f"metal_atom or ({sel} and (elem N+O+S+P+H+Cl+Br+I) within {cutoff} of metal_atom)",
    )

    cmd.select(
        "near_core",
        f"{sel} and byres (coord_core around 3.0)",
    )


def _define_scientific_colors():
    colors = {
        "sci_C_gray": [0.58, 0.58, 0.56],
        "sci_C_dark": [0.10, 0.10, 0.10],
        "sci_C_ivory": [0.82, 0.79, 0.70],
        "sci_H_white": [0.96, 0.96, 0.94],
        "sci_N_blue": [0.12, 0.26, 0.92],
        "sci_O_red": [0.92, 0.05, 0.04],
        "sci_S_yellow": [0.95, 0.72, 0.16],
        "sci_P_orange": [1.00, 0.45, 0.10],
        "sci_halogen": [0.20, 0.78, 0.30],
        "metal_gold": [0.95, 0.63, 0.20],
        "metal_rose": [0.90, 0.45, 0.75],
        "metal_silver": [0.78, 0.80, 0.84],
        "deep_black": [0.01, 0.01, 0.015],
        "neon_cyan": [0.00, 0.85, 1.00],
        "neon_magenta": [1.00, 0.05, 0.75],
        "neon_green": [0.15, 1.00, 0.35],
        "clay_carbon": [0.70, 0.64, 0.56],
        "clay_blue": [0.32, 0.45, 0.78],
        "clay_red": [0.78, 0.30, 0.25],
        "surface_sky": [0.55, 0.78, 1.00],
        "surface_warm": [1.00, 0.76, 0.45],
    }

    for name, rgb in colors.items():
        try:
            cmd.set_color(name, rgb)
        except Exception:
            pass


def _color_by_element(carbon="sci_C_gray", metal="metal_gold"):
    cmd.color(carbon, "elem C")
    cmd.color("sci_H_white", "elem H")
    cmd.color("sci_N_blue", "elem N")
    cmd.color("sci_O_red", "elem O")
    cmd.color("sci_S_yellow", "elem S")
    cmd.color("sci_P_orange", "elem P")
    cmd.color("sci_halogen", "elem F+Cl+Br+I")
    cmd.color(metal, "metal_atom")


def _base_quality():
    _safe_set("antialias", 2)
    _safe_set("ray_trace_antialias", 2)
    _safe_set("sphere_quality", 3)
    _safe_set("stick_quality", 30)
    _safe_set("two_sided_lighting", 1)
    _safe_set("use_shaders", 1)
    _safe_set("depth_cue", 1)
    _safe_set("label_connector", 0)


def _finish_style(selection):
    cmd.zoom(selection, buffer=2.0)
    cmd.orient(selection)
    cmd.refresh()


def render_editorial_minimal(selection="all"):
    """Editorial minimal white style for main-text mechanistic figures."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", "coord_core")

    cmd.set("stick_radius", 0.095)
    cmd.set("stick_radius", 0.135, "near_core")
    cmd.set("sphere_scale", 0.24, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", 0.12, "elem H")
    cmd.set("sphere_scale", 0.45, "metal_atom")
    cmd.set("sphere_scale", 0.32, "coord_core and not metal_atom")

    _color_by_element(carbon="sci_C_gray", metal="metal_gold")

    cmd.set("specular", 0.25)
    cmd.set("spec_reflect", 0.05)
    cmd.set("spec_power", 80)
    cmd.set("ambient", 0.42)
    cmd.set("direct", 0.62)
    cmd.set("reflect", 0.05)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 1)
    cmd.set("field_of_view", 25)
    cmd.set("fog_start", 0.65)
    cmd.set("ray_shadow", "on")
    _safe_ray_shadows("light")
    _base_quality()

    cmd.label("all", '""')
    _finish_style(selection)
    print("Editorial minimal white style applied.")


def render_black_gold_cover(selection="all"):
    """Black-gold journal-cover style for high-impact centerpiece figures."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", "coord_core")

    cmd.set("stick_radius", 0.12)
    cmd.set("stick_radius", 0.17, "near_core")
    cmd.set("sphere_scale", 0.26, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", 0.13, "elem H")
    cmd.set("sphere_scale", 0.58, "metal_atom")
    cmd.set("sphere_scale", 0.40, "coord_core and not metal_atom")

    cmd.color("sci_C_ivory", "elem C")
    cmd.color("sci_H_white", "elem H")
    cmd.color("sci_N_blue", "elem N")
    cmd.color("sci_O_red", "elem O")
    cmd.color("sci_S_yellow", "elem S")
    cmd.color("sci_P_orange", "elem P")
    cmd.color("sci_halogen", "elem F+Cl+Br+I")
    cmd.color("metal_gold", "metal_atom")

    cmd.set("specular", 0.90)
    cmd.set("spec_reflect", 0.72)
    cmd.set("spec_power", 320)
    cmd.set("ambient", 0.12)
    cmd.set("direct", 0.88)
    cmd.set("reflect", 0.52)
    cmd.set("shininess", 95)

    cmd.bg_color("black")
    cmd.set("ray_opaque_background", "on")
    cmd.set("orthoscopic", 0)
    cmd.set("field_of_view", 38)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.22)
    cmd.set("ray_shadow", "on")
    cmd.set("ray_trace_gain", 0.18)
    _safe_ray_shadows("light")
    _base_quality()

    cmd.label("all", '""')
    _finish_style(selection)
    print("Black-gold cover style applied.")


def render_neon_coordination_core(selection="all"):
    """Neon coordination-core style for reactive centers and catalytic pockets."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", "coord_core")

    cmd.set("stick_radius", 0.085)
    cmd.set("stick_radius", 0.18, "near_core")
    cmd.set("sphere_scale", 0.20, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", 0.11, "elem H")
    cmd.set("sphere_scale", 0.62, "metal_atom")
    cmd.set("sphere_scale", 0.42, "coord_core and not metal_atom")

    cmd.color("sci_C_dark", "elem C")
    cmd.color("sci_H_white", "elem H")
    cmd.color("neon_cyan", "elem N")
    cmd.color("neon_magenta", "elem O")
    cmd.color("neon_green", "elem S")
    cmd.color("sci_P_orange", "elem P")
    cmd.color("neon_green", "elem F+Cl+Br+I")
    cmd.color("metal_rose", "metal_atom")

    cmd.set("specular", 0.78)
    cmd.set("spec_reflect", 0.45)
    cmd.set("spec_power", 180)
    cmd.set("ambient", 0.08)
    cmd.set("direct", 0.95)
    cmd.set("reflect", 0.28)

    cmd.bg_color("black")
    cmd.set("ray_opaque_background", "on")
    cmd.set("orthoscopic", 0)
    cmd.set("field_of_view", 42)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.18)
    cmd.set("ray_shadow", "on")
    cmd.set("ray_trace_gain", 0.28)
    _safe_ray_shadows("light")
    _base_quality()

    cmd.label("all", '""')
    _finish_style(selection)
    print("Neon coordination-core style applied.")


def render_matte_clay(selection="all"):
    """Matte clay model style for soft graphical abstracts."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", "coord_core")

    cmd.set("stick_radius", 0.14)
    cmd.set("stick_radius", 0.18, "near_core")
    cmd.set("sphere_scale", 0.30, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", 0.15, "elem H")
    cmd.set("sphere_scale", 0.55, "metal_atom")
    cmd.set("sphere_scale", 0.40, "coord_core and not metal_atom")

    cmd.color("clay_carbon", "elem C")
    cmd.color("sci_H_white", "elem H")
    cmd.color("clay_blue", "elem N")
    cmd.color("clay_red", "elem O")
    cmd.color("sci_S_yellow", "elem S")
    cmd.color("sci_P_orange", "elem P")
    cmd.color("sci_halogen", "elem F+Cl+Br+I")
    cmd.color("metal_silver", "metal_atom")

    cmd.set("specular", 0.08)
    cmd.set("spec_reflect", 0.00)
    cmd.set("spec_power", 30)
    cmd.set("ambient", 0.55)
    cmd.set("direct", 0.48)
    cmd.set("reflect", 0.03)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 1)
    cmd.set("field_of_view", 28)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.50)
    cmd.set("ray_shadow", "on")
    _safe_ray_shadows("soft")
    _base_quality()

    cmd.label("all", '""')
    _finish_style(selection)
    print("Matte clay style applied.")


def render_xray_wire(selection="all"):
    """X-ray crystallography wire style for SI structure verification."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)
    cmd.show("lines", selection)
    cmd.show("sticks", "near_core")
    cmd.show("spheres", "metal_atom")

    cmd.set("line_width", 2.0)
    cmd.set("stick_radius", 0.06)
    cmd.set("stick_radius", 0.11, "near_core")
    cmd.set("sphere_scale", 0.38, "metal_atom")

    cmd.color("sci_C_dark", selection)
    cmd.color("sci_N_blue", "coord_core and elem N")
    cmd.color("sci_O_red", "coord_core and elem O")
    cmd.color("sci_S_yellow", "coord_core and elem S")
    cmd.color("metal_gold", "metal_atom")

    cmd.set("specular", 0.05)
    cmd.set("spec_reflect", 0.00)
    cmd.set("ambient", 0.70)
    cmd.set("direct", 0.35)
    cmd.set("reflect", 0.00)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 1)
    cmd.set("field_of_view", 18)
    cmd.set("ray_shadow", "off")
    cmd.set("depth_cue", 0)
    _base_quality()

    cmd.label("all", '""')
    _finish_style(selection)
    print("X-ray wire style applied.")


def render_steric_surface(selection="all"):
    """Transparent steric surface style for catalyst pockets."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)

    cmd.show("sticks", selection)
    cmd.show("spheres", "coord_core")

    cmd.set("stick_radius", 0.095)
    cmd.set("stick_radius", 0.14, "near_core")
    cmd.set("sphere_scale", 0.24, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", 0.12, "elem H")
    cmd.set("sphere_scale", 0.52, "metal_atom")
    cmd.set("sphere_scale", 0.34, "coord_core and not metal_atom")

    _color_by_element(carbon="sci_C_gray", metal="metal_gold")

    cmd.show("surface", selection)
    cmd.color("surface_sky", selection)
    cmd.set("transparency", 0.68, selection)
    cmd.set("surface_quality", 1)
    cmd.set("surface_solvent", 1)

    _color_by_element(carbon="sci_C_gray", metal="metal_gold")

    cmd.set("specular", 0.45)
    cmd.set("spec_reflect", 0.22)
    cmd.set("spec_power", 120)
    cmd.set("ambient", 0.32)
    cmd.set("direct", 0.74)
    cmd.set("reflect", 0.20)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 1)
    cmd.set("field_of_view", 30)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.45)
    cmd.set("ray_shadow", "on")
    _safe_ray_shadows("light")
    _base_quality()

    cmd.label("all", '""')
    _finish_style(selection)
    print("Transparent steric surface style applied.")


def render_quasi_chemdraw_bold(selection="all"):
    """Quasi-ChemDraw bold 3D style with formula-like clarity."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", "metal_atom or coord_core")

    cmd.hide("sticks", f"({selection}) and elem H and not coord_core")
    cmd.hide("spheres", f"({selection}) and elem H and not coord_core")

    cmd.set("stick_radius", 0.18)
    cmd.set("stick_radius", 0.24, "near_core")
    cmd.set("sphere_scale", 0.23, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", 0.13, "elem H")
    cmd.set("sphere_scale", 0.50, "metal_atom")
    cmd.set("sphere_scale", 0.34, "coord_core and not metal_atom")

    cmd.color("sci_C_dark", "elem C")
    cmd.color("sci_H_white", "elem H")
    cmd.color("sci_N_blue", "elem N")
    cmd.color("sci_O_red", "elem O")
    cmd.color("sci_S_yellow", "elem S")
    cmd.color("sci_P_orange", "elem P")
    cmd.color("sci_halogen", "elem F+Cl+Br+I")
    cmd.color("metal_gold", "metal_atom")

    cmd.set("specular", 0.18)
    cmd.set("spec_reflect", 0.02)
    cmd.set("spec_power", 70)
    cmd.set("ambient", 0.62)
    cmd.set("direct", 0.42)
    cmd.set("reflect", 0.02)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 1)
    cmd.set("field_of_view", 12)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadow", "off")
    _base_quality()

    cmd.label("all", '""')
    _finish_style(selection)
    print("Quasi-ChemDraw bold 3D style applied.")


def render_labeled_coordination_core(selection="all"):
    """Coordination-core labeled style with explicit element labels."""
    _define_scientific_colors()
    _common_select_core(selection)

    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", "coord_core")

    cmd.set("stick_radius", 0.11)
    cmd.set("stick_radius", 0.16, "near_core")
    cmd.set("sphere_scale", 0.25, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", 0.13, "elem H")
    cmd.set("sphere_scale", 0.56, "metal_atom")
    cmd.set("sphere_scale", 0.38, "coord_core and not metal_atom")

    _color_by_element(carbon="sci_C_gray", metal="metal_gold")

    cmd.set("specular", 0.55)
    cmd.set("spec_reflect", 0.28)
    cmd.set("spec_power", 180)
    cmd.set("ambient", 0.30)
    cmd.set("direct", 0.70)
    cmd.set("reflect", 0.18)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 1)
    cmd.set("field_of_view", 28)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.45)
    cmd.set("ray_shadow", "on")
    _safe_ray_shadows("light")
    _base_quality()

    cmd.label("all", '""')
    cmd.label("coord_core", "elem")
    cmd.set("label_size", 24)
    cmd.set("label_font_id", 7)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_position", [0, 0, 0])
    cmd.set("label_connector", 0)

    _finish_style(selection)
    print("Labeled coordination-core style applied.")


cmd.extend("metallic_poster_render", metallic_poster_render)
cmd.extend("metallic_poster_png", metallic_poster_png)
cmd.extend(
    "render_comic_metallic_labeled_final", render_comic_metallic_labeled_final
)
cmd.extend("comic_render", comic_render)
cmd.extend("comic_png", comic_png)
cmd.extend("render_soft_cartoon", render_soft_cartoon)
cmd.extend("soft_cartoon_render", soft_cartoon_render)
cmd.extend("soft_cartoon_png", soft_cartoon_png)
cmd.extend("render_editorial_minimal", render_editorial_minimal)
cmd.extend("render_black_gold_cover", render_black_gold_cover)
cmd.extend("render_neon_coordination_core", render_neon_coordination_core)
cmd.extend("render_matte_clay", render_matte_clay)
cmd.extend("render_xray_wire", render_xray_wire)
cmd.extend("render_steric_surface", render_steric_surface)
cmd.extend("render_quasi_chemdraw_bold", render_quasi_chemdraw_bold)
cmd.extend(
    "render_labeled_coordination_core", render_labeled_coordination_core
)
