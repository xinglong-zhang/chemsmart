"""
Comic metallic ball-and-stick PyMOL style with ray-traced outlines and labels.

ChemSmart applies this template when ``visualize -s comic`` is used.
Alias: ``-s hybrid``.

In PyMOL directly::

    run comic_style.py
    render_comic_metallic_labeled_final all
    render_comic_metallic_labeled_final all, white
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


def _metal_element_label(metal_selection):
    """Return the element symbol of the first atom in a metal selection."""
    model = cmd.get_model(metal_selection)
    if not model.atom:
        return "Mn"
    symbol = model.atom[0].symbol.strip()
    return symbol or "Mn"


def render_comic_metallic_labeled_final(selection="all", background="white"):
    """
    Apply comic metallic ball-and-stick rendering with black outlines and labels.

    Parameters
    ----------
    selection : str
        Molecule or object to style.
    background : str
        ``white`` (default) or ``dark`` for a black slide background.
    """
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


cmd.extend(
    "render_comic_metallic_labeled_final", render_comic_metallic_labeled_final
)
cmd.extend("comic_render", comic_render)
cmd.extend("comic_png", comic_png)
