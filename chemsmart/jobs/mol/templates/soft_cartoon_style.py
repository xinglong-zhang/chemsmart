"""
Soft cartoon PyMOL style with premium colors and light metallic shading.

ChemSmart applies this template when ``visualize -s soft-cartoon`` is used.
In PyMOL directly::

    run soft_cartoon_style.py
    render_soft_cartoon all
    render_soft_cartoon all, white
"""

from pymol import cmd


def _safe_set(setting, value, selection=None):
    """Set a PyMOL setting safely across different PyMOL versions."""
    try:
        if selection is None:
            cmd.set(setting, value)
        else:
            cmd.set(setting, value, selection)
    except Exception:
        pass


def render_soft_cartoon(selection="all", background="white"):
    """
    Apply soft cartoon ball-and-stick rendering with premium cover colors.

    Parameters
    ----------
    selection : str
        Molecule or object to style.
    background : str
        ``white`` (default) or ``dark`` for a black slide background.
    """
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
    cmd.util.ray_shadows("light")
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


cmd.extend("render_soft_cartoon", render_soft_cartoon)
cmd.extend("soft_cartoon_render", soft_cartoon_render)
cmd.extend("soft_cartoon_png", soft_cartoon_png)
