"""
Shared element color palette for glossy and comic_ballstick PyMOL styles.

Both ``glossy_metal_style.py`` and ``comic_ballstick_style.py`` import this
module so every element renders with identical colors in either style.
"""

from pymol import cmd

DEFAULT_METAL_COLOR = "poster_mn_gold"

METALLIC_POSTER_COLOR_RGB = {
    "poster_carbon": [0.62, 0.62, 0.62],
    "poster_hydrogen": [0.96, 0.96, 0.96],
    "poster_nitrogen": [0.12, 0.22, 0.95],
    "poster_oxygen": [0.95, 0.06, 0.04],
    "poster_sulfur": [0.95, 0.72, 0.10],
    "poster_phosphorus": [1.00, 0.48, 0.08],
    "poster_halogen": [0.20, 0.82, 0.32],
    "poster_mn_gold": [0.92, 0.60, 0.22],
    "poster_mn_pink": [0.88, 0.48, 0.78],
    "poster_metal_gray": [0.78, 0.78, 0.84],
    "poster_label_white": [1.00, 1.00, 1.00],
    "poster_label_black": [0.02, 0.02, 0.02],
}

OTHER_METALS_SELECTION = "elem Fe+Co+Ni+Cu+Zn+Ru+Rh+Pd+Ag+Ir+Pt+Au"


def define_metallic_poster_colors():
    """Register the shared poster element colors in the current PyMOL session."""
    for name, rgb in METALLIC_POSTER_COLOR_RGB.items():
        try:
            cmd.set_color(name, rgb)
        except Exception:
            pass


def apply_metallic_poster_element_colors(
    sel,
    metal,
    metal_color=DEFAULT_METAL_COLOR,
):
    """
    Apply the shared element color scheme to a selection.

    Parameters
    ----------
    sel : str
        PyMOL selection for the full molecule/object, e.g. ``(mol1)``.
    metal : str
        PyMOL selection for the primary metal center, e.g. ``comic_metal_atom``
        or ``(mol1) and (elem Mn)``.
    metal_color : str
        PyMOL color name for the primary metal (default: ``poster_mn_gold``).
    """
    define_metallic_poster_colors()

    cmd.color("poster_carbon", f"{sel} and elem C")
    cmd.color("poster_hydrogen", f"{sel} and elem H")
    cmd.color("poster_nitrogen", f"{sel} and elem N")
    cmd.color("poster_oxygen", f"{sel} and elem O")
    cmd.color("poster_sulfur", f"{sel} and elem S")
    cmd.color("poster_phosphorus", f"{sel} and elem P")
    cmd.color("poster_halogen", f"{sel} and elem F+Cl+Br+I")

    if cmd.count_atoms(metal) > 0:
        cmd.color(metal_color, metal)

    cmd.color(
        "poster_metal_gray",
        f"{sel} and {OTHER_METALS_SELECTION} and not ({metal})",
    )
