"""
PyMOL visualization styles for publication, cover, and presentation figures.

CHEMSMART applies this template for ``visualize -s`` choices including
``glossy``, ``comic``, ``soft-cartoon``, ``editorial-minimal``,
``soft-ceramic``, and other scientific styles.
In PyMOL directly::

    run scientific_styles.py
    metallic_poster_render all
    render_comic_metallic_labeled_final all
    render_soft_cartoon all
    render_editorial_minimal all
    render_soft_ceramic all
    render_matte_clay all
"""

import numpy as np
from pymol import cmd

from chemsmart.utils.geometry import get_coordinating_atoms

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
    """Register the poster element colors in the current PyMOL session."""
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
    """Apply the poster element color scheme to a selection."""
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


def _set_transparent_background():
    """Configure ray-traced PNG export with a transparent background."""
    cmd.bg_color("white")
    _safe_set("ray_opaque_background", 0)


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
    highlight_bonds="",
):
    """Apply a glossy metallic poster-style rendering to a molecular complex."""
    define_metallic_poster_colors()

    coord_sel = _normalize_none(coord_sel)
    highlight_bonds = _normalize_none(highlight_bonds)
    if coord_sel is None and highlight_bonds:
        coord_sel = _highlight_bonds_to_coord_sel(highlight_bonds)

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

    _set_transparent_background()
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

    if highlight_bonds:
        _apply_coordination_highlight_bonds(highlight_bonds, sel, metal)

    cmd.zoom(selection, buffer=2.0)
    cmd.orient(selection)
    cmd.rebuild()


def _metal_element_label(metal_selection):
    """Return the element symbol of the first atom in a metal selection."""
    model = cmd.get_model(metal_selection)
    if not model.atom:
        return "Mn"
    symbol = model.atom[0].symbol.strip()
    return symbol or "Mn"


def _parse_highlight_bond_pairs(highlight_bonds):
    """Parse CHEMSMART-encoded bond pairs such as ``1-8+1-15+1-36``."""
    highlight_bonds = _normalize_none(highlight_bonds)
    if not highlight_bonds:
        return []

    pairs = []
    for token in str(highlight_bonds).split("+"):
        token = token.strip()
        if not token:
            continue
        if "," in token:
            parts = [part.strip() for part in token.split(",")]
        elif "-" in token:
            parts = [part.strip() for part in token.split("-")]
        else:
            continue
        if len(parts) != 2:
            continue
        pairs.append((int(parts[0]), int(parts[1])))
    return pairs


def _highlight_bonds_to_coord_sel(highlight_bonds):
    """Return a PyMOL selection for all atoms referenced in bond pairs."""
    pairs = _parse_highlight_bond_pairs(highlight_bonds)
    if not pairs:
        return None
    atom_ids = sorted({atom_id for pair in pairs for atom_id in pair})
    return " or ".join(f"id {atom_id}" for atom_id in atom_ids)


def _apply_coordination_highlight_bonds(highlight_bonds, sel, metal_name):
    """Style ``-c`` bond pairs like coordination-core metal-ligand bonds."""
    pairs = _parse_highlight_bond_pairs(highlight_bonds)
    if not pairs:
        return

    for atom_a, atom_b in pairs:
        cmd.bond(f"id {atom_a}", f"id {atom_b}")

    highlight_atoms = " or ".join(
        f"(id {atom_a} or id {atom_b})" for atom_a, atom_b in pairs
    )
    _safe_set("stick_radius", 0.18, highlight_atoms)
    _safe_set(
        "sphere_scale",
        0.34,
        f"({highlight_atoms}) and not ({metal_name})",
    )
    if cmd.count_atoms(metal_name) > 0:
        _safe_set(
            "sphere_scale",
            0.42,
            f"({highlight_atoms}) and ({metal_name})",
        )


def _apply_comic_highlight_bonds(highlight_bonds, sel, metal_name):
    """Backward-compatible alias for comic highlight styling."""
    _apply_coordination_highlight_bonds(highlight_bonds, sel, metal_name)


def render_comic_metallic_labeled_final(
    selection="all", highlight_bonds="", background=None
):
    """Apply comic metallic ball-and-stick rendering with black outlines."""
    sel = f"({selection})"
    metal_name = "comic_metal_atom"

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

    _apply_comic_highlight_bonds(highlight_bonds, sel, metal_name)

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

    _set_transparent_background()

    _safe_set("orthoscopic", 1)
    _safe_set("antialias", 2)

    cmd.zoom(selection, buffer=2.0)
    cmd.orient(selection)
    cmd.refresh()
    print(
        "Comic metallic style applied. Run 'ray 1200, 1200' to generate the image."
    )


def comic_render(
    selection="all", highlight_bonds="", background=None, *_args, **_kwargs
):
    """CHEMSMART entry point for the comic style."""
    render_comic_metallic_labeled_final(
        selection=selection,
        highlight_bonds=highlight_bonds,
    )


def render_soft_cartoon(selection="all", background=None, highlight_bonds=""):
    """Apply soft cartoon ball-and-stick rendering with premium cover colors."""
    sel = f"({selection})"
    metal_name = "soft_cartoon_metal_atom"
    coord_core_name = "soft_cartoon_coord_core"
    highlight_bonds = _normalize_none(highlight_bonds)

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

    if highlight_bonds:
        coord_sel = _highlight_bonds_to_coord_sel(highlight_bonds)
        cmd.select(coord_core_name, f"{sel} and ({coord_sel})")
    else:
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

    if highlight_bonds:
        _apply_coordination_highlight_bonds(highlight_bonds, sel, metal_name)

    _set_transparent_background()

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
    selection="all", background=None, highlight_bonds="", *_args, **_kwargs
):
    """CHEMSMART entry point for the soft cartoon style."""
    render_soft_cartoon(
        selection=selection,
        background=background,
        highlight_bonds=highlight_bonds,
    )


_METAL_ELEMENTS = "elem Mn+Fe+Co+Ni+Cu+Zn+Ru+Rh+Pd+Ag+Ir+Pt+Au+Mg+Al+Ti+Zr"


def _parse_metal_symbols(metal):
    """Parse a PyMOL metal selection string into a set of element symbols."""
    text = (metal or _METAL_ELEMENTS).replace("elem", " ")
    return {
        token.strip()
        for token in text.replace("+", " ").split()
        if token.strip()
    }


def _pymol_index_selection(sel, pymol_indices):
    """Build a PyMOL selection expression from atom indices."""
    indices = sorted({int(index) for index in pymol_indices})
    if not indices:
        return "none"
    return "%s and index %s" % (
        sel,
        "+".join(str(index) for index in indices),
    )


def _get_coordinating_atoms(
    selection="all",
    prefix="coord",
    metal=None,
    include_nh_h=False,
):
    """Build named selections via covalent-radius-ratio coordination spheres.

    Extracts ChemPy coordinates, calls
    ``chemsmart.utils.geometry.get_coordinating_atoms`` for each metal center,
    and partitions the resulting primary/secondary spheres into role-based
    named selections. Applies no colors or rendering.
    """
    sel = f"({selection})"
    names = {
        "metal": "%s_metal" % prefix,
        "donor_s": "%s_donor_s" % prefix,
        "donor_n": "%s_donor_n" % prefix,
        "donor_p": "%s_donor_p" % prefix,
        "donors": "%s_donors" % prefix,
        "co_c": "%s_co_c" % prefix,
        "co_o": "%s_co_o" % prefix,
        "hydride": "%s_hydride" % prefix,
        "nh_h": "%s_nh_h" % prefix,
        "important_h": "%s_important_h" % prefix,
        "coordination_core": "%s_coordination_core" % prefix,
    }

    model = cmd.get_model(sel)
    if len(model.atom) == 0:
        for name in names.values():
            cmd.select(name, "none")
        return names

    elements = []
    pymol_indices = []
    coords = []
    for atom in model.atom:
        elements.append(str(getattr(atom, "symbol", "")).strip())
        pymol_indices.append(int(atom.index))
        coords.append([float(value) for value in atom.coord])
    coordinates = np.asarray(coords, dtype=float)

    metal_symbols = _parse_metal_symbols(metal)
    metal_local = [
        idx for idx, element in enumerate(elements) if element in metal_symbols
    ]

    primary_local = set()
    secondary_local = set()
    for metal_idx in metal_local:
        primary, secondary = get_coordinating_atoms(
            metal_idx, elements, coordinates
        )
        primary_local.update(primary)
        secondary_local.update(secondary)

    def role_indices(predicate, local_indices):
        return [
            pymol_indices[idx]
            for idx in sorted(local_indices)
            if predicate(elements[idx])
        ]

    metal_pymol = [pymol_indices[idx] for idx in metal_local]
    hydride_pymol = role_indices(lambda el: el == "H", primary_local)
    nh_h_pymol = (
        role_indices(lambda el: el == "H", secondary_local)
        if include_nh_h
        else []
    )
    donor_s_pymol = role_indices(lambda el: el == "S", primary_local)
    donor_n_pymol = role_indices(lambda el: el == "N", primary_local)
    donor_p_pymol = role_indices(lambda el: el == "P", primary_local)
    co_c_pymol = role_indices(lambda el: el == "C", primary_local)
    co_o_pymol = role_indices(
        lambda el: el == "O", primary_local | secondary_local
    )
    important_h_pymol = sorted(set(hydride_pymol) | set(nh_h_pymol))
    core_pymol = sorted(
        set(metal_pymol)
        | {pymol_indices[idx] for idx in primary_local}
        | set(co_o_pymol)
    )

    cmd.select(names["metal"], _pymol_index_selection(sel, metal_pymol))
    cmd.select(names["donor_s"], _pymol_index_selection(sel, donor_s_pymol))
    cmd.select(names["donor_n"], _pymol_index_selection(sel, donor_n_pymol))
    cmd.select(names["donor_p"], _pymol_index_selection(sel, donor_p_pymol))
    donor_parts = []
    if donor_s_pymol:
        donor_parts.append(names["donor_s"])
    if donor_n_pymol:
        donor_parts.append(names["donor_n"])
    if donor_p_pymol:
        donor_parts.append(names["donor_p"])
    if donor_parts:
        cmd.select(names["donors"], " or ".join(donor_parts))
    else:
        cmd.select(names["donors"], "none")
    cmd.select(names["co_c"], _pymol_index_selection(sel, co_c_pymol))
    cmd.select(names["co_o"], _pymol_index_selection(sel, co_o_pymol))
    cmd.select(names["hydride"], _pymol_index_selection(sel, hydride_pymol))
    cmd.select(names["nh_h"], _pymol_index_selection(sel, nh_h_pymol))
    cmd.select(
        names["important_h"], _pymol_index_selection(sel, important_h_pymol)
    )
    cmd.select(
        names["coordination_core"], _pymol_index_selection(sel, core_pymol)
    )
    return names


def _common_select_core(selection="all", cutoff=2.6, highlight_bonds=None):
    sel = f"({selection})"
    highlight_bonds = _normalize_none(highlight_bonds)

    if highlight_bonds:
        coord_sel = _highlight_bonds_to_coord_sel(highlight_bonds)
        cmd.select("coord_core", f"{sel} and ({coord_sel})")
        cmd.select("metal_atom", f"coord_core and ({_METAL_ELEMENTS})")
    else:
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


def _begin_scientific_style(selection, highlight_bonds=""):
    """Shared setup for editorial / scientific derived styles."""
    _define_scientific_colors()
    _common_select_core(selection, highlight_bonds=highlight_bonds or None)
    cmd.hide("everything", selection)


def _finish_scientific_style(selection, highlight_bonds="", message=""):
    """Apply optional bond highlighting and finish a scientific style."""
    highlight_bonds = _normalize_none(highlight_bonds)
    if highlight_bonds:
        _apply_coordination_highlight_bonds(
            highlight_bonds,
            f"({selection})",
            "metal_atom",
        )
    _end_scientific_style(selection, message)


def _apply_coord_sphere_scales(
    stick_radius,
    stick_radius_near,
    sphere_heavy,
    sphere_h,
    sphere_metal,
    sphere_coord,
):
    """Apply the common stick/sphere scale tiers used by scientific styles."""
    cmd.set("stick_radius", stick_radius)
    cmd.set("stick_radius", stick_radius_near, "near_core")
    cmd.set("sphere_scale", sphere_heavy, "elem C+N+O+S+P+F+Cl+Br+I")
    cmd.set("sphere_scale", sphere_h, "elem H")
    cmd.set("sphere_scale", sphere_metal, "metal_atom")
    cmd.set("sphere_scale", sphere_coord, "coord_core and not metal_atom")


def _apply_lighting(
    specular,
    spec_reflect,
    ambient,
    direct,
    reflect,
    spec_power=None,
    shininess=None,
):
    cmd.set("specular", specular)
    cmd.set("spec_reflect", spec_reflect)
    if spec_power is not None:
        cmd.set("spec_power", spec_power)
    cmd.set("ambient", ambient)
    cmd.set("direct", direct)
    cmd.set("reflect", reflect)
    if shininess is not None:
        cmd.set("shininess", shininess)


def _apply_view(
    orthoscopic,
    field_of_view,
    ray_shadow="on",
    depth_cue=None,
    fog_start=None,
    ray_trace_gain=None,
    ray_shadows_mode="light",
):
    """Transparent background plus camera / ray settings for scientific styles."""
    _set_transparent_background()
    cmd.set("orthoscopic", orthoscopic)
    cmd.set("field_of_view", field_of_view)
    if depth_cue is not None:
        cmd.set("depth_cue", depth_cue)
    if fog_start is not None:
        cmd.set("fog_start", fog_start)
    cmd.set("ray_shadow", ray_shadow)
    if ray_trace_gain is not None:
        cmd.set("ray_trace_gain", ray_trace_gain)
    if ray_shadows_mode is not None:
        _safe_ray_shadows(ray_shadows_mode)


def _end_scientific_style(selection, message):
    _base_quality()
    cmd.label("all", '""')
    _finish_style(selection)
    print(message)


def _show_coord_ball_and_stick(selection, spheres="coord_core"):
    cmd.show("sticks", selection)
    cmd.show("spheres", spheres)


def _define_editorial_minimal_colors():
    """Register palette used by the editorial-minimal style."""
    _define_scientific_colors()
    for name, rgb in {
        "mn_rose": [0.78, 0.38, 0.48],
        "sulfur_gold": [0.95, 0.68, 0.05],
        "hydrogen_soft": [0.88, 0.88, 0.86],
    }.items():
        try:
            cmd.set_color(name, rgb)
        except Exception:
            pass


def render_editorial_minimal(selection="all", highlight_bonds=""):
    """Editorial minimal white style for main-text mechanistic figures."""
    highlight_bonds = _normalize_none(highlight_bonds)
    sel = f"({selection})"

    _define_editorial_minimal_colors()
    cmd.hide("everything", sel)

    atoms = _get_coordinating_atoms(
        selection=selection,
        prefix="editorial",
    )
    sphere_atoms = "%s or %s or %s or %s" % (
        atoms["metal"],
        atoms["donors"],
        atoms["co_c"],
        atoms["important_h"],
    )

    cmd.show("sticks", sel)
    cmd.show("spheres", sphere_atoms)
    cmd.show("spheres", atoms["co_o"])
    cmd.hide("sticks", f"{sel} and elem H and not {atoms['important_h']}")
    cmd.hide("spheres", f"{sel} and elem H and not {atoms['important_h']}")

    _safe_set("sphere_scale", 0.60, atoms["metal"])
    _safe_set("sphere_scale", 0.36, atoms["donor_n"])
    _safe_set("sphere_scale", 0.39, atoms["donor_s"])
    _safe_set("sphere_scale", 0.26, atoms["co_c"])
    _safe_set("sphere_scale", 0.21, atoms["co_o"])
    _safe_set("sphere_scale", 0.26, atoms["important_h"])

    _safe_set("stick_radius", 0.12, sel)
    _safe_set("stick_radius", 0.15, atoms["coordination_core"])
    _safe_set("stick_radius", 0.07, f"{sel} and elem H")
    _safe_set("stick_quality", 30)

    cmd.color("mn_rose", atoms["metal"])
    cmd.color("sulfur_gold", atoms["donor_s"])
    cmd.color("sci_N_blue", atoms["donor_n"])
    cmd.color("sci_C_gray", atoms["co_c"])
    cmd.color("sci_O_red", atoms["co_o"])
    cmd.color("hydrogen_soft", atoms["important_h"])
    cmd.color(
        "sci_C_gray",
        f"{sel} and elem C and not {atoms['co_c']}",
    )

    _set_transparent_background()
    _safe_set("orthoscopic", 1)
    _safe_set("field_of_view", 45)

    _safe_set("ambient", 0.22)
    _safe_set("direct", 0.80)
    _safe_set("specular", 0.72)
    _safe_set("spec_reflect", 0.55)
    _safe_set("spec_power", 220)
    _safe_set("shininess", 75)
    _safe_set("ray_shadow", "on")
    _safe_set("ambient_occlusion_mode", 1)
    _safe_set("ambient_occlusion_scale", 15)
    _safe_set("ambient_occlusion_smooth", 10)

    _safe_set("specular", 0.85, atoms["metal"])
    _safe_set("shininess", 90, atoms["metal"])
    soft_atoms = "%s or %s or %s" % (
        atoms["donors"],
        atoms["co_c"],
        atoms["important_h"],
    )
    _safe_set("specular", 0.15, soft_atoms)
    _safe_set("shininess", 20, soft_atoms)

    if highlight_bonds:
        _apply_coordination_highlight_bonds(
            highlight_bonds,
            sel,
            atoms["metal"],
        )

    _base_quality()
    cmd.label("all", '""')
    _finish_style(selection)
    print("Editorial minimal white style applied.")


def _define_soft_ceramic_colors():
    """Register palette used by the soft-ceramic studio style."""
    for name, rgb in {
        "studio_background": [0.970, 0.965, 0.945],
        "ligand_ivory": [0.69, 0.67, 0.60],
        "carbonyl_cream": [0.76, 0.74, 0.67],
        "metal_bronze": [0.78, 0.44, 0.14],
        "sulfur_soft_gold": [0.88, 0.61, 0.08],
        "nitrogen_cobalt": [0.06, 0.17, 0.72],
        "oxygen_deep_red": [0.84, 0.04, 0.025],
        "hydrogen_warm": [0.94, 0.94, 0.91],
    }.items():
        try:
            cmd.set_color(name, rgb)
        except Exception:
            pass


def render_soft_ceramic(selection="all", highlight_bonds=""):
    """Soft ceramic / studio ball-and-stick style for coordination complexes.

    Metal centers are selected dynamically from common transition and main-group
    coordination metals rather than a hardcoded Mn-only selection.
    Optional ``-c`` bond pairs are applied via ``highlight_bonds`` without
    changing the shared CHEMSMART coordinate-highlight plumbing.
    """
    highlight_bonds = _normalize_none(highlight_bonds)
    sel = f"({selection})"

    _define_soft_ceramic_colors()
    cmd.hide("everything", sel)
    _safe_set("valence", 0)
    _safe_set("stick_ball", 0)

    atoms = _get_coordinating_atoms(
        selection=selection,
        prefix="soft_ceramic",
        include_nh_h=True,
    )

    cmd.bg_color("studio_background")
    cmd.color(
        "ligand_ivory",
        f"{sel} and (elem C+H) and not {atoms['co_c']} "
        f"and not {atoms['important_h']}",
    )
    cmd.color("carbonyl_cream", atoms["co_c"])
    cmd.color("metal_bronze", atoms["metal"])
    cmd.color("sulfur_soft_gold", atoms["donor_s"])
    cmd.color("nitrogen_cobalt", atoms["donor_n"])
    cmd.color("oxygen_deep_red", atoms["co_o"])
    cmd.color("hydrogen_warm", atoms["important_h"])

    cmd.show("sticks", sel)
    cmd.hide("sticks", f"{sel} and elem H and not {atoms['important_h']}")
    cmd.hide("spheres", f"{sel} and elem H and not {atoms['hydride']}")
    cmd.show("spheres", atoms["coordination_core"])

    _safe_set("sphere_scale", 0.56, atoms["metal"])
    _safe_set("sphere_scale", 0.40, atoms["donor_s"])
    _safe_set("sphere_scale", 0.37, atoms["donor_n"])
    _safe_set("sphere_scale", 0.28, atoms["co_c"])
    _safe_set("sphere_scale", 0.18, atoms["co_o"])
    _safe_set("sphere_scale", 0.25, atoms["hydride"])

    _safe_set("stick_radius", 0.115, sel)
    _safe_set(
        "stick_radius",
        0.155,
        "%s or %s or %s"
        % (atoms["metal"], atoms["donor_s"], atoms["donor_n"]),
    )
    _safe_set(
        "stick_radius",
        0.130,
        "%s or %s or %s" % (atoms["co_c"], atoms["co_o"], atoms["hydride"]),
    )

    _safe_set("stick_quality", 40)
    _safe_set("sphere_quality", 4)
    _safe_set("use_shaders", 1)
    _safe_set("two_sided_lighting", 1)
    _safe_set("light_count", 8)

    _safe_set("ambient", 0.30)
    _safe_set("direct", 0.70)
    _safe_set("reflect", 0.25)
    _safe_set("specular", 0.55)
    _safe_set("spec_reflect", 0.32)
    _safe_set("spec_power", 115)
    _safe_set("shininess", 58)

    _safe_set("ray_shadow", 1)
    _safe_set("ray_trace_mode", 0)
    _safe_set("ray_trace_gain", 0.06)
    _safe_set("antialias", 2)
    _safe_set("ray_trace_antialias", 2)
    _safe_set("ambient_occlusion_mode", 1)
    _safe_set("ambient_occlusion_scale", 12)
    _safe_set("ambient_occlusion_smooth", 10)
    _safe_set("ray_shadow_decay_factor", 0.25)
    _safe_set("ray_shadow_decay_range", 2.0)
    _safe_set("ray_opaque_background", 1)

    if highlight_bonds:
        _apply_coordination_highlight_bonds(
            highlight_bonds,
            sel,
            atoms["metal"],
        )

    cmd.label("all", '""')
    _finish_style(selection)
    print("Soft ceramic studio style applied.")


def render_neon_coordination_core(selection="all", highlight_bonds=""):
    """Neon coordination-core style for reactive centers and catalytic pockets."""
    _begin_scientific_style(selection, highlight_bonds=highlight_bonds)
    _show_coord_ball_and_stick(selection)
    _apply_coord_sphere_scales(0.085, 0.18, 0.20, 0.11, 0.62, 0.42)

    cmd.color("sci_C_dark", "elem C")
    cmd.color("sci_H_white", "elem H")
    cmd.color("neon_cyan", "elem N")
    cmd.color("neon_magenta", "elem O")
    cmd.color("neon_green", "elem S")
    cmd.color("sci_P_orange", "elem P")
    cmd.color("neon_green", "elem F+Cl+Br+I")
    cmd.color("metal_rose", "metal_atom")

    _apply_lighting(0.78, 0.45, 0.08, 0.95, 0.28, spec_power=180)
    _apply_view(
        orthoscopic=0,
        field_of_view=42,
        depth_cue=1,
        fog_start=0.18,
        ray_trace_gain=0.28,
    )
    _finish_scientific_style(
        selection,
        highlight_bonds=highlight_bonds,
        message="Neon coordination-core style applied.",
    )


def _define_matte_clay_colors():
    """Register muted ceramic/chalk colors used by the matte-clay style."""
    for name, rgb in {
        "mc_framework": [0.555, 0.515, 0.455],
        "mc_hydrogen": [0.875, 0.865, 0.820],
        "mc_nitrogen": [0.300, 0.395, 0.585],
        "mc_oxygen": [0.715, 0.285, 0.235],
        "mc_chalcogen": [0.790, 0.600, 0.205],
        "mc_pnictogen": [0.655, 0.420, 0.275],
        "mc_halogen": [0.455, 0.565, 0.455],
        "mc_metal_center": [0.690, 0.705, 0.745],
    }.items():
        try:
            cmd.set_color(name, rgb)
        except Exception:
            pass


def render_matte_clay(selection="all", highlight_bonds=""):
    """Matte clay style for soft graphical abstracts.

    Uses covalent-radius-ratio coordination spheres shared with soft-ceramic /
    editorial styles. Non-core hydrogens are hidden. Optional ``highlight_bonds``
    preserves CHEMSMART ``-c`` highlighting.
    """
    highlight_bonds = _normalize_none(highlight_bonds)
    sel = f"({selection})"

    _define_matte_clay_colors()
    cmd.hide("everything", sel)
    _safe_set("valence", 0)
    _safe_set("stick_ball", 0)

    atoms = _get_coordinating_atoms(
        selection=selection,
        prefix="matte_clay",
    )

    cmd.color(
        "mc_framework",
        f"{sel} and (elem C+H) and not {atoms['co_c']} "
        f"and not {atoms['hydride']}",
    )
    cmd.color("mc_framework", atoms["co_c"])
    cmd.color("mc_metal_center", atoms["metal"])
    cmd.color("mc_chalcogen", atoms["donor_s"])
    cmd.color("mc_nitrogen", atoms["donor_n"])
    cmd.color("mc_pnictogen", atoms["donor_p"])
    cmd.color("mc_oxygen", atoms["co_o"])
    cmd.color("mc_hydrogen", atoms["hydride"])
    cmd.color("mc_halogen", f"{sel} and elem F+Cl+Br+I")

    cmd.show("sticks", sel)
    cmd.hide("sticks", f"{sel} and elem H and not {atoms['hydride']}")
    cmd.hide("spheres", f"{sel} and elem H and not {atoms['hydride']}")
    cmd.show("spheres", atoms["coordination_core"])

    _safe_set("sphere_scale", 0.62, atoms["metal"])
    _safe_set("sphere_scale", 0.36, atoms["donor_s"])
    _safe_set("sphere_scale", 0.34, atoms["donor_n"])
    _safe_set("sphere_scale", 0.34, atoms["donor_p"])
    _safe_set("sphere_scale", 0.30, atoms["co_c"])
    _safe_set("sphere_scale", 0.28, atoms["co_o"])
    _safe_set("sphere_scale", 0.24, atoms["hydride"])

    _safe_set("stick_radius", 0.108, sel)
    _safe_set(
        "stick_radius",
        0.138,
        "%s or %s" % (atoms["metal"], atoms["coordination_core"]),
    )

    # Transparent background and zero-gloss clay material with soft AO.
    cmd.bg_color("white")
    _safe_set("ray_opaque_background", 0)
    _safe_set("specular", 0.0)
    _safe_set("specular_intensity", 0.0)
    _safe_set("spec_direct", 0.0)
    _safe_set("spec_reflect", 0.0)
    _safe_set("spec_power", 1.0)
    _safe_set("shininess", 0.0)
    _safe_set("reflect", 0.0)
    _safe_set("ray_transparency_specular", 0.0)
    _safe_set("ambient", 0.42)
    _safe_set("direct", 0.58)
    _safe_set("light_count", 3)
    _safe_set("light", [-0.40, -0.55, -1.00])
    _safe_set("two_sided_lighting", 1)
    _safe_set("ray_shadow", 1)
    _safe_set("ray_trace_mode", 0)
    _safe_set("ray_trace_gain", 0.0)
    _safe_set("ambient_occlusion_mode", 1)
    _safe_set("ambient_occlusion_scale", 18)
    _safe_set("ambient_occlusion_smooth", 14)
    _safe_set("ray_shadow_decay_factor", 0.28)
    _safe_set("ray_shadow_decay_range", 2.0)
    _safe_set("use_shaders", 1)
    _safe_set("antialias", 2)
    _safe_set("ray_trace_antialias", 2)
    _safe_set("sphere_quality", 4)
    _safe_set("stick_quality", 32)
    _safe_set("orthoscopic", 1)
    _safe_set("depth_cue", 0)
    _safe_set("ray_trace_fog", 0)

    if highlight_bonds:
        _apply_coordination_highlight_bonds(
            highlight_bonds,
            sel,
            atoms["metal"],
        )

    cmd.label("all", '""')
    _finish_style(selection)
    print("Matte clay style applied.")


def render_xray_wire(selection="all", highlight_bonds=""):
    """X-ray crystallography wire style for SI structure verification."""
    _begin_scientific_style(selection, highlight_bonds=highlight_bonds)
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

    _apply_lighting(0.05, 0.00, 0.70, 0.35, 0.00)
    _apply_view(
        orthoscopic=1,
        field_of_view=18,
        ray_shadow="off",
        depth_cue=0,
        ray_shadows_mode=None,
    )
    _finish_scientific_style(
        selection,
        highlight_bonds=highlight_bonds,
        message="X-ray wire style applied.",
    )


def render_steric_surface(selection="all", highlight_bonds=""):
    """Transparent steric surface style for catalyst pockets."""
    _begin_scientific_style(selection, highlight_bonds=highlight_bonds)
    _show_coord_ball_and_stick(selection)
    _apply_coord_sphere_scales(0.095, 0.14, 0.24, 0.12, 0.52, 0.34)
    _color_by_element(carbon="sci_C_gray", metal="metal_gold")

    cmd.show("surface", selection)
    cmd.color("surface_sky", selection)
    cmd.set("transparency", 0.68, selection)
    cmd.set("surface_quality", 1)
    cmd.set("surface_solvent", 1)

    # Re-apply atom colors after surface coloring paints the whole selection.
    _color_by_element(carbon="sci_C_gray", metal="metal_gold")

    _apply_lighting(0.45, 0.22, 0.32, 0.74, 0.20, spec_power=120)
    _apply_view(
        orthoscopic=1,
        field_of_view=30,
        depth_cue=1,
        fog_start=0.45,
    )
    _finish_scientific_style(
        selection,
        highlight_bonds=highlight_bonds,
        message="Transparent steric surface style applied.",
    )


def render_quasi_chemdraw_bold(selection="all", highlight_bonds=""):
    """Quasi-ChemDraw bold 3D style with formula-like clarity."""
    _begin_scientific_style(selection, highlight_bonds=highlight_bonds)
    _show_coord_ball_and_stick(selection, spheres="metal_atom or coord_core")

    cmd.hide("sticks", f"({selection}) and elem H and not coord_core")
    cmd.hide("spheres", f"({selection}) and elem H and not coord_core")

    _apply_coord_sphere_scales(0.18, 0.24, 0.23, 0.13, 0.50, 0.34)
    _color_by_element(carbon="sci_C_dark", metal="metal_gold")
    _apply_lighting(0.18, 0.02, 0.62, 0.42, 0.02, spec_power=70)
    _apply_view(
        orthoscopic=1,
        field_of_view=12,
        depth_cue=0,
        ray_shadow="off",
        ray_shadows_mode=None,
    )
    _finish_scientific_style(
        selection,
        highlight_bonds=highlight_bonds,
        message="Quasi-ChemDraw bold 3D style applied.",
    )


def render_labeled_coordination_core(selection="all", highlight_bonds=""):
    """Coordination-core labeled style with explicit element labels."""
    _begin_scientific_style(selection, highlight_bonds=highlight_bonds)
    _show_coord_ball_and_stick(selection)
    _apply_coord_sphere_scales(0.11, 0.16, 0.25, 0.13, 0.56, 0.38)
    _color_by_element(carbon="sci_C_gray", metal="metal_gold")
    _apply_lighting(0.55, 0.28, 0.30, 0.70, 0.18, spec_power=180)
    _apply_view(
        orthoscopic=1,
        field_of_view=28,
        depth_cue=1,
        fog_start=0.45,
    )

    _base_quality()
    cmd.label("all", '""')
    cmd.label("coord_core", "elem")
    cmd.set("label_size", 24)
    cmd.set("label_font_id", 7)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_position", [0, 0, 0])
    cmd.set("label_connector", 0)

    highlight_bonds = _normalize_none(highlight_bonds)
    if highlight_bonds:
        _apply_coordination_highlight_bonds(
            highlight_bonds,
            f"({selection})",
            "metal_atom",
        )

    _finish_style(selection)
    print("Labeled coordination-core style applied.")


cmd.extend("metallic_poster_render", metallic_poster_render)
cmd.extend(
    "render_comic_metallic_labeled_final", render_comic_metallic_labeled_final
)
cmd.extend("comic_render", comic_render)
cmd.extend("render_soft_cartoon", render_soft_cartoon)
cmd.extend("soft_cartoon_render", soft_cartoon_render)
cmd.extend("render_editorial_minimal", render_editorial_minimal)
cmd.extend("render_soft_ceramic", render_soft_ceramic)
cmd.extend("render_neon_coordination_core", render_neon_coordination_core)
cmd.extend("render_matte_clay", render_matte_clay)
cmd.extend("render_xray_wire", render_xray_wire)
cmd.extend("render_steric_surface", render_steric_surface)
cmd.extend("render_quasi_chemdraw_bold", render_quasi_chemdraw_bold)
cmd.extend(
    "render_labeled_coordination_core", render_labeled_coordination_core
)
