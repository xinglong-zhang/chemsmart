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

# ---------------------------------------------------------------------------
# ScientificStyle base class (shared coordination, palette, camera, lighting)
# ---------------------------------------------------------------------------


class ScientificStyle:
    """Base class for scientific PyMOL visualization styles.

    Subclasses declare ``name``, ``command``, ``prefix``, ``colors``, and
    implement ``render()`` with style-specific overrides only. Shared
    coordination selection, element categories, camera resets, and palette
    helpers live on this parent.
    """

    METAL_ELEMENTS = (
        "elem Li+Na+K+Rb+Cs+Be+Mg+Ca+Sr+Ba+Al+Ga+In+Sn+Pb+"
        "Sc+Ti+V+Cr+Mn+Fe+Co+Ni+Cu+Zn+Y+Zr+Nb+Mo+Tc+Ru+Rh+Pd+Ag+Cd+"
        "Hf+Ta+W+Re+Os+Ir+Pt+Au+Hg"
    )

    ELEMENT_CATEGORIES = {
        "carbon": "C",
        "hydrogen": "H",
        "nitrogen": "N",
        "oxygen": "O",
        "sulfur": "S",
        "phosphorus": "P",
        "C": "C",
        "H": "H",
        "N": "N",
        "O": "O",
        "S": "S",
        "P": "P",
        "N+O": "N+O",
        "S+P": "S+P",
        "heavy": "C+N+O+S+P+F+Cl+Br+I",
        "halogen": "F+Cl+Br+I",
        "Br+I": "Br+I",
        "chalcogen": "O+S+Se+Te",
        "pnictogen": "N+P+As+Sb",
        "metal": METAL_ELEMENTS,
    }

    SHARED_COLORS = {
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

    name = "scientific"
    command = None
    prefix = "coord"
    colors = {}
    include_nh_h = False
    message = "Scientific style applied."

    @staticmethod
    def safe_set(setting, value, selection=None):
        """Set a PyMOL parameter safely across different PyMOL versions."""
        try:
            if selection is None:
                cmd.set(setting, value)
            else:
                cmd.set(setting, value, selection)
        except Exception:
            pass

    @staticmethod
    def safe_ray_shadows(mode="light"):
        try:
            cmd.util.ray_shadows(mode)
        except Exception:
            pass

    @staticmethod
    def hide_distance_value_labels():
        """Hide numeric labels on ``-c`` bond-distance objects (``d1``, ``d2``, …)."""
        for name in ScientificStyle._distance_object_names():
            try:
                cmd.hide("labels", name)
            except Exception:
                pass

    @staticmethod
    def _distance_object_names():
        """Return sorted ``d1``, ``d2``, … measurement object names."""
        try:
            names = cmd.get_names("objects", enabled_only=0)
        except Exception:
            return []

        distance_names = []
        for name in names:
            suffix = name[1:] if name.startswith("d") else ""
            if suffix.isdigit():
                distance_names.append((int(suffix), name))
        distance_names.sort()
        return [name for _, name in distance_names]

    @staticmethod
    def pairs_from_distance_objects(selection="all"):
        """Return 1-based atom-index pairs stored in ``d1``, ``d2``, … objects."""
        try:
            state = cmd.get_state()
        except Exception:
            return []

        distance_names = ScientificStyle._distance_object_names()
        if not distance_names:
            return []

        xyz2idx = {}
        try:
            cmd.iterate_state(
                state,
                selection,
                "xyz2idx[(x, y, z)] = index",
                space={"xyz2idx": xyz2idx},
            )
        except Exception:
            return []

        pairs = []
        try:
            raw_objects = cmd.get_session(
                " ".join(distance_names), 1, 1, 0, 0
            )["names"]
        except Exception:
            return pairs

        for obj in raw_objects:
            try:
                points = obj[5][2][state - 1][1]
                if points is None:
                    continue
            except (KeyError, IndexError, TypeError):
                continue
            for point_index in range(0, len(points), 6):
                xyz1 = tuple(points[point_index : point_index + 3])
                xyz2 = tuple(points[point_index + 3 : point_index + 6])
                try:
                    pairs.append((xyz2idx[xyz1], xyz2idx[xyz2]))
                except KeyError:
                    continue
        return pairs

    @staticmethod
    def bond_atom_index_pairs(pairs):
        """Draw PyMOL sticks between 1-based atom-index pairs."""
        for atom_a, atom_b in pairs:
            try:
                cmd.bond(f"id {int(atom_a)}", f"id {int(atom_b)}")
            except Exception:
                pass

    @staticmethod
    def remove_distance_objects():
        """Delete ``d1``, ``d2``, … measurement objects created by ``-c``."""
        for name in ScientificStyle._distance_object_names():
            try:
                cmd.delete(name)
            except Exception:
                pass

    @classmethod
    def element_category_selection(cls, base_selection, category):
        """Build a PyMOL selection for an element category within ``base_selection``."""
        elements = cls.ELEMENT_CATEGORIES.get(category, category)
        if elements.startswith("elem "):
            return "(%s) and (%s)" % (base_selection, elements)
        return "(%s) and elem %s" % (base_selection, elements)

    @classmethod
    def apply_element_palette(cls, selection, palette, overrides=None):
        """Apply colors to element categories and optional named-selection overrides."""
        for category, color_name in palette.items():
            try:
                cmd.color(
                    color_name,
                    cls.element_category_selection(selection, category),
                )
            except Exception:
                pass
        if overrides:
            for sel_expr, color_name in overrides.items():
                try:
                    cmd.color(color_name, sel_expr)
                except Exception:
                    pass

    @staticmethod
    def selection_expr(selection):
        return "(%s)" % selection

    def define_colors(self):
        """Register ``self.colors`` in the current PyMOL session."""
        for color_name, rgb in self.colors.items():
            try:
                cmd.set_color(color_name, rgb)
            except Exception:
                pass

    @classmethod
    def define_shared_colors(cls):
        """Register the shared scientific palette (``sci_*``, ``metal_*``, etc.)."""
        for color_name, rgb in cls.SHARED_COLORS.items():
            try:
                cmd.set_color(color_name, rgb)
            except Exception:
                pass

    @classmethod
    def _parse_metal_symbols(cls, metal):
        text = (metal or cls.METAL_ELEMENTS).replace("elem", " ")
        return {
            token.strip()
            for token in text.replace("+", " ").split()
            if token.strip()
        }

    @staticmethod
    def _pymol_index_selection(sel, pymol_indices):
        indices = sorted({int(index) for index in pymol_indices})
        if not indices:
            return "none"
        return "%s and index %s" % (
            sel,
            "+".join(str(index) for index in indices),
        )

    @classmethod
    def build_coordination_atoms(
        cls,
        selection="all",
        prefix="coord",
        metal=None,
        include_nh_h=False,
    ):
        """Build named selections via covalent-radius-ratio coordination spheres."""
        sel = cls.selection_expr(selection)
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

        metal_symbols = cls._parse_metal_symbols(metal)
        metal_local = [
            idx
            for idx, element in enumerate(elements)
            if element in metal_symbols
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
            | {pymol_indices[idx] for idx in secondary_local}
        )

        cmd.select(
            names["metal"], cls._pymol_index_selection(sel, metal_pymol)
        )
        cmd.select(
            names["donor_s"], cls._pymol_index_selection(sel, donor_s_pymol)
        )
        cmd.select(
            names["donor_n"], cls._pymol_index_selection(sel, donor_n_pymol)
        )
        cmd.select(
            names["donor_p"], cls._pymol_index_selection(sel, donor_p_pymol)
        )
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
        cmd.select(names["co_c"], cls._pymol_index_selection(sel, co_c_pymol))
        cmd.select(names["co_o"], cls._pymol_index_selection(sel, co_o_pymol))
        cmd.select(
            names["hydride"], cls._pymol_index_selection(sel, hydride_pymol)
        )
        cmd.select(names["nh_h"], cls._pymol_index_selection(sel, nh_h_pymol))
        cmd.select(
            names["important_h"],
            cls._pymol_index_selection(sel, important_h_pymol),
        )
        cmd.select(
            names["coordination_core"],
            cls._pymol_index_selection(sel, core_pymol),
        )
        return names

    def select_coordination(self, selection="all"):
        """Build named coordination selections via radius-ratio helpers."""
        return self.build_coordination_atoms(
            selection=selection,
            prefix=self.prefix,
            include_nh_h=self.include_nh_h,
        )

    def first_shell(self, atoms):
        """Return selection expression for non-metal coordination-core atoms."""
        return "%s and not %s" % (
            atoms["coordination_core"],
            atoms["metal"],
        )

    @staticmethod
    def coordination_sphere_atoms(atoms, include_hydride=True):
        """Join role selections that receive explicit sphere scales."""
        parts = [atoms["metal"], atoms["donors"], atoms["co_c"], atoms["co_o"]]
        if include_hydride:
            parts.append(atoms["hydride"])
        return " or ".join(part for part in parts if part and part != "none")

    @classmethod
    def set_transparent_background(cls):
        """Configure ray-traced PNG export with a transparent background."""
        cmd.bg_color("white")
        cls.safe_set("ray_opaque_background", 0)

    @classmethod
    def apply_base_quality(cls):
        cls.safe_set("antialias", 2)
        cls.safe_set("ray_trace_antialias", 2)
        cls.safe_set("sphere_quality", 3)
        cls.safe_set("stick_quality", 30)
        cls.safe_set("two_sided_lighting", 1)
        cls.safe_set("use_shaders", 1)
        cls.safe_set("depth_cue", 1)
        cls.safe_set("label_connector", 0)

    @classmethod
    def finish_camera(cls, selection, buffer=2.0):
        cmd.zoom(selection, buffer=buffer)
        cmd.orient(selection)
        cmd.refresh()

    @staticmethod
    def apply_lighting(
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

    @classmethod
    def apply_transparent_view(
        cls,
        orthoscopic,
        field_of_view,
        ray_shadow="on",
        depth_cue=None,
        fog_start=None,
        ray_trace_gain=None,
        ray_shadows_mode="light",
    ):
        """Transparent background plus camera / ray settings for scientific styles."""
        cls.set_transparent_background()
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
            cls.safe_ray_shadows(ray_shadows_mode)

    def apply_style_palette(self, selection, palette, overrides=None):
        """Apply element-category colors from a category-to-color-name mapping."""
        self.apply_element_palette(selection, palette, overrides=overrides)

    def _safe_set(self, setting, value, selection=None, category=None):
        """Set a PyMOL parameter, optionally scoped to an element category."""
        if selection is None and category is None:
            self.safe_set(setting, value)
            return
        target = self.element_category_selection(selection, category)
        self.safe_set(setting, value, target)

    def frame(self, selection, core_name, zoom_buffer=1.4):
        """Orient / center / zoom around the coordination environment."""
        sel = self.selection_expr(selection)
        cmd.label("all", '""')
        cmd.orient(sel)
        cmd.center(core_name)
        try:
            cmd.origin(core_name)
        except Exception:
            pass
        cmd.zoom(sel, zoom_buffer)
        try:
            cmd.rebuild()
        except Exception:
            pass
        cmd.refresh()

    def finalize(self):
        """Post-render cleanup shared by all scientific styles."""
        self.hide_distance_value_labels()

    def finish_default(self, selection):
        """Zoom/orient via ``finish_camera`` and print ``self.message``."""
        self.finalize()
        self.finish_camera(selection)
        print(self.message)

    def apply_illustrated_camera(self, field_of_view=25, depth_cue=0, fog=0.0):
        """Flat orthographic camera without forcing a transparent background."""
        self.safe_set("orthoscopic", 1)
        self.safe_set("ray_orthoscopic", 1)
        self.safe_set("depth_cue", depth_cue)
        self.safe_set("fog", fog)
        self.safe_set("field_of_view", field_of_view)

    def apply_soft_shadows(self, decay_factor=0.35, decay_range=2.5):
        self.safe_set("ray_shadow", 1)
        self.safe_set("ray_shadow_decay_factor", decay_factor)
        self.safe_set("ray_shadow_decay_range", decay_range)

    def apply_ambient_occlusion(self, scale=13, smooth=15):
        self.safe_set("ambient_occlusion_mode", 1)
        self.safe_set("ambient_occlusion_scale", scale)
        self.safe_set("ambient_occlusion_smooth", smooth)

    def apply_coordination_sci_palette(self, sel, atoms):
        """Apply the shared scientific element palette to coordination roles."""
        self.apply_style_palette(
            sel,
            {
                "C": "sci_C_gray",
                "H": "sci_H_white",
                "N": "sci_N_blue",
                "O": "sci_O_red",
                "S": "sci_S_yellow",
                "P": "sci_P_orange",
                "halogen": "sci_halogen",
            },
            overrides={atoms["metal"]: "metal_gold"},
        )
        cmd.color("sci_N_blue", atoms["donor_n"])
        cmd.color("sci_S_yellow", atoms["donor_s"])
        cmd.color("sci_P_orange", atoms["donor_p"])
        cmd.color("sci_C_gray", atoms["co_c"])
        cmd.color("sci_O_red", atoms["co_o"])
        cmd.color("sci_H_white", atoms["hydride"])

    def render(self, selection="all", **kwargs):
        raise NotImplementedError(
            "%s must implement render()" % type(self).__name__
        )


# ---------------------------------------------------------------------------
# Style implementations
# ---------------------------------------------------------------------------


class MetallicPosterStyle(ScientificStyle):
    """Glossy metallic poster-style rendering for cover figures."""

    name = "metallic_poster"
    command = "metallic_poster_render"
    prefix = "poster"
    include_nh_h = True
    message = "Metallic poster style applied."
    colors = {
        "poster_carbon": [0.62, 0.62, 0.62],
        "poster_hydrogen": [0.96, 0.96, 0.96],
        "poster_nitrogen": [0.12, 0.22, 0.95],
        "poster_oxygen": [0.95, 0.06, 0.04],
        "poster_sulfur": [0.95, 0.72, 0.10],
        "poster_phosphorus": [1.00, 0.48, 0.08],
        "poster_halogen": [0.20, 0.82, 0.32],
        "poster_mn_gold": [0.92, 0.60, 0.22],
        "poster_metal_gray": [0.78, 0.78, 0.84],
        "poster_label_white": [1.00, 1.00, 1.00],
        "poster_label_black": [0.02, 0.02, 0.02],
    }

    ELEMENT_PALETTE = {
        "C": "poster_carbon",
        "H": "poster_hydrogen",
        "N": "poster_nitrogen",
        "O": "poster_oxygen",
        "S": "poster_sulfur",
        "P": "poster_phosphorus",
        "halogen": "poster_halogen",
    }

    @classmethod
    def metal_palette_overrides(cls, selection, metal):
        """Return poster/comic metal-center palette overrides."""
        sel = cls.selection_expr(selection)
        return {
            f"{sel} and ({cls.METAL_ELEMENTS}) and not ({metal})": (
                "poster_metal_gray"
            ),
            metal: "poster_mn_gold",
        }

    @staticmethod
    def metal_element_label(metal_selection):
        """Return the element symbol of the first atom in a metal selection."""
        model = cmd.get_model(metal_selection)
        if not model.atom:
            return "?"
        symbol = model.atom[0].symbol.strip()
        return symbol or "?"

    def render(self, selection="all"):
        sel = f"({selection})"
        self.define_colors()
        atoms = self.select_coordination(selection)
        shell = self.first_shell(atoms)
        core = atoms["coordination_core"]
        metal = atoms["metal"]

        self.set_transparent_background()
        self.apply_base_quality()
        self.apply_lighting(
            0.85,
            0.60,
            0.22,
            0.82,
            0.45,
            spec_power=260,
            shininess=90,
        )
        self.safe_set("fog_start", 0.60)
        self.safe_set("orthoscopic", 1)
        self.safe_set("field_of_view", 35)
        self.safe_set("depth_cue", 1)
        self.safe_set("light_count", 8)
        self.safe_set("ray_trace_mode", 0)
        self.safe_set("ray_trace_gain", 0.08)
        self.safe_set("ray_trace_disco_factor", 1)
        self.safe_ray_shadows("light")
        self.apply_soft_shadows(decay_factor=0.25, decay_range=2.0)
        self.apply_ambient_occlusion(scale=18, smooth=12)
        self.safe_set("sphere_quality", 4)

        cmd.hide("everything", sel)
        cmd.show("sticks", sel)
        self.safe_set("stick_radius", 0.13, sel)
        self.safe_set("stick_quality", 30, sel)
        self.safe_set("stick_ball", 0, sel)
        self.safe_set("valence", 0, sel)
        cmd.hide("sticks", f"{sel} and elem H and not ({shell})")
        cmd.show("spheres", core)
        self.safe_set("sphere_scale", 0.62, metal)
        self.safe_set("sphere_scale", 0.50, atoms["donors"])
        self.safe_set("sphere_scale", 0.34, atoms["important_h"])
        self.safe_set("stick_radius", 0.16, core)

        self.apply_style_palette(
            sel,
            self.ELEMENT_PALETTE,
            overrides=self.metal_palette_overrides(selection, metal),
        )

        label_prefix = "metallic_poster_labels"
        cmd.delete(label_prefix + "*")
        label_colors = {"N": "poster_label_white", "O": "poster_label_white"}
        counters = {}
        for atom in cmd.get_model(core).atom:
            elem = atom.symbol.strip() or atom.name.strip()[0]
            counters[elem] = counters.get(elem, 0) + 1
            obj_name = f"{label_prefix}_{elem}"
            cmd.pseudoatom(
                object=obj_name,
                name=f"L_{elem}_{counters[elem]}",
                pos=atom.coord,
                label=elem,
            )
        for elem in counters:
            obj_name = f"{label_prefix}_{elem}"
            cmd.hide("everything", obj_name)
            cmd.show("labels", obj_name)
            self.safe_set("label_position", [0, 0, 0], obj_name)
            self.safe_set("label_font_id", 7, obj_name)
            self.safe_set("label_size", 24.0, obj_name)
            self.safe_set("label_connector", 0, obj_name)
            self.safe_set("label_shadow_mode", 2, obj_name)
            label_color = label_colors.get(elem, "poster_label_black")
            outline_color = (
                "poster_label_black"
                if label_color == "poster_label_white"
                else "poster_label_white"
            )
            self.safe_set("label_color", label_color, obj_name)
            self.safe_set("label_outline_color", outline_color, obj_name)

        self.finish_default(selection)


class ComicMetallicStyle(ScientificStyle):
    """Comic metallic ball-and-stick rendering with black outlines."""

    name = "comic_metallic"
    command = "render_comic_metallic_labeled_final"
    prefix = "comic"
    message = "Comic metallic style applied. Run 'ray 1200, 1200' to generate the image."

    colors = MetallicPosterStyle.colors

    DONOR_ROLES = (
        ("donor_n", "N"),
        ("donor_p", "P"),
        ("donor_s", "S"),
    )

    PERIPHERAL_SCALES = (
        ("heavy", 0.25),
        ("H", 0.15),
    )

    OUTLINE_COLOR = "black"
    RAY_TRACE_GAIN = 0.6
    METAL_SPHERE_SCALE = 0.42
    STICK_RADIUS = 0.14

    def render(self, selection="all", **kwargs):
        sel = f"({selection})"
        self.define_colors()
        atoms = self.select_coordination(selection)
        metal = atoms["metal"]

        self.apply_base_quality()
        self.set_transparent_background()

        cmd.hide("everything", sel)
        cmd.show("sticks", sel)
        cmd.show("spheres", sel)
        self.safe_set("stick_radius", self.STICK_RADIUS, sel)
        for category, scale in self.PERIPHERAL_SCALES:
            self._safe_set(
                "sphere_scale", scale, selection=sel, category=category
            )
        if cmd.count_atoms(metal) > 0:
            self.safe_set("sphere_scale", self.METAL_SPHERE_SCALE, metal)

        self.apply_style_palette(
            sel,
            MetallicPosterStyle.ELEMENT_PALETTE,
            overrides=MetallicPosterStyle.metal_palette_overrides(
                selection, metal
            ),
        )

        highlight_pairs = self.pairs_from_distance_objects(selection)
        if highlight_pairs:
            self.bond_atom_index_pairs(highlight_pairs)
            self.remove_distance_objects()
        elif cmd.count_atoms(metal) > 0:
            for donor_key, _element in self.DONOR_ROLES:
                donors = atoms[donor_key]
                if donors != "none" and cmd.count_atoms(donors) > 0:
                    cmd.bond(metal, donors)

        self.apply_lighting(
            0.85,
            0.70,
            0.25,
            0.75,
            0.0,
            spec_power=300,
            shininess=90,
        )
        self.safe_set("ray_trace_mode", 1)
        self.safe_set("ray_trace_color", self.OUTLINE_COLOR)
        self.safe_set("ray_trace_gain", self.RAY_TRACE_GAIN)

        cmd.label("all", '""')
        if cmd.count_atoms(metal) > 0:
            metal_label = MetallicPosterStyle.metal_element_label(metal)
            cmd.label(metal, f'"{metal_label}"')
        for donor_key, element in self.DONOR_ROLES:
            donors = atoms[donor_key]
            if donors != "none" and cmd.count_atoms(donors) > 0:
                cmd.label(donors, f'"{element}"')
        for setting, value in (
            ("label_position", [0.0, 0.0, 2.0]),
            ("label_shadow_mode", 0),
            ("ray_label_specular", 0),
            ("label_font_id", 7),
            ("label_color", "white"),
            ("label_size", 28),
        ):
            self.safe_set(setting, value)

        self.safe_set("orthoscopic", 1)
        self.finish_default(selection)


class SoftCartoonStyle(ScientificStyle):
    """Soft cartoon ball-and-stick with muted pastels and rounded outlines.

    Uses :attr:`ScientificStyle.METAL_ELEMENTS` and radius-ratio coordination
    (plus 1.6 Å geometric partner expansion) via ``select_coordination``.
    """

    name = "soft_cartoon"
    command = "render_soft_cartoon"
    prefix = "sc"
    include_nh_h = True
    message = "Soft cartoon style applied."
    colors = {
        "sc_background": [0.965, 0.950, 0.915],
        "sc_outline": [0.18, 0.17, 0.20],
        "sc_carbon": [0.43, 0.45, 0.50],
        "sc_hydrogen": [0.91, 0.90, 0.86],
        "sc_nitrogen": [0.37, 0.49, 0.84],
        "sc_oxygen": [0.87, 0.35, 0.31],
        "sc_sulfur": [0.91, 0.66, 0.20],
        "sc_phosphorus": [0.91, 0.49, 0.27],
        "sc_halogen": [0.43, 0.68, 0.49],
        "sc_metal": [0.66, 0.50, 0.80],
    }

    def render(self, selection="all", background=None):
        sel = f"({selection})"
        self.define_colors()
        atoms = self.select_coordination(selection)
        shell = self.first_shell(atoms)
        peripheral_heavy = (
            f"{sel} and not elem H and not {atoms['coordination_core']}"
        )
        sphere_heavy = f"{sel} and not elem H"
        coordinating_h = f"({shell}) and elem H"

        cmd.bg_color("sc_background")
        cmd.hide("everything", sel)
        cmd.show("sticks", sel)
        cmd.show("spheres", sphere_heavy)
        cmd.show("spheres", coordinating_h)
        cmd.hide("sticks", f"{sel} and elem H and not ({shell})")
        cmd.hide("spheres", f"{sel} and elem H and not ({shell})")

        self.safe_set("stick_radius", 0.135, sel)
        self.safe_set("stick_h_scale", 0.72)
        self.safe_set("stick_quality", 32)
        self.safe_set("stick_ball", 0)
        self.safe_set("smooth_half_bonds", 1)
        self.safe_set("render_as_cylinders", 1)
        self.safe_set("stick_color", -1, sel)
        try:
            cmd.set_bond("stick_radius", 0.165, atoms["metal"], shell)
        except Exception:
            pass

        self.safe_set("sphere_scale", 0.215, peripheral_heavy)
        self.safe_set("sphere_scale", 0.40, atoms["metal"])
        # Shell heavy atoms default to framework scale; category rules override below.
        self.safe_set("sphere_scale", 0.215, f"({shell}) and not elem H")
        for category, scale in (
            ("N+O", 0.275),
            ("S+P", 0.255),
            ("halogen", 0.235),
        ):
            self._safe_set(
                "sphere_scale", scale, selection=shell, category=category
            )
        self.safe_set("sphere_scale", 0.175, coordinating_h)
        self.safe_set("sphere_quality", 4)
        self.safe_set(
            "sphere_color", -1, f"{sphere_heavy} or {coordinating_h}"
        )

        self.apply_style_palette(
            sel,
            {
                "C": "sc_carbon",
                "H": "sc_hydrogen",
                "N": "sc_nitrogen",
                "O": "sc_oxygen",
                "S": "sc_sulfur",
                "P": "sc_phosphorus",
                "halogen": "sc_halogen",
            },
            overrides={atoms["metal"]: "sc_metal"},
        )

        self.apply_lighting(
            0.22, 0.10, 0.42, 0.58, 0.06, spec_power=28, shininess=18
        )
        self.safe_set("light_count", 5)
        self.safe_set("two_sided_lighting", 1)
        self.safe_set("use_shaders", 1)
        self.apply_soft_shadows(decay_factor=0.35, decay_range=2.5)
        self.apply_ambient_occlusion(scale=13, smooth=15)
        self.safe_set("ray_trace_mode", 1)
        self.safe_set("ray_trace_color", "sc_outline")
        self.safe_set("ray_trace_gain", 0.025)
        self.safe_set("ray_trace_disco_factor", 0.25)
        self.apply_illustrated_camera(field_of_view=25)
        self.safe_set("antialias", 2)
        self.safe_set("ray_trace_antialias", 2)
        self.safe_set("dash_round_ends", 1)
        self.safe_set("dash_radius", 0.065)
        self.safe_set("dash_length", 0.18)
        self.safe_set("dash_gap", 0.16)
        self.safe_set("dash_color", "sc_outline")

        self.frame(selection, atoms["coordination_core"], zoom_buffer=1.6)
        self.finalize()
        print(self.message)


class EditorialMinimalStyle(ScientificStyle):
    """Editorial minimal white style for main-text mechanistic figures."""

    name = "editorial_minimal"
    command = "render_editorial_minimal"
    prefix = "editorial"
    message = "Editorial minimal white style applied."
    colors = {
        "mn_rose": [0.78, 0.38, 0.48],
        "sulfur_gold": [0.95, 0.68, 0.05],
        "hydrogen_soft": [0.88, 0.88, 0.86],
    }

    def render(self, selection="all"):
        sel = f"({selection})"

        self.define_shared_colors()
        self.define_colors()
        cmd.hide("everything", sel)

        atoms = self.select_coordination(selection)
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

        self.safe_set("sphere_scale", 0.60, atoms["metal"])
        self.safe_set("sphere_scale", 0.36, atoms["donor_n"])
        self.safe_set("sphere_scale", 0.39, atoms["donor_s"])
        self.safe_set("sphere_scale", 0.26, atoms["co_c"])
        self.safe_set("sphere_scale", 0.21, atoms["co_o"])
        self.safe_set("sphere_scale", 0.26, atoms["important_h"])

        self.safe_set("stick_radius", 0.12, sel)
        self.safe_set("stick_radius", 0.15, atoms["coordination_core"])
        self.safe_set("stick_radius", 0.07, f"{sel} and elem H")
        self.safe_set("stick_quality", 30)

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

        self.set_transparent_background()
        self.safe_set("orthoscopic", 1)
        self.safe_set("field_of_view", 45)

        self.safe_set("ambient", 0.22)
        self.safe_set("direct", 0.80)
        self.safe_set("specular", 0.72)
        self.safe_set("spec_reflect", 0.55)
        self.safe_set("spec_power", 220)
        self.safe_set("shininess", 75)
        self.safe_set("ray_shadow", "on")
        self.safe_set("ambient_occlusion_mode", 1)
        self.safe_set("ambient_occlusion_scale", 15)
        self.safe_set("ambient_occlusion_smooth", 10)

        self.safe_set("specular", 0.85, atoms["metal"])
        self.safe_set("shininess", 90, atoms["metal"])
        soft_atoms = "%s or %s or %s" % (
            atoms["donors"],
            atoms["co_c"],
            atoms["important_h"],
        )
        self.safe_set("specular", 0.15, soft_atoms)
        self.safe_set("shininess", 20, soft_atoms)

        self.apply_base_quality()
        cmd.label("all", '""')
        self.finish_default(selection)


class SoftCeramicStyle(ScientificStyle):
    """Soft ceramic / studio ball-and-stick style for coordination complexes."""

    name = "soft_ceramic"
    command = "render_soft_ceramic"
    prefix = "soft_ceramic"
    include_nh_h = True
    message = "Soft ceramic studio style applied."
    colors = {
        "studio_background": [0.970, 0.965, 0.945],
        "ligand_ivory": [0.69, 0.67, 0.60],
        "carbonyl_cream": [0.76, 0.74, 0.67],
        "metal_bronze": [0.78, 0.44, 0.14],
        "sulfur_soft_gold": [0.88, 0.61, 0.08],
        "nitrogen_cobalt": [0.06, 0.17, 0.72],
        "oxygen_deep_red": [0.84, 0.04, 0.025],
        "hydrogen_warm": [0.94, 0.94, 0.91],
    }

    def render(self, selection="all"):
        sel = f"({selection})"

        self.define_colors()
        cmd.hide("everything", sel)
        self.safe_set("valence", 0)
        self.safe_set("stick_ball", 0)

        atoms = self.select_coordination(selection)

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

        self.safe_set("sphere_scale", 0.56, atoms["metal"])
        self.safe_set("sphere_scale", 0.40, atoms["donor_s"])
        self.safe_set("sphere_scale", 0.37, atoms["donor_n"])
        self.safe_set("sphere_scale", 0.28, atoms["co_c"])
        self.safe_set("sphere_scale", 0.18, atoms["co_o"])
        self.safe_set("sphere_scale", 0.25, atoms["hydride"])

        self.safe_set("stick_radius", 0.115, sel)
        self.safe_set(
            "stick_radius",
            0.155,
            "%s or %s or %s"
            % (atoms["metal"], atoms["donor_s"], atoms["donor_n"]),
        )
        self.safe_set(
            "stick_radius",
            0.130,
            "%s or %s or %s"
            % (atoms["co_c"], atoms["co_o"], atoms["hydride"]),
        )

        self.safe_set("stick_quality", 40)
        self.safe_set("sphere_quality", 4)
        self.safe_set("use_shaders", 1)
        self.safe_set("two_sided_lighting", 1)
        self.safe_set("light_count", 8)

        self.safe_set("ambient", 0.30)
        self.safe_set("direct", 0.70)
        self.safe_set("reflect", 0.25)
        self.safe_set("specular", 0.55)
        self.safe_set("spec_reflect", 0.32)
        self.safe_set("spec_power", 115)
        self.safe_set("shininess", 58)

        self.safe_set("ray_shadow", 1)
        self.safe_set("ray_trace_mode", 0)
        self.safe_set("ray_trace_gain", 0.06)
        self.safe_set("antialias", 2)
        self.safe_set("ray_trace_antialias", 2)
        self.safe_set("ambient_occlusion_mode", 1)
        self.safe_set("ambient_occlusion_scale", 12)
        self.safe_set("ambient_occlusion_smooth", 10)
        self.safe_set("ray_shadow_decay_factor", 0.25)
        self.safe_set("ray_shadow_decay_range", 2.0)
        self.safe_set("ray_opaque_background", 1)

        cmd.label("all", '""')
        self.finish_default(selection)


class NeonCoordinationCoreStyle(ScientificStyle):
    """Neon coordination-core style for reactive centers and catalytic pockets."""

    name = "neon_coordination_core"
    command = "render_neon_coordination_core"
    prefix = "ncc"
    message = "Neon coordination-core style applied."
    colors = {
        "ncc_background": [0.008, 0.012, 0.026],
        "ncc_carbon": [0.20, 0.23, 0.30],
        "ncc_hydrogen": [0.88, 0.92, 1.00],
        "ncc_metal_c": [0.20, 1.00, 0.56],
        "ncc_nitrogen": [0.02, 0.88, 1.00],
        "ncc_oxygen": [1.00, 0.10, 0.52],
        "ncc_sulfur": [1.00, 0.76, 0.08],
        "ncc_phosphor": [1.00, 0.46, 0.08],
        "ncc_halogen": [0.76, 0.28, 1.00],
    }

    def render(self, selection="all"):
        sel = f"({selection})"

        self.define_colors()
        atoms = self.select_coordination(selection)
        first_shell = self.first_shell(atoms)

        self.safe_set("opaque_background", 0)
        self.safe_set("ray_opaque_background", 0)
        self.safe_set("show_alpha_checker", 1)

        cmd.hide("everything", sel)
        cmd.show("sticks", sel)
        cmd.show("spheres", atoms["coordination_core"])
        cmd.hide("sticks", f"{sel} and elem H and not ({first_shell})")

        self.safe_set("stick_radius", 0.10, sel)
        self.safe_set("stick_h_scale", 0.70)
        self.safe_set("stick_quality", 32)
        self.safe_set("stick_transparency", 0.0)
        self.safe_set("stick_ball", 0)
        self.safe_set("smooth_half_bonds", 1)
        try:
            cmd.set_bond("stick_radius", 0.145, atoms["metal"], first_shell)
        except Exception:
            pass

        self.safe_set("sphere_scale", 0.44, atoms["metal"])
        self.safe_set(
            "sphere_scale",
            0.28,
            f"({first_shell}) and not (elem H+Br+I)",
        )
        self._safe_set(
            "sphere_scale", 0.22, selection=first_shell, category="Br+I"
        )
        self._safe_set(
            "sphere_scale", 0.20, selection=first_shell, category="H"
        )
        self.safe_set("sphere_transparency", 0.0, atoms["coordination_core"])
        self.safe_set("sphere_quality", 4)

        self.apply_style_palette(
            sel,
            {
                "C": "ncc_carbon",
                "H": "ncc_hydrogen",
                "N": "ncc_nitrogen",
                "O": "ncc_oxygen",
                "S": "ncc_sulfur",
                "P": "ncc_phosphor",
                "halogen": "ncc_halogen",
            },
            overrides={atoms["metal"]: "ncc_metal_c"},
        )

        self.safe_set("ambient", 0.10)
        self.safe_set("direct", 0.82)
        self.safe_set("reflect", 0.20)
        self.safe_set("specular", 0.78)
        self.safe_set("spec_reflect", 0.50)
        self.safe_set("spec_power", 170)
        self.safe_set("shininess", 65)
        self.safe_set("light_count", 8)
        self.safe_set("two_sided_lighting", 1)
        self.safe_set("ray_shadow", 1)
        self.safe_set("ray_shadow_decay_factor", 0.22)
        self.safe_set("ray_shadow_decay_range", 2.0)
        self.safe_set("ambient_occlusion_mode", 1)
        self.safe_set("ambient_occlusion_scale", 18)
        self.safe_set("ambient_occlusion_smooth", 10)
        self.safe_set("ray_trace_mode", 0)
        self.safe_set("ray_trace_gain", 0.06)
        self.safe_set("depth_cue", 0)
        self.safe_set("fog", 0.0)
        self.safe_set("orthoscopic", 1)
        self.safe_set("ray_orthoscopic", 1)
        self.safe_set("use_shaders", 1)
        self.safe_set("render_as_cylinders", 1)
        self.safe_set("antialias", 2)
        self.safe_set("ray_trace_antialias", 2)

        self.frame(selection, atoms["coordination_core"], zoom_buffer=1.7)
        self.finalize()
        print(self.message)


class MatteClayStyle(ScientificStyle):
    """Matte clay style for soft graphical abstracts."""

    name = "matte_clay"
    command = "render_matte_clay"
    prefix = "matte_clay"
    message = "Matte clay style applied."
    colors = {
        "mc_framework": [0.555, 0.515, 0.455],
        "mc_hydrogen": [0.875, 0.865, 0.820],
        "mc_nitrogen": [0.300, 0.395, 0.585],
        "mc_oxygen": [0.715, 0.285, 0.235],
        "mc_chalcogen": [0.790, 0.600, 0.205],
        "mc_pnictogen": [0.655, 0.420, 0.275],
        "mc_halogen": [0.455, 0.565, 0.455],
        "mc_metal_center": [0.690, 0.705, 0.745],
    }

    def render(self, selection="all"):
        sel = f"({selection})"

        self.define_colors()
        cmd.hide("everything", sel)
        self.safe_set("valence", 0)
        self.safe_set("stick_ball", 0)

        atoms = self.select_coordination(selection)

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
        self.apply_style_palette(sel, {"halogen": "mc_halogen"})

        cmd.show("sticks", sel)
        cmd.hide("sticks", f"{sel} and elem H and not {atoms['hydride']}")
        cmd.hide("spheres", f"{sel} and elem H and not {atoms['hydride']}")
        cmd.show("spheres", atoms["coordination_core"])

        self.safe_set("sphere_scale", 0.62, atoms["metal"])
        self.safe_set("sphere_scale", 0.36, atoms["donor_s"])
        self.safe_set("sphere_scale", 0.34, atoms["donor_n"])
        self.safe_set("sphere_scale", 0.34, atoms["donor_p"])
        self.safe_set("sphere_scale", 0.30, atoms["co_c"])
        self.safe_set("sphere_scale", 0.28, atoms["co_o"])
        self.safe_set("sphere_scale", 0.24, atoms["hydride"])

        self.safe_set("stick_radius", 0.108, sel)
        self.safe_set(
            "stick_radius",
            0.138,
            "%s or %s" % (atoms["metal"], atoms["coordination_core"]),
        )

        # Transparent background and zero-gloss clay material with soft AO.
        cmd.bg_color("white")
        self.safe_set("ray_opaque_background", 0)
        self.safe_set("specular", 0.0)
        self.safe_set("specular_intensity", 0.0)
        self.safe_set("spec_direct", 0.0)
        self.safe_set("spec_reflect", 0.0)
        self.safe_set("spec_power", 1.0)
        self.safe_set("shininess", 0.0)
        self.safe_set("reflect", 0.0)
        self.safe_set("ray_transparency_specular", 0.0)
        self.safe_set("ambient", 0.42)
        self.safe_set("direct", 0.58)
        self.safe_set("light_count", 3)
        self.safe_set("light", [-0.40, -0.55, -1.00])
        self.safe_set("two_sided_lighting", 1)
        self.safe_set("ray_shadow", 1)
        self.safe_set("ray_trace_mode", 0)
        self.safe_set("ray_trace_gain", 0.0)
        self.safe_set("ambient_occlusion_mode", 1)
        self.safe_set("ambient_occlusion_scale", 18)
        self.safe_set("ambient_occlusion_smooth", 14)
        self.safe_set("ray_shadow_decay_factor", 0.28)
        self.safe_set("ray_shadow_decay_range", 2.0)
        self.safe_set("use_shaders", 1)
        self.safe_set("antialias", 2)
        self.safe_set("ray_trace_antialias", 2)
        self.safe_set("sphere_quality", 4)
        self.safe_set("stick_quality", 32)
        self.safe_set("orthoscopic", 1)
        self.safe_set("depth_cue", 0)
        self.safe_set("ray_trace_fog", 0)

        cmd.label("all", '""')
        self.finish_default(selection)


class XrayWireStyle(ScientificStyle):
    """X-ray crystallography wire style for SI structure verification."""

    name = "xray_wire"
    command = "render_xray_wire"
    prefix = "xray"
    message = "X-ray wire style applied."
    colors = {
        "xw_charcoal": [0.2, 0.2, 0.2],
        "xw_metal": [0.95, 0.55, 0.15],
        "xw_sulfur": [0.85, 0.75, 0.25],
        "xw_nitrogen": [0.2, 0.4, 0.8],
        "xw_oxygen": [0.85, 0.25, 0.25],
        "xw_hydrogen": [1.0, 1.0, 1.0],
    }

    def render(self, selection="all"):
        sel = f"({selection})"
        self.define_colors()
        atoms = self.select_coordination(selection)
        center_metal = atoms["metal"]

        cmd.hide("everything", sel)
        cmd.show("sticks", sel)
        cmd.show("spheres", center_metal)
        cmd.hide("lines", sel)
        cmd.hide("nonbonded", sel)

        self.safe_set("stick_radius", 0.18, sel)
        self.safe_set("stick_ball", 1)
        self.safe_set("stick_ball_ratio", 1.0)
        self.safe_set("sphere_scale", 0.6, center_metal)

        self.apply_style_palette(
            sel,
            {
                "C": "xw_charcoal",
                "H": "xw_hydrogen",
                "N": "xw_nitrogen",
                "O": "xw_oxygen",
                "S": "xw_sulfur",
            },
            overrides={center_metal: "xw_metal"},
        )

        self.safe_set("ray_trace_mode", 3)
        self.safe_set("ray_trace_gain", 0.05)
        self.safe_set("ambient", 0.5)
        self.safe_set("direct", 0.6)
        self.safe_set("reflect", 0.0)
        self.safe_set("spec_power", 1.0)
        self.safe_set("spec_count", 0)
        self.safe_set("specular", 0.0)
        self.safe_set("shininess", 0.0)
        self.safe_set("ray_shadow", 1)
        self.safe_set("ray_trace_fog", 0)

        self.set_transparent_background()

        self.apply_illustrated_camera(field_of_view=25, depth_cue=0, fog=0.0)
        self.safe_set("antialias", 2)
        self.safe_set("two_sided_lighting", 1)
        self.safe_set("use_shaders", 1)

        self.finish_default(selection)


class StericSurfaceStyle(ScientificStyle):
    """Transparent steric surface style for catalyst pockets."""

    name = "steric_surface"
    command = "render_steric_surface"
    prefix = "steric"
    message = "Transparent steric surface style applied."

    def render(self, selection="all"):
        sel = self.selection_expr(selection)
        self.define_shared_colors()
        cmd.hide("everything", sel)

        atoms = self.select_coordination(selection)

        cmd.show("sticks", sel)
        cmd.hide("sticks", f"{sel} and elem H and not {atoms['hydride']}")
        cmd.hide("spheres", f"{sel} and elem H and not {atoms['hydride']}")
        sphere_atoms = self.coordination_sphere_atoms(atoms)
        if sphere_atoms:
            cmd.show("spheres", sphere_atoms)

        self.safe_set("stick_radius", 0.095, sel)
        self.safe_set("stick_radius", 0.14, atoms["coordination_core"])
        self.safe_set("sphere_scale", 0.52, atoms["metal"])
        self.safe_set("sphere_scale", 0.34, atoms["donors"])
        self.safe_set("sphere_scale", 0.34, atoms["co_c"])
        self.safe_set("sphere_scale", 0.34, atoms["co_o"])
        self.safe_set("sphere_scale", 0.12, atoms["hydride"])

        self.apply_coordination_sci_palette(sel, atoms)

        cmd.show("surface", sel)
        cmd.color("surface_sky", selection)
        self.safe_set("transparency", 0.68, sel)
        self.safe_set("surface_quality", 1)
        self.safe_set("surface_solvent", 1)

        # Re-apply atom colors after surface coloring paints the whole selection.
        self.apply_coordination_sci_palette(sel, atoms)

        self.apply_lighting(0.45, 0.22, 0.32, 0.74, 0.20, spec_power=120)
        self.apply_transparent_view(
            orthoscopic=1,
            field_of_view=30,
            depth_cue=1,
            fog_start=0.45,
        )
        self.finish_default(selection)


class QuasiChemDrawBoldStyle(ScientificStyle):
    """Flat ChemDraw-like bold ball-and-stick style."""

    name = "quasi_chemdraw_bold"
    command = "render_quasi_chemdraw_bold"
    prefix = "qcd"
    message = "Quasi-ChemDraw bold 3D style applied."
    colors = {
        "qcd_bond": [0.075, 0.080, 0.090],
        "qcd_carbon": [0.075, 0.080, 0.090],
        "qcd_hydrogen": [0.92, 0.92, 0.92],
        "qcd_metal": [0.95, 0.54, 0.10],
        "qcd_nitrogen": [0.10, 0.25, 0.88],
        "qcd_oxygen": [0.90, 0.04, 0.03],
        "qcd_sulfur": [0.96, 0.66, 0.05],
        "qcd_phosphor": [0.95, 0.42, 0.08],
        "qcd_halogen": [0.20, 0.64, 0.20],
    }

    def render(self, selection="all"):
        sel = f"({selection})"
        self.define_colors()
        atoms = self.select_coordination(selection)
        shell = self.first_shell(atoms)

        cmd.hide("everything", sel)
        cmd.show("sticks", sel)
        cmd.show("spheres", atoms["coordination_core"])
        cmd.hide("sticks", f"{sel} and elem H and not ({shell})")
        cmd.hide("spheres", f"{sel} and elem H and not ({shell})")

        self.safe_set("stick_radius", 0.155, sel)
        self.safe_set("stick_h_scale", 0.70)
        self.safe_set("stick_quality", 32)
        self.safe_set("stick_ball", 0)
        self.safe_set("smooth_half_bonds", 1)
        self.safe_set("stick_color", "qcd_bond", sel)
        try:
            cmd.set_bond("stick_radius", 0.185, atoms["metal"], shell)
        except Exception:
            pass

        self.safe_set("sphere_scale", 0.34, atoms["metal"])
        for category, scale in (
            ("N+O", 0.23),
            ("S+P", 0.18),
            ("H", 0.16),
            ("halogen", 0.17),
            ("C", 0.20),
        ):
            self._safe_set(
                "sphere_scale", scale, selection=shell, category=category
            )
        self.safe_set("sphere_transparency", 0.0, atoms["coordination_core"])
        self.safe_set("sphere_quality", 4)

        self.apply_style_palette(
            sel,
            {
                "C": "qcd_carbon",
                "H": "qcd_hydrogen",
                "N": "qcd_nitrogen",
                "O": "qcd_oxygen",
                "S": "qcd_sulfur",
                "P": "qcd_phosphor",
                "halogen": "qcd_halogen",
            },
            overrides={atoms["metal"]: "qcd_metal"},
        )

        self.safe_set("ambient", 0.58)
        self.safe_set("direct", 0.48)
        self.safe_set("reflect", 0.08)
        self.safe_set("specular", 0.20)
        self.safe_set("spec_reflect", 0.12)
        self.safe_set("spec_power", 45)
        self.safe_set("shininess", 25)
        self.safe_set("light_count", 8)
        self.safe_set("two_sided_lighting", 1)
        self.safe_set("ray_shadow", 0)
        self.safe_set("ambient_occlusion_mode", 0)
        self.safe_set("depth_cue", 0)
        self.safe_set("fog", 0.0)
        self.safe_set("orthoscopic", 1)
        self.safe_set("ray_orthoscopic", 1)
        self.safe_set("use_shaders", 1)
        self.safe_set("render_as_cylinders", 1)
        self.safe_set("antialias", 2)
        self.safe_set("ray_trace_antialias", 2)

        self.frame(selection, atoms["coordination_core"], zoom_buffer=1.4)
        self.finalize()
        print(self.message)


# ---------------------------------------------------------------------------
# Public PyMOL / CHEMSMART entry points
# ---------------------------------------------------------------------------


def _make_style_wrapper(style_cls):
    """Return a thin ``render_*`` wrapper for a :class:`ScientificStyle` subclass."""

    def wrapper(selection="all", **kwargs):
        kwargs.pop("_self", None)
        return style_cls().render(selection=selection)

    wrapper.__name__ = style_cls.command or style_cls.name
    wrapper.__doc__ = style_cls.__doc__
    return wrapper


metallic_poster_render = _make_style_wrapper(MetallicPosterStyle)
render_comic_metallic_labeled_final = _make_style_wrapper(ComicMetallicStyle)
render_soft_cartoon = _make_style_wrapper(SoftCartoonStyle)
render_editorial_minimal = _make_style_wrapper(EditorialMinimalStyle)
render_soft_ceramic = _make_style_wrapper(SoftCeramicStyle)
render_neon_coordination_core = _make_style_wrapper(NeonCoordinationCoreStyle)
render_matte_clay = _make_style_wrapper(MatteClayStyle)
render_xray_wire = _make_style_wrapper(XrayWireStyle)
render_steric_surface = _make_style_wrapper(StericSurfaceStyle)
render_quasi_chemdraw_bold = _make_style_wrapper(QuasiChemDrawBoldStyle)


_PYMOL_STYLE_COMMANDS = (
    ("metallic_poster_render", metallic_poster_render),
    (
        "render_comic_metallic_labeled_final",
        render_comic_metallic_labeled_final,
    ),
    ("render_soft_cartoon", render_soft_cartoon),
    ("render_editorial_minimal", render_editorial_minimal),
    ("render_soft_ceramic", render_soft_ceramic),
    ("render_neon_coordination_core", render_neon_coordination_core),
    ("render_matte_clay", render_matte_clay),
    ("render_xray_wire", render_xray_wire),
    ("render_steric_surface", render_steric_surface),
    ("render_quasi_chemdraw_bold", render_quasi_chemdraw_bold),
)

for _command_name, _command_fn in _PYMOL_STYLE_COMMANDS:
    cmd.extend(_command_name, _command_fn)
