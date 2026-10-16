"""
Direct unit tests for chemsmart.io.xtb.input.XTBInput's numeric-parsing
fallback branches.

XTBInput parses ``$group`` / ``key=value`` sections of an xTB input
file. Nearly every numeric property follows the same pattern:

    try:
        return float(self._get_key(group, key, default))
    except ValueError:
        return default

The "valid value" path for these properties is already covered end to
end against real xTB input files in tests/test_XTBIO.py::TestXTBInput.
This file instead targets the "malformed value falls back to the
documented default" branch for every such property, which real input
files never naturally exercise. It bypasses the file-existence check
in __init__ by constructing a bare instance and injecting a fake
content_groups dict directly (content_groups is a functools.cached_property,
so assigning into the instance __dict__ pre-empties recomputation from
a real file).
"""

import pytest

from chemsmart.io.xtb.input import XTBInput

# group -> {key: "not-a-number"} for every property whose ValueError
# fallback we want to exercise.
_MALFORMED_GROUPS = {
    "chrg": ["notanumber"],
    "spin": ["notanumber"],
    "cube": [
        "step=notanumber",
        "pthr=notanumber",
        "boff=notanumber",
        "cal=notanumber",
    ],
    "embedding": ["at=notanumber"],
    "gfn": ["dispscale=notanumber"],
    "hess": [
        "sccacc=notanumber",
        "step=notanumber",
        "scale=notanumber",
    ],
    "md": [
        "temp=notanumber",
        "time=notanumber",
        "dump=notanumber",
        "skip=notanumber",
        "step=notanumber",
        "hmass=notanumber",
        "shake=notanumber",
        "sccacc=notanumber",
    ],
    "modef": [
        "n=notanumber",
        "step=notanumber",
        "updat=notanumber",
        "local=notanumber",
        "vthr=notanumber",
        "prj=notanumber",
        "mode=notanumber",
    ],
    "opt": [
        "optlevel=notanumber",
        "microcycle=notanumber",
        "maxcycle=notanumber",
        "maxdispl=notanumber",
        "hlow=notanumber",
        "s6=notanumber",
        "kstretch=notanumber",
        "kbend=notanumber",
        "ktorsion=notanumber",
        "koutofp=notanumber",
        "kvdw=notanumber",
        "kes=notanumber",
        "rcut=notanumber",
    ],
    "path": [
        "nrun=notanumber",
        "npoint=notanumber",
        "anopt=notanumber",
        "kpush=notanumber",
        "kpull=notanumber",
        "alp=notanumber",
    ],
    "scc": [
        "maxiterations=notanumber",
        "temp=notanumber",
        "broydamp=notanumber",
    ],
    "symmetry": ["desy=notanumber", "maxat=notanumber"],
    "thermo": [
        "temp=notanumber",
        "sthr=notanumber",
        "imagthr=notanumber",
        "scale=notanumber",
    ],
    "wall": [
        "alpha=notanumber",
        "beta=notanumber",
        "temp=notanumber",
        "autoscale=notanumber",
        "axisshift=notanumber",
    ],
}


def _make_malformed_input():
    """Bare XTBInput with every numeric key set to a non-numeric
    string, to force every ValueError-fallback branch at once."""
    inp = XTBInput.__new__(XTBInput)
    inp.__dict__["content_groups"] = _MALFORMED_GROUPS
    return inp


class TestXTBInputFileNotFound:
    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError, match="not found"):
            XTBInput(filename="/no/such/xtb/input/file.inp")


class TestXTBInputGetKey:
    def test_group_not_present_returns_default(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {}
        assert inp._get_key("cube", "step", "0.4") == "0.4"

    def test_key_not_present_in_group_returns_default(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {"cube": ["boff=3.0"]}
        assert inp._get_key("cube", "step", "0.4") == "0.4"


class TestXTBInputContentGroups:
    def test_lines_before_any_group_marker_are_ignored(self):
        """A stray non-'$'-prefixed line before any group has been
        opened has no current_group to append to."""
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["contents"] = [
            "stray line with no group yet",
            "$cube",
            "step=0.4",
        ]
        assert inp.content_groups == {"cube": ["step=0.4"]}


class TestXTBInputChargeAndSpin:
    def test_charge_falls_back_to_zero_on_malformed_value(self):
        assert _make_malformed_input().charge == 0

    def test_spin_falls_back_to_zero_on_malformed_value(self):
        assert _make_malformed_input().spin == 0

    def test_charge_default_when_group_absent(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {}
        assert inp.charge == 0

    def test_spin_default_when_group_absent(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {}
        assert inp.spin == 0


class TestXTBInputNoneOnInvalidRatherThanNumericDefault:
    """atom_type, modef_n, and mode_following return None (rather
    than falling back to a numeric default) when the key is present
    but not a valid integer."""

    def test_atom_type_none_on_malformed_value(self):
        assert _make_malformed_input().atom_type is None

    def test_atom_type_none_when_key_absent(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {}
        assert inp.atom_type is None

    def test_modef_n_none_on_malformed_value(self):
        assert _make_malformed_input().modef_n is None

    def test_modef_n_none_when_key_absent(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {}
        assert inp.modef_n is None

    def test_mode_following_none_on_malformed_value(self):
        assert _make_malformed_input().mode_following is None

    def test_mode_following_none_when_key_absent(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {}
        assert inp.mode_following is None


class TestXTBInputOptimizationLevelFallback:
    def test_unrecognized_string_and_non_integer_falls_back_to_normal(self):
        assert _make_malformed_input().optimization_level == "normal"

    @pytest.mark.parametrize(
        "int_level,expected",
        [
            ("-3", "crude"),
            ("-2", "sloppy"),
            ("-1", "loose"),
            ("0", "normal"),
            ("1", "tight"),
            ("2", "verytight"),
            ("3", "extreme"),
        ],
    )
    def test_integer_optlevel_maps_to_keyword(self, int_level, expected):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {"opt": [f"optlevel={int_level}"]}
        assert inp.optimization_level == expected

    def test_out_of_range_integer_falls_back_to_normal(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {"opt": ["optlevel=99"]}
        assert inp.optimization_level == "normal"


class TestXTBInputBooleanKeywordProperties:
    """isotropic_electrostatic/scc/periodic/exact_rational_function/
    average_convergence all compare `_get_key(...).lower() == "true"`,
    so any non-"true" string (not just a literal "false") yields
    False -- confirmed here rather than assumed."""

    def test_non_true_string_is_falsy(self):
        inp = XTBInput.__new__(XTBInput)
        inp.__dict__["content_groups"] = {
            "embedding": ["es=nonsense"],
            "gfn": ["scc=nonsense", "periodic=nonsense"],
            "opt": ["exact rf=nonsense", "average conv=nonsense"],
        }
        assert inp.isotropic_electrostatic is False
        assert inp.scc is False
        assert inp.periodic is False
        assert inp.exact_rational_function is False
        assert inp.average_convergence is False


@pytest.mark.parametrize(
    "attr,expected_default",
    [
        ("cube_step", 0.4),
        ("density_matrix_threshold", 0.05),
        ("boundary_offset", 3.0),
        ("cube_output", 1),
        ("dispersion_energy_scale", 1.0),
        ("hess_scc_accuracy", 0.3),
        ("hess_step", 0.005),
        ("hess_scale", 1.0),
        ("md_temperature", 298.15),
        ("md_time", 50.0),
        ("dump_structure", 50.0),
        ("skip_interval", 500),
        ("md_step", 4.0),
        ("hydrogen_mass", 4),
        ("shake_algorithm", 2),
        ("md_scc_accuracy", 2.0),
        ("modef_step", 1.0),
        ("modef_update", 0.2),
        ("modef_local", 0),
        ("modef_threshold", 0.0),
        ("projected_mode", 0),
        ("anc_microcycles", 20),
        ("max_optcycles", 0),
        ("max_displacement", 1.0),
        ("low_frequency_cutoff", 0.01),
        ("s6_in_model_hessian", 20.0),
        ("stretch_force_constant", 0.4),
        ("bend_force_constant", 0.13),
        ("torsion_force_constant", 0.0075),
        ("out_of_plane_force_constant", 0.0),
        ("additional_vdw_contribution", 0.0),
        ("electrostatic_contribution", 0.0),
        ("distance_cutoff", 8.366600265340756),
        ("pathfinder_runs", 3),
        ("path_points", 50),
        ("path_optimization_steps", 3),
        ("rmsd_push_factor", 0.05),
        ("rmsd_pull_factor", -0.04),
        ("rmsd_width", 0.7),
        ("max_iterations", 250),
        ("electronic_temperature", 300.0),
        ("broyden_damping", 0.4),
        ("symmetry_threshold", 0.1),
        ("symmetry_max_atoms", 200),
        ("thermo_temperature", 298.15),
        ("rotor_cutoff", 50.0),
        ("imaginary_frequency_cutoff", -20.0),
        ("scaling_factor", 1.0),
        ("wall_potential_exponent", 30),
        ("logfermi_bias_exponent", 6.0),
        ("wall_temperature", 300.0),
        ("auto_scale", 1.0),
        ("axis_shift", 3.5),
    ],
)
def test_numeric_property_falls_back_to_default_on_malformed_value(
    attr, expected_default
):
    inp = _make_malformed_input()
    assert getattr(inp, attr) == expected_default
