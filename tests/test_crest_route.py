"""
Direct unit tests for chemsmart.io.crest.route.CRESTRoute.

CRESTRoute is a pure string parser (no I/O), so every test simply
constructs it with a synthetic CREST command-line route string.
"""

import pytest

from chemsmart.io.crest.route import CRESTRoute


class TestCRESTRouteInit:
    def test_lowercases_and_splits_route_string(self):
        route = CRESTRoute("CREST Mol.xyz --GFN2 --Chrg 0")
        assert route.route_string == "crest mol.xyz --gfn2 --chrg 0"
        assert route.route_inputs == [
            "crest",
            "mol.xyz",
            "--gfn2",
            "--chrg",
            "0",
        ]


class TestCRESTRouteGfnVersionAndMethod:
    @pytest.mark.parametrize(
        "flag,expected",
        [
            ("--gfn1", "gfn1"),
            ("--gfn2", "gfn2"),
            ("--gfnff", "gfnff"),
            ("--gfn2//gfnff", "gfn2//gfnff"),
        ],
    )
    def test_explicit_gfn_flags(self, flag, expected):
        route = CRESTRoute(f"crest mol.xyz {flag}")
        assert route.gfn_version == expected
        assert route.method == expected

    def test_gff_alias_for_gfnff(self):
        route = CRESTRoute("crest mol.xyz --gff")
        assert route.gfn_version == "gfnff"

    @pytest.mark.parametrize(
        "num,expected",
        [("0", "gfn0"), ("1", "gfn1"), ("2", "gfn2")],
    )
    def test_gfn_int_format(self, num, expected):
        route = CRESTRoute(f"crest mol.xyz --gfn {num}")
        assert route.gfn_version == expected

    def test_gfn_int_out_of_supported_range_returns_none(self):
        route = CRESTRoute("crest mol.xyz --gfn 3")
        assert route.gfn_version is None

    def test_gfn_non_integer_value_returns_none(self):
        route = CRESTRoute("crest mol.xyz --gfn abc")
        assert route.gfn_version is None

    def test_gfn_flag_at_end_of_route_returns_none(self):
        route = CRESTRoute("crest mol.xyz --gfn")
        assert route.gfn_version is None

    def test_no_method_specified_returns_none(self):
        route = CRESTRoute("crest mol.xyz")
        assert route.gfn_version is None
        assert route.method is None


class TestCRESTRouteBasis:
    def test_basis_always_default(self):
        assert CRESTRoute("crest mol.xyz").basis == "default"


class TestCRESTRouteOptimizationLevel:
    @pytest.mark.parametrize(
        "flag",
        ["--optlev", "-optlev", "--opt", "-o", "--optlevel", "--ancopt"],
    )
    def test_each_flag_reads_following_level(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} tight")
        assert route.optimization_level == "tight"

    def test_verytight_normalized_to_vtight(self):
        route = CRESTRoute("crest mol.xyz --opt verytight")
        assert route.optimization_level == "vtight"

    def test_unknown_level_returns_none(self):
        route = CRESTRoute("crest mol.xyz --opt bogus_level")
        assert route.optimization_level is None

    def test_flag_at_end_of_route_returns_none(self):
        route = CRESTRoute("crest mol.xyz --opt")
        assert route.optimization_level is None

    def test_no_opt_flag_returns_none(self):
        assert CRESTRoute("crest mol.xyz").optimization_level is None


class TestCRESTRouteChargeAndUhf:
    @pytest.mark.parametrize("flag", ["--chrg", "-chrg"])
    def test_charge_flags(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} -1")
        assert route.charge == -1

    def test_charge_invalid_value_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --chrg notanumber")
        assert route.charge == 0

    def test_charge_flag_at_end_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --chrg")
        assert route.charge == 0

    def test_no_charge_flag_defaults_to_zero(self):
        assert CRESTRoute("crest mol.xyz").charge == 0

    @pytest.mark.parametrize("flag", ["--uhf", "-uhf"])
    def test_uhf_flags(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} 2")
        assert route.uhf == 2
        assert route.multiplicity == 3

    def test_uhf_invalid_value_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --uhf notanumber")
        assert route.uhf == 0

    def test_uhf_flag_at_end_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --uhf")
        assert route.uhf == 0

    def test_no_uhf_flag_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz")
        assert route.uhf == 0
        assert route.multiplicity == 1


class TestCRESTRouteJobtype:
    def test_jobtype_is_always_conformers(self):
        assert CRESTRoute("crest mol.xyz").jobtype == "conformers"
        assert CRESTRoute("crest mol.xyz --gfn2").get_jobtype() == "conformers"


class TestCRESTRouteSolventModelAndId:
    @pytest.mark.parametrize("model", ["gbsa", "alpb"])
    def test_each_solvent_model_flag(self, model):
        route = CRESTRoute(f"crest mol.xyz --{model} water")
        assert route.solvent_model == model
        assert route.get_solvent_model() == model

    def test_g_alias_for_gbsa(self):
        route = CRESTRoute("crest mol.xyz -g water")
        assert route.solvent_model == "gbsa"

    def test_no_solvent_model_returns_none(self):
        assert CRESTRoute("crest mol.xyz").solvent_model is None

    def test_solvent_id_extracted_when_model_present(self):
        route = CRESTRoute("crest mol.xyz --alpb toluene")
        assert route.solvent_id == "toluene"

    def test_solvent_id_none_when_no_solvent_model(self):
        route = CRESTRoute("crest mol.xyz toluene")
        assert route.solvent_id is None

    def test_solvent_id_none_when_model_present_but_no_matching_token(self):
        route = CRESTRoute("crest mol.xyz --alpb")
        assert route.solvent_id is None


class TestCRESTRouteBooleanFlags:
    @pytest.mark.parametrize(
        "attr,long_flag,short_flag",
        [
            ("nci", "--nci", "-nci"),
            ("constrained", "--cinp", "-cinp"),
            (
                "no_reference_topology_check",
                "--noreftopo",
                "-noreftopo",
            ),
            ("no_topology_check", "--notopo", "-notopo"),
        ],
    )
    def test_flag_true_for_long_and_short_forms(
        self, attr, long_flag, short_flag
    ):
        assert getattr(CRESTRoute(f"crest mol.xyz {long_flag}"), attr) is True
        assert getattr(CRESTRoute(f"crest mol.xyz {short_flag}"), attr) is True

    @pytest.mark.parametrize(
        "attr",
        [
            "nci",
            "constrained",
            "no_reference_topology_check",
            "no_topology_check",
        ],
    )
    def test_flag_false_when_absent(self, attr):
        assert getattr(CRESTRoute("crest mol.xyz"), attr) is False


class TestCRESTRouteThreads:
    @pytest.mark.parametrize("flag", ["-T", "--T"])
    def test_threads_flags_case_insensitive(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} 8")
        assert route.threads == 8

    def test_threads_invalid_value_returns_none(self):
        route = CRESTRoute("crest mol.xyz -T notanumber")
        assert route.threads is None

    def test_threads_flag_at_end_returns_none(self):
        route = CRESTRoute("crest mol.xyz -T")
        assert route.threads is None

    def test_no_threads_flag_returns_none(self):
        assert CRESTRoute("crest mol.xyz").threads is None


class TestCRESTRouteQuickMode:
    @pytest.mark.parametrize("mode", ["quick", "squick", "mquick"])
    @pytest.mark.parametrize("prefix", ["--", "-"])
    def test_quick_mode_flags(self, mode, prefix):
        route = CRESTRoute(f"crest mol.xyz {prefix}{mode}")
        assert route.quick_mode == mode

    def test_no_quick_mode_returns_none(self):
        assert CRESTRoute("crest mol.xyz").quick_mode is None


class TestCRESTRouteNumericThresholds:
    """Every 'flag REAL/INT' property sharing the same
    present/invalid/flag-at-end/absent shape."""

    _CASES = [
        ("rmsd_threshold", "--rthr", "-rthr", "0.5", 0.5),
        ("energy_threshold", "--ethr", "-ethr", "0.1", 0.1),
        ("bconst_threshold", "--bthr", "-bthr", "0.02", 0.02),
        ("population_threshold", "--pthr", "-pthr", "0.05", 0.05),
        ("temperature", "--temp", "-temp", "298.15", 298.15),
        ("shake", "--shake", "-shake", "1", 1),
        ("md_timestep", "--tstep", "-tstep", "4", 4),
        ("md_dump_step", "--mddump", "-mddump", "50", 50),
        ("vbias_dump_interval", "--vbdump", "-vbdump", "1.0", 1.0),
        ("additional_md_temperature", "--tnmd", "-tnmd", "400", 400.0),
        ("energy_window", "--ewin", "-ewin", "6.0", 6.0),
    ]

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_long_flag_parses_value(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {long_flag} {raw}")
        assert getattr(route, attr) == pytest.approx(expected)

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_short_flag_parses_value(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {short_flag} {raw}")
        assert getattr(route, attr) == pytest.approx(expected)

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_invalid_value_returns_none(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {long_flag} notanumber")
        assert getattr(route, attr) is None

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_flag_at_end_returns_none(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {long_flag}")
        assert getattr(route, attr) is None

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_absent_returns_none(
        self, attr, long_flag, short_flag, raw, expected
    ):
        assert getattr(CRESTRoute("crest mol.xyz"), attr) is None


class TestCRESTRouteMdLength:
    @pytest.mark.parametrize("flag", ["--mdlen", "-mdlen", "--len", "-len"])
    def test_numeric_value_parsed_as_float(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} 50")
        assert route.md_length == pytest.approx(50.0)

    def test_x_prefixed_multiplicative_factor_returned_as_string(self):
        route = CRESTRoute("crest mol.xyz --mdlen x2.0")
        assert route.md_length == "x2.0"

    def test_non_numeric_non_x_value_returned_as_string(self):
        route = CRESTRoute("crest mol.xyz --mdlen notanumber")
        assert route.md_length == "notanumber"

    def test_flag_at_end_returns_none(self):
        assert CRESTRoute("crest mol.xyz --mdlen").md_length is None

    def test_absent_returns_none(self):
        assert CRESTRoute("crest mol.xyz").md_length is None
