"""
Direct unit tests for chemsmart.io.xtb.route.XTBRoute.

XTBRoute is a pure string parser (no I/O), so every test simply
constructs it with a synthetic xTB command-line route string.
"""

import pytest

from chemsmart.io.xtb.route import XTBRoute


class TestXTBRouteInit:
    def test_lowercases_and_splits_route_string(self):
        route = XTBRoute("XTB Coord.xyz --GFN 2 --Opt Tight")
        assert route.route_string == "xtb coord.xyz --gfn 2 --opt tight"
        assert route.route_inputs == [
            "xtb",
            "coord.xyz",
            "--gfn",
            "2",
            "--opt",
            "tight",
        ]


class TestXTBRouteGfnVersionAndMethod:
    @pytest.mark.parametrize(
        "flag,expected",
        [
            ("--gfn0", "gfn0"),
            ("--gfn1", "gfn1"),
            ("--gfn2", "gfn2"),
            ("--gfnff", "gfnff"),
        ],
    )
    def test_explicit_gfn_flags(self, flag, expected):
        route = XTBRoute(f"xtb coord.xyz {flag}")
        assert route.gfn_version == expected
        assert route.method == expected

    def test_gff_alias_for_gfnff(self):
        route = XTBRoute("xtb coord.xyz --gff")
        assert route.gfn_version == "gfnff"

    @pytest.mark.parametrize(
        "num,expected",
        [("0", "gfn0"), ("1", "gfn1"), ("2", "gfn2")],
    )
    def test_gfn_int_format(self, num, expected):
        route = XTBRoute(f"xtb coord.xyz --gfn {num}")
        assert route.gfn_version == expected

    def test_gfn_int_out_of_supported_range_returns_none(self):
        route = XTBRoute("xtb coord.xyz --gfn 3")
        assert route.gfn_version is None

    def test_gfn_non_integer_value_returns_none(self):
        route = XTBRoute("xtb coord.xyz --gfn abc")
        assert route.gfn_version is None

    def test_gfn_flag_at_end_of_route_returns_none(self):
        route = XTBRoute("xtb coord.xyz --gfn")
        assert route.gfn_version is None

    def test_no_method_specified_returns_none(self):
        route = XTBRoute("xtb coord.xyz")
        assert route.gfn_version is None
        assert route.method is None


class TestXTBRouteBasis:
    def test_basis_always_default(self):
        assert XTBRoute("xtb coord.xyz").basis == "default"


class TestXTBRouteOptimizationLevel:
    @pytest.mark.parametrize(
        "flag", ["--opt", "-o", "--optlevel", "--metaopt"]
    )
    def test_each_flag_reads_following_level(self, flag):
        route = XTBRoute(f"xtb coord.xyz {flag} tight")
        assert route.optimization_level == "tight"

    def test_verytight_normalized_to_vtight(self):
        route = XTBRoute("xtb coord.xyz --opt verytight")
        assert route.optimization_level == "vtight"

    def test_unknown_level_returns_none(self):
        route = XTBRoute("xtb coord.xyz --opt bogus_level")
        assert route.optimization_level is None

    def test_flag_at_end_of_route_returns_none(self):
        route = XTBRoute("xtb coord.xyz --opt")
        assert route.optimization_level is None

    def test_no_opt_flag_returns_none(self):
        route = XTBRoute("xtb coord.xyz")
        assert route.optimization_level is None


class TestXTBRouteChargeAndUhf:
    @pytest.mark.parametrize("flag", ["--chrg", "-c"])
    def test_charge_flags(self, flag):
        route = XTBRoute(f"xtb coord.xyz {flag} -1")
        assert route.charge == -1

    def test_charge_invalid_value_returns_none(self):
        route = XTBRoute("xtb coord.xyz --chrg notanumber")
        assert route.charge is None

    def test_charge_flag_at_end_returns_none(self):
        route = XTBRoute("xtb coord.xyz --chrg")
        assert route.charge is None

    def test_no_charge_flag_returns_none(self):
        assert XTBRoute("xtb coord.xyz").charge is None

    @pytest.mark.parametrize("flag", ["--uhf", "-u"])
    def test_uhf_flags(self, flag):
        route = XTBRoute(f"xtb coord.xyz {flag} 2")
        assert route.uhf == 2

    def test_uhf_invalid_value_returns_none(self):
        route = XTBRoute("xtb coord.xyz --uhf notanumber")
        assert route.uhf is None

    def test_uhf_flag_at_end_returns_none(self):
        route = XTBRoute("xtb coord.xyz --uhf")
        assert route.uhf is None

    def test_no_uhf_flag_returns_none(self):
        assert XTBRoute("xtb coord.xyz").uhf is None


class TestXTBRouteJobtype:
    @pytest.mark.parametrize(
        "flag",
        ["--opt", "-o", "--omd", "--metaopt", "--ohess", "--bhess"],
    )
    def test_opt_keywords(self, flag):
        route = XTBRoute(f"xtb coord.xyz {flag}")
        assert route.jobtype == "opt"
        assert route.get_jobtype() == "opt"

    def test_hess_keyword(self):
        assert XTBRoute("xtb coord.xyz --hess").jobtype == "hess"

    @pytest.mark.parametrize("flag", ["--md", "--metadyn"])
    def test_md_keywords(self, flag):
        assert XTBRoute(f"xtb coord.xyz {flag}").jobtype == "md"

    def test_path_keyword(self):
        assert XTBRoute("xtb coord.xyz --path").jobtype == "path"

    def test_modef_keyword(self):
        assert XTBRoute("xtb coord.xyz --modef").jobtype == "modef"

    def test_default_is_sp(self):
        assert XTBRoute("xtb coord.xyz").jobtype == "sp"

    def test_opt_takes_precedence_over_hess(self):
        """--ohess/--bhess do an optimization first, so they resolve to
        'opt' even though they also imply a Hessian calculation."""
        assert XTBRoute("xtb coord.xyz --ohess").jobtype == "opt"
        assert XTBRoute("xtb coord.xyz --bhess").jobtype == "opt"


class TestXTBRouteFreqAndGrad:
    @pytest.mark.parametrize("flag", ["--hess", "--ohess", "--bhess"])
    def test_freq_true_for_hessian_flags(self, flag):
        assert XTBRoute(f"xtb coord.xyz {flag}").freq is True

    def test_freq_false_without_hessian_flags(self):
        assert XTBRoute("xtb coord.xyz --opt").freq is False

    def test_grad_true(self):
        assert XTBRoute("xtb coord.xyz --grad").grad is True

    def test_grad_false(self):
        assert XTBRoute("xtb coord.xyz").grad is False


class TestXTBRouteSolventModelAndId:
    @pytest.mark.parametrize(
        "model", ["gbsa", "alpb", "cosmo", "tmcosmo", "cpcmx"]
    )
    def test_each_solvent_model_flag(self, model):
        route = XTBRoute(f"xtb coord.xyz --{model} water")
        assert route.solvent_model == model
        assert route.get_solvent_model() == model

    def test_g_alias_for_gbsa(self):
        route = XTBRoute("xtb coord.xyz -g water")
        assert route.solvent_model == "gbsa"

    def test_no_solvent_model_returns_none(self):
        assert XTBRoute("xtb coord.xyz").solvent_model is None

    def test_solvent_id_extracted_when_model_present(self):
        route = XTBRoute("xtb coord.xyz --alpb toluene")
        assert route.solvent_id == "toluene"

    def test_solvent_id_none_when_no_solvent_model(self):
        route = XTBRoute("xtb coord.xyz toluene")
        assert route.solvent_id is None

    def test_solvent_id_none_when_model_present_but_no_matching_token(self):
        route = XTBRoute("xtb coord.xyz --alpb")
        assert route.solvent_id is None


class TestXTBRouteAccuracyAndElectronicTemperature:
    @pytest.mark.parametrize("flag", ["--acc", "-a"])
    def test_accuracy_flags(self, flag):
        route = XTBRoute(f"xtb coord.xyz {flag} 0.5")
        assert route.accuracy == pytest.approx(0.5)

    def test_accuracy_invalid_value_returns_none(self):
        route = XTBRoute("xtb coord.xyz --acc notanumber")
        assert route.accuracy is None

    def test_accuracy_flag_at_end_returns_none(self):
        assert XTBRoute("xtb coord.xyz --acc").accuracy is None

    def test_no_accuracy_flag_returns_none(self):
        assert XTBRoute("xtb coord.xyz").accuracy is None

    def test_electronic_temperature(self):
        route = XTBRoute("xtb coord.xyz --etemp 500")
        assert route.electronic_temperature == pytest.approx(500.0)

    def test_electronic_temperature_invalid_value_returns_none(self):
        route = XTBRoute("xtb coord.xyz --etemp notanumber")
        assert route.electronic_temperature is None

    def test_electronic_temperature_flag_at_end_returns_none(self):
        assert XTBRoute("xtb coord.xyz --etemp").electronic_temperature is None

    def test_no_electronic_temperature_flag_returns_none(self):
        assert XTBRoute("xtb coord.xyz").electronic_temperature is None


class TestXTBRouteBooleanFlags:
    """Every simple "--flag in route_inputs" boolean property."""

    @pytest.mark.parametrize(
        "attr,flag",
        [
            ("ptb", "--ptb"),
            ("spinpol", "--spinpol"),
            ("ceh", "--ceh"),
            ("pop", "--pop"),
            ("wbo", "--wbo"),
            ("dipole", "--dipole"),
            ("molden", "--molden"),
            ("lmo", "--lmo"),
            ("fod", "--fod"),
            ("esp", "--esp"),
            ("stm", "--stm"),
            ("vip", "--vip"),
            ("vea", "--vea"),
            ("vipea", "--vipea"),
            ("vfukui", "--vfukui"),
            ("vomega", "--vomega"),
            ("alpha", "--alpha"),
            ("cma", "--cma"),
        ],
    )
    def test_flag_true_when_present(self, attr, flag):
        route = XTBRoute(f"xtb coord.xyz {flag}")
        assert getattr(route, attr) is True

    @pytest.mark.parametrize(
        "attr",
        [
            "ptb",
            "spinpol",
            "ceh",
            "pop",
            "wbo",
            "dipole",
            "molden",
            "lmo",
            "fod",
            "esp",
            "stm",
            "vip",
            "vea",
            "vipea",
            "vfukui",
            "vomega",
            "alpha",
            "cma",
        ],
    )
    def test_flag_false_when_absent(self, attr):
        route = XTBRoute("xtb coord.xyz")
        assert getattr(route, attr) is False
