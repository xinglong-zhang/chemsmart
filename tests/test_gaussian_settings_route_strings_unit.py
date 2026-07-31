"""
Direct unit tests for route-string generation and settings-building logic
in chemsmart.jobs.gaussian.settings that lacks direct test coverage:

- GaussianIRCJobSettings._get_route_string_from_jobtype
- GaussianLinkJobSettings.link_route_string /
  _get_route_string_from_jobtype / _get_link_route_string_from_jobtype
- GaussianpKaJobSettings.build_gaussian_pka_settings
"""

import pytest

from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianpKaJobSettings,
)


class TestIRCRouteString:
    def test_jobtype_neither_ircf_nor_ircr_leaves_direction_unset(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="irc",
            predictor="LQA",
            recorrect="never",
            direction="forward",
        )
        settings._get_route_string_from_jobtype()
        # direction is untouched by the ircf/ircr-specific branches
        assert settings.direction == "forward"

    def test_predictor_and_recorrect_both_specified(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
            recorrect="never",
        )
        route = settings._get_route_string_from_jobtype()
        assert "irc(LQA,calcfc,recorrect=never" in route
        assert "forward" in route
        assert settings.direction == "forward"

    def test_ircr_sets_reverse_direction(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircr",
            predictor="LQA",
            recorrect="never",
        )
        route = settings._get_route_string_from_jobtype()
        assert "reverse" in route
        assert settings.direction == "reverse"

    def test_neither_predictor_nor_recorrect_uses_basic_route(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="ircf"
        )
        route = settings._get_route_string_from_jobtype()
        assert "irc(calcfc,recalc=" in route
        assert "LQA" not in route

    def test_only_predictor_specified_raises(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
        )
        with pytest.raises(ValueError, match="Only one of predictor"):
            settings._get_route_string_from_jobtype()

    def test_only_recorrect_specified_raises(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            recorrect="never",
        )
        with pytest.raises(ValueError, match="Only one of predictor"):
            settings._get_route_string_from_jobtype()

    def test_flat_irc_applies_defaults(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            flat_irc=True,
        )
        route = settings._get_route_string_from_jobtype()
        assert settings.predictor == "LQA"
        assert settings.recorrect == "never"
        assert settings.recalc_step == -5
        assert "irc(LQA,calcfc,recorrect=never" in route

    def test_flat_irc_does_not_override_explicit_recalc_step(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            flat_irc=True,
            recalc_step=3,
        )
        settings._get_route_string_from_jobtype()
        assert settings.recalc_step == 3

    def test_additional_route_parameters_appended_once(self):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            additional_route_parameters="scf=qc",
        )
        route = settings._get_route_string_from_jobtype()
        assert route.count("scf=qc") == 1

    def test_additional_route_parameters_not_duplicated_if_present(self):
        # The base class (via super()._get_route_string_from_jobtype())
        # unconditionally appends additional_route_parameters once
        # already; IRC's own duplicate-check at the end of the method
        # then correctly recognizes it's already present and skips
        # re-appending it, so it still only shows up once overall.
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            additional_route_parameters="iop(1/8=3)",
        )
        route = settings._get_route_string_from_jobtype()
        assert route.count("iop(1/8=3)") == 1


class TestLinkRouteString:
    def test_link_route_string_uses_custom_link_route(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            link_route="#p opt",
        )
        result = settings.link_route_string
        assert "b3lyp" in result
        assert "sto-3g" in result
        assert "geom=check" in result
        assert "guess=read" in result

    def test_link_route_string_skips_already_present_tokens(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            link_route="#p opt b3lyp sto-3g geom=check guess=read",
        )
        result = settings.link_route_string
        assert result.count("b3lyp") == 1
        assert result.count("sto-3g") == 1
        assert result.count("geom=check") == 1
        assert result.count("guess=read") == 1

    def test_link_route_string_falls_back_to_jobtype_logic(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="opt"
        )
        result = settings.link_route_string
        assert "geom=check" in result
        assert "guess=read" in result

    def test_stability_route_removes_opt_and_freq_keywords(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            freq=True,
            stable="opt",
            guess="mix",
        )
        route = settings._get_route_string_from_jobtype()
        assert "opt" not in route.replace("stable=opt", "")
        assert "freq" not in route
        assert "stable=opt" in route
        assert "guess=mix" in route

    def test_guess_with_comma_gets_parenthesized(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            guess="mix,always",
        )
        route = settings._get_route_string_from_jobtype()
        assert "guess=(mix,always)" in route

    def test_guess_with_surrounding_parens_normalized(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            guess="(mix,always)",
        )
        route = settings._get_route_string_from_jobtype()
        assert "guess=(mix,always)" in route
        assert "guess=((mix,always))" not in route

    def test_no_stable_or_guess_omits_both(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="opt",
            stable=None,
            guess=None,
        )
        route = settings._get_route_string_from_jobtype()
        assert "stable=" not in route
        assert "guess=" not in route

    def test_link_route_for_non_irc_jobtype(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="opt"
        )
        route = settings._get_link_route_string_from_jobtype()
        assert "geom=check" in route
        assert "guess=read" in route

    def test_link_route_for_irc_jobtype_delegates_to_irc_settings(self):
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
            recorrect="never",
        )
        route = settings._get_link_route_string_from_jobtype()
        assert "irc(LQA" in route
        assert "geom=check" in route
        assert "guess=read" in route

    def test_link_route_for_irc_does_not_duplicate_existing_tokens(self):
        # If the IRC-generated route already ends up containing
        # "geom=check"/"guess=read" (unlikely in practice, but the guard
        # exists), the method must not duplicate them. We simulate this
        # by checking the guard directly via the additional_route_parameters
        # channel, which is appended by the underlying IRC settings.
        settings = GaussianLinkJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ircf",
            predictor="LQA",
            recorrect="never",
            additional_route_parameters="geom=check guess=read",
        )
        route = settings._get_link_route_string_from_jobtype()
        assert route.count("geom=check") == 1
        assert route.count("guess=read") == 1


class TestBuildGaussianPkaSettings:
    def test_builds_settings_with_defaults(self):
        opt_settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", charge=0, multiplicity=1
        )
        shared = {
            "reference": "ref.log",
            "delta_g_proton": -270.3,
            "reference_color_code": "red",
            "skip_completed": True,
        }
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=1, shared=shared, opt_settings=opt_settings
        )
        assert settings.proton_index == 1
        assert settings.reference_file == "ref.log"
        assert settings.delta_G_proton == -270.3
        assert settings.functional == "b3lyp"
        assert settings.basis == "sto-3g"
        # cli_only keys must not leak into the constructed settings kwargs
        assert not hasattr(settings, "reference_color_code")
        # Defaults applied since neither pka_kwargs nor opt/sp settings
        # specify solvent info.
        assert settings.solvent_model == "SMD"
        assert settings.solvent_id == "water"

    def test_solvent_settings_prefer_explicit_shared_values(self):
        opt_settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            solvent_model="PCM",
            solvent_id="dmso",
        )
        shared = {"solvent_model": "SMD", "solvent_id": "acetonitrile"}
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=2, shared=shared, opt_settings=opt_settings
        )
        assert settings.solvent_model == "SMD"
        assert settings.solvent_id == "acetonitrile"

    def test_solvent_settings_from_opt_settings_crashes(self):
        # Documents a real bug: when opt_settings carries a non-None
        # solvent_model/solvent_id and `shared` doesn't override it, the
        # resolved value ends up in both `pka_kwargs` (via the
        # solvent_model/solvent_id fallback resolution) and `opt_kwargs`
        # (via the vars(opt_settings) passthrough, since at the time
        # opt_kwargs is built pka_kwargs doesn't contain solvent_model/
        # solvent_id yet). Passing both dicts as **kwargs to the same
        # call then raises "got multiple values for keyword argument".
        # See BUGS_FOUND.md.
        opt_settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            solvent_model="PCM",
            solvent_id="dmso",
        )
        with pytest.raises(TypeError, match="multiple values"):
            GaussianpKaJobSettings.build_gaussian_pka_settings(
                proton_index=2, shared={}, opt_settings=opt_settings
            )

    def test_solvent_settings_fall_back_to_sp_settings(self):
        opt_settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        sp_settings = GaussianJobSettings(
            functional="wb97xd",
            basis="def2-tzvp",
            solvent_model="SMD",
            solvent_id="methanol",
        )
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=3,
            shared={},
            opt_settings=opt_settings,
            sp_settings=sp_settings,
        )
        assert settings.solvent_model == "SMD"
        assert settings.solvent_id == "methanol"

    def test_opt_kwargs_do_not_override_explicit_pka_kwargs(self):
        opt_settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", charge=5
        )
        shared = {"charge": 0}
        settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=1, shared=shared, opt_settings=opt_settings
        )
        assert settings.charge == 0


class TestTDDFTRouteString:
    def test_eqsolv_none_omits_eqsolv_option(self):
        from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

        settings = GaussianTDDFTJobSettings(
            functional="cam-b3lyp",
            basis="def2svp",
            jobtype="sp",
            eqsolv=None,
        )
        route = settings._get_route_string_from_jobtype()
        assert "TD(singlets,nstates=3,root=1)" in route

    def test_valid_eqsolv_option_included(self):
        from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

        settings = GaussianTDDFTJobSettings(
            functional="cam-b3lyp",
            basis="def2svp",
            jobtype="sp",
            eqsolv="noneqsolv",
        )
        route = settings._get_route_string_from_jobtype()
        assert "TD(singlets,nstates=3,root=1,noneqsolv)" in route

    def test_invalid_eqsolv_option_raises_assertion(self):
        from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings

        settings = GaussianTDDFTJobSettings(
            functional="cam-b3lyp",
            basis="def2svp",
            jobtype="sp",
            eqsolv="bogus",
        )
        with pytest.raises(AssertionError, match="equilibrium solvation"):
            settings._get_route_string_from_jobtype()
