"""
Direct unit tests for chemsmart.jobs.gaussian.settings.GaussianQMMMJobSettings
covering legacy charge/multiplicity parameter fallbacks, ONIOM string
generation edge cases, dieze_tag/irc route branches, the "cannot
override" error paths in _get_charge_and_multiplicity, and __eq__ --
none of which had direct test coverage (tests/test_GaussianSettings.py
covers the well-formed 2-/3-layer happy paths).
"""

import pytest

from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings


class TestLegacyChargeMultiplicityFallback:
    def test_canonical_charge_intermediate_used_directly(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            charge_intermediate=5,
            mult_intermediate=6,
        )
        assert settings.charge_intermediate == 5
        assert settings.mult_intermediate == 6
        assert settings.int_charge == 5
        assert settings.int_multiplicity == 6


class TestRouteStringDiezeTagAndIrc:
    def test_dieze_tag_prepended(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            parent_jobtype="sp",
            dieze_tag="p",
        )
        assert settings.route_string.startswith("#p")

    def test_irc_parent_jobtype_adds_irc_keyword(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            parent_jobtype="irc",
        )
        assert " irc" in settings.route_string


class TestOniomStringBranches:
    def test_medium_only_level_specified(self):
        settings = GaussianQMMMJobSettings(
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            charge_total=0,
            mult_total=1,
        )
        oniom_string = settings._get_oniom_string()
        assert oniom_string == " oniom(b3lyp/6-31g(d))"

    def test_low_only_level_specified(self):
        settings = GaussianQMMMJobSettings(
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        oniom_string = settings._get_oniom_string()
        assert oniom_string == " oniom(uff)"

    def test_no_levels_specified_bare_oniom(self):
        settings = GaussianQMMMJobSettings(charge_total=0, mult_total=1)
        assert settings._get_oniom_string() == " oniom"


class TestGetQmmmLevelOfTheoryString:
    def test_warns_when_jobtype_none(self, caplog):
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            jobtype=None,
            parent_jobtype="sp",
        )
        with caplog.at_level("WARNING"):
            settings.get_qmmm_level_of_theory_string()
        assert "Job type not specified" in caplog.text


class TestChargeAndMultiplicityErrors:
    def test_two_layer_inconsistent_specification_raises(self):
        # Only model_charge given (not model_multiplicity) for a 2-layer
        # job (no medium level) -- doesn't match any of the "fill from
        # a single reference level" patterns, so it must raise.
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            model_charge=-1,
        )
        with pytest.raises(ValueError, match="cannot override"):
            settings.charge_and_multiplicity_string

    def test_three_layer_fills_from_intermediate_when_only_real_and_int_set(
        self,
    ):
        # This is the only reachable "partial fill" branch other than
        # "only real specified": see BUGS_FOUND.md for why the
        # "int_med only" and "model_high only" fill branches are dead
        # code (int_med/int_low and model_high/model_med/model_low are
        # always equal, since they all derive from the same single
        # int_charge/int_multiplicity and model_charge/model_multiplicity
        # constructor parameters respectively).
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            int_charge=-1,
            int_multiplicity=2,
        )
        assert (
            settings.charge_and_multiplicity_string
            == "0 1 -1 2 -1 2 -1 2 -1 2 -1 2"
        )

    def test_three_layer_inconsistent_specification_raises(self):
        # int_charge set without int_multiplicity, and model_charge set
        # without model_multiplicity: doesn't match any single-level
        # fill pattern, so it must raise.
        settings = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
            int_charge=1,
            model_charge=0,
        )
        with pytest.raises(ValueError, match="cannot override"):
            settings.charge_and_multiplicity_string


class TestEquality:
    def test_different_types_returns_not_implemented(self):
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        assert settings.__eq__(object()) is NotImplemented

    def test_equal_settings(self):
        kwargs = dict(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        s1 = GaussianQMMMJobSettings(**kwargs)
        s2 = GaussianQMMMJobSettings(**kwargs)
        assert s1 == s2

    def test_unequal_settings_logs_diff(self, caplog):
        s1 = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        s2 = GaussianQMMMJobSettings(
            high_level_functional="wb97xd",
            high_level_basis="def2-tzvp",
            low_level_force_field="uff",
            charge_total=0,
            mult_total=1,
        )
        with caplog.at_level("INFO"):
            assert not (s1 == s2)
        assert "not equal" in caplog.text
