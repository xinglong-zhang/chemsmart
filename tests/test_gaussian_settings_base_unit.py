"""
Direct unit tests for chemsmart.jobs.gaussian.settings.GaussianJobSettings
(the base class) covering construction validation, dunder methods,
alternate constructors (from_outfile, from_filepath, PBC fallback in
from_logfile), route-string level-of-theory branches, GenECP helpers,
and solvent model validation that lacked direct test coverage.
"""

from unittest.mock import patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


def _mol(symbols=("H", "H"), positions=None):
    if positions is None:
        positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]
    return Molecule(symbols=list(symbols), positions=positions)


class TestInitValidation:
    def test_expands_user_home_in_gen_genecp_file(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            gen_genecp_file="~/genecp.txt",
        )
        assert "~" not in settings.gen_genecp_file

    def test_forces_and_freq_both_true_raises(self):
        with pytest.raises(ValueError, match="Frequency and Force"):
            GaussianJobSettings(
                functional="b3lyp", basis="sto-3g", forces=True, freq=True
            )


class TestMerge:
    def test_merge_with_none_keywords_merges_everything(self):
        base = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        other = GaussianJobSettings(functional="wb97xd", basis="def2-tzvp")
        merged = base.merge(other, keywords=None)
        assert merged.functional == "wb97xd"
        assert merged.basis == "def2-tzvp"


class TestDunderMethods:
    def test_getitem_returns_attribute(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert settings["functional"] == "b3lyp"

    def test_eq_different_types_returns_not_implemented(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert settings.__eq__(object()) is NotImplemented

    def test_eq_equal_settings(self):
        s1 = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        s2 = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert s1 == s2

    def test_eq_unequal_settings_logs_diff(self, caplog):
        s1 = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        s2 = GaussianJobSettings(functional="wb97xd", basis="def2-tzvp")
        with caplog.at_level("INFO"):
            assert not (s1 == s2)
        assert "not equal" in caplog.text


class TestAlternateConstructors:
    def test_from_outfile_merges_orca_settings(self):
        mock_orca_settings = GaussianJobSettings(
            functional="wb97xd", basis="def2-tzvp"
        )
        with patch(
            "chemsmart.io.orca.output.ORCAOutput.read_settings",
            return_value=mock_orca_settings,
        ):
            settings = GaussianJobSettings.from_outfile("fake.out")
        assert settings.functional == "wb97xd"
        assert settings.basis == "def2-tzvp"

    def test_from_logfile_falls_back_to_pbc_on_value_error(self):
        mock_pbc_settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g"
        )
        with (
            patch(
                "chemsmart.io.gaussian.output.Gaussian16Output.read_settings",
                side_effect=ValueError("not a valid non-PBC log"),
            ),
            patch(
                "chemsmart.io.gaussian.output.Gaussian16OutputWithPBC.read_settings",
                return_value=mock_pbc_settings,
            ),
        ):
            settings = GaussianJobSettings.from_logfile("fake.log")
        assert settings.functional == "b3lyp"

    def test_from_filepath_inp_extension(self):
        from unittest.mock import MagicMock

        mock_orca_settings = GaussianJobSettings(
            functional="wb97xd", basis="def2-tzvp"
        )
        mock_orca_input = MagicMock()
        mock_orca_input.read_settings.return_value = mock_orca_settings
        with patch(
            "chemsmart.io.orca.input.ORCAInput",
            return_value=mock_orca_input,
        ):
            settings = GaussianJobSettings.from_filepath("fake.inp")
        assert settings.functional == "wb97xd"

    def test_from_filepath_unsupported_extension_raises(self):
        with pytest.raises(ValueError, match="Could not create"):
            GaussianJobSettings.from_filepath("fake.xyz")


class TestRouteStringSetterAndUserInput:
    def test_route_string_setter(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        settings.route_string = "# custom route"
        assert settings._route_string == "# custom route"

    def test_user_input_route_without_dieze_tag(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            route_to_be_written="opt freq",
        )
        assert settings.route_string == "# opt freq"

    def test_user_input_route_with_dieze_tag(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            dieze_tag="p",
            route_to_be_written="opt freq",
        )
        assert settings.route_string == "#p opt freq"

    def test_user_input_route_already_prefixed(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            route_to_be_written="#p opt freq",
        )
        assert settings.route_string == "#p opt freq"


class TestGetLightElements:
    def test_returns_none_without_heavy_elements(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        assert settings.get_light_elements(_mol()) is None


class TestDiezeTagAndOptFreq:
    def test_custom_dieze_tag(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="sp", dieze_tag="p"
        )
        assert settings._get_dieze_tag().startswith("#p")

    @pytest.mark.parametrize(
        "jobtype,expected",
        [
            ("opt", "opt=(maxstep=5)"),
            ("modred", "opt=(modredundant,maxstep=5)"),
            ("scan", "opt=(modredundant,maxstep=5)"),
            ("sp", ""),
        ],
    )
    def test_additional_opt_options_by_jobtype(self, jobtype, expected):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype=jobtype,
            additional_opt_options_in_route="maxstep=5",
        )
        route = settings._get_dieze_tag()
        if expected:
            assert expected in route
        if jobtype in ("modred",):
            assert settings.freq is True
        if jobtype in ("scan", "sp"):
            assert settings.freq is False

    def test_ts_additional_opt_options_without_calcall(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ts",
            additional_opt_options_in_route="maxstep=5",
        )
        route = settings._get_dieze_tag()
        assert "opt=(ts,calcfc,noeigentest,maxstep=5)" in route

    def test_ts_additional_opt_options_with_calcall(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="ts",
            additional_opt_options_in_route="calcall",
        )
        route = settings._get_dieze_tag()
        assert "opt=(ts,noeigentest,calcall)" in route

    def test_jobtype_not_matching_any_branch_falls_through(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="nci",
            additional_opt_options_in_route=None,
        )
        route = settings._get_dieze_tag()
        assert "opt" not in route

    def test_jobtype_not_matching_any_branch_with_additional_opt_options(
        self,
    ):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="nci",
            additional_opt_options_in_route="maxstep=5",
        )
        route = settings._get_dieze_tag()
        assert "opt" not in route

    def test_numfreq_without_freq(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
            freq=False,
            numfreq=True,
        )
        route = settings._get_dieze_tag()
        assert "freq=numer" in route

    def test_freq_and_numfreq_both_true_raises(self):
        # jobtype="sp" would itself force freq=False as a side effect of
        # _get_dieze_tag(), so use a jobtype that doesn't touch freq/numfreq
        # to let the manually-set combination reach the final check.
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="nci"
        )
        settings.freq = True
        settings.numfreq = True
        with pytest.raises(ValueError, match="cannot be True at the same"):
            settings._get_dieze_tag()


class TestLevelOfTheoryString:
    def test_semiempirical_with_basis_warns(self, caplog):
        settings = GaussianJobSettings(
            semiempirical="PM6", basis="sto-3g", jobtype="sp"
        )
        with caplog.at_level("WARNING"):
            route = settings._get_level_of_theory_string()
        assert "PM6" in route
        assert "not required" in caplog.text

    def test_semiempirical_without_basis_no_warning(self, caplog):
        settings = GaussianJobSettings(semiempirical="PM6", jobtype="sp")
        with caplog.at_level("WARNING"):
            route = settings._get_level_of_theory_string()
        assert "PM6" in route
        assert "not required" not in caplog.text

    def test_ab_initio_with_basis_succeeds(self):
        settings = GaussianJobSettings(
            ab_initio="MP2", basis="sto-3g", jobtype="sp"
        )
        route = settings._get_level_of_theory_string()
        assert "MP2 sto-3g" in route

    def test_ab_initio_without_basis_raises(self):
        settings = GaussianJobSettings(ab_initio="MP2", jobtype="sp")
        with pytest.raises(ValueError, match="ab initio methods"):
            settings._get_level_of_theory_string()

    def test_functional_without_basis_raises(self):
        settings = GaussianJobSettings(functional="b3lyp", jobtype="sp")
        with pytest.raises(ValueError, match="DFT methods"):
            settings._get_level_of_theory_string()

    def test_both_ab_initio_and_functional_raises(self):
        settings = GaussianJobSettings(
            ab_initio="MP2",
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
        )
        with pytest.raises(ValueError, match="Both ab initio and DFT"):
            settings._get_level_of_theory_string()

    def test_no_method_provided_raises(self):
        settings = GaussianJobSettings(jobtype="sp")
        with pytest.raises(ValueError, match="No computational method"):
            settings._get_level_of_theory_string()

    def test_forces_keyword_added(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="sp", forces=True
        )
        route = settings._get_level_of_theory_string()
        assert " force" in route

    def test_custom_solvent_without_model_or_id_uses_pcm_default(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
            custom_solvent="SolventName=custom",
        )
        route = settings._get_level_of_theory_string()
        assert "scrf=(pcm,read" in route

    def test_incomplete_solvation_spec_raises(self):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            jobtype="sp",
            solvent_model="SMD",
        )
        with pytest.raises(ValueError, match="Both solvent model"):
            settings._get_level_of_theory_string()

    def test_nci_jobtype_adds_output_wfn(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="nci"
        )
        route = settings._get_route_string_from_jobtype()
        assert "output=wfn" in route

    def test_wbi_jobtype_adds_pop_nboread(self):
        settings = GaussianJobSettings(
            functional="b3lyp", basis="sto-3g", jobtype="wbi"
        )
        route = settings._get_route_string_from_jobtype()
        assert "pop=nboread" in route


class TestGenecpSection:
    def test_raises_when_neither_elements_nor_file_specified(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="genecp")
        with pytest.raises(ValueError, match="Could not get GenECPSection"):
            settings.get_genecp_section(_mol())


class TestDetermineBasisKeyword:
    def test_no_heavy_elements_in_structure_without_light_basis_returns_basis(
        self,
    ):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="genecp",
            heavy_elements=["Fe"],
            light_elements_basis=None,
        )
        # Molecule contains no Fe, so no heavy elements are present.
        result = settings.determine_basis_keyword(_mol(symbols=["H", "H"]))
        assert result == "genecp"


class TestCheckSolvent:
    def test_unsupported_solvent_model_raises(self):
        settings = GaussianJobSettings(functional="b3lyp", basis="sto-3g")
        with pytest.raises(ValueError, match="not in"):
            settings._check_solvent("totally_bogus_solvent_model")
