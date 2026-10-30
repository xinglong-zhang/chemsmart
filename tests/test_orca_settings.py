"""Tests for the ORCA job-settings classes in ``chemsmart.jobs.orca.settings``:
``ORCAJobSettings`` and its job-type subclasses (pKa, IRC, TS, QMMM, NEB),
covering route-string construction, solvent handling, merging, and the
per-job-type block options that feed ``ORCAInputWriter``.
"""

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAJobSettings,
    ORCANEBJobSettings,
    ORCApKaJobSettings,
    ORCAQMMMJobSettings,
    ORCATSJobSettings,
)
from chemsmart.jobs.orca.writer import ORCAInputWriter


class TestORCAJobSettingsInit:
    def test_default_construction(self):
        settings = ORCAJobSettings()
        assert settings.freq is False
        assert settings.gbw is True

    def test_forces_and_freq_conflict_raises(self):
        with pytest.raises(ValueError, match="Frequency and Force"):
            ORCAJobSettings(forces=True, freq=True)

    def test_forces_and_numfreq_conflict_raises(self):
        with pytest.raises(ValueError, match="Frequency and Force"):
            ORCAJobSettings(forces=True, numfreq=True)

    def test_forces_alone_is_fine(self):
        settings = ORCAJobSettings(forces=True)
        assert settings.forces is True


class TestMergeCopyGetitem:
    def test_merge_with_keywords_only_updates_those(self):
        base = ORCAJobSettings(charge=0, multiplicity=1, basis="sto-3g")
        merged = base.merge({"charge": 1, "multiplicity": 2, "basis": "x"})
        assert merged.charge == 1
        assert merged.multiplicity == 2
        assert merged.basis == "sto-3g"  # not in default keywords, untouched

    def test_merge_all_updates_everything_present(self):
        base = ORCAJobSettings(charge=0, multiplicity=1, basis="sto-3g")
        merged = base.merge({"basis": "def2-svp"}, merge_all=True)
        assert merged.basis == "def2-svp"

    def test_merge_with_settings_object(self):
        base = ORCAJobSettings(charge=0)
        other = ORCAJobSettings(charge=5)
        merged = base.merge(other, keywords=("charge",))
        assert merged.charge == 5

    def test_merge_returns_new_object(self):
        base = ORCAJobSettings(charge=0)
        merged = base.merge({"charge": 1})
        assert merged is not base

    def test_merge_with_keywords_none_uses_full_other_dict(self):
        base = ORCAJobSettings(charge=0, multiplicity=1, basis="sto-3g")
        merged = base.merge({"basis": "def2-svp"}, keywords=None)
        assert merged.basis == "def2-svp"

    def test_copy_is_independent_deep_copy(self):
        base = ORCAJobSettings(heavy_elements=["Fe"])
        copied = base.copy()
        copied.heavy_elements.append("Zn")
        assert base.heavy_elements == ["Fe"]

    def test_getitem(self):
        settings = ORCAJobSettings(charge=3)
        assert settings["charge"] == 3


class TestEqualityBug:
    """`ORCAJobSettings.__eq__` unconditionally pops
    "append_additional_info" from both operands' __dict__, but that
    attribute is a Gaussian-only field that's never set on ORCA
    settings objects -- so equality comparisons crash. See
    BUGS_FOUND.md (bug #11)."""

    def test_equality_crashes_with_key_error(self):
        a = ORCAJobSettings.default()
        b = ORCAJobSettings.default()
        with pytest.raises(KeyError, match="append_additional_info"):
            a == b

    def test_equality_crashes_even_for_identical_object(self):
        a = ORCAJobSettings.default()
        with pytest.raises(KeyError):
            a == a

    def test_different_types_return_notimplemented_without_crashing(self):
        a = ORCAJobSettings.default()
        # `!=` against an unrelated type falls back to identity via
        # NotImplemented on both sides, never reaching the buggy pop.
        assert (a == "not a settings object") is False

    def test_comparison_body_once_attribute_present(self):
        # Demonstrates the comparison logic itself is correct; the bug
        # is purely the unconditional, default-less `.pop()`.
        a = ORCAJobSettings.default()
        b = ORCAJobSettings.default()
        a.append_additional_info = None
        b.append_additional_info = None
        assert a == b
        b.charge = 99
        a.append_additional_info = None
        b.append_additional_info = None
        assert a != b


class TestFromFileClassmethods:
    def test_from_comfile_merges_gaussian_settings(self, tmp_path):
        com_path = tmp_path / "job.com"
        com_path.write_text("dummy")
        fake_gaussian_settings = SimpleNamespace(basis="sto-3g")
        with patch(
            "chemsmart.io.gaussian.input.Gaussian16Input"
        ) as mocked_cls:
            mocked_cls.return_value.read_settings.return_value = (
                fake_gaussian_settings
            )
            settings = ORCAJobSettings.from_comfile(str(com_path))
        assert settings.basis == "sto-3g"

    def test_from_logfile_merges_and_normalizes_def2_basis(self, tmp_path):
        log_path = tmp_path / "job.log"
        log_path.write_text("dummy")
        fake_gaussian_settings = SimpleNamespace(basis="def2svp")
        with patch(
            "chemsmart.io.gaussian.output.Gaussian16Output"
        ) as mocked_cls:
            mocked_cls.return_value.read_settings.return_value = (
                fake_gaussian_settings
            )
            settings = ORCAJobSettings.from_logfile(str(log_path))
        assert settings.basis == "def2-svp"

    def test_from_logfile_leaves_non_def2_basis_untouched(self, tmp_path):
        log_path = tmp_path / "job.log"
        log_path.write_text("dummy")
        fake_gaussian_settings = SimpleNamespace(basis="sto-3g")
        with patch(
            "chemsmart.io.gaussian.output.Gaussian16Output"
        ) as mocked_cls:
            mocked_cls.return_value.read_settings.return_value = (
                fake_gaussian_settings
            )
            settings = ORCAJobSettings.from_logfile(str(log_path))
        assert settings.basis == "sto-3g"

    def test_from_inpfile_returns_parsed_settings(self, tmp_path):
        inp_path = tmp_path / "job.inp"
        inp_path.write_text("dummy")
        expected = ORCAJobSettings.default()
        with patch("chemsmart.io.orca.input.ORCAInput") as mocked_cls:
            mocked_cls.return_value.read_settings.return_value = expected
            result = ORCAJobSettings.from_inpfile(str(inp_path))
        assert result is expected

    def test_from_outfile_returns_parsed_settings(self, tmp_path):
        out_path = tmp_path / "job.out"
        out_path.write_text("dummy")
        expected = ORCAJobSettings.default()
        with patch("chemsmart.io.orca.output.ORCAOutput") as mocked_cls:
            mocked_cls.return_value.read_settings.return_value = expected
            result = ORCAJobSettings.from_outfile(str(out_path))
        assert result is expected

    def test_from_xyzfile_returns_default(self):
        settings = ORCAJobSettings.from_xyzfile()
        assert settings.freq is True  # ORCAJobSettings.default() sets freq

    def test_from_filepath_dispatches_by_extension(self, tmp_path):
        with (
            patch.object(
                ORCAJobSettings, "from_comfile", return_value="COM"
            ) as com,
            patch.object(
                ORCAJobSettings, "from_logfile", return_value="LOG"
            ) as log,
            patch.object(
                ORCAJobSettings, "from_inpfile", return_value="INP"
            ) as inp,
            patch.object(
                ORCAJobSettings, "from_outfile", return_value="OUT"
            ) as out,
            patch.object(
                ORCAJobSettings, "from_xyzfile", return_value="XYZ"
            ) as xyz,
        ):
            assert ORCAJobSettings.from_filepath("a.com") == "COM"
            assert ORCAJobSettings.from_filepath("a.gjf") == "COM"
            assert ORCAJobSettings.from_filepath("a.log") == "LOG"
            assert ORCAJobSettings.from_filepath("a.inp") == "INP"
            assert ORCAJobSettings.from_filepath("a.out") == "OUT"
            assert ORCAJobSettings.from_filepath("a.xyz") == "XYZ"
        assert com.call_count == 2  # .com and .gjf
        log.assert_called_once()
        inp.assert_called_once()
        out.assert_called_once()
        xyz.assert_called_once()

    def test_from_filepath_unknown_extension_returns_none(self):
        assert ORCAJobSettings.from_filepath("a.pdb") is None


class TestRouteStringUserInput:
    def test_user_route_string_gets_bang_prefix(self):
        settings = ORCAJobSettings(route_to_be_written="B3LYP def2-SVP")
        assert settings.route_string == "! B3LYP def2-SVP"

    def test_user_route_string_keeps_existing_bang_prefix(self):
        settings = ORCAJobSettings(route_to_be_written="! B3LYP def2-SVP")
        assert settings.route_string == "! B3LYP def2-SVP"


class TestRouteStringFromJobtype:
    def test_opt_jobtype(self):
        settings = ORCAJobSettings(
            jobtype="opt", functional="B3LYP", basis="def2-SVP"
        )
        assert "Opt" in settings.route_string
        assert "B3LYP def2-SVP" in settings.route_string

    def test_ts_jobtype_uses_optts_keyword(self):
        settings = ORCAJobSettings(
            jobtype="ts", functional="B3LYP", basis="def2-SVP"
        )
        assert "OptTS" in settings.route_string

    def test_irc_jobtype(self):
        settings = ORCAJobSettings(
            jobtype="irc", functional="B3LYP", basis="def2-SVP"
        )
        assert "IRC" in settings.route_string

    def test_sp_jobtype_has_no_opt_keyword(self):
        settings = ORCAJobSettings(
            jobtype="sp", functional="B3LYP", basis="def2-SVP"
        )
        assert "Opt" not in settings.route_string

    def test_freq_and_numfreq_both_true_raises(self):
        settings = ORCAJobSettings(
            functional="B3LYP", basis="def2-SVP", freq=True, numfreq=True
        )
        with pytest.raises(ValueError, match="Cannot specify both"):
            _ = settings.route_string

    def test_freq_appends_freq_keyword(self):
        settings = ORCAJobSettings(
            functional="B3LYP", basis="def2-SVP", freq=True
        )
        assert " Freq" in settings.route_string

    def test_numfreq_appends_numfreq_keyword(self):
        settings = ORCAJobSettings(
            functional="B3LYP", basis="def2-SVP", numfreq=True
        )
        assert " NumFreq" in settings.route_string

    def test_defgrid_appended(self):
        settings = ORCAJobSettings(
            functional="B3LYP", basis="def2-SVP", defgrid="defgrid3"
        )
        assert "defgrid3" in settings.route_string

    def test_scf_tol_gets_scf_suffix(self):
        settings = ORCAJobSettings(
            functional="B3LYP", basis="def2-SVP", scf_tol="tight"
        )
        assert "tightSCF" in settings.route_string

    def test_scf_tol_with_existing_scf_suffix_untouched(self):
        settings = ORCAJobSettings(
            functional="B3LYP", basis="def2-SVP", scf_tol="TightSCF"
        )
        assert "TightSCF" in settings.route_string

    def test_scf_algorithm_appended(self):
        settings = ORCAJobSettings(
            functional="B3LYP", basis="def2-SVP", scf_algorithm="KDIIS"
        )
        assert "KDIIS" in settings.route_string

    def test_custom_solvent_with_solvent_id(self):
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            custom_solvent="Epsilon 16.7\n",
            solvent_id="chloroform",
        )
        assert "CPCM(chloroform)" in settings.route_string

    def test_custom_solvent_without_solvent_id(self):
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            custom_solvent="Epsilon 16.7\n",
        )
        assert settings.route_string.strip().endswith("CPCM")

    def test_solvent_model_and_id_both_set(self):
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="smd",
            solvent_id="water",
        )
        assert "SMD(water)" in settings.route_string

    def test_solvent_model_without_id(self):
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="cpcm",
        )
        assert settings.route_string.strip().endswith("CPCM")

    def test_solvent_id_without_model_warns_and_defaults_to_cpcm(self):
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            solvent_id="water",
        )
        assert "CPCM(water)" in settings.route_string

    def test_neither_solvent_model_nor_id(self):
        settings = ORCAJobSettings(functional="B3LYP", basis="def2-SVP")
        assert "CPCM" not in settings.route_string

    def test_solvent_model_deduplicated_when_repeated(self):
        # Contrived: use the solvent keyword itself as the "functional"
        # so the level-of-theory string and the solvent route keyword
        # both contain a standalone "CPCM" token, forcing the
        # deduplication branch to run.
        settings = ORCAJobSettings(
            functional="CPCM",
            basis="def2-SVP",
            solvent_model="cpcm",
            solvent_id="water",
        )
        import re as _re

        route = settings.route_string
        occurrences = len(_re.findall(r"\bcpcm\b", route, _re.IGNORECASE))
        assert occurrences == 1
        assert "CPCM(water)" in route


class TestGetLevelOfTheory:
    def test_both_ab_initio_and_functional_raises(self):
        settings = ORCAJobSettings(
            ab_initio="MP2", functional="B3LYP", basis="def2-SVP"
        )
        with pytest.raises(ValueError, match="both ab initio and DFT"):
            settings._get_level_of_theory()

    def test_neither_ab_initio_nor_functional_raises(self):
        settings = ORCAJobSettings(basis="def2-SVP")
        with pytest.raises(ValueError, match="neither ab initio nor DFT"):
            settings._get_level_of_theory()

    def test_neither_method_allowed_for_qmmm_jobtype(self):
        settings = ORCAJobSettings(jobtype="qmmm", basis="def2-SVP")
        # Should not raise despite missing ab_initio/functional.
        assert settings._get_level_of_theory() == " def2-SVP"

    def test_missing_basis_raises(self):
        settings = ORCAJobSettings(functional="B3LYP")
        with pytest.raises(ValueError, match="basis is missing"):
            settings._get_level_of_theory()

    def test_missing_basis_allowed_for_qmmm_jobtype(self):
        settings = ORCAJobSettings(jobtype="qmmm", functional="B3LYP")
        assert settings._get_level_of_theory() == "B3LYP"

    def test_aux_and_extrapolation_basis_appended(self):
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            aux_basis="def2/J",
            extrapolation_basis="extrapolate(2/3,def2)",
        )
        theory = settings._get_level_of_theory()
        assert theory == "B3LYP def2-SVP def2/J extrapolate(2/3,def2)"

    def test_ab_initio_only_used_when_functional_not_set(self):
        settings = ORCAJobSettings(ab_initio="MP2", basis="def2-SVP")
        assert settings._get_level_of_theory() == "MP2 def2-SVP"


class TestWriteGeometry:
    def test_writes_charge_multiplicity_and_coordinates(
        self, ethanol_molecule, tmp_path
    ):
        settings = ORCAJobSettings(charge=0, multiplicity=1)
        outfile = tmp_path / "geom.inp"
        with open(outfile, "w") as f:
            settings._write_geometry(f, ethanol_molecule)
        content = outfile.read_text()
        assert content.startswith("* xyz 0 1\n")
        assert content.strip().endswith("*")
        assert "O" in content

    def test_missing_charge_raises(self, ethanol_molecule, tmp_path):
        settings = ORCAJobSettings(multiplicity=1)
        outfile = tmp_path / "geom.inp"
        with open(outfile, "w") as f:
            with pytest.raises(AssertionError, match="No charge"):
                settings._write_geometry(f, ethanol_molecule)

    def test_missing_multiplicity_raises(self, ethanol_molecule, tmp_path):
        settings = ORCAJobSettings(charge=0)
        outfile = tmp_path / "geom.inp"
        with open(outfile, "w") as f:
            with pytest.raises(AssertionError, match="No multiplicity"):
                settings._write_geometry(f, ethanol_molecule)

    def test_missing_geometry_raises(self, tmp_path):
        settings = ORCAJobSettings(charge=0, multiplicity=1)
        outfile = tmp_path / "geom.inp"
        with open(outfile, "w") as f:
            with pytest.raises(AssertionError, match="No molecular geometry"):
                settings._write_geometry(f, None)


class TestSolventHelpers:
    @pytest.mark.parametrize(
        "model,expected",
        [
            ("cpcm", "CPCM"),
            ("cpcmc", "CPCMC"),
            ("smd", "SMD"),
            ("SMD", "SMD"),
            ("cosmors", "COSMORS"),
            (None, "CPCM"),
            ("unknown", "CPCM"),
        ],
    )
    def test_get_solvent_route_keyword(self, model, expected):
        settings = ORCAJobSettings(solvent_model=model)
        assert settings._get_solvent_route_keyword() == expected

    def test_check_solvent_accepts_valid_model(self):
        settings = ORCAJobSettings()
        settings._check_solvent("cpcm")  # must not raise

    def test_check_solvent_rejects_invalid_model(self):
        settings = ORCAJobSettings()
        with pytest.raises(ValueError, match="not in"):
            settings._check_solvent("not_a_model")


class TestDefault:
    def test_default_settings_values(self):
        settings = ORCAJobSettings.default()
        assert settings.freq is True
        assert settings.gbw is True
        assert settings.forces is False


ETHANOL_XYZ = """9

O -1.1712 0.2997 0.0000
C -0.0463 -0.5665 0.0000
C 1.2175 0.2668 0.0000
H -0.0958 -1.2120 0.8819
H -0.0952 -1.1938 -0.8946
H 2.1050 -0.3720 -0.0177
H 1.2426 0.9307 -0.8704
H 1.2616 0.9052 0.8886
H -1.1291 0.8364 0.8099
"""


@pytest.fixture()
def reference_xyz_file(tmp_path):
    path = tmp_path / "href.xyz"
    path.write_text(ETHANOL_XYZ)
    return str(path)


@pytest.fixture()
def pka_settings():
    return ORCApKaJobSettings(proton_index=4, charge=0, multiplicity=1)


@pytest.fixture()
def pka_settings_with_reference(reference_xyz_file):
    return ORCApKaJobSettings(
        proton_index=4,
        charge=0,
        multiplicity=1,
        reference_file=reference_xyz_file,
        reference_proton_index=4,
        reference_charge=0,
        reference_multiplicity=1,
    )


class TestORCApKaJobSettingsInit:
    def test_default_delta_g_proton(self, pka_settings):
        assert (
            pka_settings.delta_G_proton
            == ORCApKaJobSettings.DEFAULT_DELTA_G_PROTON
        )

    def test_custom_delta_g_proton(self):
        settings = ORCApKaJobSettings(proton_index=1, delta_G_proton=-100.0)
        assert settings.delta_G_proton == -100.0

    def test_reference_fields_none_without_reference_file(self, pka_settings):
        assert pka_settings.reference_file is None
        assert pka_settings.reference_proton_index is None

    def test_reference_fields_set_with_reference_file(
        self, pka_settings_with_reference
    ):
        assert pka_settings_with_reference.reference_file is not None
        assert pka_settings_with_reference.reference_proton_index == 4

    def test_deprecated_thermodynamic_cycle_kwarg_sets_scheme(self, caplog):
        import logging

        with caplog.at_level(logging.WARNING):
            settings = ORCApKaJobSettings(
                proton_index=1, thermodynamic_cycle="direct"
            )
        assert settings.scheme == "direct"
        assert "deprecated" in caplog.text

    def test_default_title_used_when_not_provided(self, pka_settings):
        assert pka_settings.title == "ORCA pKa calculation job"

    def test_custom_title_not_overridden(self):
        settings = ORCApKaJobSettings(proton_index=1, title="My custom title")
        assert settings.title == "My custom title"


class TestBuildOrcaPkaSettings:
    def test_uses_shared_solvent_when_provided(self):
        shared = {
            "solvent_model": "smd",
            "solvent_id": "water",
            "scheme": "direct",
            "reference": None,
            "reference_proton_index": None,
            "reference_charge": None,
            "reference_multiplicity": None,
            "reference_conjugate_base_charge": None,
            "reference_conjugate_base_multiplicity": None,
            "delta_g_proton": None,
            "conjugate_base_charge": None,
            "conjugate_base_multiplicity": None,
            "temperature": 298.15,
            "concentration": 1.0,
            "pressure": 1.0,
            "cutoff_entropy_grimme": 100.0,
            "cutoff_enthalpy": 100.0,
        }
        opt_settings = ORCAJobSettings(
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
        )
        settings = ORCApKaJobSettings.build_orca_pka_settings(
            proton_index=4, shared=shared, opt_settings=opt_settings
        )
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.charge == 0
        assert settings.functional == "B3LYP"

    def test_defaults_solvent_when_not_provided_anywhere(self):
        shared = {
            "solvent_model": None,
            "solvent_id": None,
            "scheme": "direct",
            "reference": None,
            "reference_proton_index": None,
            "reference_charge": None,
            "reference_multiplicity": None,
            "reference_conjugate_base_charge": None,
            "reference_conjugate_base_multiplicity": None,
            "delta_g_proton": None,
            "conjugate_base_charge": None,
            "conjugate_base_multiplicity": None,
            "temperature": 298.15,
            "concentration": 1.0,
            "pressure": 1.0,
            "cutoff_entropy_grimme": 100.0,
            "cutoff_enthalpy": 100.0,
        }
        opt_settings = ORCAJobSettings(
            charge=0, multiplicity=1, functional="B3LYP", basis="def2-SVP"
        )
        settings = ORCApKaJobSettings.build_orca_pka_settings(
            proton_index=4, shared=shared, opt_settings=opt_settings
        )
        assert settings.solvent_model == "CPCM"
        assert settings.solvent_id == "water"

    def test_falls_back_to_opt_settings_solvent(self):
        shared = {
            "solvent_model": None,
            "solvent_id": None,
            "scheme": "direct",
            "reference": None,
            "reference_proton_index": None,
            "reference_charge": None,
            "reference_multiplicity": None,
            "reference_conjugate_base_charge": None,
            "reference_conjugate_base_multiplicity": None,
            "delta_g_proton": None,
            "conjugate_base_charge": None,
            "conjugate_base_multiplicity": None,
            "temperature": 298.15,
            "concentration": 1.0,
            "pressure": 1.0,
            "cutoff_entropy_grimme": 100.0,
            "cutoff_enthalpy": 100.0,
        }
        opt_settings = ORCAJobSettings(
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="cpcmc",
            solvent_id="chloroform",
        )
        settings = ORCApKaJobSettings.build_orca_pka_settings(
            proton_index=4, shared=shared, opt_settings=opt_settings
        )
        assert settings.solvent_model == "cpcmc"
        assert settings.solvent_id == "chloroform"

    def test_handles_opt_settings_without_solvent_attrs(self):
        shared = {
            "solvent_model": None,
            "solvent_id": None,
            "scheme": "direct",
            "reference": None,
            "reference_proton_index": None,
            "reference_charge": None,
            "reference_multiplicity": None,
            "reference_conjugate_base_charge": None,
            "reference_conjugate_base_multiplicity": None,
            "delta_g_proton": None,
            "conjugate_base_charge": None,
            "conjugate_base_multiplicity": None,
            "temperature": 298.15,
            "concentration": 1.0,
            "pressure": 1.0,
            "cutoff_entropy_grimme": 100.0,
            "cutoff_enthalpy": 100.0,
        }
        opt_settings = SimpleNamespace(
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
            ab_initio=None,
            dispersion=None,
            aux_basis=None,
            defgrid=None,
            semiempirical=None,
            additional_route_parameters=None,
            gen_genecp_file=None,
            heavy_elements=None,
            heavy_elements_basis=None,
            light_elements_basis=None,
        )
        settings = ORCApKaJobSettings.build_orca_pka_settings(
            proton_index=4, shared=shared, opt_settings=opt_settings
        )
        assert settings.solvent_model == "CPCM"
        assert settings.solvent_id == "water"


class TestHasReferenceFileAndValidation:
    def test_has_reference_file_false_by_default(self, pka_settings):
        assert pka_settings.has_reference_file is False

    def test_has_reference_file_true_when_set(
        self, pka_settings_with_reference
    ):
        assert pka_settings_with_reference.has_reference_file is True

    def test_has_reference_file_false_for_direct_scheme(
        self, reference_xyz_file
    ):
        settings = ORCApKaJobSettings(
            proton_index=4, scheme="direct", reference_file=reference_xyz_file
        )
        assert settings.has_reference_file is False

    def test_validate_reference_settings_raises_without_file(
        self, pka_settings
    ):
        with pytest.raises(ValueError, match="Reference acid file"):
            pka_settings.validate_reference_settings()

    def test_validate_reference_settings_raises_when_missing_fields(
        self, reference_xyz_file
    ):
        settings = ORCApKaJobSettings(
            proton_index=4, reference_file=reference_xyz_file
        )
        with pytest.raises(ValueError, match="Missing required"):
            settings.validate_reference_settings()

    def test_validate_reference_settings_passes_when_complete(
        self, pka_settings_with_reference
    ):
        pka_settings_with_reference.validate_reference_settings()  # no raise


class TestReferenceMoleculeConstruction:
    def test_get_reference_molecule_sets_charge_and_multiplicity(
        self, pka_settings_with_reference
    ):
        ref_mol = pka_settings_with_reference.get_reference_molecule()
        assert ref_mol.charge == 0
        assert ref_mol.multiplicity == 1

    def test_get_reference_conjugate_base_molecule(
        self, pka_settings_with_reference
    ):
        ref_cb_mol = (
            pka_settings_with_reference.get_reference_conjugate_base_molecule()
        )
        assert len(ref_cb_mol) == 8  # one H removed

    def test_create_reference_conjugate_base_proton_index_none_raises(
        self, pka_settings_with_reference
    ):
        ref_mol = pka_settings_with_reference.get_reference_molecule()
        pka_settings_with_reference.reference_proton_index = None
        with pytest.raises(ValueError, match="reference_proton_index"):
            pka_settings_with_reference._create_reference_conjugate_base_molecule(
                ref_mol
            )

    def test_create_reference_conjugate_base_out_of_range_raises(
        self, pka_settings_with_reference
    ):
        pka_settings_with_reference.reference_proton_index = 99
        ref_mol = pka_settings_with_reference.get_reference_molecule()
        with pytest.raises(ValueError, match="out of range"):
            pka_settings_with_reference._create_reference_conjugate_base_molecule(
                ref_mol
            )

    def test_create_reference_conjugate_base_non_hydrogen_raises(
        self, pka_settings_with_reference
    ):
        pka_settings_with_reference.reference_proton_index = 1  # oxygen
        ref_mol = pka_settings_with_reference.get_reference_molecule()
        with pytest.raises(ValueError, match="not hydrogen"):
            pka_settings_with_reference._create_reference_conjugate_base_molecule(
                ref_mol
            )

    def test_create_reference_conjugate_base_charge_mult_defaults(
        self, pka_settings_with_reference
    ):
        ref_mol = pka_settings_with_reference.get_reference_molecule()
        ref_cb_mol = pka_settings_with_reference._create_reference_conjugate_base_molecule(
            ref_mol
        )
        assert ref_cb_mol.charge == -1  # reference_charge(0) - 1
        assert ref_cb_mol.multiplicity == 1  # reference_multiplicity

    def test_create_reference_conjugate_base_charge_mult_overrides(
        self, reference_xyz_file
    ):
        settings = ORCApKaJobSettings(
            proton_index=4,
            reference_file=reference_xyz_file,
            reference_proton_index=4,
            reference_charge=0,
            reference_multiplicity=1,
            reference_conjugate_base_charge=-2,
            reference_conjugate_base_multiplicity=3,
        )
        ref_mol = settings.get_reference_molecule()
        ref_cb_mol = settings._create_reference_conjugate_base_molecule(
            ref_mol
        )
        assert ref_cb_mol.charge == -2
        assert ref_cb_mol.multiplicity == 3


class TestProtonatedChargeMultiplicityAliases:
    def test_protonated_charge_getter(self, pka_settings):
        assert pka_settings.protonated_charge == 0

    def test_protonated_charge_setter(self, pka_settings):
        pka_settings.protonated_charge = 5
        assert pka_settings.charge == 5

    def test_protonated_multiplicity_getter(self, pka_settings):
        assert pka_settings.protonated_multiplicity == 1

    def test_protonated_multiplicity_setter(self, pka_settings):
        pka_settings.protonated_multiplicity = 3
        assert pka_settings.multiplicity == 3


class TestConjugateBaseMoleculeCreation:
    def test_create_conjugate_base_proton_index_none_raises(
        self, ethanol_molecule
    ):
        settings = ORCApKaJobSettings()
        with pytest.raises(ValueError, match="proton_index"):
            settings._create_conjugate_base_molecule(ethanol_molecule)

    def test_create_conjugate_base_out_of_range_raises(self, ethanol_molecule):
        settings = ORCApKaJobSettings(proton_index=99)
        with pytest.raises(ValueError, match="out of range"):
            settings._create_conjugate_base_molecule(ethanol_molecule)

    def test_create_conjugate_base_non_hydrogen_raises(self, ethanol_molecule):
        settings = ORCApKaJobSettings(proton_index=1)  # oxygen
        with pytest.raises(ValueError, match="not hydrogen"):
            settings._create_conjugate_base_molecule(ethanol_molecule)

    def test_create_conjugate_base_charge_mult_default_from_molecule(
        self, ethanol_molecule
    ):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(proton_index=4)
        cb_mol = settings._create_conjugate_base_molecule(ethanol_molecule)
        assert cb_mol.charge == -1
        assert cb_mol.multiplicity == 1
        assert len(cb_mol) == 8

    def test_create_conjugate_base_charge_mult_default_when_molecule_unset(
        self, ethanol_molecule
    ):
        ethanol_molecule.charge = None
        ethanol_molecule.multiplicity = None
        settings = ORCApKaJobSettings(proton_index=4)
        cb_mol = settings._create_conjugate_base_molecule(ethanol_molecule)
        assert cb_mol.charge == -1  # falls back to 0 - 1
        assert cb_mol.multiplicity == 1

    def test_create_conjugate_base_charge_mult_settings_override(
        self, ethanol_molecule
    ):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(
            proton_index=4,
            conjugate_base_charge=-3,
            conjugate_base_multiplicity=2,
        )
        cb_mol = settings._create_conjugate_base_molecule(ethanol_molecule)
        assert cb_mol.charge == -3
        assert cb_mol.multiplicity == 2

    def test_conjugate_base_molecule_delegates(self, ethanol_molecule):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(proton_index=4)
        cb_mol = settings.conjugate_base_molecule(ethanol_molecule)
        assert len(cb_mol) == 8

    def test_conjugate_pair_molecules(self, ethanol_molecule):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(proton_index=4)
        prot_mol, conj_mol = settings.conjugate_pair_molecules(
            ethanol_molecule
        )
        assert prot_mol is ethanol_molecule
        assert len(conj_mol) == 8


class TestOrcaSettingsKwargs:
    def test_returns_shared_orca_fields(self):
        settings = ORCApKaJobSettings(
            proton_index=4,
            basis="def2-SVP",
            functional="B3LYP",
            scf_tol="tight",
        )
        kwargs = settings._orca_settings_kwargs()
        assert kwargs["basis"] == "def2-SVP"
        assert kwargs["functional"] == "B3LYP"
        assert kwargs["scf_tol"] == "tight"


class TestSubJobSettingsFactories:
    def test_create_gas_phase_job_settings(self, ethanol_molecule):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(
            proton_index=4,
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="smd",
            solvent_id="water",
        )
        prot_settings, conj_settings = settings._create_gas_phase_job_settings(
            ethanol_molecule
        )
        assert prot_settings.jobtype == "opt"
        assert prot_settings.freq is True
        assert prot_settings.solvent_model is None  # gas phase: no solvent
        assert prot_settings.charge == 0
        assert conj_settings.charge == -1  # conjugate base default

    def test_create_gas_phase_job_settings_infers_charge_from_molecule(
        self, ethanol_molecule
    ):
        ethanol_molecule.charge = 2
        ethanol_molecule.multiplicity = 3
        settings = ORCApKaJobSettings(
            proton_index=4, functional="B3LYP", basis="def2-SVP"
        )
        prot_settings, _ = settings._create_gas_phase_job_settings(
            ethanol_molecule
        )
        assert prot_settings.charge == 2
        assert prot_settings.multiplicity == 3

    def test_create_gas_phase_job_settings_conjugate_base_overrides(
        self, ethanol_molecule
    ):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(
            proton_index=4,
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
            conjugate_base_charge=-5,
            conjugate_base_multiplicity=4,
        )
        _, conj_settings = settings._create_gas_phase_job_settings(
            ethanol_molecule
        )
        assert conj_settings.charge == -5
        assert conj_settings.multiplicity == 4

    def test_create_solution_phase_sp_settings(self, ethanol_molecule):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(
            proton_index=4,
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="smd",
            solvent_id="water",
        )
        prot_sp, conj_sp = settings._create_solution_phase_sp_settings(
            ethanol_molecule
        )
        assert prot_sp.jobtype == "sp"
        assert prot_sp.freq is False
        assert prot_sp.solvent_model == "smd"  # solution phase: keeps solvent
        assert conj_sp.charge == -1

    def test_create_reference_gas_phase_job_settings(
        self, pka_settings_with_reference
    ):
        ref_settings, ref_cb_settings = (
            pka_settings_with_reference._create_reference_gas_phase_job_settings()
        )
        assert ref_settings.charge == 0
        assert ref_settings.jobtype == "opt"
        assert ref_settings.freq is True
        assert ref_cb_settings.charge == -1  # reference_charge - 1

    def test_create_reference_gas_phase_job_settings_requires_reference(
        self, pka_settings
    ):
        with pytest.raises(ValueError, match="Reference acid file"):
            pka_settings._create_reference_gas_phase_job_settings()

    def test_create_reference_solution_phase_sp_settings(
        self, pka_settings_with_reference
    ):
        pka_settings_with_reference.solvent_model = "smd"
        pka_settings_with_reference.solvent_id = "water"
        ref_sp, ref_cb_sp = (
            pka_settings_with_reference._create_reference_solution_phase_sp_settings()
        )
        assert ref_sp.jobtype == "sp"
        assert ref_sp.solvent_model == "smd"
        assert ref_cb_sp.charge == -1

    def test_create_reference_solution_phase_sp_settings_conjugate_base_overrides(
        self, reference_xyz_file
    ):
        settings = ORCApKaJobSettings(
            proton_index=4,
            reference_file=reference_xyz_file,
            reference_proton_index=4,
            reference_charge=0,
            reference_multiplicity=1,
            reference_conjugate_base_charge=-9,
            reference_conjugate_base_multiplicity=5,
        )
        _, ref_cb_sp = settings._create_reference_solution_phase_sp_settings()
        assert ref_cb_sp.charge == -9
        assert ref_cb_sp.multiplicity == 5


class TestPairDelegationMethods:
    def test_conjugate_pair_job_settings_delegates(self, ethanol_molecule):
        ethanol_molecule.charge = 0
        ethanol_molecule.multiplicity = 1
        settings = ORCApKaJobSettings(
            proton_index=4, functional="B3LYP", basis="def2-SVP"
        )
        with patch.object(
            settings, "_create_gas_phase_job_settings", return_value="X"
        ) as mocked:
            result = settings.conjugate_pair_job_settings(ethanol_molecule)
        mocked.assert_called_once_with(ethanol_molecule)
        assert result == "X"

    def test_conjugate_pair_sp_job_settings_delegates(self, ethanol_molecule):
        settings = ORCApKaJobSettings(proton_index=4)
        with patch.object(
            settings,
            "_create_solution_phase_sp_settings",
            return_value="Y",
        ) as mocked:
            result = settings.conjugate_pair_sp_job_settings(ethanol_molecule)
        mocked.assert_called_once_with(ethanol_molecule)
        assert result == "Y"

    def test_reference_pair_molecules_delegates(
        self, pka_settings_with_reference
    ):
        ref_mol, ref_cb_mol = (
            pka_settings_with_reference.reference_pair_molecules()
        )
        assert ref_mol.charge == 0
        assert len(ref_cb_mol) == 8

    def test_reference_pair_job_settings_delegates(
        self, pka_settings_with_reference
    ):
        with patch.object(
            pka_settings_with_reference,
            "_create_reference_gas_phase_job_settings",
            return_value="REF_OPT",
        ) as mocked:
            result = pka_settings_with_reference.reference_pair_job_settings()
        mocked.assert_called_once()
        assert result == "REF_OPT"

    def test_reference_pair_sp_job_settings_delegates(
        self, pka_settings_with_reference
    ):
        with patch.object(
            pka_settings_with_reference,
            "_create_reference_solution_phase_sp_settings",
            return_value="REF_SP",
        ) as mocked:
            result = (
                pka_settings_with_reference.reference_pair_sp_job_settings()
            )
        mocked.assert_called_once()
        assert result == "REF_SP"


class TestORCATSJobSettings:
    def test_defaults(self):
        settings = ORCATSJobSettings()
        assert settings.tssearch_type == "optts"
        assert settings.recalc_hess == 5
        assert settings.full_scan is False

    def test_route_string_uses_optts_by_default(self):
        settings = ORCATSJobSettings(functional="B3LYP", basis="def2-SVP")
        assert "OptTS" in settings.route_string

    def test_route_string_uses_scants_when_configured(self):
        settings = ORCATSJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            tssearch_type="scants",
        )
        assert "ScanTS" in settings.route_string
        assert "OptTS" not in settings.route_string

    def test_route_string_sets_jobtype_to_ts(self):
        settings = ORCATSJobSettings(functional="B3LYP", basis="def2-SVP")
        _ = settings.route_string
        assert settings.jobtype == "ts"


class TestORCAIRCJobSettings:
    def test_defaults(self):
        settings = ORCAIRCJobSettings()
        assert settings.maxiter is None
        assert settings.direction is None

    def test_route_string_sets_jobtype_and_strips_freq(self):
        settings = ORCAIRCJobSettings(
            functional="B3LYP", basis="def2-SVP", freq=True
        )
        route = settings.route_string
        assert settings.jobtype == "irc"
        assert "freq" not in route.lower()
        assert "IRC" in route

    def test_route_string_without_freq_skips_substitution(self):
        """When neither freq nor numfreq is set, the route string never
        contains "freq" in the first place, so the re.sub substitution
        branch is simply skipped."""
        settings = ORCAIRCJobSettings(functional="B3LYP", basis="def2-SVP")
        route = settings.route_string
        assert "freq" not in route.lower()
        assert "IRC" in route

    def test_write_irc_block_noop_when_all_options_none(self, tmp_path):
        settings = ORCAIRCJobSettings()
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            settings._write_irc_block(f)
        assert outfile.read_text() == ""

    def test_write_irc_block_writes_simple_options(self, tmp_path):
        settings = ORCAIRCJobSettings(maxiter=20, printlevel=1)
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            settings._write_irc_block(f)
        content = outfile.read_text()
        assert content.startswith("%irc\n")
        assert "maxiter 20" in content
        assert "printlevel 1" in content
        assert content.strip().endswith("end")

    def test_write_irc_block_inithess_read_requires_hess_file(self, tmp_path):
        settings = ORCAIRCJobSettings(inithess="read")
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            with pytest.raises(AssertionError, match="No Hessian file"):
                settings._write_irc_block(f)

    def test_write_irc_block_inithess_read_requires_existing_file(
        self, tmp_path
    ):
        settings = ORCAIRCJobSettings(
            inithess="read", hess_filename=str(tmp_path / "missing.hess")
        )
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            with pytest.raises(AssertionError, match="is not found"):
                settings._write_irc_block(f)

    def test_write_irc_block_inithess_read_writes_hess_filename(
        self, tmp_path
    ):
        hess_file = tmp_path / "job.hess"
        hess_file.write_text("dummy")
        settings = ORCAIRCJobSettings(
            inithess="read", hess_filename=str(hess_file)
        )
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            settings._write_irc_block(f)
        content = outfile.read_text()
        assert "inithess read" in content
        assert f'Hess_Filename "{hess_file}"' in content
        # hess_filename itself is not written as a separate key/value line
        assert f"hess_filename {hess_file}" not in content

    def test_write_irc_block_monitor_internals_true_writes_modred(
        self, tmp_path
    ):
        settings = ORCAIRCJobSettings(
            monitor_internals=True, internal_modred=[[1, 2]]
        )
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            settings._write_irc_block(f)
        content = outfile.read_text()
        assert "monitor_internals\n" in content
        assert "end\n" in content

    def test_write_irc_block_monitor_internals_true_requires_modred(
        self, tmp_path
    ):
        settings = ORCAIRCJobSettings(monitor_internals=True)
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            with pytest.raises(AssertionError, match="internal modred"):
                settings._write_irc_block(f)

    def test_write_irc_block_monitor_internals_false_not_written(
        self, tmp_path
    ):
        settings = ORCAIRCJobSettings(monitor_internals=False, maxiter=10)
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            settings._write_irc_block(f)
        content = outfile.read_text()
        assert "monitor_internals" not in content
        assert "maxiter 10" in content

    def test_write_irc_block_inithess_not_read_skips_hess_filename(
        self, tmp_path
    ):
        """When inithess is set to something other than 'read' (e.g. one of
        the computed-Hessian options), the Hess_Filename line and its
        existence assertions are skipped entirely; the loop just moves on
        to the next option."""
        settings = ORCAIRCJobSettings(inithess="calc_anfreq", maxiter=5)
        outfile = tmp_path / "irc.inp"
        with open(outfile, "w") as f:
            settings._write_irc_block(f)
        content = outfile.read_text()
        assert "inithess calc_anfreq" in content
        assert "Hess_Filename" not in content
        assert "maxiter 5" in content


class TestQMMMChargeMultiplicityFallback:
    def test_high_level_takes_priority(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_high=0,
            mult_high=1,
            charge_intermediate=1,
            mult_intermediate=2,
            charge_total=2,
            mult_total=3,
        )
        assert s.charge == 0
        assert s.multiplicity == 1

    def test_intermediate_used_when_no_high(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=1,
            mult_intermediate=2,
            charge_total=2,
            mult_total=3,
        )
        assert s.charge == 1
        assert s.multiplicity == 2

    def test_total_used_when_no_high_or_intermediate(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=2,
            mult_total=3,
        )
        assert s.charge == 2
        assert s.multiplicity == 3

    def test_none_when_nothing_specified(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        assert s.charge is None
        assert s.multiplicity is None

    def test_functional_and_basis_aliased_from_high_level(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        assert s.functional == "B3LYP"
        assert s.basis == "def2-SVP"


class TestQMMMIntermediateValidation:
    def test_qm_qm2_without_intermediate_params_raises(self):
        with pytest.raises(ValueError, match="requires QM2"):
            ORCAQMMMJobSettings(
                jobtype="QM/QM2",
                high_level_functional="B3LYP",
                high_level_basis="def2-SVP",
                charge_high=0,
                mult_high=1,
            )

    def test_qm_qm2_mm_without_intermediate_params_raises(self):
        with pytest.raises(ValueError, match="requires QM2"):
            ORCAQMMMJobSettings(
                jobtype="QM/QM2/MM",
                low_level_method="ff.prms",
                charge_total=0,
                mult_total=1,
            )

    def test_intermediate_method_without_atoms_but_with_low_level_raises(
        self,
    ):
        with pytest.raises(ValueError, match="intermediate_level_atoms"):
            ORCAQMMMJobSettings(
                jobtype="QM/QM2/MM",
                intermediate_level_method="XTB2",
                low_level_method="ff.prms",
                charge_total=0,
                mult_total=1,
            )

    def test_intermediate_method_without_atoms_or_low_level_is_allowed(self):
        # The atoms-required guard is only triggered when a low-level
        # method is also supplied; a pure QM/QM2 job with just an
        # intermediate method and no MM layer sails through.
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            charge_high=0,
            mult_high=1,
        )
        assert s.intermediate_level_atoms is None

    def test_intermediate_params_for_non_qm2_jobtype_with_low_level_raises(
        self,
    ):
        with pytest.raises(ValueError, match="does not require QM2"):
            ORCAQMMMJobSettings(
                jobtype="QMMM",
                intermediate_level_functional="B3LYP",
                intermediate_level_basis="def2-SVP",
                intermediate_level_atoms=[1, 2],
                low_level_method="ff.prms",
                charge_total=0,
                mult_total=1,
            )

    def test_all_provided_intermediate_params_are_listed_when_no_low_level(
        self,
    ):
        """No low_level_method means no error is raised, but the
        param-collection loop that builds the (unused-in-this-case)
        error message must still run over every intermediate-level field
        so that each optional attribute is exercised."""
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            intermediate_level_method="XTB2",
            charge_intermediate=0,
            mult_intermediate=1,
            intermediate_level_solvation="cpcm(water)",
        )
        assert s.intermediate_level_method == "XTB2"
        assert s.charge_intermediate == 0
        assert s.mult_intermediate == 1
        assert s.intermediate_level_solvation == "cpcm(water)"

    def test_intermediate_params_for_non_qm2_jobtype_without_low_level_is_silently_accepted(
        self,
    ):
        """Documents a validation gap (see BUGS_FOUND.md #11): the
        "does not require QM2 but intermediate params were provided"
        check is only reached when ``low_level_method`` is also set, so a
        plain ``QMMM`` job can carry a fully-specified (and meaningless)
        intermediate QM2 layer without any error being raised."""
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            intermediate_level_functional="B3LYP",
            intermediate_level_basis="def2-SVP",
            charge_total=0,
            mult_total=1,
        )
        assert s.intermediate_level_functional == "B3LYP"

    def test_re_init_and_validate_recomputes_charge_and_reruns_validation(
        self,
    ):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_high=0,
            mult_high=1,
        )
        s.charge_high = 1
        s.mult_high = 3
        s.high_level_functional = "PBE0"
        s.re_init_and_validate()
        assert s.charge == 1
        assert s.multiplicity == 3
        assert s.functional == "PBE0"

        s.intermediate_level_method = None
        s.intermediate_level_functional = None
        s.intermediate_level_basis = None
        with pytest.raises(ValueError, match="requires QM2"):
            s.re_init_and_validate()

    def test_re_init_and_validate_falls_back_to_intermediate_charge(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=1,
            mult_intermediate=2,
        )
        assert s.charge == 1
        assert s.multiplicity == 2
        # Clear high-level charge/mult (already unset) and re-run: should
        # still fall back to the intermediate charge/mult pair.
        s.high_level_functional = "PBE0"
        s.re_init_and_validate()
        assert s.charge == 1
        assert s.multiplicity == 2
        assert s.functional == "PBE0"

    def test_re_init_and_validate_falls_back_to_total_charge(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            charge_total=2,
            mult_total=3,
        )
        assert s.charge == 2
        assert s.multiplicity == 3
        s.high_level_functional = "PBE0"
        s.re_init_and_validate()
        assert s.charge == 2
        assert s.multiplicity == 3
        assert s.functional == "PBE0"


class TestQMMMLevelOfTheoryString:
    def test_additive_qmmm_route_string(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            parent_jobtype="opt",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert route.startswith("! Opt QMMM")
        assert "B3LYP def2-SVP" in route

    def test_qm_qm2_route_string(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            parent_jobtype="sp",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_high=0,
            mult_high=1,
        )
        route = s._get_level_of_theory_string()
        assert "QM/XTB2" in route
        assert "MM" not in route.replace("QM/XTB2", "")

    def test_qm_qm2_mm_route_string_includes_low_level(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2/MM",
            parent_jobtype="ts",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert route.startswith("! OptTS")
        assert "QM/XTB2/MM" in route

    def test_irc_parent_jobtype(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            parent_jobtype="irc",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            freq=True,
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert route.startswith("! IRC Freq")

    def test_crystal_jobtype_route_string_is_just_the_jobtype(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
        )
        assert s._get_level_of_theory_string() == "! MOL-CRYSTAL-QMMM"

    def test_qmmm_route_string_property_delegates(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
        )
        assert s.qmmm_route_string == s._get_level_of_theory_string()

    def test_unrecognized_parent_jobtype_appends_no_opt_keyword(self):
        """A parent_jobtype that doesn't match any of opt/modred/scan/ts/
        irc/sp simply falls through without adding any job-type keyword."""
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            parent_jobtype="something-else",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert route.startswith("! QMMM")

    def test_solvent_model_appended_to_qmmm_route(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            parent_jobtype="opt",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            solvent_model="cpcm",
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert "cpcm" in route

    def test_intermediate_level_solvation_alpb_requires_xtb(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            parent_jobtype="sp",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB",
            intermediate_level_atoms=[1, 2],
            intermediate_level_solvation="alpb(water)",
            charge_high=0,
            mult_high=1,
        )
        route = s._get_level_of_theory_string()
        assert "alpb(water)" in route

    def test_intermediate_level_solvation_alpb_raises_when_not_xtb(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            parent_jobtype="sp",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_functional="B3LYP",
            intermediate_level_basis="def2-SVP",
            intermediate_level_atoms=[1, 2],
            intermediate_level_solvation="alpb(water)",
            charge_high=0,
            mult_high=1,
        )
        with pytest.raises(AssertionError, match="only compatible"):
            s._get_level_of_theory_string()

    def test_intermediate_level_solvation_unrecognized_value_appends_nothing(
        self,
    ):
        """A truthy intermediate_level_solvation that matches neither the
        ALPB/DDCOSMO/CPCMX list nor 'cpcm(water)' falls through both
        branches without appending anything."""
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            parent_jobtype="sp",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            intermediate_level_solvation="some_other_model",
            charge_high=0,
            mult_high=1,
        )
        route = s._get_level_of_theory_string()
        assert "some_other_model" not in route

    def test_intermediate_level_solvation_cpcm_water(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            parent_jobtype="sp",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            intermediate_level_solvation="cpcm(water)",
            charge_high=0,
            mult_high=1,
        )
        route = s._get_level_of_theory_string()
        assert "cpcm(water)" in route


class TestValidateAndAssignLevel:
    def test_functional_and_basis_and_builtin_together_raises(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        with pytest.raises(ValueError, match="not both"):
            s.validate_and_assign_level(
                "B3LYP", "def2-SVP", "XTB2", level_name="high_level"
            )

    def test_builtin_method_for_intermediate_level(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            None, None, "XTB2", level_name="intermediate_level"
        )
        assert result == "XTB2"

    def test_unknown_builtin_method_for_intermediate_level_returns_empty(
        self,
    ):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            None, None, "NOT-A-METHOD", level_name="intermediate_level"
        )
        assert result == ""

    def test_functional_and_basis_for_intermediate_level_is_qm2(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            "B3LYP", "def2-SVP", None, level_name="intermediate_level"
        )
        assert result == "QM2"

    def test_functional_and_basis_for_high_level(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            "B3LYP", "def2-SVP", None, level_name="high_level"
        )
        assert result == "B3LYP def2-SVP"

    def test_low_level_always_mm(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            None, None, None, level_name="low_level"
        )
        assert result == "MM"


class TestCheckCrystalQmmm:
    def test_mol_crystal_qmmm_requires_no_multiplicity(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
            mult_total=1,
        )
        with pytest.raises(AssertionError, match="Multiplicity"):
            s.check_crystal_qmmm()

    def test_mol_crystal_qmmm_sets_multiplicity_to_zero(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
        )
        s.check_crystal_qmmm()
        assert s.multiplicity == 0

    def test_mol_crystal_qmmm_requires_n_unit_cell_atoms(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM", low_level_method="ff.prms"
        )
        with pytest.raises(AssertionError, match="number of atoms"):
            s.check_crystal_qmmm()

    def test_ionic_crystal_qmmm_requires_ecp_layer_ecp(self):
        s = ORCAQMMMJobSettings(
            jobtype="IONIC-CRYSTAL-QMMM", low_level_method="ff.prms"
        )
        with pytest.raises(AssertionError, match="cECPs"):
            s.check_crystal_qmmm()

    def test_ionic_crystal_qmmm_rejects_n_unit_cell_atoms(self):
        s = ORCAQMMMJobSettings(
            jobtype="IONIC-CRYSTAL-QMMM",
            ecp_layer_ecp="def2-ECP",
            n_unit_cell_atoms=5,
            low_level_method="ff.prms",
        )
        with pytest.raises(AssertionError, match="only applicable"):
            s.check_crystal_qmmm()

    def test_mol_crystal_qmmm_conv_charges_false_with_low_level_method_passes(
        self,
    ):
        """When conv_charges is False, a force-field file must be given
        via low_level_method; if it is, the assertion simply passes."""
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
            conv_charges=False,
        )
        s.check_crystal_qmmm()  # must not raise
        assert s.multiplicity == 0

    def test_non_crystal_jobtype_is_a_no_op(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        s.check_crystal_qmmm()
        assert s.multiplicity == 1


class TestFormattedPartitionStrings:
    def _settings(self):
        return ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")

    def test_none_returns_none(self):
        s = self._settings()
        assert s._get_formatted_partition_strings(None) is None

    def test_empty_list_returns_none(self):
        s = self._settings()
        assert s._get_formatted_partition_strings([]) is None

    def test_list_of_ints_compresses_ranges_and_shifts_to_zero_index(self):
        s = self._settings()
        result = s._get_formatted_partition_strings([1, 2, 3, 5, 7, 8, 9])
        assert result == "0:2 4 6:8"

    def test_string_range_is_parsed_and_shifted(self):
        s = self._settings()
        result = s._get_formatted_partition_strings("1-15,37,39")
        assert result == "0:14 36 38"

    def test_duplicates_are_deduplicated(self):
        s = self._settings()
        result = s._get_formatted_partition_strings([1, 1, 2, 2, 3])
        assert result == "0:2"

    def test_unsupported_type_raises_type_error(self):
        s = self._settings()
        with pytest.raises(TypeError):
            s._get_formatted_partition_strings(object())

    def test_fallback_parser_handles_dash_colon_and_plain_tokens(self):
        """`get_list_from_string_range` cannot parse ':'-style tokens, so
        mixing a ':' token into the string forces the whole call to raise
        and fall back to the local regex-based parser, which understands
        '-', ':', and bare-integer tokens."""
        s = self._settings()
        result = s._get_formatted_partition_strings("1-3,1:5,7")
        # union of {1,2,3} (dash), {1,2,3,4,5} (colon), {7} (plain) =
        # {1,2,3,4,5,7}; shifted to 0-index and compressed.
        assert result == "0:4 6"

    def test_min_atom_already_zero_indexed_is_not_shifted(self):
        s = self._settings()
        result = s._get_formatted_partition_strings([0, 1, 2])
        assert result == "0:2"


class TestQmmmBlockGeneration:
    def test_additive_qmmm_full_block(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            high_level_atoms=[1, 2, 3],
            low_level_method="system.prms",
            charge_total=0,
            mult_total=1,
            embedding_type="electronic",
            delete_la_double_counting=True,
        )
        block = s._write_qmmm_block()
        assert block.startswith("%qmmm\n")
        assert block.endswith("end\n")
        assert "QMAtoms {0:2} end" in block
        assert "Charge_Total 0" in block
        assert "Mult_Total 1" in block
        assert 'ORCAFFFilename "system.prms"' in block
        assert "DeleteLADoubleCounting true" in block
        # Electronic embedding is ORCA's default and is omitted from
        # the block rather than written explicitly (only mechanical
        # embedding needs an explicit "Embedding Mechanical" line).
        assert "Embedding" not in block

    def test_qmmm_block_property_delegates(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="system.prms",
            charge_total=0,
            mult_total=1,
        )
        assert s.qmmm_block == s._write_qmmm_block()

    def test_qmmm_requires_force_field_file(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", charge_total=0, mult_total=1)
        with pytest.raises(AssertionError, match="Force field file"):
            s._write_qmmm_block()

    def test_intermediate_atoms_uses_intermediate_charge_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=1,
            mult_intermediate=2,
        )
        block = s._write_qmmm_block()
        assert "Charge_Intermediate 1" in block
        assert "Mult_Intermediate 2" in block

    def test_none_charge_lines_are_omitted(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        block = s._write_qmmm_block()
        assert "Charge_Total None" not in block
        assert "None" not in block

    def test_custom_intermediate_method_file_written_when_not_builtin(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_method="my_custom_method",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=0,
            mult_intermediate=1,
        )
        block = s._write_qmmm_block()
        assert 'QM2CustomFile "my_custom_method"' in block

    def test_custom_functional_basis_written_when_no_builtin_method(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_functional="B3LYP",
            intermediate_level_basis="def2-SVP",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=0,
            mult_intermediate=1,
        )
        block = s._write_qmmm_block()
        assert 'QM2CUSTOMMETHOD "B3LYP"' in block
        assert 'QM2CUSTOMBASIS "def2-SVP"' in block

    def test_use_active_info_from_pbc_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            use_active_info_from_pbc=True,
        )
        block = s._write_qmmm_block()
        assert "Use_Active_InfoFromPDB true" in block

    def test_optregion_fixed_atoms_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            optregion_fixed_atoms=[1, 2, 3],
        )
        block = s._write_qmmm_block()
        assert "OptRegion_FixedAtoms {0:2} end" in block

    def test_crystal_subblock_included(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
            conv_charges=False,
        )
        block = s._write_qmmm_block()
        assert "Conv_Charges False" in block
        assert "NumUnitCellAtoms 12" in block

    def test_intermediate_solv_scheme_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            intermediate_solv_scheme="cds",
            charge_intermediate=0,
            mult_intermediate=1,
        )
        block = s._write_qmmm_block()
        assert "solv_scheme cds" in block

    def test_active_atoms_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            active_atoms=[1, 2, 3],
            charge_total=0,
            mult_total=1,
        )
        block = s._write_qmmm_block()
        assert "ActiveAtoms {0:2} end" in block

    def test_optregion_fixed_atoms_that_format_to_none_skips_line(self):
        """A truthy but content-free optregion_fixed_atoms value (e.g. a
        whitespace-only string) reaches the elif branch but formats down
        to None, so no OptRegion_FixedAtoms line is written."""
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            optregion_fixed_atoms=" ",
            charge_total=0,
            mult_total=1,
        )
        block = s._write_qmmm_block()
        assert "OptRegion_FixedAtoms" not in block

    def test_h_bond_length_dict_included_in_full_block(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            high_level_h_bond_length={("C", "H"): 1.09},
            charge_total=0,
            mult_total=1,
        )
        block = s._write_qmmm_block()
        assert "Dist_C_HLA 1.09" in block

    def test_empty_dict_h_bond_length_produces_falsy_block_and_is_skipped(
        self,
    ):
        """An empty dict is `is not None` (so the outer guard is entered)
        but `_get_h_bond_length` then returns an empty string for it,
        which is falsy, so nothing is appended to the qmmm block."""
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            high_level_h_bond_length={},
            charge_total=0,
            mult_total=1,
        )
        block = s._write_qmmm_block()
        assert "Dist_" not in block

    def test_delete_la_bond_double_counting_atoms_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            delete_la_bond_double_counting_atoms=True,
            charge_total=0,
            mult_total=1,
        )
        block = s._write_qmmm_block()
        assert "DeleteLABondDoubleCounting true" in block


class TestCrystalQmmmSubblockOptionalFields:
    def test_all_optional_fields_written(self):
        s = ORCAQMMMJobSettings(
            jobtype="IONIC-CRYSTAL-QMMM",
            low_level_method="ff.prms",
            ecp_layer_ecp="def2-ECP",
            ecp_layer=2,
            conv_charges_max_n_cycles=10,
            conv_charges_conv_thresh=1e-4,
            scale_formal_charge_mm_atom=0.5,
            scale_formal_charge_ecp_atom=0.75,
        )
        block = s._write_crystal_qmmm_subblock()
        assert "Conv_Charges_MaxNCycles 10" in block
        assert "Conv_Charges_ConvThresh 0.0001" in block
        assert "Scale_FormalCharge_MMAtom 0.5" in block
        assert "cECPs def2-ECP" in block
        assert "ECPLayers 2" in block
        assert "Scale_FormalCharge_ECPAtom 0.75" in block


class TestHBondLength:
    def _settings(self, **kwargs):
        return ORCAQMMMJobSettings(
            jobtype="QMMM", low_level_method="ff.prms", **kwargs
        )

    def test_dict_form(self):
        s = self._settings(
            high_level_h_bond_length={("C", "H"): 1.09, ("O", "H"): 0.98}
        )
        result = s._get_h_bond_length()
        assert "Dist_C_HLA 1.09" in result
        assert "Dist_O_HLA 0.98" in result

    def test_string_form_requires_existing_file(self, tmp_path):
        missing = tmp_path / "does_not_exist.txt"
        s = self._settings(high_level_h_bond_length=str(missing))
        with pytest.raises(AssertionError, match="does not exist"):
            s._get_h_bond_length()

    def test_string_form_with_existing_file(self, tmp_path):
        present = tmp_path / "dist.txt"
        present.write_text("Dist_C_HLA 1.09\n")
        s = self._settings(high_level_h_bond_length=str(present))
        result = s._get_h_bond_length()
        assert result == f'H_Dist_FileName "{present}"'

    def test_neither_dict_nor_str_returns_none_implicitly(self):
        """Only dict/str forms are handled; any other type (bypassing the
        `is not None` guard in `_write_qmmm_block` via a direct call)
        falls through to an implicit `None` return."""
        s = self._settings()
        s.high_level_h_bond_length = 123
        assert s._get_h_bond_length() is None


class TestEmbeddingType:
    def test_electronic(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            embedding_type="electronic",
        )
        # Electronic embedding is ORCA's default and is represented by
        # omitting the line (see reference fixtures under
        # tests/data/ORCATests/inputs/dna_qmmm*.inp, which never write
        # an "Embedding Electronic" line).
        assert s._get_embedding_type() is None

    def test_mechanical(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            embedding_type="Mechanical",
        )
        assert s._get_embedding_type() == "Embedding Mechanical"

    def test_invalid_raises(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            embedding_type="bogus",
        )
        with pytest.raises(ValueError, match="Invalid embedding type"):
            s._get_embedding_type()


class TestQMMMEquality:
    def test_equal_settings_are_equal(self):
        kwargs = dict(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        assert ORCAQMMMJobSettings(**kwargs) == ORCAQMMMJobSettings(**kwargs)

    def test_different_settings_are_not_equal(self):
        s1 = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        s2 = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=1,
            mult_total=2,
        )
        assert s1 != s2

    def test_different_type_returns_not_implemented(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        assert s.__eq__(5) is NotImplemented
        assert s != 5


class TestORCANEBJobSettings:
    def test_basic_construction(self):
        s = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
        )
        assert s.joboption == "NEB-TS"
        assert s.nimages == 8
        assert s.ending_xyzfile == "product.xyz"
        assert s.preopt_ends is False

    def test_route_string_with_dft(self):
        s = ORCANEBJobSettings(
            joboption="NEB-CI",
            functional="B3LYP",
            basis="def2-SVP",
        )
        route = s.route_string
        assert "B3LYP def2-SVP NEB-CI" in route

    def test_route_string_with_semiempirical(self):
        s = ORCANEBJobSettings(joboption="NEB-TS", semiempirical="XTB2")
        route = s.route_string
        assert "XTB2 NEB-TS" in route

    def test_route_string_with_freq(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            freq=True,
        )
        assert "Freq" in s.route_string

    def test_route_string_freq_and_numfreq_conflict_raises(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
        )
        s.freq = True
        s.numfreq = True
        with pytest.raises(ValueError, match="Cannot specify both"):
            s.route_string

    def test_route_string_appends_solvent_with_id(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="CPCM",
            solvent_id="water",
        )
        assert "CPCM(water)" in s.route_string

    def test_route_string_with_numfreq(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            numfreq=True,
        )
        assert "NumFreq" in s.route_string

    def test_route_string_scf_tol_and_algorithm(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            scf_tol="tight",
            scf_algorithm="KDIIS",
        )
        route = s.route_string
        assert "tightSCF" in route
        assert "KDIIS" in route

    def test_route_string_scf_tol_with_existing_scf_suffix_untouched(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            scf_tol="TightSCF",
        )
        assert "TightSCF" in s.route_string

    def test_route_string_custom_solvent_with_id(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            custom_solvent="Epsilon 16.7\n",
            solvent_id="chloroform",
        )
        assert "CPCM(chloroform)" in s.route_string

    def test_route_string_custom_solvent_without_id(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            custom_solvent="Epsilon 16.7\n",
        )
        assert s.route_string.strip().endswith("CPCM")

    def test_route_string_solvent_model_without_id(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="cpcm",
        )
        assert s.route_string.strip().endswith("CPCM")

    def test_route_string_solvent_id_without_model_defaults_to_cpcm(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            solvent_id="water",
        )
        assert "CPCM(water)" in s.route_string

    def test_route_string_solvent_model_deduplicated_when_repeated(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            semiempirical="CPCM",
            solvent_model="cpcm",
            solvent_id="water",
        )
        import re as _re

        route = s.route_string
        occurrences = len(_re.findall(r"\bcpcm\b", route, _re.IGNORECASE))
        assert occurrences == 1
        assert "CPCM(water)" in route

    def test_equal_settings_are_equal(self):
        kwargs = dict(joboption="NEB-TS", nimages=8)
        assert ORCANEBJobSettings(**kwargs) == ORCANEBJobSettings(**kwargs)

    def test_different_type_returns_not_implemented(self):
        s = ORCANEBJobSettings(joboption="NEB-TS", nimages=8)
        assert s.__eq__("not a settings object") is NotImplemented
        assert s != "not a settings object"

    def test_init_default(self):
        """Test default initialization."""
        settings = ORCANEBJobSettings()
        assert settings.joboption is None
        assert settings.nimages is None
        assert settings.preopt_ends is False

    def test_init_with_parameters(self):
        """Test initialization with parameters."""
        settings = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=8, semiempirical="XTB2"
        )
        assert settings.joboption == "NEB-TS"
        assert settings.nimages == 8
        assert settings.semiempirical == "XTB2"

    def test_route_string_generation(self):
        """Test route string generation."""
        settings = ORCANEBJobSettings(joboption="NEB-CI", semiempirical="XTB2")
        assert settings.route_string == "!  XTB2 NEB-CI"

    def test_neb_block_basic(self, tmpdir):
        """Test basic NEB block generation via ORCAInputWriter."""
        from unittest.mock import MagicMock

        settings = ORCANEBJobSettings(
            nimages=5, starting_xyz="start.xyz", ending_xyzfile="end.xyz"
        )
        # Create a mock job with the settings
        mock_job = MagicMock()
        mock_job.settings = settings

        writer = ORCAInputWriter(job=mock_job)
        neb_block = writer.neb_block
        assert "%NEB" in neb_block
        assert "NImages 5" in neb_block
        assert 'NEB_END_XYZFILE "end.xyz"' in neb_block

    def test_inheritance(self):
        """Test inheritance from ORCAJobSettings."""
        settings = ORCANEBJobSettings(functional="B3LYP", basis="def2-SVP")
        assert isinstance(settings, ORCANEBJobSettings)
        assert settings.functional == "B3LYP"
        assert settings.basis == "def2-SVP"

    def test_validation_errors(self):
        """Test validation raises appropriate errors."""
        from unittest.mock import MagicMock

        settings = ORCANEBJobSettings(starting_xyz="start.xyz")
        mock_job = MagicMock()
        mock_job.settings = settings
        writer = ORCAInputWriter(job=mock_job)

        with pytest.raises(
            AssertionError, match="The number of images is missing"
        ):
            _ = writer.neb_block

    def test_equality_identical_settings(self):
        """Test that identical NEB settings are equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts_guess.xyz",
            preopt_ends=True,
            semiempirical="XTB2",
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts_guess.xyz",
            preopt_ends=True,
            semiempirical="XTB2",
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        assert settings1 == settings2

    def test_equality_different_joboption(self):
        """Test that settings with different joboption are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=8, ending_xyzfile="product.xyz"
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-CI", nimages=8, ending_xyzfile="product.xyz"
        )
        assert settings1 != settings2

    def test_equality_different_nimages(self):
        """Test that settings with different nimages are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=8, ending_xyzfile="product.xyz"
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS", nimages=12, ending_xyzfile="product.xyz"
        )
        assert settings1 != settings2

    def test_equality_different_intermediate_xyzfile(self):
        """Test that settings with different intermediate file are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts1.xyz",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts2.xyz",
        )
        assert settings1 != settings2

    def test_equality_different_restarting_xyzfile(self):
        """Test that settings with different restart file are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            restarting_xyzfile="restart1.allxyz",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            restarting_xyzfile="restart2.allxyz",
        )
        assert settings1 != settings2

    def test_equality_different_preopt_ends(self):
        """Test that settings with different preopt_ends are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            preopt_ends=True,
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            preopt_ends=False,
        )
        assert settings1 != settings2

    def test_equality_different_semiempirical(self):
        """Test that settings with different semiempirical method are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            semiempirical="XTB2",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            semiempirical="XTB1",
        )
        assert settings1 != settings2

    def test_equality_different_parent_attributes(self):
        """Test that settings with different parent class attributes are not equal."""
        settings1 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            functional="B3LYP",
            basis="def2-SVP",
        )
        settings2 = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            functional="PBE0",
            basis="def2-SVP",
        )
        assert settings1 != settings2

    def test_equality_includes_all_neb_attributes(self):
        """Test that all 7 NEB-specific attributes are included in equality check."""
        # Create two settings that differ only in each NEB attribute
        base_kwargs = {
            "joboption": "NEB-TS",
            "nimages": 8,
            "ending_xyzfile": "product.xyz",
            "intermediate_xyzfile": "ts.xyz",
            "restarting_xyzfile": "restart.allxyz",
            "preopt_ends": True,
            "semiempirical": "XTB2",
        }

        # Test each attribute individually
        for attr in [
            "joboption",
            "nimages",
            "ending_xyzfile",
            "intermediate_xyzfile",
            "restarting_xyzfile",
            "preopt_ends",
            "semiempirical",
        ]:
            settings1 = ORCANEBJobSettings(**base_kwargs)
            modified_kwargs = base_kwargs.copy()

            # Change the attribute to a different value
            if attr == "joboption":
                modified_kwargs[attr] = "NEB-CI"
            elif attr == "nimages":
                modified_kwargs[attr] = 12
            elif attr == "preopt_ends":
                modified_kwargs[attr] = False
            elif attr == "semiempirical":
                modified_kwargs[attr] = "XTB1"
            else:  # file paths
                modified_kwargs[attr] = "different_file.xyz"

            settings2 = ORCANEBJobSettings(**modified_kwargs)
            assert (
                settings1 != settings2
            ), f"Equality failed for attribute: {attr}"


class TestORCAQMMMJobSettings:
    def test_partition_string_single_and_list_input(self):
        """Partition string should accept
        '1-15,37,39' or list and compress ranges."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings()
        s.high_level_atoms = "1-15,37,39"
        out = s._get_partition_string()
        assert out.strip() == "QMAtoms {0:14 36 38} end"

        # list input should produce the same output
        s.high_level_atoms = [1, *range(2, 16), 37, 39]
        out2 = s._get_partition_string()
        assert out2.strip() == "QMAtoms {0:14 36 38} end"

    def test_partition_string_qm_and_qm2(self):
        """When both high_level_atoms and medium_level_atoms
        provided, both lines should be returned."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings()
        s.high_level_atoms = "1-3,5"
        s.intermediate_level_atoms = "7-9,12"
        out = s._get_partition_string()
        # order: QMAtoms then QM2Atoms
        lines = [ln for ln in out.splitlines() if ln.strip()]
        assert lines[0].strip() == "QMAtoms {0:2 4} end"
        assert lines[1].strip() == "QM2Atoms {6:8 11} end"

    def test_charge_and_multiplicity_population(self):
        """ORCAQMMMJobSettings should populate .charge and
        .multiplicity from intermediate or high fields."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        # when both high and intermediate
        # specified, high-region takes precedence
        s1 = ORCAQMMMJobSettings(
            charge_intermediate=0,
            mult_intermediate=1,
            charge_high=2,
            mult_high=3,
        )
        assert s1.charge == 2
        assert s1.multiplicity == 3

        # intermediate missing -> fall back to high
        s2 = ORCAQMMMJobSettings(charge_high=-1, mult_high=2)
        assert s2.charge == -1
        assert s2.multiplicity == 2

        # re_init keeps high-first priority after overrides
        s1.charge_high = 5
        s1.mult_high = 1
        s1.re_init_and_validate()
        assert s1.charge == 5
        assert s1.multiplicity == 1

    def test_qmmm_xyz_line_uses_charge_high_not_total(self):
        from io import StringIO
        from unittest.mock import MagicMock

        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        settings = ORCAQMMMJobSettings(
            charge_high=1,
            mult_high=2,
            charge_total=0,
            mult_total=1,
        )
        job = MagicMock()
        job.settings = settings
        job.jobrunner = MagicMock()
        buf = StringIO()
        ORCAInputWriter(job=job)._write_charge_and_multiplicity(buf)
        assert buf.getvalue() == "* xyz 1 2\n"

    def test_qmmm_xyz_line_does_not_fall_back_to_total_charge(self):
        from io import StringIO
        from unittest.mock import MagicMock

        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        settings = ORCAQMMMJobSettings(
            charge_total=0,
            mult_total=1,
        )
        assert settings.charge == 0
        assert settings.charge_high is None
        job = MagicMock()
        job.settings = settings
        job.jobrunner = MagicMock()
        with pytest.raises(
            AssertionError, match="Charge and multiplicity must be specified"
        ):
            ORCAInputWriter(job=job)._write_charge_and_multiplicity(StringIO())

    def test_route_string_without_parent_jobtype(self):
        """Unset parent_jobtype should not crash route generation."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings(
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="XTB",
            jobtype="QMMM",
            charge_high=0,
            mult_high=1,
            charge_total=0,
            mult_total=1,
            high_level_atoms="1-5",
            parent_jobtype=None,
        )
        route = s._get_level_of_theory_string()
        assert "QMMM" in route
        assert "B3LYP" in route.upper() or "b3lyp" in route.lower()

    def test_partition_string_empty_and_none(self):
        """Empty string or None should return empty partition block."""
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings()
        s.high_level_atoms = ""
        assert s._get_partition_string() == ""

        s.high_level_atoms = None
        assert s._get_partition_string() == ""

    def test_qmmm_block_uses_orca_keywords(self):
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            high_level_functional="MN15",
            high_level_basis="def2-SVP",
            intermediate_level_functional="B3LYP",
            intermediate_level_basis="def2-SVP",
            high_level_atoms="1-3",
            charge_high=0,
            mult_high=2,
            charge_total=0,
            mult_total=2,
            high_level_h_bond_length={("C", "H"): 1.09, "N_H": 1.01},
            delete_la_double_counting=True,
            delete_la_bond_double_counting_atoms=True,
            embedding_type="electronic",
            parent_jobtype="opt",
        )
        block = s.qmmm_block
        assert 'QM2CUSTOMMETHOD "B3LYP"' in block
        assert 'QM2CUSTOMMETHOD "B3LYP" end' not in block
        assert 'QM2CUSTOMBASIS "def2-SVP"' in block
        assert 'QM2CUSTOMBASIS "def2-SVP" end' not in block
        assert "Dist_C_HLA 1.09" in block
        assert "Dist_N_HLA 1.01" in block
        assert "DeleteLADoubleCounting true" in block
        assert "Delete_LA_Double_Counting" not in block
        assert "Embedding" not in block
        assert "ORCAFFFilename" not in block
        route = s.qmmm_route_string
        assert "QM/QM2" in route

        custom_file = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            high_level_functional="MN15",
            high_level_basis="def2-SVP",
            intermediate_level_method="myQM2Method.txt",
            high_level_atoms="1-3",
            charge_high=0,
            mult_high=2,
            charge_total=0,
            mult_total=2,
        )
        custom_block = custom_file.qmmm_block
        assert 'QM2CustomFile "myQM2Method.txt"' in custom_block
        assert 'QM2CustomFile "myQM2Method.txt" end' not in custom_block

    def test_embedding_type_helpers(self, tmp_path):
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        mechanical = ORCAQMMMJobSettings(
            jobtype="QMMM",
            embedding_type="mechanical",
            high_level_atoms="1-3",
            charge_high=0,
            mult_high=1,
            charge_total=0,
            mult_total=1,
            low_level_method="system.ORCAFF.prms",
            high_level_h_bond_length={("C", "X"): 1.2, ("O", "HLA"): 0.98},
        )
        assert mechanical._get_embedding_type() == "Embedding Mechanical"
        block = mechanical.qmmm_block
        assert "Embedding Mechanical" in block
        assert "Dist_C_X 1.2" in block
        assert "Dist_O_HLA 0.98" in block

        electrostatic = ORCAQMMMJobSettings(embedding_type="electrostatic")
        assert electrostatic._get_embedding_type() is None

        invalid = ORCAQMMMJobSettings(embedding_type="invalid")
        with pytest.raises(ValueError, match="Invalid embedding type"):
            invalid._get_embedding_type()

        h_file = tmp_path / "QM_H_dist.txt"
        h_file.write_text("C 1.09\n")
        from_file = ORCAQMMMJobSettings(high_level_h_bond_length=str(h_file))
        assert from_file._get_h_bond_length() == (
            f'H_Dist_FileName "{h_file}"'
        )

        missing = ORCAQMMMJobSettings(
            high_level_h_bond_length=str(tmp_path / "missing.txt")
        )
        with pytest.raises(AssertionError, match="does not exist"):
            missing._get_h_bond_length()
