"""
Direct unit tests for ``chemsmart.jobs.orca.settings``.

Covers ``ORCAJobSettings`` (constructor validation, merge/copy/eq,
the various ``from_*file`` classmethods, ``route_string`` generation,
level-of-theory/solvent helpers), ``ORCApKaJobSettings`` (reference/
conjugate-base molecule construction and sub-job settings factories),
``ORCATSJobSettings``, and ``ORCAIRCJobSettings``.

Most of this is pure logic (string building, molecule bookkeeping) so
these tests construct real settings/molecule objects rather than
mocking. The ``from_comfile``/``from_logfile``/``from_inpfile``/
``from_outfile`` classmethods delegate to file-parser classes imported
locally inside each method; those parser classes are mocked since we
only care about how their results get merged into settings.
"""

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAJobSettings,
    ORCApKaJobSettings,
    ORCATSJobSettings,
)


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
