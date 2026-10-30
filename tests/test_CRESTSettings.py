import logging
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.settings.crest import CRESTProjectSettings


class TestCRESTJobSettingsConstructorValidation:
    def test_unknown_values_are_warned_but_stored_lowercased(self, caplog):
        with caplog.at_level(logging.WARNING):
            settings = CRESTJobSettings(
                gfn_version="GFN99",
                optimization_level="SUPERTIGHT",
                jobtype="BOGUS",
                solvent_model="MYSOLV",
                solvent_id="MYID",
                quick_mode="ULTRA",
            )
        assert settings.gfn_version == "gfn99"
        assert settings.optimization_level == "supertight"
        assert settings.jobtype == "bogus"
        assert settings.solvent_model == "mysolv"
        assert settings.solvent_id == "myid"
        assert settings.quick_mode == "ultra"
        assert "is not in the known CREST values" in caplog.text
        assert caplog.text.count("is not in the known CREST values") == 6

    def test_known_values_do_not_warn(self, caplog):
        with caplog.at_level(logging.WARNING):
            CRESTJobSettings(
                gfn_version="GFN2",
                optimization_level="TIGHT",
                jobtype="CONFORMERS",
                solvent_model="GBSA",
                solvent_id="water",
                quick_mode="QUICK",
            )
        assert "is not in the known CREST values" not in caplog.text

    def test_none_values_do_not_warn_and_are_not_lowercased(self, caplog):
        with caplog.at_level(logging.WARNING):
            settings = CRESTJobSettings()
        assert settings.optimization_level is None
        assert settings.jobtype is None
        assert settings.solvent_model is None
        assert settings.solvent_id is None
        assert settings.quick_mode is None
        assert "is not in the known CREST values" not in caplog.text

    def test_default_gfn_version_is_gfn2(self):
        assert CRESTJobSettings().gfn_version == "gfn2"

    def test_explicit_none_gfn_version_is_not_lowercased_or_warned(
        self, caplog
    ):
        with caplog.at_level(logging.WARNING):
            settings = CRESTJobSettings(gfn_version=None)
        assert settings.gfn_version is None
        assert "is not in the known CREST values" not in caplog.text

    def test_extra_kwargs_are_accepted_and_ignored(self):
        # **kwargs sink lets merge()/from_dict() round-trip a superset of
        # fields (e.g. from another settings class) without raising.
        settings = CRESTJobSettings(some_unrelated_field="x")
        assert not hasattr(settings, "some_unrelated_field")


class TestDefaultAndFromDict:
    def test_default_returns_settings_instance(self):
        settings = CRESTJobSettings.default()
        assert isinstance(settings, CRESTJobSettings)
        assert settings.charge is None

    def test_from_dict_constructs_settings(self):
        settings = CRESTJobSettings.from_dict(
            {"gfn_version": "gfnff", "charge": -1}
        )
        assert settings.gfn_version == "gfnff"
        assert settings.charge == -1


class TestCopyAndMerge:
    def test_copy_produces_independent_deep_copy(self):
        settings = CRESTJobSettings(constraints=["a"])
        copied = settings.copy()
        copied.constraints.append("b")
        assert settings.constraints == ["a"]
        assert copied.constraints == ["a", "b"]

    def test_merge_default_keywords_only_merges_charge_and_multiplicity(self):
        settings = CRESTJobSettings(
            charge=0, multiplicity=1, gfn_version="gfn2"
        )
        other = CRESTJobSettings(
            charge=-1, multiplicity=3, gfn_version="gfnff"
        )
        merged = settings.merge(other)
        assert merged.charge == -1
        assert merged.multiplicity == 3
        assert merged.gfn_version == "gfn2"

    def test_merge_with_explicit_keywords(self):
        settings = CRESTJobSettings(gfn_version="gfn2")
        other = CRESTJobSettings(gfn_version="gfnff", charge=5)
        merged = settings.merge(other, keywords=("gfn_version",))
        assert merged.gfn_version == "gfnff"
        assert merged.charge is None

    def test_merge_all_true_merges_every_field(self):
        settings = CRESTJobSettings(gfn_version="gfn2")
        other = CRESTJobSettings(gfn_version="gfnff", charge=5, nci=True)
        merged = settings.merge(other, merge_all=True)
        assert merged.gfn_version == "gfnff"
        assert merged.charge == 5
        assert merged.nci is True

    def test_merge_accepts_plain_dict_as_other(self):
        settings = CRESTJobSettings()
        merged = settings.merge({"charge": 7, "multiplicity": 2})
        assert merged.charge == 7
        assert merged.multiplicity == 2

    def test_merge_keywords_missing_from_other_are_skipped(self):
        settings = CRESTJobSettings(charge=1, multiplicity=1)
        merged = settings.merge({"charge": 9})
        assert merged.charge == 9
        assert merged.multiplicity == 1

    def test_merge_keywords_none_merges_every_field_like_merge_all(self):
        # keywords=None skips the keyword-filtering step entirely, so the
        # full other_dict is merged in -- same practical effect as
        # merge_all=True, but through the other branch of the `if
        # keywords is not None` check.
        settings = CRESTJobSettings(gfn_version="gfn2")
        other = CRESTJobSettings(gfn_version="gfnff", charge=5)
        merged = settings.merge(other, keywords=None)
        assert merged.gfn_version == "gfnff"
        assert merged.charge == 5


class TestSolventHelpers:
    def test_remove_solvent_clears_both_fields(self):
        settings = CRESTJobSettings(solvent_model="gbsa", solvent_id="water")
        settings.remove_solvent()
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_update_solvent_updates_and_lowercases_given_fields(self):
        settings = CRESTJobSettings()
        settings.update_solvent(solvent_model="GBSA", solvent_id="WATER")
        assert settings.solvent_model == "gbsa"
        assert settings.solvent_id == "water"

    def test_update_solvent_warns_on_unknown_values(self, caplog):
        settings = CRESTJobSettings()
        with caplog.at_level(logging.WARNING):
            settings.update_solvent(solvent_model="notreal")
        assert "is not in the known CREST values" in caplog.text
        assert settings.solvent_model == "notreal"

    def test_update_solvent_no_op_when_nothing_given(self):
        settings = CRESTJobSettings(solvent_model="gbsa", solvent_id="water")
        settings.update_solvent()
        assert settings.solvent_model == "gbsa"
        assert settings.solvent_id == "water"

    def test_modify_solvent_remove_flag_removes(self):
        settings = CRESTJobSettings(solvent_model="gbsa", solvent_id="water")
        settings.modify_solvent(remove_solvent=True)
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_modify_solvent_without_remove_updates(self):
        settings = CRESTJobSettings()
        settings.modify_solvent(
            remove_solvent=False, solvent_model="alpb", solvent_id="toluene"
        )
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "toluene"


class TestFromDatabase:
    def test_raises_when_file_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Database file not"):
            CRESTJobSettings.from_database(filepath=str(tmp_path / "nope.db"))

    def test_structure_id_and_record_selector_mutually_exclusive(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        with patch("chemsmart.database.database.Database"):
            with pytest.raises(ValueError, match="not both"):
                CRESTJobSettings.from_database(
                    filepath=str(db_file),
                    structure_id="abc123",
                    record_index=1,
                )

    def test_structure_id_not_found_raises(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        mock_db = MagicMock()
        mock_db.get_structure_by_partial_id.return_value = "full-sid"
        mock_db.get_structure.return_value = None
        with patch(
            "chemsmart.database.database.Database", return_value=mock_db
        ):
            with pytest.raises(ValueError, match="No structure found"):
                CRESTJobSettings.from_database(
                    filepath=str(db_file), structure_id="abc"
                )

    def test_structure_id_populates_charge_and_multiplicity_only(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        mock_db = MagicMock()
        mock_db.get_structure_by_partial_id.return_value = "full-sid"
        mock_db.get_structure.return_value = {"charge": 1, "multiplicity": 2}
        with patch(
            "chemsmart.database.database.Database", return_value=mock_db
        ):
            settings = CRESTJobSettings.from_database(
                filepath=str(db_file), structure_id="abc"
            )
        assert settings.charge == 1
        assert settings.multiplicity == 2
        assert "chemsmart database" in settings.title
        # structure_id path uses cls.default() and returns immediately,
        # without touching gfn_version/solvent from any record metadata.
        assert settings.gfn_version == "gfn2"

    def test_record_not_found_returns_none(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record", return_value=None
            ),
        ):
            result = CRESTJobSettings.from_database(
                filepath=str(db_file), record_index=1
            )
        assert result is None

    def test_record_selector_populates_settings_from_meta_and_structure(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {
            "meta": {
                "method": "GFNFF",
                "solvent_model": "gbsa",
                "solvent_id": "water",
            },
            "molecules": [{"charge": 0, "multiplicity": 1}],
        }
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            settings = CRESTJobSettings.from_database(
                filepath=str(db_file), record_index=1, structure_index="-1"
            )
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.gfn_version == "gfnff"
        assert settings.solvent_model == "gbsa"
        assert settings.solvent_id == "water"
        assert "chemsmart database" in settings.title

    def test_record_method_not_in_known_methods_leaves_default_gfn_version(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {
            "meta": {"method": "not-a-crest-method"},
            "molecules": [{"charge": 0, "multiplicity": 1}],
        }
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            settings = CRESTJobSettings.from_database(
                filepath=str(db_file), record_index=1
            )
        assert settings.gfn_version == "gfn2"

    def test_record_with_no_meta_or_molecules_uses_empty_structure(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {}
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            settings = CRESTJobSettings.from_database(
                filepath=str(db_file), record_index=1
            )
        assert settings.charge is None
        assert settings.multiplicity is None

    def test_record_with_slice_structure_index_raises(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {"meta": {}, "molecules": [{}]}
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            with pytest.raises(ValueError, match="one structure at a time"):
                CRESTJobSettings.from_database(
                    filepath=str(db_file),
                    record_index=1,
                    structure_index=":",
                )

    def test_record_structure_index_out_of_range_raises(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        fake_record = {"meta": {}, "molecules": [{"charge": 0}]}
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            with pytest.raises(ValueError, match="out of range"):
                CRESTJobSettings.from_database(
                    filepath=str(db_file),
                    record_index=1,
                    structure_index="5",
                )


class TestFromFilepath:
    def test_com_extension_dispatches_to_gaussian_input(self, tmp_path):
        filepath = tmp_path / "mol.com"
        filepath.write_text("")
        fake_settings = SimpleNamespace(**{"charge": 1, "multiplicity": 2})
        with patch("chemsmart.io.gaussian.input.Gaussian16Input") as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = CRESTJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 1
        assert settings.multiplicity == 2

    def test_gjf_extension_also_dispatches_to_gaussian_input(self, tmp_path):
        filepath = tmp_path / "mol.gjf"
        filepath.write_text("")
        fake_settings = SimpleNamespace(**{"charge": 0, "multiplicity": 1})
        with patch("chemsmart.io.gaussian.input.Gaussian16Input") as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = CRESTJobSettings.from_filepath(str(filepath))
        assert settings.charge == 0

    def test_log_extension_dispatches_to_gaussian_output(self, tmp_path):
        filepath = tmp_path / "mol.log"
        filepath.write_text("")
        fake_settings = SimpleNamespace(**{"charge": -1, "multiplicity": 1})
        with patch(
            "chemsmart.io.gaussian.output.Gaussian16Output"
        ) as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = CRESTJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == -1

    def test_inp_extension_dispatches_to_orca_input(self, tmp_path):
        filepath = tmp_path / "mol.inp"
        filepath.write_text("")
        fake_settings = SimpleNamespace(**{"charge": 2, "multiplicity": 3})
        with patch("chemsmart.io.orca.input.ORCAInput") as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = CRESTJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 2
        assert settings.multiplicity == 3

    def test_out_extension_gaussian_program_dispatches_to_gaussian_output(
        self, tmp_path
    ):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        fake_settings = SimpleNamespace(**{"charge": 0, "multiplicity": 1})
        with (
            patch(
                "chemsmart.utils.io.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.io.gaussian.output.Gaussian16Output") as mock_cls,
        ):
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = CRESTJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 0

    def test_out_extension_orca_program_dispatches_to_orca_output(
        self, tmp_path
    ):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        fake_settings = SimpleNamespace(**{"charge": 0, "multiplicity": 1})
        with (
            patch(
                "chemsmart.utils.io.get_program_type_from_file",
                return_value="orca",
            ),
            patch("chemsmart.io.orca.output.ORCAOutput") as mock_cls,
        ):
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = CRESTJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 0

    def test_out_extension_xtb_program_dispatches_to_xtb_settings(
        self, tmp_path
    ):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        fake_settings = SimpleNamespace(**{"charge": 0, "multiplicity": 1})
        with (
            patch(
                "chemsmart.utils.io.get_program_type_from_file",
                return_value="xtb",
            ),
            patch(
                "chemsmart.jobs.xtb.settings.XTBJobSettings.from_outfile",
                return_value=fake_settings,
            ) as mock_from_outfile,
        ):
            settings = CRESTJobSettings.from_filepath(str(filepath))
        mock_from_outfile.assert_called_once_with(str(filepath))
        assert settings.charge == 0

    def test_out_extension_unsupported_program_raises(self, tmp_path):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        with patch(
            "chemsmart.utils.io.get_program_type_from_file",
            return_value="mopac",
        ):
            with pytest.raises(ValueError, match="Unsupported .out file"):
                CRESTJobSettings.from_filepath(str(filepath))

    def test_unrecognized_extension_returns_defaults(self, tmp_path):
        filepath = tmp_path / "mol.xyz"
        filepath.write_text("")
        settings = CRESTJobSettings.from_filepath(str(filepath))
        assert isinstance(settings, CRESTJobSettings)
        assert settings.charge is None
        assert settings.gfn_version == "gfn2"


class TestEq:
    def test_equal_settings_are_equal(self):
        assert CRESTJobSettings(charge=0) == CRESTJobSettings(charge=0)

    def test_different_field_values_are_not_equal(self):
        assert CRESTJobSettings(charge=0) != CRESTJobSettings(charge=1)

    def test_different_type_returns_not_implemented(self):
        # __eq__ returning NotImplemented lets Python fall back to the
        # other type's comparison (or default identity inequality) rather
        # than crashing when compared against an unrelated type.
        assert CRESTJobSettings().__eq__(object()) is NotImplemented
        assert CRESTJobSettings() != object()


class TestCRESTJobSettings:
    def test_merge_dict(self):
        settings = CRESTJobSettings.default()
        settings.charge = None
        settings.multiplicity = None

        merged_settings = settings.merge({"charge": -2, "multiplicity": 3})
        assert merged_settings.charge == -2
        assert merged_settings.multiplicity == 3

    def test_get_settings_from_yaml(
        self, crest_yaml_settings_gas_project_name
    ):
        project_settings = CRESTProjectSettings.from_project(
            crest_yaml_settings_gas_project_name
        )
        conformers = project_settings.conformer_settings()

        assert isinstance(conformers, CRESTJobSettings)
        assert conformers.jobtype == "conformers"
        assert conformers.gfn_version == "gfn2"
        assert conformers.energy_window == 6.0
        assert conformers.optimization_level == "vtight"
        assert conformers.nci is False
        assert conformers.solvent_model is None
        assert conformers.solvent_id is None

    def test_get_settings_from_yaml_solv(
        self, crest_yaml_settings_solv_project_name
    ):
        project_settings = CRESTProjectSettings.from_project(
            crest_yaml_settings_solv_project_name
        )
        conformers = project_settings.conformer_settings()

        assert isinstance(conformers, CRESTJobSettings)
        assert conformers.jobtype == "conformers"
        assert conformers.gfn_version == "gfnff"
        assert conformers.optimization_level == "tight"
        assert conformers.solvent_model == "gbsa"
        assert conformers.solvent_id == "water"
