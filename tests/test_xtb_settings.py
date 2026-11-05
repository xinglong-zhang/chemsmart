import logging
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.xtb.settings import XTBJobSettings


class TestXTBJobSettingsConstructorValidation:
    def test_unknown_values_are_warned_but_stored_lowercased(self, caplog):
        with caplog.at_level(logging.WARNING):
            settings = XTBJobSettings(
                gfn_version="GFN99",
                optimization_level="SUPERTIGHT",
                jobtype="BOGUS",
                solvent_model="MYSOLV",
                solvent_id="MYID",
            )
        assert settings.gfn_version == "gfn99"
        assert settings.optimization_level == "supertight"
        assert settings.jobtype == "bogus"
        assert settings.solvent_model == "mysolv"
        assert settings.solvent_id == "myid"
        assert "is not in the known xTB values" in caplog.text
        assert caplog.text.count("is not in the known xTB values") == 5

    def test_known_values_do_not_warn(self, caplog):
        with caplog.at_level(logging.WARNING):
            XTBJobSettings(
                gfn_version="GFN2",
                optimization_level="TIGHT",
                jobtype="OPT",
                solvent_model="GBSA",
                solvent_id="water",
            )
        assert "is not in the known xTB values" not in caplog.text

    def test_none_values_do_not_warn_and_are_not_lowercased(self, caplog):
        with caplog.at_level(logging.WARNING):
            settings = XTBJobSettings()
        assert settings.optimization_level is None
        assert settings.jobtype is None
        assert settings.solvent_model is None
        assert settings.solvent_id is None
        assert "is not in the known xTB values" not in caplog.text

    def test_default_gfn_version_is_gfn2(self):
        assert XTBJobSettings().gfn_version == "gfn2"

    def test_explicit_none_gfn_version_is_not_lowercased_or_warned(
        self, caplog
    ):
        with caplog.at_level(logging.WARNING):
            settings = XTBJobSettings(gfn_version=None)
        assert settings.gfn_version is None
        assert "is not in the known xTB values" not in caplog.text

    def test_defaults_for_freq_grad_and_title(self):
        settings = XTBJobSettings()
        assert settings.freq is False
        assert settings.grad is False
        assert settings.title is None

    def test_extra_kwargs_are_accepted_and_ignored(self):
        # **kwargs sink lets merge()/from_dict() round-trip a superset of
        # fields (e.g. from another settings class) without raising.
        settings = XTBJobSettings(some_unrelated_field="x")
        assert not hasattr(settings, "some_unrelated_field")


class TestDefaultAndFromDict:
    def test_default_returns_settings_instance(self):
        settings = XTBJobSettings.default()
        assert isinstance(settings, XTBJobSettings)
        assert settings.charge is None

    def test_from_dict_constructs_settings(self):
        settings = XTBJobSettings.from_dict(
            {"gfn_version": "gfnff", "charge": -1}
        )
        assert settings.gfn_version == "gfnff"
        assert settings.charge == -1


class TestCopyAndMerge:
    def test_copy_produces_independent_deep_copy(self):
        settings = XTBJobSettings(constraints=["a"])
        copied = settings.copy()
        copied.constraints.append("b")
        assert settings.constraints == ["a"]
        assert copied.constraints == ["a", "b"]

    def test_merge_default_keywords_merge_charge_multiplicity_and_title(self):
        settings = XTBJobSettings(
            charge=0, multiplicity=1, title="old", gfn_version="gfn2"
        )
        other = XTBJobSettings(
            charge=-1, multiplicity=3, title="new", gfn_version="gfnff"
        )
        merged = settings.merge(other)
        assert merged.charge == -1
        assert merged.multiplicity == 3
        assert merged.title == "new"
        assert merged.gfn_version == "gfn2"

    def test_merge_with_explicit_keywords(self):
        settings = XTBJobSettings(gfn_version="gfn2")
        other = XTBJobSettings(gfn_version="gfnff", charge=5)
        merged = settings.merge(other, keywords=("gfn_version",))
        assert merged.gfn_version == "gfnff"
        assert merged.charge is None

    def test_merge_keywords_none_merges_every_field_like_merge_all(self):
        settings = XTBJobSettings(gfn_version="gfn2")
        other = XTBJobSettings(gfn_version="gfnff", charge=5)
        merged = settings.merge(other, keywords=None)
        assert merged.gfn_version == "gfnff"
        assert merged.charge == 5

    def test_merge_all_true_merges_every_field(self):
        settings = XTBJobSettings(gfn_version="gfn2")
        other = XTBJobSettings(gfn_version="gfnff", charge=5, freq=True)
        merged = settings.merge(other, merge_all=True)
        assert merged.gfn_version == "gfnff"
        assert merged.charge == 5
        assert merged.freq is True

    def test_merge_accepts_plain_dict_as_other(self):
        settings = XTBJobSettings()
        merged = settings.merge({"charge": 7, "multiplicity": 2})
        assert merged.charge == 7
        assert merged.multiplicity == 2

    def test_merge_keywords_missing_from_other_are_skipped(self):
        settings = XTBJobSettings(charge=1, multiplicity=1)
        merged = settings.merge({"charge": 9})
        assert merged.charge == 9
        assert merged.multiplicity == 1


class TestSolventHelpers:
    def test_remove_solvent_clears_both_fields(self):
        settings = XTBJobSettings(solvent_model="gbsa", solvent_id="water")
        settings.remove_solvent()
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_update_solvent_updates_and_lowercases_given_fields(self):
        settings = XTBJobSettings()
        settings.update_solvent(solvent_model="GBSA", solvent_id="WATER")
        assert settings.solvent_model == "gbsa"
        assert settings.solvent_id == "water"

    def test_update_solvent_warns_on_unknown_values(self, caplog):
        settings = XTBJobSettings()
        with caplog.at_level(logging.WARNING):
            settings.update_solvent(solvent_model="notreal")
        assert "is not in the known xTB values" in caplog.text
        assert settings.solvent_model == "notreal"

    def test_update_solvent_no_op_when_nothing_given(self):
        settings = XTBJobSettings(solvent_model="gbsa", solvent_id="water")
        settings.update_solvent()
        assert settings.solvent_model == "gbsa"
        assert settings.solvent_id == "water"

    def test_modify_solvent_remove_flag_removes(self):
        settings = XTBJobSettings(solvent_model="gbsa", solvent_id="water")
        settings.modify_solvent(remove_solvent=True)
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_modify_solvent_without_remove_updates(self):
        settings = XTBJobSettings()
        settings.modify_solvent(
            remove_solvent=False, solvent_model="alpb", solvent_id="toluene"
        )
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "toluene"


class TestFromDatabase:
    def test_raises_when_file_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Database file not"):
            XTBJobSettings.from_database(filepath=str(tmp_path / "nope.db"))

    def test_structure_id_and_record_selector_mutually_exclusive(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        with patch("chemsmart.database.database.Database"):
            with pytest.raises(ValueError, match="not both"):
                XTBJobSettings.from_database(
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
                XTBJobSettings.from_database(
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
            settings = XTBJobSettings.from_database(
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
            result = XTBJobSettings.from_database(
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
            settings = XTBJobSettings.from_database(
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
            "meta": {"method": "not-an-xtb-method"},
            "molecules": [{"charge": 0, "multiplicity": 1}],
        }
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record",
                return_value=fake_record,
            ),
        ):
            settings = XTBJobSettings.from_database(
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
            settings = XTBJobSettings.from_database(
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
                XTBJobSettings.from_database(
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
                XTBJobSettings.from_database(
                    filepath=str(db_file),
                    record_index=1,
                    structure_index="5",
                )


class TestFromOutfile:
    def test_from_outfile_maps_all_fields_from_xtb_main_out(self, tmp_path):
        filepath = tmp_path / "job.out"
        filepath.write_text("")
        fake_output = MagicMock()
        fake_output.method = "gfn2"
        fake_output.optimization_level = "tight"
        fake_output.net_charge = 0
        fake_output.multiplicity = 1
        fake_output.jobtype = "opt"
        fake_output.route_object.freq = True
        fake_output.route_object.grad = False
        fake_output.solvent_model = "gbsa"
        fake_output.solvent_id = "water"
        with patch(
            "chemsmart.io.xtb.file.XTBMainOut", return_value=fake_output
        ) as mock_cls:
            settings = XTBJobSettings.from_outfile(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.gfn_version == "gfn2"
        assert settings.optimization_level == "tight"
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.jobtype == "opt"
        assert settings.freq is True
        assert settings.grad is False
        assert settings.solvent_model == "gbsa"
        assert settings.solvent_id == "water"
        assert "chemsmart database" not in settings.title
        assert "job.out" in settings.title


class TestFromFilepath:
    def test_com_extension_dispatches_to_gaussian_input(self, tmp_path):
        filepath = tmp_path / "mol.com"
        filepath.write_text("")
        fake_settings = SimpleNamespace(charge=1, multiplicity=2, title="t")
        with patch("chemsmart.io.gaussian.input.Gaussian16Input") as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = XTBJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 1
        assert settings.multiplicity == 2
        assert settings.title == "t"

    def test_gjf_extension_also_dispatches_to_gaussian_input(self, tmp_path):
        filepath = tmp_path / "mol.gjf"
        filepath.write_text("")
        fake_settings = SimpleNamespace(charge=0, multiplicity=1)
        with patch("chemsmart.io.gaussian.input.Gaussian16Input") as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = XTBJobSettings.from_filepath(str(filepath))
        assert settings.charge == 0

    def test_log_extension_dispatches_to_gaussian_output(self, tmp_path):
        filepath = tmp_path / "mol.log"
        filepath.write_text("")
        fake_settings = SimpleNamespace(charge=-1, multiplicity=1)
        with patch(
            "chemsmart.io.gaussian.output.Gaussian16Output"
        ) as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = XTBJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == -1

    def test_inp_extension_dispatches_to_orca_input(self, tmp_path):
        filepath = tmp_path / "mol.inp"
        filepath.write_text("")
        fake_settings = SimpleNamespace(charge=2, multiplicity=3)
        with patch("chemsmart.io.orca.input.ORCAInput") as mock_cls:
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = XTBJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 2
        assert settings.multiplicity == 3

    def test_out_extension_xtb_program_dispatches_to_from_outfile(
        self, tmp_path
    ):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        fake_settings = XTBJobSettings(charge=0, multiplicity=1)
        with (
            patch(
                "chemsmart.utils.io.get_program_type_from_file",
                return_value="xtb",
            ),
            patch.object(
                XTBJobSettings,
                "from_outfile",
                return_value=fake_settings,
            ) as mock_from_outfile,
        ):
            settings = XTBJobSettings.from_filepath(str(filepath))
        mock_from_outfile.assert_called_once_with(str(filepath))
        assert settings is fake_settings

    def test_out_extension_gaussian_program_dispatches_to_gaussian_output(
        self, tmp_path
    ):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        fake_settings = SimpleNamespace(charge=0, multiplicity=1)
        with (
            patch(
                "chemsmart.utils.io.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.io.gaussian.output.Gaussian16Output") as mock_cls,
        ):
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = XTBJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 0

    def test_out_extension_orca_program_dispatches_to_orca_output(
        self, tmp_path
    ):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        fake_settings = SimpleNamespace(charge=0, multiplicity=1)
        with (
            patch(
                "chemsmart.utils.io.get_program_type_from_file",
                return_value="orca",
            ),
            patch("chemsmart.io.orca.output.ORCAOutput") as mock_cls,
        ):
            mock_cls.return_value.read_settings.return_value = fake_settings
            settings = XTBJobSettings.from_filepath(str(filepath))
        mock_cls.assert_called_once_with(filename=str(filepath))
        assert settings.charge == 0

    def test_out_extension_unsupported_program_raises(self, tmp_path):
        filepath = tmp_path / "mol.out"
        filepath.write_text("")
        with patch(
            "chemsmart.utils.io.get_program_type_from_file",
            return_value="mopac",
        ):
            with pytest.raises(ValueError, match="Unsupported .out file"):
                XTBJobSettings.from_filepath(str(filepath))

    def test_unrecognized_extension_returns_defaults(self, tmp_path):
        filepath = tmp_path / "mol.xyz"
        filepath.write_text("")
        settings = XTBJobSettings.from_filepath(str(filepath))
        assert isinstance(settings, XTBJobSettings)
        assert settings.charge is None
        assert settings.gfn_version == "gfn2"


class TestEq:
    def test_equal_settings_are_equal(self):
        assert XTBJobSettings(charge=0) == XTBJobSettings(charge=0)

    def test_different_field_values_are_not_equal(self):
        assert XTBJobSettings(charge=0) != XTBJobSettings(charge=1)

    def test_different_type_returns_not_implemented(self):
        assert XTBJobSettings().__eq__(object()) is NotImplemented
        assert XTBJobSettings() != object()
