"""
Direct unit tests for :class:`chemsmart.jobs.settings.MolecularJobSettings`
(the shared base for ``GaussianJobSettings``/``ORCAJobSettings``).

Solvent handling and custom-solvent-from-file are exercised via
``GaussianJobSettings`` as a concrete stand-in. ``from_database`` is
tested with the ``Database``/``resolve_record`` collaborators mocked
out, since its actual database-reading behavior is covered elsewhere
(test_io_database.py, test_database.py).
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.gaussian.settings import GaussianJobSettings


class TestCustomSolventValidation:
    def test_rejects_non_string_custom_solvent(self):
        with pytest.raises(ValueError, match="must be a string"):
            GaussianJobSettings(custom_solvent=123)

    def test_custom_solvent_free_text_stored_with_trailing_newline(self):
        settings = GaussianJobSettings(custom_solvent="eps=78.4")
        assert settings.custom_solvent == "eps=78.4\n"

    def test_custom_solvent_free_text_preserves_existing_newline(self):
        settings = GaussianJobSettings(custom_solvent="eps=78.4\n")
        assert settings.custom_solvent == "eps=78.4\n"

    def test_custom_solvent_from_file(self, tmp_path):
        solvent_file = tmp_path / "solvent.txt"
        solvent_file.write_text("eps=78.4\nrsolv=1.0\n")
        settings = GaussianJobSettings(custom_solvent=str(solvent_file))
        assert settings.custom_solvent == "eps=78.4\nrsolv=1.0\n"

    def test_set_custom_solvent_via_file_missing_raises(self):
        settings = GaussianJobSettings()
        with pytest.raises(ValueError, match="does not exist"):
            settings.set_custom_solvent_via_file("/no/such/file.txt")

    def test_custom_solvent_defaults_to_none(self):
        settings = GaussianJobSettings()
        assert settings.custom_solvent is None


class TestSolventHelpers:
    def test_remove_solvent_clears_all_solvent_fields(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.custom_solvent = "eps=78.4\n"
        settings.remove_solvent()
        assert settings.solvent_model is None
        assert settings.solvent_id is None
        assert settings.custom_solvent is None

    def test_update_solvent_only_updates_given_fields(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.update_solvent(solvent_id="toluene")
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_update_solvent_no_op_when_nothing_given(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.update_solvent()
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_modify_solvent_remove_flag_removes(self):
        settings = GaussianJobSettings(solvent_model="smd", solvent_id="water")
        settings.modify_solvent(remove_solvent=True)
        assert settings.solvent_model is None

    def test_modify_solvent_without_remove_updates(self):
        settings = GaussianJobSettings()
        settings.modify_solvent(
            remove_solvent=False, solvent_model="cpcm", solvent_id="dmso"
        )
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "dmso"


class TestFromDict:
    def test_from_dict_constructs_settings(self):
        settings = GaussianJobSettings.from_dict(
            {"functional": "b3lyp", "basis": "6-31g(d)"}
        )
        assert settings.functional == "b3lyp"
        assert settings.basis == "6-31g(d)"


class TestFromDatabase:
    def test_raises_when_file_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Database file not"):
            GaussianJobSettings.from_database(
                filepath=str(tmp_path / "nope.db")
            )

    def test_structure_id_and_record_selector_mutually_exclusive(
        self, tmp_path
    ):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        with patch("chemsmart.database.database.Database"):
            with pytest.raises(ValueError, match="not both"):
                GaussianJobSettings.from_database(
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
                GaussianJobSettings.from_database(
                    filepath=str(db_file), structure_id="abc"
                )

    def test_structure_id_populates_charge_and_multiplicity(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        mock_db = MagicMock()
        mock_db.get_structure_by_partial_id.return_value = "full-sid"
        mock_db.get_structure.return_value = {
            "charge": 1,
            "multiplicity": 2,
        }
        with patch(
            "chemsmart.database.database.Database", return_value=mock_db
        ):
            settings = GaussianJobSettings.from_database(
                filepath=str(db_file), structure_id="abc"
            )
        assert settings.charge == 1
        assert settings.multiplicity == 2
        assert "chemsmart database" in settings.title

    def test_record_not_found_returns_none(self, tmp_path):
        db_file = tmp_path / "test.db"
        db_file.write_text("")
        with (
            patch("chemsmart.database.database.Database"),
            patch(
                "chemsmart.database.utils.resolve_record", return_value=None
            ),
        ):
            result = GaussianJobSettings.from_database(
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
                "method": "b3lyp",
                "basis": "6-31g(d)",
                "jobtype": "opt",
                "solvent_model": "smd",
                "solvent_id": "water",
                "custom_solvent": None,
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
            settings = GaussianJobSettings.from_database(
                filepath=str(db_file), record_index=1, structure_index="-1"
            )

        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.functional == "b3lyp"
        assert settings.basis == "6-31g(d)"
        assert settings.jobtype == "opt"
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

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
                GaussianJobSettings.from_database(
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
                GaussianJobSettings.from_database(
                    filepath=str(db_file),
                    record_index=1,
                    structure_index="5",
                )
