"""Tests for chemsmart.cli.database (assemble, export, inspect, query)."""

import importlib
from unittest.mock import MagicMock

import click
import pytest
from click.testing import CliRunner

from chemsmart.cli.database import database

# `chemsmart.cli.database`'s __init__.py rebinds the package attribute
# ``database`` to the Click group, shadowing the submodule of the same
# name. Look submodules up directly via the module cache so we can patch
# names inside them without hitting the Click group instead.
assemble_module = importlib.import_module("chemsmart.cli.database.assemble")
export_module = importlib.import_module("chemsmart.cli.database.export")
inspect_module = importlib.import_module("chemsmart.cli.database.inspect")
query_module = importlib.import_module("chemsmart.cli.database.query")


def invoke(args):
    runner = CliRunner()
    return runner.invoke(database, args, catch_exceptions=False)


def invoke_allow_exceptions(args):
    """Like invoke(), but for paths that raise a raw (non-Click) exception."""
    runner = CliRunner()
    return runner.invoke(database, args, catch_exceptions=True)


class TestAssembleCommand:
    def test_appends_db_extension_when_missing(self, mocker, tmp_path):
        mocker.patch.object(
            assemble_module,
            "BaseFolder",
            return_value=MagicMock(
                **{
                    "get_all_output_files_in_current_folder_and_subfolders_by_program.return_value": [
                        str(tmp_path / "a.log")
                    ]
                }
            ),
        )
        mock_assembler_cls = mocker.patch.object(
            assemble_module, "SingleFileAssembler"
        )
        mock_assembler_cls.return_value.assemble_data = {"foo": "bar"}
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        mock_db = mock_db_cls.return_value
        mock_db.insert_records.return_value = 1
        mock_db.count_records.return_value = 1

        result = invoke(
            [
                "assemble",
                "-d",
                str(tmp_path),
                "-p",
                "gaussian",
                "-o",
                "mydb",
            ]
        )
        assert result.exit_code == 0, result.output
        mock_db_cls.assert_called_once_with(db_file="mydb.db")

    def test_directory_defaults_to_cwd(self, mocker, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mock_folder_cls = mocker.patch.object(
            assemble_module,
            "BaseFolder",
            return_value=MagicMock(
                **{
                    "get_all_output_files_in_current_folder_and_subfolders_by_program.return_value": []
                }
            ),
        )
        result = invoke(["assemble"])
        assert result.exit_code == 0, result.output
        called_folder = mock_folder_cls.call_args[1]["folder"]
        assert called_folder == str(tmp_path)

    def test_nonexistent_directory_raises(self, tmp_path):
        result = invoke_allow_exceptions(
            ["assemble", "-d", str(tmp_path / "does_not_exist")]
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, FileNotFoundError)

    def test_unsupported_program_raises(self, tmp_path):
        result = invoke_allow_exceptions(
            ["assemble", "-d", str(tmp_path), "-p", "bogus"]
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)

    def test_no_program_uses_both_supported_programs(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = (
            []
        )
        mocker.patch.object(
            assemble_module, "BaseFolder", return_value=mock_folder
        )
        result = invoke(["assemble", "-d", str(tmp_path)])
        assert result.exit_code == 0, result.output
        assert (
            mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.call_count
            == 3
        )

    def test_no_files_found_returns_none(self, mocker, tmp_path):
        mocker.patch.object(
            assemble_module,
            "BaseFolder",
            return_value=MagicMock(
                **{
                    "get_all_output_files_in_current_folder_and_subfolders_by_program.return_value": []
                }
            ),
        )
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        result = invoke(["assemble", "-d", str(tmp_path), "-p", "orca"])
        assert result.exit_code == 0, result.output
        mock_db_cls.assert_not_called()

    def test_parse_failure_is_logged_and_skipped(self, mocker, tmp_path):
        mocker.patch.object(
            assemble_module,
            "BaseFolder",
            return_value=MagicMock(
                **{
                    "get_all_output_files_in_current_folder_and_subfolders_by_program.return_value": [
                        str(tmp_path / "bad.log")
                    ]
                }
            ),
        )
        mocker.patch.object(
            assemble_module,
            "SingleFileAssembler",
            side_effect=RuntimeError("parse failed"),
        )
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        result = invoke(["assemble", "-d", str(tmp_path), "-p", "gaussian"])
        assert result.exit_code == 0, result.output
        mock_db_cls.assert_not_called()

    def test_falsy_assemble_data_is_skipped(self, mocker, tmp_path):
        mocker.patch.object(
            assemble_module,
            "BaseFolder",
            return_value=MagicMock(
                **{
                    "get_all_output_files_in_current_folder_and_subfolders_by_program.return_value": [
                        str(tmp_path / "empty.log")
                    ]
                }
            ),
        )
        mock_assembler_cls = mocker.patch.object(
            assemble_module, "SingleFileAssembler"
        )
        mock_assembler_cls.return_value.assemble_data = None
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        result = invoke(["assemble", "-d", str(tmp_path), "-p", "gaussian"])
        assert result.exit_code == 0, result.output
        mock_db_cls.assert_not_called()

    def test_duplicate_records_logs_warning(self, mocker, tmp_path):
        mocker.patch.object(
            assemble_module,
            "BaseFolder",
            return_value=MagicMock(
                **{
                    "get_all_output_files_in_current_folder_and_subfolders_by_program.return_value": [
                        str(tmp_path / "a.log"),
                        str(tmp_path / "b.log"),
                    ]
                }
            ),
        )
        mock_assembler_cls = mocker.patch.object(
            assemble_module, "SingleFileAssembler"
        )
        mock_assembler_cls.return_value.assemble_data = {"foo": "bar"}
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        mock_db = mock_db_cls.return_value
        mock_db.insert_records.return_value = 2
        mock_db.count_records.return_value = 1

        result = invoke(["assemble", "-d", str(tmp_path), "-p", "gaussian"])
        assert result.exit_code == 0, result.output
        mock_db.create.assert_called_once()

    def test_include_failed_flag_forwarded(self, mocker, tmp_path):
        mocker.patch.object(
            assemble_module,
            "BaseFolder",
            return_value=MagicMock(
                **{
                    "get_all_output_files_in_current_folder_and_subfolders_by_program.return_value": [
                        str(tmp_path / "a.log")
                    ]
                }
            ),
        )
        mock_assembler_cls = mocker.patch.object(
            assemble_module, "SingleFileAssembler"
        )
        mock_assembler_cls.return_value.assemble_data = {"foo": "bar"}
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        mock_db_cls.return_value.insert_records.return_value = 1
        mock_db_cls.return_value.count_records.return_value = 1

        result = invoke(
            [
                "assemble",
                "-d",
                str(tmp_path),
                "-p",
                "gaussian",
                "--include-failed",
            ]
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_assembler_cls.call_args
        assert call_kwargs["include_failed"] is True

    def test_xtb_program_assembles_valid_calculation_directory(
        self, mocker, tmp_path
    ):
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = [
            str(tmp_path / "xtbdir" / "xtb.out")
        ]
        mocker.patch.object(
            assemble_module, "BaseFolder", return_value=mock_folder
        )
        mock_xtb_cls = mocker.patch.object(assemble_module, "XTBFolder")
        mock_xtb_cls.return_value.is_xtb_calculation_directory = True
        mock_assembler_cls = mocker.patch.object(
            assemble_module, "SingleFolderAssembler"
        )
        mock_assembler_cls.return_value.assemble_data = {"foo": "bar"}
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        mock_db_cls.return_value.insert_records.return_value = 1
        mock_db_cls.return_value.count_records.return_value = 1

        result = invoke(["assemble", "-d", str(tmp_path), "-p", "xtb"])

        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_assembler_cls.call_args
        assert call_kwargs["folder"] == str(tmp_path / "xtbdir")
        mock_db_cls.assert_called_once()

    def test_xtb_invalid_calculation_directory_is_skipped(
        self, mocker, tmp_path
    ):
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = [
            str(tmp_path / "xtbdir" / "xtb.out")
        ]
        mocker.patch.object(
            assemble_module, "BaseFolder", return_value=mock_folder
        )
        mock_xtb_cls = mocker.patch.object(assemble_module, "XTBFolder")
        type(mock_xtb_cls.return_value).is_xtb_calculation_directory = (
            mocker.PropertyMock(side_effect=ValueError("bad xtb dir"))
        )
        mock_db_cls = mocker.patch.object(assemble_module, "Database")

        result = invoke(["assemble", "-d", str(tmp_path), "-p", "xtb"])

        assert result.exit_code == 0, result.output
        mock_db_cls.assert_not_called()

    def test_xtb_program_no_valid_folders_returns_none(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = [
            str(tmp_path / "xtbdir" / "xtb.out")
        ]
        mocker.patch.object(
            assemble_module, "BaseFolder", return_value=mock_folder
        )
        mock_xtb_cls = mocker.patch.object(assemble_module, "XTBFolder")
        mock_xtb_cls.return_value.is_xtb_calculation_directory = False
        mock_db_cls = mocker.patch.object(assemble_module, "Database")

        result = invoke(["assemble", "-d", str(tmp_path), "-p", "xtb"])

        assert result.exit_code == 0, result.output
        mock_db_cls.assert_not_called()

    def test_xtb_folder_parse_failure_is_logged_and_skipped(
        self, mocker, tmp_path
    ):
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = [
            str(tmp_path / "xtbdir" / "xtb.out")
        ]
        mocker.patch.object(
            assemble_module, "BaseFolder", return_value=mock_folder
        )
        mock_xtb_cls = mocker.patch.object(assemble_module, "XTBFolder")
        mock_xtb_cls.return_value.is_xtb_calculation_directory = True
        mocker.patch.object(
            assemble_module,
            "SingleFolderAssembler",
            side_effect=RuntimeError("boom"),
        )
        mock_db_cls = mocker.patch.object(assemble_module, "Database")

        result = invoke(["assemble", "-d", str(tmp_path), "-p", "xtb"])

        assert result.exit_code == 0, result.output
        mock_db_cls.assert_not_called()

    def test_no_program_reports_mixed_files_and_folders(
        self, mocker, tmp_path
    ):
        """Covers the program-is-None log-message branch, exercising
        both the file-based and folder-based collection loops together."""

        def files_by_program(program):
            if program == "xtb":
                return [str(tmp_path / "xtbdir" / "xtb.out")]
            return [str(tmp_path / "a.log")]

        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.side_effect = (
            files_by_program
        )
        mocker.patch.object(
            assemble_module, "BaseFolder", return_value=mock_folder
        )
        mock_xtb_cls = mocker.patch.object(assemble_module, "XTBFolder")
        mock_xtb_cls.return_value.is_xtb_calculation_directory = True
        mock_file_assembler_cls = mocker.patch.object(
            assemble_module, "SingleFileAssembler"
        )
        mock_file_assembler_cls.return_value.assemble_data = {"a": 1}
        mock_folder_assembler_cls = mocker.patch.object(
            assemble_module, "SingleFolderAssembler"
        )
        mock_folder_assembler_cls.return_value.assemble_data = {"b": 2}
        mock_db_cls = mocker.patch.object(assemble_module, "Database")
        mock_db_cls.return_value.insert_records.return_value = 2
        mock_db_cls.return_value.count_records.return_value = 2

        result = invoke(["assemble", "-d", str(tmp_path)])

        assert result.exit_code == 0, result.output
        mock_file_assembler_cls.assert_called_once()
        mock_folder_assembler_cls.assert_called_once()


class TestExportCommandValidation:
    def _mock_db_checks(self, mocker, tmp_path, valid=True, schema_error=None):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        mocker.patch(
            "chemsmart.database.utils.is_chemsmart_database",
            return_value=valid,
        )
        if schema_error is not None:
            mocker.patch(
                "chemsmart.database.utils.check_schema_version",
                side_effect=schema_error,
            )
        else:
            mocker.patch(
                "chemsmart.database.utils.check_schema_version",
                return_value=None,
            )
        return str(db_file)

    def test_file_not_found_raises(self, tmp_path):
        result = invoke(
            [
                "export",
                "-f",
                str(tmp_path / "missing.db"),
                "-o",
                str(tmp_path / "out.json"),
            ]
        )
        assert result.exit_code != 0
        assert "not found" in result.output

    def test_invalid_database_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path, valid=False)
        result = invoke(
            ["export", "-f", db_file, "-o", str(tmp_path / "out.json")]
        )
        assert result.exit_code != 0
        assert "not a valid chemsmart database" in result.output

    def test_schema_version_error_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(
            mocker, tmp_path, schema_error=RuntimeError("bad schema")
        )
        result = invoke(
            ["export", "-f", db_file, "-o", str(tmp_path / "out.json")]
        )
        assert result.exit_code != 0
        assert "bad schema" in result.output

    def test_json_with_selector_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "--rid",
                "abc123",
                "-o",
                str(tmp_path / "out.json"),
            ]
        )
        assert result.exit_code != 0
        assert "not allowed" in result.output

    def test_json_with_keys_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "-k",
                "total_energy",
                "-o",
                str(tmp_path / "out.json"),
            ]
        )
        assert result.exit_code != 0
        assert "only valid for .csv" in result.output

    def test_csv_with_method_basis_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "-x",
                "MN15/def2tzvp",
                "-o",
                str(tmp_path / "out.csv"),
            ]
        )
        assert result.exit_code != 0
        assert "only valid for .xyz/.extxyz" in result.output

    def test_xyz_with_keys_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "-k",
                "total_energy",
                "--rid",
                "abc",
                "-o",
                str(tmp_path / "out.xyz"),
            ]
        )
        assert result.exit_code != 0
        assert "only valid for .csv" in result.output

    def test_xyz_without_selector_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            ["export", "-f", db_file, "-o", str(tmp_path / "out.xyz")]
        )
        assert result.exit_code != 0
        assert "requires exactly one" in result.output

    def test_xyz_multiple_selectors_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "--rid",
                "abc",
                "--sid",
                "def",
                "-o",
                str(tmp_path / "out.xyz"),
            ]
        )
        assert result.exit_code != 0
        assert "mutually" in result.output

    def test_xyz_structure_index_without_record_selector_raises(
        self, mocker, tmp_path
    ):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "--sid",
                "def",
                "--si",
                "2",
                "-o",
                str(tmp_path / "out.xyz"),
            ]
        )
        assert result.exit_code != 0
        assert "--si/--structure-index can only be used" in result.output

    def test_method_basis_without_sid_or_mid_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "--rid",
                "abc",
                "-x",
                "MN15/def2tzvp",
                "-o",
                str(tmp_path / "out.xyz"),
            ]
        )
        assert result.exit_code != 0
        assert "-x/--method-basis can only be used" in result.output

    def test_unsupported_output_extension_skips_cli_validation(
        self, mocker, tmp_path
    ):
        # An extension outside {.json, .csv, .xyz, .extxyz} matches
        # neither the `if`/`elif` branch in the CLI command, so none of
        # its selector/keys/method-basis validation runs; only
        # DatabaseExporter's own _infer_format() rejects it -- and that
        # constructor call isn't wrapped in the try/except ValueError
        # block (only exporter.export() is), so it surfaces as a raw
        # ValueError rather than a clean click.UsageError.
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke_allow_exceptions(
            [
                "export",
                "-f",
                db_file,
                "--rid",
                "abc",
                "-o",
                str(tmp_path / "out.txt"),
            ]
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "Unsupported output format" in str(result.exception)


class TestExportParseMethodBasis:
    def test_no_slash_raises(self, mocker, tmp_path):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        with pytest.raises(click.UsageError, match="method/basis"):
            export_module._parse_method_basis(str(db_file), "MN15only")

    def test_empty_parts_raises(self, tmp_path):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        with pytest.raises(click.UsageError, match="non-empty"):
            export_module._parse_method_basis(str(db_file), "MN15/")

    def test_unresolved_raises(self, mocker, tmp_path):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        mocker.patch(
            "chemsmart.database.database.Database.resolve_method_basis",
            return_value=None,
        )
        with pytest.raises(click.UsageError, match="not present"):
            export_module._parse_method_basis(str(db_file), "MN15/def2tzvp")

    def test_success_returns_tuple(self, mocker, tmp_path):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        mocker.patch(
            "chemsmart.database.database.Database.resolve_method_basis",
            return_value=("MN15", "def2tzvp"),
        )
        result = export_module._parse_method_basis(
            str(db_file), "MN15/def2tzvp"
        )
        assert result == ("MN15", "def2tzvp")


class TestExportCommandSuccess:
    def _mock_db_checks(self, mocker, tmp_path):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        mocker.patch(
            "chemsmart.database.utils.is_chemsmart_database",
            return_value=True,
        )
        mocker.patch(
            "chemsmart.database.utils.check_schema_version",
            return_value=None,
        )
        return str(db_file)

    def test_successful_json_export(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_exporter_cls = mocker.patch.object(
            export_module, "DatabaseExporter"
        )
        result = invoke(
            ["export", "-f", db_file, "-o", str(tmp_path / "out.json")]
        )
        assert result.exit_code == 0, result.output
        mock_exporter_cls.return_value.export.assert_called_once()

    def test_export_value_error_raises_click_exception(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_exporter_cls = mocker.patch.object(
            export_module, "DatabaseExporter"
        )
        mock_exporter_cls.return_value.export.side_effect = ValueError(
            "export failed"
        )
        result = invoke(
            ["export", "-f", db_file, "-o", str(tmp_path / "out.json")]
        )
        assert result.exit_code != 0
        assert "export failed" in result.output

    def test_successful_xyz_export_with_method_basis(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mocker.patch(
            "chemsmart.database.database.Database.resolve_method_basis",
            return_value=("MN15", "def2tzvp"),
        )
        mock_exporter_cls = mocker.patch.object(
            export_module, "DatabaseExporter"
        )
        result = invoke(
            [
                "export",
                "-f",
                db_file,
                "--sid",
                "abc123",
                "-x",
                "MN15/def2tzvp",
                "-o",
                str(tmp_path / "out.xyz"),
            ]
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_exporter_cls.call_args
        assert call_kwargs["method"] == "MN15"
        assert call_kwargs["basis"] == "def2tzvp"


class TestInspectCommand:
    def _mock_db_checks(self, mocker, tmp_path, valid=True, schema_error=None):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        mocker.patch(
            "chemsmart.database.utils.is_chemsmart_database",
            return_value=valid,
        )
        if schema_error is not None:
            mocker.patch(
                "chemsmart.database.utils.check_schema_version",
                side_effect=schema_error,
            )
        else:
            mocker.patch(
                "chemsmart.database.utils.check_schema_version",
                return_value=None,
            )
        return str(db_file)

    def test_file_not_found_raises(self, tmp_path):
        result = invoke(["inspect", "-f", str(tmp_path / "missing.db")])
        assert result.exit_code != 0
        assert "not found" in result.output

    def test_invalid_database_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path, valid=False)
        result = invoke(["inspect", "-f", db_file])
        assert result.exit_code != 0
        assert "not a valid chemsmart database" in result.output

    def test_schema_version_error_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(
            mocker, tmp_path, schema_error=RuntimeError("bad schema")
        )
        result = invoke(["inspect", "-f", db_file])
        assert result.exit_code != 0
        assert "bad schema" in result.output

    def test_multiple_entity_options_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(
            ["inspect", "-f", db_file, "--ri", "1", "--mid", "abc"]
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_structure_index_without_record_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(["inspect", "-f", db_file, "--si", "1"])
        assert result.exit_code != 0
        assert "requires --ri" in result.output

    def test_non_integer_structure_index_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(["inspect", "-f", db_file, "--ri", "1", "--si", "abc"])
        assert result.exit_code != 0
        assert "must be an integer index" in result.output

    def test_molecule_detail_path(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_inspector_cls = mocker.patch.object(
            inspect_module, "DatabaseInspector"
        )
        mock_inspector_cls.return_value.format_molecule_detail.return_value = (
            "MOLECULE DETAIL"
        )
        result = invoke(["inspect", "-f", db_file, "--mid", "abc123"])
        assert result.exit_code == 0, result.output
        assert "MOLECULE DETAIL" in result.output

    def test_structure_detail_standalone_path(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_inspector_cls = mocker.patch.object(
            inspect_module, "DatabaseInspector"
        )
        mock_inspector_cls.return_value.format_standalone_structure_detail.return_value = (
            "STRUCTURE DETAIL"
        )
        result = invoke(["inspect", "-f", db_file, "--sid", "def456"])
        assert result.exit_code == 0, result.output
        assert "STRUCTURE DETAIL" in result.output

    def test_overview_path(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_inspector_cls = mocker.patch.object(
            inspect_module, "DatabaseInspector"
        )
        mock_inspector_cls.return_value.format_overview.return_value = (
            "OVERVIEW"
        )
        result = invoke(["inspect", "-f", db_file])
        assert result.exit_code == 0, result.output
        assert "OVERVIEW" in result.output

    def test_record_detail_path_by_index(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_inspector_cls = mocker.patch.object(
            inspect_module, "DatabaseInspector"
        )
        mock_inspector_cls.return_value.format_record_detail.return_value = (
            "RECORD DETAIL"
        )
        result = invoke(["inspect", "-f", db_file, "--ri", "2"])
        assert result.exit_code == 0, result.output
        assert "RECORD DETAIL" in result.output

    def test_record_detail_path_by_id(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_inspector_cls = mocker.patch.object(
            inspect_module, "DatabaseInspector"
        )
        mock_inspector_cls.return_value.format_record_detail.return_value = (
            "RECORD DETAIL"
        )
        result = invoke(["inspect", "-f", db_file, "--rid", "abc123"])
        assert result.exit_code == 0, result.output
        assert "RECORD DETAIL" in result.output

    def test_structure_detail_path_by_record_index(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_inspector_cls = mocker.patch.object(
            inspect_module, "DatabaseInspector"
        )
        mock_inspector_cls.return_value.format_structure_detail.return_value = (
            "STRUCTURE IN RECORD"
        )
        result = invoke(["inspect", "-f", db_file, "--ri", "2", "--si", "3"])
        assert result.exit_code == 0, result.output
        assert "STRUCTURE IN RECORD" in result.output

    def test_structure_detail_path_by_record_id(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_inspector_cls = mocker.patch.object(
            inspect_module, "DatabaseInspector"
        )
        mock_inspector_cls.return_value.format_structure_detail.return_value = (
            "STRUCTURE IN RECORD"
        )
        result = invoke(
            ["inspect", "-f", db_file, "--rid", "abc123", "--si", "3"]
        )
        assert result.exit_code == 0, result.output
        assert "STRUCTURE IN RECORD" in result.output


class TestQueryCommand:
    def _mock_db_checks(self, mocker, tmp_path, valid=True, schema_error=None):
        db_file = tmp_path / "my.db"
        db_file.write_text("")
        mocker.patch(
            "chemsmart.database.utils.is_chemsmart_database",
            return_value=valid,
        )
        if schema_error is not None:
            mocker.patch(
                "chemsmart.database.utils.check_schema_version",
                side_effect=schema_error,
            )
        else:
            mocker.patch(
                "chemsmart.database.utils.check_schema_version",
                return_value=None,
            )
        return str(db_file)

    def test_file_not_found_raises(self, tmp_path):
        result = invoke_allow_exceptions(
            ["query", "-f", str(tmp_path / "missing.db")]
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, FileNotFoundError)

    def test_invalid_database_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path, valid=False)
        result = invoke(["query", "-f", db_file])
        assert result.exit_code != 0
        assert "not a valid chemsmart database" in result.output

    def test_schema_version_error_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(
            mocker, tmp_path, schema_error=RuntimeError("bad schema")
        )
        result = invoke(["query", "-f", db_file])
        assert result.exit_code != 0
        assert "bad schema" in result.output

    def test_non_positive_limit_raises(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        result = invoke(["query", "-f", db_file, "-l", "0"])
        assert result.exit_code != 0
        assert "positive integer" in result.output

    def test_invalid_query_logs_error_and_returns(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_query_cls = mocker.patch.object(query_module, "DatabaseQuery")
        mock_query_cls.return_value.query_summaries.side_effect = ValueError(
            "bad query"
        )
        result = invoke(["query", "-f", db_file, "-q", "bogus"])
        assert result.exit_code == 0, result.output

    def test_successful_query_prints_summary(self, mocker, tmp_path):
        db_file = self._mock_db_checks(mocker, tmp_path)
        mock_query_cls = mocker.patch.object(query_module, "DatabaseQuery")
        mock_dq = mock_query_cls.return_value
        mock_dq.query_summaries.return_value = ["summary1"]
        mock_dq._config = {"entity_name": "record"}
        mock_dq.count_total.return_value = 5
        mock_dq.format_summary.return_value = "FORMATTED SUMMARY"
        result = invoke(["query", "-f", db_file])
        assert result.exit_code == 0, result.output
        assert "FORMATTED SUMMARY" in result.output
