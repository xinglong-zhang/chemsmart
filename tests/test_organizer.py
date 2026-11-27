import os

import pytest

from chemsmart.io.organizer import FileOrganizer


class TestOrganizer:

    def test_file_organizer(self, tmpdir, excel_file):
        file_organizer = FileOrganizer(
            filename=excel_file, sheetname="co2", cols="B:D", skip=2, row=100
        )
        df = file_organizer.load_data()
        assert df.shape == (98, 3)
        folders, new_filenames, old_filenames = file_organizer.clean_data()
        assert len(folders) == len(new_filenames) == len(old_filenames) == 31


class TestFileOrganizerInit:
    def test_missing_filename_raises_assertion_error(self, tmp_path):
        with pytest.raises(AssertionError, match="does not exist"):
            FileOrganizer(
                directory=str(tmp_path),
                filename=str(tmp_path / "missing.xlsx"),
                sheetname="co2",
            )

    def test_default_directory_is_cwd(self, excel_file, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        file_organizer = FileOrganizer(filename=excel_file, sheetname="co2")
        assert file_organizer.directory == str(tmp_path)


class TestCreateDirectories:
    def test_creates_new_directories(self, excel_file, tmp_path):
        file_organizer = FileOrganizer(
            directory=str(tmp_path), filename=excel_file, sheetname="co2"
        )
        file_organizer.create_directories(["new_folder_a", "new_folder_b"])
        assert (tmp_path / "new_folder_a").is_dir()
        assert (tmp_path / "new_folder_b").is_dir()

    def test_skips_non_string_entries(self, excel_file, tmp_path):
        file_organizer = FileOrganizer(
            directory=str(tmp_path), filename=excel_file, sheetname="co2"
        )
        # Should not raise even with a non-string (e.g. NaN-like) entry.
        file_organizer.create_directories([float("nan"), "valid_folder"])
        assert (tmp_path / "valid_folder").is_dir()

    def test_skips_existing_relative_directory(
        self, excel_file, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        os.makedirs("already_here")
        file_organizer = FileOrganizer(
            directory=str(tmp_path), filename=excel_file, sheetname="co2"
        )
        # Should not raise; existing relative dir is skipped entirely.
        file_organizer.create_directories(["already_here"])
        assert (tmp_path / "already_here").is_dir()


class TestCopyAndRename:
    def test_copies_and_renames_found_file(self, mocker, excel_file, tmp_path):
        target_folder = "pathway1"
        os.makedirs(tmp_path / target_folder)
        source_file = tmp_path / "runtime_name.log"
        source_file.write_text("dummy log content")

        mocker.patch(
            "chemsmart.io.organizer.search_file",
            return_value=(str(source_file), str(tmp_path)),
        )

        file_organizer = FileOrganizer(
            directory=str(tmp_path), filename=excel_file, sheetname="co2"
        )
        file_organizer.copy_and_rename(
            target_folder, "final_name", "runtime_name"
        )

        copied_file = tmp_path / target_folder / "final_name.log"
        assert copied_file.exists()
        assert copied_file.read_text() == "dummy log content"

    def test_does_nothing_when_file_not_found(
        self, mocker, excel_file, tmp_path
    ):
        mocker.patch(
            "chemsmart.io.organizer.search_file",
            return_value=(None, None),
        )
        file_organizer = FileOrganizer(
            directory=str(tmp_path), filename=excel_file, sheetname="co2"
        )
        result = file_organizer.copy_and_rename(
            "pathway1", "final_name", "runtime_name"
        )
        assert result is None

    def test_does_not_double_append_extension(
        self, mocker, excel_file, tmp_path
    ):
        target_folder = "pathway1"
        os.makedirs(tmp_path / target_folder)
        source_file = tmp_path / "runtime_name.log"
        source_file.write_text("dummy log content")

        mock_search = mocker.patch(
            "chemsmart.io.organizer.search_file",
            return_value=(str(source_file), str(tmp_path)),
        )

        file_organizer = FileOrganizer(
            directory=str(tmp_path), filename=excel_file, sheetname="co2"
        )
        file_organizer.copy_and_rename(
            target_folder, "final_name.log", "runtime_name.log"
        )
        mock_search.assert_called_once_with("runtime_name.log")


class TestOrganizeFiles:
    def test_organize_files_end_to_end(self, mocker, excel_file, tmp_path):
        file_organizer = FileOrganizer(
            directory=str(tmp_path),
            filename=excel_file,
            sheetname="co2",
            cols="B:D",
            skip=2,
            row=100,
        )
        folders, new_filenames, old_filenames = file_organizer.clean_data()

        source_file = tmp_path / "source.log"
        source_file.write_text("dummy")
        mocker.patch(
            "chemsmart.io.organizer.search_file",
            return_value=(str(source_file), str(tmp_path)),
        )

        result = file_organizer.organize_files()
        assert result is None
        for folder in set(folders):
            if isinstance(folder, str):
                assert (tmp_path / folder).is_dir()
