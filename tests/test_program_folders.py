"""Tests for chemsmart.io.gaussian.folder and chemsmart.io.orca.folder."""

import pytest

from chemsmart.io.gaussian.folder import (
    GaussianInputFolder,
    GaussianOutputFolder,
)
from chemsmart.io.orca.folder import ORCAInputFolder, ORCAOutputFolder


@pytest.fixture
def make_files(tmp_path):
    def _make(*names):
        for name in names:
            (tmp_path / name).write_text("dummy content\n")
        return tmp_path

    return _make


class TestGaussianInputFolder:
    def test_all_com_files(self, make_files):
        folder_path = make_files("a.com", "b.com", "c.txt")
        folder = GaussianInputFolder(str(folder_path))
        assert sorted(folder.all_com_files) == sorted(
            str(folder_path / f) for f in ("a.com", "b.com")
        )

    def test_all_gjf_files(self, make_files):
        folder_path = make_files("a.gjf", "b.txt")
        folder = GaussianInputFolder(str(folder_path))
        assert folder.all_gjf_files == [str(folder_path / "a.gjf")]

    def test_all_input_files_combines_com_and_gjf(self, make_files):
        folder_path = make_files("a.com", "b.gjf")
        folder = GaussianInputFolder(str(folder_path))
        assert sorted(folder.all_input_files) == sorted(
            str(folder_path / f) for f in ("a.com", "b.gjf")
        )

    def test_all_com_files_in_current_folder(self, make_files):
        folder_path = make_files("a.com")
        folder = GaussianInputFolder(str(folder_path))
        assert folder.all_com_files_in_current_folder == [
            str(folder_path / "a.com")
        ]

    def test_all_gjf_files_in_current_folder(self, make_files):
        folder_path = make_files("a.gjf")
        folder = GaussianInputFolder(str(folder_path))
        assert folder.all_gjf_files_in_current_folder == [
            str(folder_path / "a.gjf")
        ]

    def test_all_input_files_in_current_folder(self, make_files):
        folder_path = make_files("a.com", "b.gjf")
        folder = GaussianInputFolder(str(folder_path))
        assert sorted(folder.all_input_files_in_current_folder) == sorted(
            str(folder_path / f) for f in ("a.com", "b.gjf")
        )


class TestGaussianOutputFolder:
    def test_all_log_files(self, make_files):
        folder_path = make_files("a.log", "b.txt")
        folder = GaussianOutputFolder(str(folder_path))
        assert folder.all_log_files == [str(folder_path / "a.log")]

    def test_all_out_files(self, make_files):
        folder_path = make_files("a.out")
        folder = GaussianOutputFolder(str(folder_path))
        assert folder.all_out_files == [str(folder_path / "a.out")]

    def test_all_log_files_in_current_folder(self, make_files):
        folder_path = make_files("a.log")
        folder = GaussianOutputFolder(str(folder_path))
        assert folder.all_log_files_in_current_folder == [
            str(folder_path / "a.log")
        ]

    def test_all_out_files_in_current_folder(self, make_files):
        folder_path = make_files("a.out")
        folder = GaussianOutputFolder(str(folder_path))
        assert folder.all_out_files_in_current_folder == [
            str(folder_path / "a.out")
        ]

    def test_all_output_files(self, mocker, make_files):
        folder_path = make_files("a.log")
        folder = GaussianOutputFolder(str(folder_path))
        mocker.patch.object(
            folder,
            "get_all_output_files_in_current_folder_and_subfolders_by_program",
            return_value=[str(folder_path / "a.log")],
        )
        assert folder.all_output_files == [str(folder_path / "a.log")]

    def test_all_output_files_in_current_folder(self, mocker, make_files):
        folder_path = make_files("a.log")
        folder = GaussianOutputFolder(str(folder_path))
        mocker.patch.object(
            folder,
            "get_all_output_files_in_current_folder_by_program",
            return_value=[str(folder_path / "a.log")],
        )
        assert folder.all_output_files_in_current_folder == [
            str(folder_path / "a.log")
        ]

    def test_all_molecules(self, mocker, make_files):
        folder_path = make_files("a.log")
        folder = GaussianOutputFolder(str(folder_path))
        mocker.patch.object(
            type(folder),
            "all_output_files",
            new_callable=mocker.PropertyMock,
            return_value=[str(folder_path / "a.log")],
        )
        mock_output_cls = mocker.patch(
            "chemsmart.io.gaussian.folder.Gaussian16Output"
        )
        mock_output_cls.return_value.molecule = "fake_molecule"
        assert folder.all_molecules == ["fake_molecule"]

    def test_total_service_units(self, mocker, make_files):
        folder_path = make_files("a.log", "b.log")
        folder = GaussianOutputFolder(str(folder_path))
        mocker.patch.object(
            type(folder),
            "all_output_files",
            new_callable=mocker.PropertyMock,
            return_value=[
                str(folder_path / "a.log"),
                str(folder_path / "b.log"),
            ],
        )
        mock_output_cls = mocker.patch(
            "chemsmart.io.gaussian.folder.Gaussian16Output"
        )
        mock_output_cls.return_value.total_core_hours = 1.5
        assert folder.total_service_units == 3.0

    def test_write_job_runtime(self, mocker, make_files, tmp_path):
        folder_path = make_files("a.log")
        folder = GaussianOutputFolder(str(folder_path))
        mocker.patch.object(
            type(folder),
            "all_output_files",
            new_callable=mocker.PropertyMock,
            return_value=[str(folder_path / "a.log")],
        )
        mock_output_cls = mocker.patch(
            "chemsmart.io.gaussian.folder.Gaussian16Output"
        )
        mock_output_cls.return_value.total_core_hours = 2.0

        runtime_file = tmp_path / "runtime.txt"
        folder.write_job_runtime(job_runtime_file=str(runtime_file))

        content = runtime_file.read_text()
        assert "Total time" in content
        assert "TOTAL core-hours" in content


class TestORCAInputFolder:
    def test_all_inp_files(self, make_files):
        folder_path = make_files("a.inp", "b.txt")
        folder = ORCAInputFolder(str(folder_path))
        assert folder.all_inp_files == [str(folder_path / "a.inp")]

    def test_all_input_files(self, make_files):
        folder_path = make_files("a.inp")
        folder = ORCAInputFolder(str(folder_path))
        assert folder.all_input_files == [str(folder_path / "a.inp")]

    def test_all_inp_files_in_current_folder(self, make_files):
        folder_path = make_files("a.inp")
        folder = ORCAInputFolder(str(folder_path))
        assert folder.all_inp_files_in_current_folder == [
            str(folder_path / "a.inp")
        ]

    def test_all_input_files_in_current_folder(self, make_files):
        folder_path = make_files("a.inp")
        folder = ORCAInputFolder(str(folder_path))
        assert folder.all_input_files_in_current_folder == [
            str(folder_path / "a.inp")
        ]


class TestORCAOutputFolder:
    def test_all_out_files(self, make_files):
        folder_path = make_files("a.out")
        folder = ORCAOutputFolder(str(folder_path))
        assert folder.all_out_files == [str(folder_path / "a.out")]

    def test_all_out_files_in_current_folder(self, make_files):
        folder_path = make_files("a.out")
        folder = ORCAOutputFolder(str(folder_path))
        assert folder.all_out_files_in_current_folder == [
            str(folder_path / "a.out")
        ]

    def test_all_output_files(self, mocker, make_files):
        folder_path = make_files("a.out")
        folder = ORCAOutputFolder(str(folder_path))
        mocker.patch.object(
            folder,
            "get_all_output_files_in_current_folder_and_subfolders_by_program",
            return_value=[str(folder_path / "a.out")],
        )
        assert folder.all_output_files == [str(folder_path / "a.out")]

    def test_all_output_files_in_current_folder(self, mocker, make_files):
        folder_path = make_files("a.out")
        folder = ORCAOutputFolder(str(folder_path))
        mocker.patch.object(
            folder,
            "get_all_output_files_in_current_folder_by_program",
            return_value=[str(folder_path / "a.out")],
        )
        assert folder.all_output_files_in_current_folder == [
            str(folder_path / "a.out")
        ]

    def test_total_service_units(self, mocker, make_files):
        folder_path = make_files("a.out", "b.out")
        folder = ORCAOutputFolder(str(folder_path))
        mocker.patch.object(
            type(folder),
            "all_output_files",
            new_callable=mocker.PropertyMock,
            return_value=[
                str(folder_path / "a.out"),
                str(folder_path / "b.out"),
            ],
        )
        mock_output_cls = mocker.patch("chemsmart.io.orca.folder.ORCAOutput")
        mock_output_cls.return_value.total_core_hours = 1.5
        assert folder.total_service_units == 3.0

    def test_write_job_runtime(self, mocker, make_files, tmp_path):
        folder_path = make_files("a.out")
        folder = ORCAOutputFolder(str(folder_path))
        mocker.patch.object(
            type(folder),
            "all_output_files",
            new_callable=mocker.PropertyMock,
            return_value=[str(folder_path / "a.out")],
        )
        mock_output_cls = mocker.patch("chemsmart.io.orca.folder.ORCAOutput")
        mock_output_cls.return_value.total_core_hours = 2.0

        runtime_file = tmp_path / "runtime.txt"
        folder.write_job_runtime(job_runtime_file=str(runtime_file))

        content = runtime_file.read_text()
        assert "Total time" in content
        assert "TOTAL core-hours" in content
