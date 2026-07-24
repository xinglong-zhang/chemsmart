"""
Tests for the standalone CLI scripts under ``chemsmart/scripts``.

These scripts are thin click-command wrappers around already-tested
library classes (``CubeFileOperator``, ``FileConverter``, ``FileOrganizer``,
output parsers, etc). Each test mocks the underlying class so the CLI
argument-parsing and dispatch logic can be exercised without touching
real files or performing real computation.
"""

from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner


class TestCubeOperationScript:
    def test_entry_point_invokes_cube_operator(self):
        from chemsmart.scripts.cube_operation import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.cube_operation.CubeFileOperator"
        ) as mock_cls:
            mock_operator = MagicMock()
            mock_cls.return_value = mock_operator
            result = runner.invoke(
                entry_point,
                [
                    "-c1",
                    "a.cube",
                    "-c2",
                    "b.cube",
                    "-x",
                    "add",
                    "-o",
                    "out.cube",
                ],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_cls.assert_called_once_with(
            cubefile1="a.cube",
            cubefile2="b.cube",
            operation="add",
            output_cubefile="out.cube",
        )
        mock_operator.write_results.assert_called_once()

    def test_requires_both_cube_files(self):
        from chemsmart.scripts.cube_operation import entry_point

        runner = CliRunner()
        result = runner.invoke(entry_point, ["-c1", "a.cube"])
        assert result.exit_code != 0


class TestFileConverterScript:
    def test_entry_point_converts_single_file(self):
        from chemsmart.scripts.file_converter import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.file_converter.FileConverter"
        ) as mock_cls:
            mock_converter = MagicMock()
            mock_cls.return_value = mock_converter
            result = runner.invoke(
                entry_point,
                ["-f", "test.xyz", "-o", "sdf"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_cls.assert_called_once_with(
            filename="test.xyz",
            output_filetype="sdf",
            include_intermediate_structures=False,
        )
        mock_converter.convert_files.assert_called_once()

    def test_entry_point_converts_directory(self):
        from chemsmart.scripts.file_converter import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.file_converter.FileConverter"
        ) as mock_cls:
            mock_converter = MagicMock()
            mock_cls.return_value = mock_converter
            result = runner.invoke(
                entry_point,
                ["-d", "some_dir", "-t", "log"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_cls.assert_called_once_with(
            directory="some_dir",
            type="log",
            program=None,
            output_filetype="xyz",
            include_intermediate_structures=False,
        )
        mock_converter.convert_files.assert_called_once()

    def test_directory_requires_filetype(self):
        from chemsmart.scripts.file_converter import entry_point

        runner = CliRunner()
        with pytest.raises(AssertionError, match="filetype"):
            runner.invoke(
                entry_point,
                ["-d", "some_dir"],
                catch_exceptions=False,
            )

    def test_requires_filename_when_no_directory(self):
        from chemsmart.scripts.file_converter import entry_point

        runner = CliRunner()
        with pytest.raises(AssertionError, match="Filename must be"):
            runner.invoke(entry_point, [], catch_exceptions=False)


class TestFileOrganizerScript:
    def test_entry_point_invokes_organizer(self):
        from chemsmart.scripts.file_organizer import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.file_organizer.FileOrganizer"
        ) as mock_cls:
            mock_organizer = MagicMock()
            mock_cls.return_value = mock_organizer
            result = runner.invoke(
                entry_point,
                ["-f", "data.xlsx", "-n", "sheet1"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_cls.assert_called_once_with(
            directory=".",
            filename="data.xlsx",
            sheetname="sheet1",
            filetype="log",
            cols=None,
            skip=2,
            row=100,
            keep_default_na=False,
        )
        mock_organizer.organize_files.assert_called_once()

    def test_requires_filename_and_name(self):
        from chemsmart.scripts.file_organizer import entry_point

        runner = CliRunner()
        result = runner.invoke(entry_point, [])
        assert result.exit_code != 0


class TestSubmitJobsScript:
    def test_entry_point_submits_each_file(self, tmp_path):
        from chemsmart.scripts.submit_jobs import entry_point

        list_file = tmp_path / "jobs.txt"
        list_file.write_text("job1.com\njob2.com\n")

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.submit_jobs.subprocess.Popen"
        ) as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"", b"")
            mock_popen.return_value = mock_process
            result = runner.invoke(
                entry_point,
                [
                    "-f",
                    str(list_file),
                    "-c",
                    "chemsmart sub gaussian -p test -f file sp",
                ],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock_popen.call_count == 2
        called_cmds = [c.args[0] for c in mock_popen.call_args_list]
        assert (
            called_cmds[0][-2:] == ["-f", "job1.com"]
            or "job1.com" in called_cmds[0]
        )
        assert any("job2.com" in c for c in called_cmds[1])


class TestWbiAnalysisScript:
    def test_entry_point_logs_natural_charges(self):
        from chemsmart.scripts.wbi_analysis import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.wbi_analysis.Gaussian16WBIOutput"
        ) as mock_cls:
            mock_output = MagicMock()
            mock_output.natural_charges = {"1Fe": -0.5, "2O": 0.3}
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.log", "-n", "1"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_cls.assert_called_once_with(filename="some.log")


class TestWriteXyzScript:
    def test_single_structure_writes_single_file(self, tmp_path):
        from chemsmart.scripts.write_xyz import entry_point

        input_file = tmp_path / "mymol.xyz"
        input_file.write_text("dummy")

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.write_xyz.Molecule.from_filepath"
        ) as mock_from_filepath:
            mol = MagicMock()
            mock_from_filepath.return_value = [mol]
            result = runner.invoke(
                entry_point,
                ["-f", str(input_file)],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mol.write_xyz.assert_called_once_with(
            str(tmp_path / "mymol_single.xyz")
        )

    def test_multiple_structures_single_file_mode(self, tmp_path):
        from chemsmart.scripts.write_xyz import entry_point

        input_file = tmp_path / "traj.xyz"
        input_file.write_text("dummy")

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.write_xyz.Molecule.from_filepath"
        ) as mock_from_filepath:
            mols = [MagicMock(), MagicMock()]
            mock_from_filepath.return_value = mols
            result = runner.invoke(
                entry_point,
                ["-f", str(input_file), "-i", ":"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        expected_path = str(tmp_path / "traj_all.xyz")
        for mol in mols:
            mol.write_xyz.assert_called_once_with(expected_path, mode="a")

    def test_multiple_structures_multi_file_mode(self, tmp_path):
        from chemsmart.scripts.write_xyz import entry_point

        input_file = tmp_path / "traj.xyz"
        input_file.write_text("dummy")

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.write_xyz.Molecule.from_filepath"
        ) as mock_from_filepath:
            mols = [MagicMock(), MagicMock()]
            mock_from_filepath.return_value = mols
            result = runner.invoke(
                entry_point,
                ["-f", str(input_file), "-i", ":", "--no-single-files"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mols[0].write_xyz.assert_called_once_with(str(tmp_path / "traj_1.xyz"))
        mols[1].write_xyz.assert_called_once_with(str(tmp_path / "traj_2.xyz"))

    def test_missing_file_logs_error_and_returns(self, tmp_path):
        from chemsmart.scripts.write_xyz import entry_point

        runner = CliRunner()
        missing_file = str(tmp_path / "does_not_exist.xyz")
        with patch(
            "chemsmart.scripts.write_xyz.Molecule.from_filepath",
            side_effect=FileNotFoundError("nope"),
        ):
            result = runner.invoke(
                entry_point,
                ["-f", missing_file],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output


class TestHirshfeldScript:
    def test_gaussian_output_dispatch(self):
        from chemsmart.scripts.hirshfeld import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.hirshfeld.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.scripts.hirshfeld.Gaussian16Output") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.hirshfeld_charges = {"1C": 0.1}
            mock_output.hirshfeld_spin_densities = {"1C": 0.0}
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.log"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_cls.assert_called_once_with(filename="some.log")

    def test_orca_output_dispatch(self):
        from chemsmart.scripts.hirshfeld import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.hirshfeld.get_program_type_from_file",
                return_value="orca",
            ),
            patch("chemsmart.scripts.hirshfeld.ORCAOutput") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.hirshfeld_charges = {"1C": 0.1}
            mock_output.hirshfeld_spin_densities = {"1C": 0.0}
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.out"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_cls.assert_called_once_with(filename="some.out")

    def test_unknown_program_raises(self):
        from chemsmart.scripts.hirshfeld import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.hirshfeld.get_program_type_from_file",
            return_value="unknown",
        ):
            with pytest.raises(TypeError, match="unknown filetype"):
                runner.invoke(
                    entry_point,
                    ["-f", "some.txt"],
                    catch_exceptions=False,
                )


class TestMullikenScript:
    def test_unrestricted_spin_reports_densities(self):
        from chemsmart.scripts.mulliken import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.mulliken.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.scripts.mulliken.Gaussian16Output") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.spin = "unrestricted"
            mock_output.mulliken_atomic_charges = {"1C": 0.1}
            mock_output.mulliken_spin_densities = {"1C": 0.5}
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.log", "-n", "1"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output

    def test_restricted_spin_skips_densities(self):
        from chemsmart.scripts.mulliken import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.mulliken.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.scripts.mulliken.Gaussian16Output") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.spin = "restricted"
            mock_output.mulliken_atomic_charges = {"1C": 0.1}
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.log"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output

    def test_missing_spin_raises(self):
        from chemsmart.scripts.mulliken import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.mulliken.get_program_type_from_file",
                return_value="orca",
            ),
            patch("chemsmart.scripts.mulliken.ORCAOutput") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.spin = None
            mock_cls.return_value = mock_output
            with pytest.raises(ValueError, match="No spin information"):
                runner.invoke(
                    entry_point,
                    ["-f", "some.out"],
                    catch_exceptions=False,
                )
