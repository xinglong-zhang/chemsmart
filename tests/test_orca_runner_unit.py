"""
Direct unit tests for the real :class:`ORCAJobRunner` (not
``FakeORCAJobRunner``, which is already exercised indirectly via the
``orca_jobrunner_no_scratch``/``orca_jobrunner_scratch`` fixtures used
throughout the suite).

Covers scratch vs. job-directory path setup, XYZ/extra-file reference
scanning and copying (used by NEB and COSMO-RS jobs), command
generation, subprocess creation, and postrun file copying.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.runner import ORCAJobRunner
from chemsmart.jobs.orca.settings import ORCAJobSettings


@pytest.fixture()
def orca_job(single_molecule_xyz_file, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = Molecule.from_filepath(single_molecule_xyz_file)
    settings = ORCAJobSettings.default()
    return ORCAJob(molecule=molecule, settings=settings, label="orcatest")


class TestORCAJobRunnerSetup:
    def test_defaults_scratch_to_true(self, pbs_server):
        assert ORCAJobRunner.SCRATCH is True
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        assert runner.scratch is False

    def test_scratch_none_resolves_same_as_explicit_true(self, pbs_server):
        # Exercises `if scratch is None: scratch = self.SCRATCH` (True):
        # the base JobRunner may still turn scratch off afterward if it
        # can't determine a scratch_dir, so just confirm scratch=None
        # produces the exact same resolved value as scratch=True,
        # rather than asserting a specific final True/False.
        runner_none = ORCAJobRunner(server=pbs_server, scratch=None)
        runner_true = ORCAJobRunner(server=pbs_server, scratch=True)
        assert runner_none.scratch == runner_true.scratch

    def test_assign_variables_job_directory_paths(self, pbs_server, orca_job):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(orca_job)

        assert runner.running_directory == orca_job.folder
        assert runner.job_inputfile == os.path.abspath(orca_job.inputfile)
        assert runner.job_gbwfile == os.path.abspath(orca_job.gbwfile)
        assert runner.job_errfile == os.path.abspath(orca_job.errfile)
        assert runner.job_outputfile == os.path.abspath(orca_job.outputfile)

    def test_assign_variables_scratch_paths(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)

        expected_dir = os.path.join(str(scratch_dir), orca_job.label)
        assert runner.running_directory == expected_dir
        assert os.path.isdir(expected_dir)
        assert runner.job_inputfile.endswith("orcatest.inp")
        assert runner.job_gbwfile.endswith("orcatest.gbw")
        assert runner.job_errfile.endswith("orcatest.err")
        assert runner.job_outputfile.endswith("orcatest.out")

    def test_prerun_job_directory_skips_xyz_copy(self, pbs_server, orca_job):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        with (
            patch.object(runner, "_copy_over_xyz_files") as mock_copy_xyz,
            patch.object(runner, "_copy_over_extra_files") as mock_copy_extra,
        ):
            runner._prerun(orca_job)

        mock_copy_xyz.assert_not_called()
        mock_copy_extra.assert_not_called()

    def test_prerun_scratch_copies_when_inputfile_exists(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        os.makedirs(os.path.dirname(orca_job.inputfile), exist_ok=True)
        with open(orca_job.inputfile, "w") as f:
            f.write("! opt\n")

        with (
            patch.object(runner, "_copy_over_xyz_files") as mock_copy_xyz,
            patch.object(runner, "_copy_over_extra_files") as mock_copy_extra,
        ):
            runner._prerun(orca_job)

        mock_copy_xyz.assert_called_once_with(orca_job)
        mock_copy_extra.assert_called_once_with(orca_job)


class TestORCAJobRunnerCopyXyzFiles:
    def test_restart_allxyz_marker_without_valid_filename_is_ignored(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)
        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            # Contains the marker text but no parseable .allxyz filename.
            f.write("Restart_ALLXYZFile\n")

        # Should not raise (no match found, nothing to copy).
        runner._copy_over_xyz_files(orca_job)

    def test_copies_referenced_xyz_file_to_scratch(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)

        referenced_xyz = os.path.join(orca_job.folder, "product.xyz")
        with open(referenced_xyz, "w") as f:
            f.write("1\ncomment\nH 0 0 0\n")

        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            f.write('NEB_End_XYZFile "product.xyz"\n')

        runner._copy_over_xyz_files(orca_job)

        assert os.path.exists(
            os.path.join(runner.running_directory, "product.xyz")
        )

    def test_raises_when_referenced_xyz_file_missing(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)
        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            f.write('NEB_End_XYZFile "missing.xyz"\n')

        with pytest.raises(FileNotFoundError):
            runner._copy_over_xyz_files(orca_job)

    def test_no_op_when_no_xyz_reference_in_input(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)
        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            f.write("! opt freq\n")

        # Should not raise.
        runner._copy_over_xyz_files(orca_job)

    def test_supports_allxyz_restart_reference(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)

        restart_file = os.path.join(orca_job.folder, "restart.allxyz")
        with open(restart_file, "w") as f:
            f.write("dummy allxyz content")

        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            f.write("Restart_ALLXYZFile restart.allxyz\n")

        runner._copy_over_xyz_files(orca_job)

        assert os.path.exists(
            os.path.join(runner.running_directory, "restart.allxyz")
        )

    def test_scratch_job_dir_already_exists_is_reused(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        expected_dir = os.path.join(str(scratch_dir), orca_job.label)
        os.makedirs(expected_dir)
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        assert runner.running_directory == expected_dir

    def test_absolute_xyz_path_used_directly(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)

        abs_xyz = tmp_path / "abs_product.xyz"
        abs_xyz.write_text("1\ncomment\nH 0 0 0\n")

        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            f.write(f'NEB_End_XYZFile "{abs_xyz}"\n')

        runner._copy_over_xyz_files(orca_job)
        assert os.path.exists(
            os.path.join(runner.running_directory, "abs_product.xyz")
        )

    def test_no_copy_when_not_using_scratch(
        self, pbs_server, orca_job, tmp_path
    ):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(orca_job)
        os.makedirs(os.path.dirname(orca_job.inputfile), exist_ok=True)

        referenced_xyz = os.path.join(orca_job.folder, "product.xyz")
        with open(referenced_xyz, "w") as f:
            f.write("1\ncomment\nH 0 0 0\n")
        with open(orca_job.inputfile, "w") as f:
            f.write('NEB_End_XYZFile "product.xyz"\n')

        # Should not raise; nothing gets copied since scratch is off.
        runner._copy_over_xyz_files(orca_job)


class TestORCAJobRunnerCopyExtraFiles:
    def test_no_op_without_scratch(self, pbs_server, orca_job):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(orca_job)
        # Should not raise even without scratch_dir.
        runner._copy_over_extra_files(orca_job)

    def test_no_op_when_input_file_does_not_exist(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)
        # Neither job_inputfile nor job.inputfile exist on disk.
        runner._copy_over_extra_files(orca_job)  # should not raise

    def test_copies_referenced_cosmorsxyz_file(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)

        cosmors_file = os.path.join(orca_job.folder, "water.cosmorsxyz")
        with open(cosmors_file, "w") as f:
            f.write("dummy cosmors content")

        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            f.write('  solventfilename "water"\n')

        runner._copy_over_extra_files(orca_job)

        assert os.path.exists(
            os.path.join(runner.running_directory, "water.cosmorsxyz")
        )

    def test_skips_when_referenced_file_not_found(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)
        os.makedirs(os.path.dirname(runner.job_inputfile), exist_ok=True)
        with open(runner.job_inputfile, "w") as f:
            f.write('  solventfilename "nonexistent"\n')

        # Should not raise, just log and skip.
        runner._copy_over_extra_files(orca_job)


class TestORCAJobRunnerCommand:
    def test_get_command_uses_executable_and_inputfile(
        self, pbs_server, orca_job
    ):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(orca_job)

        with patch.object(
            runner, "_get_executable", return_value="/path/to/orca"
        ):
            command = runner._get_command(orca_job)

        assert command == f"/path/to/orca {runner.job_inputfile}"

    def test_create_process_invokes_subprocess_popen(
        self, pbs_server, orca_job
    ):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)

        mock_executable = MagicMock()
        mock_executable.env = {}
        with (
            patch.object(type(runner), "executable", new=mock_executable),
            patch("chemsmart.jobs.orca.runner.subprocess.Popen") as mock_popen,
        ):
            mock_popen.return_value = MagicMock()
            result = runner._create_process(orca_job, "orca input.inp", env={})

        mock_popen.assert_called_once()
        call_args = mock_popen.call_args
        assert call_args.args[0] == ["orca", "input.inp"]
        assert call_args.kwargs["cwd"] == runner.running_directory
        assert result is mock_popen.return_value


class TestORCAJobRunnerPostrun:
    def test_postrun_no_op_without_scratch(self, pbs_server, orca_job):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(orca_job)
        # Should not raise.
        runner._postrun(orca_job)

    def test_postrun_copies_label_prefixed_files_skips_tmp(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)
        os.makedirs(orca_job.folder, exist_ok=True)

        keep_file = os.path.join(runner.running_directory, "orcatest.out")
        with open(keep_file, "w") as f:
            f.write("output")
        tmp_file = os.path.join(runner.running_directory, "orcatest.tmp")
        with open(tmp_file, "w") as f:
            f.write("temp")

        runner._postrun(orca_job)

        assert os.path.exists(os.path.join(orca_job.folder, "orcatest.out"))
        assert not os.path.exists(
            os.path.join(orca_job.folder, "orcatest.tmp")
        )

    def test_postrun_logs_error_when_copy_fails(
        self, pbs_server, orca_job, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)
        os.makedirs(orca_job.folder, exist_ok=True)

        keep_file = os.path.join(runner.running_directory, "orcatest.out")
        with open(keep_file, "w") as f:
            f.write("output")

        with patch(
            "chemsmart.jobs.orca.runner.copy",
            side_effect=OSError("disk full"),
        ):
            runner._postrun(orca_job)  # should not raise


class TestORCAJobRunnerExecutable:
    def test_executable_raises_file_not_found_with_context(self, pbs_server):
        from chemsmart.settings.executable import ORCAExecutable

        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        type(runner).executable.fget.cache_clear()
        with patch.object(
            ORCAExecutable,
            "from_servername",
            side_effect=FileNotFoundError("no such server config"),
        ):
            with pytest.raises(FileNotFoundError):
                _ = runner.executable
        type(runner).executable.fget.cache_clear()

    def test_get_executable_delegates_to_executable(
        self, pbs_server, orca_job
    ):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        mock_executable = MagicMock()
        mock_executable.get_executable.return_value = "/opt/orca/orca"
        with patch.object(type(runner), "executable", new=mock_executable):
            exe = runner._get_executable()
        assert exe == "/opt/orca/orca"

    def test_assign_variables_sets_local_run_from_executable(
        self, pbs_server, orca_job
    ):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        mock_executable = MagicMock()
        mock_executable.local_run = True
        with patch.object(type(runner), "executable", new=mock_executable):
            runner._assign_variables(orca_job)
        assert orca_job.local is True

    def test_assign_variables_skips_local_run_when_none(
        self, pbs_server, orca_job
    ):
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        mock_executable = MagicMock()
        mock_executable.local_run = None
        with patch.object(type(runner), "executable", new=mock_executable):
            runner._assign_variables(orca_job)
        # job.local should remain at whatever default it started with,
        # i.e. untouched by _assign_variables.
        assert not hasattr(orca_job, "local") or orca_job.local is not True


class TestWriteInput:
    def test_write_input_writes_and_copies_extra_files(
        self, pbs_server, orca_job, tmp_path
    ):
        orca_job.settings.functional = "b3lyp"
        orca_job.settings.basis = "sto-3g"
        orca_job.settings.charge = 0
        orca_job.settings.multiplicity = 1
        scratch_dir = tmp_path / "scratch"
        runner = ORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        orca_job.jobrunner = runner
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)

        with patch.object(runner, "_copy_over_extra_files") as mock_copy_extra:
            runner._write_input(orca_job)

        assert os.path.exists(runner.job_inputfile)
        mock_copy_extra.assert_called_once_with(orca_job)

    def test_write_input_skips_extra_file_copy_without_scratch(
        self, pbs_server, orca_job
    ):
        orca_job.settings.functional = "b3lyp"
        orca_job.settings.basis = "sto-3g"
        orca_job.settings.charge = 0
        orca_job.settings.multiplicity = 1
        runner = ORCAJobRunner(server=pbs_server, scratch=False)
        orca_job.jobrunner = runner
        runner._assign_variables(orca_job)
        os.makedirs(runner.running_directory, exist_ok=True)

        with patch.object(runner, "_copy_over_extra_files") as mock_copy_extra:
            runner._write_input(orca_job)

        mock_copy_extra.assert_not_called()


class TestFakeORCAJobRunnerRun:
    def test_run_executes_full_fake_pipeline(
        self, pbs_server, orca_job, tmp_path
    ):
        from chemsmart.jobs.orca.runner import FakeORCAJobRunner

        orca_job.settings.functional = "b3lyp"
        orca_job.settings.basis = "sto-3g"
        orca_job.settings.charge = 0
        orca_job.settings.multiplicity = 1
        scratch_dir = tmp_path / "scratch"
        runner = FakeORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        orca_job.jobrunner = runner

        returncode = runner.run(orca_job)

        assert returncode is None or returncode == 0
        # orca_job.label has already had "_fake" appended in place by
        # _append_suffix_to_job_label during the run.
        assert orca_job.label.endswith("_fake")
        assert os.path.exists(
            os.path.join(orca_job.folder, orca_job.label + ".out")
        )

    def test_scratch_job_dir_already_exists_is_reused(
        self, pbs_server, orca_job, tmp_path
    ):
        from chemsmart.jobs.orca.runner import FakeORCAJobRunner

        scratch_dir = tmp_path / "scratch"
        expected_dir = os.path.join(str(scratch_dir), orca_job.label)
        os.makedirs(expected_dir)
        runner = FakeORCAJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(orca_job)
        assert runner.running_directory == expected_dir
