import io
import os
from shutil import copy
from types import SimpleNamespace

import pytest

from chemsmart.jobs.nciplot import NCIPLOTJob
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings
from chemsmart.jobs.nciplot.writer import NCIPLOTInputWriter


class TestNCIPLOTInputWriter:
    def test_write_nci_for_one_file(
        self,
        tmpdir,
        monkeypatch,
        single_molecule_xyz_file,
        nciplot_jobrunner_no_scratch,
    ):
        job_settings = NCIPLOTJobSettings(
            rthres=10.0,
            ligand_file_number=1,
            ligand_radius=1.0,
            radius_positions="(1.0, 1.1, 1.2)",
            radius_r=1.5,
            intercut1=0.5,
            intercut2=0.6,
            increments="0.1,0.1,0.1",
            fragments={1: [1, 2, 3], 2: [4, 5, 6]},
            cutoff_density_dat=0.01,
            cutoff_rdg_dat=0.02,
            cutoff_density_cube=0.03,
            cutoff_rdg_cube=0.04,
            dgrid=True,
            integrate=True,
            ranges=[[-0.1, -0.02], [-0.02, 0.02], [0.02, 0.1]],
            grid_quality="coarse",
        )

        # copy file to tmpdir
        tmpdir_xyz_file = os.path.join(tmpdir, "single_molecule.xyz")
        copy(single_molecule_xyz_file, tmpdir_xyz_file)

        # change to tmpdir and run test
        monkeypatch.chdir(tmpdir)

        # create nciplot job
        job = NCIPLOTJob(
            filenames=("single_molecule.xyz",),
            settings=job_settings,
            jobrunner=nciplot_jobrunner_no_scratch,
        )
        # set up correct variables by calling prerun()
        nciplot_jobrunner_no_scratch._prerun(job)
        nciplot_writer = NCIPLOTInputWriter(job=job)

        # write input file
        nciplot_writer.write(target_directory=tmpdir)

        # if no job label, then base filename is used
        nci_file = os.path.join(tmpdir, "single_molecule.nci")
        assert os.path.exists(nci_file)
        lines = open(nci_file, "r").readlines()
        assert lines[0] == "1\n"
        assert lines[1] == "single_molecule.xyz\n"
        assert lines[2] == "RTHRES 10.0\n"
        assert lines[3] == "LIGAND 1 1.0\n"
        assert lines[4] == "RADIUS 1.0  1.1  1.2 1.5\n"
        assert lines[5] == "INTERMOLECULAR\n"
        assert lines[6] == "INTERCUT 0.5 0.6\n"
        assert lines[7] == "INCREMENTS 0.1 0.1 0.1\n"
        assert lines[8] == "FRAGMENTS\n"
        assert lines[9] == "  1 1 2 3\n"
        assert lines[10] == "  2 4 5 6\n"
        assert lines[11] == "END\n"
        assert lines[12] == "CUTOFFS 0.01 0.02\n"
        assert lines[13] == "CUTPLOT 0.03 0.04\n"
        assert lines[14] == "DGRID\n"
        assert lines[15] == "INTEGRATE\n"
        assert lines[16] == "RANGE 3\n"
        assert lines[17] == "-0.1 -0.02\n"
        assert lines[18] == "-0.02 0.02\n"
        assert lines[19] == "0.02 0.1\n"
        assert lines[20] == "COARSE\n"

    def test_write_nci_for_two_files(
        self,
        tmpdir,
        monkeypatch,
        single_molecule_xyz_file,
        nciplot_jobrunner_no_scratch,
    ):
        job_settings = NCIPLOTJobSettings(
            rthres=10.0,
            ligand_file_number=1,
            ligand_radius=1.0,
            radius_positions="(1.0, 1.1, 1.2)",
            radius_r=1.5,
            intercut1=0.5,
            intercut2=0.6,
            increments="0.1,0.1,0.1",
            fragments={1: [1, 2, 3], 2: [4, 5, 6]},
            cutoff_density_dat=0.01,
            cutoff_rdg_dat=0.02,
            cutoff_density_cube=0.03,
            cutoff_rdg_cube=0.04,
            dgrid=True,
            integrate=True,
            ranges=[[-0.1, -0.02], [-0.02, 0.02], [0.02, 0.1]],
            grid_quality="ultrafine",
        )

        # copy file to tmpdir
        tmpdir_xyz_file = os.path.join(tmpdir, "single_molecule.xyz")
        copy(single_molecule_xyz_file, tmpdir_xyz_file)

        tmpdir_xyz_file2 = os.path.join(tmpdir, "single_molecule2.xyz")
        copy(single_molecule_xyz_file, tmpdir_xyz_file2)

        # change to tmpdir and run test
        monkeypatch.chdir(tmpdir)

        # create nciplot job
        job = NCIPLOTJob(
            filenames=("single_molecule.xyz", "single_molecule2.xyz"),
            settings=job_settings,
            label="nci_two_files",
            jobrunner=nciplot_jobrunner_no_scratch,
        )

        # set up correct variables by calling prerun()
        nciplot_jobrunner_no_scratch._prerun(job)

        nciplot_writer = NCIPLOTInputWriter(job=job)

        # write input file
        nciplot_writer.write(target_directory=tmpdir)

        # if job label is given, then filename takes job label
        # NCIJobSettings.label is not used
        nci_file = os.path.join(tmpdir, "nci_two_files.nci")
        assert os.path.exists(nci_file)
        lines = open(nci_file, "r").readlines()
        assert lines[0] == "2\n"
        assert lines[1] == "single_molecule.xyz\n"
        assert lines[2] == "single_molecule2.xyz\n"
        assert lines[3] == "RTHRES 10.0\n"
        assert lines[4] == "LIGAND 1 1.0\n"
        assert lines[5] == "RADIUS 1.0  1.1  1.2 1.5\n"
        assert lines[6] == "INTERMOLECULAR\n"
        assert lines[7] == "INTERCUT 0.5 0.6\n"
        assert lines[8] == "INCREMENTS 0.1 0.1 0.1\n"
        assert lines[9] == "FRAGMENTS\n"
        assert lines[10] == "  1 1 2 3\n"
        assert lines[11] == "  2 4 5 6\n"
        assert lines[12] == "END\n"
        assert lines[13] == "CUTOFFS 0.01 0.02\n"
        assert lines[14] == "CUTPLOT 0.03 0.04\n"
        assert lines[15] == "DGRID\n"
        assert lines[16] == "INTEGRATE\n"
        assert lines[17] == "RANGE 3\n"
        assert lines[18] == "-0.1 -0.02\n"
        assert lines[19] == "-0.02 0.02\n"
        assert lines[20] == "0.02 0.1\n"
        assert lines[21] == "ULTRAFINE\n"

        # job run will result in the job being run and
        # the output file copied back to run folder
        # job.run()
        # assert job.is_complete()

    def test_write_nci_promolecular_with_scratch(
        self,
        tmpdir,
        monkeypatch,
        single_molecule_xyz_file,
        nciplot_jobrunner_scratch,
    ):
        """Test that promolecular density files are found in scratch directory.

        This test reproduces the bug where the writer looks for files in the
        wrong directory when using scratch and promolecular density.
        """
        job_settings = NCIPLOTJobSettings()

        # copy file to tmpdir to simulate original location
        tmpdir_xyz_file = os.path.join(tmpdir, "test_molecule.xyz")
        copy(single_molecule_xyz_file, tmpdir_xyz_file)

        # change to tmpdir
        monkeypatch.chdir(tmpdir)

        # create nciplot job with promolecular suffix in label
        # This simulates what happens in the CLI when processing .xyz files
        job = NCIPLOTJob(
            filenames=("test_molecule.xyz",),
            settings=job_settings,
            label="test_molecule_promolecular",
            jobrunner=nciplot_jobrunner_scratch,
        )

        # Simulate the runner's prerun which
        # creates scratch dir and copies files
        nciplot_jobrunner_scratch._prerun(job)

        # Now write the input file - this should not raise FileNotFoundError
        nciplot_writer = NCIPLOTInputWriter(job=job)
        nciplot_writer.write(
            target_directory=nciplot_jobrunner_scratch.running_directory
        )

        # Verify the input file was created correctly
        nci_file = os.path.join(
            nciplot_jobrunner_scratch.running_directory,
            "test_molecule_promolecular.nci",
        )
        assert os.path.exists(nci_file)

        # Check that the filename is written correctly
        with open(nci_file, "r") as f:
            lines = f.readlines()
        assert lines[0] == "1\n"
        assert lines[1] == "test_molecule.xyz\n"


class TestNCIPLOTPathHandling:
    """Test NCIPLOT handling of full file paths."""

    def test_write_nci_with_full_path_in_scratch(
        self,
        tmpdir,
        monkeypatch,
        single_molecule_xyz_file,
        nciplot_jobrunner_scratch,
    ):
        """Test that NCIPLOT correctly handles full file paths.

        This test verifies the fix for a bug where providing a full path to a
        .xyz file caused a FileNotFoundError because the writer tried to use
        the full path instead of just the basename when validating file
        existence in scratch.
        """
        job_settings = NCIPLOTJobSettings()

        # Create a source directory to simulate
        # files being in a different location
        source_dir = os.path.join(tmpdir, "source")
        os.makedirs(source_dir)

        # Copy file to source directory
        source_xyz_file = os.path.join(source_dir, "test_molecule.xyz")
        copy(single_molecule_xyz_file, source_xyz_file)

        # Change to a different working directory
        work_dir = os.path.join(tmpdir, "work")
        os.makedirs(work_dir)
        monkeypatch.chdir(work_dir)

        # Create job with FULL PATH to the file (simulating real-world usage)
        # and with promolecular suffix in label (as the CLI would add)
        job = NCIPLOTJob(
            filenames=(source_xyz_file,),  # Full path!
            settings=job_settings,
            label="test_molecule_promolecular",
            jobrunner=nciplot_jobrunner_scratch,
        )

        # Set up scratch directory (with promolecular suffix from label)
        nciplot_jobrunner_scratch._prerun(job)

        # This should NOT raise FileNotFoundError
        # The bug was that the writer would try to find the file at:
        # <scratch_dir>/<full_path> instead of <scratch_dir>/<basename>
        nciplot_writer = NCIPLOTInputWriter(job=job)
        nciplot_writer.write(
            target_directory=nciplot_jobrunner_scratch.running_directory
        )

        # Verify the input file was created correctly
        nci_file = os.path.join(
            nciplot_jobrunner_scratch.running_directory,
            "test_molecule_promolecular.nci",
        )
        assert os.path.exists(nci_file)

        # Check that the filename written is
        # just the basename, not the full path
        with open(nci_file, "r") as f:
            lines = f.readlines()
        assert lines[0] == "1\n"
        assert lines[1] == "test_molecule.xyz\n"  # Should be basename only

    def test_write_nci_converted_file_with_full_path(
        self,
        tmpdir,
        monkeypatch,
        nciplot_jobrunner_scratch,
    ):
        """Test NCIPLOT correctly handles non-.xyz files with full paths.

        This test verifies that when a non-supported file (e.g., .log) is
        converted to .xyz, the writer correctly uses the basename with
        _promolecular suffix for both validation and writing to the .nci file.
        """
        job_settings = NCIPLOTJobSettings()

        # Create a source directory with a fake .log file
        source_dir = os.path.join(tmpdir, "source")
        os.makedirs(source_dir)
        source_log_file = os.path.join(source_dir, "test_calculation.log")

        # Create a minimal fake log file (won't
        # actually be converted in this test)
        with open(source_log_file, "w") as f:
            f.write("Fake log file content\n")

        # Change to work directory
        work_dir = os.path.join(tmpdir, "work")
        os.makedirs(work_dir)
        monkeypatch.chdir(work_dir)

        # Create the converted .xyz file that the runner would have created
        # (simulating what _write_xyz_from_input_files does)
        converted_xyz = os.path.join(
            source_dir, "test_calculation_promolecular.xyz"
        )
        with open(converted_xyz, "w") as f:
            f.write("1\nConverted XYZ\nC 0.0 0.0 0.0\n")

        # Create job with full path to .log file
        job = NCIPLOTJob(
            filenames=(source_log_file,),  # Full path to .log
            settings=job_settings,
            label="test_calculation_promolecular",
            jobrunner=nciplot_jobrunner_scratch,
        )

        # Manually set up scratch and copy the converted file
        # (normally done by runner._write_xyz_from_input_files)
        nciplot_jobrunner_scratch._assign_variables(job)
        copy(
            converted_xyz,
            os.path.join(
                nciplot_jobrunner_scratch.running_directory,
                "test_calculation_promolecular.xyz",
            ),
        )

        # Write input file - should handle the basename correctly
        nciplot_writer = NCIPLOTInputWriter(job=job)
        nciplot_writer.write(
            target_directory=nciplot_jobrunner_scratch.running_directory
        )

        # Verify correct filename was written
        nci_file = os.path.join(
            nciplot_jobrunner_scratch.running_directory,
            "test_calculation_promolecular.nci",
        )
        assert os.path.exists(nci_file)

        with open(nci_file, "r") as f:
            lines = f.readlines()
        assert lines[0] == "1\n"
        # Should be basename with _promolecular suffix
        assert lines[1] == "test_calculation_promolecular.xyz\n"


def _make_writer(tmp_path, **setting_overrides):
    """Direct-construction helper bypassing NCIPLOTJob/JobRunner
    entirely, since most of the writer's validation branches are pure
    functions of self.job.filenames/self.settings and don't need a
    real job or running directory."""
    writer = NCIPLOTInputWriter.__new__(NCIPLOTInputWriter)
    writer.job = SimpleNamespace(
        filenames=None, label="test", folder=str(tmp_path)
    )
    writer.settings = NCIPLOTJobSettings()
    writer.jobrunner = SimpleNamespace(running_directory=str(tmp_path))
    for key, value in setting_overrides.items():
        setattr(writer.settings, key, value)
    return writer


class TestNCIPLOTWriterValidationBranches:
    """Direct unit coverage for the writer's error/edge branches, none
    of which the full-job success-path tests above reach."""

    def test_write_creates_missing_target_directory(self, tmp_path):
        writer = _make_writer(tmp_path)
        new_dir = tmp_path / "newsub"
        writer._write(target_directory=str(new_dir))
        assert new_dir.is_dir()
        assert (new_dir / "test.nci").exists()

    def test_write_uses_job_folder_when_no_target_directory(self, tmp_path):
        writer = _make_writer(tmp_path)
        writer._write(target_directory=None)
        assert (tmp_path / "test.nci").exists()

    def test_write_filenames_empty_list_raises(self, tmp_path):
        writer = _make_writer(tmp_path)
        writer.job.filenames = []
        with pytest.raises(ValueError, match="No filenames provided"):
            writer._write_filenames(io.StringIO())

    def test_write_filenames_missing_file_raises(self, tmp_path):
        writer = _make_writer(tmp_path)
        writer.job.filenames = ["missing.xyz"]
        with pytest.raises(FileNotFoundError, match="does not exist"):
            writer._write_filenames(io.StringIO())

    def test_write_ligand_only_one_value_raises(self, tmp_path):
        writer = _make_writer(
            tmp_path, ligand_file_number=1, ligand_radius=None
        )
        with pytest.raises(ValueError, match="must be provided or"):
            writer._write_ligand(io.StringIO())

    def test_write_radius_wrong_coordinate_count_raises(self, tmp_path):
        writer = _make_writer(
            tmp_path, radius_positions="1.0,2.0", radius_r=1.5
        )
        with pytest.raises(ValueError, match="exactly 3 coordinates"):
            writer._write_radius(io.StringIO())

    def test_write_radius_invalid_coordinate_value_raises(self, tmp_path):
        writer = _make_writer(
            tmp_path, radius_positions="1.0,x,2.0", radius_r=1.5
        )
        with pytest.raises(ValueError, match="Invalid coordinate value"):
            writer._write_radius(io.StringIO())

    def test_write_radius_only_one_value_raises(self, tmp_path):
        writer = _make_writer(
            tmp_path, radius_positions="1.0,2.0,3.0", radius_r=None
        )
        with pytest.raises(ValueError, match="must be provided or"):
            writer._write_radius(io.StringIO())

    def test_write_intermolecular_negative_value_raises(self, tmp_path):
        writer = _make_writer(tmp_path, intercut1=-0.5)
        with pytest.raises(ValueError, match="must be positive"):
            writer._write_intermolecular(io.StringIO())

    def test_write_increments_invalid_value_raises(self, tmp_path):
        writer = _make_writer(tmp_path, increments="0.1,x,0.1")
        with pytest.raises(ValueError, match="Invalid increment value"):
            writer._write_increments(io.StringIO())

    def test_write_fragments_value_not_list_or_tuple_raises(self, tmp_path):
        writer = _make_writer(tmp_path, fragments={1: "notalist"})
        with pytest.raises(ValueError, match="must be a list or tuple"):
            writer._write_fragments(io.StringIO())

    def test_write_fragments_not_a_dict_raises(self, tmp_path):
        writer = _make_writer(tmp_path, fragments=["a", "b"])
        with pytest.raises(ValueError, match="must be a dictionary"):
            writer._write_fragments(io.StringIO())

    def test_write_cutoffs_negative_value_raises(self, tmp_path):
        writer = _make_writer(tmp_path, cutoff_density_dat=-0.5)
        with pytest.raises(ValueError, match="must be positive"):
            writer._write_cutoffs(io.StringIO())

    def test_write_cutplot_filenames_not_list_or_tuple_raises(self, tmp_path):
        writer = _make_writer(tmp_path)
        writer.job.filenames = "notalist"
        with pytest.raises(TypeError, match="list or tuple"):
            writer._write_cutplot(io.StringIO())

    def test_write_cutplot_empty_filenames_raises(self, tmp_path):
        writer = _make_writer(tmp_path)
        writer.job.filenames = []
        with pytest.raises(ValueError, match="No filenames provided"):
            writer._write_cutplot(io.StringIO())

    def test_write_cutplot_wfn_uses_scf_defaults(self, tmp_path):
        writer = _make_writer(
            tmp_path, cutoff_density_cube=None, cutoff_rdg_cube=0.6
        )
        writer.job.filenames = ["mol.wfn"]
        buf = io.StringIO()
        writer._write_cutplot(buf)
        # SCF default for r1 is 0.05 (vs. promolecular's 0.07)
        assert "CUTPLOT 0.05" in buf.getvalue()

    def test_write_cutplot_negative_value_raises(self, tmp_path):
        writer = _make_writer(
            tmp_path, cutoff_density_cube=-0.1, cutoff_rdg_cube=0.5
        )
        writer.job.filenames = ["mol.xyz"]
        with pytest.raises(ValueError, match="must be positive"):
            writer._write_cutplot(io.StringIO())

    def test_write_ranges_invalid_value_raises(self, tmp_path):
        writer = _make_writer(tmp_path, ranges=[["a", 1.0]])
        with pytest.raises(ValueError, match="Invalid range values"):
            writer._write_ranges(io.StringIO())

    def test_write_ranges_wrong_shape_raises(self, tmp_path):
        writer = _make_writer(tmp_path, ranges=[[1.0, 2.0, 3.0]])
        with pytest.raises(ValueError, match="list or tuple of two"):
            writer._write_ranges(io.StringIO())
