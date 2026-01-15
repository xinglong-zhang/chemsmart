"""
Test for NCIPLOT file path handling bug fix.

This test verifies that NCIPLOT correctly handles filenames with full paths
when using scratch directories and promolecular density.
"""

import os
import tempfile
from shutil import copy

import pytest

from chemsmart.jobs.nciplot import NCIPLOTJob
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings
from chemsmart.jobs.nciplot.writer import NCIPLOTInputWriter


class TestNCIPLOTPathHandling:
    """Test NCIPLOT handling of full file paths."""

    def test_write_nci_with_full_path_in_scratch(
        self,
        tmpdir,
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

        # Create a source directory to simulate files being in a different location
        source_dir = os.path.join(tmpdir, "source")
        os.makedirs(source_dir)

        # Copy file to source directory
        source_xyz_file = os.path.join(source_dir, "test_molecule.xyz")
        copy(single_molecule_xyz_file, source_xyz_file)

        # Change to a different working directory
        work_dir = os.path.join(tmpdir, "work")
        os.makedirs(work_dir)
        os.chdir(work_dir)

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

        # Check that the filename written is just the basename, not the full path
        with open(nci_file, "r") as f:
            lines = f.readlines()
        assert lines[0] == "1\n"
        assert lines[1] == "test_molecule.xyz\n"  # Should be basename only

    def test_write_nci_converted_file_with_full_path(
        self,
        tmpdir,
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

        # Create a minimal fake log file (won't actually be converted in this test)
        with open(source_log_file, "w") as f:
            f.write("Fake log file content\n")

        # Change to work directory
        work_dir = os.path.join(tmpdir, "work")
        os.makedirs(work_dir)
        os.chdir(work_dir)

        # Create the converted .xyz file that the runner would have created
        # (simulating what _write_xyz_from_input_files does)
        converted_xyz = os.path.join(source_dir, "test_calculation_promolecular.xyz")
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
