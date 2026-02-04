"""Tests for BaseFolder and folder operations."""

import os
import tempfile

from chemsmart.io.folder import BaseFolder


class TestFolderOutputFileDiscovery:
    """Tests for the BaseFolder output file discovery methods."""

    def test_find_gaussian_outputs(self, gaussian_outputs_test_directory):
        """Test finding Gaussian output files."""
        gaussian_output_files = BaseFolder(
            folder=gaussian_outputs_test_directory
        ).get_all_output_files_in_current_folder_by_program(program="gaussian")
        assert len(gaussian_output_files) > 0
        for file in gaussian_output_files:
            assert file.endswith((".log", ".out"))

    def test_find_orca_outputs(self, orca_outputs_directory):
        """Test finding ORCA output files."""
        orca_output_files = BaseFolder(
            folder=orca_outputs_directory
        ).get_all_output_files_in_current_folder_by_program(program="orca")
        assert len(orca_output_files) > 0
        for file in orca_output_files:
            assert file.endswith(".out")

    def test_find_xtb_outputs_in_parent_folder(self, xtb_outputs_directory):
        """Test that non-recursive search finds no xTB files in parent directory."""
        xtb_output_files = BaseFolder(
            folder=xtb_outputs_directory
        ).get_all_output_files_in_current_folder_by_program(program="xtb")
        assert not xtb_output_files  # No files in the top-level folder

    def test_find_xtb_outputs_in_subfolders(self, xtb_outputs_directory):
        """Test finding xTB output files with recursive search."""
        xtb_output_files = BaseFolder(
            folder=xtb_outputs_directory
        ).get_all_output_files_in_current_folder_and_subfolders_by_program(
            program="xtb"
        )
        # Recursive should find files in subdirectories
        assert len(xtb_output_files) > 0
        for file in xtb_output_files:
            assert file.endswith(".out")


class TestFolderIsProgramCalculationDirectory:
    """Tests for the is_program_calculation_directory method."""

    def test_xtb_output_directory(self, xtb_outputs_directory):
        """Test that parent directory is not detected as xTB calculation directory."""
        assert not BaseFolder(
            folder=xtb_outputs_directory
        ).is_program_calculation_directory("xtb")

    def test_xtb_co2_calculation_directory(self, xtb_co2_outfolder):
        """Test detection of xTB calculation directory for CO2."""
        assert BaseFolder(
            folder=xtb_co2_outfolder
        ).is_program_calculation_directory("xtb")

    def test_not_xtb_directory(self, gaussian_outputs_test_directory):
        """Test that Gaussian directory is not detected as xTB."""
        assert not BaseFolder(
            folder=gaussian_outputs_test_directory
        ).is_program_calculation_directory("xtb")


class TestFolderGetProgramType:
    """Tests for the get_program_type_from_folder method."""

    def test_xtb_output_detection(self, xtb_water_outfolder):
        """Test detection of xTB calculation folder."""
        result = BaseFolder(
            folder=xtb_water_outfolder
        ).get_program_type_from_folder()
        assert result == "xtb"

    def test_unknown_folder(self):
        """Test unknown folder type."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create some random files
            with open(os.path.join(temp_dir, "random.txt"), "w") as f:
                f.write("Random content\n")
            result = BaseFolder(folder=temp_dir).get_program_type_from_folder()
            assert result == "unknown"
