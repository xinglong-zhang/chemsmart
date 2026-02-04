"""Tests for BaseFolder and folder operations."""

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
