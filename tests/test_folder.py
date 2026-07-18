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


class TestClickFolderOptionsDecorator:
    """``click_folder_options`` registers the expected CLI options."""

    def test_directory_option_present_in_help(self, invoke_folder_command):
        """``-d/--directory`` appears in the command help text."""
        result = invoke_folder_command(["--help"])
        assert result.exit_code == 0, result.output
        assert "--directory" in result.output

    def test_filetype_option_present_in_help(self, invoke_folder_command):
        """``-t/--filetype`` appears in the command help text."""
        result = invoke_folder_command(["--help"])
        assert result.exit_code == 0, result.output
        assert "--filetype" in result.output

    def test_program_option_present_in_help(self, invoke_folder_command):
        """``-p/--program`` appears in the command help text."""
        result = invoke_folder_command(["--help"])
        assert result.exit_code == 0, result.output
        assert "--program" in result.output

    def test_directory_option_received_by_callback(
        self, invoke_folder_command
    ):
        """Passing ``-d /some/path`` propagates the value to the callback."""
        result = invoke_folder_command(["-d", "/some/path"])
        assert result.exit_code == 0, result.output

    def test_filetype_option_received_by_callback(self, invoke_folder_command):
        """Passing ``-t log`` propagates the value to the callback."""
        result = invoke_folder_command(["-t", "log"])
        assert result.exit_code == 0, result.output

    def test_program_option_received_by_callback(self, invoke_folder_command):
        """Passing ``-p gaussian`` propagates the value to the callback."""
        result = invoke_folder_command(["-p", "gaussian"])
        assert result.exit_code == 0, result.output
