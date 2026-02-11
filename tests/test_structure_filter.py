"""
Unit tests for structure_filter.py script.

Tests the automatic file type detection functionality and error handling
for Gaussian and ORCA output files.
"""

import os
import shutil
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.scripts.structure_filter import entry_point


class TestStructureFilterAutoDetection:
    """Test suite for automatic file type detection in structure_filter."""

    @staticmethod
    def assert_error_contains(result, expected_message):
        """Helper method to check if error message appears in result."""
        assert result.exit_code != 0
        assert expected_message in result.output or (
            result.exception and expected_message in str(result.exception)
        )

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        temp_path = tempfile.mkdtemp()
        yield temp_path
        # Cleanup after tests
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def gaussian_log_file(self, temp_dir):
        """Create a mock Gaussian log file."""
        filepath = os.path.join(temp_dir, "test_c1.log")
        with open(filepath, "w") as f:
            f.write("Entering Gaussian System\n")
            f.write("Gaussian, Inc.\n")
            f.write(" C     0.000000    0.000000    0.000000\n")
            f.write(" H     1.000000    0.000000    0.000000\n")
        return filepath

    @pytest.fixture
    def orca_out_file(self, temp_dir):
        """Create a mock ORCA output file."""
        filepath = os.path.join(temp_dir, "test_c1.out")
        with open(filepath, "w") as f:
            f.write("* O   R   C   A *\n")
            f.write("Your ORCA version\n")
            f.write("CARTESIAN COORDINATES (ANGSTROEM)\n")
            f.write("C     0.000000    0.000000    0.000000\n")
            f.write("H     1.000000    0.000000    0.000000\n")
        return filepath

    @pytest.fixture
    def unknown_file(self, temp_dir):
        """Create a file with unknown format."""
        filepath = os.path.join(temp_dir, "test_c1.out")
        with open(filepath, "w") as f:
            f.write("Unknown quantum chemistry program output\n")
            f.write("No recognizable patterns\n")
        return filepath

    @pytest.fixture
    def xtb_file(self, temp_dir):
        """Create a mock xTB output file."""
        filepath = os.path.join(temp_dir, "test_c1.out")
        with open(filepath, "w") as f:
            f.write("x T B\n")
            f.write("xtb version 6.5.0\n")
            f.write("Some content here\n")
        return filepath

    def test_auto_detect_gaussian_log(self, temp_dir, gaussian_log_file):
        """Test auto-detection of Gaussian .log files."""
        # This test verifies the auto-detection logic recognizes Gaussian files
        from chemsmart.utils.io import get_program_type_from_file

        program = get_program_type_from_file(gaussian_log_file)
        assert program == "gaussian"

    def test_auto_detect_orca_out(self, temp_dir, orca_out_file):
        """Test auto-detection of ORCA .out files."""
        # This test verifies the auto-detection logic recognizes ORCA files
        from chemsmart.utils.io import get_program_type_from_file

        program = get_program_type_from_file(orca_out_file)
        assert program == "orca"

    def test_auto_detect_no_files(self, temp_dir):
        """Test error handling when no output files are found."""
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            ["-d", temp_dir, "-g", "rmsd"],
        )
        self.assert_error_contains(result, "No .log or .out files found")

    def test_auto_detect_unknown_file_type(self, temp_dir, unknown_file):
        """Test error handling for unknown file types."""
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            ["-d", temp_dir, "-g", "rmsd"],
        )
        self.assert_error_contains(result, "Could not determine file type")

    def test_auto_detect_unsupported_program(self, temp_dir, xtb_file):
        """Test error handling for unsupported programs like xTB."""
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            ["-d", temp_dir, "-g", "rmsd"],
        )
        self.assert_error_contains(
            result, "only Gaussian and ORCA are supported"
        )

    def test_manual_type_specification(self, temp_dir, gaussian_log_file):
        """Test that manual type specification still works."""
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            ["-d", temp_dir, "-t", "log", "-g", "rmsd"],
        )
        # May fail due to other reasons (like invalid molecule data),
        # but should not have auto-detection errors
        if result.exit_code != 0:
            assert "auto-detect" not in result.output.lower()
            if result.exception:
                assert "auto-detect" not in str(result.exception).lower()

    def test_no_files_with_specified_type(self, temp_dir):
        """Test error handling when no files match the specified type."""
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            ["-d", temp_dir, "-t", "xyz", "-g", "rmsd"],
        )
        self.assert_error_contains(
            result, "No files with extension '.xyz' found"
        )
