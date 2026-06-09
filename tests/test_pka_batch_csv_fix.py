"""
Test suite for pKa CSV batch processing with ORCA and Gaussian backends.

Tests the fix for handling CSV input files in pKa workflows. CSV files
require explicit subcommand specification (batch or submit) and should not
attempt to call resolve_proton_index() without a proton_index parameter.
"""

import os
import tempfile

import pytest
from click.testing import CliRunner


@pytest.fixture
def temp_xyz_file():
    """Create a temporary XYZ molecule file."""
    content = """3

H 0 0 0
H 0 0 0.74
O 0 0 0.37
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".xyz", delete=False
    ) as f:
        f.write(content)
        temp_path = f.name
    yield temp_path
    if os.path.exists(temp_path):
        os.unlink(temp_path)


@pytest.fixture
def temp_csv_file(temp_xyz_file):
    """Create a temporary CSV table file pointing to the XYZ file."""
    content = f"""filepath,proton_index,charge,multiplicity
{temp_xyz_file},1,0,1
{temp_xyz_file},2,0,1
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", delete=False
    ) as f:
        f.write(content)
        temp_path = f.name
    yield temp_path
    if os.path.exists(temp_path):
        os.unlink(temp_path)


class TestGaussianPkaCSVHandling:
    """Test Gaussian pKa CSV batch processing."""

    def test_csv_without_subcommand_raises_error(self, temp_csv_file):
        """CSV files without explicit subcommand should raise helpful error."""
        from chemsmart.cli.gaussian.gaussian import gaussian

        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-p",
                "test_project",
                "-f",
                temp_csv_file,
                "pka",
                "-s",
                "direct",
                # NOTE: No explicit subcommand (no 'batch' or 'submit')
            ],
        )

        assert result.exit_code != 0, "Should fail without explicit subcommand"
        assert (
            "CSV input files require an explicit subcommand" in result.output
        ), f"Expected helpful error message. Got:\n{result.output}"

    def test_csv_batch_help(self, temp_csv_file):
        """CSV with batch subcommand should display help without error."""
        from chemsmart.cli.gaussian.gaussian import gaussian

        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-p",
                "test_project",
                "-f",
                temp_csv_file,
                "pka",
                "-s",
                "direct",
                "batch",
                "--help",
            ],
        )

        # Help should show batch command documentation
        # Exit code 0 is expected for --help
        assert (
            result.exit_code == 0
        ), f"Batch help failed with exit code {result.exit_code}:\n{result.output}"
        assert (
            "batch" in result.output.lower()
        ), f"Expected 'batch' in help output:\n{result.output}"

    def test_csv_submit_help(self, temp_csv_file):
        """CSV with explicit submit subcommand should work (though not typical)."""
        from chemsmart.cli.gaussian.gaussian import gaussian

        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-p",
                "test_project",
                "-f",
                temp_csv_file,
                "pka",
                "-s",
                "direct",
                "submit",
                "--help",
            ],
        )

        # This is allowed - submit can be explicitly specified
        # But it will likely fail during actual execution (no proton_index)
        # For now we just check that the help works
        assert (
            result.exit_code == 0
        ), f"Submit help failed with exit code {result.exit_code}:\n{result.output}"


class TestOrcaPkaCSVHandling:
    """Test ORCA pKa CSV batch processing."""

    def test_csv_without_subcommand_raises_error(self, temp_csv_file):
        """CSV files without explicit subcommand should raise helpful error."""
        from chemsmart.cli.orca.orca import orca

        runner = CliRunner()
        result = runner.invoke(
            orca,
            [
                "-p",
                "test_project",
                "-f",
                temp_csv_file,
                "pka",
                "-s",
                "direct",
                # NOTE: No explicit subcommand (no 'batch' or 'submit')
            ],
        )

        assert result.exit_code != 0, "Should fail without explicit subcommand"
        assert (
            "CSV input files require an explicit subcommand" in result.output
        ), f"Expected helpful error message. Got:\n{result.output}"

    def test_csv_batch_help(self, temp_csv_file):
        """CSV with batch subcommand should display help without error."""
        from chemsmart.cli.orca.orca import orca

        runner = CliRunner()
        result = runner.invoke(
            orca,
            [
                "-p",
                "test_project",
                "-f",
                temp_csv_file,
                "pka",
                "-s",
                "direct",
                "batch",
                "--help",
            ],
        )

        # Help should show batch command documentation
        # Exit code 0 is expected for --help
        assert (
            result.exit_code == 0
        ), f"Batch help failed with exit code {result.exit_code}:\n{result.output}"
        assert (
            "batch" in result.output.lower()
        ), f"Expected 'batch' in help output:\n{result.output}"

    def test_csv_submit_help(self, temp_csv_file):
        """CSV with explicit submit subcommand should work (though not typical)."""
        from chemsmart.cli.orca.orca import orca

        runner = CliRunner()
        result = runner.invoke(
            orca,
            [
                "-p",
                "test_project",
                "-f",
                temp_csv_file,
                "pka",
                "-s",
                "direct",
                "submit",
                "--help",
            ],
        )

        # This is allowed - submit can be explicitly specified
        # But it will likely fail during actual execution (no proton_index)
        # For now we just check that the help works
        assert (
            result.exit_code == 0
        ), f"Submit help failed with exit code {result.exit_code}:\n{result.output}"


class TestNonCSVFilesStillWork:
    """Ensure non-CSV workflows are not affected by the fix."""

    def test_xyz_file_without_subcommand(self, temp_xyz_file):
        """XYZ files without subcommand should auto-invoke submit."""
        from chemsmart.cli.gaussian.gaussian import gaussian

        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-p",
                "test_project",
                "-f",
                temp_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-pi",
                "1",
                "-s",
                "direct",
                "--help",
            ],
        )

        # XYZ files should NOT raise the CSV error
        # --help on submit command should work
        if result.exit_code != 0:
            # If error, it should NOT be about CSV requirement
            assert (
                "CSV input files require an explicit subcommand"
                not in result.output
            ), f"XYZ file incorrectly rejected as CSV:\n{result.output}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
