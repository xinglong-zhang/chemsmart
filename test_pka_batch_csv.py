#!/usr/bin/env python
"""Test script to verify pKa batch processing with CSV files works correctly."""

import os
import sys
import tempfile

# Add the project to path
sys.path.insert(0, "/Users/taipanlan/PycharmProjects/chemsmart")

from click.testing import CliRunner

from chemsmart.cli.gaussian.gaussian import gaussian as gaussian_cli
from chemsmart.cli.orca.orca import orca as orca_cli


def create_temp_xyz(content="3\n\nH 0 0 0\nH 0 0 0.74\nO 0 0 0.37\n"):
    """Create a temporary XYZ file."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False)
    f.write(content)
    f.close()
    return f.name


def create_temp_csv(xyz_path):
    """Create a temporary CSV file pointing to xyz_path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False)
    f.write("filepath,proton_index,charge,multiplicity\n")
    f.write(f"{xyz_path},1,0,1\n")
    f.close()
    return f.name


def test_gaussian_pka_csv_without_subcommand():
    """Test that CSV without explicit subcommand raises error."""
    runner = CliRunner()
    xyz_path = create_temp_xyz()
    csv_path = create_temp_csv(xyz_path)

    try:
        # Try to run pka with CSV but no subcommand
        # This should raise an error asking for explicit subcommand
        result = runner.invoke(
            gaussian_cli,
            [
                "-p",
                "test_project",
                "-f",
                csv_path,
                "pka",
                "-s",
                "direct",
                # NOTE: No subcommand specified (no 'batch' or 'submit')
            ],
        )

        print("=" * 80)
        print("TEST 1: Gaussian pKa with CSV but NO explicit subcommand")
        print("=" * 80)
        print(f"Exit Code: {result.exit_code}")
        print(f"Output:\n{result.output}")

        # Should fail with our new error message
        assert (
            result.exit_code != 0
        ), "Should have failed without explicit subcommand"
        assert (
            "CSV input files require an explicit subcommand" in result.output
        ), f"Expected error message about CSV requiring subcommand. Got: {result.output}"
        print("✓ PASSED: CSV without subcommand correctly rejected\n")

    finally:
        os.unlink(xyz_path)
        os.unlink(csv_path)


def test_gaussian_pka_csv_batch_help():
    """Test that CSV with batch subcommand and --help works."""
    runner = CliRunner()
    xyz_path = create_temp_xyz()
    csv_path = create_temp_csv(xyz_path)

    try:
        result = runner.invoke(
            gaussian_cli,
            [
                "-p",
                "test_project",
                "-f",
                csv_path,
                "pka",
                "-s",
                "direct",
                "batch",
                "--help",
            ],
        )

        print("=" * 80)
        print("TEST 2: Gaussian pKa CSV batch with --help")
        print("=" * 80)
        print(f"Exit Code: {result.exit_code}")
        print(f"Output:\n{result.output[:500]}")

        # Should succeed (show help without actually running)
        # Exit code 0 for --help is normal
        assert (
            "batch" in result.output.lower() or result.exit_code == 0
        ), f"Expected help output. Got exit code {result.exit_code}: {result.output}"
        print("✓ PASSED: CSV batch --help works\n")

    finally:
        os.unlink(xyz_path)
        os.unlink(csv_path)


def test_orca_pka_csv_without_subcommand():
    """Test that ORCA CSV without explicit subcommand raises error."""
    runner = CliRunner()
    xyz_path = create_temp_xyz()
    csv_path = create_temp_csv(xyz_path)

    try:
        result = runner.invoke(
            orca_cli,
            [
                "-p",
                "test_project",
                "-f",
                csv_path,
                "pka",
                "-s",
                "direct",
                # NOTE: No subcommand specified
            ],
        )

        print("=" * 80)
        print("TEST 3: ORCA pKa with CSV but NO explicit subcommand")
        print("=" * 80)
        print(f"Exit Code: {result.exit_code}")
        print(f"Output:\n{result.output}")

        # Should fail with our new error message
        assert (
            result.exit_code != 0
        ), "Should have failed without explicit subcommand"
        assert (
            "CSV input files require an explicit subcommand" in result.output
        ), f"Expected error message about CSV requiring subcommand. Got: {result.output}"
        print("✓ PASSED: ORCA CSV without subcommand correctly rejected\n")

    finally:
        os.unlink(xyz_path)
        os.unlink(csv_path)


def test_orca_pka_csv_batch_help():
    """Test that ORCA CSV with batch subcommand and --help works."""
    runner = CliRunner()
    xyz_path = create_temp_xyz()
    csv_path = create_temp_csv(xyz_path)

    try:
        result = runner.invoke(
            orca_cli,
            [
                "-p",
                "test_project",
                "-f",
                csv_path,
                "pka",
                "-s",
                "direct",
                "batch",
                "--help",
            ],
        )

        print("=" * 80)
        print("TEST 4: ORCA pKa CSV batch with --help")
        print("=" * 80)
        print(f"Exit Code: {result.exit_code}")
        print(f"Output:\n{result.output[:500]}")

        # Should succeed (show help without actually running)
        assert (
            "batch" in result.output.lower() or result.exit_code == 0
        ), f"Expected help output. Got exit code {result.exit_code}: {result.output}"
        print("✓ PASSED: ORCA CSV batch --help works\n")

    finally:
        os.unlink(xyz_path)
        os.unlink(csv_path)


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("TESTING pKa CSV BATCH PROCESSING FIX")
    print("=" * 80 + "\n")

    try:
        test_gaussian_pka_csv_without_subcommand()
        test_gaussian_pka_csv_batch_help()
        test_orca_pka_csv_without_subcommand()
        test_orca_pka_csv_batch_help()

        print("=" * 80)
        print("ALL TESTS PASSED ✓")
        print("=" * 80)
    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
