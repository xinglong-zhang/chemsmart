import tempfile
from pathlib import Path

from click.testing import CliRunner

from chemsmart.cli.update import (
    DependencyUpdater,
    VersionUpdater,
    update,
)


def test_dependency_cli_delegates_to_dependency_updater(monkeypatch):
    called = False

    class FakeDependencyUpdater:
        def update(self):
            nonlocal called
            called = True

    monkeypatch.setattr(
        "chemsmart.cli.update.DependencyUpdater", FakeDependencyUpdater
    )

    result = CliRunner().invoke(update, ["deps"])

    assert result.exit_code == 0, result.output
    assert called


def test_version_cli_constructs_version_updater_with_cli_value(monkeypatch):
    seen_version = None
    called = False

    class FakeVersionUpdater:
        def __init__(self, version_number):
            nonlocal seen_version
            seen_version = version_number

        def update(self):
            nonlocal called
            called = True

    monkeypatch.setattr(
        "chemsmart.cli.update.VersionUpdater", FakeVersionUpdater
    )

    result = CliRunner().invoke(
        update, ["version", "--version-number", "9.8.7"]
    )

    assert result.exit_code == 0, result.output
    assert seen_version == "9.8.7"
    assert called


def test_version_updater_updates_version_files_in_order(tmp_path):
    version_file = tmp_path / "VERSION"
    pyproject_file = tmp_path / "pyproject.toml"
    docs_conf_file = tmp_path / "docs" / "source" / "conf.py"
    docs_conf_file.parent.mkdir(parents=True)

    version_file.write_text("0.1.9\n", encoding="utf-8")
    pyproject_file.write_text(
        'name = "chemsmart"\nversion = "0.1.9"  # keep pyproject comment\n',
        encoding="utf-8",
    )
    docs_conf_file.write_text(
        'project = "CHEMSMART"\nrelease = "0.1.9"  # keep docs comment\n',
        encoding="utf-8",
    )

    updater = VersionUpdater("4.5.6")
    updater._version_file_path = version_file
    updater._pyproject_path = pyproject_file
    updater._docs_conf_file_path = docs_conf_file

    updater.update()

    assert version_file.read_text(encoding="utf-8") == "4.5.6\n"
    assert (
        pyproject_file.read_text(encoding="utf-8")
        == 'name = "chemsmart"\nversion = "4.5.6"  # keep pyproject comment\n'
    )
    assert docs_conf_file.read_text(encoding="utf-8") == (
        'project = "CHEMSMART"\nrelease = "4.5.6"  # keep docs comment\n'
    )


class TestUpdater:
    """Tests for the Updater class focusing on dependency normalization."""

    def test_extract_pkg_name_normalization(self):
        """Test that package names are normalized correctly (PEP 503)."""
        updater = DependencyUpdater()

        # Create a temporary requirements file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".txt"
        ) as tmp:
            tmp.write("scikit_learn>=1.6.1\n")
            tmp.write("pandas>=2.0.0\n")
            tmp.write("some-package>=1.0.0\n")
            tmp_path = Path(tmp.name)

        try:
            # Mock the existing dependencies with dash notation
            existing_deps = {
                "scikit-learn>=1.6.0",  # Different version, dash notation
                "pandas>=2.0.0",  # Same version
                "some_package>=0.9.0",  # Different version, underscore notation
            }

            # Temporarily replace the method to return our mock dependencies
            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps

            try:
                # Get missing dependencies
                missing = updater._get_missing_dependencies(tmp_path)

                # Assert no duplicates are found
                # scikit_learn should match scikit-learn
                # pandas should match pandas
                # some-package should match some_package
                assert (
                    len(missing) == 0
                ), f"Expected no missing deps, got: {missing}"
            finally:
                # Restore original method
                updater._get_existing_dependencies = original_method
        finally:
            # Always cleanup temporary file
            tmp_path.unlink()

    def test_underscore_dash_equivalence(self):
        """Test that underscore and dash package
        names are treated as equivalent."""
        updater = DependencyUpdater()

        # Create a temporary requirements file with underscore
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".txt"
        ) as tmp:
            tmp.write("pytest_mock>=3.14.0\n")
            tmp_path = Path(tmp.name)

        try:
            # Existing dependency with dash
            existing_deps = {"pytest-mock>=3.0.0"}

            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps

            try:
                missing = updater._get_missing_dependencies(tmp_path)

                # Should not consider pytest_mock as
                # missing since pytest-mock exists
                assert (
                    len(missing) == 0
                ), f"Expected no missing deps, got: {missing}"
            finally:
                updater._get_existing_dependencies = original_method
        finally:
            tmp_path.unlink()

    def test_dash_underscore_equivalence(self):
        """Test that dash and underscore package
        names are treated as equivalent."""
        updater = DependencyUpdater()

        # Create a temporary requirements file with dash
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".txt"
        ) as tmp:
            tmp.write("scikit-learn>=1.6.1\n")
            tmp_path = Path(tmp.name)

        try:
            # Existing dependency with underscore
            existing_deps = {"scikit_learn>=1.5.0"}

            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps

            try:
                missing = updater._get_missing_dependencies(tmp_path)

                # Should not consider scikit-learn as
                # missing since scikit_learn exists
                assert (
                    len(missing) == 0
                ), f"Expected no missing deps, got: {missing}"
            finally:
                updater._get_existing_dependencies = original_method
        finally:
            tmp_path.unlink()

    def test_new_dependency_detected(self):
        """Test that genuinely new dependencies are still detected."""
        updater = DependencyUpdater()

        # Create a temporary requirements file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".txt"
        ) as tmp:
            tmp.write("new-package>=1.0.0\n")
            tmp.write("scikit-learn>=1.6.1\n")
            tmp_path = Path(tmp.name)

        try:
            # Existing dependencies without new-package
            existing_deps = {"scikit-learn>=1.6.0"}

            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps

            try:
                missing = updater._get_missing_dependencies(tmp_path)

                # Should detect new-package as missing, but not scikit-learn
                assert len(missing) == 1
                assert "new-package>=1.0.0" in missing
            finally:
                updater._get_existing_dependencies = original_method
        finally:
            tmp_path.unlink()
