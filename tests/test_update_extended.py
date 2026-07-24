"""
Extended direct unit tests for :class:`chemsmart.cli.update.Updater` and
the ``update`` CLI group, beyond what ``test_update.py`` already covers
(package-name normalization in ``_get_missing_dependencies``).

Every test that needs a pyproject.toml/environment.yml/VERSION/conf.py
redirects the ``Updater`` instance's private path attributes to files
under ``tmp_path`` — never the real project files.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.update import Updater, update


@pytest.fixture()
def updater(tmp_path):
    u = Updater()
    u._package_path = tmp_path
    u._pyproject_path = tmp_path / "pyproject.toml"
    u._version_file_path = tmp_path / "VERSION"
    u._docs_conf_file_path = tmp_path / "conf.py"
    return u


class TestUpdaterPathProperties:
    def test_properties_expose_private_attrs(self, updater, tmp_path):
        assert updater.package_path == tmp_path
        assert updater.pyproject_path == tmp_path / "pyproject.toml"
        assert updater.version_file_path == tmp_path / "VERSION"
        assert updater.docs_conf_file_path == tmp_path / "conf.py"


class TestGetExistingDependenciesFromPyproject:
    def test_missing_file_returns_empty_set(self, updater):
        assert updater._get_existing_dependencies_from_pyproject() == set()

    def test_reads_dependencies_from_project_table(self, updater):
        updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy>=1.0", "click"]\n'
        )
        deps = updater._get_existing_dependencies_from_pyproject()
        assert deps == {"numpy>=1.0", "click"}

    def test_missing_dependencies_key_returns_empty_set(self, updater):
        updater.pyproject_path.write_text('[project]\nname = "chemsmart"\n')
        assert updater._get_existing_dependencies_from_pyproject() == set()


class TestGetExistingDependenciesFromEnvironmentYml:
    def test_missing_file_returns_empty_set(
        self, updater, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        assert (
            updater._get_existing_dependencies_from_environment_yml() == set()
        )

    def test_reads_plain_string_deps(self, updater, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "environment.yml").write_text(
            "dependencies:\n  - numpy\n  - pandas\n"
        )
        deps = updater._get_existing_dependencies_from_environment_yml()
        assert deps == {"numpy", "pandas"}

    def test_reads_pip_subsection_deps(self, updater, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "environment.yml").write_text(
            "dependencies:\n"
            "  - numpy\n"
            "  - pip:\n"
            "    - click>=8.0\n"
            "    - requests\n"
        )
        deps = updater._get_existing_dependencies_from_environment_yml()
        assert deps == {"numpy", "click>=8.0", "requests"}


class TestGetExistingDependencies:
    def test_combines_pyproject_and_environment_yml(
        self, updater, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy"]\n'
        )
        (tmp_path / "environment.yml").write_text(
            "dependencies:\n  - pandas\n"
        )
        deps = updater._get_existing_dependencies()
        assert deps == {"numpy", "pandas"}


class TestGetMissingDependenciesFileHandling:
    def test_missing_requirements_file_returns_empty_set(
        self, updater, tmp_path
    ):
        missing = updater._get_missing_dependencies(
            tmp_path / "does_not_exist.txt"
        )
        assert missing == set()


class TestUpdateToml:
    def test_no_missing_dependencies_is_noop(self, updater, tmp_path):
        updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy"]\n'
        )
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("numpy\n")

        original_text = updater.pyproject_path.read_text()
        updater._update_toml(req_file)
        assert updater.pyproject_path.read_text() == original_text

    def test_appends_missing_dependency_with_gte_normalized(
        self, updater, tmp_path
    ):
        updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy"]\n'
        )
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("Click==8.1.0\n")

        updater._update_toml(req_file)

        new_text = updater.pyproject_path.read_text()
        assert "click>=8.1.0" in new_text
        assert "numpy" in new_text


class TestGenerateRequirements:
    def test_returns_none_on_pipreqs_failure(self, updater):
        mock_process = MagicMock()
        mock_process.communicate.return_value = (b"", b"pipreqs error")
        mock_process.returncode = 1
        with patch(
            "chemsmart.cli.update.subprocess.Popen",
            return_value=mock_process,
        ):
            result = updater._generate_requirements()
        assert result is None

    def test_returns_tmp_path_on_success(self, updater):
        mock_process = MagicMock()
        mock_process.communicate.return_value = (b"done", b"")
        mock_process.returncode = 0
        with patch(
            "chemsmart.cli.update.subprocess.Popen",
            return_value=mock_process,
        ):
            result = updater._generate_requirements()
        assert isinstance(result, Path)


class TestUpdatePyprojectToml:
    def test_no_requirements_path_skips_toml_update(self, updater):
        with (
            patch.object(updater, "_generate_requirements", return_value=None),
            patch.object(updater, "_update_toml") as mock_update_toml,
        ):
            updater.update_pyproject_toml()
        mock_update_toml.assert_not_called()

    def test_requirements_path_triggers_toml_update(self, updater, tmp_path):
        req_path = tmp_path / "req.txt"
        with (
            patch.object(
                updater, "_generate_requirements", return_value=req_path
            ),
            patch.object(updater, "_update_toml") as mock_update_toml,
        ):
            updater.update_pyproject_toml()
        mock_update_toml.assert_called_once_with(req_path)


class TestUpdateVersionNumber:
    def test_missing_version_file_logs_error_but_continues(
        self, updater, caplog
    ):
        # version_file_path doesn't exist; pyproject/conf also don't
        # exist, so this should log errors/warnings but not raise.
        with caplog.at_level("ERROR"):
            updater.update_version_number("9.9.9")
        assert "Version file not found" in caplog.text

    def test_writes_version_file(self, updater):
        updater.version_file_path.write_text("0.0.1\n")
        updater.update_version_number("1.2.3")
        assert updater.version_file_path.read_text() == "1.2.3\n"

    def test_updates_pyproject_version_field(self, updater):
        updater.pyproject_path.write_text(
            '[project]\nname = "chemsmart"\nversion = "0.1.0"\n'
        )
        updater.update_version_number("2.0.0")
        new_text = updater.pyproject_path.read_text()
        assert 'version = "2.0.0"' in new_text

    def test_pyproject_without_version_field_warns(self, updater, caplog):
        updater.pyproject_path.write_text('[project]\nname = "chemsmart"\n')
        with caplog.at_level("WARNING"):
            updater.update_version_number("2.0.0")
        assert "not updated" in caplog.text

    def test_updates_docs_conf_release_field(self, updater):
        updater.docs_conf_file_path.write_text('release = "0.1.0"\n')
        updater.update_version_number("2.0.0")
        new_text = updater.docs_conf_file_path.read_text()
        assert 'release = "2.0.0"' in new_text

    def test_missing_docs_conf_warns_but_continues(self, updater, caplog):
        with caplog.at_level("WARNING"):
            updater.update_version_number("2.0.0")
        assert "skipping docs version update" in caplog.text


class TestUpdateCLI:
    def test_deps_command_calls_update_pyproject_toml(self):
        runner = CliRunner()
        with patch("chemsmart.cli.update.Updater") as mock_updater_cls:
            mock_updater = MagicMock()
            mock_updater_cls.return_value = mock_updater
            result = runner.invoke(update, ["deps"], obj={})

        assert result.exit_code == 0, result.output
        mock_updater.update_pyproject_toml.assert_called_once()

    def test_version_command_requires_version_number(self):
        runner = CliRunner()
        result = runner.invoke(update, ["version"], obj={})
        assert result.exit_code != 0

    def test_version_command_calls_update_version_number(self):
        runner = CliRunner()
        with patch("chemsmart.cli.update.Updater") as mock_updater_cls:
            mock_updater = MagicMock()
            mock_updater_cls.return_value = mock_updater
            result = runner.invoke(update, ["version", "-v", "3.2.1"], obj={})

        assert result.exit_code == 0, result.output
        mock_updater.update_version_number.assert_called_once_with("3.2.1")
