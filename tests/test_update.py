import builtins
import os
import tempfile
from collections.abc import Mapping
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

import chemsmart.cli.update as update_module
from chemsmart.cli.update import (
    ConfigUpdateError,
    ConfigUpdateReport,
    ConfigurationUpdater,
    DependencyUpdater,
    ProjectsUpdater,
    ProjectUpdateError,
    Updater,
    VersionUpdater,
    _PreparedConfig,
    update,
)


class _FakeYamlParser:
    """Stand-in for the real ruamel.yaml.YAML() parser instance returned
    by ConfigurationUpdater._yaml(), for exercising load()'s error and
    not-a-Mapping branches without needing malformed real files."""

    def __init__(self, load_result=None, load_error=None):
        self._result = load_result
        self._error = load_error

    def load(self, text):
        if self._error is not None:
            raise self._error
        return self._result


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


class TestConfigsUpdater:
    def test_update_config_updates_selected_files_by_missing_program(
        self,
        tmp_path,
        monkeypatch,
        write_server_yaml,
        invoke_update_configs,
        read_yaml_file,
    ):
        fake_home = tmp_path / "home"
        monkeypatch.setenv("HOME", str(fake_home))
        template = {
            "SERVER": {"TEMPLATE_VALUE": "should-not-be-copied"},
            "FAKE_PROGRAM_A": {"EXEFOLDER": None, "LOCAL_RUN": True},
            "FAKE_PROGRAM_B": {"EXEFOLDER": None, "LOCAL_RUN": True},
            "FAKE_PROGRAM_C": {"EXEFOLDER": None, "LOCAL_RUN": True},
        }
        monkeypatch.setattr(
            ConfigurationUpdater,
            "_load_server_template",
            lambda self: template,
        )
        monkeypatch.setattr(
            ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
        )

        a_yaml = write_server_yaml(
            tmp_path,
            "A.yaml",
            """
            SERVER:
                USER_VALUE: preserved
            FAKE_PROGRAM_B:
                EXEFOLDER: /old/a-fake-b
            FAKE_PROGRAM_C:
                EXEFOLDER: /old/a-fake-c
            """,
        )
        b_yaml = write_server_yaml(
            tmp_path,
            "B.yaml",
            """
            SERVER:
                USER_VALUE: preserved
            FAKE_PROGRAM_A:
                EXEFOLDER: /old/b-fake-a
            FAKE_PROGRAM_C:
                EXEFOLDER: /old/b-fake-c
            """,
        )
        c_yaml = write_server_yaml(
            tmp_path,
            "C.yaml",
            """
            SERVER:
                USER_VALUE: preserved
            """,
        )
        unselected = write_server_yaml(
            tmp_path,
            "unselected.yaml",
            """
            SERVER:
                USER_VALUE: preserved
            """,
        )
        unselected_before = unselected.read_text(encoding="utf-8")

        prompts = []
        answers = iter(["~/path/to/fake-a", "nUlL", ""])

        def prompt_for_program(text, **kwargs):
            prompts.append(text)
            return next(answers)

        monkeypatch.setattr(update_module.click, "prompt", prompt_for_program)

        result = invoke_update_configs(
            tmp_path, ["-s", "A", "-s", "B.yaml", "-s", "C"]
        )

        assert result.exit_code == 0, result.output
        assert "A.yaml:" in result.output
        assert "B.yaml:" in result.output
        assert "C.yaml:" in result.output
        assert "unselected.yaml:" not in result.output
        assert len(prompts) == 3
        assert prompts[0].startswith("FAKE_PROGRAM_A EXEFOLDER")
        assert prompts[1].startswith("FAKE_PROGRAM_B EXEFOLDER")
        assert prompts[2].startswith("FAKE_PROGRAM_C EXEFOLDER")
        assert not any(
            filename in prompt
            for prompt in prompts
            for filename in ("A.yaml", "B.yaml", "C.yaml")
        )

        a_data = read_yaml_file(a_yaml)
        b_data = read_yaml_file(b_yaml)
        c_data = read_yaml_file(c_yaml)
        expanded_fake_a = str(fake_home / "path" / "to" / "fake-a")
        assert a_data["SERVER"] == {"USER_VALUE": "preserved"}
        assert a_data["FAKE_PROGRAM_A"] == {
            "EXEFOLDER": expanded_fake_a,
            "LOCAL_RUN": True,
        }
        assert a_data["FAKE_PROGRAM_B"] == {"EXEFOLDER": "/old/a-fake-b"}
        assert a_data["FAKE_PROGRAM_C"] == {"EXEFOLDER": "/old/a-fake-c"}
        assert b_data["FAKE_PROGRAM_A"] == {"EXEFOLDER": "/old/b-fake-a"}
        assert b_data["FAKE_PROGRAM_B"] == {
            "EXEFOLDER": None,
            "LOCAL_RUN": True,
        }
        assert b_data["FAKE_PROGRAM_C"] == {"EXEFOLDER": "/old/b-fake-c"}
        assert c_data["FAKE_PROGRAM_A"] == {
            "EXEFOLDER": expanded_fake_a,
            "LOCAL_RUN": True,
        }
        assert c_data["FAKE_PROGRAM_B"] == {
            "EXEFOLDER": None,
            "LOCAL_RUN": True,
        }
        assert c_data["FAKE_PROGRAM_C"] == {
            "EXEFOLDER": None,
            "LOCAL_RUN": True,
        }
        assert "EXEFOLDER: Null" in b_yaml.read_text(encoding="utf-8")
        c_text = c_yaml.read_text(encoding="utf-8")
        assert c_text.count("EXEFOLDER: Null") == 2
        assert "&id" not in c_text
        assert "*id" not in c_text
        assert unselected.read_text(encoding="utf-8") == unselected_before

    def test_update_config_missing_ruamel_shows_install_hint(
        self, tmp_path, monkeypatch, invoke_update_configs
    ):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "ruamel.yaml":
                raise ModuleNotFoundError(name="ruamel")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)

        result = invoke_update_configs(tmp_path)

        assert result.exit_code != 0
        assert "ruamel.yaml is required" in result.output
        assert "pip install" in result.output
        assert "make env" in result.output
        assert "Traceback" not in result.output


class TestProjectsUpdater:
    def test_update_projects_adds_missing_template_content_without_overwriting(
        self, tmp_path, monkeypatch, write_text_file, invoke_update_projects
    ):
        assert ProjectsUpdater._template_config_dir().is_dir()

        template_root = tmp_path / "templates" / ".chemsmart"
        config_root = tmp_path / "user" / ".chemsmart"

        write_text_file(
            template_root / "fake_a" / "existing.yaml", "template existing\n"
        )
        write_text_file(
            template_root / "fake_a" / "missing.yaml", "template missing\n"
        )
        write_text_file(
            template_root / "fake_a" / "nested" / "nested.yaml",
            "template nested\n",
        )
        write_text_file(
            template_root / "fake_b" / "new.yaml", "template new\n"
        )
        (template_root / "fake_empty").mkdir(parents=True)
        write_text_file(
            template_root / "server" / "ignored.yaml", "server template\n"
        )
        write_text_file(
            template_root / "usersettings.yaml", "user settings template\n"
        )

        existing = write_text_file(
            config_root / "fake_a" / "existing.yaml", "user existing\n"
        )
        custom = write_text_file(
            config_root / "fake_a" / "custom.yaml", "user custom\n"
        )
        personal = write_text_file(
            config_root / "custom_project" / "personal.yaml", "user personal\n"
        )

        monkeypatch.setattr(
            ProjectsUpdater,
            "_template_config_dir",
            staticmethod(lambda: template_root),
        )

        result = invoke_update_projects(config_root)

        assert result.exit_code == 0, result.output
        assert "added missing content: fake_a" in result.output
        assert "added missing content: fake_b" in result.output
        assert "added missing content: fake_empty" in result.output
        assert existing.read_text(encoding="utf-8") == "user existing\n"
        assert custom.read_text(encoding="utf-8") == "user custom\n"
        assert personal.read_text(encoding="utf-8") == "user personal\n"
        assert (config_root / "fake_a" / "missing.yaml").read_text(
            encoding="utf-8"
        ) == "template missing\n"
        assert (config_root / "fake_a" / "nested" / "nested.yaml").read_text(
            encoding="utf-8"
        ) == "template nested\n"
        assert (config_root / "fake_b" / "new.yaml").read_text(
            encoding="utf-8"
        ) == "template new\n"
        assert (config_root / "fake_empty").is_dir()
        assert not (config_root / "server").exists()
        assert not (config_root / "usersettings.yaml").exists()

        files_before = {
            path.relative_to(config_root): path.read_bytes()
            for path in config_root.rglob("*")
            if path.is_file()
        }
        repeated = invoke_update_projects(config_root)
        files_after = {
            path.relative_to(config_root): path.read_bytes()
            for path in config_root.rglob("*")
            if path.is_file()
        }

        assert repeated.exit_code == 0, repeated.output
        assert "already up to date" in repeated.output
        assert files_after == files_before

    def test_update_projects_reports_errors_without_damaging_user_content(
        self, tmp_path, monkeypatch, write_text_file, invoke_update_projects
    ):
        template_root = tmp_path / "templates" / ".chemsmart"
        config_root = tmp_path / "user" / ".chemsmart"
        write_text_file(
            template_root / "fake_a" / "settings.yaml", "template value\n"
        )
        conflicting_path = write_text_file(
            config_root / "fake_a", "user value\n"
        )

        monkeypatch.setattr(
            ProjectsUpdater,
            "_template_config_dir",
            staticmethod(lambda: template_root),
        )

        result = invoke_update_projects(config_root)

        assert result.exit_code != 0
        assert "Expected a directory" in result.output
        assert "Traceback" not in result.output
        assert conflicting_path.read_text(encoding="utf-8") == "user value\n"

        copy_config_root = tmp_path / "copy-user" / ".chemsmart"
        custom = write_text_file(
            copy_config_root / "custom.yaml", "user custom\n"
        )

        def fail_during_copy(source_file, destination_file):
            destination_file.write(b"partial content")
            raise OSError("copy failed")

        monkeypatch.setattr(
            update_module.shutil, "copyfileobj", fail_during_copy
        )

        copy_result = invoke_update_projects(copy_config_root)

        assert copy_result.exit_code != 0
        assert "copy failed" in copy_result.output
        assert not (copy_config_root / "fake_a" / "settings.yaml").exists()
        assert custom.read_text(encoding="utf-8") == "user custom\n"


class TestUpdaterAbstractBase:
    def test_base_update_is_a_no_op_pass(self, dependency_updater):
        # Updater.update() is @abstractmethod with a `pass` body; concrete
        # subclasses always override it, so this only runs by calling the
        # base implementation directly against a subclass instance.
        assert Updater.update(dependency_updater) is None


class TestUpdaterPathProperties:
    def test_dependency_updater_properties_expose_private_attrs(
        self, dependency_updater, tmp_path
    ):
        assert dependency_updater.package_path == tmp_path
        assert dependency_updater.pyproject_path == tmp_path / "pyproject.toml"

    def test_version_updater_properties_expose_private_attrs(
        self, version_updater, tmp_path
    ):
        assert version_updater.pyproject_path == tmp_path / "pyproject.toml"
        assert version_updater.version_file_path == tmp_path / "VERSION"
        assert version_updater.docs_conf_file_path == tmp_path / "conf.py"


class TestGetExistingDependenciesFromPyproject:
    def test_missing_file_returns_empty_set(self, dependency_updater):
        assert (
            dependency_updater._get_existing_dependencies_from_pyproject()
            == set()
        )

    def test_reads_dependencies_from_project_table(self, dependency_updater):
        dependency_updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy>=1.0", "click"]\n'
        )
        deps = dependency_updater._get_existing_dependencies_from_pyproject()
        assert deps == {"numpy>=1.0", "click"}

    def test_missing_dependencies_key_returns_empty_set(
        self, dependency_updater
    ):
        dependency_updater.pyproject_path.write_text(
            '[project]\nname = "chemsmart"\n'
        )
        assert (
            dependency_updater._get_existing_dependencies_from_pyproject()
            == set()
        )


class TestGetExistingDependenciesFromEnvironmentYml:
    def test_missing_file_returns_empty_set(
        self, dependency_updater, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        assert (
            dependency_updater._get_existing_dependencies_from_environment_yml()
            == set()
        )

    def test_reads_plain_string_deps(
        self, dependency_updater, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "environment.yml").write_text(
            "dependencies:\n  - numpy\n  - pandas\n"
        )
        deps = (
            dependency_updater._get_existing_dependencies_from_environment_yml()
        )
        assert deps == {"numpy", "pandas"}

    def test_reads_pip_subsection_deps(
        self, dependency_updater, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "environment.yml").write_text(
            "dependencies:\n"
            "  - numpy\n"
            "  - pip:\n"
            "    - click>=8.0\n"
            "    - requests\n"
        )
        deps = (
            dependency_updater._get_existing_dependencies_from_environment_yml()
        )
        assert deps == {"numpy", "click>=8.0", "requests"}


class TestGetExistingDependencies:
    def test_combines_pyproject_and_environment_yml(
        self, dependency_updater, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        dependency_updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy"]\n'
        )
        (tmp_path / "environment.yml").write_text(
            "dependencies:\n  - pandas\n"
        )
        deps = dependency_updater._get_existing_dependencies()
        assert deps == {"numpy", "pandas"}


class TestGetMissingDependenciesFileHandling:
    def test_missing_requirements_file_returns_empty_set(
        self, dependency_updater, tmp_path
    ):
        missing = dependency_updater._get_missing_dependencies(
            tmp_path / "does_not_exist.txt"
        )
        assert missing == set()


class TestUpdateToml:
    def test_no_missing_dependencies_is_noop(
        self, dependency_updater, tmp_path
    ):
        dependency_updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy"]\n'
        )
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("numpy\n")

        original_text = dependency_updater.pyproject_path.read_text()
        dependency_updater._update_toml(req_file)
        assert dependency_updater.pyproject_path.read_text() == original_text

    def test_appends_missing_dependency_with_gte_normalized(
        self, dependency_updater, tmp_path
    ):
        dependency_updater.pyproject_path.write_text(
            '[project]\ndependencies = ["numpy"]\n'
        )
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("Click==8.1.0\n")

        dependency_updater._update_toml(req_file)

        new_text = dependency_updater.pyproject_path.read_text()
        assert "click>=8.1.0" in new_text
        assert "numpy" in new_text


class TestGenerateRequirements:
    def test_returns_none_on_pipreqs_failure(self, dependency_updater):
        mock_process = MagicMock()
        mock_process.communicate.return_value = (b"", b"pipreqs error")
        mock_process.returncode = 1
        with patch(
            "chemsmart.cli.update.subprocess.Popen",
            return_value=mock_process,
        ):
            result = dependency_updater._generate_requirements()
        assert result is None

    def test_returns_tmp_path_on_success(self, dependency_updater):
        mock_process = MagicMock()
        mock_process.communicate.return_value = (b"done", b"")
        mock_process.returncode = 0
        with patch(
            "chemsmart.cli.update.subprocess.Popen",
            return_value=mock_process,
        ):
            result = dependency_updater._generate_requirements()
        assert isinstance(result, Path)


class TestDependencyUpdaterUpdateEntryPoint:
    def test_update_delegates_to_update_pyproject_toml(
        self, dependency_updater
    ):
        with patch.object(
            dependency_updater, "update_pyproject_toml"
        ) as mock_update_pyproject:
            dependency_updater.update()
        mock_update_pyproject.assert_called_once_with()


class TestUpdateTomlNoVersionSpecifier:
    def test_dependency_without_version_specifier_normalized_lowercase(
        self, dependency_updater, tmp_path
    ):
        dependency_updater.pyproject_path.write_text(
            "[project]\ndependencies = []\n"
        )
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("SomePackage\n")

        dependency_updater._update_toml(req_file)

        new_text = dependency_updater.pyproject_path.read_text()
        assert "somepackage" in new_text


class TestUpdatePyprojectToml:
    def test_no_requirements_path_skips_toml_update(self, dependency_updater):
        with (
            patch.object(
                dependency_updater, "_generate_requirements", return_value=None
            ),
            patch.object(
                dependency_updater, "_update_toml"
            ) as mock_update_toml,
        ):
            dependency_updater.update_pyproject_toml()
        mock_update_toml.assert_not_called()

    def test_requirements_path_triggers_toml_update(
        self, dependency_updater, tmp_path
    ):
        req_path = tmp_path / "req.txt"
        with (
            patch.object(
                dependency_updater,
                "_generate_requirements",
                return_value=req_path,
            ),
            patch.object(
                dependency_updater, "_update_toml"
            ) as mock_update_toml,
        ):
            dependency_updater.update_pyproject_toml()
        mock_update_toml.assert_called_once_with(req_path)


class TestUpdateVersionNumber:
    def test_missing_version_file_logs_error_but_continues(
        self, version_updater, caplog
    ):
        # version_file_path doesn't exist; pyproject/conf also don't
        # exist, so this should log errors/warnings but not raise.
        with caplog.at_level("ERROR"):
            version_updater.update_version_number("9.9.9")
        assert "Version file not found" in caplog.text

    def test_writes_version_file(self, version_updater):
        version_updater.version_file_path.write_text("0.0.1\n")
        version_updater.update_version_number("1.2.3")
        assert version_updater.version_file_path.read_text() == "1.2.3\n"

    def test_updates_pyproject_version_field(self, version_updater):
        version_updater.pyproject_path.write_text(
            '[project]\nname = "chemsmart"\nversion = "0.1.0"\n'
        )
        version_updater.update_version_number("2.0.0")
        new_text = version_updater.pyproject_path.read_text()
        assert 'version = "2.0.0"' in new_text

    def test_pyproject_without_version_field_warns(
        self, version_updater, caplog
    ):
        version_updater.pyproject_path.write_text(
            '[project]\nname = "chemsmart"\n'
        )
        with caplog.at_level("WARNING"):
            version_updater.update_version_number("2.0.0")
        assert "not updated" in caplog.text

    def test_updates_docs_conf_release_field(self, version_updater):
        version_updater.docs_conf_file_path.write_text('release = "0.1.0"\n')
        version_updater.update_version_number("2.0.0")
        new_text = version_updater.docs_conf_file_path.read_text()
        assert 'release = "2.0.0"' in new_text

    def test_missing_docs_conf_warns_but_continues(
        self, version_updater, caplog
    ):
        with caplog.at_level("WARNING"):
            version_updater.update_version_number("2.0.0")
        assert "skipping docs version update" in caplog.text

    def test_docs_conf_without_release_field_warns(
        self, version_updater, caplog
    ):
        version_updater.docs_conf_file_path.write_text(
            'project = "CHEMSMART"\n'
        )
        with caplog.at_level("WARNING"):
            version_updater.update_version_number("2.0.0")
        assert "Sphinx version not updated" in caplog.text


class TestProjectsUpdaterSymlinkAndTypeErrors:
    def test_ensure_directory_rejects_symlinked_config_dir(self, tmp_path):
        real_dir = tmp_path / "real"
        real_dir.mkdir()
        link = tmp_path / "link"
        link.symlink_to(real_dir, target_is_directory=True)

        with pytest.raises(ProjectUpdateError, match="Refusing to update"):
            ProjectsUpdater._ensure_directory(link)

    def test_copy_missing_rejects_symlinked_source(self, tmp_path):
        real_file = tmp_path / "real.yaml"
        real_file.write_text("content\n")
        source_link = tmp_path / "source_link.yaml"
        source_link.symlink_to(real_file)
        destination = tmp_path / "dest" / "source_link.yaml"

        with pytest.raises(ProjectUpdateError, match="Refusing to update"):
            ProjectsUpdater()._copy_missing(source_link, destination)

    def test_copy_missing_rejects_symlinked_destination(self, tmp_path):
        source = tmp_path / "source.yaml"
        source.write_text("content\n")
        real_dest_target = tmp_path / "elsewhere.yaml"
        real_dest_target.write_text("other\n")
        destination = tmp_path / "dest_link.yaml"
        destination.symlink_to(real_dest_target)

        with pytest.raises(ProjectUpdateError, match="Refusing to update"):
            ProjectsUpdater()._copy_missing(source, destination)

    def test_copy_missing_rejects_unsupported_source_entry(
        self, tmp_path, monkeypatch
    ):
        source = tmp_path / "fifo_like"
        source.write_text("")  # exists, but we'll fake it as neither
        destination = tmp_path / "dest" / "fifo_like"

        monkeypatch.setattr(Path, "is_dir", lambda self: False)
        monkeypatch.setattr(Path, "is_file", lambda self: False)

        with pytest.raises(
            ProjectUpdateError, match="Unsupported project template entry"
        ):
            ProjectsUpdater()._copy_missing(source, destination)

    def test_copy_missing_rejects_destination_that_is_a_directory(
        self, tmp_path
    ):
        source = tmp_path / "source.yaml"
        source.write_text("content\n")
        destination = tmp_path / "dest.yaml"
        destination.mkdir()

        with pytest.raises(ProjectUpdateError, match="Expected a file"):
            ProjectsUpdater()._copy_missing(source, destination)

    def test_copy_missing_cleans_up_partial_file_when_open_fails_early(
        self, tmp_path, monkeypatch
    ):
        """When destination.open('xb') itself fails (before any bytes are
        written), `created` is still False, so the except handler skips
        straight to re-raising without attempting an unlink."""
        source = tmp_path / "source.yaml"
        source.write_text("content\n")
        destination = tmp_path / "dest.yaml"

        real_open = Path.open

        def failing_open(self, mode="r", *args, **kwargs):
            if "x" in mode:
                raise OSError("cannot create destination")
            return real_open(self, mode, *args, **kwargs)

        monkeypatch.setattr(Path, "open", failing_open)

        with pytest.raises(OSError, match="cannot create destination"):
            ProjectsUpdater()._copy_missing(source, destination)
        assert not destination.exists()

    def test_copy_missing_ignores_filenotfound_when_cleaning_up(
        self, tmp_path, monkeypatch
    ):
        """If the partially-written destination is already gone by the
        time the except handler tries to unlink it, that FileNotFoundError
        is swallowed and the original error still propagates."""
        source = tmp_path / "source.yaml"
        source.write_text("content\n")
        destination = tmp_path / "dest.yaml"

        def failing_copystat(src, dst, *args, **kwargs):
            os.unlink(dst)
            raise OSError("copystat failed")

        monkeypatch.setattr(update_module.shutil, "copystat", failing_copystat)

        with pytest.raises(OSError, match="copystat failed"):
            ProjectsUpdater()._copy_missing(source, destination)
        assert not destination.exists()


class TestConfigurationUpdaterInit:
    def test_none_server_gives_empty_tuple(self):
        assert ConfigurationUpdater(server=None).servers == ()

    def test_string_server_becomes_single_element_tuple(self):
        assert ConfigurationUpdater(server="A.yaml").servers == ("A.yaml",)

    def test_sequence_server_becomes_tuple(self):
        updater = ConfigurationUpdater(server=["A.yaml", "B.yaml"])
        assert updater.servers == ("A.yaml", "B.yaml")


class TestSelectServerYamlFiles:
    def test_duplicate_server_names_are_deduplicated(self, tmp_path):
        server_dir = tmp_path / "server"
        server_dir.mkdir()
        (server_dir / "A.yaml").write_text("SERVER: {}\n")

        updater = ConfigurationUpdater(server=["A.yaml", "A.yaml", "A"])
        selected = updater._select_server_yaml_files(
            server_dir, updater.servers
        )
        assert selected == [server_dir / "A.yaml"]

    def test_missing_named_server_raises(self, tmp_path):
        server_dir = tmp_path / "server"
        server_dir.mkdir()

        updater = ConfigurationUpdater(server="missing")
        with pytest.raises(ConfigUpdateError, match="not found"):
            updater._select_server_yaml_files(server_dir, updater.servers)

    def test_missing_server_directory_raises(self, tmp_path):
        server_dir = tmp_path / "does_not_exist"
        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="not found"):
            updater._select_server_yaml_files(server_dir, ())

    def test_glob_oserror_wrapped_in_config_update_error(
        self, tmp_path, monkeypatch
    ):
        server_dir = tmp_path / "server"
        server_dir.mkdir()

        def failing_glob(self, pattern):
            raise OSError("glob failed")

        monkeypatch.setattr(Path, "glob", failing_glob)

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="Could not read"):
            updater._select_server_yaml_files(server_dir, ())


class TestNormalizeServerFilename:
    def test_empty_server_raises(self):
        with pytest.raises(ConfigUpdateError, match="must not be empty"):
            ConfigurationUpdater._normalize_server_filename("")

    def test_path_with_directory_component_raises(self):
        with pytest.raises(ConfigUpdateError, match="not a path"):
            ConfigurationUpdater._normalize_server_filename("sub/A.yaml")

    def test_wrong_suffix_raises(self):
        with pytest.raises(ConfigUpdateError, match="must name a .yaml"):
            ConfigurationUpdater._normalize_server_filename("A.txt")

    def test_appends_yaml_suffix_when_absent(self):
        assert ConfigurationUpdater._normalize_server_filename("A") == (
            "A.yaml"
        )

    def test_leaves_yaml_suffix_as_is(self):
        assert ConfigurationUpdater._normalize_server_filename("A.yaml") == (
            "A.yaml"
        )


class TestLoadServerTemplate:
    def test_loads_real_bundled_template_successfully(self):
        data = ConfigurationUpdater()._load_server_template()
        assert isinstance(data, Mapping)
        assert len(data) > 0

    def test_yaml_error_wrapped_in_config_update_error(self, monkeypatch):
        from ruamel.yaml.error import YAMLError

        updater = ConfigurationUpdater()
        monkeypatch.setattr(
            updater,
            "_yaml",
            lambda: _FakeYamlParser(load_error=YAMLError("bad yaml")),
        )
        with pytest.raises(ConfigUpdateError, match="not valid YAML"):
            updater._load_server_template()

    def test_read_text_oserror_wrapped_in_config_update_error(
        self, monkeypatch
    ):
        class _FakeTemplatePath:
            def __truediv__(self, other):
                return self

            def read_text(self, encoding="utf-8"):
                raise OSError("cannot read bundled template")

            @property
            def name(self):
                return "server.yaml"

        monkeypatch.setattr(
            update_module.resources, "files", lambda pkg: _FakeTemplatePath()
        )

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="Could not read"):
            updater._load_server_template()

    def test_non_mapping_template_raises(self, monkeypatch):
        updater = ConfigurationUpdater()
        monkeypatch.setattr(
            updater, "_yaml", lambda: _FakeYamlParser(load_result=["a", "b"])
        )
        with pytest.raises(ConfigUpdateError, match="must contain a mapping"):
            updater._load_server_template()


class TestLoadUserYaml:
    def test_malformed_yaml_wrapped_in_config_update_error(self, tmp_path):
        bad_yaml = tmp_path / "bad.yaml"
        bad_yaml.write_text("key: [unclosed\n")

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="not valid YAML"):
            updater._load_user_yaml(bad_yaml)

    def test_read_oserror_wrapped_in_config_update_error(
        self, tmp_path, monkeypatch
    ):
        yaml_file = tmp_path / "A.yaml"
        yaml_file.write_text("SERVER: {}\n")

        def failing_read_text(self, encoding=None):
            raise OSError("cannot read")

        monkeypatch.setattr(Path, "read_text", failing_read_text)

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="Could not read"):
            updater._load_user_yaml(yaml_file)

    def test_non_mapping_yaml_raises(self, tmp_path):
        list_yaml = tmp_path / "list.yaml"
        list_yaml.write_text("- a\n- b\n")

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="YAML mapping"):
            updater._load_user_yaml(list_yaml)


class TestIsInteractive:
    def test_returns_false_under_pytest(self):
        # stdin is never a real tty during a pytest run.
        assert ConfigurationUpdater._is_interactive() is False

    def test_prompt_program_exefolders_returns_empty_when_not_interactive(
        self,
    ):
        updater = ConfigurationUpdater()
        template_data = {"FAKE": {"EXEFOLDER": None}}
        prepared = _PreparedConfig(
            path=Path("A.yaml"),
            original_text="SERVER: {}\n",
            data={},
            template_data=template_data,
            missing_programs=("FAKE",),
        )
        assert updater._prompt_program_exefolders([prepared]) == {}


class TestApplyUpdatePlanExefolderBranches:
    def test_missing_program_section_without_exefolder_key(self):
        updater = ConfigurationUpdater()
        template_data = {"FAKE": {"LOCAL_RUN": True}}
        prepared = _PreparedConfig(
            path=Path("A.yaml"),
            original_text="SERVER: {}\n",
            data={},
            template_data=template_data,
            missing_programs=("FAKE",),
        )
        reports = updater._apply_update_plan([prepared], {})
        assert prepared.data["FAKE"] == {"LOCAL_RUN": True}
        assert reports[0].added_programs == ("FAKE",)

    def test_missing_program_section_with_non_none_exefolder_kept_as_is(
        self,
    ):
        updater = ConfigurationUpdater()
        template_data = {"FAKE": {"EXEFOLDER": "/opt/fake"}}
        prepared = _PreparedConfig(
            path=Path("A.yaml"),
            original_text="SERVER: {}\n",
            data={},
            template_data=template_data,
            missing_programs=("FAKE",),
        )
        updater._apply_update_plan([prepared], {})
        assert prepared.data["FAKE"]["EXEFOLDER"] == "/opt/fake"


class TestWriteChangedFilesAndAtomicWrite:
    def test_rejects_symlinked_server_yaml(self, tmp_path):
        real_file = tmp_path / "real.yaml"
        real_file.write_text("SERVER: {}\n")
        link = tmp_path / "link.yaml"
        link.symlink_to(real_file)

        updater = ConfigurationUpdater()
        prepared = _PreparedConfig(
            path=link,
            original_text="SERVER: {}\n",
            data={"SERVER": {}},
            template_data={},
        )
        report = ConfigUpdateReport(
            path=link, added_programs=("FAKE",), changed=True
        )
        with pytest.raises(ConfigUpdateError, match="symlinked"):
            updater._write_changed_files([prepared], [report])

    def test_dump_failure_wrapped_in_config_update_error(self, tmp_path):
        yaml_file = tmp_path / "A.yaml"
        yaml_file.write_text("SERVER: {}\n")

        updater = ConfigurationUpdater()

        def failing_dump_yaml(data):
            raise OSError("dump failed")

        updater._dump_yaml = failing_dump_yaml
        prepared = _PreparedConfig(
            path=yaml_file,
            original_text="SERVER: {}\n",
            data={"SERVER": {}},
            template_data={},
        )
        report = ConfigUpdateReport(
            path=yaml_file, added_programs=("FAKE",), changed=True
        )
        with pytest.raises(ConfigUpdateError, match="Could not update"):
            updater._write_changed_files([prepared], [report])

    def test_no_actual_text_change_skips_write(self, tmp_path):
        """A report marked changed=True but whose dumped YAML happens to
        equal the original text is not written out."""
        yaml_file = tmp_path / "A.yaml"
        original_text = "SERVER: {}\n"
        yaml_file.write_text(original_text)
        before = yaml_file.stat().st_mtime_ns

        updater = ConfigurationUpdater()
        updater._dump_yaml = lambda data: original_text
        prepared = _PreparedConfig(
            path=yaml_file,
            original_text=original_text,
            data={"SERVER": {}},
            template_data={},
        )
        report = ConfigUpdateReport(
            path=yaml_file, added_programs=("FAKE",), changed=True
        )
        updater._write_changed_files([prepared], [report])

        assert yaml_file.stat().st_mtime_ns == before
        assert yaml_file.read_text() == original_text

    def test_atomic_write_text_rejects_symlink(self, tmp_path):
        real_file = tmp_path / "real.yaml"
        real_file.write_text("content\n")
        link = tmp_path / "link.yaml"
        link.symlink_to(real_file)

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="symlinked"):
            updater._atomic_write_text(link, "new content\n")

    def test_atomic_write_text_cleans_up_tmp_file_on_oserror(
        self, tmp_path, monkeypatch
    ):
        target = tmp_path / "A.yaml"
        target.write_text("old content\n")

        def failing_replace(src, dst):
            raise OSError("replace failed")

        monkeypatch.setattr(update_module.os, "replace", failing_replace)

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="Could not update"):
            updater._atomic_write_text(target, "new content\n")

        # the temporary file created via tempfile.mkstemp must have been
        # cleaned up, leaving only the original target file behind.
        leftovers = [p for p in tmp_path.iterdir() if p.name != "A.yaml"]
        assert leftovers == []
        assert target.read_text() == "old content\n"

    def test_atomic_write_text_failure_before_tmp_file_created(
        self, tmp_path, monkeypatch
    ):
        """When the failure happens before tempfile.mkstemp() runs (e.g.
        path.stat() itself fails), tmp_name is still None, so the cleanup
        branch is skipped entirely."""
        target = tmp_path / "A.yaml"
        target.write_text("old content\n")

        real_stat = Path.stat

        def failing_stat(self, *args, **kwargs):
            # is_symlink() internally calls stat(follow_symlinks=False);
            # only fail the real "get file mode" call at line ~841.
            if kwargs.get("follow_symlinks", True) is False:
                return real_stat(self, *args, **kwargs)
            raise OSError("stat failed")

        monkeypatch.setattr(Path, "stat", failing_stat)

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="Could not update"):
            updater._atomic_write_text(target, "new content\n")

    def test_atomic_write_text_ignores_filenotfound_during_tmp_cleanup(
        self, tmp_path, monkeypatch
    ):
        """If the temp file created via mkstemp is already gone by the
        time the except handler tries to unlink it, that
        FileNotFoundError is swallowed and the original error still
        propagates."""
        target = tmp_path / "A.yaml"
        target.write_text("old content\n")

        def failing_replace(src, dst):
            os.unlink(src)
            raise OSError("replace failed")

        monkeypatch.setattr(update_module.os, "replace", failing_replace)

        updater = ConfigurationUpdater()
        with pytest.raises(ConfigUpdateError, match="Could not update"):
            updater._atomic_write_text(target, "new content\n")

        leftovers = [p for p in tmp_path.iterdir() if p.name != "A.yaml"]
        assert leftovers == []
        assert target.read_text() == "old content\n"


class TestFormatError:
    def test_oserror_with_strerror_returns_strerror(self):
        exc = OSError(2, "No such file or directory")
        assert (
            ConfigurationUpdater._format_error(exc)
            == "No such file or directory"
        )

    def test_non_oserror_returns_str(self):
        exc = ValueError("something else")
        assert ConfigurationUpdater._format_error(exc) == "something else"

    def test_oserror_without_strerror_returns_str(self):
        exc = OSError("plain message")
        exc.strerror = None
        assert ConfigurationUpdater._format_error(exc) == str(exc)


class TestConfigsCliConfigUpdateError:
    def test_config_update_error_becomes_click_exception(
        self, tmp_path, invoke_update_configs
    ):
        # No server yaml files exist at all under tmp_path, so
        # ConfigurationUpdater.update() raises ConfigUpdateError for the
        # missing server config directory, which the CLI must catch and
        # convert to a ClickException rather than letting it propagate
        # as a raw traceback.
        result = invoke_update_configs(tmp_path)

        assert result.exit_code != 0
        assert "Traceback" not in result.output
        assert "not found" in result.output


class TestConfigsCliModuleNotFoundError:
    def test_unrelated_module_not_found_error_propagates(
        self, tmp_path, monkeypatch, invoke_update_configs
    ):
        def raise_unrelated(self):
            raise ModuleNotFoundError(name="some_other_module")

        monkeypatch.setattr(
            update_module.ConfigurationUpdater, "update", raise_unrelated
        )

        result = invoke_update_configs(tmp_path)

        assert result.exit_code != 0
        assert isinstance(result.exception, ModuleNotFoundError)
        assert result.exception.name == "some_other_module"


class TestConfigsCliUpToDateReport:
    def test_already_up_to_date_server_prints_message(
        self,
        tmp_path,
        monkeypatch,
        write_server_yaml,
        invoke_update_configs,
    ):
        template = {"FAKE_PROGRAM": {"LOCAL_RUN": True}}
        monkeypatch.setattr(
            ConfigurationUpdater,
            "_load_server_template",
            lambda self: template,
        )
        write_server_yaml(
            tmp_path,
            "A.yaml",
            """
            SERVER:
                USER_VALUE: preserved
            FAKE_PROGRAM:
                LOCAL_RUN: true
            """,
        )

        result = invoke_update_configs(tmp_path, ["-s", "A"])

        assert result.exit_code == 0, result.output
        assert "A.yaml:" in result.output
        assert "already up to date" in result.output


class TestUpdateCLI:
    def test_deps_command_calls_update_pyproject_toml(self):
        runner = CliRunner()
        with patch(
            "chemsmart.cli.update.DependencyUpdater"
        ) as mock_updater_cls:
            mock_updater = MagicMock()
            mock_updater_cls.return_value = mock_updater
            result = runner.invoke(update, ["deps"], obj={})

        assert result.exit_code == 0, result.output
        mock_updater.update.assert_called_once()

    def test_version_command_requires_version_number(self):
        runner = CliRunner()
        result = runner.invoke(update, ["version"], obj={})
        assert result.exit_code != 0

    def test_version_command_calls_update_version_number(self):
        runner = CliRunner()
        with patch("chemsmart.cli.update.VersionUpdater") as mock_updater_cls:
            mock_updater = MagicMock()
            mock_updater_cls.return_value = mock_updater
            result = runner.invoke(update, ["version", "-v", "3.2.1"], obj={})

        assert result.exit_code == 0, result.output
        mock_updater_cls.assert_called_once_with("3.2.1")
        mock_updater.update.assert_called_once()
