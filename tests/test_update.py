import builtins
import tempfile
from pathlib import Path

from click.testing import CliRunner

import chemsmart.cli.update as update_module
from chemsmart.cli.update import (
    ConfigurationUpdater,
    DependencyUpdater,
    ProjectsUpdater,
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
