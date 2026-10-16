"""
Extended direct unit tests for :class:`chemsmart.cli.update.DependencyUpdater`
and :class:`chemsmart.cli.update.VersionUpdater`, and the ``update`` CLI
group, beyond what ``test_update.py`` already covers (package-name
normalization in ``_get_missing_dependencies``).

Every test that needs a pyproject.toml/environment.yml/VERSION/conf.py
redirects the updater instance's private path attributes to files
under ``tmp_path`` — never the real project files.
"""

import os
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


@pytest.fixture()
def updater(tmp_path):
    u = DependencyUpdater()
    u._package_path = tmp_path
    u._pyproject_path = tmp_path / "pyproject.toml"
    return u


@pytest.fixture()
def version_updater(tmp_path):
    u = VersionUpdater("0.0.0")
    u._pyproject_path = tmp_path / "pyproject.toml"
    u._version_file_path = tmp_path / "VERSION"
    u._docs_conf_file_path = tmp_path / "conf.py"
    return u


class TestUpdaterAbstractBase:
    def test_base_update_is_a_no_op_pass(self, updater):
        # Updater.update() is @abstractmethod with a `pass` body; concrete
        # subclasses always override it, so this only runs by calling the
        # base implementation directly against a subclass instance.
        assert Updater.update(updater) is None


class TestUpdaterPathProperties:
    def test_dependency_updater_properties_expose_private_attrs(
        self, updater, tmp_path
    ):
        assert updater.package_path == tmp_path
        assert updater.pyproject_path == tmp_path / "pyproject.toml"

    def test_version_updater_properties_expose_private_attrs(
        self, version_updater, tmp_path
    ):
        assert version_updater.pyproject_path == tmp_path / "pyproject.toml"
        assert version_updater.version_file_path == tmp_path / "VERSION"
        assert version_updater.docs_conf_file_path == tmp_path / "conf.py"


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


class TestDependencyUpdaterUpdateEntryPoint:
    def test_update_delegates_to_update_pyproject_toml(self, updater):
        with patch.object(
            updater, "update_pyproject_toml"
        ) as mock_update_pyproject:
            updater.update()
        mock_update_pyproject.assert_called_once_with()


class TestUpdateTomlNoVersionSpecifier:
    def test_dependency_without_version_specifier_normalized_lowercase(
        self, updater, tmp_path
    ):
        updater.pyproject_path.write_text("[project]\ndependencies = []\n")
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("SomePackage\n")

        updater._update_toml(req_file)

        new_text = updater.pyproject_path.read_text()
        assert "somepackage" in new_text


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
