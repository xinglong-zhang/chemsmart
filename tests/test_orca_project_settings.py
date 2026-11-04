"""
Direct unit tests for chemsmart.settings.orca.ORCAProjectSettings (the
plain, non-YAML base class) and its supporting manager/loader classes.

tests/test_ORCAWriter.py exercises ORCAProjectSettings.from_project(),
but that always resolves to a YamlORCAProjectSettings, whose job-type
methods override the base class entirely. So the base class's own
main_settings()/opt_settings()/.../qmmm_settings() computations (driven
by the functional/small_basis/large_basis class attributes) were never
directly exercised.
"""

from unittest.mock import patch

import pytest

from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCANEBJobSettings,
    ORCAQMMMJobSettings,
    ORCATSJobSettings,
)
from chemsmart.settings.orca import (
    ORCAProjectSettings,
    ORCAProjectSettingsManager,
    YamlORCAProjectSettings,
)


class _MyProjectSettings(ORCAProjectSettings):
    PROJECT_NAME = "my_project"
    functional = "b3lyp"
    small_basis = "sto-3g"
    large_basis = "def2-tzvp"


class TestMainSettings:
    def test_main_settings_uses_functional_and_small_basis(self):
        settings = _MyProjectSettings().main_settings()
        assert settings.functional == "b3lyp"
        assert settings.basis == "sto-3g"


class TestOptAndModredSettings:
    def test_opt_settings(self):
        settings = _MyProjectSettings().opt_settings()
        assert settings.jobtype == "opt"
        assert settings.functional == "b3lyp"

    def test_modred_settings(self):
        settings = _MyProjectSettings().modred_settings()
        assert settings.jobtype == "modred"


class TestTsSettings:
    def test_ts_settings_returns_ts_job_settings_type(self):
        settings = _MyProjectSettings().ts_settings()
        assert isinstance(settings, ORCATSJobSettings)
        assert settings.jobtype == "ts"
        assert settings.functional == "b3lyp"


class TestIrcSettings:
    def test_irc_settings_disables_freq(self):
        settings = _MyProjectSettings().irc_settings()
        assert isinstance(settings, ORCAIRCJobSettings)
        assert settings.jobtype == "irc"
        assert settings.freq is False


class TestScanSettings:
    def test_scan_settings_disables_freq(self):
        settings = _MyProjectSettings().scan_settings()
        assert settings.jobtype == "scan"
        assert settings.freq is False


class TestNciSettings:
    def test_nci_settings_disables_freq(self):
        settings = _MyProjectSettings().nci_settings()
        assert settings.jobtype == "nci"
        assert settings.freq is False


class TestWbiSettings:
    def test_wbi_settings_disables_freq(self):
        settings = _MyProjectSettings().wbi_settings()
        assert settings.jobtype == "wbi"
        assert settings.freq is False


class TestSpSettings:
    def test_sp_settings_uses_large_basis_and_disables_freq(self):
        settings = _MyProjectSettings().sp_settings()
        assert settings.jobtype == "sp"
        assert settings.freq is False
        assert settings.basis == "def2-tzvp"


class TestNebSettings:
    def test_neb_settings_returns_neb_job_settings_type(self):
        settings = _MyProjectSettings().neb_settings()
        assert isinstance(settings, ORCANEBJobSettings)
        assert settings.jobtype == "neb"
        assert settings.freq is False


class TestQmmmSettings:
    def test_qmmm_settings_returns_qmmm_job_settings_type(self):
        settings = _MyProjectSettings().qmmm_settings()
        assert isinstance(settings, ORCAQMMMJobSettings)
        assert settings.jobtype == "qmmm"
        assert settings.freq is False


class TestFromProject:
    def test_raises_when_project_not_found_anywhere(self):
        with pytest.raises(
            FileNotFoundError, match="No project settings implemented"
        ):
            ORCAProjectSettings.from_project(
                "definitely_not_a_real_project_xyz"
            )

    def test_falls_back_to_chemsmart_test_projects_when_user_absent(self):
        # "orca" corresponds to tests/data/ORCATests/project_yaml/orca.yaml,
        # a bare project name that _from_chemsmart_test_projects resolves
        # (unlike the fixtures in test_ORCAWriter.py, which pass absolute
        # paths as the project name and so never exercise this fallback).
        with patch.object(
            ORCAProjectSettings,
            "_from_user_project_name",
            return_value=None,
        ):
            settings = ORCAProjectSettings.from_project("orca")
        assert isinstance(settings, YamlORCAProjectSettings)

    def test_uses_user_project_when_present(self):
        sentinel = object()
        with patch.object(
            ORCAProjectSettings,
            "_from_user_project_name",
            return_value=sentinel,
        ):
            result = ORCAProjectSettings.from_project("anything")
        assert result is sentinel


class TestFromProjectsManager:
    def test_returns_none_on_file_not_found(self):
        class _FailingManager:
            def create(self):
                raise FileNotFoundError("no such file")

        result = ORCAProjectSettings._from_projects_manager(_FailingManager())
        assert result is None

    def test_returns_manager_result_on_success(self):
        sentinel = object()

        class _SucceedingManager:
            def create(self):
                return sentinel

        result = ORCAProjectSettings._from_projects_manager(
            _SucceedingManager()
        )
        assert result is sentinel


class TestOrcaProjectSettingsManager:
    def test_rejects_none_filename(self):
        with pytest.raises(ValueError, match="filename is not specified"):
            ORCAProjectSettingsManager(filename=None)

    def test_stores_absolute_path(self, tmp_path):
        rel_path = "project.yaml"
        manager = ORCAProjectSettingsManager(filename=rel_path)
        assert manager.filename == __import__("os").path.abspath(rel_path)


class TestYamlORCAProjectSettingsGetters:
    def _make(self, **overrides):
        base = {
            "opt_settings": "opt",
            "modred_settings": "modred",
            "ts_settings": "ts",
            "irc_settings": "irc",
            "scan_settings": "scan",
            "nci_settings": "nci",
            "sp_settings": "sp",
            "td_settings": "td",
            "wbi_settings": "wbi",
            "qmmm_settings": "qmmm",
            "neb_settings": "neb",
        }
        base.update(overrides)
        return YamlORCAProjectSettings(**base)

    def test_all_simple_getters_return_stored_settings(self):
        settings = self._make()
        assert settings.irc_settings() == "irc"
        assert settings.scan_settings() == "scan"
        assert settings.nci_settings() == "nci"
        assert settings.sp_settings() == "sp"
        assert settings.td_settings() == "td"
        assert settings.wbi_settings() == "wbi"
        assert settings.qmmm_settings() == "qmmm"

    def test_neb_settings_falls_back_to_opt_when_none(self):
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        opt = ORCAJobSettings.default()
        opt.functional = "b3lyp"
        settings = self._make(opt_settings=opt, neb_settings=None)
        neb = settings.neb_settings()
        assert isinstance(neb, ORCANEBJobSettings)
        assert neb.functional == "b3lyp"
        # Cached on second call rather than rebuilt.
        assert settings.neb_settings() is neb


class TestProjectSettingsForJobKeyError:
    def test_key_error_from_settings_class_is_wrapped_in_runtime_error(self):
        from chemsmart.settings.orca import YamlORCAProjectSettingsBuilder

        builder = YamlORCAProjectSettingsBuilder(filename="fake.yaml")
        with (
            patch.object(
                builder,
                "_read_config",
                return_value={"opt": {"functional": "b3lyp"}},
            ),
            patch.object(
                __import__(
                    "chemsmart.jobs.orca.settings",
                    fromlist=["ORCAJobSettings"],
                ).ORCAJobSettings,
                "from_dict",
                side_effect=KeyError("missing_field"),
            ),
        ):
            with pytest.raises(
                RuntimeError, match="ORCA settings for job opt cannot"
            ):
                builder._project_settings_for_job(jobtype="opt")
