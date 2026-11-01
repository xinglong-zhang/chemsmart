"""
Direct unit tests for chemsmart.settings.gaussian, covering code paths
that the extensive existing Gaussian test suite never reaches:

- GaussianProjectSettings' own default job-type-settings methods.
  Every real project in this codebase is YAML-backed
  (YamlGaussianProjectSettings), which overrides all of these methods,
  so the base implementation is otherwise never exercised.
- from_project's FileNotFoundError when a project exists in neither
  the user's nor the test-fixtures project directory.
- GaussianProjectSettingsManager's missing-filename validation.
- YamlGaussianProjectSettingsBuilder's KeyError -> RuntimeError
  wrapping in _project_settings_for_job (only reachable if
  Settings.from_dict itself raises KeyError, since dict.get() never
  does).
"""

from unittest.mock import patch

import pytest

from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
    GaussianQMMMJobSettings,
)
from chemsmart.settings.gaussian import (
    GaussianProjectSettings,
    GaussianProjectSettingsManager,
    YamlGaussianProjectSettingsBuilder,
)


class _ConcreteProjectSettings(GaussianProjectSettings):
    """A hardcoded (non-YAML) project settings subclass, exercising
    the base class's own default method implementations."""

    PROJECT_NAME = "concrete_test_project"
    functional = "b3lyp"
    small_basis = "sto-3g"
    large_basis = "def2tzvp"


@pytest.fixture()
def concrete_project():
    return _ConcreteProjectSettings()


class TestGaussianProjectSettingsDefaults:
    def test_main_settings_uses_functional_and_small_basis(
        self, concrete_project
    ):
        settings = concrete_project.main_settings()
        assert isinstance(settings, GaussianJobSettings)
        assert settings.functional == "b3lyp"
        assert settings.basis == "sto-3g"

    def test_opt_settings(self, concrete_project):
        settings = concrete_project.opt_settings()
        assert settings.jobtype == "opt"
        assert settings.functional == "b3lyp"
        assert settings.basis == "sto-3g"

    def test_modred_settings(self, concrete_project):
        settings = concrete_project.modred_settings()
        assert settings.jobtype == "modred"

    def test_ts_settings(self, concrete_project):
        settings = concrete_project.ts_settings()
        assert settings.jobtype == "ts"

    def test_irc_settings_returns_irc_specific_class_with_freq_disabled(
        self, concrete_project
    ):
        settings = concrete_project.irc_settings()
        assert isinstance(settings, GaussianIRCJobSettings)
        assert settings.jobtype == "irc"
        assert settings.freq is False

    def test_scan_settings_disables_freq(self, concrete_project):
        settings = concrete_project.scan_settings()
        assert settings.jobtype == "scan"
        assert settings.freq is False

    def test_nci_settings_disables_freq(self, concrete_project):
        settings = concrete_project.nci_settings()
        assert settings.jobtype == "nci"
        assert settings.freq is False

    def test_wbi_settings_disables_freq(self, concrete_project):
        settings = concrete_project.wbi_settings()
        assert settings.jobtype == "wbi"
        assert settings.freq is False

    def test_sp_settings_uses_large_basis_and_disables_freq(
        self, concrete_project
    ):
        settings = concrete_project.sp_settings()
        assert settings.jobtype == "sp"
        assert settings.freq is False
        assert settings.basis == "def2tzvp"

    def test_qmmm_settings_returns_qmmm_specific_class_with_freq_disabled(
        self, concrete_project
    ):
        settings = concrete_project.qmmm_settings()
        assert isinstance(settings, GaussianQMMMJobSettings)
        assert settings.jobtype == "qmmm"
        assert settings.freq is False


class TestFromProjectNotFound:
    def test_raises_filenotfounderror_when_project_missing_everywhere(self):
        with pytest.raises(
            FileNotFoundError, match="No project settings implemented"
        ):
            GaussianProjectSettings.from_project(
                "totally_nonexistent_project_xyz"
            )

    def test_from_chemsmart_test_projects_returns_none_when_missing(self):
        result = GaussianProjectSettings._from_chemsmart_test_projects(
            "totally_nonexistent_project_xyz"
        )
        assert result is None


class TestGaussianProjectSettingsManager:
    def test_none_filename_raises_value_error(self):
        with pytest.raises(ValueError, match="filename is not specified"):
            GaussianProjectSettingsManager(filename=None)


class TestYamlGaussianProjectSettingsBuilderKeyError:
    def test_from_dict_keyerror_wrapped_as_runtime_error(self):
        """dict.get() itself never raises KeyError for a missing
        jobtype, so the only way to reach this wrapping is for
        Settings.from_dict to raise KeyError on a present-but-malformed
        job config."""
        builder = YamlGaussianProjectSettingsBuilder(
            filename="tests/data/GaussianTests/project_yaml/gas_solv.yaml"
        )
        with patch.object(
            GaussianJobSettings, "from_dict", side_effect=KeyError("bad")
        ):
            with pytest.raises(RuntimeError, match="cannot be found"):
                builder._project_settings_for_job("opt")
