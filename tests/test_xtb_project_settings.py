"""
Direct unit tests for chemsmart.settings.xtb.XTBProjectSettings (the
plain, non-YAML base class) and its supporting manager/loader classes.

tests/test_xtb_cli.py exercises XTBProjectSettings.from_project() via
real YAML fixtures, but those always resolve to a YamlXTBProjectSettings
(whose sp/opt/hess_settings() override the base class entirely), so the
base class's own main_settings()/sp_settings()/opt_settings()/
hess_settings() computation and several from_project()/manager branches
were never directly exercised.
"""

from unittest.mock import patch

import pytest

from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.settings.xtb import (
    XTBProjectSettings,
    XTBProjectSettingsManager,
    YamlXTBProjectSettings,
    YamlXTBProjectSettingsBuilder,
)


class _MyProjectSettings(XTBProjectSettings):
    PROJECT_NAME = "my_project"
    gfn_version = "gfnff"


class TestMainSettings:
    def test_main_settings_uses_gfn_version(self):
        settings = _MyProjectSettings().main_settings()
        assert isinstance(settings, XTBJobSettings)
        assert settings.gfn_version == "gfnff"

    def test_default_gfn_version_is_gfn2(self):
        assert XTBProjectSettings().main_settings().gfn_version == "gfn2"


class TestJobtypeSettings:
    def test_sp_settings_disables_freq(self):
        settings = _MyProjectSettings().sp_settings()
        assert settings.jobtype == "sp"
        assert settings.freq is False
        assert settings.gfn_version == "gfnff"

    def test_opt_settings_enables_freq(self):
        settings = _MyProjectSettings().opt_settings()
        assert settings.jobtype == "opt"
        assert settings.freq is True

    def test_hess_settings_enables_freq(self):
        settings = _MyProjectSettings().hess_settings()
        assert settings.jobtype == "hess"
        assert settings.freq is True

    def test_each_jobtype_returns_independent_copy(self):
        project = _MyProjectSettings()
        a = project.sp_settings()
        b = project.sp_settings()
        a.charge = 5
        assert b.charge is None


class TestFromProject:
    def test_raises_when_project_not_found_anywhere(self):
        with pytest.raises(
            FileNotFoundError, match="No xTB project settings implemented"
        ):
            XTBProjectSettings.from_project(
                "definitely_not_a_real_project_xyz"
            )

    def test_falls_back_to_packaged_project_when_user_absent(self):
        # "test" corresponds to the packaged .chemsmart/xtb/test.yaml
        # template.
        with patch.object(
            XTBProjectSettings, "_from_user_project_name", return_value=None
        ):
            settings = XTBProjectSettings.from_project("test")
        assert isinstance(settings, YamlXTBProjectSettings)

    def test_uses_user_project_when_present(self):
        sentinel = object()
        with patch.object(
            XTBProjectSettings,
            "_from_user_project_name",
            return_value=sentinel,
        ):
            result = XTBProjectSettings.from_project("anything")
        assert result is sentinel

    def test_none_project_name_falls_through_to_not_found(self):
        with pytest.raises(FileNotFoundError):
            XTBProjectSettings.from_project(None)


class TestFromProjectsManager:
    def test_returns_none_on_file_not_found(self):
        class _FailingManager:
            def create(self):
                raise FileNotFoundError("no such file")

        result = XTBProjectSettings._from_projects_manager(_FailingManager())
        assert result is None

    def test_returns_manager_result_on_success(self):
        sentinel = object()

        class _SucceedingManager:
            def create(self):
                return sentinel

        result = XTBProjectSettings._from_projects_manager(
            _SucceedingManager()
        )
        assert result is sentinel


class TestXTBProjectSettingsManager:
    def test_rejects_none_filename(self):
        with pytest.raises(ValueError, match="filename is not specified"):
            XTBProjectSettingsManager(filename=None)

    def test_stores_absolute_path(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        manager = XTBProjectSettingsManager(filename="project.yaml")
        assert manager.filename == str(tmp_path / "project.yaml")

    def test_create_delegates_to_yaml_xtb_project_settings(self, tmp_path):
        yaml_file = tmp_path / "myproj.yaml"
        yaml_file.write_text("sp:\n  gfn_version: gfnff\n")
        manager = XTBProjectSettingsManager(filename=str(yaml_file))
        settings = manager.create()
        assert isinstance(settings, YamlXTBProjectSettings)
        assert settings.PROJECT_NAME == "myproj"


class TestYamlXTBProjectSettingsBuilder:
    def test_build_parses_project_name_from_filename(self, tmp_path):
        yaml_file = tmp_path / "gas_solv.yaml"
        yaml_file.write_text("sp:\n  gfn_version: gfn2\n")
        builder = YamlXTBProjectSettingsBuilder(filename=str(yaml_file))
        project_settings = builder.build()
        assert project_settings.PROJECT_NAME == "gas_solv"

    def test_build_reads_all_three_jobtypes_from_yaml(self, tmp_path):
        yaml_file = tmp_path / "proj.yaml"
        yaml_file.write_text(
            "sp:\n  gfn_version: gfnff\n"
            "opt:\n  optimization_level: tight\n"
            "hess:\n  charge: -1\n"
        )
        builder = YamlXTBProjectSettingsBuilder(filename=str(yaml_file))
        project_settings = builder.build()

        sp = project_settings.sp_settings()
        assert sp.jobtype == "sp"
        assert sp.gfn_version == "gfnff"
        assert sp.freq is False

        opt = project_settings.opt_settings()
        assert opt.jobtype == "opt"
        assert opt.optimization_level == "tight"
        assert opt.freq is True

        hess = project_settings.hess_settings()
        assert hess.jobtype == "hess"
        assert hess.charge == -1
        assert hess.freq is True

    def test_empty_yaml_file_returns_empty_config(self, tmp_path):
        yaml_file = tmp_path / "empty.yaml"
        yaml_file.write_text("")
        builder = YamlXTBProjectSettingsBuilder(filename=str(yaml_file))
        assert builder._read_config() == {}

    def test_missing_jobtype_section_uses_defaults(self, tmp_path):
        yaml_file = tmp_path / "sp_only.yaml"
        yaml_file.write_text("sp:\n  gfn_version: gfnff\n")
        builder = YamlXTBProjectSettingsBuilder(filename=str(yaml_file))
        opt_settings = builder._settings_for_job("opt")
        assert opt_settings.jobtype == "opt"
        assert opt_settings.gfn_version == "gfn2"
        assert opt_settings.freq is True

    def test_unsupported_jobtype_raises(self, tmp_path):
        yaml_file = tmp_path / "proj.yaml"
        yaml_file.write_text("sp:\n  gfn_version: gfn2\n")
        builder = YamlXTBProjectSettingsBuilder(filename=str(yaml_file))
        with pytest.raises(ValueError, match="Unsupported xTB job type"):
            builder._settings_for_job("bogus")


class TestYamlXTBProjectSettings:
    def test_each_jobtype_returns_independent_copy(self):
        base_sp = XTBJobSettings(jobtype="sp")
        base_opt = XTBJobSettings(jobtype="opt")
        base_hess = XTBJobSettings(jobtype="hess")
        project_settings = YamlXTBProjectSettings(
            sp_settings=base_sp, opt_settings=base_opt, hess_settings=base_hess
        )
        a = project_settings.sp_settings()
        b = project_settings.sp_settings()
        assert a is not b
        a.charge = 3
        assert b.charge is None

    def test_from_yaml_builds_via_builder(self, tmp_path):
        yaml_file = tmp_path / "proj.yaml"
        yaml_file.write_text("sp:\n  gfn_version: gfn2\n")
        project_settings = YamlXTBProjectSettings.from_yaml(str(yaml_file))
        assert isinstance(project_settings, YamlXTBProjectSettings)
        assert project_settings.PROJECT_NAME == "proj"
