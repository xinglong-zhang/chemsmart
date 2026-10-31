from unittest.mock import patch

import pytest

from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.settings.crest import (
    CRESTProjectSettings,
    CRESTProjectSettingsManager,
    YamlCRESTProjectSettings,
    YamlCRESTProjectSettingsBuilder,
)


class TestMainSettings:
    def test_main_settings_uses_gfn_version(self, my_crest_project_settings):
        settings = my_crest_project_settings.main_settings()
        assert isinstance(settings, CRESTJobSettings)
        assert settings.gfn_version == "gfnff"

    def test_default_gfn_version_is_gfn2(self):
        settings = CRESTProjectSettings().main_settings()
        assert settings.gfn_version == "gfn2"


class TestConformerSettings:
    def test_conformer_settings_sets_jobtype(self, my_crest_project_settings):
        settings = my_crest_project_settings.conformer_settings()
        assert settings.jobtype == "conformers"
        assert settings.gfn_version == "gfnff"

    def test_conformer_settings_is_independent_copy(
        self, my_crest_project_settings
    ):
        project = my_crest_project_settings
        settings_a = project.conformer_settings()
        settings_b = project.conformer_settings()
        settings_a.charge = 5
        assert settings_b.charge is None


class TestFromProject:
    def test_raises_when_project_not_found_anywhere(self):
        with pytest.raises(
            FileNotFoundError,
            match="No CREST project settings implemented",
        ):
            CRESTProjectSettings.from_project(
                "definitely_not_a_real_project_xyz"
            )

    def test_falls_back_to_packaged_project_when_user_absent(self):
        # "test" corresponds to the packaged
        # .chemsmart/crest/test.yaml template.
        with patch.object(
            CRESTProjectSettings,
            "_from_user_project_name",
            return_value=None,
        ):
            settings = CRESTProjectSettings.from_project("test")
        assert isinstance(settings, YamlCRESTProjectSettings)

    def test_uses_user_project_when_present(self):
        sentinel = object()
        with patch.object(
            CRESTProjectSettings,
            "_from_user_project_name",
            return_value=sentinel,
        ):
            result = CRESTProjectSettings.from_project("anything")
        assert result is sentinel

    def test_none_project_name_falls_through_to_not_found(self):
        with pytest.raises(FileNotFoundError):
            CRESTProjectSettings.from_project(None)


class TestFromProjectsManager:
    def test_returns_none_on_file_not_found(self):
        class _FailingManager:
            def create(self):
                raise FileNotFoundError("no such file")

        result = CRESTProjectSettings._from_projects_manager(_FailingManager())
        assert result is None

    def test_returns_manager_result_on_success(self):
        sentinel = object()

        class _SucceedingManager:
            def create(self):
                return sentinel

        result = CRESTProjectSettings._from_projects_manager(
            _SucceedingManager()
        )
        assert result is sentinel


class TestCRESTProjectSettingsManager:
    def test_rejects_none_filename(self):
        with pytest.raises(ValueError, match="filename is not specified"):
            CRESTProjectSettingsManager(filename=None)

    def test_stores_absolute_path(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        manager = CRESTProjectSettingsManager(filename="project.yaml")
        assert manager.filename == str(tmp_path / "project.yaml")

    def test_create_delegates_to_yaml_crest_project_settings(self, tmp_path):
        yaml_file = tmp_path / "myproj.yaml"
        yaml_file.write_text(
            "conformers:\n  gfn_version: gfnff\n  charge: -1\n"
        )
        manager = CRESTProjectSettingsManager(filename=str(yaml_file))
        settings = manager.create()
        assert isinstance(settings, YamlCRESTProjectSettings)
        assert settings.PROJECT_NAME == "myproj"


class TestYamlCRESTProjectSettingsBuilder:
    def test_build_parses_project_name_from_filename(self, tmp_path):
        yaml_file = tmp_path / "gas_solv.yaml"
        yaml_file.write_text("conformers:\n  gfn_version: gfn2\n")
        builder = YamlCRESTProjectSettingsBuilder(filename=str(yaml_file))
        project_settings = builder.build()
        assert project_settings.PROJECT_NAME == "gas_solv"

    def test_build_reads_conformer_settings_from_yaml(self, tmp_path):
        yaml_file = tmp_path / "proj.yaml"
        yaml_file.write_text(
            "conformers:\n  gfn_version: gfnff\n  energy_window: 4.0\n"
        )
        builder = YamlCRESTProjectSettingsBuilder(filename=str(yaml_file))
        project_settings = builder.build()
        settings = project_settings.conformer_settings()
        assert settings.gfn_version == "gfnff"
        assert settings.energy_window == 4.0
        assert settings.jobtype == "conformers"

    def test_empty_yaml_file_returns_empty_config(self, tmp_path):
        yaml_file = tmp_path / "empty.yaml"
        yaml_file.write_text("")
        builder = YamlCRESTProjectSettingsBuilder(filename=str(yaml_file))
        assert builder._read_config() == {}

    def test_missing_jobtype_section_warns_and_uses_defaults(
        self, tmp_path, caplog
    ):
        import logging

        yaml_file = tmp_path / "no_conformers.yaml"
        yaml_file.write_text("other_job:\n  gfn_version: gfnff\n")
        builder = YamlCRESTProjectSettingsBuilder(filename=str(yaml_file))
        with caplog.at_level(logging.WARNING):
            settings = builder._settings_for_job("conformers")
        assert "No CREST configuration found for conformers" in caplog.text
        assert settings.jobtype == "conformers"
        assert settings.gfn_version == "gfn2"

    def test_settings_for_job_sets_default_jobtype_when_absent_from_config(
        self, tmp_path
    ):
        yaml_file = tmp_path / "proj.yaml"
        yaml_file.write_text("conformers:\n  gfn_version: gfnff\n")
        builder = YamlCRESTProjectSettingsBuilder(filename=str(yaml_file))
        settings = builder._settings_for_job("conformers")
        assert settings.jobtype == "conformers"

    def test_config_explicit_jobtype_is_preserved(self, tmp_path):
        yaml_file = tmp_path / "proj.yaml"
        yaml_file.write_text("conformers:\n  jobtype: conformers\n")
        builder = YamlCRESTProjectSettingsBuilder(filename=str(yaml_file))
        settings = builder._settings_for_job("conformers")
        assert settings.jobtype == "conformers"


class TestYamlCRESTProjectSettings:
    def test_conformer_settings_returns_independent_copy(self):
        base = CRESTJobSettings(gfn_version="gfnff")
        project_settings = YamlCRESTProjectSettings(conformer_settings=base)
        settings_a = project_settings.conformer_settings()
        settings_b = project_settings.conformer_settings()
        assert settings_a is not settings_b
        settings_a.charge = 3
        assert settings_b.charge is None

    def test_from_yaml_builds_via_builder(self, tmp_path):
        yaml_file = tmp_path / "proj.yaml"
        yaml_file.write_text("conformers:\n  gfn_version: gfn2\n")
        project_settings = YamlCRESTProjectSettings.from_yaml(str(yaml_file))
        assert isinstance(project_settings, YamlCRESTProjectSettings)
        assert project_settings.PROJECT_NAME == "proj"


class TestCRESTJobSettings:
    def test_merge_dict(self):
        settings = CRESTJobSettings.default()
        settings.charge = None
        settings.multiplicity = None

        merged_settings = settings.merge({"charge": -2, "multiplicity": 3})
        assert merged_settings.charge == -2
        assert merged_settings.multiplicity == 3

    def test_get_settings_from_yaml(
        self, crest_yaml_settings_gas_project_name
    ):
        project_settings = CRESTProjectSettings.from_project(
            crest_yaml_settings_gas_project_name
        )
        conformers = project_settings.conformer_settings()

        assert isinstance(conformers, CRESTJobSettings)
        assert conformers.jobtype == "conformers"
        assert conformers.gfn_version == "gfn2"
        assert conformers.energy_window == 6.0
        assert conformers.optimization_level == "vtight"
        assert conformers.nci is False
        assert conformers.solvent_model is None
        assert conformers.solvent_id is None

    def test_get_settings_from_yaml_solv(
        self, crest_yaml_settings_solv_project_name
    ):
        project_settings = CRESTProjectSettings.from_project(
            crest_yaml_settings_solv_project_name
        )
        conformers = project_settings.conformer_settings()

        assert isinstance(conformers, CRESTJobSettings)
        assert conformers.jobtype == "conformers"
        assert conformers.gfn_version == "gfnff"
        assert conformers.optimization_level == "tight"
        assert conformers.solvent_model == "gbsa"
        assert conformers.solvent_id == "water"
