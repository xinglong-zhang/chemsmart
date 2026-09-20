import os

import pytest

from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.user import CHEMSMARTUserSettings


@pytest.fixture()
def chain_tests_directory():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_directory, "data", "ChainTests")


@pytest.fixture()
def chain1_yaml_stem(chain_tests_directory):
    return os.path.join(chain_tests_directory, "chain1")


@pytest.fixture()
def isolated_config_dir(tmp_path, monkeypatch):
    config_dir = tmp_path / "chemsmart_config"
    config_dir.mkdir()
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_dir))
    monkeypatch.setattr(
        CHEMSMARTUserSettings,
        "USER_CONFIG_DIR",
        str(config_dir),
    )
    return config_dir


class TestCHEMSMARTUserChainSettings:
    def test_user_chain_settings_dir(self, isolated_config_dir):
        settings = CHEMSMARTUserSettings()
        assert settings.user_chain_settings_dir == os.path.join(
            str(isolated_config_dir), "chain"
        )

    def test_all_available_chain_projects_empty(self, isolated_config_dir):
        settings = CHEMSMARTUserSettings()
        assert settings.all_available_chain_projects == []

    def test_all_available_chain_projects(self, isolated_config_dir):
        chain_dir = isolated_config_dir / "chain"
        chain_dir.mkdir()
        (chain_dir / "chain1.yaml").write_text("gaussian: test\n")
        (chain_dir / "test.yaml").write_text("orca: test\n")
        settings = CHEMSMARTUserSettings()
        assert sorted(settings.all_available_chain_projects) == [
            "chain1",
            "test",
        ]


class TestChainProjectSettings:
    def test_aliases_from_yaml(self, chain1_yaml_stem):
        settings = ChainProjectSettings.from_project(chain1_yaml_stem)
        assert settings.PROJECT_NAME == "chain1"
        assert settings.project_for("crest") == "crest_project1"
        assert settings.project_for("xtb") == "xtb_project1"
        assert settings.project_for("gaussian") == "gaussian_project2"
        assert settings.project_for("orca") == "orca_project3"

    def test_aliases_only_yaml(self, chain_tests_directory):
        settings = ChainProjectSettings.from_project(
            os.path.join(chain_tests_directory, "aliases_only")
        )
        assert settings.project_for("gaussian") == "gaussian_project2"

    def test_unknown_steps_key_errors(self, chain_tests_directory):
        with pytest.raises(
            ValueError, match="Unknown chain project keys: steps"
        ):
            ChainProjectSettings.from_project(
                os.path.join(chain_tests_directory, "missing_step_alias")
            )

    def test_unknown_top_level_keys_error(self, chain_tests_directory):
        with pytest.raises(
            ValueError, match="Unknown chain project keys: extra"
        ):
            ChainProjectSettings.from_project(
                os.path.join(chain_tests_directory, "unknown_key")
            )

    def test_missing_program_section_error(self, chain_tests_directory):
        settings = ChainProjectSettings.from_project(
            os.path.join(chain_tests_directory, "missing_program_section")
        )
        assert settings.project_for("crest") == "crest_project1"
        with pytest.raises(
            ValueError,
            match="has no gaussian section",
        ):
            settings.project_for("gaussian")

    def test_missing_target_project_yaml_uses_program_from_project(
        self, chain_tests_directory, isolated_config_dir
    ):
        settings = ChainProjectSettings.from_project(
            os.path.join(chain_tests_directory, "missing_target")
        )
        assert settings.project_for("gaussian") == "missing_gaussian_project"
        with pytest.raises(
            FileNotFoundError, match="missing_gaussian_project"
        ):
            GaussianProjectSettings.from_project(
                settings.project_for("gaussian")
            )

    def test_from_project_uses_packaged_chain1(self, isolated_config_dir):
        settings = ChainProjectSettings.from_project("chain1")
        assert settings.PROJECT_NAME == "chain1"
        assert settings.project_for("crest") == "test"
        assert settings.project_for("xtb") == "test"
        assert settings.project_for("gaussian") == "test"
        assert settings.project_for("orca") == "test"

    def test_from_project_uses_packaged_test(self, isolated_config_dir):
        settings = ChainProjectSettings.from_project("test")
        assert settings.PROJECT_NAME == "test"
        assert settings.project_for("gaussian") == "test"
        with pytest.raises(ValueError, match="has no crest section"):
            settings.project_for("crest")

    def test_from_project_user_dir_overrides_packaged(
        self, isolated_config_dir
    ):
        chain_dir = isolated_config_dir / "chain"
        chain_dir.mkdir()
        (chain_dir / "chain1.yaml").write_text("gaussian: user_gaussian\n")
        settings = ChainProjectSettings.from_project("chain1")
        assert settings.project_for("gaussian") == "user_gaussian"
        with pytest.raises(ValueError, match="has no crest section"):
            settings.project_for("crest")

    def test_from_project_loads_user_chain_test(self, isolated_config_dir):
        chain_dir = isolated_config_dir / "chain"
        chain_dir.mkdir()
        (chain_dir / "test.yaml").write_text("gaussian: user_gaussian\n")
        settings = ChainProjectSettings.from_project("test")
        assert settings.PROJECT_NAME == "test"
        assert settings.project_for("gaussian") == "user_gaussian"

    def test_from_project_missing_chain_file(self, isolated_config_dir):
        with pytest.raises(
            FileNotFoundError,
            match="No chain project settings implemented for missing",
        ):
            ChainProjectSettings.from_project("missing")

    def test_unknown_program_error(self, chain1_yaml_stem):
        settings = ChainProjectSettings.from_project(chain1_yaml_stem)
        with pytest.raises(ValueError, match="Unknown chain program 'foo'"):
            settings.project_for("foo")
