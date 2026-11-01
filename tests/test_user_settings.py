"""Tests for chemsmart.settings.user.CHEMSMARTUserSettings."""

from chemsmart.settings.user import CHEMSMARTUserSettings
from tests.helpers import write_text_file


class TestResolveConfigDir:
    def test_uses_env_override_when_set(self, tmp_path, monkeypatch):
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
        assert CHEMSMARTUserSettings.resolve_config_dir() == str(tmp_path)

    def test_falls_back_to_default_when_unset(self, monkeypatch):
        monkeypatch.delenv("CHEMSMART_CONFIG_DIR", raising=False)
        assert (
            CHEMSMARTUserSettings.resolve_config_dir()
            == CHEMSMARTUserSettings.USER_CONFIG_DIR
        )


class TestInit:
    def test_empty_data_when_yaml_missing(self, user_settings):
        assert user_settings.data == {}

    def test_loads_yaml_data_when_present(
        self, tmp_path, monkeypatch, write_tmp_file
    ):
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
        write_tmp_file(
            CHEMSMARTUserSettings.USER_YAML_FILE,
            "SCRATCH: /tmp/scratch\nEMAIL: user@example.com\n",
        )

        settings = CHEMSMARTUserSettings()
        assert settings.data["SCRATCH"] == "/tmp/scratch"
        assert settings.data["EMAIL"] == "user@example.com"


class TestPathProperties:
    def test_user_server_dir(self, user_settings, tmp_path):
        assert user_settings.user_server_dir == str(tmp_path / "server")

    def test_user_gaussian_settings_dir(self, user_settings, tmp_path):
        assert user_settings.user_gaussian_settings_dir == str(
            tmp_path / "gaussian"
        )

    def test_user_gaussian_envars(self, user_settings, tmp_path):
        assert user_settings.user_gaussian_envars == str(
            tmp_path / "gaussian" / "gaussian.envars"
        )

    def test_user_gaussian_modules(self, user_settings, tmp_path):
        assert user_settings.user_gaussian_modules == str(
            tmp_path / "gaussian" / "gaussian.modules"
        )

    def test_user_gaussian_script(self, user_settings, tmp_path):
        assert user_settings.user_gaussian_script == str(
            tmp_path / "gaussian" / "gaussian.sh"
        )

    def test_user_orca_settings_dir(self, user_settings, tmp_path):
        assert user_settings.user_orca_settings_dir == str(tmp_path / "orca")

    def test_user_orca_envars(self, user_settings, tmp_path):
        assert user_settings.user_orca_envars == str(
            tmp_path / "orca" / "orca.envars"
        )

    def test_user_orca_modules(self, user_settings, tmp_path):
        assert user_settings.user_orca_modules == str(
            tmp_path / "orca" / "orca.modules"
        )

    def test_user_orca_script(self, user_settings, tmp_path):
        assert user_settings.user_orca_script == str(
            tmp_path / "orca" / "orca.sh"
        )


class TestDataBackedCachedProperties:
    def test_scratch_present(self, tmp_path, monkeypatch, write_tmp_file):
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
        write_tmp_file(
            CHEMSMARTUserSettings.USER_YAML_FILE, "SCRATCH: /tmp/scratch\n"
        )
        settings = CHEMSMARTUserSettings()
        assert settings.scratch == "/tmp/scratch"

    def test_scratch_absent_returns_none(self, user_settings):
        assert user_settings.scratch is None

    def test_email_present(self, tmp_path, monkeypatch, write_tmp_file):
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
        write_tmp_file(
            CHEMSMARTUserSettings.USER_YAML_FILE, "EMAIL: user@example.com\n"
        )
        settings = CHEMSMARTUserSettings()
        assert settings.email == "user@example.com"

    def test_email_absent_returns_none(self, user_settings):
        assert user_settings.email is None

    def test_project_present(self, tmp_path, monkeypatch, write_tmp_file):
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
        write_tmp_file(
            CHEMSMARTUserSettings.USER_YAML_FILE, "PROJECT: myproject\n"
        )
        settings = CHEMSMARTUserSettings()
        assert settings.project == "myproject"

    def test_project_absent_returns_none(self, user_settings):
        assert user_settings.project is None


class TestGlobBasedLists:
    def test_server_yaml_files(self, user_settings, tmp_path):
        server_dir = tmp_path / "server"
        write_text_file(server_dir / "cluster1.yaml", "")
        write_text_file(server_dir / "cluster2.yaml", "")
        write_text_file(server_dir / "notes.txt", "")

        assert sorted(user_settings.server_yaml_files) == sorted(
            [
                str(server_dir / "cluster1.yaml"),
                str(server_dir / "cluster2.yaml"),
            ]
        )

    def test_server_yaml_files_empty_when_dir_missing(self, user_settings):
        assert user_settings.server_yaml_files == []

    def test_gaussian_project_yaml_files(self, user_settings, tmp_path):
        gaussian_dir = tmp_path / "gaussian"
        write_text_file(gaussian_dir / "b3lyp.yaml", "")
        write_text_file(gaussian_dir / "m062x.yaml", "")

        assert sorted(user_settings.gaussian_project_yaml_files) == sorted(
            [
                str(gaussian_dir / "b3lyp.yaml"),
                str(gaussian_dir / "m062x.yaml"),
            ]
        )

    def test_orca_project_yaml_files(self, user_settings, tmp_path):
        orca_dir = tmp_path / "orca"
        write_text_file(orca_dir / "dlpno.yaml", "")

        assert user_settings.orca_project_yaml_files == [
            str(orca_dir / "dlpno.yaml")
        ]

    def test_all_available_servers(self, user_settings, tmp_path):
        server_dir = tmp_path / "server"
        write_text_file(server_dir / "cluster1.yaml", "")

        assert user_settings.all_available_servers == ["cluster1"]

    def test_all_available_gaussian_projects(self, user_settings, tmp_path):
        gaussian_dir = tmp_path / "gaussian"
        write_text_file(gaussian_dir / "m062x.yaml", "")

        assert user_settings.all_available_gaussian_projects == ["m062x"]

    def test_all_available_orca_projects(self, user_settings, tmp_path):
        orca_dir = tmp_path / "orca"
        write_text_file(orca_dir / "dlpno.yaml", "")

        assert user_settings.all_available_orca_projects == ["dlpno"]
